import modal
import os
import uuid
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from loguru import logger
import httpx
import asyncio
from Bio.Seq import Seq
from Bio.Data import CodonTable

# --- Guided Warfare Doctrine Constants ---
# BEAM_WIDTH = 5  # Obsolete under new doctrine
# CHUNK_SIZE = 30 # Obsolete
# MAX_LENGTH = 120 # Obsolete

# --- Service URLs (Kill-Chain Coordinates) ---
# As per the Zeta Engineering Doctrine, we use direct HTTP strikes, not the false god of `lookup()`.
ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run"
EVO_SERVICE_URL = "https://crispro--evo-service-evoservice-api.modal.run"

# --- Image Definition ---
# DOCTRINE: This is the battle-tested image definition from the proven Zeta Oracle.
# It correctly handles all dependencies and build complexities for Evo2.
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.12"
    )
    .apt_install(
        ["build-essential", "cmake", "ninja-build",
            "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++", "gfortran", "llvm"]
    )
    # Create a symlink to satisfy the hardcoded path in the numpy build script.
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/llvm-ar /tools/llvm/bin/llvm-ar",
        "ln -s /usr/bin/llvm-ranlib /tools/llvm/bin/llvm-ranlib"
    )
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
    })
    .run_commands(
        "git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .",
        "pip uninstall -y transformer-engine transformer_engine",
        "pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation"
    )
    .pip_install(
        "fastapi",
        "uvicorn[standard]",
        "loguru",
        "pydantic",
        "httpx",
        "biopython",
        "scikit-learn>=1.5.0",
        "flash-attn"
    )
)

# --- App Definition ---
app = modal.App(
    "zeta-forge-v1", image=evo2_image, secrets=[modal.Secret.from_dotenv()]
)

# --- SINGLE, SECURED FASTAPI APP ---
# All endpoints will be attached to this single, globally-defined app.
fastapi_app = FastAPI(
    title="Zeta Forge",
    description="Generates novel biologic inhibitors using the Evo2 model.",
)

origins = [
    "http://localhost:5173",  # Local development frontend
    "http://localhost:5174",
    "https://crispro.io",  # Production domain
]

fastapi_app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# --- Pydantic Models (Data Contracts) ---
class ForgeRequest(BaseModel):
    # --- DOCTRINE: THE AMBUSH ---
    # The API is now simplified to only support the Ambush doctrine.
    bait_sequence: str = Field(..., description="The reverse complement of the critical target motif for DNA-based inhibitor generation.")
    
    # --- Shared Generation Parameters ---
    num_candidates_per_temp: int = Field(5, description="Number of inhibitors to generate for each temperature setting.")
    temperatures: list[float] = Field([0.2, 0.7, 1.2], description="List of temperature settings for generation to balance fidelity and creativity.")
    generation_length: int = Field(120, description="The desired length for the generated inhibitor sequence (amino acids for protein, nucleotides for DNA).")

class JobSubmitResponse(BaseModel):
    job_id: str

class ForgeStatusResponse(BaseModel):
    status: str
    message: str | None = None 
    # --- DOCTRINE: REVISED BATTLE REPORT ---
    # The Forge now only returns a list of raw DNA candidates.
    # The CommandCenter is responsible for validation and ranking.
    candidates: list[str] = []
    error: str | None = None

# --- Shared State ---
# A modal.Dict is used to share state (job status) between the submission endpoint
# and the background function running the computationally expensive forge loop.
job_status_dict = modal.Dict.from_name("zeta-forge-job-status", create_if_missing=True)

# --- Global FastAPI Endpoints ---
@fastapi_app.post("/generate_inhibitor", status_code=202, response_model=JobSubmitResponse)
def submit_forge_job(request: ForgeRequest):
    """
    Submits a job to the forge. Immediately returns a job ID.
    The actual work is spawned in the background.
    """
    job_id = str(uuid.uuid4())
    logger.info(f"Received forge submission. Assigning job ID: {job_id}")
    job_status_dict[job_id] = {"status": "pending"}
    
    # Spawn the long-running process in the background.
    run_forge_in_background.spawn(job_id, request)
    
    return JobSubmitResponse(job_id=job_id)

@fastapi_app.get("/status/{job_id}", response_model=ForgeStatusResponse)
def get_forge_status(job_id: str):
    """Poll this endpoint to get the status of a forge job."""
    logger.info(f"Received status query for job ID: {job_id}")
    if job_id not in job_status_dict:
        raise HTTPException(status_code=404, detail="Job ID not found.")
    
    status_data = job_status_dict.get(job_id, {})
    return ForgeStatusResponse(**status_data)

# --- Mount the global FastAPI app ---
@app.function()
@modal.asgi_app()
def api():
    """Serve the FastAPI app."""
    return fastapi_app


# --- Helper and Background Functions ---
async def _make_request(client, method, url, **kwargs):
    """A helper to make async HTTP requests with error handling."""
    try:
        response = await client.request(method, url, **kwargs)
        response.raise_for_status()
        return response.json()
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error for {e.request.url}: {e.response.status_code} - {e.response.text}")
        return {"error": f"HTTP {e.response.status_code}", "details": e.response.text}
    except httpx.RequestError as e:
        logger.error(f"Request error for {e.request.url}: {e}")
        return {"error": "Request failed", "details": str(e)}

def translate_dna_to_protein(dna_sequence: str) -> str:
    """
    Translates a DNA sequence into a protein sequence, trying multiple genetic
    codes to find the most plausible translation (fewest stop codons).
    """
    try:
        # Ensure the sequence length is a multiple of 3
        remainder = len(dna_sequence) % 3
        if remainder != 0:
            dna_sequence = dna_sequence[:-remainder]
            
        seq = Seq(dna_sequence)
        
        best_translation = ""
        min_stop_codons = float('inf')

        # --- DOCTRINE: Multi-Code Translation (Corrected) ---
        # Iterate through the dictionary of available unambiguous DNA codon tables.
        for table_id in CodonTable.unambiguous_dna_by_id:
            try:
                protein_sequence = str(seq.translate(table=table_id, to_stop=False))
                
                stop_codon_count = protein_sequence.count('*')
                
                if stop_codon_count < min_stop_codons:
                    min_stop_codons = stop_codon_count
                    best_translation = protein_sequence
                
                # If we find a perfect translation (no stops), we can stop early.
                if min_stop_codons == 0:
                    break
            except Exception:
                # Some sequences might not be translatable with some tables, which is fine.
                continue

        if not best_translation:
             # Fallback to standard table if no other translation worked
             return str(seq.translate(to_stop=False))

        return best_translation

    except Exception as e:
        logger.error(f"Biopython translation failed: {e}")
        return "TRANSLATION_ERROR"

@app.function(
    gpu="H100:2",
    timeout=3600,
    image=evo2_image, # Ensure the function uses our defined image
)
async def run_forge_in_background(job_id: str, request: ForgeRequest):
    """
    This is the core of the Zeta Forge, re-aligned with the Ambush doctrine.
    """
    try:
        from evo2 import Evo2
        import torch

        logger.info(f"Job {job_id}: Executing AMBUSH (DNA) doctrine.")
        # Format: "PREPEND_TOKENS<|completion|>APPEND_TOKENS"
        # We provide the bait as the prepend sequence to guide generation.
        prompt = f"{request.bait_sequence}<|completion|>"
        job_status_dict[job_id] = {
            "status": "running", 
            "message": f"Competitive Generation started. Generating {request.num_candidates_per_temp * len(request.temperatures)} total DNA candidates.",
            "candidates": [],
        }

        # Load the Evo2 model onto the GPU
        logger.info(f"Job {job_id}: Loading Evo2 40B model...")
        device = "cuda"
        evo_model = Evo2('evo2_40b')
        logger.info(f"Job {job_id}: Evo2 model loaded successfully.")

        all_generated_candidates = set()
        
        # --- DOCTRINE: DIRECT, IN-PROCESS GENERATION ---
        for temp in request.temperatures:
            loop_msg = f"Generating {request.num_candidates_per_temp} candidates with temperature {temp}..."
            logger.info(f"Job {job_id}: {loop_msg}")
            current_status = job_status_dict.get(job_id, {})
            current_status["message"] = loop_msg
            job_status_dict[job_id] = current_status

            # Directly generate candidates using the loaded model
            generation_output = evo_model.generate(
                [prompt],
                n_samples=request.num_candidates_per_temp,
                max_new_tokens=request.generation_length,
                temperature=temp,
                top_k=4,
                top_p=0.9,
                return_prompt=False,
            )
            
            if generation_output and generation_output[0]:
                 all_generated_candidates.update(generation_output[0])
            else:
                 logger.warning(f"Job {job_id}: Generation at temp {temp} produced no output.")

        if not all_generated_candidates:
            raise RuntimeError("The Zeta Forge failed to generate any candidates after all attempts.")
        
        unique_candidates = list(all_generated_candidates)
        logger.success(f"Job {job_id}: Competitive Generation complete. Generated {len(unique_candidates)} unique candidates.")

        final_status = {
            "status": "complete",
            "message": f"Competitive Generation complete. Returning {len(unique_candidates)} unique candidates.",
            "candidates": unique_candidates,
        }
        job_status_dict[job_id] = final_status

    except Exception as e:
        logger.error(f"Job {job_id}: A critical error occurred in the forge: {e}")
        import traceback
        logger.error(traceback.format_exc())
        job_status_dict[job_id] = {"status": "failed", "message": str(e), "candidates": []} 