# FORCE REDEPLOY: v2.0 (Phoenix Protocol)
import modal
import os
import re
import uuid
from fastapi import FastAPI, APIRouter, HTTPException
from pydantic import BaseModel
from loguru import logger

# --- Image Definition ---
# This is the battle-tested image definition from the proven Zeta Oracle.
# It correctly handles all dependencies and build complexities for Evo2.
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11"
    )
    .apt_install(
        ["build-essential", "cmake", "ninja-build",
            "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"]
    )
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
    })
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic")
)

# --- App Definition ---
app = modal.App("evo-service", image=evo2_image)

# --- Shared State for Jobs ---
job_status_dict = modal.Dict.from_name("evo-job-status", create_if_missing=True)

# --- Persistent Volume for Caching ---
# This ensures the model is downloaded only once.
volume = modal.Volume.from_name("evo-model-cache", create_if_missing=True)
EVO_MODEL_DIR = "/root/.cache/huggingface"

# --- Pydantic Models ---
class GuideRNARequest(BaseModel):
    target_sequence: str
    num_guides: int = 5

class JobSubmitResponse(BaseModel):
    job_id: str

# --- Main Service Class ---
@app.cls(
    gpu="H100:2", 
    volumes={EVO_MODEL_DIR: volume}, 
    scaledown_window=300, 
    timeout=1800
)
class EvoService:
    @modal.enter()
    def load_model(self):
        """Load the Evo2 40B model into memory when the container starts."""
        from evo2 import Evo2
        logger.info("🚀 Loading Evo2 40B model for guide RNA generation...")
        # The HF_HOME env var is set by Modal from the volume mount.
        # Evo2 will automatically use this directory.
        self.model = Evo2('evo2_40b')
        logger.info("🎉 Evo2 40B model loaded successfully! WE ARE OPERATIONAL. 💥")

    @modal.method()
    def generate_in_background(self, job_id: str, target_sequence: str, num_guides: int):
        """This runs the actual model inference in the background."""
        logger.info(f"Job {job_id}: Starting background generation...")
        job_status_dict[job_id] = {"status": "running"}
        try:
            # The high-level .generate() method handles tokenization.
            # The result is an object with a .sequences attribute.
            generation_result = self.model.generate(
                prompt_seqs=[target_sequence],
                n_tokens=num_guides * 24, # A rough estimate
            )
            decoded_guides = generation_result.sequences
            logger.info(f"Job {job_id}: Successfully generated {len(decoded_guides)} candidates.")
            job_status_dict[job_id] = {"status": "complete", "result": decoded_guides}
        except Exception as e:
            logger.error(f"Job {job_id}: Error during guide RNA generation: {e}")
            job_status_dict[job_id] = {"status": "failed", "error": str(e)}

# --- FastAPI Endpoints ---
fastapi_app = FastAPI(title="Evo2 Guide RNA Generation Service")

@fastapi_app.post("/generate", status_code=202, response_model=JobSubmitResponse)
def submit_generation(request: GuideRNARequest):
    """Accepts a sequence and submits a background job to generate guide RNAs."""
    job_id = str(uuid.uuid4())
    logger.info(f"Received submission request. Assigning job ID: {job_id}")
    job_status_dict[job_id] = {"status": "pending"}
    EvoService().generate_in_background.spawn(job_id, request.target_sequence, request.num_guides)
    return {"job_id": job_id}

@fastapi_app.get("/status/{job_id}")
def get_status(job_id: str):
    """Polls for the status and result of a guide RNA generation job."""
    logger.info(f"Received status query for job ID: {job_id}")
    if job_id not in job_status_dict:
        raise HTTPException(status_code=404, detail="Job ID not found.")
    return job_status_dict[job_id]

@app.asgi_app()
def api():
    return fastapi_app

@app.local_entrypoint()
def main():
    """A local test function for smoke testing."""
    print("--- ⚔️ LOCAL EVO-SERVICE TEST ⚔️ ---")
    model = EvoService()
    
    # This is a mock test, as the async part doesn't run locally this way.
    # It mainly serves to check for syntax errors.
    print("Local entrypoint for syntax validation. Does not run generation.")
    print("✅ Syntax validation passed.")