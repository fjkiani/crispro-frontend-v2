# FORCE REDEPLOY: v2.0 (Phoenix Protocol)
import modal
import os
import uuid
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from loguru import logger

# --- Image Definition ---
# This is the battle-tested image definition from the proven Zeta Oracle.
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
    # --- Fix for numpy 1.22.0 build ---
    # As per Zeta Doctrine (evo2-deployment-guide.mdc), this is a critical fix
    # to address a build bug in this specific numpy version.
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar",
    )
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0")
)

# --- App Definition ---
app = modal.App("evo-service", image=evo2_image)

# --- Shared State & Volume ---
job_status_dict = modal.Dict.from_name("evo-job-status", create_if_missing=True)
volume = modal.Volume.from_name("evo-model-cache", create_if_missing=True)
EVO_MODEL_DIR = "/root/.cache/huggingface"

# --- Pydantic Models ---
class JobSubmitResponse(BaseModel):
    job_id: str

class GenerationRequest(BaseModel):
    prompt: str = Field(..., description="A detailed natural language prompt for the generative model.")
    n_tokens: int = Field(100, description="The number of tokens to generate.")

# --- Main Service Class ---
@app.cls(
    gpu="H100:2", 
    volumes={EVO_MODEL_DIR: volume}, 
    scaledown_window=300, 
    timeout=1800
)
class EvoService:
    @modal.enter()
    def load_model_and_api(self):
        """Load model and initialize the FastAPI app within the container."""
        from evo2 import Evo2
        logger.info("🚀 Loading Evo2 40B model for GENERAL PURPOSE generation...")
        self.model = Evo2('evo2_40b')
        logger.info("🎉 Evo2 40B model loaded successfully!")

        self.fastapi_app = FastAPI(title="Evo2 General-Purpose Generation Service")

        @self.fastapi_app.post("/generate", status_code=202, response_model=JobSubmitResponse)
        def submit_generation(request: GenerationRequest):
            job_id = str(uuid.uuid4())
            logger.info(f"Received submission request. Assigning job ID: {job_id}")
            job_status_dict[job_id] = {"status": "pending"}
            self.generate_in_background.spawn(job_id, request.prompt, request.n_tokens)
            return {"job_id": job_id}

        @self.fastapi_app.get("/status/{job_id}")
        def get_status(job_id: str):
            logger.info(f"Received status query for job ID: {job_id}")
            if job_id not in job_status_dict:
                raise HTTPException(status_code=404, detail="Job ID not found.")
            return job_status_dict[job_id]

    @modal.method()
    def generate_in_background(self, job_id: str, prompt: str, n_tokens: int):
        """Runs the model inference in the background with a general-purpose prompt."""
        logger.info(f"Job {job_id}: Starting background generation with prompt: '{prompt[:50]}...'")
        job_status_dict[job_id] = {"status": "running"}
        try:
            generation_result = self.model.generate(
                prompt_seqs=[prompt],
                n_tokens=n_tokens,
            )
            # We now return a structured result, not a simple list.
            # This is critical for reliable parsing by the client (the Zeta Forge).
            decoded_sequence = generation_result.sequences[0]
            logger.info(f"Job {job_id}: Successfully generated sequence.")
            job_status_dict[job_id] = {"status": "complete", "result": {"sequence": decoded_sequence}}
        except Exception as e:
            logger.error(f"Job {job_id}: Error during generation: {e}")
            job_status_dict[job_id] = {"status": "failed", "error": str(e)}

    @modal.asgi_app()
    def api(self):
        """Serve the FastAPI app."""
        return self.fastapi_app

# --- Local Entrypoint for Testing ---
@app.local_entrypoint()
def main():
    print("--- ⚔️ LOCAL EVO-SERVICE TEST ⚔️ ---")
    print("Local entrypoint for syntax validation. Does not run generation.")
    print("✅ Syntax validation passed.")