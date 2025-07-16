# FORCE REDEPLOY: v2.0 (Phoenix Protocol)
import modal
import os
import uuid
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
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
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install(
        "fastapi", 
        "uvicorn[standard]", 
        "loguru", 
        "pydantic"
    )
    # Force the correct numpy version AFTER all other dependencies are installed
    # This overwrites any incompatible version brought in by evo2's setup.py
    .run_commands("pip install numpy==1.23.5 --force-reinstall")
)

# --- App Definition ---
app = modal.App("evo-service", image=evo2_image)

# --- Shared State & Volume ---
job_status_dict = modal.Dict.from_name("evo-job-status", create_if_missing=True)
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
    def __init__(self):
        """
        Initialize the FastAPI app and define the routes.
        This ensures the routes are part of the class instance.
        """
        self.fastapi_app = FastAPI(title="Evo2 Guide RNA Generation Service")
        self.fastapi_app.add_api_route(
            "/generate", self.submit_generation, methods=["POST"], status_code=202, response_model=JobSubmitResponse
        )
        self.fastapi_app.add_api_route("/status/{job_id}", self.get_status, methods=["GET"])

    @modal.enter()
    def load_model(self):
        """Load the model into memory."""
        from evo2 import Evo2
        logger.info("🚀 Loading Evo2 40B model for guide RNA generation...")
        self.model = Evo2('evo2_40b')
        logger.info("🎉 Evo2 40B model loaded successfully!")

    def submit_generation(self, request: GuideRNARequest):
        """Endpoint to submit a guide generation job."""
        job_id = str(uuid.uuid4())
        logger.info(f"Received submission request. Assigning job ID: {job_id}")
        job_status_dict[job_id] = {"status": "pending"}
        # We call the background generation method using .spawn()
        self.generate_in_background.spawn(job_id, request.target_sequence, request.num_guides)
        return {"job_id": job_id}

    def get_status(self, job_id: str):
        """Endpoint to get the status of a job."""
        logger.info(f"Received status query for job ID: {job_id}")
        if job_id not in job_status_dict:
            raise HTTPException(status_code=404, detail="Job ID not found.")
        return job_status_dict[job_id]

    @modal.method()
    def generate_in_background(self, job_id: str, target_sequence: str, num_guides: int):
        """Runs the model inference in the background."""
        logger.info(f"Job {job_id}: Starting background generation...")
        job_status_dict[job_id] = {"status": "running"}
        try:
            generation_result = self.model.generate(
                prompt_seqs=[target_sequence],
                n_tokens=num_guides * 24,
            )
            decoded_guides = generation_result.sequences
            logger.info(f"Job {job_id}: Successfully generated {len(decoded_guides)} candidates.")
            # Align status and key with what the command_center expects
            job_status_dict[job_id] = {"status": "completed", "guides": decoded_guides}
        except Exception as e:
            logger.error(f"Job {job_id}: Error during guide RNA generation: {e}")
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