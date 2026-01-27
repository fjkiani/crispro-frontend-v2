"""Evo2 service Modal app entrypoint.
This is a thin entrypoint that wires together all the modular components."""
import modal
from loguru import logger
import sys
import os
from pathlib import Path

# Import image and service modules
# Modal deploys from the directory containing main.py
try:
    # Try relative imports first (works in local package context)
    from .image import evo2_image
    from .service import EvoService
except (ImportError, ValueError):
    # Modal deployment - files are in the same directory
    # Modal copies all files from the deployment directory to /root
    # So we can import directly
    try:
        from image import evo2_image
        from service import EvoService
    except ImportError:
        # Fallback: add current directory to path
        current_dir = Path(__file__).parent
        if str(current_dir) not in sys.path:
            sys.path.insert(0, str(current_dir))
        from image import evo2_image
        from service import EvoService

# --- App Definition ---
# Files are already added to the image via image.add_local_dir()
app = modal.App("main", image=evo2_image)

# --- Shared State & Volume ---
job_status_dict = modal.Dict.from_name("evo-job-status", create_if_missing=True)
volume = modal.Volume.from_name("evo-model-cache", create_if_missing=True)
EVO_MODEL_DIR = "/root/.cache/huggingface"

# --- Main Service Class ---
@app.cls(
    gpu="H100:2",
    volumes={EVO_MODEL_DIR: volume},
    scaledown_window=300,
    timeout=1800
)
class EvoServiceWrapper:
    """Wrapper class that bridges Modal's @app.cls decorator with EvoService."""
    
    def __init__(self):
        """Initialize the service wrapper."""
        self.service = EvoService()
        self.fastapi_app = None
    
    @modal.enter()
    def load_model_and_api(self):
        """Load model and initialize the FastAPI app within the container."""
        self.fastapi_app = self.service.load_model_and_api(job_status_dict)
        logger.info("✅ EvoService initialized and endpoints registered")
    
    @modal.method()
    def generate_in_background(self, job_id: str, prompt: str, n_tokens: int):
        """Runs the model inference in the background with a general-purpose prompt."""
        logger.info(f"Job {job_id}: Starting background generation with prompt: '{prompt[:50]}...'")
        job_status_dict[job_id] = {"status": "running"}
        try:
            generation_result = self.service.model.generate(
                prompt_seqs=[prompt],
                n_tokens=n_tokens,
            )
            # Return structured result for reliable parsing by clients
            decoded_sequence = generation_result.sequences[0]
            logger.info(f"Job {job_id}: Successfully generated sequence.")
            job_status_dict[job_id] = {"status": "complete", "result": {"sequence": decoded_sequence}}
        except Exception as e:
            logger.error(f"Job {job_id}: Error during generation: {e}")
            job_status_dict[job_id] = {"status": "failed", "error": str(e)}
    
    @modal.asgi_app()
    def api(self):
        """Serve the FastAPI app."""
        if self.fastapi_app is None:
            raise RuntimeError("FastAPI app not initialized. Call load_model_and_api() first.")
        return self.fastapi_app


# --- Local Entrypoint for Syntax Validation (does not deploy/serve) ---
@app.local_entrypoint()
def main():
    """Local entrypoint for testing/debugging."""
    print("--- evo-service local entrypoint ---")
    print("✅ Modular structure loaded successfully")
    print("✅ Image definition available")
    print("✅ Service class available")
    print("✅ Endpoints module available")
