import modal
import os
import torch
from fastapi import FastAPI

# --- Step 1: Add all dependencies from the main service ---
# This will test for installation or import-time conflicts.
debug_image = (
    modal.Image.debian_slim(python_version="3.11")
    # CRITICAL FIX: Add build-essential to provide compilers (g++, etc.)
    # that may be needed by pip to build underlying C/C++ dependencies.
    .apt_install("build-essential")
    .pip_install(
        "fastapi",
        "pydantic",
        "pandas",
        "pyarrow",
        "fair-esm", 
        "torch",
        "requests",
        "numpy"
    )
    # Step 2: Add the tools directory, just like the main service
    .add_local_dir("./tools", remote_path="/root/tools")
)

# --- Step 3: Add the client imports ---
# We are only importing them, not initializing them. This tests for import-time errors.
import sys
sys.path.append("/root")
from tools.alphamissense_client import AlphaMissenseClient
from tools.esm_client import ESMClient


# --- Volume Definition ---
alphamissense_volume = modal.Volume.from_name("alphamissense-data")

# --- App Definition ---
app = modal.App("fusion-engine-debug-v2", image=debug_image)

# --- Path to check ---
PARQUET_PATH = "/data/AlphaMissense_hg38.parquet"

@app.cls(volumes={"/data": alphamissense_volume}, gpu="any")
class SanityCheck:
    @modal.enter()
    def check_resources(self):
        """
        This method runs on container startup and checks if the resources are available.
        """
        print("--- SANITY CHECK (v2) ---")
        print("✅ Libraries installed and imported successfully.")
        
        # Check for GPU
        gpu_available = torch.cuda.is_available()
        print(f"✅ GPU Available: {gpu_available}")
            
        # Check for Data Volume file
        file_exists = os.path.exists(PARQUET_PATH)
        print(f"✅ Parquet File Exists at {PARQUET_PATH}: {file_exists}")

        print("--- SANITY CHECK (v2) COMPLETE ---")

    @modal.asgi_app()
    def api(self):
        """Defines the FastAPI web service."""
        fastapi_app = FastAPI(title="Sanity Check Service (v2)")

        @fastapi_app.get("/health")
        def health_check():
            """A simple health check endpoint."""
            return {"status": "ok", "message": "Sanity check v2 container is running."}
        
        return fastapi_app

@app.local_entrypoint()
def main():
    print("Deploying the sanity check app (v2)...")
    print("Run: modal deploy services/fusion_engine/debug_main.py") 