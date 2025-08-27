# Mission: Fusion Engine Service
# This service is a lightweight, decoupled engine for orchestrating variant scoring
# using AlphaMissense and ESM, deliberately excluding the heavy Evo2 dependencies.
import modal
import os
import re
from fastapi import FastAPI
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Dict, Any
import asyncio
import logging

# It's critical to be able to import from the tools directory.
import sys
# CORRECTED PATHING LOGIC: Add the project's root `src` directory to the path
# to ensure that the `tools` module can be found during the Modal build process.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from tools.alphamissense_client import AlphaMissenseClient
from tools.esm_client import ESMClient

# --- Lightweight Image Definition ---
# This image contains all necessary dependencies for both clients.
fusion_image = (
    modal.Image.debian_slim(python_version="3.11")
    # CRITICAL FIX: Add build-essential to provide compilers (g++, etc.)
    # that are needed by some of torch's dependencies (e.g., triton).
    .apt_install("build-essential")
    .pip_install(
        "fastapi",
        "pydantic",
        "pandas",
        "pyarrow",
        "fair-esm", 
        "torch",
        "requests",
        "numpy" # Explicitly add numpy to resolve dependency issue.
    )
    # Copy the 'tools' directory into the container image at build time.
    .add_local_dir("src/tools", remote_path="/root/tools")
)

# Define the volume where our large data files are stored.
alphamissense_volume = modal.Volume.from_name("alphamissense-data")

# --- App Definition ---
app = modal.App("fusion-engine", image=fusion_image)

# --- Pydantic Models for API Validation (Corrected Structure) ---
class VariantInput(BaseModel):
    variant_id: str
    hgvs: str
    alphamissense_variant_str: str 

class ScoringRequest(BaseModel):
    protein_sequence: str
    variants: List[VariantInput]

class ScoredVariant(BaseModel):
    variant_id: str
    alphamissense_score: float
    esm_score: float
    zeta_score: float # The fused score

class ScoringResponse(BaseModel):
    scored_variants: List[ScoredVariant]

# With the correct data, a generic GPU is sufficient.
@app.cls(cpu=8, memory=65536, timeout=600, volumes={"/data": alphamissense_volume}, gpu="any")
class FusionEngine:
    @modal.enter()
    def load_clients(self):
        """Initializes the clients when the container starts."""
        print("🚀 Initializing Fusion Engine clients...")
        # CRITICAL FIX: Point the client to the correct path in the mounted volume.
        self.am_client = AlphaMissenseClient(data_path="/data/AlphaMissense_hg38.parquet")
        self.esm_client = ESMClient()
        print("✅ Clients initialized. Fusion Engine is OPERATIONAL.")

    @modal.asgi_app()
    def api(self):
        """Defines the FastAPI web service."""
        fastapi_app = FastAPI(
            title="Fusion Engine",
            description="A lightweight service to score variants using AlphaMissense and ESM.",
            version="1.0.0"
        )

        @fastapi_app.post("/score_variants", response_model=ScoringResponse)
        async def score_variants(request: ScoringRequest):
            """
            This is the main fusion engine endpoint.
            It orchestrates calls to AlphaMissense and ESM in parallel.
            """
            try:
                # --- Step 1: Dispatch All Tasks in Parallel ---
                
                # Only run ESM task if the model loaded correctly
                if self.esm_client and self.esm_client.model:
                    esm_task = asyncio.to_thread(
                        self.esm_client.get_scores,
                        request.protein_sequence,
                        [v.hgvs for v in request.variants]
                    )
                else:
                    logging.warning("ESM model not loaded. Skipping ESM scoring.")
                    # Create a dummy task that returns an empty dict
                    async def dummy_esm_task():
                        return {}
                    esm_task = asyncio.create_task(dummy_esm_task())

                # CORRECT LOGIC: AlphaMissense client takes one variant at a time.
                # We create a separate, concurrent task for each one.
                am_tasks = [
                    asyncio.to_thread(self.am_client.get_score, v.alphamissense_variant_str)
                    for v in request.variants
                ]
                
                # Gather results from ESM and all individual AlphaMissense calls
                all_tasks = [esm_task] + am_tasks
                results = await asyncio.gather(*all_tasks, return_exceptions=True)
                
                # --- Step 2: Unpack and Fuse Results ---
                esm_results = results[0] if not isinstance(results[0], Exception) else {}
                am_results_list = results[1:]
                
                final_scores = []
                for i, v in enumerate(request.variants):
                    # Unpack AlphaMissense result for this variant
                    am_result = am_results_list[i] if i < len(am_results_list) and isinstance(am_results_list[i], dict) else {}
                    am_score = am_result.get("score", -999.0)

                    # Unpack ESM result for this variant
                    esm_result = esm_results.get(v.hgvs, {}) if isinstance(esm_results, dict) else {}
                    esm_score = esm_result.get("esm_score", -999.0)
                    
                    # Fuse the scores. For now, a simple average, ignoring errors.
                    valid_scores = [s for s in [am_score, esm_score] if s not in [-999.0, -998.0, -997.0]]
                    zeta_score = sum(valid_scores) / len(valid_scores) if valid_scores else -999.0
                    
                    final_scores.append(ScoredVariant(
                        variant_id=v.variant_id,
                        alphamissense_score=round(am_score, 4),
                        esm_score=round(esm_score, 4),
                        zeta_score=round(zeta_score, 4)
                    ))
                
                return ScoringResponse(scored_variants=final_scores)
            except Exception as e:
                logging.error(f"💥 FUSION ENGINE CRITICAL ERROR: {e}", exc_info=True)
                return JSONResponse(
                    content={"status": "error", "message": f"Internal server error: {str(e)}"},
                    status_code=500
                )
        
        @fastapi_app.get("/health")
        def health_check():
            """A simple health check endpoint."""
            return {"status": "healthy", "service": "fusion-engine"}
        
        return fastapi_app

@app.local_entrypoint()
def main():
    """A local test function to verify the fusion engine is working."""
    print("--- ⚔️ LOCAL FUSION ENGINE TEST ⚔️ ---")
    print("   This confirms the service code is syntactically valid.")
    print("   To perform a full end-to-end test, deploy the service with:")
    print("   modal deploy services/fusion_engine/main.py") 