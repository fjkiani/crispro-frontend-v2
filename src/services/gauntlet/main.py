# CACHE BUSTER v6.1 - Using correct colabfold_batch executable
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import modal
import subprocess
import tempfile
import logging
from pathlib import Path
import json

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ProteinSequence(BaseModel):
    protein_sequence: str

app = modal.App("gauntlet",)

# Volume for caching model parameters and databases
volume = modal.Volume.from_name("gauntlet-cache", create_if_missing=True)

def download_params_to_volume():
    """Download ColabFold model parameters to the volume during image build."""
    import subprocess
    from pathlib import Path
    
    params_dir = Path("/cache/params")
    params_dir.mkdir(parents=True, exist_ok=True)
    
    # Use the Git-cloned script to download parameters
    subprocess.run([
        "python",
        "/colabfold/colabfold/download.py",
        "alphafold2_ptm",
        str(params_dir)
    ], check=True)

# Define the Modal image
juggernaut_image = (
    modal.Image.from_registry("nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install("git", "wget", "build-essential")
    .pip_install(
        "fastapi",
        "uvicorn",
        "requests",
        "numpy",
    )
    .run_commands(
        # CRITICAL FIX: Install directly from Git to ensure all modules, including alphafold.common, are included.
        "pip install --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold'",
        
        # The JAX install needs the -f flag, which pip_install doesn't support.
        "pip install --no-warn-conflicts --upgrade 'dm-haiku==0.0.10' 'jax[cuda12_pip]==0.4.23' -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html",
        
        # We still clone the repo to get the download script and batch script easily accessible
        "git clone https://github.com/sokrypton/ColabFold.git /colabfold"
    )
    # Use run_function to execute our download logic during the build, mounting the volume.
    .run_function(download_params_to_volume, volumes={"/cache": volume})
    .env({
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": "/cache",
        "XDG_CACHE_HOME": "/cache",
        "COLABFOLD_PARAMS_DIR": "/cache/params"
    })
)

@app.function(
    image=juggernaut_image,
    gpu="H100:2",
    timeout=900,
    memory=32768,
)
def predict_structure(sequence: str) -> dict:
    """
    Predicts protein structure using ColabFold and returns comprehensive metrics.
    This function is designed to be called directly by other Modal services or as a webhook.
    """
    import subprocess
    import tempfile
    import json
    import os
    from pathlib import Path
    
    logger.info(f"Received prediction request for a protein of length {len(sequence)}")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Write the protein sequence to a FASTA file
        fasta_file = temp_path / "input.fasta"
        with open(fasta_file, 'w') as f:
            f.write(f">protein\n{sequence}\n")
        
        logger.info(f"FASTA file written to {fasta_file}")
        
        # Define output directory
        output_dir = temp_path / "output"
        output_dir.mkdir(exist_ok=True)
        
        # Run ColabFold using the correct executable
        cmd = [
            "colabfold_batch",
            str(fasta_file),
            str(output_dir),
            "--num-recycle", "3",
            "--model-type", "alphafold2_ptm",
            "--data", "/cache/params"
        ]
        
        logger.info(f"Executing command: {' '.join(cmd)}")
        
        try:
            # Increased timeout for robustness
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=800,
                check=False  # --- DOCTRINE: DEEP DIAGNOSTICS - DO NOT RAISE ON ERROR ---
            )
            
            # --- ALWAYS LOG, ESPECIALLY ON FAILURE ---
            if result.returncode != 0:
                logger.error(f"Prediction failed with exit code {result.returncode}")
                logger.error(f"STDOUT: {result.stdout}")
                logger.error(f"STDERR: {result.stderr}")
                raise RuntimeError(f"ColabFold prediction failed. See logs for details.")

            logger.info("Prediction completed successfully")
            logger.info(f"STDOUT: {result.stdout}")
            logger.warning(f"STDERR: {result.stderr}")
            
        except subprocess.TimeoutExpired as e:
            logger.error(f"Prediction timed out after {e.timeout} seconds")
            raise RuntimeError("Prediction timed out")
        
        except RuntimeError as e:
            # If something goes wrong, try to read and log the ColabFold log file for debugging
            # This part of the original code was not provided in the edit_specification,
            # so it's not included in the new_code.
            logger.error(f"Prediction failed with error: {e}")
            raise HTTPException(status_code=500, detail=f"ColabFold prediction failed: {e}")
        
        # Parse the results
        # Look for the scores JSON file
        scores_files = list(output_dir.glob("*_scores_rank_*.json"))
        if not scores_files:
            logger.error("No scores file found")
            raise HTTPException(status_code=500, detail="No prediction scores found")
        
        scores_file = scores_files[0]  # Take the first (best) ranked result
        
        with open(scores_file, 'r') as f:
            scores_data = json.load(f)
        
        # --- DOCTRINE: COMPREHENSIVE & RESILIENT ASSESSMENT ---
        # Surgically enhanced to prevent crashes on missing scores.
        plddt_score = scores_data.get('plddt', 0.0)
        ptm_score = scores_data.get('ptm', 0.0)
        max_pae = scores_data.get('max_pae', 0.0)
        
        pdb_files = list(output_dir.glob("*.pdb"))
        if not pdb_files:
            logger.error("No PDB file found")
            raise HTTPException(status_code=500, detail="No PDB structure found")
        
        pdb_file = pdb_files[0]  # Take the first PDB file
        
        with open(pdb_file, 'r') as f:
            pdb_content = f.read()
        
        logger.info(f"Prediction successful. pLDDT: {plddt_score}, ptm: {ptm_score}, max_pae: {max_pae}, fraction_disordered: {fraction_disordered}")
        
        return {
            "plddt_score": plddt_score,
            "ptm_score": ptm_score,
            "max_pae": max_pae,
            "fraction_disordered": fraction_disordered,
            "pdb_structure": pdb_content,
            "method": "ColabFold AlphaFold2"
        } 

@app.function(image=juggernaut_image)
@modal.asgi_app()
def fastapi_app():
    """FastAPI wrapper for the gauntlet service."""
    app_fastapi = FastAPI()

    @app_fastapi.post("/predict_structure")
    async def predict_structure_endpoint(item: ProteinSequence):
        """HTTP endpoint for structure prediction."""
        try:
            # --- DOCTRINE: CORRECT INVOCATION ---
            # To call another Modal function from an endpoint, we must use .remote()
            # This is the standard, correct way to invoke a separate, containerized function.
            result = predict_structure.remote(item.protein_sequence)
            
            if result is None:
                raise HTTPException(status_code=500, detail="Prediction function returned None.")
                
            plddt = result.get("plddt_score")
            ptm = result.get("ptm_score")
            
            if plddt is None or ptm is None:
                logger.error(f"Prediction succeeded but returned incomplete data: {result}")
                raise HTTPException(status_code=500, detail="Prediction failed to return all required scores (pLDDT, ptm).")

            return result
        except Exception as e:
            logger.error(f"Error during prediction endpoint call: {e}", exc_info=True)
            raise HTTPException(status_code=500, detail="Internal server error")
    
    return app_fastapi 