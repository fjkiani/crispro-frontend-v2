import modal
import os
import subprocess
import yaml
import tempfile
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from loguru import logger
import time
import json
import glob

# --- Image Definition ---
# The Boltz-2 model has complex dependencies. We define a dedicated image.
boltz_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "git+https://github.com/exscientia/boltz.git@9d2b45f47a6154b83030b4578b36873b221a71b4#egg=boltz[inference]",
        "torch==2.1.0",
        "pyyaml"
    )
    .apt_install("git", "clang")
)

# --- App Definition ---
app = modal.App("boltz-service-reforged", image=boltz_image)
fastapi_app = FastAPI(title="Boltz Service: Reforged Gauntlet")

# --- Pydantic Models ---

# Original model for binding affinity
class BoltzAffinityRequest(BaseModel):
    entities: list[dict] = Field(..., description="List of entities (proteins, ligands) for the simulation.")

class BoltzAffinityResponse(BaseModel):
    status: str
    affinity_score: float | None = None
    reason: str | None = None

# New models for structural integrity prediction (The Gauntlet)
class StructuralIntegrityRequest(BaseModel):
    protein_sequence: str = Field(..., description="The single protein sequence to be validated.")
    target_name: str = Field(..., description="A unique name for the prediction job, e.g., 'candidate_xyz_validation'.")

class StructuralIntegrityResponse(BaseModel):
    status: str
    plddt: float | None = None
    reason: str | None = None


# --- Main Service Logic ---

@app.function(
    gpu="A100", # Boltz requires significant GPU resources
    timeout=1800 # 30-minute timeout for the simulation
)
def run_boltz_affinity_job(req: BoltzAffinityRequest):
    """Original function to predict binding affinity between two entities."""
    logger.info("Starting Boltz-2 BINDING AFFINITY job.")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        config_path = os.path.join(temp_dir, "config.yml")
        fasta_path = os.path.join(temp_dir, "candidate.fasta")
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)

        # Write the FASTA file for the candidate
        with open(fasta_path, "w") as f:
            f.write(">candidate_inhibitor\n")
            f.write(req.entities[1]['sequence']) # Assuming second entity is the inhibitor

        # Write the YAML config
        config_data = {"entities": req.entities}
        with open(config_path, "w") as f:
            yaml.dump(config_data, f)
            
        command = [
            'boltz', 'predict',
            '--config_path', config_path,
            '--fasta_paths', fasta_path,
            '--output_dir', output_dir,
            '--random_seed', '0'
        ]
        
        logger.info(f"Running Boltz command: {' '.join(command)}")
        
        try:
            subprocess.run(command, capture_output=True, text=True, check=True)
            
            affinity_file = os.path.join(output_dir, "affinity.json")
            if os.path.exists(affinity_file):
                with open(affinity_file, "r") as f:
                    affinity_data = json.load(f)
                affinity_score = affinity_data.get("affinity", -999.0)
                logger.success(f"Simulation complete. Affinity score: {affinity_score}")
                return {"status": "complete", "affinity_score": affinity_score}
            else:
                logger.error("Simulation finished, but affinity.json not found.")
                return {"status": "error", "reason": "affinity.json not found"}

        except subprocess.CalledProcessError as e:
            logger.error(f"Boltz-2 simulation failed.")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            return {"status": "failed", "reason": e.stderr}

@app.function(
    gpu="A100",
    timeout=1800
)
def run_structure_prediction_job(req: StructuralIntegrityRequest):
    """
    New function for structural integrity, using the "self-docking" trick.
    """
    logger.info(f"Starting Boltz-2 STRUCTURAL INTEGRITY job for '{req.target_name}'.")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        config_path = os.path.join(temp_dir, "config.yml")
        fasta_path = os.path.join(temp_dir, "protein.fasta")
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)

        # Write the FASTA file for our single protein
        with open(fasta_path, "w") as f:
            f.write(f">{req.target_name}\n")
            f.write(req.protein_sequence)

        # Create the "self-docking" YAML config
        config_data = {
            "entities": [
                {"name": "self_target", "type": "protein", "sequence": req.protein_sequence},
                {"name": "self_candidate", "type": "protein", "sequence": req.protein_sequence}
            ]
        }
        with open(config_path, "w") as f:
            yaml.dump(config_data, f)
            
        command = [
            'boltz', 'predict',
            '--config_path', config_path,
            '--fasta_paths', fasta_path,
            '--output_dir', output_dir,
            '--random_seed', '0'
        ]
        
        logger.info(f"Running Boltz 'self-docking' command: {' '.join(command)}")
        
        try:
            # Execute the command
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            logger.debug(f"Boltz STDOUT: {result.stdout}")

            # --- Extract pLDDT from confidence file ---
            # Instead of affinity.json, we find the confidence JSON file.
            confidence_files = glob.glob(os.path.join(output_dir, "predictions", "*", "confidence_*.json"))
            if not confidence_files:
                logger.error("Simulation finished, but no confidence JSON file was found.")
                return {"status": "error", "reason": "confidence.json not found"}

            confidence_file_path = confidence_files[0]
            logger.info(f"Found confidence file: {confidence_file_path}")
            with open(confidence_file_path, "r") as f:
                confidence_data = json.load(f)
            
            # The key is typically 'complex_plddt' or just 'plddt'
            plddt_score = confidence_data.get("complex_plddt", confidence_data.get("plddt"))
            
            if plddt_score is None:
                logger.error("Confidence file found, but it does not contain a 'plddt' or 'complex_plddt' score.")
                return {"status": "error", "reason": "pLDDT score not in confidence file"}

            # The score is 0-1, we scale it to 0-100
            plddt_score_scaled = plddt_score * 100
            
            logger.success(f"Structural integrity job complete. pLDDT score: {plddt_score_scaled:.2f}")
            return {"status": "complete", "plddt": plddt_score_scaled}

        except subprocess.CalledProcessError as e:
            logger.error("Boltz structural integrity job failed.")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            return {"status": "failed", "reason": e.stderr}


# --- API Endpoints ---
@fastapi_app.post("/v1/predict_affinity", response_model=BoltzAffinityResponse)
def predict_affinity(request: BoltzAffinityRequest):
    """(Original Endpoint) Submits a job to predict binding affinity."""
    try:
        result = run_boltz_affinity_job.remote(request)
        return BoltzAffinityResponse(**result)
    except Exception as e:
        logger.error(f"Error submitting affinity job to Boltz: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@fastapi_app.post("/v1/predict_structure", response_model=StructuralIntegrityResponse)
def predict_structure(request: StructuralIntegrityRequest):
    """(New Endpoint) Submits a job to predict structural integrity (pLDDT)."""
    try:
        result = run_structure_prediction_job.remote(request)
        return StructuralIntegrityResponse(**result)
    except Exception as e:
        logger.error(f"Error submitting structure job to Boltz: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.function()
@modal.asgi_app()
def api():
    return fastapi_app 