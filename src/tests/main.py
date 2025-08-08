import modal
import os
import sys
from io import BytesIO
from pathlib import Path # Path is a standard library and safe to import globally.

# --- SETUP: PATHS & SCHEMAS ---
# This sys.path.append is for LOCAL execution and testing only.
# The container's PYTHONPATH is set in the image definition.
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# --- MODAL APP CONFIGURATION ---
app = modal.App(name="boltz-service")

# --- IMAGE DEFINITIONS ---
# This image installs boltz, which brings in PyYAML as a dependency.
boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .run_commands("pip install uv")
    .run_commands("uv pip install --system --compile-bytecode boltz==2.1.1 pydantic rdkit")
    .env({"PYTHONPATH": "/root"})
    .add_local_dir("src", "/root/src")
)

web_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("fastapi[standard]")
    .env({"PYTHONPATH": "/root"})
    .add_local_dir("src", "/root/src")
)

# --- VOLUME DEFINITION: FOR PERSISTENT MODEL WEIGHTS ---
boltz_model_volume = modal.Volume.from_name("boltz-models", create_if_missing=True)
models_dir = Path("/models/boltz")

# --- MODEL DOWNLOADER: A DEDICATED, ONE-OFF FUNCTION ---
downloader_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("huggingface_hub[hf_transfer]==0.26.3")
    .env({"HF_HUB_ENABLE_HF_TRANSFER": "1"})
)

@app.function(
    image=downloader_image,
    volumes={
        str(models_dir): boltz_model_volume,
    },
    timeout=1200,
)
def download_model(force_download: bool = False):
    """
    Downloads the Boltz-2 model weights from Hugging Face into the modal.Volume.
    This only needs to be run once, or when the model is updated.
    """
    from huggingface_hub import snapshot_download
    print("Downloading Boltz-2 model weights...")
    snapshot_download(
        repo_id="boltz-community/boltz-2",
        local_dir=models_dir,
        force_download=force_download,
    )
    boltz_model_volume.commit()
    print("✅ Model download complete.")

# --- INFERENCE FUNCTION: THE CORE WEAPON ---
@app.function(
    image=boltz_image,
    volumes={str(models_dir): boltz_model_volume},
    gpu="H100", # Boltz-2 requires a powerful GPU
    timeout=1800, # Allow ample time for batch processing
)
def run_boltz_inference(target_sequence: str, candidate_sequences: list, job_id: str) -> list:
    """
    Runs Boltz-2 inference on a target and a list of candidate sequences.
    This function now implements the "False Flag V2" doctrine, derived
    from direct source code analysis of the Boltz-2 tool.
    """
    # --- DOCTRINE: LAZY IMPORTS ---
    import yaml
    import subprocess
    import json
    from src.services.boltz_service.schemas import RankedCandidate
    from rdkit import Chem

    results = []
    
    # We must process each candidate as a separate job because the `affinity`
    # property only supports one binder at a time.
    for i, candidate_seq in enumerate(candidate_sequences):
        complex_name = f'complex_{i}'
        candidate_id = f'CANDIDATE_{i}'

        # --- "FALSE FLAG V2" PAYLOAD CONSTRUCTION ---
        # This is the one true schema dictated by `boltz/data/parse/schema.py`.
        # We declare the target as a protein and falsely classify our candidate
        # nanobody as a 'ligand' to trigger the interaction modeling pipeline.
        # RDKit is used to convert the peptide sequence to a SMILES string.
        
        mol = Chem.MolFromSequence(candidate_seq)
        smiles = Chem.MolToSmiles(mol) if mol else ""

        if not smiles:
            print(f"Could not convert candidate sequence {candidate_seq} to SMILES. Skipping.")
            continue

        final_input_data = {
            'version': 1,
            'sequences': [
                {'protein': {'id': 'TARGET', 'sequence': target_sequence}},
                {'ligand': {'id': candidate_id, 'smiles': smiles}}
            ],
            'properties': [
                {'affinity': {'binder': candidate_id}}
            ]
        }

        input_path = Path(f"/tmp/input_{job_id}_{complex_name}.yaml")
        with open(input_path, 'w') as f:
            yaml.dump(final_input_data, f)

        # Run the Boltz-2 command line tool for this specific complex
        output_dir = f"boltz_results_{job_id}"
        cmd = [
            "boltz", "predict", str(input_path),
            "--use_msa_server",
            "--cache", str(models_dir),
            "--out_dir", output_dir
        ]
        print(f"Running command for {complex_name}: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            print(f"--- BOLTZ STDOUT for {complex_name} ---")
            print(result.stdout)
            print(f"--- BOLTZ STDERR for {complex_name} ---")
            print(result.stderr)
            result.check_returncode() # Manually raise exception if failed
        except subprocess.CalledProcessError as e:
            print(f"Error processing {complex_name}. Skipping.")
            continue

        # --- PARSE THE TRUE SIGNAL: IPTM ---
        # The affinity score is irrelevant. We hunt for the iptm score in the
        # confidence file, which proves the interaction was modeled.
        confidence_file = Path(output_dir) / "predictions" / complex_name / f"confidence_{complex_name}_model_0.json"
        
        print(f"Searching for confidence file at: {confidence_file}")
        
        binding_affinity = -999.0 # Default to a low score
        if confidence_file.exists():
            print(f"Found confidence file for {complex_name}. Parsing for iptm score...")
            with open(confidence_file, 'r') as f:
                confidence_data = json.load(f)
                binding_affinity = confidence_data.get('iptm', -999.0)
        else:
            print(f"Confidence file NOT FOUND for {complex_name}.")

        commentary = "High Confidence Interface" if binding_affinity > 0.7 else "Moderate Confidence Interface"
        results.append(
            RankedCandidate(
                sequence=candidate_seq,
                binding_affinity=binding_affinity,
                commentary=commentary
            )
        )
    
    results.sort(key=lambda x: x.binding_affinity, reverse=True)
    return results

# --- WEB ENDPOINT: THE PUBLIC-FACING API ---
@app.function(image=web_image, timeout=1810)
@modal.asgi_app()
def fastapi_app():
    from fastapi import FastAPI
    from src.services.boltz_service.schemas import InteractionRequest, InteractionResponse
    
    web_app = FastAPI(title="Boltz-2 Interaction Prediction Service")

    @web_app.post("/", response_model=InteractionResponse)
    async def predict(request: InteractionRequest):
        print(f"Received screening request for job {request.job_id}")
        ranked_candidates = run_boltz_inference.remote(
            request.target_sequence,
            request.candidate_sequences,
            request.job_id
        )
        return InteractionResponse(job_id=request.job_id, ranked_candidates=ranked_candidates)
    
    return web_app

# --- CLI ENTRYPOINT FOR DOWNLOADING THE MODEL ---
@app.local_entrypoint()
def download(force: bool = False):
    """A local CLI entrypoint to trigger the model download."""
    print("Triggering remote model download...")
    download_model.remote(force_download=force)
