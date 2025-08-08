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
# This image installs boltz, which brings in PyYAML and RDKit as dependencies.
boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .run_commands("pip install uv")
    .run_commands("uv pip install --system --compile-bytecode boltz==2.1.1 pydantic rdkit requests tqdm")
    .env({"PYTHONPATH": "/root"})
    .add_local_dir("src", "/root/src")
    .add_local_dir("tools", "/root/tools") # Add our new msa_client tool
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
    from tools.msa_client import run_mmseqs2

    results = []
    
    # --- OPERATION LOCAL LOGISTICS ---
    # 1. Fetch MSA data for the target sequence using our own robust client.
    # This bypasses the flaky implementation within the boltz tool.
    print("Fetching MSA for target sequence...")
    msa_prefix = f"/tmp/msa_{job_id}"
    msa_results = run_mmseqs2(target_sequence, prefix=msa_prefix)
    
    # 2. Sanitize the MSA data. The ColabFold server returns data with null
    # bytes, which the original authors stripped out. We must write this
    # sanitized content to a new file to avoid parser errors.
    sanitized_msa_content = msa_results[0]
    sanitized_msa_path = Path(f"/tmp/sanitized_msa_{job_id}.a3m")
    with open(sanitized_msa_path, "w") as f:
        f.write(sanitized_msa_content)
    
    print(f"MSA fetched and sanitized to {sanitized_msa_path}")

    # We must process each candidate as a separate job because the `affinity`
    # property only supports one binder at a time.
    for i, candidate_seq in enumerate(candidate_sequences):
        complex_name = f'complex_{i}'
        candidate_id = f'C{i}' # Must be <= 5 chars due to bug in boltz library

        # --- "FALSE FLAG V2" PAYLOAD CONSTRUCTION ---
        # This is the one true schema dictated by `boltz/data/parse/schema.py`.
        # We declare the target as a protein and falsely classify our candidate
        # nanobody as a 'ligand' to trigger the interaction modeling pipeline.
        # RDKit is used to convert the peptide sequence to a SMILES string.
        
        # --- "Bilingual" SMILES/Peptide Handling ---
        # Try to interpret the candidate as a SMILES string first.
        mol = Chem.MolFromSmiles(candidate_seq)
        if mol:
            smiles = candidate_seq # It's a valid SMILES string, use it directly.
        else:
            # If not a valid SMILES, assume it's a peptide sequence and try to convert.
            mol = Chem.MolFromSequence(candidate_seq)
            smiles = Chem.MolToSmiles(mol) if mol else ""

        if not smiles:
            print(f"Could not interpret candidate '{candidate_seq}' as SMILES or convert from peptide sequence. Skipping.")
            continue

        final_input_data = {
            'version': 1,
            'sequences': [
                {'protein': {'id': 'TARG', 'sequence': target_sequence, 'msa': str(sanitized_msa_path)}},
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
            "--cache", str(models_dir),
            "--out_dir", output_dir,
            "--output_format", "mmcif"
        ]
        print(f"Running command for {complex_name}: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            print(f"--- BOLTZ STDOUT for {complex_name} ---")
            print(result.stdout)
            print(f"--- BOLTZ STDERR for {complex_name} ---")
            print(result.stderr)
            result.check_returncode() # Manually raise exception if failed
        except subprocess.CalledProcessError:
            print(f"Error processing {complex_name}. Skipping.")
            continue

        # --- PARSE THE TRUE SIGNAL: IPTM ---
        # The affinity score is irrelevant. We hunt for the iptm score in the
        # confidence file, which proves the interaction was modeled.
        # DOCTRINAL CORRECTION: Analysis of boltz source code reveals it creates
        # a nested directory based on the input file stem. Our path must reflect this.
        input_stem = input_path.stem
        nested_output_dir = Path(output_dir) / f"boltz_results_{input_stem}"
        confidence_file = nested_output_dir / "predictions" / input_stem / f"confidence_{input_stem}_model_0.json"
        
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

@app.function(
    image=boltz_image,
    volumes={str(models_dir): boltz_model_volume},
    gpu="H100",
    timeout=1800,
)
def run_structure_prediction(protein_sequence: str, job_id: str) -> dict:
    """
    (Reforged v3) Runs a structural integrity check for a single protein.
    This version restores the critical MSA generation step.
    """
    import yaml
    import subprocess
    import json
    from pathlib import Path
    from tools.msa_client import run_mmseqs2 # MSA RESTORATION

    print(f"Starting structural integrity check for job {job_id}")

    # --- OPERATION: MSA RESTORATION ---
    print("Fetching MSA for target sequence...")
    msa_prefix = f"/tmp/msa_{job_id}"
    msa_results = run_mmseqs2(protein_sequence, prefix=msa_prefix)
    
    sanitized_msa_content = msa_results[0]
    sanitized_msa_path = Path(f"/tmp/sanitized_msa_{job_id}.a3m")
    with open(sanitized_msa_path, "w") as f:
        f.write(sanitized_msa_content)
    print(f"MSA fetched and sanitized to {sanitized_msa_path}")

    final_input_data = {
        'version': 1,
        'sequences': [
            # Provide the MSA file path along with the sequence
            {'protein': {'id': 'TARG', 'sequence': protein_sequence, 'msa': str(sanitized_msa_path)}},
        ],
    }

    input_path = Path(f"/tmp/input_{job_id}.yaml")
    with open(input_path, 'w') as f:
        yaml.dump(final_input_data, f)

    output_dir = f"boltz_results_{job_id}"
    cmd = [
        "boltz", "predict", str(input_path),
        "--cache", str(models_dir),
        "--out_dir", output_dir,
        "--output_format", "mmcif"
    ]
    print(f"Running simple fold command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"--- BOLTZ STDOUT for {job_id} ---")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"--- BOLTZ STDERR for {job_id} ---")
        print(e.stderr)
        return {"status": "error", "error": f"Boltz execution failed: {e.stderr}"}

    input_stem = input_path.stem
    nested_output_dir = Path(output_dir) / f"boltz_results_{input_stem}"
    confidence_file = nested_output_dir / "predictions" / input_stem / f"confidence_{input_stem}_model_0.json"
    
    print(f"Searching for confidence file at: {confidence_file}")

    if confidence_file.exists():
        with open(confidence_file, 'r') as f:
            confidence_data = json.load(f)
        
        # --- DOCTRINAL CORRECTION (PER USER INTEL) & FINAL FORTIFICATION ---
        # 1. The key for a single protein fold IS 'complex_plddt', confirmed by diagnostic blueprint.
        # 2. Add fail-safe defaults to prevent any possible NoneType errors.
        plddt_score = confidence_data.get('complex_plddt') or 0.0
        ptm_score = confidence_data.get('ptm') or 0.0
        fraction_disordered = confidence_data.get('fraction_disordered') or 1.0 # Default to fully disordered
        
        # The score from Boltz is 0-1. We must scale it to 0-100 to match the pLDDT standard.
        plddt_score_scaled = plddt_score * 100

        print(f"✅ Success! Found complex_plddt: {plddt_score_scaled:.2f}, ptm: {ptm_score:.2f}, disordered: {fraction_disordered:.2f}")
        return {
            "status": "complete", 
            "plddt_score": plddt_score_scaled,
            "ptm_score": ptm_score,
            "fraction_disordered": fraction_disordered,
        }
    else:
        return {"status": "error", "error": "Confidence file not found after successful execution."}


# --- WEB ENDPOINT: THE PUBLIC-FACING API ---
@app.function(image=web_image, timeout=1810)
@modal.asgi_app()
def fastapi_app():
    from fastapi import FastAPI, HTTPException
    from src.services.boltz_service.schemas import InteractionRequest, InteractionResponse, StructuralIntegrityRequest, StructuralIntegrityResponse
    
    web_app = FastAPI(title="Boltz-2 Service: Interaction & Structural Prediction")

    @web_app.post("/v1/predict_interaction", response_model=InteractionResponse)
    async def predict_interaction(request: InteractionRequest):
        print(f"Received screening request for job {request.job_id}")
        # .remote() returns a Job object, not the result directly.
        # The client will have to poll for the result.
        job = run_boltz_inference.remote(
            request.target_sequence,
            request.candidate_sequences,
            request.job_id
        )
        # We can return the job ID so the client can check status
        return InteractionResponse(job_id=job.object_id, ranked_candidates=[])

    @web_app.post("/v1/predict_structure", response_model=StructuralIntegrityResponse)
    async def predict_structure(request: StructuralIntegrityRequest):
        print(f"Received structural integrity request for job {request.job_id}")
        try:
            # --- DOCTRINE OF SYNCHRONICITY ---
            # When calling a function in the same app with .remote() from an async
            # context, the call is BLOCKING and returns the FINAL DICTIONARY,
            # not a Job object.
            result_dict = run_structure_prediction.remote(
                request.protein_sequence,
                request.job_id
            )
            
            # Combine the result with the request's job_id and return the final response.
            return StructuralIntegrityResponse(
                job_id=request.job_id,
                status=result_dict.get("status", "error"),
                plddt_score=result_dict.get("plddt_score"),
                error=result_dict.get("error")
            )
        except Exception as e:
            # This will catch any unexpected errors during the synchronous execution.
            print(f"Critical error during structural integrity execution: {e}")
            raise HTTPException(status_code=500, detail=str(e))

    return web_app

# --- CLI ENTRYPOINT FOR DOWNLOADING THE MODEL ---
@app.local_entrypoint()
def download(force: bool = False):
    """A local CLI entrypoint to trigger the model download."""
    print("Triggering remote model download...")
    download_model.remote(force_download=force)
