import sys
import os
from loguru import logger
import requests
import time
import json

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

# --- CONFIGURATION ---
# CORRECTED: Pointing to the specialized Zeta Forge, not the generic Evo service.
ZETA_FORGE_URL = "https://crispro--zeta-forge-v1-api.modal.run"
TARGET_GENE = "VEGFA"
# These coordinates are from our successful reconnaissance run
TARGET_START_AMINO_ACID = 1
TARGET_END_AMINO_ACID = 395

def run_point_blank_test():
    """
    This test validates the full "Hunt, Slice, and Forge" workflow.
    1. Fetches the full DNA sequence for a gene.
    2. Slices out the vulnerable subdomain based on recon coordinates.
    3. Submits the subsequence to the Forge to generate a weapon.
    """
    logger.info("--- 🔥 POINT-BLANK FORGE STRIKE TEST 🔥 ---")

    # --- 1. Fetch Full DNA Sequence ---
    client = NCBIClient(threat_matrix={}) # Threat matrix not needed for this test
    logger.info(f"Fetching full DNA sequence for target gene: {TARGET_GENE}")
    full_sequence = client.get_gene_sequence(TARGET_GENE)

    if not full_sequence:
        logger.error("TEST FAILED: Could not retrieve DNA sequence. Aborting.")
        return

    # --- 2. Slice the Vulnerable Subdomain ---
    # Convert amino acid coordinates to nucleotide coordinates (1 AA = 3 bases)
    start_dna = (TARGET_START_AMINO_ACID - 1) * 3
    end_dna = TARGET_END_AMINO_ACID * 3
    target_subsequence = full_sequence[start_dna:end_dna]
    logger.info(f"Subdomain sliced. DNA coordinates: {start_dna}-{end_dna}. Length: {len(target_subsequence)} bp.")

    if len(target_subsequence) == 0:
        logger.error("TEST FAILED: Sliced subsequence is empty. Check coordinates.")
        return

    # --- 3. Submit to the Forge ---
    logger.info(f"Submitting {len(target_subsequence)}bp subsequence to ZetaForge at {ZETA_FORGE_URL}")
    # CORRECTED: Use the correct '/generate_inhibitor' endpoint for the Zeta Forge
    forge_endpoint = f"{ZETA_FORGE_URL}/generate_inhibitor"
    payload = {
        # CORRECTED: The Zeta Forge expects 'target_gene_sequence' and a 'design_goal'
        "target_gene_sequence": target_subsequence,
        "design_goal": f"Point-blank test strike against {TARGET_GENE} subdomain.",
        "optimization_loops": 3 # Use fewer loops for a quick test
    }

    try:
        response = requests.post(forge_endpoint, json=payload, timeout=60)
        response.raise_for_status()
        job_data = response.json()
        job_id = job_data.get("job_id")
        logger.success(f"✅ Successfully submitted job to ZetaForge. Job ID: {job_id}")
        logger.info("Polling for results... (this may take several minutes)")

        # --- 4. Poll for Results ---
        status_endpoint = f"{ZETA_FORGE_URL}/status/{job_id}"
        while True:
            time.sleep(30) # Poll every 30 seconds
            status_response = requests.get(status_endpoint)
            status_data = status_response.json()
            status = status_data.get("status")
            message = status_data.get("message", "Polling...") # 'message' may not exist, provide default
            logger.info(f"Forge status: {status} - {message}")

            if status == 'complete':
                logger.success("✅🔥 POINT-BLANK STRIKE SUCCESSFUL! 🔥✅")
                print("\n" + "="*50)
                print("          FORGED WEAPON BLUEPRINT")
                print("="*50)
                # The result is now at the top level of the response
                print(json.dumps(status_data, indent=2))
                print("="*50 + "\n")
                break
            elif status == 'failed':
                logger.error("❌ FORGE JOB FAILED ❌")
                logger.error(f"Failure reason: {status_data.get('error', 'Unknown error')}")
                break

    except Exception as e:
        logger.error(f"TEST FAILED during communication with ZetaForge: {e}")

if __name__ == "__main__":
    run_point_blank_test() 