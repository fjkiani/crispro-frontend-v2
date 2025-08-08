import requests
import os
import json
from loguru import logger
import time

# --- CONFIGURATION ---
BOLTZ_SERVICE_URL = "https://crispro--boltz-service-fastapi-app.modal.run"

def print_header(title):
    logger.info("=" * 80)
    logger.info(f"   {title.upper()}")
    logger.info("=" * 80)

def print_boltz_response(response):
    print_header("BOLTZ-2 DIRECT FIRE MISSION REPORT")
    if not response.get("ranked_candidates"):
        logger.warning("No candidates were returned from the Boltz service.")
        return

    for i, candidate in enumerate(response["ranked_candidates"]):
        logger.info(
            f"  Rank {i+1:02d}: {candidate['sequence']} "
            f"| Affinity: {candidate['binding_affinity']:.4f} "
            f"| Commentary: {candidate['commentary']}"
        )
    logger.info("-" * 80)

def run_direct_boltz_test():
    """
    Executes a direct-fire test on the boltz_service endpoint.
    """
    try:
        print_header("DIRECT-FIRE TEST INITIATED: OPERATION CANARY (FINAL)")
        
        # --- DOCTRINE: CANARY IN A COAL MINE ---
        # We have a perfectly running service that produces no output file.
        # This final test uses a simple, canonical ligand (Benzene) to determine
        # if the tool is capable of writing *any* confidence file at all.
        # This will prove whether the issue is our complex candidates or the tool itself.
        payload = {
            "target_sequence": "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVS",
            "candidate_sequences": [
                "c1ccccc1" # SMILES string for Benzene
            ],
            "job_id": f"direct_fire_{int(time.time())}"
        }

        logger.info(f"Dispatching FINAL CANARY command to Boltz service with payload: {payload}")
        
        response = requests.post(
            BOLTZ_SERVICE_URL, 
            json=payload, 
            timeout=1800 # Long timeout for the first run
        )
        response.raise_for_status()
        
        response_data = response.json()
        print_boltz_response(response_data)
        
        logger.success("--- ✅ MISSION COMPLETE: BOLTZ-2 RESPONDED TO DIRECT FIRE. ---")

    except requests.exceptions.HTTPError as e:
        logger.error(f"A fatal HTTP error occurred: {e.response.status_code} - {e.response.text}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    run_direct_boltz_test() 