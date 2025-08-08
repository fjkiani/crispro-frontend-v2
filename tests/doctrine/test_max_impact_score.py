import sys
import os
import requests
from loguru import logger
import json

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

# --- CONFIGURATION ---
# As per Zeta Doctrine, we use direct HTTP strikes.
ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run"
TARGET_GENE = "VEGFA"

def calculate_max_impact_score():
    """
    This test calculates the theoretical "Maximum Impact Score" for a gene.
    It does this by comparing the gene's full sequence against a total knockout (empty sequence).
    The resulting Zeta Score is our benchmark for a perfect weapon.
    """
    logger.info("--- 💥 OPERATION: MAXIMUM IMPACT 💥 ---")
    logger.info(f"Calculating theoretical kill-shot score for target: {TARGET_GENE}")

    # 1. Acquire Target Sequence
    # We instantiate the NCBIClient with an empty threat matrix because we don't need
    # domain analysis for this operation, only sequence fetching.
    client = NCBIClient(threat_matrix={})
    logger.info(f"Fetching full DNA sequence for {TARGET_GENE}...")
    reference_sequence = client.get_gene_sequence(TARGET_GENE)

    if not reference_sequence:
        logger.error("TEST FAILED: Could not retrieve reference DNA sequence. Aborting.")
        return

    logger.success(f"Acquired reference sequence (Length: {len(reference_sequence)} bp).")

    # 2. Define Annihilation & Strike Oracle
    # The ultimate perturbation is a frameshift mutation at the start of the sequence.
    # This scrambles every subsequent codon, resulting in a non-functional protein.
    # We use this instead of an empty string, which the Oracle rejects.
    perturbed_sequence = "G" + reference_sequence
    
    logger.info(f"Striking Zeta Oracle at {ZETA_ORACLE_URL} to calculate impact score...")
    
    oracle_endpoint = f"{ZETA_ORACLE_URL}/invoke"
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": reference_sequence,
            "alternate_sequence": perturbed_sequence
        }
    }

    try:
        with requests.post(oracle_endpoint, json=payload, timeout=900, stream=True) as response:
            response.raise_for_status()
            
            # Since the response might be large or streamed, iterate over lines
            for line in response.iter_lines():
                if line:
                    try:
                        data = json.loads(line)
                        zeta_score = data.get("zeta_score")

                        if zeta_score is not None:
                            logger.success("✅--- MAXIMUM IMPACT SCORE CALCULATED ---✅")
                            print("\n" + "="*60)
                            print(f"  Theoretical Max Impact Zeta Score for {TARGET_GENE}: {zeta_score}")
                            print("="*60 + "\n")
                            logger.info("This score is now our benchmark for weapon effectiveness.")
                            return
                    except json.JSONDecodeError:
                        logger.warning(f"Received non-JSON line from Oracle: {line}")
            
            logger.error("TEST FAILED: Oracle response did not contain a valid Zeta Score.")

    except requests.exceptions.RequestException as e:
        logger.error(f"TEST FAILED: Error communicating with Zeta Oracle: {e}")

if __name__ == "__main__":
    calculate_max_impact_score() 