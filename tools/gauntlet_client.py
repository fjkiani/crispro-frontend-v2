import requests
import sys
import json
from loguru import logger

# --- Configuration ---
# This URL should point to your newly deployed Gauntlet service
GAUNTLET_URL = "https://crispro--gauntlet-service-juggernaut-fastapi-app.modal.run/predict_structure"

# A simple, well-behaved protein sequence for testing (a small part of human insulin)
# This is a known-good sequence that should fold without issues.
TEST_PROTEIN_SEQUENCE = "GIVEQCCTSICSLYQLENYCN"

def run_gauntlet_test():
    """
    Performs a direct, isolated test of the Gauntlet service.
    """
    logger.remove()
    logger.add(sys.stderr, level="INFO")

    logger.info(f"--- 🔬 ISOLATED TEST: STRIKING GAUNTLET AT {GAUNTLET_URL} 🔬 ---")
    
    payload = {"protein_sequence": TEST_PROTEIN_SEQUENCE}
    
    try:
        logger.info(f"Sending payload: {payload}")
        response = requests.post(GAUNTLET_URL, json=payload, timeout=900) # Long timeout for cold starts
        response.raise_for_status()
        
        result = response.json()
        
        logger.info("--- ✅ GAUNTLET RESPONDED SUCCESSFULLY ---")
        print(json.dumps(result, indent=2))
        
        # Final verification of the key metrics
        plddt = result.get("plddt_score")
        ptm = result.get("ptm_score")
        
        if plddt is not None and ptm is not None and plddt > 0:
            logger.success("Verification PASSED: Received valid pLDDT and PTM scores.")
        else:
            logger.error("Verification FAILED: Response did not contain valid scores.")

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"--- 💥 HTTP ERROR 💥 ---")
        logger.error(f"Status Code: {http_err.response.status_code}")
        logger.error(f"Response Body: {http_err.response.text}")
    except Exception as e:
        logger.error(f"--- 💥 CRITICAL FAILURE 💥 ---")
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)

if __name__ == "__main__":
    run_gauntlet_test() 
 