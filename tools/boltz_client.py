import requests
import sys
import json
import uuid
from loguru import logger

# --- Configuration ---
# This URL should point to your newly deployed Boltz service
BOLTZ_URL = "https://crispro--boltz-service-fastapi-app.modal.run/v1/predict_structure"

# Using our robust, 35-amino-acid Villin Headpiece Subdomain for the test
TEST_PROTEIN_SEQUENCE = "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF"

def run_boltz_test():
    """
    Performs a direct, isolated test of the Boltz service.
    """
    logger.remove()
    logger.add(sys.stderr, level="INFO")

    job_id = f"boltz_isolated_test_{uuid.uuid4().hex[:8]}"
    logger.info(f"--- 🔬 ISOLATED TEST: STRIKING BOLTZ AT {BOLTZ_URL} 🔬 ---")
    logger.info(f"Using Job ID: {job_id}")
    
    payload = {
        "protein_sequence": TEST_PROTEIN_SEQUENCE,
        "job_id": job_id
    }
    
    try:
        logger.info(f"Sending payload: {json.dumps(payload)}")
        response = requests.post(BOLTZ_URL, json=payload, timeout=1800) # Long timeout for cold starts and GPU work
        response.raise_for_status()
        
        result = response.json()
        
        logger.info("--- ✅ BOLTZ RESPONDED SUCCESSFULLY ---")
        print(json.dumps(result, indent=2))
        
        # Final verification of the key metrics
        status = result.get("status")
        plddt = result.get("plddt_score")
        
        if status == "complete" and plddt is not None and plddt > 0:
            logger.success(f"Verification PASSED: Received status '{status}' and valid pLDDT score: {plddt:.2f}")
        else:
            logger.error(f"Verification FAILED: Status was '{status}' or scores were invalid.")
            logger.error(f"Full response: {result}")

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"--- 💥 HTTP ERROR 💥 ---")
        logger.error(f"Status Code: {http_err.response.status_code}")
        logger.error(f"Response Body: {http_err.response.text}")
    except Exception as e:
        logger.error(f"--- 💥 CRITICAL FAILURE 💥 ---")
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)

if __name__ == "__main__":
    run_boltz_test() 