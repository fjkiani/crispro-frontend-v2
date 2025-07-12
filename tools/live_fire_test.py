import httpx
from loguru import logger
import time

# NOTE: You may need to update this if your app name changes in the modal deploy
COMMAND_CENTER_URL = "https://crispro--crispr-assistant-command-center-v2-commandcenter-api.modal.run"

def run_test():
    """
    Runs an end-to-end test of the Command Center and Evo Service.
    1. Creates a new Battle Plan for a synthetic patient.
    2. Requests guide RNA design for that Battle Plan.
    """
    patient_id = f"live_fire_{int(time.time())}"
    logger.info(f"STARTING LIVE FIRE TEST for patient: {patient_id}")
    logger.info(f"Targeting Command Center at: {COMMAND_CENTER_URL}")

    # Step 1: Formulate a battle plan to get a valid battle_plan_id
    # We need a valid ID to associate the generated guides with.
    logger.info("--- STEP 1: Formulating New Battle Plan ---")
    plan_payload = {
        "patient_identifier": patient_id, # Corrected field name
        "mutation_hgvs_c": "c.1799T>A", # Dummy data
        "mutation_hgvs_p": "p.Val600Glu", # Dummy data
        "gene": "BRAF", # Dummy data
        "sequence_for_perplexity": "GCTGAAACTTTGATGAAGAGACCTC",
        "protein_sequence": "LTEEGE", # Dummy data
        "transcript_id": "ENST00000288602", # Dummy data
        "bam_file_path": "data/vcf/dummy.bam"
    }
    
    battle_plan_id = None
    try:
        with httpx.Client(timeout=300.0) as client:
            response = client.post(f"{COMMAND_CENTER_URL}/v2/workflow/formulate_battle_plan", json=plan_payload, follow_redirects=True)
            response.raise_for_status()
            battle_plan = response.json()
            battle_plan_id = battle_plan.get("id")
            if not battle_plan_id:
                raise ValueError("Response did not contain a battle plan ID.")
            logger.success(f"Battle Plan CREATED. ID: {battle_plan_id}")
    except Exception as e:
        logger.error(f"CRITICAL FAILURE during battle plan creation: {e}")
        if hasattr(e, 'response'):
            logger.error(f"Response content: {e.response.text}")
        return

    # Step 2: Request guide RNA design for the created battle plan
    logger.info("--- STEP 2: Requesting Guide RNA Design from Zeta Forge ---")
    target_sequence_for_design = "GCTGAAACTTTGATGAAGAGACCTC"

    design_payload = {
        "battle_plan_id": battle_plan_id,
        "target_sequence": target_sequence_for_design
    }

    try:
        with httpx.Client(timeout=300.0) as client:
            response = client.post(f"{COMMAND_CENTER_URL}/v2/design/guide_rna", json=design_payload, follow_redirects=True)
            response.raise_for_status()
            design_result = response.json()
            logger.success("Guide RNA design request SUCCESSFUL!")
            logger.info("Response from Zeta Forge:")
            print(design_result)
            logger.info("--- TEST COMPLETE: SUCCESS ---")
    except Exception as e:
        logger.error(f"CRITICAL FAILURE during guide RNA design request: {e}")
        if hasattr(e, 'response'):
            logger.error(f"Response content: {e.response.text}")
        return

if __name__ == "__main__":
    run_test() 