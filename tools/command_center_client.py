import httpx
import time
import os
from loguru import logger

# --- Configuration ---
COMMAND_CENTER_URL = "https://crispro--command-center-v8-override-fix-web-app.modal.run"
EXECUTE_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/execute"
STATUS_ENDPOINT = f"{COMMAND_CENTER_URL}/status"

def load_brca1_bait_sequence() -> str:
    """Loads the BRCA1 bait sequence from the local FASTA file."""
    try:
        # Construct the path relative to this script's location
        script_dir = os.path.dirname(__file__)
        fasta_path = os.path.join(script_dir, '..', 'data', 'gene_database', 'BRCA1.fasta')
        
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            # Skip the header line and join the rest
            return "".join(line.strip() for line in lines[1:])
    except FileNotFoundError:
        logger.error(f"FATAL: BRCA1.fasta not found at {fasta_path}. Run tools/brca1_recon.py to generate it.")
        return None
    except Exception as e:
        logger.error(f"An error occurred loading the BRCA1 sequence: {e}")
        return None

def run_client():
    """
    Main function to run the command center client.
    --- DOCTRINE: MANUAL TARGETING OVERRIDE ---
    This client now bypasses the Hunter-Analyst by providing a high-quality
    bait sequence directly, aligning with the "BRCA1 Gauntlet" and "Precision Forge" protocols.
    """
    logger.info("--- 🚀 Initiating Structural Integrity Protocol (Live Fire Exercise) 🚀 ---")

    brca1_sequence = load_brca1_bait_sequence()
    if not brca1_sequence:
        return

    # Using the first 2000bp as a high-context prompt. The Oracle's limit is ~4k.
    # This ensures we provide rich context without exceeding payload limits.
    payload = {
        "bait_sequence": brca1_sequence[:2000]
    }

    try:
        # --- 1. Execute Workflow ---
        logger.info(f"Sending request to {EXECUTE_ENDPOINT} with manual BRCA1 bait sequence...")
        with httpx.Client() as client:
            response = client.post(EXECUTE_ENDPOINT, json=payload, timeout=60.0)
            response.raise_for_status()
            data = response.json()
            workflow_id = data.get("workflow_id")
            if not workflow_id:
                logger.error("❌ Failed to initiate workflow. No workflow_id received.")
                return
            logger.success(f"✅ Workflow initiated successfully. Workflow ID: {workflow_id}")

        # --- 2. Poll for Status ---
        status_url = f"{STATUS_ENDPOINT}/{workflow_id}"
        logger.info(f"Polling for status at {status_url}...")
        final_status = None
        timeout = 600  # 10 minutes timeout for the entire workflow
        start_time = time.time()

        while time.time() - start_time < timeout:
            try:
                with httpx.Client() as client:
                    response = client.get(status_url, timeout=30.0)
                if response.status_code == 200:
                    status_data = response.json()
                    current_status = status_data.get("status")
                    message = status_data.get("message")
                    logger.info(f"[Status: {current_status}] {message}")

                    if current_status in ["complete", "failed"]:
                        final_status = status_data
                        break
                else:
                    logger.warning(f"Received status {response.status_code} while polling.")
                
                time.sleep(5) # Wait 5 seconds between polls

            except httpx.RequestError as e:
                logger.error(f"An error occurred while polling: {e}")
                time.sleep(10) # Wait longer on network errors

        # --- 3. Display Final Report ---
        if final_status:
            logger.info("\n--- 🏁 FINAL REPORT 🏁 ---")
            logger.info(f"Workflow ID: {final_status.get('job_id')}")
            logger.info(f"Final Status: {final_status.get('status')}")
            logger.info(f"Final Message: {final_status.get('message')}")
            
            events = final_status.get("events", [])
            logger.info("\n--- 📜 Event Log 📜 ---")
            for event in events:
                logger.info(f"  [{event.get('timestamp')}] {event.get('type')}: {event.get('message')}")

            results = final_status.get("result", [])
            if results:
                logger.info("\n--- 🏆 Validated Weapons 🏆 ---")
                for i, weapon in enumerate(results):
                    logger.success(f"  Weapon {i+1}:")
                    logger.success(f"    pLDDT: {weapon.get('plddt_score'):.2f}")
                    logger.success(f"    Delta Score: {weapon.get('delta_score'):.4f}")
                    logger.success(f"    Protein (first 30aa): {weapon.get('protein_sequence', '')[:30]}...")
                    logger.success(f"    DNA (first 30bp): {weapon.get('dna_sequence', '')[:30]}...")
            else:
                logger.warning("\nNo weapons survived the full protocol.")
                
        else:
            logger.error("❌ Workflow timed out. No final status received.")

    except httpx.HTTPStatusError as e:
        logger.error(f"❌ A critical HTTP error occurred: {e.response.status_code} - {e.response.text}")
    except httpx.RequestError as e:
        logger.error(f"❌ A critical network error occurred: {e}")

if __name__ == "__main__":
    run_client() 
 
 
 
 