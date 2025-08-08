import sys
import os

# This ensures the project's src directory is on the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

import requests
import time
import json
from datetime import datetime

from tools.candidate_triage import CandidateTriage

# --- Configuration ---
ZETA_FORGE_URL = "https://crispro--zeta-forge-v1-api.modal.run"
# This is the full VEGFA protein sequence, a more robust target
TARGET_PROTEIN_SEQUENCE = "MNFLLSWVHWSLALLLYLHHAKWSQAAPMAEGGGQNHHEVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDESLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPKKDRARQENCDKPRR"


# Normally, this would be fetched from a database or a centralized config
# For this test script, we define it directly.
THREAT_MATRIX = {
    'high_value_keywords': [
        "binding pocket", "active site", "allosteric site", 
        "protein-protein interface", "catalytic domain", "receptor binding domain",
        "variable loop", "complementarity-determining region (CDR)", "paratope"
    ],
    'low_value_keywords': [
        "repetitive sequence", "low complexity", "dna sequence", 
        "artifact", "non-amino acid characters", "gtp-binding", "atp-binding"
    ]
}


def run_forge_and_triage(target_protein_id: str, num_candidates: int, output_dir: str = "results/forge_triage"):
    """
    Runs the full workflow by correctly interacting with the async Zeta-Forge API.
    1. Submits a job to generate candidate sequences.
    2. Polls the status endpoint until the job is complete.
    3. Retrieves the candidates and triages them.
    4. Saves the viable candidates and the full triage report.
    """
    print(f"--- Starting Operation: Forge & Triage for target: {target_protein_id} ---")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # === STAGE 1: SUBMIT FORGE JOB ===
    job_id = None
    try:
        print("\n==================================================")
        print("           SUBMITTING JOB TO ZETA-FORGE")
        print("==================================================")
        submit_response = requests.post(
            f"{ZETA_FORGE_URL}/generate_inhibitor",
            json={
                # --- NEW DOCTRINE: Conceptual Prompting ---
                "target_protein_id": target_protein_id, 
                # "target_protein_sequence": TARGET_PROTEIN_SEQUENCE, # This is now legacy
                "num_candidates_per_temp": num_candidates // 3 or 1, # Distribute among temps
                "temperatures": [0.7, 1.0, 1.2],
                "generation_length": 150 
            },
            timeout=120
        )
        submit_response.raise_for_status()
        job_id = submit_response.json()["job_id"]
        print(f"SUCCESS: Job submitted. Job ID: {job_id}")

    except requests.exceptions.RequestException as e:
        print(f"FATAL: Failed to submit job to Zeta-Forge: {e}")
        return

    # === STAGE 2: POLL FOR JOB COMPLETION ===
    candidates = []
    start_time = time.time()
    try:
        print("\n==================================================")
        print("           AWAITING FORGE COMPLETION")
        print("==================================================")
        while True:
            status_response = requests.get(f"{ZETA_FORGE_URL}/status/{job_id}", timeout=60)
            status_response.raise_for_status()
            status_data = status_response.json()
            
            current_status = status_data.get("status")
            message = status_data.get("message", "No message.")
            print(f"[{datetime.now().strftime('%H:%M:%S')}] Status: {current_status} - {message}")

            if current_status == "complete":
                candidates = status_data.get("candidates", [])
                print(f"\nSUCCESS: Forge operation complete. Retrieved {len(candidates)} candidates.")
                break
            elif current_status == "failed":
                raise RuntimeError(f"Forge job failed: {status_data.get('error', 'Unknown error')}")

            time.sleep(15) # Wait 15 seconds before polling again

    except requests.exceptions.RequestException as e:
        print(f"FATAL: Failed to get job status from Zeta-Forge: {e}")
        return
    except RuntimeError as e:
        print(f"FATAL: {e}")
        return

    forge_duration = time.time() - start_time
    print(f"Total forge time: {forge_duration:.2f} seconds.")

    # === STAGE 3: TRIAGE CANDIDATES ===
    if not candidates:
        print("\nWARNING: No candidates were generated. Skipping triage.")
        return
        
    print("\n==================================================")
    print("           INITIATING TRIAGE PROTOCOL")
    print("==================================================")
    start_triage_time = time.time()

    triage_tool = CandidateTriage(threat_matrix=THREAT_MATRIX)
    triage_report = triage_tool.run_triage(candidates)
    
    triage_duration = time.time() - start_triage_time
    print(f"\nTriage protocol completed in {triage_duration:.2f} seconds.")

    # === STAGE 3: SAVE RESULTS ===
    viable_candidates = [
        item for item in triage_report if item.get('rating') == 'VIABLE'
    ]

    print("\n==================================================")
    print("           TRIAGE SUMMARY")
    print("==================================================")
    print(f"Total Candidates Generated: {len(candidates)}")
    print(f"Candidates Passing Low-Complexity Filter: {len([c for c in triage_report])}")
    print(f"Candidates Rated VIABLE by AI Analyst: {len(viable_candidates)}")
    print(f"Candidates REJECTED/ERROR: {len(triage_report) - len(viable_candidates)}")
    
    # Save the full triage report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_filename = os.path.join(output_dir, f"triage_report_{target_protein_id}_{timestamp}.json")
    
    with open(report_filename, 'w') as f:
        json.dump(triage_report, f, indent=4)
        
    print(f"\nFull triage report saved to: {report_filename}")

    # Save just the viable candidates' sequences for easy access
    if viable_candidates:
        viable_filename = os.path.join(output_dir, f"viable_candidates_{target_protein_id}_{timestamp}.json")
        with open(viable_filename, 'w') as f:
            # Storing them in a simple list within a json object
            json.dump({"viable_sequences": [vc['candidate'] for vc in viable_candidates]}, f, indent=4)
        print(f"Viable candidates saved to: {viable_filename}")

    print(f"\n--- Operation Forge & Triage Complete ---")


if __name__ == "__main__":
    # Example command to run the test
    TARGET_PROTEIN_ID = "VEGFA" 
    NUM_CANDIDATES_TO_GENERATE = 9 # Use a number divisible by 3 for even distribution

    run_forge_and_triage(
        target_protein_id=TARGET_PROTEIN_ID,
        num_candidates=NUM_CANDIDATES_TO_GENERATE
    ) 