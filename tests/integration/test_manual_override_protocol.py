import pytest
import httpx
import time
import os
import json
from datetime import datetime

# Target the new V8 "Scorched Earth" deployment
BASE_URL = "https://crispro--command-center-v8-override-fix.modal.run"
POLLING_INTERVAL = 20  # seconds
MAX_POLLS = 60  # 20s * 60 = 1200s or 20 minutes max wait time
LOG_DIR = "results/integration_logs"

# --- OPERATION: BRCA1 GAUNTLET ---
# This new bait sequence targets a critical splice site in BRCA1 (Exon 11 donor).
# It is a high-value, doctrinally-sound target that plays to Evo2's strengths.
BRCA1_SPLICE_SITE_BAIT = "TGGATACTGAAACACATTAGAAGAAATGACTGAACTGTCTAAAGCTGTAGTAAATATACAATGCACTTTCTTCCAAAATGTAGAAAAGGCTGAATTCATCTAAAAAGAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCT"

@pytest.mark.integration
def test_full_brca1_gauntlet_workflow():
    """
    Tests the full end-to-end "Focused Forge" protocol against a high-value
    BRCA1 splice site target.
    """
    # --- 1. Setup Logging ---
    os.makedirs(LOG_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file_path = os.path.join(LOG_DIR, f"brca1_gauntlet_run_{timestamp}.json")

    start_payload = {"bait_sequence": BRCA1_SPLICE_SITE_BAIT}
    final_status_data = {}

    try:
        # Increased timeout to 300s (5 mins) to account for Modal cold starts.
        with httpx.Client(timeout=300.0) as client:
            # --- 2. Initiate the Workflow ---
            print(f"Initiating BRCA1 Gauntlet workflow...")
            start_response = client.post(f"{BASE_URL}/workflow/execute", json=start_payload)
            assert start_response.status_code == 200
            workflow_id = start_response.json()["workflow_id"]
            print(f"Workflow {workflow_id} initiated successfully. Log file: {log_file_path}")

            # --- 3. Poll for Status ---
            for i in range(MAX_POLLS):
                print(f"Polling attempt {i+1}/{MAX_POLLS} for workflow {workflow_id}...")
                status_response = client.get(f"{BASE_URL}/status/{workflow_id}")
                assert status_response.status_code == 200, f"Failed to get status (HTTP {status_response.status_code})"
                status_data = status_response.json()
                final_status_data = status_data # Continuously update the final log data
                
                # More verbose console logging
                print(json.dumps(status_data, indent=2))

                if status_data["status"] == "complete":
                    # --- 4. Final Validation ---
                    print("Workflow complete! Validating results...")
                    final_result = status_data.get("result")
                    assert final_result is not None, "Result is missing from the final status."
                    assert isinstance(final_result, list), "Result should be a list."
                    assert len(final_result) > 0, "Expected at least one validated weapon."
                    
                    first_weapon = final_result[0]
                    assert "structural_confidence_plddt" in first_weapon
                    plddt = first_weapon["structural_confidence_plddt"]
                    assert plddt >= 70.0, f"Weapon pLDDT score {plddt} is below the 70.0 threshold."
                    
                    print(f"✅ SUCCESS: Workflow {workflow_id} produced a validated weapon with pLDDT: {plddt:.2f}")
                    return # Test passed

                elif status_data["status"] == "failed":
                    pytest.fail(f"Workflow {workflow_id} failed with message: {status_data['message']}")
                
                time.sleep(POLLING_INTERVAL)

            pytest.fail(f"Workflow {workflow_id} timed out after {MAX_POLLS * POLLING_INTERVAL} seconds.")

    finally:
        # --- 5. Write to Log File ---
        print(f"Writing final status to log file: {log_file_path}")
        with open(log_file_path, "w") as f:
            json.dump(final_status_data, f, indent=2) 
 
 
 
 
 
 