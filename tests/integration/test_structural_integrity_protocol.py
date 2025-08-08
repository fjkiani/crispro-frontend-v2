import pytest
import httpx
import time
import os

# --- Test Configuration ---
# It's better to fetch the URL from an environment variable if available,
# allowing for flexibility between local and staging environments.
BASE_URL = "https://crispro--command-center-v6-orchestrator-web-app.modal.run"
POLLING_INTERVAL = 15  # seconds
MAX_POLLS = 40  # 15s * 40 = 600s or 10 minutes max wait time

@pytest.mark.integration
def test_full_structural_integrity_workflow():
    """
    Tests the full end-to-end Structural Integrity Protocol.
    It initiates the workflow, polls for completion, and validates the output.
    This test is designed to catch "wet noodle" failures by asserting
    on the final pLDDT score.
    """
    target_gene = "BRAF"  # TACTICAL PIVOT: Using a gene we have in our local database
    start_payload = {"target_gene_symbol": target_gene}
    
    with httpx.Client(timeout=60.0) as client:
        # --- 1. Initiate the Protocol ---
        start_url = f"{BASE_URL}/workflow/execute_structural_integrity_protocol"
        print(f"\n[TEST] Initiating workflow for gene: {target_gene} at {start_url}")
        
        try:
            start_response = client.post(start_url, json=start_payload)
            start_response.raise_for_status()
            start_data = start_response.json()
            job_id = start_data.get("job_id")
            
            assert job_id is not None, "Failed to get a job_id from the command center."
            print(f"[TEST] Workflow initiated successfully. Job ID: {job_id}")

        except httpx.HTTPStatusError as e:
            pytest.fail(f"HTTP error during workflow initiation: {e.response.status_code} - {e.response.text}")
        except httpx.RequestError as e:
            pytest.fail(f"Request error during workflow initiation. Is the service running at {BASE_URL}? Error: {e}")

        # --- 2. Poll for Status ---
        status_url = f"{BASE_URL}/status/{job_id}"
        for i in range(MAX_POLLS):
            print(f"[TEST] Polling status ({i+1}/{MAX_POLLS})...")
            try:
                status_response = client.get(status_url)
                status_response.raise_for_status()
                status_data = status_response.json()
                
                status = status_data.get("status")
                message = status_data.get("message")
                print(f"[TEST]   Status: {status} - {message}")

                if status == "awaiting_manual_gauntlet":
                    # --- 3. The Final Verdict ---
                    print("[TEST] Workflow complete. Verifying qualified candidates.")
                    final_result = status_data.get("result")
                    
                    assert final_result is not None, "Workflow completed but returned no result."
                    assert isinstance(final_result, list), "Result should be a list of candidates."
                    assert len(final_result) > 0, "Workflow completed but returned an empty list of candidates."
                    
                    top_candidate = final_result[0]
                    stability_score = top_candidate.get("stability_score")
                    
                    assert stability_score is not None, "Top candidate is missing a stability_score."
                    assert "protein_sequence" in top_candidate, "Top candidate is missing a protein_sequence."
                    assert len(top_candidate["protein_sequence"]) > 0, "Top candidate has an empty protein_sequence."
                    
                    # THE CRITICAL ASSERTION: The Sieve is working correctly.
                    assert stability_score >= -5.0, f"Validation Failed: Top candidate stability score ({stability_score}) is below the threshold of -5.0."
                    
                    print(f"✅ VALIDATION SUCCESS: Top candidate passed the Sieve with stability score: {stability_score}")
                    return # Test success

                elif status == "failed":
                    pytest.fail(f"Workflow failed with message: {message}")

            except httpx.HTTPStatusError as e:
                pytest.fail(f"HTTP error during status polling: {e.response.status_code} - {e.response.text}")
            except httpx.RequestError as e:
                pytest.fail(f"Request error during polling. Service might have crashed. Error: {e}")

            time.sleep(POLLING_INTERVAL)

        pytest.fail("Polling timed out. The workflow took too long to complete.")

if __name__ == "__main__":
    # Allows running this test directly for debugging.
    # Make sure the CommandCenter service is running first.
    # e.g., `modal serve src/services/command_center/main.py`
    test_full_structural_integrity_workflow() 