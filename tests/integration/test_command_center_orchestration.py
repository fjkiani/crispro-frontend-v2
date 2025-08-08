import pytest
import httpx
import time
import os

COMMAND_CENTER_URL = "https://crispro--command-center-v6-orchestrator-web-app.modal.run"

@pytest.mark.integration
def test_full_predator_protocol_orchestration():
    """
    This is the ultimate integration test. It validates the entire kill chain
    by issuing a high-level command to the CommandCenter and observing the result.
    It confirms that our Triumvirate architecture works in concert.
    """
    target_gene = "MMP9" # Using a known good target for validation
    
    with httpx.Client(timeout=900.0) as client:
        # 1. Initiate the workflow
        print(f"\n--- Dispatching command to CommandCenter: Engage target {target_gene} ---")
        initiate_response = client.post(
            f"{COMMAND_CENTER_URL}/workflow/execute_predator_protocol",
            json={"target_gene_symbol": target_gene}
        )
        assert initiate_response.status_code == 200
        job_id = initiate_response.json().get("job_id")
        assert job_id is not None
        print(f"--- CommandCenter Acknowledged. Job ID: {job_id} ---")

        # 2. Poll for the final result
        while True:
            time.sleep(10)
            status_response = client.get(f"{COMMAND_CENTER_URL}/status/{job_id}")
            assert status_response.status_code == 200, f"Status check failed: {status_response.text}"
            
            job_data = status_response.json()
            status = job_data.get("status")
            message = job_data.get("message")
            
            print(f"POLLING... STATUS: {status.upper()} | MESSAGE: {message}")
            
            if status == "complete":
                print("--- MISSION COMPLETE ---")
                final_result = job_data.get("result")
                assert final_result is not None
                assert isinstance(final_result, list)
                assert len(final_result) > 0
                
                # Validate the structure of the final blueprint
                blueprint = final_result[0]
                assert "inhibitor_sequence" in blueprint
                assert "stability_score" in blueprint
                assert "commentary" in blueprint
                
                print("--- Final Weapon Blueprint Validated ---")
                print(final_result)
                break
            
            elif status == "failed":
                pytest.fail(f"The orchestrated workflow failed with message: {message}")
        
        print(f"--- Full Triumvirate Orchestration for {target_gene} VERIFIED ---") 