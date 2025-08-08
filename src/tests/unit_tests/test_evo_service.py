import pytest
import httpx
import time
import asyncio
import json

# NOTE: This is a live integration test and requires the evo-service to be deployed on Modal.
EVO_SERVICE_URL = "https://crispro--evo-service-serve.modal.run"
BRAF_V600E_TARGET_SEQUENCE = "AGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATGGTA"

@pytest.mark.asyncio
@pytest.mark.integration
async def test_generate_guide_rna_endpoint():
    """
    Tests the full asynchronous guide RNA generation workflow:
    1. Submits a job to /generate_guide_rna.
    2. Polls the /jobs/{job_id} endpoint until completion.
    3. Verifies the final result.
    """
    async with httpx.AsyncClient() as client:
        # Step 1: Submit the job
        submit_payload = {
            "target_sequence": BRAF_V600E_TARGET_SEQUENCE,
            "num_guides": 3,
        }
        submit_response = await client.post(
            f"{EVO_SERVICE_URL}/generate",
            json=submit_payload,
            timeout=60.0,
        )
        assert submit_response.status_code == 202
        job_id = submit_response.json()["job_id"]
        print(f"\n✅ Job submitted successfully. Job ID: {job_id}")

        # Step 2: Poll for the result
        start_time = time.time()
        timeout_seconds = 1800  # 30-minute overall timeout for the first run
        poll_interval_seconds = 15

        while time.time() - start_time < timeout_seconds:
            status_response = await client.get(f"{EVO_SERVICE_URL}/status/{job_id}")
            
            if status_response.status_code == 200:
                status_data = status_response.json()
                current_status = status_data.get("status")
                print(f"Polling... current status: {current_status}")

                if current_status == "complete":
                    # Step 3: Verify the final result
                    assert "result" in status_data
                    guides = status_data["result"]
                    assert isinstance(guides, list)
                    assert len(guides) > 0
                    assert all(isinstance(g, str) for g in guides)
                    print(f"✅ SUCCESS: Job {job_id} completed. Received {len(guides)} guide RNA candidates.")
                    return  # Test success
                
                elif current_status == "failed":
                    pytest.fail(f"Job {job_id} failed with error: {status_data.get('error')}")

            await asyncio.sleep(poll_interval_seconds)

        pytest.fail(f"Polling timed out after {timeout_seconds} seconds. The job did not complete.") 