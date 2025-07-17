import pytest
import os
import sys
import httpx
import asyncio
import logging

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

# --- Test Configuration ---
# This is the real, live URL for our DEPLOYED V3 service
COMMAND_CENTER_URL = "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run"
CLIENT_TIMEOUT = 900.0  # Increased timeout for the full, multi-step workflow
POLLING_INTERVAL = 30   # Seconds to wait between status checks
MAX_POLLS = 30          # Max number of polls (30 * 30s = 15 minutes)

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@pytest.mark.integration
# The test is parameterized to run with asyncio, which is all we need.
# The trio backend is not configured in this environment and causes a failure.
@pytest.mark.asyncio
async def test_v3_full_workflow_integration():
    """
    This is a true end-to-end integration test for the V3 asynchronous workflow.
    It performs the following steps:
    1. Calls POST /v3/workflow/formulate_battle_plan to create a new plan.
    2. Calls POST /v3/workflow/execute_guide_design_campaign to start the async process.
    3. Polls GET /v3/battle_plan/{plan_id} until the status is 'interventions_designed'.
    
    NOTE: This test requires both CommandCenter and EvoService to be deployed and reachable.
    It will be slow due to model loading and computation.
    """
    async with httpx.AsyncClient(timeout=CLIENT_TIMEOUT, verify=False) as client:
        # --- Step 1: Formulate the Battle Plan ---
        formulate_url = f"{COMMAND_CENTER_URL}/v3/workflow/formulate_battle_plan"
        # Using RUNX1 test data
        request_payload = {
            "patient_identifier": "integration_test_patient_v3_01",
            "gene": "RUNX1",
            "mutation_hgvs_c": "c.587G>A",
            "mutation_hgvs_p": "p.R196Q",
            "transcript_id": "NM_001754.4",
            "bam_file_path": "/dev/null",
            "protein_sequence": "MSDSILAEYEYPHDFVNSEEEQRAKRSPIRVADPVAFVADGERAREGSPLLQLPGQPAPGALSKMIKRPRTOPV",
            "sequence_for_perplexity": "CCTGCCTGGGCTCCCCAGCCCCAGCAGCCCTGCAGCCCAGCTCTGGCCTGCCGGAGGCCACGCGC"
        }

        logger.info(f"\n🚀 STEP 1: Calling Formulate Battle Plan -> {formulate_url}")
        response = await client.post(formulate_url, json=request_payload)
        logger.info(f"📦 Response status: {response.status_code}")
        logger.info(f"📄 Response JSON: {response.json()}")
        
        response.raise_for_status()
        assert response.status_code == 201
        
        formulate_data = response.json()
        plan_id = formulate_data.get("plan_id")
        assert plan_id is not None
        logger.info(f"✅ STEP 1 COMPLETE: Battle Plan created with ID: {plan_id}")

        # --- Step 2: Execute Guide Design Campaign ---
        execute_url = f"{COMMAND_CENTER_URL}/v3/workflow/execute_guide_design_campaign/{plan_id}"
        logger.info(f"\n🚀 STEP 2: Calling Execute Guide Design -> {execute_url}")
        response = await client.post(execute_url)
        logger.info(f"📦 Response status: {response.status_code}")
        logger.info(f"📄 Response JSON: {response.json()}")

        response.raise_for_status()
        assert response.status_code == 202
        logger.info(f"✅ STEP 2 COMPLETE: Guide design campaign initiated.")

        # --- Step 3: Poll for Results ---
        status_url = f"{COMMAND_CENTER_URL}/v3/battle_plan/{plan_id}"
        logger.info(f"\n🚀 STEP 3: Polling for results -> {status_url}")
        
        final_status = ""
        for i in range(MAX_POLLS):
            logger.info(f"  - Poll attempt {i+1}/{MAX_POLLS}...")
            await asyncio.sleep(POLLING_INTERVAL)
            response = await client.get(status_url)
            
            if response.status_code == 200:
                status_data = response.json()
                current_status = status_data.get("status")
                logger.info(f"  - Current status: '{current_status}'")
                if current_status == "interventions_designed":
                    final_status = current_status
                    final_results = status_data.get("results")
                    logger.info(f"  - SUCCESS! Results are in: {final_results}")
                    # --- ADDING MORE STRINGENT CHECK ---
                    # We must verify that the BLAST check actually ran and didn't return a dummy value.
                    assert final_results["guides"][0]["off_target_count"] != 999, "BLAST service failed! Off-target count is a dummy value."
                    break
                elif current_status == "design_failed":
                    final_status = current_status
                    logger.error(f"  - FAILURE! Design process failed: {status_data.get('results')}")
                    break
            else:
                logger.warning(f"  - Received non-200 status code: {response.status_code}")
        
        assert final_status == "interventions_designed", f"Workflow did not complete successfully. Final status was '{final_status}'"
        logger.info(f"\n✅ STEP 3 COMPLETE: Workflow finished with status 'interventions_designed'.")
        logger.info("🎉 Integration test passed! 🎉") 