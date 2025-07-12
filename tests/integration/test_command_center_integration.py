import pytest
import os
import sys
import httpx

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from services.command_center.schemas import PatientDataPacket, BattlePlanResponse

# This is the real, live URL for our deployed service
COMMAND_CENTER_URL = "https://crispro--crispr-assistant-command-center-v2-commandcenter-api.modal.run"

# We need a long timeout because these are real, cold-starting services
CLIENT_TIMEOUT = 300.0

@pytest.mark.integration
@pytest.mark.anyio
async def test_reconnaissance_phase_integration():
    """
    This is a true integration test. It calls the live, deployed CommandCenter 
    endpoint, which in turn calls the live, deployed ZetaOracle.
    
    NOTE: This test requires both services to be deployed and reachable.
    It will be slow on the first run due to model loading.
    """
    endpoint_url = f"{COMMAND_CENTER_URL}/v2/workflow/formulate_battle_plan"

    # We need to provide real sequences for the oracle to score.
    # These are for BRAF V600E, a well-known mutation.
    # In a real app, we'd fetch these from a database like RefSeq.
    reference_dna_sequence = "CTGTAGCTAGCAGAAATCT" # Simplified sequence around the mutation
    alternate_dna_sequence = "CTGTAGCTAGCagAAATCT" # V600E is a T>A substitution

    request_payload = {
        "patient_identifier": "integration_test_patient_01",
        "gene": "BRAF",
        "mutation_hgvs_c": "c.1799T>A",
        "mutation_hgvs_p": "p.Val600Glu",
        "transcript_id": "NM_004333.3",
        "bam_file_path": "/dev/null", # Ignored by mock
        "protein_sequence": alternate_dna_sequence, # Using this field for alternate sequence
        "sequence_for_perplexity": reference_dna_sequence, # Using this field for reference sequence
    }

    async with httpx.AsyncClient(timeout=CLIENT_TIMEOUT) as client:
        print(f"\\n🚀 Calling live endpoint: {endpoint_url}")
        response = await client.post(endpoint_url, json=request_payload)

        print(f"📦 Response status: {response.status_code}")
        print(f"📄 Response JSON: {response.json()}")
        
        response.raise_for_status()
        
        # Assertions to validate the response structure
        response_data = response.json()
        assert "plan_id" in response_data
        assert "message" in response_data
        assert response_data["message"] == "Battle plan formulated successfully."

    # Future step: We could add a subsequent check to query the DB via another
    # endpoint to confirm the zeta_score was stored correctly. For now, a 200 OK
    # from the full pipeline is our success metric. 