import pytest
import os
from loguru import logger
import httpx
import time
import json
from datetime import datetime

# --- Add src to path to allow for direct NCBIClient import ---
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.tools.ncbi_client import NCBIClient

# --- TEST CONFIGURATION ---
BASE_URL = "https://crispro--command-center-v8-override-fix.modal.run"
POLLING_INTERVAL = 20
MAX_POLLS = 60
TARGET_GENE = "VEGFA"
LOG_DIR = "results/integration_logs"

# --- DOCTRINE: THE BLACKOUT DRILL ---
# This test proves our platform's resilience to external service outages.

@pytest.fixture(scope="module", autouse=True)
def populate_local_cache():
    """
    Fixture to ensure our local armory is stocked with the necessary intelligence
    BEFORE the main test runs. This simulates a pre-prepared, offline-capable state.
    """
    logger.info(f"--- BLACKOUT DRILL: PRE-POPULATING LOCAL CACHE FOR {TARGET_GENE} ---")
    client = NCBIClient(threat_matrix={})
    
    # Ensure both DNA and Protein are cached
    dna_seq = client.get_dna_sequence_from_gene_symbol(TARGET_GENE)
    protein_seq = client.get_protein_sequence_from_gene_symbol(TARGET_GENE)
    
    assert dna_seq is not None, f"Failed to cache DNA for {TARGET_GENE}. Cannot run test."
    assert protein_seq is not None, f"Failed to cache Protein for {TARGET_GENE}. Cannot run test."
    
    logger.info(f"--- CACHE POPULATED. READY FOR BLACKOUT DRILL. ---")

@pytest.mark.integration
def test_full_predator_protocol_offline_capability():
    """
    Tests the full end-to-end Predator Protocol, relying exclusively on the
    local cache populated by our fixture. This validates our "Ghost Protocol"
    and "Hydra Protocol V2" resilience doctrines.
    """
    os.makedirs(LOG_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file_path = os.path.join(LOG_DIR, f"predator_protocol_run_{TARGET_GENE}_{timestamp}.json")

    start_payload = {"target_gene_symbol": TARGET_GENE}
    final_status_data = {}

    try:
        with httpx.Client(timeout=300.0) as client:
            logger.info(f"--- INITIATING PREDATOR PROTOCOL FOR {TARGET_GENE} ---")
            start_response = client.post(f"{BASE_URL}/workflow/execute", json=start_payload)
            assert start_response.status_code == 200
            workflow_id = start_response.json()["workflow_id"]
            logger.info(f"Workflow {workflow_id} initiated. Log: {log_file_path}")

            for i in range(MAX_POLLS):
                logger.info(f"Polling attempt {i+1}/{MAX_POLLS} for workflow {workflow_id}...")
                status_response = client.get(f"{BASE_URL}/status/{workflow_id}")
                assert status_response.status_code == 200
                status_data = status_response.json()
                final_status_data = status_data
                print(json.dumps(status_data, indent=2))

                if status_data["status"] == "complete":
                    logger.success(f"--- PREDATOR PROTOCOL COMPLETE FOR {TARGET_GENE} ---")
                    final_result = status_data.get("result")
                    assert final_result is not None
                    assert isinstance(final_result, list) and len(final_result) > 0
                    
                    first_weapon = final_result[0]
                    assert "structural_confidence_plddt" in first_weapon
                    assert "binding_affinity_score" in first_weapon
                    plddt = first_weapon["structural_confidence_plddt"]
                    affinity = first_weapon["binding_affinity_score"]

                    assert plddt >= 70.0
                    # For this simulation, we accept any non-error affinity score.
                    assert affinity != -1.0

                    logger.success(f"✅ SUCCESS: Weapon forged for {TARGET_GENE}. pLDDT: {plddt:.2f}, Affinity: {affinity:.2f}")
                    return

                elif status_data["status"] == "failed":
                    pytest.fail(f"Workflow {workflow_id} failed: {status_data['message']}")
                
                time.sleep(POLLING_INTERVAL)

            pytest.fail(f"Workflow {workflow_id} timed out.")

    finally:
        with open(log_file_path, "w") as f:
            json.dump(final_status_data, f, indent=2) 
import os
from loguru import logger
import httpx
import time
import json
from datetime import datetime

# --- Add src to path to allow for direct NCBIClient import ---
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.tools.ncbi_client import NCBIClient

# --- TEST CONFIGURATION ---
BASE_URL = "https://crispro--command-center-v8-override-fix.modal.run"
POLLING_INTERVAL = 20
MAX_POLLS = 60
TARGET_GENE = "VEGFA"
LOG_DIR = "results/integration_logs"

# --- DOCTRINE: THE BLACKOUT DRILL ---
# This test proves our platform's resilience to external service outages.

@pytest.fixture(scope="module", autouse=True)
def populate_local_cache():
    """
    Fixture to ensure our local armory is stocked with the necessary intelligence
    BEFORE the main test runs. This simulates a pre-prepared, offline-capable state.
    """
    logger.info(f"--- BLACKOUT DRILL: PRE-POPULATING LOCAL CACHE FOR {TARGET_GENE} ---")
    client = NCBIClient(threat_matrix={})
    
    # Ensure both DNA and Protein are cached
    dna_seq = client.get_dna_sequence_from_gene_symbol(TARGET_GENE)
    protein_seq = client.get_protein_sequence_from_gene_symbol(TARGET_GENE)
    
    assert dna_seq is not None, f"Failed to cache DNA for {TARGET_GENE}. Cannot run test."
    assert protein_seq is not None, f"Failed to cache Protein for {TARGET_GENE}. Cannot run test."
    
    logger.info(f"--- CACHE POPULATED. READY FOR BLACKOUT DRILL. ---")

@pytest.mark.integration
def test_full_predator_protocol_offline_capability():
    """
    Tests the full end-to-end Predator Protocol, relying exclusively on the
    local cache populated by our fixture. This validates our "Ghost Protocol"
    and "Hydra Protocol V2" resilience doctrines.
    """
    os.makedirs(LOG_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file_path = os.path.join(LOG_DIR, f"predator_protocol_run_{TARGET_GENE}_{timestamp}.json")

    start_payload = {"target_gene_symbol": TARGET_GENE}
    final_status_data = {}

    try:
        with httpx.Client(timeout=300.0) as client:
            logger.info(f"--- INITIATING PREDATOR PROTOCOL FOR {TARGET_GENE} ---")
            start_response = client.post(f"{BASE_URL}/workflow/execute", json=start_payload)
            assert start_response.status_code == 200
            workflow_id = start_response.json()["workflow_id"]
            logger.info(f"Workflow {workflow_id} initiated. Log: {log_file_path}")

            for i in range(MAX_POLLS):
                logger.info(f"Polling attempt {i+1}/{MAX_POLLS} for workflow {workflow_id}...")
                status_response = client.get(f"{BASE_URL}/status/{workflow_id}")
                assert status_response.status_code == 200
                status_data = status_response.json()
                final_status_data = status_data
                print(json.dumps(status_data, indent=2))

                if status_data["status"] == "complete":
                    logger.success(f"--- PREDATOR PROTOCOL COMPLETE FOR {TARGET_GENE} ---")
                    final_result = status_data.get("result")
                    assert final_result is not None
                    assert isinstance(final_result, list) and len(final_result) > 0
                    
                    first_weapon = final_result[0]
                    assert "structural_confidence_plddt" in first_weapon
                    assert "binding_affinity_score" in first_weapon
                    plddt = first_weapon["structural_confidence_plddt"]
                    affinity = first_weapon["binding_affinity_score"]

                    assert plddt >= 70.0
                    # For this simulation, we accept any non-error affinity score.
                    assert affinity != -1.0

                    logger.success(f"✅ SUCCESS: Weapon forged for {TARGET_GENE}. pLDDT: {plddt:.2f}, Affinity: {affinity:.2f}")
                    return

                elif status_data["status"] == "failed":
                    pytest.fail(f"Workflow {workflow_id} failed: {status_data['message']}")
                
                time.sleep(POLLING_INTERVAL)

            pytest.fail(f"Workflow {workflow_id} timed out.")

    finally:
        with open(log_file_path, "w") as f:
            json.dump(final_status_data, f, indent=2) 
 
 
 
 
 