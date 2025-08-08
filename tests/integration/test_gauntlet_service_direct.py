import pytest
import requests
import os
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# The URL for the deployed Gauntlet service.
# This must be set as an environment variable.
# Find the URL on your Modal deployment page: https://modal.com/apps
# Example: export GAUNTLET_SERVICE_URL="https://<your-modal-user>--gauntlet-service-gauntlet-fastapi-app.modal.run"
SERVICE_URL = "https://crispro--gauntlet-service-gauntlet-fastapi-app.modal.run"

# Using the same test sequence from the ColabFold notebook for consistency
TEST_PROTEIN_SEQUENCE = "PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASKK"

# Short test sequence for debugging (only 10 amino acids)
TEST_SHORT_SEQUENCE = "MKLLISGLQG"

@pytest.mark.skipif(not SERVICE_URL, reason="GAUNTLET_SERVICE_URL environment variable not set")
@pytest.mark.integration
def test_gauntlet_debug_short():
    """
    Debug test with a very short protein sequence.
    """
    endpoint_url = f"{SERVICE_URL}/predict_structure"
    logger.info(f"Testing endpoint: {endpoint_url}")

    payload = {
        "protein_sequence": TEST_SHORT_SEQUENCE
    }

    # Shorter timeout for debugging
    response = requests.post(endpoint_url, json=payload, timeout=300)

    logger.info(f"Response status: {response.status_code}")
    logger.info(f"Response text: {response.text}")
    
    assert response.status_code == 200, f"Request failed with status {response.status_code}: {response.text}" 
    # Set a very long timeout because folding can take several minutes
    response = requests.post(endpoint_url, json=payload, timeout=1800)
    
    assert response.status_code == 200, f"Request failed with status {response.status_code}: {response.text}"
    
    data = response.json()
    
    assert "plddt_score" in data
    assert "pdb_string" in data
    
    assert isinstance(data["plddt_score"], float)
    assert isinstance(data["pdb_string"], str)
    
    assert data["plddt_score"] > 0.0, "pLDDT score should be greater than 0"
    assert "ATOM" in data["pdb_string"], "Response should contain a valid PDB string"

    logger.info(f"Test successful! Received pLDDT score: {data['plddt_score']}") 