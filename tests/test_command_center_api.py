import pytest
import requests
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

@pytest.fixture
def command_center_url():
    """Fixture to provide the CommandCenter URL from environment variables."""
    url = os.getenv("COMMAND_CENTER_URL")
    if not url:
        pytest.fail("COMMAND_CENTER_URL environment variable is not set. Cannot run integration tests.")
    return url.rstrip('/')

def test_health_check(command_center_url):
    """Tests if the CommandCenter health check endpoint is responsive."""
    response = requests.post(f"{command_center_url}/health")
    assert response.status_code == 200
    assert response.json() == {"status": "healthy", "service": "command-center"}

def test_assess_threat_endpoint_frameshift(command_center_url):
    """
    Tests the /workflow/assess_threat endpoint with a known frameshift mutation.
    This should trigger the 'Truncation Sieve'.
    """
    endpoint = "/workflow/assess_threat"
    payload = {
        "gene_symbol": "ASXL1",
        "protein_change": "p.Gly646fs"
    }
    response = requests.post(f"{command_center_url}{endpoint}", json=payload)
    
    assert response.status_code == 200
    data = response.json()
    
    assert data['request']['gene_symbol'] == "ASXL1"
    assert data['assessment']['source'] == "Truncation Sieve"
    assert data['assessment']['is_truncating'] is True
    assert "vep" in data['external_intelligence']

def test_assess_threat_endpoint_missense(command_center_url):
    """
    Tests the /workflow/assess_threat endpoint with a known missense mutation.
    This should trigger the 'Zeta Oracle'.
    """
    endpoint = "/workflow/assess_threat"
    payload = {
        "gene_symbol": "BRAF",
        "protein_change": "p.Val600Glu" # V600E
    }
    response = requests.post(f"{command_center_url}{endpoint}", json=payload)
    
    assert response.status_code == 200
    data = response.json()
    
    assert data['request']['gene_symbol'] == "BRAF"
    assert data['assessment']['source'] == "Zeta Oracle"
    assert data['assessment']['is_truncating'] is False
    assert data['assessment']['zeta_score'] is not None
    assert "vep" in data['external_intelligence']

@pytest.mark.skip(reason="Seed & Soil workflow has placeholder components and is too slow for unit tests.")
def test_seed_and_soil_full_workflow(command_center_url):
    """
    An integration test to validate the full Seed & Soil workflow chain.
    This test is currently skipped as it's long-running and depends on workflows
    that may have their own placeholders (e.g., essentiality prediction).
    """
    # This would be a more complex test simulating the UI's behavior.
    # 1. Call run_gene_essentiality_prediction
    # 2. Call simulate_environmental_stress with the result of step 1
    # 3. Call run_threat_report_generation with the results of steps 1 and 2
    # 4. Assert the final report is valid
    pass 