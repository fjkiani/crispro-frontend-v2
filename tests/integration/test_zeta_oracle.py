import pytest
import requests
import os

# Get the base URL for the ZetaOracle from an environment variable, with a fallback
# NOTE: This assumes the oracle will be deployed with this app name.
# This may need to be updated after the first deployment.
BASE_URL = os.environ.get("ZETA_ORACLE_URL", "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run")

@pytest.fixture(scope="module")
def zeta_oracle_service():
    """
    A pytest fixture to check that the ZetaOracle service is online before running tests.
    """
    try:
        # A simple GET to the root of a FastAPI app should not be a valid endpoint,
        # so we expect a 404, which still proves the service is running and responsive.
        response = requests.get(BASE_URL)
        assert response.status_code == 404, f"Expected 404 from root, but got {response.status_code}"
    except requests.exceptions.RequestException as e:
        pytest.fail(f"ZetaOracle service is not responding at {BASE_URL}. Error: {e}")
    
    return BASE_URL

def test_score_variant_ensemble_success(zeta_oracle_service):
    """
    Tests the /score_variant_ensemble endpoint with a valid batch request.
    Verifies the structure of the response and the fused zeta_score.
    """
    base_url = zeta_oracle_service
    endpoint = f"{base_url}/score_variant_ensemble"

    # A batch of two variants to score
    payload = {
        "variants": [
            {
                "variant_id": "TEST_VAR_1",
                "alphamissense_variant_str": "chr7:g.140753336A>T", # BRAF V600E like
                "mutation_hgvs": "V600E",
                "protein_sequence": "PROTEINSEQ1"
            },
            {
                "variant_id": "TEST_VAR_2",
                "alphamissense_variant_str": "chr1:g.12345A>G", # Benign-like
                "mutation_hgvs": "A123T",
                "protein_sequence": "PROTEINSEQ2"
            }
        ]
    }

    print(f"Submitting batch scoring request to {endpoint}...")
    response = requests.post(endpoint, json=payload)
    
    assert response.status_code == 200, f"Request failed: {response.status_code} - {response.text}"
    
    data = response.json()
    assert "scored_variants" in data
    
    scored_variants = data["scored_variants"]
    assert isinstance(scored_variants, list)
    assert len(scored_variants) == 2

    # Validate the structure of the first scored variant
    first_variant = scored_variants[0]
    assert first_variant["variant_id"] == "TEST_VAR_1"
    assert "alphamissense_score" in first_variant
    assert "esm_score" in first_variant
    assert "evo2_score" in first_variant
    assert "zeta_score" in first_variant
    
    zeta_score = first_variant["zeta_score"]
    assert isinstance(zeta_score, (float, int))

    print("Successfully received and validated the ensemble response. ✅")

def test_score_variant_ensemble_empty_request(zeta_oracle_service):
    """
    Tests that the endpoint handles a request with an empty list of variants.
    """
    base_url = zeta_oracle_service
    endpoint = f"{base_url}/score_variant_ensemble"
    
    payload = {"variants": []}
    
    response = requests.post(endpoint, json=payload)
    assert response.status_code == 200
    
    data = response.json()
    assert "scored_variants" in data
    assert len(data["scored_variants"]) == 0

def test_score_variant_ensemble_bad_payload(zeta_oracle_service):
    """
    Tests that the endpoint returns a 422 Unprocessable Entity for a malformed payload.
    """
    base_url = zeta_oracle_service
    endpoint = f"{base_url}/score_variant_ensemble"
    
    # Malformed payload (e.g., missing 'variants' key)
    payload = {"not_variants": []}
    
    response = requests.post(endpoint, json=payload)
    assert response.status_code == 422 