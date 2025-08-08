import pytest
import httpx
import os

# The URL for our deployed service
HUNTER_ANALYST_URL = "https://crispro--hunter-analyst-service-hunteranalyst-hunt.modal.run"

@pytest.mark.integration
def test_hunter_analyst_hunt_success():
    """
    Tests the /hunt endpoint of the Hunter-Analyst service with a valid gene symbol.
    This validates the "happy path" for our reconnaissance service.
    """
    target_gene = "MMP9"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # 1. Assert successful response
        assert response.status_code == 200, f"Expected status 200, but got {response.status_code}: {response.text}"
        
        # 2. Assert the response payload structure
        data = response.json()
        expected_keys = ["gene_symbol", "domain_name", "target_motif_dna", "start_aa", "end_aa"]
        for key in expected_keys:
            assert key in data, f"Response missing expected key: {key}"
        
        # 3. Assert data integrity
        assert data["gene_symbol"] == target_gene
        assert isinstance(data["target_motif_dna"], str)
        assert len(data["target_motif_dna"]) > 0
        assert isinstance(data["start_aa"], int)
        assert isinstance(data["end_aa"], int)

@pytest.mark.integration
def test_hunter_analyst_hunt_failure_not_found():
    """
    Tests the /hunt endpoint with a gene symbol that is unlikely to be found.
    This validates the service's ability to handle bad input gracefully.
    """
    target_gene = "FAKEGENE123XYZ"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # The service should return a 404 Not Found
        assert response.status_code == 404, f"Expected status 404, but got {response.status_code}: {response.text}"
        data = response.json()
        assert "detail" in data
        assert "Could not retrieve DNA sequence" in data["detail"] 
import httpx
import os

# The URL for our deployed service
HUNTER_ANALYST_URL = "https://crispro--hunter-analyst-service-hunteranalyst-hunt.modal.run"

@pytest.mark.integration
def test_hunter_analyst_hunt_success():
    """
    Tests the /hunt endpoint of the Hunter-Analyst service with a valid gene symbol.
    This validates the "happy path" for our reconnaissance service.
    """
    target_gene = "MMP9"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # 1. Assert successful response
        assert response.status_code == 200, f"Expected status 200, but got {response.status_code}: {response.text}"
        
        # 2. Assert the response payload structure
        data = response.json()
        expected_keys = ["gene_symbol", "domain_name", "target_motif_dna", "start_aa", "end_aa"]
        for key in expected_keys:
            assert key in data, f"Response missing expected key: {key}"
        
        # 3. Assert data integrity
        assert data["gene_symbol"] == target_gene
        assert isinstance(data["target_motif_dna"], str)
        assert len(data["target_motif_dna"]) > 0
        assert isinstance(data["start_aa"], int)
        assert isinstance(data["end_aa"], int)

@pytest.mark.integration
def test_hunter_analyst_hunt_failure_not_found():
    """
    Tests the /hunt endpoint with a gene symbol that is unlikely to be found.
    This validates the service's ability to handle bad input gracefully.
    """
    target_gene = "FAKEGENE123XYZ"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # The service should return a 404 Not Found
        assert response.status_code == 404, f"Expected status 404, but got {response.status_code}: {response.text}"
        data = response.json()
        assert "detail" in data
        assert "Could not retrieve DNA sequence" in data["detail"] 
import httpx
import os

# The URL for our deployed service
HUNTER_ANALYST_URL = "https://crispro--hunter-analyst-service-hunteranalyst-hunt.modal.run"

@pytest.mark.integration
def test_hunter_analyst_hunt_success():
    """
    Tests the /hunt endpoint of the Hunter-Analyst service with a valid gene symbol.
    This validates the "happy path" for our reconnaissance service.
    """
    target_gene = "MMP9"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # 1. Assert successful response
        assert response.status_code == 200, f"Expected status 200, but got {response.status_code}: {response.text}"
        
        # 2. Assert the response payload structure
        data = response.json()
        expected_keys = ["gene_symbol", "domain_name", "target_motif_dna", "start_aa", "end_aa"]
        for key in expected_keys:
            assert key in data, f"Response missing expected key: {key}"
        
        # 3. Assert data integrity
        assert data["gene_symbol"] == target_gene
        assert isinstance(data["target_motif_dna"], str)
        assert len(data["target_motif_dna"]) > 0
        assert isinstance(data["start_aa"], int)
        assert isinstance(data["end_aa"], int)

@pytest.mark.integration
def test_hunter_analyst_hunt_failure_not_found():
    """
    Tests the /hunt endpoint with a gene symbol that is unlikely to be found.
    This validates the service's ability to handle bad input gracefully.
    """
    target_gene = "FAKEGENE123XYZ"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # The service should return a 404 Not Found
        assert response.status_code == 404, f"Expected status 404, but got {response.status_code}: {response.text}"
        data = response.json()
        assert "detail" in data
        assert "Could not retrieve DNA sequence" in data["detail"] 
import httpx
import os

# The URL for our deployed service
HUNTER_ANALYST_URL = "https://crispro--hunter-analyst-service-hunteranalyst-hunt.modal.run"

@pytest.mark.integration
def test_hunter_analyst_hunt_success():
    """
    Tests the /hunt endpoint of the Hunter-Analyst service with a valid gene symbol.
    This validates the "happy path" for our reconnaissance service.
    """
    target_gene = "MMP9"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # 1. Assert successful response
        assert response.status_code == 200, f"Expected status 200, but got {response.status_code}: {response.text}"
        
        # 2. Assert the response payload structure
        data = response.json()
        expected_keys = ["gene_symbol", "domain_name", "target_motif_dna", "start_aa", "end_aa"]
        for key in expected_keys:
            assert key in data, f"Response missing expected key: {key}"
        
        # 3. Assert data integrity
        assert data["gene_symbol"] == target_gene
        assert isinstance(data["target_motif_dna"], str)
        assert len(data["target_motif_dna"]) > 0
        assert isinstance(data["start_aa"], int)
        assert isinstance(data["end_aa"], int)

@pytest.mark.integration
def test_hunter_analyst_hunt_failure_not_found():
    """
    Tests the /hunt endpoint with a gene symbol that is unlikely to be found.
    This validates the service's ability to handle bad input gracefully.
    """
    target_gene = "FAKEGENE123XYZ"
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(HUNTER_ANALYST_URL, json={"gene_symbol": target_gene})
        
        # The service should return a 404 Not Found
        assert response.status_code == 404, f"Expected status 404, but got {response.status_code}: {response.text}"
        data = response.json()
        assert "detail" in data
        assert "Could not retrieve DNA sequence" in data["detail"] 