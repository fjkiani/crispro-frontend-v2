import pytest
import httpx

UNIFIED_ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

@pytest.mark.integration
def test_unified_oracle_generate_action():
    """
    Tests the 'generate' action of the Unified Oracle.
    Validates that our foundational weapon can forge new sequences.
    """
    payload = {
        "action": "generate",
        "params": {
            "prompt": "ATG",
            "gen_params": {"n_tokens": 10}
        }
    }
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "completion" in data
        assert isinstance(data["completion"], str)
        assert len(data["completion"]) > 3

@pytest.mark.integration
def test_unified_oracle_score_action():
    """
    Tests the 'score' action of the Unified Oracle.
    Validates that our foundational weapon can judge the impact of a variant.
    """
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": "ATGC",
            "alternate_sequence": "ATGG"
        }
    }

    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "delta_score" in data
        assert isinstance(data["delta_score"], float)
        assert "interpretation" in data 
import httpx

UNIFIED_ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

@pytest.mark.integration
def test_unified_oracle_generate_action():
    """
    Tests the 'generate' action of the Unified Oracle.
    Validates that our foundational weapon can forge new sequences.
    """
    payload = {
        "action": "generate",
        "params": {
            "prompt": "ATG",
            "gen_params": {"n_tokens": 10}
        }
    }
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "completion" in data
        assert isinstance(data["completion"], str)
        assert len(data["completion"]) > 3

@pytest.mark.integration
def test_unified_oracle_score_action():
    """
    Tests the 'score' action of the Unified Oracle.
    Validates that our foundational weapon can judge the impact of a variant.
    """
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": "ATGC",
            "alternate_sequence": "ATGG"
        }
    }

    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "delta_score" in data
        assert isinstance(data["delta_score"], float)
        assert "interpretation" in data 
import httpx

UNIFIED_ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

@pytest.mark.integration
def test_unified_oracle_generate_action():
    """
    Tests the 'generate' action of the Unified Oracle.
    Validates that our foundational weapon can forge new sequences.
    """
    payload = {
        "action": "generate",
        "params": {
            "prompt": "ATG",
            "gen_params": {"n_tokens": 10}
        }
    }
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "completion" in data
        assert isinstance(data["completion"], str)
        assert len(data["completion"]) > 3

@pytest.mark.integration
def test_unified_oracle_score_action():
    """
    Tests the 'score' action of the Unified Oracle.
    Validates that our foundational weapon can judge the impact of a variant.
    """
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": "ATGC",
            "alternate_sequence": "ATGG"
        }
    }

    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "delta_score" in data
        assert isinstance(data["delta_score"], float)
        assert "interpretation" in data 
import httpx

UNIFIED_ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

@pytest.mark.integration
def test_unified_oracle_generate_action():
    """
    Tests the 'generate' action of the Unified Oracle.
    Validates that our foundational weapon can forge new sequences.
    """
    payload = {
        "action": "generate",
        "params": {
            "prompt": "ATG",
            "gen_params": {"n_tokens": 10}
        }
    }
    
    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "completion" in data
        assert isinstance(data["completion"], str)
        assert len(data["completion"]) > 3

@pytest.mark.integration
def test_unified_oracle_score_action():
    """
    Tests the 'score' action of the Unified Oracle.
    Validates that our foundational weapon can judge the impact of a variant.
    """
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": "ATGC",
            "alternate_sequence": "ATGG"
        }
    }

    with httpx.Client(timeout=60.0) as client:
        response = client.post(UNIFIED_ORACLE_URL, json=payload)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "delta_score" in data
        assert isinstance(data["delta_score"], float)
        assert "interpretation" in data 