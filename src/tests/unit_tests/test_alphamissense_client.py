import pytest
from unittest.mock import patch, MagicMock
from tools.alphamissense_client import AlphaMissenseClient
import requests

@pytest.fixture
def client():
    """Provides an AlphaMissenseClient instance for testing."""
    return AlphaMissenseClient(retries=2)

@patch('requests.Session.get')
def test_score_variant_success(mock_get, client):
    """
    Tests a successful API call where a variant is found and scored.
    """
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"pathogenicity_score": 0.95}
    mock_get.return_value = mock_response

    variant_str = "chr7:g.140753336A>T"
    result = client.score_variant(variant_str)

    mock_get.assert_called_once()
    assert result is not None
    assert result["pathogenicity_score"] == 0.95

@patch('requests.Session.get')
def test_score_variant_not_found(mock_get, client):
    """
    Tests the scenario where the API returns a 404, indicating the variant
    was not in the precomputed database. This is a valid, non-error case.
    """
    mock_response = MagicMock()
    mock_response.status_code = 404
    mock_get.return_value = mock_response

    variant_str = "chr1:g.12345A>G"
    result = client.score_variant(variant_str)

    mock_get.assert_called_once()
    assert result is not None
    assert result["error"] == "not_found"
    assert result["pathogenicity_score"] == 0.0

@patch('requests.Session.get')
def test_score_variant_server_error_with_retry(mock_get, client):
    """
    Tests the retry logic when the API returns a server error (e.g., 500)
    on the first attempt and succeeds on the second.
    """
    # Simulate a server error, then a success
    mock_error_response = MagicMock()
    mock_error_response.status_code = 500
    mock_error_response.raise_for_status.side_effect = requests.exceptions.HTTPError

    mock_success_response = MagicMock()
    mock_success_response.status_code = 200
    mock_success_response.json.return_value = {"pathogenicity_score": 0.88}

    mock_get.side_effect = [mock_error_response, mock_success_response]

    variant_str = "chrX:g.100000C>T"
    result = client.score_variant(variant_str)

    assert mock_get.call_count == 2
    assert result is not None
    assert result["pathogenicity_score"] == 0.88

@patch('requests.Session.get')
def test_score_variant_all_retries_fail(mock_get, client):
    """
    Tests the scenario where the API call fails for all retry attempts.
    """
    mock_error_response = MagicMock()
    mock_error_response.status_code = 503
    mock_error_response.raise_for_status.side_effect = requests.exceptions.HTTPError
    mock_get.return_value = mock_error_response

    variant_str = "chr2:g.200000G>A"
    result = client.score_variant(variant_str)

    assert mock_get.call_count == 2  # Initial call + 1 retry = 2 calls for retries=2
    assert result is None

def test_score_variant_empty_input(client):
    """
    Tests that the client handles an empty variant string gracefully.
    """
    result = client.score_variant("")
    assert result is None 
 