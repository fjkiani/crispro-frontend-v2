import pytest
from tools.esm_client import ESMClient

@pytest.fixture
def client():
    """Provides an ESMClient instance for testing."""
    return ESMClient()

def test_score_missense_variant_mock(client):
    """
    Tests that the mock client returns the expected score for a missense variant.
    """
    protein_seq = "PROTEINSEQUENCE"
    mutation = "V600E"
    result = client.score_variant(protein_seq, mutation)
    
    assert result is not None
    assert result["mutation"] == mutation
    assert "log_likelihood_score" in result
    assert result["log_likelihood_score"] == -2.5 # Default mock score

def test_score_frameshift_variant_mock(client):
    """
    Tests that the mock client returns a more severe score for a frameshift.
    """
    protein_seq = "PROTEINSEQUENCE"
    mutation = "G646fs"
    result = client.score_variant(protein_seq, mutation)

    assert result is not None
    assert result["mutation"] == mutation
    assert "log_likelihood_score" in result
    assert result["log_likelihood_score"] == -15.0

def test_score_nonsense_variant_mock(client):
    """
    Tests that the mock client returns a more severe score for a nonsense mutation.
    """
    protein_seq = "PROTEINSEQUENCE"
    mutation = "Q50*"
    result = client.score_variant(protein_seq, mutation)

    assert result is not None
    assert result["mutation"] == mutation
    assert "log_likelihood_score" in result
    assert result["log_likelihood_score"] == -15.0

def test_empty_sequence_input(client):
    """
    Tests that the client returns None if the protein sequence is empty.
    """
    result = client.score_variant("", "V600E")
    assert result is None

def test_empty_mutation_input(client):
    """
    Tests that the client returns None if the mutation string is empty.
    """
    result = client.score_variant("PROTEINSEQUENCE", "")
    assert result is None 
 
from tools.esm_client import ESMClient

@pytest.fixture
def client():
    """Provides an ESMClient instance for testing."""
    return ESMClient()

def test_score_missense_variant_mock(client):
    """
    Tests that the mock client returns the expected score for a missense variant.
    """
    protein_seq = "PROTEINSEQUENCE"
    mutation = "V600E"
    result = client.score_variant(protein_seq, mutation)
    
    assert result is not None
    assert result["mutation"] == mutation
    assert "log_likelihood_score" in result
    assert result["log_likelihood_score"] == -2.5 # Default mock score

def test_score_frameshift_variant_mock(client):
    """
    Tests that the mock client returns a more severe score for a frameshift.
    """
    protein_seq = "PROTEINSEQUENCE"
    mutation = "G646fs"
    result = client.score_variant(protein_seq, mutation)

    assert result is not None
    assert result["mutation"] == mutation
    assert "log_likelihood_score" in result
    assert result["log_likelihood_score"] == -15.0

def test_score_nonsense_variant_mock(client):
    """
    Tests that the mock client returns a more severe score for a nonsense mutation.
    """
    protein_seq = "PROTEINSEQUENCE"
    mutation = "Q50*"
    result = client.score_variant(protein_seq, mutation)

    assert result is not None
    assert result["mutation"] == mutation
    assert "log_likelihood_score" in result
    assert result["log_likelihood_score"] == -15.0

def test_empty_sequence_input(client):
    """
    Tests that the client returns None if the protein sequence is empty.
    """
    result = client.score_variant("", "V600E")
    assert result is None

def test_empty_mutation_input(client):
    """
    Tests that the client returns None if the mutation string is empty.
    """
    result = client.score_variant("PROTEINSEQUENCE", "")
    assert result is None 