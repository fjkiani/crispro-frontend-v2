import pytest
import sys
import os
import tempfile
import unittest

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from fastapi.testclient import TestClient
from fastapi import Depends
from sqlalchemy.orm import Session

from services.command_center import models
from services.command_center.main import fastapi_app, get_db
# Import the new initializer
from services.command_center.database import initialize_database, Base

from unittest.mock import patch

# --- Test Database Setup ---
# Use a temporary file-based SQLite database instead of in-memory
# This ensures both test fixtures and FastAPI app access the same database
def get_test_db_url():
    """Create a temporary SQLite database file for testing."""
    db_fd, db_path = tempfile.mkstemp(suffix='.db')
    os.close(db_fd)  # Close the file descriptor, we only need the path
    return f"sqlite:///{db_path}", db_path

# Global variables for test database
test_db_url, test_db_path = get_test_db_url()

# Initialize the database with the test URL *before* creating the engine
initialize_database(test_db_url)

# Now create the engine and session from the initialized test config
from services.command_center.database import engine as test_engine
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)


# --- Pytest Fixtures ---
@pytest.fixture(scope="function")
def dbsession():
    """
    Creates a clean database session for a test.
    This is the master session that will be shared.
    """
    Base.metadata.create_all(bind=test_engine)
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()
        Base.metadata.drop_all(bind=test_engine)

@pytest.fixture(scope="function")
def client(dbsession):
    """
    Provides a test client that correctly initializes the CommandCenter
    and its database dependencies.
    """
    # Instead of trying to instantiate the Modal class (which fails outside Modal runtime),
    # we manually perform the same setup that the Modal class would do.
    
    # The database is already initialized with the test URL globally.
    # We just need to ensure the app uses the test session.
    def override_get_db_for_test():
        yield dbsession
    
    fastapi_app.dependency_overrides[get_db] = override_get_db_for_test
    
    yield TestClient(fastapi_app)
    
    # Clean up the override
    fastapi_app.dependency_overrides.clear()


# Cleanup the temporary database file after all tests
def pytest_sessionfinish(session, exitstatus):
    """Clean up the temporary database file after all tests."""
    if os.path.exists(test_db_path):
        os.unlink(test_db_path)

# --- Tests ---
def test_database_connection_and_model_creation(dbsession):
    """
    Tests that the tables are created in the test database.
    This test now uses the dbsession fixture directly.
    """
    with test_engine.connect() as connection:
        inspector = connection.connection.cursor().connection.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = inspector.fetchall()
        table_names = [table[0] for table in tables]
        assert "patients" in table_names
        assert "battle_plans" in table_names

def test_create_patient(dbsession):
    """Tests creating a patient directly in the database."""
    patient_identifier = "test_patient_01"
    patient = models.Patient(patient_identifier=patient_identifier)
    dbsession.add(patient)
    dbsession.commit()
    retrieved_patient = dbsession.query(models.Patient).filter(models.Patient.patient_identifier == patient_identifier).first()
    assert retrieved_patient is not None
    assert retrieved_patient.patient_identifier == patient_identifier

def test_battle_plan_can_store_interventions(dbsession):
    """Tests that the proposed_interventions JSON field can be written to and read from."""
    # 1. Create prerequisite patient and mutation
    patient = models.Patient(patient_identifier="intervention_patient_01")
    dbsession.add(patient)
    dbsession.commit()

    mutation = models.Mutation(
        patient_id=patient.id,
        gene="BRAF",
        hgvs_c="c.1799T>A",
        hgvs_p="p.Val600Glu"
    )
    dbsession.add(mutation)
    dbsession.commit()

    # 2. Create BattlePlan with intervention data
    guide_data = {"guides": ["ACGTACGTACGTACGTACGT", "TCGATCGATCGATCGATCGA"]}
    battle_plan = models.BattlePlan(
        mutation_id=mutation.id,
        proposed_interventions=guide_data
    )
    dbsession.add(battle_plan)
    dbsession.commit()
    dbsession.refresh(battle_plan)

    # 3. Retrieve and verify
    retrieved_plan = dbsession.query(models.BattlePlan).filter(models.BattlePlan.id == battle_plan.id).first()
    assert retrieved_plan is not None
    assert retrieved_plan.proposed_interventions == guide_data
    assert retrieved_plan.proposed_interventions["guides"][0] == "ACGTACGTACGTACGTACGT"


@pytest.mark.anyio(backend='asyncio')
@patch('tools.sv_scanner.scan_for_svs')
@patch('tools.alphamissense_client.AlphaMissenseClient.get_score')
@patch('tools.esm_client.ESMClient.get_log_likelihood_ratio')
@patch('httpx.AsyncClient.post')
async def test_formulate_battle_plan_endpoint(
    mock_httpx_post, mock_esm, mock_am, mock_sv, client, dbsession
):
    """
    This test now creates the data directly to be more robust and focuses
    on the endpoint logic itself.
    """
    # This test is now massively simplified because the client fixture handles setup.
    # We just need to set up mocks and make the request.
    from services.command_center.main import BattlePlanRequest
    
    # Setup mocks
    mock_sv.return_value = {"sv_type": "DEL", "length": 500}
    mock_am.return_value = 0.95
    mock_esm.return_value = -12.34

    async def get_mock_response(*args, **kwargs):
        # A simple mock that returns a successful response for the perplexity score
        class MockResponse:
            def raise_for_status(self):
                pass
            def json(self):
                return {"mean_perplexity": 0.8}
        return MockResponse()

    mock_httpx_post.side_effect = get_mock_response
    
    # Make the request
    response = client.post(
        "/v2/workflow/formulate_battle_plan",
        json={
            "patient_identifier": "test_patient_02",
            "mutation_hgvs_c": "c.1799T>A",
            "mutation_hgvs_p": "p.Val600Glu",
            "gene": "BRAF",
            "bam_file_path": "/fake/path.bam",
            "sequence_for_perplexity": "ATGC...",
            "protein_sequence": "MVL...",
            "transcript_id": "NM_004333"
        }
    )
    
    assert response.status_code == 200
    data = response.json()
    assert "plan_id" in data
    assert data["message"] == "Battle plan formulated successfully."

    # Verify data in the database
    plan_id = data["plan_id"]
    battle_plan = dbsession.query(models.BattlePlan).filter(models.BattlePlan.id == plan_id).first()
    assert battle_plan is not None
    assert battle_plan.perplexity_score == 0.8
    assert battle_plan.ensemble_scores["alphamissense"] == 0.95


@pytest.mark.anyio(backend='asyncio')
@patch('httpx.AsyncClient.post')
async def test_design_guide_rna_endpoint(mock_httpx_post, client, dbsession):
    """
    Tests the guide RNA design endpoint. This also becomes much simpler.
    """
    from services.command_center.main import GuideDesignRequest
    
    # 1. Create a battle plan to design guides for
    patient = models.Patient(patient_identifier="guide_patient_01")
    dbsession.add(patient); dbsession.commit()
    mutation = models.Mutation(patient_id=patient.id, gene="BRAF", hgvs_c="c.1799T>A", hgvs_p="p.Val600Glu")
    dbsession.add(mutation); dbsession.commit()
    battle_plan = models.BattlePlan(mutation_id=mutation.id)
    dbsession.add(battle_plan); dbsession.commit()

    # 2. Setup mock for the evo_service call
    async def get_mock_response(*args, **kwargs):
        class MockResponse:
            def raise_for_status(self):
                pass
            def json(self):
                return {"generated_guides": ["guide1", "guide2"]}
        return MockResponse()

    mock_httpx_post.side_effect = get_mock_response

    # 3. Make the request
    response = client.post(
        "/v2/design/guide_rna",
        json={
            "battle_plan_id": battle_plan.id,
            "target_sequence": "ATGC...",
            "num_guides": 2
        }
    )

    assert response.status_code == 200
    data = response.json()
    assert data["generated_guides"] == ["guide1", "guide2"]

    # 4. Verify the guides were stored
    dbsession.refresh(battle_plan)
    assert battle_plan.proposed_interventions["guides"] == ["guide1", "guide2"] 