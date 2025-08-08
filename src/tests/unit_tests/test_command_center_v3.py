import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from unittest.mock import patch, MagicMock, AsyncMock

# Adjust the path to import from the service
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from services.command_center.main import fastapi_app, get_db
from services.command_center.models import Base, BattlePlan as BattlePlanModel
from services.command_center.schemas import PatientDataPacket

# --- Test Database Setup ---
SQLALCHEMY_DATABASE_URL = "sqlite:///./test.db"
engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base.metadata.create_all(bind=engine)

def override_get_db():
    try:
        db = TestingSessionLocal()
        yield db
    finally:
        db.close()

fastapi_app.dependency_overrides[get_db] = override_get_db

client = TestClient(fastapi_app)

# --- Fixtures ---
@pytest.fixture(scope="function")
def db_session():
    """Create a new database session for a test."""
    Base.metadata.drop_all(bind=engine)
    Base.metadata.create_all(bind=engine)
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()

# --- Unit Tests ---
def test_formulate_battle_plan_saves_sequence(db_session):
    """
    Tests if the /v2/formulate_battle_plan endpoint correctly saves the
    target_sequence to the database.
    """
    # Mock the call to Zeta Oracle
    with patch('httpx.AsyncClient.post', new_callable=AsyncMock) as mock_post:
        mock_post.return_value.status_code = 200
        mock_post.return_value.json.return_value = {"zeta_score": -2.5}

        patient_data = {
            "patient_identifier": "TEST_PATIENT_01",
            "gene": "RUNX1",
            "mutation_hgvs_c": "c.123A>G",
            "mutation_hgvs_p": "p.Lys41Arg",
            "transcript_id": "ENSTFAKE123",
            "bam_file_path": "/fake/path.bam",
            "protein_sequence": "FAKEPROTEIN",
            "sequence_for_perplexity": "AGCTAGCTAGCT", # This is the sequence we want to save
        }

        response = client.post("/v2/workflow/formulate_battle_plan", json=patient_data)
        
        assert response.status_code == 201
        response_data = response.json()
        assert "plan_id" in response_data

        # Verify the data was saved correctly in the database
        plan_id = int(response_data["plan_id"])
        db_plan = db_session.query(BattlePlanModel).filter(BattlePlanModel.id == plan_id).first()
        
        assert db_plan is not None
        assert db_plan.target_sequence == "AGCTAGCTAGCT"

@patch('services.command_center.main.run_guide_design_campaign_background', new_callable=MagicMock)
def test_execute_guide_design_campaign_starts_background_task(mock_background_task, db_session):
    """
    Tests if the /v3/execute_guide_design_campaign endpoint correctly
    kicks off the background task.
    """
    # First, create a battle plan to work with
    test_plan = BattlePlanModel(
        patient_identifier="TEST_PATIENT_02",
        gene="RUNX1",
        mutation_hgvs_p="p.Arg135fs",
        status="pending_review",
        target_sequence="GATTACAGATTACA"
    )
    db_session.add(test_plan)
    db_session.commit()
    db_session.refresh(test_plan)

    # Call the endpoint to trigger the background task
    response = client.post(f"/v3/workflow/execute_guide_design_campaign/{test_plan.id}")

    assert response.status_code == 202
    assert response.json()["message"] == "Guide design campaign started."
    
    # Check that the battle plan status was updated
    db_plan = db_session.query(BattlePlanModel).filter(BattlePlanModel.id == test_plan.id).first()
    assert db_plan.status == "design_in_progress"

    # We cannot directly test if the background task ran, as it's out of band.
    # The `patch` on the function itself is the standard way to test this.
    # Here we can at least confirm our endpoint behaves as expected.

def test_get_battle_plan_status_and_results(db_session):
    """
    Tests the polling endpoint /v3/battle_plan/{battle_plan_id}.
    """
    # Create a completed battle plan in the DB
    results_payload = {"guides": [{"guide_sequence": "ACGT", "assassin_score": 0.95}]}
    test_plan = BattlePlanModel(
        patient_identifier="TEST_PATIENT_03",
        gene="RUNX1",
        mutation_hgvs_p="p.Gly123Asp",
        status="interventions_designed",
        target_sequence="GATTACAGATTACA",
        proposed_interventions=results_payload
    )
    db_session.add(test_plan)
    db_session.commit()
    db_session.refresh(test_plan)

    # Poll the endpoint
    response = client.get(f"/v3/battle_plan/{test_plan.id}")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "interventions_designed"
    assert data["results"] == results_payload 
 
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from unittest.mock import patch, MagicMock, AsyncMock

# Adjust the path to import from the service
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from services.command_center.main import fastapi_app, get_db
from services.command_center.models import Base, BattlePlan as BattlePlanModel
from services.command_center.schemas import PatientDataPacket

# --- Test Database Setup ---
SQLALCHEMY_DATABASE_URL = "sqlite:///./test.db"
engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base.metadata.create_all(bind=engine)

def override_get_db():
    try:
        db = TestingSessionLocal()
        yield db
    finally:
        db.close()

fastapi_app.dependency_overrides[get_db] = override_get_db

client = TestClient(fastapi_app)

# --- Fixtures ---
@pytest.fixture(scope="function")
def db_session():
    """Create a new database session for a test."""
    Base.metadata.drop_all(bind=engine)
    Base.metadata.create_all(bind=engine)
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()

# --- Unit Tests ---
def test_formulate_battle_plan_saves_sequence(db_session):
    """
    Tests if the /v2/formulate_battle_plan endpoint correctly saves the
    target_sequence to the database.
    """
    # Mock the call to Zeta Oracle
    with patch('httpx.AsyncClient.post', new_callable=AsyncMock) as mock_post:
        mock_post.return_value.status_code = 200
        mock_post.return_value.json.return_value = {"zeta_score": -2.5}

        patient_data = {
            "patient_identifier": "TEST_PATIENT_01",
            "gene": "RUNX1",
            "mutation_hgvs_c": "c.123A>G",
            "mutation_hgvs_p": "p.Lys41Arg",
            "transcript_id": "ENSTFAKE123",
            "bam_file_path": "/fake/path.bam",
            "protein_sequence": "FAKEPROTEIN",
            "sequence_for_perplexity": "AGCTAGCTAGCT", # This is the sequence we want to save
        }

        response = client.post("/v2/workflow/formulate_battle_plan", json=patient_data)
        
        assert response.status_code == 201
        response_data = response.json()
        assert "plan_id" in response_data

        # Verify the data was saved correctly in the database
        plan_id = int(response_data["plan_id"])
        db_plan = db_session.query(BattlePlanModel).filter(BattlePlanModel.id == plan_id).first()
        
        assert db_plan is not None
        assert db_plan.target_sequence == "AGCTAGCTAGCT"

@patch('services.command_center.main.run_guide_design_campaign_background', new_callable=MagicMock)
def test_execute_guide_design_campaign_starts_background_task(mock_background_task, db_session):
    """
    Tests if the /v3/execute_guide_design_campaign endpoint correctly
    kicks off the background task.
    """
    # First, create a battle plan to work with
    test_plan = BattlePlanModel(
        patient_identifier="TEST_PATIENT_02",
        gene="RUNX1",
        mutation_hgvs_p="p.Arg135fs",
        status="pending_review",
        target_sequence="GATTACAGATTACA"
    )
    db_session.add(test_plan)
    db_session.commit()
    db_session.refresh(test_plan)

    # Call the endpoint to trigger the background task
    response = client.post(f"/v3/workflow/execute_guide_design_campaign/{test_plan.id}")

    assert response.status_code == 202
    assert response.json()["message"] == "Guide design campaign started."
    
    # Check that the battle plan status was updated
    db_plan = db_session.query(BattlePlanModel).filter(BattlePlanModel.id == test_plan.id).first()
    assert db_plan.status == "design_in_progress"

    # We cannot directly test if the background task ran, as it's out of band.
    # The `patch` on the function itself is the standard way to test this.
    # Here we can at least confirm our endpoint behaves as expected.

def test_get_battle_plan_status_and_results(db_session):
    """
    Tests the polling endpoint /v3/battle_plan/{battle_plan_id}.
    """
    # Create a completed battle plan in the DB
    results_payload = {"guides": [{"guide_sequence": "ACGT", "assassin_score": 0.95}]}
    test_plan = BattlePlanModel(
        patient_identifier="TEST_PATIENT_03",
        gene="RUNX1",
        mutation_hgvs_p="p.Gly123Asp",
        status="interventions_designed",
        target_sequence="GATTACAGATTACA",
        proposed_interventions=results_payload
    )
    db_session.add(test_plan)
    db_session.commit()
    db_session.refresh(test_plan)

    # Poll the endpoint
    response = client.get(f"/v3/battle_plan/{test_plan.id}")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "interventions_designed"
    assert data["results"] == results_payload 