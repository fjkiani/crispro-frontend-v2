import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.pool import StaticPool
from unittest.mock import patch, MagicMock, AsyncMock

# Import the single, definitive Base and key modules
from services.command_center import database, models, main
from services.command_center.main import fastapi_app
from services.command_center.models import Base, Patient, Mutation, BattlePlan


@pytest.fixture(scope="function")
def session_override_fixture():
    """
    This is the definitive, all-in-one fixture for integration testing.
    It creates a single, shared in-memory database connection for the entire test
    by using StaticPool, which is essential for testing with SQLite in-memory databases.
    """
    # 1. Create the in-memory database engine with StaticPool
    engine = create_engine(
        "sqlite:///:memory:",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    # 2. Create all tables using the imported Base metadata
    Base.metadata.create_all(bind=engine)

    # 3. Create a sessionmaker bound to this engine
    SessionTesting = sessionmaker(autocommit=False, autoflush=False, bind=engine)

    # 4. Override the app's dependency
    def override_get_db():
        db = SessionTesting()
        try:
            yield db
        finally:
            db.close()

    fastapi_app.dependency_overrides[database.get_db] = override_get_db

    # 5. Yield the client for making requests
    yield TestClient(fastapi_app)

    # 6. Clean up after the test
    fastapi_app.dependency_overrides.clear()
    Base.metadata.drop_all(bind=engine)


def get_session_for_setup(client: TestClient) -> Session:
    """Helper to get a db session for setting up data before a test."""
    override_func = client.app.dependency_overrides[database.get_db]
    return next(override_func())


def test_full_dossier_integration(session_override_fixture: TestClient):
    """
    An end-to-end integration test for the dossier generation workflow.
    """
    db_session = get_session_for_setup(session_override_fixture)

    patient = Patient(patient_identifier="test-patient-001")
    mutation = Mutation(
        gene="BRAF", hgvs_p="p.Val600Glu", locus="chr7:140753336", patient=patient
    )
    battle_plan = BattlePlan(
        mutation=mutation, status="pending", target_sequence="ATGCG..."
    )
    db_session.add(patient)
    db_session.commit()

    mock_cmd_center_instance = MagicMock()
    mock_cmd_center_instance.score_guides_with_oracle = AsyncMock(
        return_value=[
            {"guide_sequence": "ACGT", "zeta_score": 0.95, "confidence": 0.99}
        ]
    )

    # Mock the threat matrix query to return a predictable baseline
    with patch(
        "services.command_center.database.query_threat_matrix",
        return_value={"pathogenic_max": 0.9, "benign_min": 0.1},
    ) as mock_threat_matrix:
        with patch(
            "services.command_center.main.CommandCenter",
            return_value=mock_cmd_center_instance,
        ) as mock_CommandCenter_class:
            response = session_override_fixture.post(
                f"/workflow/generate_intelligence_dossier/{battle_plan.id}"
            )

    assert response.status_code == 200

    db_session_verify = get_session_for_setup(session_override_fixture)
    battle_plan_from_db = db_session_verify.query(BattlePlan).get(battle_plan.id)

    assert battle_plan_from_db.status == "dossier_generated"
    results = battle_plan_from_db.results
    assert results is not None
    assert results["threat_assessment"]["verdict"] == "PATHOGENIC"
    assert "intelligence_dossier" in results


def test_dossier_battle_plan_not_found(session_override_fixture: TestClient):
    """
    Tests the failure case where a non-existent battle_plan_id is provided.
    """
    response = session_override_fixture.post(
        "/workflow/generate_intelligence_dossier/9999"
    )
    assert response.status_code == 404
    # Correct the expected error message to match the application's actual response
    assert response.json() == {"detail": "Battle Plan not found"}
