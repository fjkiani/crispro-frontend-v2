from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base, Session
import os
import asyncio
from typing import Callable, Any
from . import models
from . import schemas
from .models import BattlePlan, Patient, Mutation


# --- Database Configuration ---
# This path is relative to the container's root when running in Modal.
# For local testing, this will be overridden.
db_path = "/root/data/databases"
# This check is for the Modal container environment, not local.
if "MODAL_ENVIRONMENT" in os.environ:
    os.makedirs(db_path, exist_ok=True)
SQLALCHEMY_DATABASE_URL = f"sqlite:///{db_path}/threat_matrix.db"

# Engine is no longer a global singleton, but created by a function
# This allows tests to inject a different database URL.
engine = None
SessionLocal = None
Base = declarative_base()

# Dependency to get a DB session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

def initialize_database(db_url: str = SQLALCHEMY_DATABASE_URL):
    """Initializes the database engine and session."""
    global engine, SessionLocal
    
    # Import Base here to avoid circular dependencies at the module level.
    from .models import Base
    
    engine = create_engine(
        db_url,
        connect_args={"check_same_thread": False}, # Needed for SQLite
        echo=False # Set to True to log SQL queries
    )
    
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    
    # Create all tables in the database using the correct Base from models
    Base.metadata.create_all(bind=engine)

# --- Data Access Functions ---

def create_full_battle_plan(data: schemas.PatientDataPacket, db: Session):
    """
    Handles the full logic of finding/creating patient, mutation, and battle plan.
    """
    # 1. Find or create the Patient
    patient = db.query(Patient).filter(Patient.patient_identifier == data.patient_identifier).first()
    if not patient:
        patient = Patient(patient_identifier=data.patient_identifier)
        db.add(patient)
        db.commit()
        db.refresh(patient)

    # 2. Find or create the Mutation for this Patient
    mutation = db.query(Mutation).filter(Mutation.hgvs_p == data.mutation_hgvs_p).first()
    if not mutation:
        mutation = Mutation(patient_id=patient.id, gene=data.gene, hgvs_p=data.mutation_hgvs_p)
        db.add(mutation)
        db.commit()
        db.refresh(mutation)

    # 3. Create the Battle Plan and link it to the mutation
    db_battle_plan = BattlePlan(
        mutation_id=mutation.id,
        status="pending_design",
        target_sequence=data.sequence_for_perplexity
    )
    db.add(db_battle_plan)
    db.commit()
    db.refresh(db_battle_plan)
    return db_battle_plan

# --- Utility Functions ---
async def run_in_threadpool(sync_func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """
    Runs a synchronous function in a separate thread to avoid blocking the
    asyncio event loop.
    """
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(None, lambda: sync_func(*args, **kwargs)) 