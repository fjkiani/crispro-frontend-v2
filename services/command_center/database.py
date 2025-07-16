from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
import os
import asyncio
from typing import Callable, Any

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
    
    engine = create_engine(
        db_url,
        connect_args={"check_same_thread": False}, # Needed for SQLite
        echo=False # Set to True to log SQL queries
    )
    
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    
    # Create all tables in the database
    Base.metadata.create_all(bind=engine)

# --- Utility Functions ---
async def run_in_threadpool(sync_func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """
    Runs a synchronous function in a separate thread to avoid blocking the
    asyncio event loop.
    """
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(None, lambda: sync_func(*args, **kwargs)) 