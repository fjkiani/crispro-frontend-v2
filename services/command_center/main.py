import modal
import os
import requests
import re
import xml.etree.ElementTree as ET
from fastapi import FastAPI, BackgroundTasks, HTTPException, Depends
from pydantic import BaseModel
import time
import random
import pandas as pd
import sqlite3
from typing import Optional
from fastapi.responses import JSONResponse
from Bio import SeqIO
import logging
import httpx
from modal import App, asgi_app, Image
from sqlalchemy.orm import Session
from contextlib import asynccontextmanager

# The sys.path hack is removed. Imports will now be absolute.

from services.command_center.database import SessionLocal, engine, get_db
from services.command_center.models import Base, BattlePlan as BattlePlanModel
from services.command_center.schemas import BattlePlan, PatientDataPacket, BattlePlanResponse, GuideDesignResponse


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Lifespan Manager & DB Setup ---
# NOTE: The imports for models and db setup are now lazy-loaded inside the class
# to avoid top-level import errors in Modal.

@asynccontextmanager
async def lifespan(app: FastAPI):
    # This context manager can be used for other startup/shutdown logic if needed,
    # but the primary db setup is now handled in the class @modal.enter()
    yield

fastapi_app = FastAPI(
    title="CRISPR Assistant - AI General Command Center",
    version="2.0",
    description="The central nervous system for orchestrating genomic threat analysis and countermeasure design.",
    lifespan=lifespan
)

# Pydantic models remain the same...
class BattlePlanRequest(BaseModel):
    patient_identifier: str
    mutation_hgvs_c: str
    mutation_hgvs_p: str
    gene: str
    bam_file_path: str 
    sequence_for_perplexity: str 
    protein_sequence: str
    transcript_id: str

class BattlePlanResponse(BaseModel):
    plan_id: int
    message: str

class GuideDesignRequest(BaseModel):
    battle_plan_id: int
    target_sequence: str
    num_guides: int = 5

class GuideDesignResponse(BaseModel):
    battle_plan_id: int
    message: str
    generated_guides: list

# --- API Endpoints ---
# Dependency placeholder, will be properly defined in the class
def get_db():
    pass 

# The endpoint implementations are now defined at the top level,
# not inside the class, to avoid the 'self' parameter issue.

ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-api.modal.run/invoke"


@fastapi_app.post("/v2/workflow/formulate_battle_plan", status_code=201, response_model=BattlePlanResponse)
async def formulate_battle_plan_v2(
    data: PatientDataPacket, db: Session = Depends(get_db)
):
    """
    Asynchronous endpoint to formulate a new battle plan.
    This version integrates with the Zeta Oracle for live variant scoring.
    """
    logger.info(f"Received battle plan formulation request for patient {data.patient_identifier}")

    # Use the existing synchronous logic but run it in a thread to keep the endpoint async
    # This is a common pattern to reuse synchronous code in an async context.
    
    # Critical Change: We are replacing the generic HTTP client with Modal's client
    # to correctly handle the service-to-service communication.
    try:
        # Step 1: Use HTTP client instead of Modal function lookup
        oracle_url = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run/invoke"
        
        # Step 2: Prepare the payload for the oracle
        ref_sequence = data.sequence_for_perplexity # Using this field for reference seq
        alt_sequence = data.protein_sequence      # Using this field for alternate seq
        
        oracle_payload = {
            "action": "score",
            "params": {
                "reference_sequence": ref_sequence,
                "alternate_sequence": alt_sequence,
            },
        }

        # Step 3: Make HTTP request to the oracle
        logger.info("🧬 Calling Zeta Oracle for variant impact score via HTTP...")
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.post(oracle_url, json=oracle_payload)
            response.raise_for_status()
            oracle_results = response.json()
        
        logger.info(f"✅ Zeta Oracle responded: {oracle_results}")

        zeta_score = oracle_results.get("zeta_score", -999.0)
        ensemble_scores = {"zeta_score": zeta_score} # Store the score

    except Exception as e:
        logger.error(f"💥 Failed to get variant score from Zeta Oracle: {e}")
        # Fail gracefully if the oracle is down
        zeta_score = -999.0
        ensemble_scores = {"zeta_score": zeta_score, "error": str(e)}


    # Create a new battle plan entry in the database
    db_battle_plan = BattlePlanModel(
        patient_identifier=data.patient_identifier,
        gene=data.gene,
        mutation_hgvs_p=data.mutation_hgvs_p,
        status="pending_review",
        # We now store the live score from the oracle
        ensemble_scores=ensemble_scores,
        perplexity_score=zeta_score, # Using perplexity field to store the main score for now
        # ... other fields
        structural_variants={}, 
        evidence_features={},
        proposed_interventions={},
    )
    db.add(db_battle_plan)
    db.commit()
    db.refresh(db_battle_plan)

    logger.info(f"Battle Plan {db_battle_plan.id} created and stored.")
    return BattlePlanResponse(plan_id=str(db_battle_plan.id), message="Battle plan formulated successfully.")

@fastapi_app.post("/v2/design/guide_rna", response_model=GuideDesignResponse)
async def design_guide_rna(request: GuideDesignRequest, db: Session = Depends(get_db)):
    """
    Designs guide RNAs for a specific target sequence related to a battle plan.
    """
    logger.info(f"Received guide RNA design request for Battle Plan ID: {request.battle_plan_id}")
    from services.command_center import models

    # 1. Find the battle plan
    battle_plan = db.query(models.BattlePlan).filter(models.BattlePlan.id == request.battle_plan_id).first()
    if not battle_plan:
        raise HTTPException(status_code=404, detail="Battle Plan not found")

    # 2. Call the Evo2 service to generate guides
    logger.info("Requesting guide RNA generation from Evo 2 service...")
    evo_service_url = "https://crispro--evo2-perplexity-service-evo2service-web-app.modal.run/generate_guide_rna"
    generated_guides = []
    
    async with httpx.AsyncClient(follow_redirects=True) as client:
        try:
            response = await client.post(
                evo_service_url,
                json={"target_sequence": request.target_sequence, "num_guides": request.num_guides},
                timeout=300.0
            )
            response.raise_for_status()
            generated_guides = response.json().get("generated_guides", [])
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error calling Evo2 service for guide generation: {e.response.text}")
            raise HTTPException(status_code=500, detail="Failed to generate guide RNAs from Evo2 service")
        except httpx.RequestError as e:
            logger.error(f"Request error calling Evo2 service for guide generation: {e}")
            raise HTTPException(status_code=500, detail="Failed to connect to Evo2 service")

    # 3. Store the generated guides in the battle plan
    if generated_guides:
        battle_plan.proposed_interventions = {"guides": generated_guides}
        db.commit()
        db.refresh(battle_plan)
        logger.info(f"Stored {len(generated_guides)} generated guides in Battle Plan {request.battle_plan_id}")
    
    return GuideDesignResponse(
        battle_plan_id=request.battle_plan_id,
        message="Guide RNA design process completed.",
        generated_guides=generated_guides
    )

# --- Modal App Definition ---
app = App("crispr-assistant-command-center-v2")

image = Image.debian_slim(python_version="3.11").pip_install(
    "fastapi", "uvicorn", "sqlalchemy", "pydantic", "requests", "biopython", "pysam", "pandas", "httpx"
).add_local_dir("tools", remote_path="/root/tools").add_local_dir("services", remote_path="/root/services").add_local_dir("data", remote_path="/root/data")

@app.cls(image=image, timeout=1200)
class CommandCenter:
    @modal.enter()
    def setup(self):
        """This runs once when the container starts."""
        import sys
        sys.path.append('/root')
        
        # Lazy import and initialize database components
        global get_db_real, models
        from services.command_center import models
        from services.command_center import database
        
        # CRITICAL: Initialize the database first
        database.initialize_database()
        
        # Now access the initialized engine and SessionLocal from the module
        models.Base.metadata.create_all(bind=database.engine)
        
        def get_db_real():
            db = database.SessionLocal()
            try:
                yield db
            finally:
                db.close()
        
        # Override the placeholder dependency
        fastapi_app.dependency_overrides[get_db] = get_db_real
        logger.info("Database tables created and dependency overridden.")

    @modal.asgi_app()
    def api(self):
        return fastapi_app 