import modal
import os
import time
import logging
import asyncio
import uuid
import httpx

from fastapi import FastAPI, BackgroundTasks, HTTPException, Depends
from pydantic import BaseModel
from sqlalchemy.orm import Session
from contextlib import asynccontextmanager

# Local application imports
from services.command_center import database
from services.command_center.models import BattlePlan as BattlePlanModel, Patient, Mutation
from services.command_center.schemas import PatientDataPacket, BattlePlanResponse
# Import our new, dedicated client
from tools.command_center_client import CommandCenterClient

# --- App Setup ---
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

fastapi_app = FastAPI(
    title="CRISPR Assistant - AI General Command Center",
    version="3.0",
    description="Central nervous system for orchestrating genomic threat analysis and countermeasure design (V3 Async Workflow)."
)

# --- Modal App Definition ---
# NOTE: The app name here dictates the final deployment URL
app = modal.App("crispr-assistant-command-center-v3") 

# FORCE REBUILD v3.1: Added a comment to invalidate the cache and force a rebuild
# with the latest version of the tools directory.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "fastapi", "uvicorn", "sqlalchemy", "pydantic", "requests", 
        "biopython", "pandas", "httpx", "psycopg2-binary", "cowsay"
    )
    .add_local_dir("services/command_center", remote_path="/root/services/command_center")
    # Add the tools directory so the new client can be imported
    .add_local_dir("tools", remote_path="/root/tools")
)

# --- API Endpoints ---

@fastapi_app.post("/v3/workflow/formulate_battle_plan", status_code=201, response_model=BattlePlanResponse)
async def formulate_battle_plan_v3(data: PatientDataPacket, db: Session = Depends(database.get_db)):
    logger.info(f"V3: Received battle plan formulation for patient {data.patient_identifier}")

    # 1. Find or create the Patient
    patient = db.query(Patient).filter(Patient.patient_identifier == data.patient_identifier).first()
    if not patient:
        patient = Patient(patient_identifier=data.patient_identifier)
        db.add(patient); db.commit(); db.refresh(patient)

    # 2. Find or create the Mutation for this Patient
    mutation = db.query(Mutation).filter(Mutation.hgvs_p == data.mutation_hgvs_p).first()
    if not mutation:
        mutation = Mutation(patient_id=patient.id, gene=data.gene, hgvs_p=data.mutation_hgvs_p)
        db.add(mutation); db.commit(); db.refresh(mutation)

    # 3. Create the Battle Plan and link it to the mutation
    db_battle_plan = BattlePlanModel(
        mutation_id=mutation.id,
        status="pending_design",
        target_sequence=data.sequence_for_perplexity
    )
    db.add(db_battle_plan); db.commit(); db.refresh(db_battle_plan)

    logger.info(f"V3: Battle Plan {db_battle_plan.id} for Mutation {mutation.id} created.")
    return BattlePlanResponse(plan_id=str(db_battle_plan.id), message="Battle plan formulated successfully.")

async def run_guide_design_campaign_background(battle_plan_id: int, db: Session):
    logger.info(f"V3 BACKGROUND TASK: Started for Battle Plan ID: {battle_plan_id}")
    
    battle_plan = db.query(BattlePlanModel).filter(BattlePlanModel.id == battle_plan_id).first()
    if not battle_plan:
        logger.error(f"[BP {battle_plan_id}] Battle Plan not found.")
        db.close(); return

    battle_plan.status = "design_in_progress"; db.commit()

    client = CommandCenterClient()
    
    try:
        # Step 1: Forge Candidates with ZetaForge (Live)
        guides = await client.forge_guides(battle_plan.target_sequence)
        if "error" in guides:
            raise Exception(f"ZetaForge failed: {guides['details']}")
        
        if not guides:
            raise Exception("ZetaForge returned no guide candidates.")

        # Step 2: Score and Validate (using a mix of live and mock clients)
        # Note: We run these in parallel for efficiency
        logger.info(f"[BP {battle_plan_id}] Running parallel validation for {len(guides)} candidates...")
        oracle_task = client.score_guides_with_oracle(guides, battle_plan.target_sequence)
        safety_task = client.check_off_targets_with_blast(guides)
        immuno_task = client.check_immunogenicity(guides)

        scored_guides, safety_data, immuno_data = await asyncio.gather(
            oracle_task, safety_task, immuno_task
        )
        
        # Step 3: Consolidate and Rank (mock ranking for now)
        logger.info(f"[BP {battle_plan_id}] Consolidating results into guide dossiers...")
        guide_dossiers = []
        # Simple consolidation by matching guide sequence
        for guide in guides:
            dossier = {"guide_sequence": guide}
            # Find corresponding data from other services
            dossier.update(next((item for item in scored_guides if item["guide_sequence"] == guide), {}))
            dossier.update(next((item for item in safety_data if item["guide_sequence"] == guide), {}))
            dossier.update(next((item for item in immuno_data if item["guide_sequence"] == guide), {}))
            
            # Mock "Assassin Score" calculation
            zeta_score = dossier.get("zeta_score", 0)
            off_targets = dossier.get("off_target_count", 1) # avoid division by zero
            immuno_score = dossier.get("immunogenicity_score", 1)
            dossier["assassin_score"] = (zeta_score * 0.6) + ((1 - off_targets) * 0.3) + ((1 - immuno_score) * 0.1)
            
            guide_dossiers.append(dossier)

        # Sort by the mock score
        guide_dossiers.sort(key=lambda x: x["assassin_score"], reverse=True)

        battle_plan.results = {"guides": guide_dossiers}
        battle_plan.status = "interventions_designed"
        logger.info(f"✅ Stored {len(guide_dossiers)} designs for Battle Plan {battle_plan_id}")

    except Exception as e:
        logger.error(f"💥 V3 guide design FAILED for BP {battle_plan_id}: {e}", exc_info=True)
        battle_plan.status = "design_failed"
        battle_plan.results = {"error": str(e)}
    
    finally:
        db.commit(); db.close()

@fastapi_app.post("/v3/workflow/execute_guide_design_campaign/{battle_plan_id}", status_code=202)
async def execute_guide_design_campaign(battle_plan_id: int, background_tasks: BackgroundTasks, db: Session = Depends(database.get_db)):
    battle_plan = db.query(BattlePlanModel).filter(BattlePlanModel.id == battle_plan_id).first()
    if not (battle_plan and battle_plan.target_sequence):
        raise HTTPException(status_code=404, detail="Battle Plan or target sequence not found")

    db_session_for_task = database.SessionLocal()
    background_tasks.add_task(run_guide_design_campaign_background, battle_plan_id, db_session_for_task)
    return {"message": "Guide design campaign started.", "battle_plan_id": battle_plan_id}

@fastapi_app.get("/v3/battle_plan/{battle_plan_id}")
async def get_battle_plan_status_and_results(battle_plan_id: int, db: Session = Depends(database.get_db)):
    battle_plan = db.query(BattlePlanModel).filter(BattlePlanModel.id == battle_plan_id).first()
    if not battle_plan: raise HTTPException(status_code=404, detail="Battle Plan not found")
    return {"battle_plan_id": battle_plan.id, "status": battle_plan.status, "results": battle_plan.results}

# --- Modal Class Definition ---
@app.cls(image=image, timeout=2000)
class CommandCenter:
    @modal.enter()
    def setup(self):
        logger.info("--- COMMAND CENTER DEPLOYMENT V3.2 (LIVE BLAST) ---")
        database.initialize_database()

    @modal.asgi_app()
    def api(self): return fastapi_app

@fastapi_app.get("/")
def read_root(): return {"message": "Command Center V3 is operational."} 