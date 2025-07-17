# FORCE REBUILD v3.3: Added a comment to invalidate the Modal cache.
import modal
import os
import time
import logging
import asyncio
import uuid
import httpx
from xml.etree import ElementTree as ET


from fastapi import FastAPI, BackgroundTasks, HTTPException, Depends
from pydantic import BaseModel
from sqlalchemy.orm import Session
from contextlib import asynccontextmanager

# Local application imports
from services.command_center import database
from services.command_center.models import BattlePlan as BattlePlanModel, Patient, Mutation
from services.command_center.schemas import PatientDataPacket, BattlePlanResponse
from tools.command_center_client import CommandCenterClient

# --- App Setup ---
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- BEGIN INLINED COMMAND CENTER CLIENT ---

# --- Service URLs ---
# These should match your live Modal deployments.
COMMAND_CENTER_URL = os.environ.get("COMMAND_CENTER_URL", "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run")
ZETA_FORGE_URL = "https://crispro--evo-service-evoservice-api.modal.run" # Corrected from deployment logs
ZETA_ORACLE_URL = "https://crispro--zetascorer-api.modal.run" # Placeholder - Needs to be deployed
BLAST_SERVICE_URL = "https://crispro--blast-service-human-blastservice-web-app.modal.run" # This is the Modal app name for direct calls
IMMUNOGENICITY_SERVICE_URL = "mock" # Placeholder

class CommandCenterClient:
    """
    An asynchronous client to orchestrate guide design campaigns 
    through the CRISPR Assistant's microservices.
    """

    async def _make_request(self, client, method, url, **kwargs):
        """A helper to make async HTTP requests with error handling."""
        try:
            response = await client.request(method, url, **kwargs)
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error for {e.request.url}: {e.response.status_code} - {e.response.text}")
            return {"error": f"HTTP {e.response.status_code}", "details": e.response.text}
        except httpx.RequestError as e:
            logger.error(f"Request error for {e.request.url}: {e}")
            return {"error": "Request failed", "details": str(e)}

    async def forge_guides(self, target_sequence: str, num_guides: int = 10):
        """
        Calls ZetaForge (evo-service) to generate guide RNA candidates.
        """
        async with httpx.AsyncClient(timeout=600.0) as client:
            logger.info(f"Calling ZetaForge to generate {num_guides} guides for sequence.")
            
            # Step 1: Start the generation job
            start_payload = {"target_sequence": target_sequence, "num_guides": num_guides}
            start_url = f"{ZETA_FORGE_URL}/generate"
            start_response = await self._make_request(client, "POST", start_url, json=start_payload)
            
            if "error" in start_response:
                return start_response
            
            job_id = start_response.get("job_id")
            if not job_id:
                return {"error": "Failed to get job_id from ZetaForge"}

            logger.info(f"ZetaForge job initiated with ID: {job_id}. Polling for results...")

            # Step 2: Poll for the results
            status_url = f"{ZETA_FORGE_URL}/status/{job_id}"
            for _ in range(20): # Poll for up to 10 minutes
                await asyncio.sleep(30)
                status_response = await self._make_request(client, "GET", status_url)
                if "error" in status_response:
                    # Log the polling error but continue, as it might be transient
                    logger.warning(f"Polling ZetaForge status failed: {status_response}. Retrying...")
                    continue
                
                status = status_response.get("status")
                logger.info(f"Polling ZetaForge... current status: '{status}'")
                
                if status == "completed":
                    logger.info("ZetaForge job complete.")
                    return status_response.get("guides", [])
                elif status == "failed":
                    logger.error(f"ZetaForge job failed: {status_response.get('error')}")
                    return {"error": "ZetaForge generation failed", "details": status_response.get('error')}
            
            return {"error": "ZetaForge polling timed out."}

    async def score_guides_with_oracle(self, guides: list[str], context_sequence: str):
        """
        Calls ZetaOracle to score guides for on-target efficacy.
        (Placeholder - returns mock scores)
        """
        logger.warning("ZetaOracle client is a placeholder. Returning mock scores.")
        await asyncio.sleep(1) # Simulate network latency
        mock_scores = []
        for i, guide in enumerate(guides):
            mock_scores.append({
                "guide_sequence": guide,
                "zeta_score": 0.85 + (i * 0.01), # Mock score
                "confidence": 0.95
            })
        return mock_scores

    def _parse_blast_xml(self, xml_string: str, query_len: int, mismatch_threshold: int = 2) -> int:
        """
        Parses BLAST XML output to count significant off-target hits.
        A 'significant' hit is one that is not the query itself and has few mismatches.
        """
        if not xml_string:
            return 999 # Return a high number if BLAST fails
        
        try:
            root = ET.fromstring(xml_string)
            off_target_count = 0
            
            # We look for 'Hit' elements in the XML
            for hit in root.findall('.//Hit'):
                # Each hit has one or more 'Hsp' (High-scoring Pair)
                for hsp in hit.findall('.//Hsp'):
                    # We are interested in the number of mismatches
                    mismatch_node = hsp.find('Hsp_mismatch')
                    identity_node = hsp.find('Hsp_identity')
                    
                    # Handle missing XML elements gracefully
                    mismatches = int(mismatch_node.text) if mismatch_node is not None else 0
                    identity = int(identity_node.text) if identity_node is not None else 0
                    
                    # The on-target hit will have perfect identity and 0 mismatches
                    if identity == query_len and mismatches == 0:
                        continue # This is the on-target sequence, ignore it
                    
                    if mismatches <= mismatch_threshold:
                        off_target_count += 1
                        
            return off_target_count
        except ET.ParseError as e:
            logger.error(f"Failed to parse BLAST XML: {e}")
            return 999

    async def check_off_targets_with_blast(self, guides: list[str]):
        """
        Calls the live BLAST service to check for potential off-targets.
        """
        logger.info(f"Connecting to live BLAST service '{BLAST_SERVICE_URL}'...")

        async def run_search(client, guide):
            logger.info(f"  BLASTing guide: {guide}")
            payload = {"query_sequence": guide}
            blast_url = f"{BLAST_SERVICE_URL}/search"
            result = await self._make_request(client, "POST", blast_url, json=payload, timeout=300.0)

            if "error" in result:
                logger.error(f"BLAST search for {guide} failed: {result.get('details', result['error'])}")
                off_targets = 999
            else:
                off_targets = self._parse_blast_xml(result.get("raw_blast_xml", ""), len(guide))

            return {
                "guide_sequence": guide,
                "off_target_count": off_targets,
                "confidence": 0.99
            }

        async with httpx.AsyncClient(follow_redirects=True) as client:
            tasks = [run_search(client, guide) for guide in guides]
            safety_data = await asyncio.gather(*tasks)
        
        logger.info("BLAST off-target analysis complete.")
        return safety_data

    async def check_immunogenicity(self, guides: list[str]):
        """
        Calls the Immunogenicity service to predict immune response.
        (Placeholder - returns mock data)
        """
        logger.warning("Immunogenicity client is a placeholder. Returning mock immunogenicity data.")
        await asyncio.sleep(1) # Simulate network latency
        mock_immuno_data = []
        for guide in guides:
            mock_immuno_data.append({
                "guide_sequence": guide,
                "immunogenicity_score": 0.05, # Low is good
                "confidence": 0.90
            })
        return mock_immuno_data 

# --- END INLINED COMMAND CENTER CLIENT ---


fastapi_app = FastAPI(
    title="CRISPR Assistant - AI General Command Center",
    version="3.0",
    description="Central nervous system for orchestrating genomic threat analysis and countermeasure design (V3 Async Workflow)."
)

# --- Modal App Definition ---
app = modal.App("crispr-assistant-command-center-v3") 

# FORCE REBUILD v3.1: Added a comment to invalidate the cache and force a rebuild
# with the latest version of the tools directory.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "fastapi",
        "uvicorn",
        "sqlalchemy",
        "psycopg2-binary",
        "httpx",
    )
    # Set up the environment *before* adding local files.
    .env({"PYTHONPATH": "/root"})
    .add_local_dir("tools", remote_path="/root/tools")
    .add_local_dir("services", remote_path="/root/services")
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