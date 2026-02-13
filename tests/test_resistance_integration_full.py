
import asyncio
import sys
import logging
from unittest.mock import MagicMock, AsyncMock

# Add backend to path
sys.path.insert(0, './oncology-coPilot/oncology-backend-minimal')

# Configure logging
logging.basicConfig(level=logging.WARN)
logger = logging.getLogger("IntegrationTest")

# Mock dependencies to avoid DB/Network calls
sys.modules['api.services.supabase_service'] = MagicMock()
sys.modules['api.services.neo4j_connection'] = MagicMock()
sys.modules['api.services.llm_provider.google_genai'] = MagicMock()

# Import Schema (Critical for proper arg passing)
try:
    from api.services.ayesha_care_plan.schemas import CompleteCareV2Request
    from api.services.ayesha_care_plan.orchestrator import AyeshaCarePlanOrchestrator
except ImportError as e:
    print(f"Failed to import: {e}")
    sys.exit(1)

async def test_full_resistance_loop():
    print("\n⚔️ TEST: Full Resistance Prophecy Loop (Detection -> Action)")
    
    orchestrator = AyeshaCarePlanOrchestrator()
    
    # --- MOCKING ---
    # 1. Mock Resistance Service (The Prophet)
    orchestrator.resistance_service = AsyncMock()
    # Scenario: High Risk of Repair Restoration (Resistance to PARP)
    orchestrator.resistance_service.get_resistance_prediction.return_value = {
        "risk_level": "HIGH",
        "signals_detected": [
            {
                "signal_type": "DNA_REPAIR_RESTORATION",
                "detected": True,
                "probability": 0.95,
                "rationale": "Secondary reversion mutation in BRCA1"
            }
        ],
        "integrated_confidence": 0.90
    }
    
    # 2. Mock SAE Service (Phase 2 Computation)
    # We KEEP the real 'get_next_test_recommendations' method by NOT mocking the whole service object
    # But we mock the compute methods attached to it?
    # Orchestrator calls: self.sae_service.compute_sae_features(...)
    # self.sae_service is an instance of AyeshaSAEService.
    # We can mock specific methods on the instance.
    orchestrator.sae_service.compute_sae_features = AsyncMock(return_value={
        "dna_repair_capacity": 0.5, # Moderate, doesn't trigger SAE rule alone
        "hotspot_mutation": False
    })
    orchestrator.sae_service.detect_resistance = MagicMock(return_value={})
    orchestrator.sae_service.get_hint_tiles = MagicMock(return_value={"hint_tiles": []})
    
    # Note: orchestrator.sae_service.get_next_test_recommendations is REAL.
    
    # 3. Mock other services to prevent crashes
    orchestrator.drug_efficacy_service = AsyncMock()
    orchestrator.drug_efficacy_service.get_drug_efficacy.return_value = {"provenance": {}}
    orchestrator.ca125_service = MagicMock()
    orchestrator.ca125_service.get_ca125_intelligence.return_value = {}
    orchestrator.soc_service = MagicMock()
    orchestrator.soc_service.get_soc_recommendation.return_value = {}
    orchestrator.foo_service = AsyncMock() # Typo in orchestrator? 'food_service'
    orchestrator.food_service = AsyncMock()
    
    # --- REQUEST ---
    # Construct pydantic model or mock it if validation is strict
    # Using Pydantic model needed.
    request = CompleteCareV2Request(
        patient_id="AYESHA",
        germline_status="positive", # MBD4 trigger
        ca125_value=35.0,
        stage="IVB",
        treatment_history=[
            {"name": "carboplatin"}, 
            {"name": "paclitaxel"}, 
            {"name": "olaparib"}
        ],
        tumor_context={"hrd_score": 58},
        include_resistance_prediction=True,
        include_wiwfm=False
    )
    
    # --- EXECUTE ---
    try:
        response = await orchestrator.get_complete_care_plan(request)
        
        # --- VERIFICATION ---
        
        # 1. Verify Resistance Prediction was called
        # Check calling args to verify MBD4 baseline injection
        call_args = orchestrator.resistance_service.get_resistance_prediction.call_args
        _, kwargs = call_args
        baseline = kwargs.get('baseline_sae_features')
        
        if baseline and baseline.get('dna_repair_deficiency') == 0.70:
            print("✅ CHECK 1: MBD4 Baseline Injection -> CONFIRMED (Deficiency 0.70)")
        else:
            print(f"❌ CHECK 1 FAILED: Baseline incorrect: {baseline}")

        # 2. Verify Next Test Recommendation updated
        # We look for 'ctDNA' promoted to Priority 1 due to 'Prophet-ALERT'
        recs = response.next_test_recommender.get("recommendations", [])
        
        ctdna_rec = next((r for r in recs if "ctDNA" in r.get("test_name")), None)
        
        if not ctdna_rec:
            print("❌ CHECK 2 FAILED: ctDNA recommendation missing entirely.")
        else:
            priority = ctdna_rec.get("priority")
            rationale = ctdna_rec.get("rationale", "")
            
            if priority == 1 and "[Prophet-ALERT]" in rationale:
                print("✅ CHECK 2: Actionable Intelligence -> CONFIRMED")
                print(f"   - ctDNA promoted to Priority {priority}")
                print(f"   - Rationale: {rationale[:60]}...")
            else:
                print(f"❌ CHECK 2 FAILED: ctDNA not promoted correctly.")
                print(f"   Priority: {priority} (Expected 1)")
                print(f"   Rationale: {rationale}")
                
    except Exception as e:
        print(f"❌ CRASH: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(test_full_resistance_loop())
