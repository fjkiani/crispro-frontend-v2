
import asyncio
import sys
import json
from pathlib import Path

# Add project root to path
BACKEND_ROOT = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal")
sys.path.insert(0, str(BACKEND_ROOT))

# Mock LLM API to avoid import errors if not needed for ranking or use real one
REPO_ROOT = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main")
sys.path.insert(0, str(REPO_ROOT / "src" / "tools"))

from api.services.ayesha_care_plan.trial_service import AyeshaTrialService
from scripts.data.ayesha_patient_profile import get_ayesha_complete_profile

async def audit_top_trials():
    print("🔍 AUDITING TOP 5 TRIALS FOR AYESHA FIT...")
    
    # 1. Initialize Service
    service = AyeshaTrialService()
    
    # 2. Get Ayesha Profile & Construct Mock Request
    ayesha = get_ayesha_complete_profile()
    
    # Mock request object matching what service expects
    class MockRequest:
        def __init__(self):
            self.max_trials = 1200 
            self.age = 40
            self.location_state = "NY"
            self.germline_status = "positive" 
            self.treatment_history = []
            self.tumor_context = {
                "disease_type": "ovarian_cancer_hgs",
                "somatic_mutations": [
                    {"gene": "TP53", "variant": "mutant"},
                    {"gene": "MBD4", "variant": "loss_of_function"} 
                ],
                "biomarkers": {
                    "p53_status": "MUTANT_TYPE",
                    "pd_l1_status": "POSITIVE",
                    "mmr_status": "PRESERVED" 
                }
            }
            
    req = MockRequest()
    
    # 3. Get Trials
    # MATCH generator script vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    mech_vector = [0.88, 0.12, 0.15, 0.10, 0.05, 0.2, 0.0]
    
    print("   Calling service.get_trials()...")
    results = await service.get_trials(req, mech_vector)
    trials = results.get("trials", [])
    
    print(f"✅ Retrieved {len(trials)} trials.")
    
    print("\n🏆 TOP 5 TRIALS:")
    for i, t in enumerate(trials[:5], 1):
        print(f"\n{i}. {t.get('nct_id')} - Score: {t.get('holistic_score', 0):.4f}")
        print(f"   Title: {t.get('title')}")
        print(f"   Phase: {t.get('phase')}")
        print(f"   Status: {t.get('status')}")
        print(f"   Rationale: {t.get('rationale')}")
        
        # Check relevance keywords
        text = (t.get('title', '') + t.get('brief_summary', '') + t.get('eligibility_text', '')).lower()
        print(f"   Relevant Keywords Found:")
        if 'ovarian' in text: print("   - Ovarian")
        if 'platinum' in text: print("   - Platinum")
        if 'brca' in text: print("   - BRCA")
        if 'recu' in text: print("   - Recurrent/Recurring")
        
        # Check if dossier exists
        dossier_path = Path(f"/Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal/.cursor/ayesha/zo_fresh_dossiers/INTELLIGENCE_{t.get('nct_id')}_TOP_TIER.md")
        if dossier_path.exists():
            print(f"   ✅ Dossier Found: {dossier_path.name}")
            # Quick content peek
            content = dossier_path.read_text()
            if "TACTICAL FIT ANALYSIS" in content:
                print("      - Contains 'TACTICAL FIT ANALYSIS'")
            else:
                print("      ❌ MISSING 'TACTICAL FIT ANALYSIS' Header")
        else:
            print(f"   ❌ Dossier MISSING")

if __name__ == "__main__":
    asyncio.run(audit_top_trials())
