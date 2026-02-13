
import asyncio
import logging
import sys
import os
from types import SimpleNamespace
from typing import List, Dict, Any

# Ensure we can import from api
sys.path.append(os.getcwd())

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Mock the Request Object
class MockCompleteCareRequest:
    def __init__(self):
        self.max_trials = 20
        self.age = 40
        self.location_state = "NY"
        self.germline_status = "positive" # MBD4
        self.treatment_history = [] # Naive
        self.tumor_context = {
            "disease_type": "ovarian_cancer_hgs",
            "biomarkers": {
                "p53_status": "MUTANT_TYPE",
                "pd_l1_status": "POSITIVE",
                "mmr_status": "PRESERVED"
            },
            "somatic_mutations": [
                {"gene": "TP53", "variant": "mutant"},
                {"gene": "MBD4", "variant": "loss_of_function"} # This is the key driver
            ]
        }

async def run_simulation():
    print("--- STARTING AYSEHHA RANKING SIMULATION ---")
    
    try:
        from api.services.ayesha_care_plan.trial_service import get_ayesha_trial_service
        
        service = get_ayesha_trial_service()
        request = MockCompleteCareRequest()
        
        # Ayesha's Mechanism Vector (DDR High)
        # [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
        # MBD4 + TP53 -> High DDR
        mechanism_vector = [0.88, 0.12, 0.15, 0.10, 0.05, 0.2, 0.0]
        
        print(f"Requesting top {request.max_trials} trials...")
        response = await service.get_trials(request, mechanism_vector)
        
        trials = response.get("trials", [])
        print(f"Received {len(trials)} trials.")
        
        print("\n--- TOP 10 TRIALS FOR AYESHA ---")
        for i, t in enumerate(trials[:10]):
            print(f"\nRANK {i+1}: {t.get('nct_id')} - {t.get('title')}")
            print(f"   HOLISTIC SCORE: {t.get('holistic_score'):.3f} ({t.get('holistic_interpretation')})")
            print(f"   BREAKDOWN:")
            print(f"      - Mechanism Fit: {t.get('mechanism_fit_score'):.3f}")
            print(f"      - Eligibility:   {t.get('eligibility_score'):.3f}")
            print(f"      - PGx Safety:    {t.get('pgx_safety_score'):.3f}")
            print(f"      - Resist. Risk:  {t.get('resistance_risk_score'):.3f} (Lower Resistance=Higher Score)")
            
            if t.get('holistic_caveats'):
                print(f"   CAVEATS: {t.get('holistic_caveats')}")
                
    except Exception as e:
        logger.exception("Simulation failed")

if __name__ == "__main__":
    asyncio.run(run_simulation())
