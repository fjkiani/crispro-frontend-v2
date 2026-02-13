
import asyncio
import sys
import logging
from unittest.mock import MagicMock, AsyncMock

# Add backend to path
sys.path.insert(0, './oncology-coPilot/oncology-backend-minimal')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("SimulationTest")

# Mock dependencies
sys.modules['api.services.supabase_service'] = MagicMock()
sys.modules['api.services.neo4j_connection'] = MagicMock()
sys.modules['api.services.llm_provider.google_genai'] = MagicMock()

from api.services.ayesha_care_plan.resistance_service import AyeshaResistanceService
from api.services.resistance_prophet.schemas import ResistancePrediction, ResistanceRiskLevel, UrgencyLevel, ResistanceSignalData, ResistanceSignal

async def test_simulation():
    print("\n⚔️ TEST: Resistance Lab Simulation Engine")
    
    # 1. Setup Service and Logic
    service = AyeshaResistanceService()
    service._prophet_service = AsyncMock()
    
    # Scenario: MBD4+ Patient with Restoration (Deficiency dropped from 0.70 to 0.30)
    # This simulates High Probability Restoration
    
    # Define Mock Prediction Object
    mock_prediction = ResistancePrediction(
        risk_level=ResistanceRiskLevel.HIGH,
        probability=0.92,
        confidence=0.85,
        signals_detected=[
            ResistanceSignalData(
                signal_type=ResistanceSignal.DNA_REPAIR_RESTORATION,
                detected=True,
                probability=0.92,
                confidence=0.90,
                rationale="Deficiency drop > 0.30 indicates restoration",
                provenance={}
            )
        ],
        signal_count=1,
        urgency=UrgencyLevel.CRITICAL,
        recommended_actions=[],
        next_line_options=[],
        rationale=[],
        provenance={},
        warnings=["INSUFFICIENT_CA-125_DATA"],
        baseline_source="simulated",
        baseline_penalty_applied=False
    )
    
    service._prophet_service.predict_resistance.return_value = mock_prediction
    
    # 2. Run Simulation
    params = {
        "simulate_germline": "positive",
        "simulate_hrd": 70, # Implies deficiency 0.3 (Proficient)
        "simulate_repair_capacity": 0.70, # Explicitly Proficient
        "simulate_treatment": ["carboplatin", "olaparib"]
    }
    
    result = await service.simulate_resistance_scenarios(params)
    
    # 3. Validation
    # A. Check Structure
    assert "prediction" in result
    assert "recommended_tests" in result
    assert "logic_steps" in result
    
    tests = result['recommended_tests']
    steps = result['logic_steps']
    
    print(f"\n✅ Simulation Result Generated")
    print(f"   - Tests Recommended: {len(tests)}")
    print(f"   - Logic Steps: {len(steps)}")
    
    # B. Validate Test Orders
    # Expect: 
    # 1. CA125_SERIES (due to warning)
    # 2. LIQUID_BIOPSY (due to Restoration)
    # 3. TISSUE_HRD (due to Restoration Prob > 0.8)
    # 4. IMAGING_ESCALATE (due to High Risk)
    
    test_ids = [t['test_id'] for t in tests]
    print(f"   - Test IDs: {test_ids}")
    
    expected_ids = ["CA125_SERIES", "LIQUID_BIOPSY", "TISSUE_HRD", "IMAGING_ESCALATE"]
    all_present = all(ex in test_ids for ex in expected_ids)
    
    if all_present:
        print("✅ ALL EXPECTED ORDERS PRESENT")
    else:
        print(f"❌ MISSING ORDERS. Expected {expected_ids}, got {test_ids}")
        
    # C. Validate Logic Log
    # Should see "GAP_DETECTED" and "SIGNAL_ACTION"
    events = [s['event'] for s in steps]
    print(f"   - Events: {events}")
    
    if "SIGNAL_ACTION" in events:
        print("✅ Logic Stream captured SIGNAL_ACTION")
    else:
        print("❌ Logic Stream missing events")

    # D. Validate Interval Logic
    # CA125 should be WEEKLY x3
    ca125 = next(t for t in tests if t['test_id'] == "CA125_SERIES")
    if ca125['frequency'] == "WEEKLY" and ca125['duration'] == "x3":
         print("✅ CA-125 Interval Logic Correct (WEEKLY x3)")
    else:
         print(f"❌ CA-125 Interval Incorrect: {ca125}")

if __name__ == "__main__":
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(test_simulation())
