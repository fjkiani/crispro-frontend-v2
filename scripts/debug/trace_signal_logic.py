import asyncio
import sys
import os
import json
from datetime import datetime

# Adjust path
sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))

from api.services.resistance_prophet_service import ResistanceProphetService
from api.services.resistance_prophet.schemas import ResistanceRiskLevel

async def trace_logic():
    print("=== Signal Logic Trace (OV-045 Replay) ===")
    
    # 1. Instantiate Service WITHOUT CA125 Service (Simulate current state)
    service = ResistanceProphetService(
        sae_service=None,
        ca125_service=None, # Explicitly None
        treatment_line_service=None,
        resistance_playbook_service=None
    )
    
    # 2. Prepare Input (Day 469: 65.0)
    # Context: Recurrence Event. Previous CA125 was 14.0.
    simulated_context = {
        "ca125_series": [
            {"date": "Day 362", "value": 14.0},
            {"date": "Day 469", "value": 65.0} # SPIKE
        ],
        "active_drugs": ["Chemotherapy"]
    }
    
    # 3. Predict
    result = await service.predict_resistance(
        patient_id="OV-045",
        disease_type="ovarian",
        current_sae_features={"dna_repair_deficiency": 0.5}, # Neutral
        baseline_sae_features={"dna_repair_deficiency": 0.5},
        treatment_history=["Chemotherapy"], 
        current_drug_class="Chemotherapy",
        return_object=True
    )
    
    # 4. Analyze Internal State
    # Inspect which signals were actually detected vs evaluated
    trace = {
        "patient": "OV-045",
        "day": 469,
        "ca125_input": 65.0,
        "service_state": {
            "ca125_service": "None",
            "sae_service": "None"
        },
        "prediction": {
            "risk_level": str(result.risk_level),
            "probability": result.probability
        },
        "signals_detected": [s.signal_type.value for s in result.signals_detected if s.detected],
        "all_signals_evaluated": [s.signal_type.value for s in result.signals_detected]
    }
    
    # Check if CA125_KINETICS is in "all_signals_evaluated"
    trace["logic_verdict"] = "CA125_Logic_Missing"
    if "CA125_KINETICS" in trace["all_signals_evaluated"]:
        trace["logic_verdict"] = "CA125_Logic_Present_But_Negative"
    
    print(json.dumps(trace, indent=2))
    
    # Write to report
    os.makedirs("reports", exist_ok=True)
    with open("reports/signal_trace_patient_OV-045_day_469.json", "w") as f:
        json.dump(trace, f, indent=2)

if __name__ == "__main__":
    asyncio.run(trace_logic())
