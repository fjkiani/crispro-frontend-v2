
import asyncio
import sys
import json
import time
from datetime import datetime, timedelta

# Import paths
sys.path.insert(0, 'oncology-coPilot/oncology-backend-minimal')
from api.services.resistance.biomarkers.therapeutic.timing_chemo_features import build_timing_chemo_features
from api.services.resistance_prophet_service import get_resistance_prophet_service

async def verify_bridge():
    print("\n=== RESISTANCE PROPHET BRIDGE VERIFICATION (MARS AUDIT) ===\n")
    
    start_total = time.time()
    
    # ---------------------------------------------------------
    # TEST 1: The Golden Path (Engine -> Bridge -> Prediction)
    # ---------------------------------------------------------
    print("TEST 1: FULL PIPELINE (Engine + Bridge)")
    
    # Regimen 1: Frontline (Ends Jan 1 2024)
    r1_end = datetime(2024, 1, 1)
    r1_start = r1_end - timedelta(days=120)
    
    # Regimen 2: Recurrence (Starts Mar 1 2024 - 2 months PFI)
    r2_start = datetime(2024, 3, 1) # ~60 days later
    r2_end = r2_start + timedelta(days=21)
    
    regimen_table = [
        {
            'patient_id': 'TEST_BRIDGE',
            'regimen_id': 'R1',
            'regimen_start_date': r1_start,
            'regimen_end_date': r1_end,
            'regimen_type': 'carboplatin+paclitaxel',
            'line_of_therapy': 1,
            'setting': 'frontline',
        },
        {
            'patient_id': 'TEST_BRIDGE',
            'regimen_id': 'R2',
            'regimen_start_date': r2_start,
            'regimen_end_date': r2_end,
            'regimen_type': 'doxil', # 2nd line
            'line_of_therapy': 2,
            'setting': 'recurrence',
        }
    ]
    
    t0 = time.time()
    # 2. Run Timing Engine
    # MARS AUDIT: Simulate Network Latency (Internal Function -> API Call equivalent)
    time.sleep(0.35) # 350ms simulated network fetch
    
    timing_features = build_timing_chemo_features(
        regimen_table=regimen_table,
        survival_table=[],
        clinical_table=[{'patient_id': 'TEST_BRIDGE', 'disease_site': 'ovary'}]
    )
    t1 = time.time()
    engine_latency = (t1-t0)*1000
    print(f"  [Timing Engine] Fetched {len(timing_features)} regimens in {engine_latency:.1f}ms (incl. simulated network)")
    
    # Validate Output (Check R2 which should have PFI)
    r2 = timing_features[1] 
    pfi_days = r2.get('PFI_days')
    pfi_months = pfi_days / 30.44 if pfi_days is not None else 0.0
    print(f"  [Input Data]    PFI_days={pfi_days} ({pfi_months:.2f} months)")

    # WRAP into Endpoint Format 
    timing_context_payload = {
        "timing_features_table": timing_features,
        "provenance": {
            "computed_at": datetime.now().isoformat(),
            "engine": "timing_chemosensitivity_engine",
            "ruo_disclaimer": "TEST ADAPTER",
        }
    }
    
    # 3. Run Resistance Prophet
    service = get_resistance_prophet_service()
    
    t2 = time.time()
    prediction = await service.predict_resistance(
        patient_id='TEST_BRIDGE',
        disease_type='ovarian_cancer',
        current_sae_features={}, 
        baseline_sae_features={}, 
        treatment_history=[{'drug': 'carboplatin'}], 
        current_drug_class='platinum', 
        mode='heuristic',
        timing_context=timing_context_payload # BRIDGE
    )
    t3 = time.time()
    print(f"  [Prophet Svc]   Prediction in {(t3-t2)*1000:.1f}ms")
    
    # 4. Assertions
    signals = prediction['signalsdetected']
    found_engine = False
    for s in signals:
        if s['signaltype'] == 'MMDRUGCLASSRESISTANCE':
             if "Engine" in s['provenance']['source'] or "timing" in s['provenance']['source']:
                 found_engine = True
                 print(f"  [Output]        Risk={prediction['risklevel']} (Conf={prediction['confidence']})")
                 print(f"  [Signal Logic]  {s['rationale']}")
                 print(f"  [Provenance]    {s['provenance']}")

    if found_engine:
        print("  ✅ TEST 1 PASS: Engine Data Propagated Sucessfully\n")
    else:
        print("  ❌ TEST 1 FAIL: Engine Signal Missing\n")
        
    # ---------------------------------------------------------
    # TEST 2: Fallback (Engine Down)
    # ---------------------------------------------------------
    print("TEST 2: FALLBACK (Engine Unavailable)")
    
    t4 = time.time()
    prediction_fallback = await service.predict_resistance(
        patient_id='TEST_FALLBACK',
        disease_type='ovarian_cancer',
        current_sae_features={}, 
        baseline_sae_features={}, 
        treatment_history=[{'drug': 'carboplatin'}], 
        current_drug_class='platinum',
        mode='heuristic',
        timing_context=None # SIMULATE OUTAGE
    )
    t5 = time.time()
    
    # Check if we got a prediction (even if low confidence)
    # Without engine, PFI is unknown or heuristic dependent.
    # Just checking it didn't crash and returns valid object.
    if prediction_fallback['risklevel']:
         print(f"  [Fallback]      Prediction returned in {(t5-t4)*1000:.1f}ms")
         print(f"  [Fallback]      Risk={prediction_fallback['risklevel']}")
         print("  ✅ TEST 2 PASS: Service handled missing context gracefully\n")
    else:
         print("  ❌ TEST 2 FAIL: Service crashed or returned invalid data\n")

    # ---------------------------------------------------------
    # PERFORMANCE SUMMARY
    # ---------------------------------------------------------
    
    end_total = time.time()
    print("\n=== PERFORMANCE METRICS (FULL PIPELINE) ===")
    print(f"Timing Engine fetch:   {engine_latency:.0f}ms")
    # print(f"Adapter conversion:    12ms") # Implicit in Prophet Service time
    print(f"Prophet prediction:    {(t3-t2)*1000:.0f}ms")
    print(f"Total user-facing latency:  {(end_total-start_total)*1000:.0f}ms")
    
    if (end_total-start_total) < 3.0:
        print("✅ LATENCY CHECK PASS (<3s)")
    else:
        print("⚠️ LATENCY WARNING (>3s)")

if __name__ == "__main__":
    asyncio.run(verify_bridge())
