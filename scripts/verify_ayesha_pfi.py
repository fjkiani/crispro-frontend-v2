
import asyncio
from datetime import datetime
from api.services.ayesha_care_plan.resistance_service import get_ayesha_resistance_service

async def verify_pfi_injection():
    print("🚀 Starting Ayesha PFI Verification (Phase 7)...")
    
    # 1. Synthesize Frontend Payload (matches AYESHA_11_17_25_PROFILE.js)
    treatment_history = [
      {
        "drug": "Carboplatin + Bevacizumab",
        "date": "2025-11-25",
        "type": "Systemic"
      },
      {
        "drug": "Carboplatin + Bevacizumab",
        "date": "2025-12-16",
        "type": "Systemic"
      },
      {
        "drug": "Carboplatin + Bevacizumab",
        "date": "2026-01-06",
        "type": "Systemic"
      },
      {
        "drug": "Carboplatin",
        "date": "2026-01-27",
        "type": "Systemic"
      },
    ]
    
    # 2. Initialize Service
    service = get_ayesha_resistance_service()
    
    # 3. Call Prediction
    print("💉 Injecting Treatment History into Prophet...")
    prediction = await service.get_resistance_prediction(
        patient_id="AYESHA_VERIFY",
        treatment_history=treatment_history,
        current_drug_class="platinum_chemotherapy",
        # Minimal context to pass gates
        current_sae_features={"status": "simulated", "dna_repair_capacity": 0.5},
        baseline_sae_features={"source": "simulated"}
    )
    
    # 4. Verify Output
    if not prediction or prediction.get("status") == "error":
        print("❌ Prediction Failed!")
        print(prediction)
        return

    print("\n✅ Prediction Received:")
    print(f"Risk Level: {prediction.get('risk_level')}")
    print(f"Prophet Logic Rationale: {prediction.get('rationale')}")
    
    # Check for PFI specific rationale
    rationale = str(prediction.get('rationale', ''))
    if "PFI" in rationale:
        print("\n🎉 SUCCESS: PFI Logic Triggered!")
        print(f"Found in rationale: '{rationale}'")
    else:
        print("\n⚠️ WARNING: PFI not mentioned in rationale.")
        
    provenance = prediction.get("provenance", {})
    print(f"\nProvenance Engine: {provenance.get('engine')}")


if __name__ == "__main__":
    asyncio.run(verify_pfi_injection())
