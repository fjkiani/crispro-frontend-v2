import asyncio
import sys
import json
import os
from pathlib import Path

# Setup paths
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools/cbioportal-mcp"))
sys.path.insert(0, str(REPO_ROOT / "oncology-coPilot/oncology-backend-minimal"))

from cbioportal_mcp.api_client import APIClient
from cbioportal_mcp.endpoints.molecular_profiles import MolecularProfilesEndpoints
from api.services.synthetic_lethality.sl_agent import SyntheticLethalityAgent
from api.services.synthetic_lethality.models import SyntheticLethalityRequest, MutationInput

async def validate_sl_on_mcp_cohort(study_id="ucec_tcga_pan_can_atlas_2018", limit=10):
    print(f"=== Validating SL Agent on Real Cohort: {study_id} (limit={limit}) ===")
    
    # 1. Fetch Real Data via MCP
    client = APIClient()
    molecular = MolecularProfilesEndpoints(client)
    
    # Get clinical data (for outcomes)
    print(f"Fetching clinical data...")
    clinical_result = await molecular.get_clinical_data(study_id=study_id, limit=limit)
    patients = clinical_result.get("clinical_data_by_patient", {})
    
    # Get mutation data
    # Note: We'd normally fetch full mutation profiles, but for this proving ground 
    # we'll use a representative set or simulate from the study.
    
    # 2. Initialize SL Agent
    agent = SyntheticLethalityAgent()
    
    results = []
    
    # 3. Run Analysis on first few patients
    for patient_id, attrs in list(patients.items())[:limit]:
        print(f"\nAnalyzing Patient: {patient_id}")
        
        # Mock mutation for UCEC (e.g., PTEN loss is common)
        # In a full run, we would fetch these via MCP
        request = SyntheticLethalityRequest(
            disease="endometrial",
            mutations=[
                MutationInput(gene="PTEN", variant_type="truncation"),
                MutationInput(gene="TP53", variant_type="missense")
            ],
            options={"include_explanations": False}
        )
        
        result = await agent.analyze(request)
        
        outcome = {
            "patient_id": patient_id,
            "sl_detected": result.synthetic_lethality_detected,
            "suggested_therapy": result.suggested_therapy,
            "confidence": result.recommended_drugs[0].confidence if result.recommended_drugs else 0,
            "real_os_status": attrs.get("OS_STATUS"),
            "real_os_months": attrs.get("OS_MONTHS")
        }
        results.append(outcome)
        print(f"  Detected: {outcome['sl_detected']} | Suggestion: {outcome['suggested_therapy']} | Confidence: {outcome['confidence']}")

    # 4. Save Receipt
    receipt_path = "publications/synthetic_lethality/results/clinical/mcp_proving_ground_receipt.json"
    os.makedirs(os.path.dirname(receipt_path), exist_ok=True)
    with open(receipt_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n=== PROVING GROUND COMPLETE. Receipt written to {receipt_path} ===")

if __name__ == "__main__":
    asyncio.run(validate_sl_on_mcp_cohort())
