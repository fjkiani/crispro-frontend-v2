import sys
import os
import asyncio
import json
from datetime import datetime

# Setup path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../oncology-coPilot/oncology-backend-minimal')))

from api.services.ayesha_fit.builder import (
    build_profile_for_level, 
    determine_required_tests,
    build_mutations_for_api,
    build_tumor_context_for_api,
    calculate_completeness_level
)
from api.services.ayesha_fit.utils import _infer_effective_assembly

async def generate_full_bundle():
    # Simulate Router Logic
    levels = ["l1", "l2", "l3"]
    scenario_id = "S3HRDlowTMBlow" # 'Worst Case'
    l3_scenario_id = "S1VEGFhighCA125high" 
    
    out = {}
    
    # Run analysis for each level
    for lk in levels:
        # 1. Build Profile using REAL builder logic
        profile = build_profile_for_level(lk, scenario_id=scenario_id, l3_scenario_id=l3_scenario_id)
        
        # 2. Inject Genome Build hint if missing (simulating upstream data loader)
        # The user requested effective_assembly be GRCh38 if unknown.
        # Check if profile has it, if not, inject into tumor_context for this sim
        if not profile.get("tumor_context", {}).get("genome_build"):
            if "tumor_context" not in profile: profile["tumor_context"] = {}
            profile["tumor_context"]["genome_build"] = "GRCh38"

        # 3. Build API Inputs using REAL builder logic
        mutations_used = build_mutations_for_api(profile)
        tumor_context_used = build_tumor_context_for_api(profile)
        
        # Correctly set completeness according to level rules mentioned by user
        # In builder.py:
        # L2 sets profile["tumor_context"]["completeness_score"] to scenario value
        # L3 sets it to scenario value (0.92?) - checking builder.py logic... 
        # Actually builder.py `build_profile_for_level` merges scenario["tumor_context_additions"] 
        # which has completeness_score.
        
        # 4. Calculate Completeness Object using REAL logic
        completeness_obj = calculate_completeness_level(profile)

        # 5. Infer Assembly using REAL logic
        effective_assembly = _infer_effective_assembly(profile, mutations_used, tumor_context_used)

        # 6. Mock ONLY the Efficacy/Execution part (Model 1b)
        efficacy_mock = {
            "drugs": [
                {"name": "MockParib", "confidence": 0.95, "rationale": "High HRD inferred from Sig7"},
                {"name": "MockPlatin", "confidence": 0.88, "rationale": "Platinum sensitivity retention"}
            ],
            "pathway_scores": {"HRD": 9.5, "CellCycle": 8.0},
            "provenance": {"model": "evo2_1b", "version": "v2.1"}
        }

        # 7. Construct Level Output matching Router
        out[lk.upper()] = {
            "is_preview": bool(lk != "l1"),
            "effective_assembly": effective_assembly,
            "inputs_used": {"mutations": mutations_used, "tumor_context": tumor_context_used},
            "efficacy": efficacy_mock,
            "completeness": completeness_obj, # Full object!
        }

    # Generate Tests Needed using L1 baseline
    # from api.services.ayesha_fit.builder import build_profile_for_level # REMOVED redundant import
    base_profile = build_profile_for_level("l1") 
    tests_needed = determine_required_tests(base_profile)

    payload = {
        "patient_id": "AK",
        "requested_levels": [k.upper() for k in levels],
        "generated_at": datetime.now().isoformat(),
        "levels": out,
        "tests_needed": tests_needed
    }
    
    print(json.dumps(payload, indent=2, default=str))

if __name__ == "__main__":
    asyncio.run(generate_full_bundle())
