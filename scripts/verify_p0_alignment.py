import sys
import os
import asyncio
import json

# Setup path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../oncology-coPilot/oncology-backend-minimal')))

from api.services.ayesha_fit.scenarios import L2_SCENARIOS, L3_SCENARIOS
from api.services.ayesha_fit.builder import build_profile_for_level
from api.services.ayesha_fit.service import run_therapy_fit_analysis

async def check_1_scenarios_endpoint():
    print("\n--- Check 1: Scenarios Endpoint (Simulated) ---")
    # Simulation of @router.get("/scenarios")
    l2_cards = []
    for sid, scn in L2_SCENARIOS.items():
        l2_cards.append({
            "id": sid,
            "meta": scn
        })
    
    # Check if Sig7 is in meta
    s1 = next((s for s in l2_cards if s["id"] == "S1HRDhighTMBhigh"), None)
    if s1:
        sig7 = s1["meta"].get("tumor_context_additions", {}).get("sig7_exposure")
        print(f"S1HRDhighTMBhigh Canonical ID found: YES")
        print(f"S1HRDhighTMBhigh Sig7 in Meta: {sig7}")
        if sig7 is None: print("FAIL: Sig7 missing from S1 metadata")
    else:
        print("FAIL: S1HRDhighTMBhigh not found (Canonical ID mismatch?)")

async def check_2_bundle_endpoint():
    print("\n--- Check 2: Bundle Endpoint (Simulated) ---")
    # Simulation of @router.post("/bundle") for L2 S3
    # Simulation of @router.post("/bundle") for L2 S3
    sid = "S3HRDlowTMBlow"
    try:
        # We only need to check what BUILDER produces, not what Service returns (which just echoes inputs anyway)
        # Service calls build_profile_for_level -> build_tumor_context_for_api
        
        # 1. Build Profile
        profile = build_profile_for_level("l2", scenario_id=sid)
        
        # 2. Extract Context (Simulating what Service sees)
        tumor_ctx = profile.get("tumor_context", {})
        sig7 = tumor_ctx.get("sig7_exposure")
        
        print(f"Bundle Builder Output for S3HRDlowTMBlow:")
        print(f"  - HRD: {tumor_ctx.get('hrd_score')}")
        print(f"  - Sig7: {sig7}")
        
        if sig7 is None:
            print("FAIL: Sig7 NOT passed through to 'tumor_context_used'")
        elif sig7 == 0.88:
            print("PASS: Sig7 (0.88) correctly passed through")
        else:
            print(f"WARN: Sig7 present but value mismatch: {sig7}")
            
    except Exception as e:
        print(f"FAIL: Bundle simulation crashed: {e}")

if __name__ == "__main__":
    asyncio.run(check_1_scenarios_endpoint())
    asyncio.run(check_2_bundle_endpoint())
