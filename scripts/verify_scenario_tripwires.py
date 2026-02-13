import sys
import os
import asyncio

# Setup path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../oncology-coPilot/oncology-backend-minimal')))

from api.services.ayesha_fit.validation import validate_demo_scenarios
from api.services.ayesha_fit.scenarios import L2_SCENARIOS, L3_SCENARIOS

def verify_sig7_presence():
    print("--- Verifying Sig7 Injection ---")
    
    # Check L2
    s1_l2 = L2_SCENARIOS.get("S1HRDhighTMBhigh")
    if s1_l2 and "sig7_exposure" in s1_l2["tumor_context_additions"]:
        print(f"PASS: S1HRDhighTMBhigh has sig7_exposure: {s1_l2['tumor_context_additions']['sig7_exposure']}")
    else:
        print("FAIL: S1HRDhighTMBhigh missing sig7_exposure")

    # Check L3
    s1_l3 = L3_SCENARIOS.get("S1VEGFhighCA125high")
    if s1_l3 and "sig7_exposure" in s1_l3["tumor_context_additions"]:
        print(f"PASS: S1VEGFhighCA125high has sig7_exposure: {s1_l3['tumor_context_additions']['sig7_exposure']}")
    else:
        print("FAIL: S1VEGFhighCA125high missing sig7_exposure")

async def main():
    print("--- Running Validator Tripwires ---")
    try:
        # validate_demo_scenarios is synchronous in the file I read? Or async?
        # Checking imports... it was imported from validation.py
        # Let's check validation.py content if needed, but assuming sync for now based on name.
        # Wait, usually service calls are async. Let's try running it.
        # Actually view_file of validation.py would confirm. 
        # But let's assume it checks the dicts directly.
        
        results = validate_demo_scenarios()
        
        failed_l2 = results["l2"]["failed"]
        failed_l3 = results["l3"]["failed"]
        
        if not failed_l2 and not failed_l3:
            print(f"PASS: All scenarios validated (L2: {len(results['l2']['ok'])}, L3: {len(results['l3']['ok'])})")
        else:
            print("FAIL: Scenario validation errors detected")
            if failed_l2: print(f"L2 Failures: {failed_l2}")
            if failed_l3: print(f"L3 Failures: {failed_l3}")
            
    except Exception as e:
        print(f"CRITICAL FAIL: {e}")
        import traceback
        traceback.print_exc()

    verify_sig7_presence()

if __name__ == "__main__":
    asyncio.run(main())
