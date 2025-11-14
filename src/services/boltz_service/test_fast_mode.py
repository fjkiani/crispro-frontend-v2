#!/usr/bin/env python3
"""
Boltz Fast-Mode Smoke Test
Tests msa='empty' single-sequence inference (expect 2-5 min)
"""
import modal
import time

# Lookup the deployed app
app = modal.App.lookup("boltz-service", create_if_missing=False)
run_structure_prediction = modal.Function.lookup("boltz-service", "run_structure_prediction")

# Test sequence: tiny protein (48 aa from test_tiny.fasta)
TEST_SEQUENCE = "MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV"

if __name__ == "__main__":
    print("⚔️ BOLTZ FAST-MODE SMOKE TEST")
    print("=" * 80)
    print(f"Sequence: {TEST_SEQUENCE} ({len(TEST_SEQUENCE)} aa)")
    print("Expected: 2-5 min (vs 60+ min with full MSA)")
    print("=" * 80)
    
    start_time = time.time()
    print(f"🚀 Starting at {time.ctime()}")
    
    try:
        # Call deployed function (no local timeout - server handles it)
        result = run_structure_prediction.remote(TEST_SEQUENCE, job_id="fast_smoke_001")
        
        elapsed = time.time() - start_time
        print("\n" + "=" * 80)
        print(f"✅ COMPLETED IN {elapsed/60:.1f} MIN")
        print("=" * 80)
        
        if result.get('status') == 'complete':
            print(f"✅ Success: Structure predicted!")
            print(f"📊 pLDDT Score: {result.get('plddt_score', 0.0):.2f}")
            print(f"📊 PTM Score: {result.get('ptm_score', 0.0):.2f}")
            print(f"📊 Fraction Disordered: {result.get('fraction_disordered', 1.0):.2f}")
            print(f"⏱️  Runtime: {elapsed:.1f}s ({elapsed/60:.1f} min)")
            
            # Assessment
            plddt = result.get('plddt_score', 0.0)
            if plddt >= 70:
                print("✅ PASS: pLDDT ≥ 70 (good confidence)")
            elif plddt >= 50:
                print("⚠️  ACCEPTABLE: pLDDT 50-70 (expected for single-sequence mode)")
            else:
                print("❌ LOW: pLDDT < 50 (unreliable structure)")
            
        else:
            print(f"❌ Failed: {result.get('error', 'Unknown error')}")
            print(f"Status: {result.get('status', 'unknown')}")
            print(f"Full result: {result}")
            
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n❌ EXCEPTION after {elapsed/60:.1f} min:")
        print(str(e))
        raise

