#!/usr/bin/env python3
"""
Test Boltz fast-mode via HTTP endpoint
"""
import requests
import time

BOLTZ_URL = "https://crispro--boltz-service-fastapi-app.modal.run"
TEST_SEQUENCE = "MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV"

if __name__ == "__main__":
    print("⚔️ BOLTZ FAST-MODE HTTP TEST")
    print("=" * 80)
    print(f"URL: {BOLTZ_URL}")
    print(f"Sequence: {TEST_SEQUENCE} ({len(TEST_SEQUENCE)} aa)")
    print("=" * 80)
    
    # Test structural integrity endpoint
    payload = {
        "target_sequence": TEST_SEQUENCE,
        "job_id": f"http_test_{int(time.time())}"
    }
    
    print(f"🚀 Starting at {time.ctime()}")
    start_time = time.time()
    
    try:
        response = requests.post(
            f"{BOLTZ_URL}/v1/predict_structural_integrity",
            json=payload,
            timeout=600  # 10 min max
        )
        
        elapsed = time.time() - start_time
        print(f"\n✅ Response received in {elapsed:.1f}s ({elapsed/60:.1f} min)")
        print("=" * 80)
        
        if response.status_code == 200:
            result = response.json()
            print(f"Status: {result.get('status', 'unknown')}")
            print(f"📊 pLDDT Score: {result.get('plddt_score', 0.0):.2f}")
            print(f"📊 PTM Score: {result.get('ptm_score', 0.0):.2f}")
            print(f"📊 Fraction Disordered: {result.get('fraction_disordered', 1.0):.2f}")
            
            plddt = result.get('plddt_score', 0.0)
            if plddt >= 70:
                print("\n✅ PASS: pLDDT ≥ 70 (good confidence)")
            elif plddt >= 50:
                print("\n⚠️  ACCEPTABLE: pLDDT 50-70 (expected for single-sequence mode)")
            else:
                print("\n❌ LOW: pLDDT < 50 (unreliable structure)")
        else:
            print(f"❌ HTTP {response.status_code}")
            print(response.text[:500])
            
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n❌ EXCEPTION after {elapsed/60:.1f} min:")
        print(str(e))


