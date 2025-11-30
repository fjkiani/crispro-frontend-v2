#!/usr/bin/env python3
"""
🔍 SAE Backend Services Health Check

Purpose: Verify SAE backend endpoints are operational
Agent: Zo (Lead Commander)
Date: January 20, 2025
"""

import asyncio
import httpx
import sys

async def check_backend_services():
    """Verify SAE backend endpoints are operational."""
    
    print("=" * 80)
    print("HEALTH CHECK: SAE Backend Services")
    print("=" * 80)
    
    base_url = "http://localhost:8000"
    
    endpoints_to_check = [
        ("/api/sae/compute_features", "POST", {
            "pathway_scores": {"ddr": 0.8, "mapk": 0.2, "pi3k": 0.1, "vegf": 0.3, "her2": 0.0},
            "insights_bundle": {"functionality": 0.6, "chromatin": 0.5, "essentiality": 0.7, "regulatory": 0.4},
            "tumor_context": {"hrd_score": 0.75, "tmb_score": 25.0}
        }),
        ("/healthz", "GET", None),
    ]
    
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            for endpoint, method, payload in endpoints_to_check:
                print(f"\n[{endpoints_to_check.index((endpoint, method, payload)) + 1}/{len(endpoints_to_check)}] Checking {method} {endpoint}...")
                
                try:
                    if method == "GET":
                        response = await client.get(f"{base_url}{endpoint}")
                    else:
                        response = await client.post(f"{base_url}{endpoint}", json=payload)
                    
                    if response.status_code == 200:
                        print(f"✅ PASS: {endpoint} returned 200")
                        if endpoint == "/api/sae/compute_features":
                            data = response.json()
                            if "dna_repair_capacity" in data:
                                print(f"   DNA repair capacity: {data['dna_repair_capacity']:.3f}")
                            if "mechanism_vector" in data:
                                print(f"   Mechanism vector: {len(data['mechanism_vector'])}D")
                    else:
                        print(f"❌ FAIL: {endpoint} returned {response.status_code}")
                        print(f"   Response: {response.text[:200]}")
                        return False
                        
                except httpx.TimeoutException:
                    print(f"❌ FAIL: {endpoint} timed out")
                    return False
                except Exception as e:
                    print(f"❌ FAIL: {endpoint} error: {e}")
                    return False
            
            print("\n✅ All backend services operational")
            return True
            
    except Exception as e:
        print(f"❌ FAIL: Backend check failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = asyncio.run(check_backend_services())
    sys.exit(0 if success else 1)

