#!/usr/bin/env python3
"""
🔍 Pathway Score Validation Health Check

Purpose: Verify pathway scores match known biology (BRCA1→DDR, KRAS→MAPK, etc.)
Agent: Zo (Lead Commander)
Date: January 20, 2025
"""

import asyncio
import httpx
import sys

async def check_pathway_validation():
    """Verify pathway scores match known biology."""
    
    print("=" * 80)
    print("HEALTH CHECK: Pathway Score Validation")
    print("=" * 80)
    
    base_url = "http://localhost:8000"
    
    # Known biology test cases
    test_cases = [
        {
            "name": "BRCA1 biallelic",
            "mutations": [{"gene": "BRCA1", "hgvs_p": "p.Q1756*", "chrom": "17", "pos": 43044295, "ref": "T", "alt": "G", "build": "GRCh38"}],
            "expected_pathway": "ddr",
            "expected_min": 0.70  # BRCA1 loss → high DDR
        },
        {
            "name": "KRAS G12D",
            "mutations": [{"gene": "KRAS", "hgvs_p": "p.G12D", "chrom": "12", "pos": 25398284, "ref": "G", "alt": "A", "build": "GRCh38"}],
            "expected_pathway": "mapk",
            "expected_min": 0.70  # KRAS mutation → high MAPK
        },
        {
            "name": "MBD4+TP53",
            "mutations": [
                {"gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2", "chrom": "3", "pos": 129430456, "ref": "A", "alt": "", "build": "GRCh38"},
                {"gene": "TP53", "hgvs_p": "p.R175H", "chrom": "17", "pos": 7577120, "ref": "G", "alt": "A", "build": "GRCh38"}
            ],
            "expected_pathway": "ddr",
            "expected_min": 0.70  # MBD4 BER + TP53 checkpoint → high DDR
        }
    ]
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            passed = 0
            failed = 0
            
            for i, test_case in enumerate(test_cases, 1):
                print(f"\n[{i}/{len(test_cases)}] Testing {test_case['name']}...")
                
                try:
                    response = await client.post(
                        f"{base_url}/api/efficacy/predict",
                        json={
                            "model_id": "evo2_1b",
                            "mutations": test_case["mutations"],
                            "options": {"adaptive": True}
                        }
                    )
                    
                    if response.status_code != 200:
                        print(f"❌ FAIL: Efficacy returned {response.status_code}")
                        failed += 1
                        continue
                    
                    data = response.json()
                    
                    # Extract pathway score
                    pathway_score = 0.0
                    if "provenance" in data and "confidence_breakdown" in data["provenance"]:
                        pathway_disruption = data["provenance"]["confidence_breakdown"].get("pathway_disruption", {})
                        pathway_score = pathway_disruption.get(test_case["expected_pathway"], 0.0)
                    
                    # Validate against expected
                    if pathway_score >= test_case["expected_min"]:
                        print(f"✅ PASS: {test_case['expected_pathway'].upper()} = {pathway_score:.3f} (expected ≥{test_case['expected_min']})")
                        passed += 1
                    else:
                        print(f"⚠️  WARNING: {test_case['expected_pathway'].upper()} = {pathway_score:.3f} (expected ≥{test_case['expected_min']})")
                        print(f"   Known biology: {test_case['name']} should have high {test_case['expected_pathway'].upper()}")
                        failed += 1
                        
                except Exception as e:
                    print(f"❌ FAIL: Test case error: {e}")
                    failed += 1
            
            print(f"\n{'=' * 80}")
            print(f"RESULTS: {passed} passed, {failed} failed out of {len(test_cases)}")
            print(f"{'=' * 80}")
            
            return failed == 0
            
    except Exception as e:
        print(f"❌ FAIL: Pathway validation failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = asyncio.run(check_pathway_validation())
    sys.exit(0 if success else 1)

