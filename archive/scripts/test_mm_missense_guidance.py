#!/usr/bin/env python3
"""
MM Missense Guidance Test Suite
Test key MM mutations against appropriate drug classes to demonstrate fused S scoring
"""
import json
import httpx
import asyncio
from typing import Dict, Any, List

# Test cases: mutation -> drug class mapping
TEST_CASES = [
    {
        "name": "BRAF V600E → BRAF inhibitor",
        "mutation": {"gene": "BRAF", "hgvs_p": "V600E", "chrom": "7", "pos": 140453136, "ref": "T", "alt": "A"},
        "drug_class": "BRAF inhibitor",
        "expected_coverage": "AlphaMissense covered"
    },
    {
        "name": "KRAS G12D → MEK inhibitor", 
        "mutation": {"gene": "KRAS", "hgvs_p": "G12D", "chrom": "12", "pos": 25398284, "ref": "G", "alt": "A"},
        "drug_class": "MEK inhibitor",
        "expected_coverage": "AlphaMissense covered"
    },
    {
        "name": "NRAS Q61R → MEK inhibitor",
        "mutation": {"gene": "NRAS", "hgvs_p": "Q61R", "chrom": "1", "pos": 115256529, "ref": "A", "alt": "G"},
        "drug_class": "MEK inhibitor", 
        "expected_coverage": "AlphaMissense covered"
    },
    {
        "name": "TP53 R248W → Proteasome inhibitor",
        "mutation": {"gene": "TP53", "hgvs_p": "R248W", "chrom": "17", "pos": 7577120, "ref": "C", "alt": "T"},
        "drug_class": "proteasome inhibitor",
        "expected_coverage": "AlphaMissense covered"
    }
]

async def test_guidance(case: Dict[str, Any], api_base: str = "http://127.0.0.1:8000") -> Dict[str, Any]:
    """Test a single guidance case"""
    payload = {
        "disease": "multiple_myeloma",
        "drug_or_class": case["drug_class"],
        "mutations": [case["mutation"]],
        "api_base": api_base
    }
    
    async with httpx.AsyncClient(timeout=120.0) as client:
        response = await client.post(f"{api_base}/api/guidance/chemo", json=payload)
        if response.status_code != 200:
            return {
                "case": case["name"],
                "error": f"HTTP {response.status_code}: {response.text}",
                "success": False
            }
        
        data = response.json()
        return {
            "case": case["name"],
            "drug_class": case["drug_class"],
            "efficacy_score": data.get("efficacy_score"),
            "confidence": data.get("confidence"),
            "tier": data.get("tier"),
            "evidence_tier": data.get("evidence_tier"),
            "insights": data.get("insights", {}),
            "provenance": data.get("provenance", {}),
            "expected_coverage": case["expected_coverage"],
            "success": True
        }

async def test_fusion_score(mutation: Dict[str, Any], fusion_url: str = "https://crispro--fusion-engine-fusionengine-api.modal.run") -> Dict[str, Any]:
    """Test fusion engine for a mutation to check AM coverage"""
    variant_key = f"chr{mutation['chrom']}:{mutation['pos']}:{mutation['ref']}:{mutation['alt']}"
    
    async with httpx.AsyncClient(timeout=60.0) as client:
        response = await client.post(f"{fusion_url}/score_variants", json={"variants": [variant_key]})
        if response.status_code != 200:
            return {"variant": variant_key, "fusion_error": response.text}
        
        data = response.json()
        return {
            "variant": variant_key,
            "fusion_result": data.get("results", [{}])[0] if data.get("results") else {}
        }

async def main():
    print("🧬 MM Missense Guidance Test Suite")
    print("=" * 60)
    
    # Test all guidance cases
    guidance_results = []
    for case in TEST_CASES:
        print(f"\n🔬 Testing: {case['name']}")
        result = await test_guidance(case)
        guidance_results.append(result)
        
        if result["success"]:
            print(f"  ✅ Efficacy: {result['efficacy_score']:.3f}")
            print(f"  ✅ Confidence: {result['confidence']:.3f}")
            print(f"  ✅ Tier: {result['tier']} ({result['evidence_tier']})")
            print(f"  ✅ Insights: F={result['insights'].get('functionality', 0):.2f}, C={result['insights'].get('chromatin', 0):.2f}, E={result['insights'].get('essentiality', 0):.2f}, R={result['insights'].get('regulatory', 0):.2f}")
        else:
            print(f"  ❌ Error: {result['error']}")
    
    # Test fusion engine coverage
    print(f"\n🔗 Testing Fusion Engine Coverage")
    print("-" * 40)
    fusion_results = []
    for case in TEST_CASES:
        print(f"\n🧪 Fusion test: {case['name']}")
        fusion_result = await test_fusion_score(case["mutation"])
        fusion_results.append(fusion_result)
        
        if "fusion_error" not in fusion_result:
            fusion_data = fusion_result["fusion_result"]
            am_score = fusion_data.get("alphamissense_score", -999)
            evo_score = fusion_data.get("evo_score", -999)
            fused_score = fusion_data.get("fused_score", -999)
            
            print(f"  🎯 AlphaMissense: {am_score}")
            print(f"  🎯 Evo: {evo_score}")
            print(f"  🎯 Fused: {fused_score}")
            
            if am_score > -900:
                print(f"  ✅ AlphaMissense coverage: YES")
            else:
                print(f"  ❌ AlphaMissense coverage: NO")
        else:
            print(f"  ❌ Fusion error: {fusion_result['fusion_error']}")
    
    # Summary table
    print(f"\n📊 SUMMARY TABLE")
    print("=" * 80)
    print(f"{'Case':<30} {'Efficacy':<10} {'Confidence':<12} {'Tier':<8} {'AM Coverage':<12}")
    print("-" * 80)
    
    for i, result in enumerate(guidance_results):
        if result["success"]:
            fusion_data = fusion_results[i].get("fusion_result", {})
            am_covered = "YES" if fusion_data.get("alphamissense_score", -999) > -900 else "NO"
            
            print(f"{result['case']:<30} {result['efficacy_score']:<10.3f} {result['confidence']:<12.3f} {result['tier']:<8} {am_covered:<12}")
        else:
            print(f"{result['case']:<30} {'ERROR':<10} {'ERROR':<12} {'ERROR':<8} {'ERROR':<12}")

if __name__ == "__main__":
    asyncio.run(main())
