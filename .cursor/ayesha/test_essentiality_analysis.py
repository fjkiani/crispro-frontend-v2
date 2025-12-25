#!/usr/bin/env python3
"""
Test Essentiality Analysis for Ayesha (MBD4 + TP53)
Tests the essentiality endpoint and displays results
"""

import json
import httpx
import asyncio
from typing import Dict, Any

API_BASE = "http://127.0.0.1:8000"

async def test_gene_essentiality(gene: str, variants: list, model_id: str = "evo2_1b") -> Dict[str, Any]:
    """Test essentiality endpoint for a gene"""
    url = f"{API_BASE}/api/insights/predict_gene_essentiality"
    
    payload = {
        "gene": gene,
        "variants": variants,
        "model_id": model_id
    }
    
    print(f"\n{'='*60}")
    print(f"Testing Essentiality: {gene}")
    print(f"{'='*60}")
    print(f"Payload: {json.dumps(payload, indent=2)}")
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.post(url, json=payload)
            
            if response.status_code == 200:
                result = response.json()
                print(f"\n✅ SUCCESS")
                print(f"Essentiality Score: {result.get('essentiality_score', 'N/A')}")
                print(f"Flags: {result.get('flags', {})}")
                print(f"Rationale: {result.get('rationale', 'N/A')}")
                print(f"Confidence: {result.get('confidence', 'N/A')}")
                
                if 'provenance' in result:
                    print(f"\nProvenance:")
                    print(f"  Model: {result['provenance'].get('model_id', 'N/A')}")
                    if 'calibration' in result['provenance']:
                        calib = result['provenance']['calibration']
                        print(f"  Calibration: {calib}")
                
                return result
            else:
                print(f"\n❌ ERROR: Status {response.status_code}")
                print(f"Response: {response.text}")
                return {"error": response.text, "status_code": response.status_code}
                
    except Exception as e:
        print(f"\n❌ EXCEPTION: {str(e)}")
        return {"error": str(e)}

async def test_synthetic_lethality(mutations: list, disease: str = "ovarian_cancer") -> Dict[str, Any]:
    """Test synthetic lethality endpoint"""
    url = f"{API_BASE}/api/guidance/synthetic_lethality"
    
    payload = {
        "mutations": mutations,
        "disease": disease,
        "model_id": "evo2_1b"
    }
    
    print(f"\n{'='*60}")
    print(f"Testing Synthetic Lethality (Combined MBD4 + TP53)")
    print(f"{'='*60}")
    print(f"Payload: {json.dumps(payload, indent=2)}")
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.post(url, json=payload)
            
            if response.status_code == 200:
                result = response.json()
                print(f"\n✅ SUCCESS")
                print(f"Result: {json.dumps(result, indent=2)}")
                return result
            else:
                print(f"\n❌ ERROR: Status {response.status_code}")
                print(f"Response: {response.text}")
                return {"error": response.text, "status_code": response.status_code}
                
    except Exception as e:
        print(f"\n❌ EXCEPTION: {str(e)}")
        return {"error": str(e)}

async def main():
    """Run all essentiality tests"""
    
    print("\n" + "="*60)
    print("AYESHA ESSENTIALITY ANALYSIS TEST")
    print("="*60)
    
    # MBD4 variant (frameshift)
    mbd4_variant = {
        "gene": "MBD4",
        "chrom": "3",
        "pos": 129430456,
        "ref": "A",
        "alt": "",
        "consequence": "frameshift_variant",
        "hgvs_p": "p.Ile413Serfs*2"
    }
    
    # TP53 variant (hotspot)
    tp53_variant = {
        "gene": "TP53",
        "chrom": "17",
        "pos": 7577120,
        "ref": "G",
        "alt": "A",
        "consequence": "missense_variant",
        "hgvs_p": "p.Arg175His"
    }
    
    # Test MBD4
    mbd4_result = await test_gene_essentiality(
        gene="MBD4",
        variants=[mbd4_variant],
        model_id="evo2_1b"
    )
    
    # Test TP53
    tp53_result = await test_gene_essentiality(
        gene="TP53",
        variants=[tp53_variant],
        model_id="evo2_1b"
    )
    
    # Test Synthetic Lethality (combined)
    synthetic_lethality_result = await test_synthetic_lethality(
        mutations=[mbd4_variant, tp53_variant],
        disease="ovarian_cancer"
    )
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"\nMBD4 Essentiality Score: {mbd4_result.get('essentiality_score', 'N/A')}")
    print(f"TP53 Essentiality Score: {tp53_result.get('essentiality_score', 'N/A')}")
    
    if 'essentiality_score' in mbd4_result and 'essentiality_score' in tp53_result:
        mbd4_score = mbd4_result['essentiality_score']
        tp53_score = tp53_result['essentiality_score']
        
        print(f"\n📊 Analysis:")
        print(f"  - MBD4 frameshift → Essentiality: {mbd4_score}")
        print(f"  - TP53 hotspot → Essentiality: {tp53_score}")
        
        if mbd4_score >= 0.7:
            print(f"  ✅ MBD4 high essentiality → BER pathway non-functional")
        if tp53_score >= 0.7:
            print(f"  ✅ TP53 high essentiality → Checkpoint pathway bypassed")
        
        if mbd4_score >= 0.7 and tp53_score >= 0.7:
            print(f"\n🎯 Combined Effect:")
            print(f"  - Double-hit essentiality → Synthetic lethality expected")
            print(f"  - PARP sensitivity: HIGH (HR pathway becomes essential)")
            print(f"  - ATR sensitivity: HIGH (backup checkpoint essential)")
    
    print(f"\n{'='*60}\n")

if __name__ == "__main__":
    asyncio.run(main())




