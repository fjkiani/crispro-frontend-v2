#!/usr/bin/env python3
"""
Validate Ayesha Drug Confidence Scores - PRODUCTION API
========================================================

Purpose: Replace ALL assumed confidence scores with actual API outputs
API: https://crispro--evo-service-evoservice1b-api-1b.modal.run
Timeline: 30 minutes
Status: EXECUTION READY

Validates:
- Olaparib: assumed 70% → actual ??%
- Niraparib: assumed 65% → actual ??%
- Pembrolizumab: assumed 65% → actual ??%
- Bevacizumab: assumed 60% → actual ??%
- Carboplatin: assumed ??% → actual ??%
- Imatinib: assumed 35% → actual ??%
"""

import asyncio
import httpx
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any

# Configuration - PRODUCTION API
API_BASE_URL = "https://crispro--evo-service-evoservice1b-api-1b.modal.run"
OUTPUT_DIR = Path("results/ayesha_validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Ayesha's Profile
AYESHA_MUTATIONS = [
    {
        "gene": "MBD4",
        "hgvs_p": "p.K431Nfs*54",
        "chrom": "3",
        "pos": 129149435,
        "ref": "A",
        "alt": "",
        "consequence": "frameshift_variant"
    },
    {
        "gene": "TP53",
        "hgvs_p": "p.R273H",
        "chrom": "17",
        "pos": 7577120,
        "ref": "G",
        "alt": "A",
        "consequence": "missense_variant"
    },
    {
        "gene": "PDGFRA",
        "hgvs_p": "p.S755P",
        "chrom": "4",
        "pos": 54280422,
        "ref": "T",
        "alt": "C",
        "consequence": "missense_variant"
    }
]

# Drugs to validate
DRUGS_TO_VALIDATE = [
    {"name": "olaparib", "assumed": 0.70, "rationale": "MBD4→PARP via SL"},
    {"name": "niraparib", "assumed": 0.65, "rationale": "MBD4→PARP via SL (alt)"},
    {"name": "pembrolizumab", "assumed": 0.65, "rationale": "TMB-HIGH + PD-L1+"},
    {"name": "bevacizumab", "assumed": 0.60, "rationale": "Stage IVB + PDGFRA"},
    {"name": "carboplatin", "assumed": None, "rationale": "Standard of care"},
    {"name": "imatinib", "assumed": 0.35, "rationale": "PDGFRA VUS"},
]


async def call_efficacy_api(client: httpx.AsyncClient, drug: str) -> Dict[str, Any]:
    """Call production efficacy API for single drug"""
    
    payload = {
        "model_id": "evo2_1b",
        "mutations": AYESHA_MUTATIONS,
        "disease": "ovarian_cancer",
        "drugs": [drug]
    }
    
    print(f"\n🔍 Calling API for {drug}...")
    
    try:
        response = await client.post(
            f"{API_BASE_URL}/api/efficacy/predict",
            json=payload,
            timeout=120.0
        )
        response.raise_for_status()
        result = response.json()
        
        # Extract drug data from response (it's in "drugs" array)
        drugs_list = result.get("drugs", [])
        drug_data = next((d for d in drugs_list if d.get("name") == drug), None)
        
        if drug_data:
            confidence = drug_data.get("confidence")
            efficacy = drug_data.get("efficacy_score")
            tier = drug_data.get("evidence_tier")
            badges = drug_data.get("badges", [])
            
            # Check provenance for cache hits
            provenance = result.get("provenance", {})
            cache_status = provenance.get("cache", "unknown")
            
            print(f"   ✅ Confidence: {confidence:.3f}" if confidence else "   ⚠️  No confidence")
            print(f"   Cache: {cache_status}")
            
            return {
                "drug": drug,
                "confidence": confidence,
                "efficacy_score": efficacy,
                "evidence_tier": tier,
                "badges": badges,
                "cache_status": cache_status,
                "full_response": result
            }
        else:
            print(f"   ❌ Drug not found in response")
            return {"drug": drug, "error": "Drug not in response"}
            
    except Exception as e:
        print(f"   ❌ Error: {str(e)}")
        return {"drug": drug, "error": str(e)}


async def validate_all_drugs():
    """Validate all Ayesha drugs"""
    
    print("=" * 80)
    print("AYESHA DRUG CONFIDENCE VALIDATION - PRODUCTION API")
    print("=" * 80)
    print(f"Started: {datetime.now().isoformat()}")
    print(f"API: {API_BASE_URL}")
    print(f"Mutations: MBD4, TP53, PDGFRA")
    print("=" * 80)
    
    results = {
        "timestamp": datetime.now().isoformat(),
        "api_endpoint": API_BASE_URL,
        "mutations": AYESHA_MUTATIONS,
        "validations": []
    }
    
    async with httpx.AsyncClient() as client:
        for drug_info in DRUGS_TO_VALIDATE:
            drug_name = drug_info["name"]
            assumed = drug_info["assumed"]
            
            print(f"\n{'='*80}")
            print(f"DRUG: {drug_name.upper()}")
            print(f"Assumed: {assumed if assumed else 'N/A'}")
            print(f"Rationale: {drug_info['rationale']}")
            print(f"{'='*80}")
            
            api_result = await call_efficacy_api(client, drug_name)
            
            validation = {
                "drug": drug_name,
                "assumed_confidence": assumed,
                "actual_confidence": api_result.get("confidence"),
                "efficacy_score": api_result.get("efficacy_score"),
                "evidence_tier": api_result.get("evidence_tier"),
                "badges": api_result.get("badges", []),
                "cache_status": api_result.get("cache_status"),
                "rationale": drug_info["rationale"]
            }
            
            # Calculate validation status
            if "error" in api_result:
                validation["status"] = f"❌ ERROR: {api_result['error']}"
                validation["difference"] = None
            elif validation["actual_confidence"] is None:
                validation["status"] = "⚠️  NO CONFIDENCE RETURNED"
                validation["difference"] = None
            elif assumed is None:
                validation["status"] = "✅ BASELINE (no assumption)"
                validation["difference"] = None
            else:
                diff = validation["actual_confidence"] - assumed
                validation["difference"] = round(diff, 3)
                
                if abs(diff) <= 0.05:
                    validation["status"] = "✅ VALIDATED (within 5%)"
                elif abs(diff) <= 0.10:
                    validation["status"] = "⚠️  CLOSE (within 10%)"
                else:
                    validation["status"] = f"❌ MISMATCH ({diff:+.1%})"
                
                print(f"\n📊 RESULT:")
                print(f"   Actual: {validation['actual_confidence']:.3f}")
                print(f"   Difference: {diff:+.3f} ({diff:+.1%})")
                print(f"   Status: {validation['status']}")
            
            results["validations"].append(validation)
    
    # Summary
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    
    validated = sum(1 for v in results["validations"] if "✅ VALIDATED" in v["status"])
    close = sum(1 for v in results["validations"] if "⚠️  CLOSE" in v["status"])
    mismatch = sum(1 for v in results["validations"] if "❌ MISMATCH" in v["status"])
    
    print(f"\n✅ Validated (within 5%): {validated}/{len(DRUGS_TO_VALIDATE)}")
    print(f"⚠️  Close (within 10%): {close}/{len(DRUGS_TO_VALIDATE)}")
    print(f"❌ Mismatch (>10%): {mismatch}/{len(DRUGS_TO_VALIDATE)}")
    
    # Detailed table
    print("\n" + "=" * 80)
    print("DETAILED RESULTS")
    print("=" * 80)
    print(f"{'Drug':<15} {'Assumed':<10} {'Actual':<10} {'Diff':<10} {'Status':<30}")
    print("-" * 80)
    
    for v in results["validations"]:
        assumed_str = f"{v['assumed_confidence']:.2f}" if v['assumed_confidence'] else "N/A"
        actual_str = f"{v['actual_confidence']:.2f}" if v['actual_confidence'] else "N/A"
        diff_str = f"{v['difference']:+.2f}" if v['difference'] else "N/A"
        print(f"{v['drug']:<15} {assumed_str:<10} {actual_str:<10} {diff_str:<10} {v['status']:<30}")
    
    # Save results
    output_file = OUTPUT_DIR / f"ayesha_validation_PRODUCTION_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 80)
    print(f"📁 Results saved to: {output_file}")
    print(f"Completed: {datetime.now().isoformat()}")
    print("=" * 80)
    
    return results


if __name__ == "__main__":
    asyncio.run(validate_all_drugs())
