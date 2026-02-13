#!/usr/bin/env python3
"""
Ayesha Confidence Validation - LOCAL BACKEND
=============================================
Using localhost:8000 (the backend that has /api/efficacy/predict)
Production Evo2 service is for scoring only, not efficacy prediction
"""

import asyncio
import httpx
import json
from datetime import datetime
from pathlib import Path

API_BASE_URL = "http://127.0.0.1:8000"
OUTPUT_DIR = Path("results/ayesha_validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

AYESHA_MUTATIONS = [
    {"gene": "MBD4", "hgvs_p": "p.K431Nfs*54", "chrom": "3", "pos": 129149435, "ref": "A", "alt": "", "consequence": "frameshift_variant"},
    {"gene": "TP53", "hgvs_p": "p.R273H", "chrom": "17", "pos": 7577120, "ref": "G", "alt": "A", "consequence": "missense_variant"},
    {"gene": "PDGFRA", "hgvs_p": "p.S755P", "chrom": "4", "pos": 54280422, "ref": "T", "alt": "C", "consequence": "missense_variant"}
]

DRUGS = [
    {"name": "olaparib", "assumed": 0.70},
    {"name": "niraparib", "assumed": 0.65},
    {"name": "pembrolizumab", "assumed": 0.65},
    {"name": "bevacizumab", "assumed": 0.60},
    {"name": "carboplatin", "assumed": None},
    {"name": "imatinib", "assumed": 0.35}
]

async def validate():
    print("="*80)
    print("AYESHA CONFIDENCE VALIDATION - FINAL")
    print("="*80)
    print(f"API: {API_BASE_URL}")
    print(f"Mutations: MBD4, TP53, PDGFRA")
    print("="*80)
    
    results = {"timestamp": datetime.now().isoformat(), "validations": []}
    
    async with httpx.AsyncClient() as client:
        for drug_info in DRUGS:
            drug = drug_info["name"]
            assumed = drug_info["assumed"]
            
            print(f"\n{'='*80}")
            print(f"VALIDATING: {drug.upper()}")
            print(f"Assumed: {assumed if assumed else 'N/A'}")
            
            try:
                response = await client.post(
                    f"{API_BASE_URL}/api/efficacy/predict",
                    json={"model_id": "evo2_1b", "mutations": AYESHA_MUTATIONS, "disease": "ovarian_cancer", "drugs": [drug]},
                    timeout=120.0
                )
                response.raise_for_status()
                result = response.json()
                
                drugs_list = result.get("drugs", [])
                drug_data = next((d for d in drugs_list if d.get("name") == drug), None)
                
                if drug_data:
                    actual = drug_data.get("confidence")
                    efficacy = drug_data.get("efficacy_score")
                    tier = drug_data.get("evidence_tier")
                    badges = drug_data.get("badges", [])
                    
                    print(f"✅ Actual: {actual:.3f}")
                    print(f"   Efficacy: {efficacy:.3f}")
                    print(f"   Tier: {tier}")
                    print(f"   Badges: {badges}")
                    
                    validation = {
                        "drug": drug,
                        "assumed": assumed,
                        "actual": actual,
                        "efficacy_score": efficacy,
                        "tier": tier,
                        "badges": badges
                    }
                    
                    if assumed:
                        diff = actual - assumed
                        validation["difference"] = round(diff, 3)
                        validation["diff_pct"] = round(diff * 100, 1)
                        
                        if abs(diff) <= 0.05:
                            validation["status"] = "✅ VALIDATED"
                        elif abs(diff) <= 0.10:
                            validation["status"] = "⚠️  CLOSE"
                        else:
                            validation["status"] = "❌ MISMATCH"
                        
                        print(f"   Difference: {diff:+.3f} ({diff:+.1%})")
                        print(f"   Status: {validation['status']}")
                    else:
                        validation["status"] = "✅ BASELINE"
                    
                    results["validations"].append(validation)
                    
            except Exception as e:
                print(f"❌ Error: {e}")
                results["validations"].append({"drug": drug, "error": str(e)})
    
    # Summary
    print("\n" + "="*80)
    print("FINAL RESULTS")
    print("="*80)
    print(f"\n{'Drug':<15} {'Assumed':<10} {'Actual':<10} {'Diff':<10} {'Status':<20}")
    print("-"*80)
    
    for v in results["validations"]:
        if "error" not in v:
            assumed_str = f"{v['assumed']:.2f}" if v['assumed'] else "N/A"
            actual_str = f"{v['actual']:.2f}"
            diff_str = f"{v.get('difference', 0):+.2f}" if v.get('difference') else "N/A"
            print(f"{v['drug']:<15} {assumed_str:<10} {actual_str:<10} {diff_str:<10} {v['status']:<20}")
    
    # Save
    output_file = OUTPUT_DIR / f"FINAL_ayesha_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n📁 Saved: {output_file}")
    print("="*80)
    
    return results

if __name__ == "__main__":
    asyncio.run(validate())
