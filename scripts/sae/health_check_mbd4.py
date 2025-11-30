#!/usr/bin/env python3
"""
🔍 MBD4+TP53 Analysis Health Check

Purpose: Verify MBD4+TP53 analysis can run end-to-end with proxy SAE
Agent: Zo (Lead Commander)
Date: January 20, 2025
"""

import asyncio
import httpx
import json
import sys

async def check_mbd4_analysis():
    """Verify MBD4+TP53 analysis pipeline works."""
    
    print("=" * 80)
    print("HEALTH CHECK: MBD4+TP53 Analysis Pipeline")
    print("=" * 80)
    
    base_url = "http://localhost:8000"
    
    # MBD4+TP53 test case
    mutations = [
        {"gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2", "chrom": "3", "pos": 129430456, "ref": "A", "alt": "", "build": "GRCh38"},
        {"gene": "TP53", "hgvs_p": "p.R175H", "chrom": "17", "pos": 7577120, "ref": "G", "alt": "A", "build": "GRCh38"}
    ]
    
    tumor_context = {
        "disease": "ovarian_cancer",
        "hrd_score": 0.75,
        "tmb_score": 25.0,
        "msi_status": "MSS"
    }
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Step 1: Efficacy prediction (gets pathway scores)
            print("\n[1/3] Calling /api/efficacy/predict...")
            efficacy_response = await client.post(
                f"{base_url}/api/efficacy/predict",
                json={
                    "model_id": "evo2_1b",
                    "mutations": mutations,
                    "tumor_context": tumor_context,
                    "options": {"adaptive": True}
                }
            )
            
            if efficacy_response.status_code != 200:
                print(f"❌ FAIL: Efficacy returned {efficacy_response.status_code}")
                return False
            
            efficacy_data = efficacy_response.json()
            
            # Extract pathway scores (proxy SAE source)
            pathway_scores = {}
            if "provenance" in efficacy_data and "confidence_breakdown" in efficacy_data["provenance"]:
                pathway_disruption = efficacy_data["provenance"]["confidence_breakdown"].get("pathway_disruption", {})
                pathway_scores = {
                    "ddr": pathway_disruption.get("ddr", 0.0),
                    "mapk": pathway_disruption.get("mapk", 0.0),
                    "pi3k": pathway_disruption.get("pi3k", 0.0),
                    "vegf": pathway_disruption.get("vegf", 0.0),
                    "her2": pathway_disruption.get("her2", 0.0)
                }
            
            print(f"✅ Efficacy prediction complete")
            print(f"   Pathway scores: DDR={pathway_scores.get('ddr', 0):.2f}, MAPK={pathway_scores.get('mapk', 0):.2f}")
            
            # Step 2: SAE features computation
            print("\n[2/3] Calling /api/sae/compute_features...")
            
            # Extract insights bundle from efficacy response
            insights_bundle = {}
            if "drugs" in efficacy_data and len(efficacy_data["drugs"]) > 0:
                first_drug = efficacy_data["drugs"][0]
                if "insights" in first_drug:
                    insights_bundle = first_drug["insights"]
            
            sae_response = await client.post(
                f"{base_url}/api/sae/compute_features",
                json={
                    "pathway_scores": pathway_scores,
                    "insights_bundle": insights_bundle,
                    "tumor_context": tumor_context
                }
            )
            
            if sae_response.status_code != 200:
                print(f"❌ FAIL: SAE features returned {sae_response.status_code}")
                return False
            
            sae_data = sae_response.json()
            
            print(f"✅ SAE features computed")
            if "dna_repair_capacity" in sae_data:
                print(f"   DNA repair capacity: {sae_data['dna_repair_capacity']:.3f}")
            if "mechanism_vector" in sae_data:
                print(f"   Mechanism vector: {sae_data['mechanism_vector']}")
            
            # Step 3: Trial matching (optional, may not be available)
            print("\n[3/3] Checking trial matching endpoint...")
            try:
                trial_response = await client.post(
                    f"{base_url}/api/trials/agent/search",
                    json={
                        "mutations": mutations,
                        "disease": "ovarian_cancer",
                        "tumor_context": tumor_context
                    },
                    timeout=30.0
                )
                
                if trial_response.status_code == 200:
                    print(f"✅ Trial matching available")
                else:
                    print(f"⚠️  Trial matching returned {trial_response.status_code} (optional)")
            except Exception as e:
                print(f"⚠️  Trial matching not available: {e} (optional)")
            
            print("\n✅ MBD4+TP53 analysis pipeline operational")
            return True
            
    except httpx.TimeoutException:
        print("❌ FAIL: Request timed out")
        return False
    except Exception as e:
        print(f"❌ FAIL: Analysis check failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = asyncio.run(check_mbd4_analysis())
    sys.exit(0 if success else 1)

