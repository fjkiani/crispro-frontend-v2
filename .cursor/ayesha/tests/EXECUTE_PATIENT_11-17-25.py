#!/usr/bin/env python3
"""
Execute Patient 11-17-25 Through AK Framework

This script demonstrates running a NEW patient through the existing Ayesha framework
WITHOUT any hardcoding. Uses the same services that work for Ayesha.

Patient Profile:
- Metastatic high-grade serous carcinoma (ovarian/peritoneal)
- PD-L1 positive (CPS 10)
- p53 mutant
- ER weakly positive (50%)
- MMR preserved
- HER-2 negative
- FOLR1 negative
- NTRK negative

Uses:
1. Complete Care v2 Orchestrator - Unified care plan
2. Research Intelligence - Answer mechanism questions
3. All existing Ayesha services (CA-125, Trials, Drug Ranking, etc.)
"""

import asyncio
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional

# Load .env files
from dotenv import load_dotenv
project_root = Path(__file__).parent.parent.parent.parent
backend_env = project_root / "oncology-coPilot" / "oncology-backend-minimal" / ".env"
root_env = project_root / ".env"

if backend_env.exists():
    load_dotenv(backend_env, override=False)
if root_env.exists():
    load_dotenv(root_env, override=True)

sys.path.insert(0, str(project_root / "oncology-coPilot" / "oncology-backend-minimal"))

import httpx

# Patient Profile (extracted from pathology report)
PATIENT_PROFILE = {
    "ca125_value": 0,  # Not available in pathology report
    "stage": "IVB",
    "treatment_line": "first-line",
    "germline_status": "negative",
    "has_ascites": True,  # Multiple nodules in omentum
    "has_peritoneal_disease": True,  # Omental biopsies
    "location_state": "NY",
    "ecog_status": None,  # Not provided
    "tumor_context": {
        "somatic_mutations": [
            {
                "gene": "TP53",
                "hgvs_p": "p.Arg175His",
                "variant_type": "missense",
                "pathway": "DDR",
                "hotspot": True
            }
        ],
        "pd_l1": {
            "cps": 10,
            "status": "POSITIVE"
        },
        "hrd": {
            "status": "UNKNOWN"
        },
        "tmb": {
            "classification": "UNKNOWN"
        },
        "msi": {
            "status": "MSS",
            "classification": "STABLE"
        },
        "er": {
            "status": "WEAKLY_POSITIVE",
            "percent": 50
        },
        "pr": {
            "status": "NEGATIVE",
            "percent": 1
        }
    }
}

# Research Intelligence Questions
RESEARCH_QUESTIONS = [
    {
        "priority": 1,
        "question": "What are the mechanisms of action for PD-L1 positive ovarian cancer immunotherapy?",
        "context": {
            "compound": "checkpoint inhibitor",
            "disease": "ovarian cancer",
            "biomarkers": {"PD-L1": "POSITIVE", "CPS": 10}
        }
    },
    {
        "priority": 2,
        "question": "Should patients with p53 mutant ovarian cancer receive PARP inhibitors even if MMR is preserved?",
        "context": {
            "compound": "PARP inhibitor",
            "disease": "ovarian cancer",
            "biomarkers": {"p53": "MUTANT", "MMR": "PRESERVED"}
        }
    },
    {
        "priority": 3,
        "question": "What are the treatment options for ER weakly positive ovarian cancer?",
        "context": {
            "compound": "hormone therapy",
            "disease": "ovarian cancer",
            "biomarkers": {"ER": "WEAKLY_POSITIVE", "percent": 50}
        }
    }
]

API_BASE = "http://localhost:8000"

async def call_complete_care_v2(profile: Dict[str, Any]) -> Dict[str, Any]:
    """Call Complete Care v2 Orchestrator."""
    print("\n" + "=" * 70)
    print("STEP 1: COMPLETE CARE v2 ORCHESTRATOR")
    print("=" * 70)
    
    request = {
        "ca125_value": profile["ca125_value"],
        "stage": profile["stage"],
        "treatment_line": profile["treatment_line"],
        "germline_status": profile["germline_status"],
        "has_ascites": profile["has_ascites"],
        "has_peritoneal_disease": profile["has_peritoneal_disease"],
        "location_state": profile["location_state"],
        "ecog_status": profile.get("ecog_status"),
        "tumor_context": profile["tumor_context"],
        "include_trials": True,
        "include_soc": True,
        "include_ca125": True,
        "include_wiwfm": True,
        "include_food": False,
        "include_resistance": True,
        "include_resistance_prediction": False,
        "max_trials": 10
    }
    
    try:
        async with httpx.AsyncClient(timeout=120.0) as client:
            response = await client.post(
                f"{API_BASE}/api/ayesha/complete_care_v2",
                json=request,
                headers={"Content-Type": "application/json"}
            )
            
            if response.status_code == 200:
                data = response.json()
                print(f"✅ Complete Care v2 successful")
                
                # Extract key metrics
                trials = data.get("trials", {}).get("trials", [])
                soc = data.get("soc_recommendation")
                ca125 = data.get("ca125_intelligence")
                wiwfm = data.get("wiwfm")
                resistance = data.get("resistance_playbook")
                
                print(f"  - Trials: {len(trials)} found")
                print(f"  - SOC: {'✅' if soc else '❌'}")
                print(f"  - CA-125: {'✅' if ca125 else '❌'}")
                print(f"  - Drug Ranking: {'✅' if wiwfm else '❌'}")
                print(f"  - Resistance Playbook: {'✅' if resistance else '❌'}")
                
                return data
            else:
                print(f"❌ Complete Care v2 failed: {response.status_code}")
                print(f"   {response.text[:200]}")
                return {}
    except Exception as e:
        print(f"❌ Complete Care v2 error: {e}")
        import traceback
        traceback.print_exc()
        return {}

async def call_research_intelligence(question_data: Dict[str, Any]) -> Dict[str, Any]:
    """Call Research Intelligence on a question."""
    from api.services.research_intelligence.orchestrator import ResearchIntelligenceOrchestrator
    
    orchestrator = ResearchIntelligenceOrchestrator()
    if not orchestrator.is_available():
        return {"error": "Orchestrator not available"}
    
    try:
        result = await orchestrator.research_question(
            question=question_data["question"],
            context=question_data["context"]
        )
        
        synthesized = result.get("synthesized_findings", {})
        mechanisms = synthesized.get("mechanisms", [])
        evidence_tier = synthesized.get("evidence_tier", "Unknown")
        confidence = synthesized.get("overall_confidence", 0.0)
        
        return {
            "question": question_data["question"],
            "mechanisms": mechanisms[:5],
            "evidence_tier": evidence_tier,
            "confidence": confidence,
            "pathways": result.get("moat_analysis", {}).get("pathways", [])
        }
    except Exception as e:
        return {"error": str(e)}

async def main():
    """Run complete patient analysis."""
    print("\n" + "=" * 70)
    print("PATIENT 11-17-25: COMPLETE FRAMEWORK EXECUTION")
    print("Metastatic High-Grade Serous Carcinoma (Ovarian/Peritoneal)")
    print("=" * 70)
    
    print(f"\n📋 Patient Profile:")
    print(f"  Disease: Ovarian Cancer (HGS)")
    print(f"  Stage: {PATIENT_PROFILE['stage']}")
    print(f"  Treatment Line: {PATIENT_PROFILE['treatment_line']}")
    print(f"  Key Biomarkers:")
    print(f"    - PD-L1: POSITIVE (CPS {PATIENT_PROFILE['tumor_context']['pd_l1']['cps']})")
    print(f"    - p53: MUTANT")
    print(f"    - ER: WEAKLY POSITIVE ({PATIENT_PROFILE['tumor_context']['er']['percent']}%)")
    print(f"    - MMR: PRESERVED")
    
    # Step 1: Complete Care v2
    care_plan = await call_complete_care_v2(PATIENT_PROFILE)
    
    # Step 2: Research Intelligence
    print("\n" + "=" * 70)
    print("STEP 2: RESEARCH INTELLIGENCE")
    print("=" * 70)
    
    research_results = {}
    for i, query_data in enumerate(RESEARCH_QUESTIONS, 1):
        print(f"\n[{i}/{len(RESEARCH_QUESTIONS)}] {query_data['question'][:60]}...")
        result = await call_research_intelligence(query_data)
        research_results[query_data['priority']] = result
        
        if "error" not in result:
            print(f"  ✅ Found {len(result.get('mechanisms', []))} mechanisms")
            print(f"     Evidence: {result.get('evidence_tier')}, Confidence: {result.get('confidence', 0):.2f}")
        else:
            print(f"  ❌ Error: {result['error']}")
        
        # Delay between queries
        if i < len(RESEARCH_QUESTIONS):
            await asyncio.sleep(3)
    
    # Save results
    output_file = Path(__file__).parent / "patient_11-17-25_results.json"
    output = {
        "patient_profile": PATIENT_PROFILE,
        "care_plan": care_plan,
        "research_intelligence": research_results,
        "timestamp": datetime.now().isoformat()
    }
    
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    if care_plan:
        trials_count = len(care_plan.get("trials", {}).get("trials", []))
        print(f"✅ Complete Care Plan Generated")
        print(f"   - Trials: {trials_count}")
        print(f"   - SOC: {'✅' if care_plan.get('soc_recommendation') else '❌'}")
        print(f"   - CA-125: {'✅' if care_plan.get('ca125_intelligence') else '❌'}")
        print(f"   - Drug Ranking: {'✅' if care_plan.get('wiwfm') else '❌'}")
        print(f"   - Resistance Playbook: {'✅' if care_plan.get('resistance_playbook') else '❌'}")
    
    total_mechanisms = sum(len(r.get("mechanisms", [])) for r in research_results.values() if "error" not in r)
    print(f"\n✅ Research Intelligence")
    print(f"   - Questions answered: {len([r for r in research_results.values() if 'error' not in r])}/{len(RESEARCH_QUESTIONS)}")
    print(f"   - Total mechanisms: {total_mechanisms}")
    
    print("\n✅ Patient analysis complete!")
    print("\nNext Steps:")
    print("  1. Review Complete Care Plan results")
    print("  2. Review Research Intelligence findings")
    print("  3. Order NGS tests (HRD, TMB, MSI) for full S/P/E scoring")
    print("  4. Order CA-125 baseline for monitoring")

if __name__ == "__main__":
    asyncio.run(main())

