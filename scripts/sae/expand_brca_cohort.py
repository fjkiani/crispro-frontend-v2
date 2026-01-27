#!/usr/bin/env python3
"""
Expand BRCA Cohort - Step 2
===========================

Mission: Extract SAE features for additional TCGA-BRCA patients
Objective: Increase sample size from 173 to 300+ patients (50+ events)
Timeline: 8-12 hours (SAE extraction is slow)
Status: EXECUTION READY
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
import asyncio
import httpx
import os

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OUTPUT_DIR = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
EXPANDED_CHECKPOINT = OUTPUT_DIR / "BRCA_TCGA_TRUE_SAE_cohort_expanded.json"

# API Configuration
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")
SAE_FEATURES_ENDPOINT = f"{BACKEND_URL}/api/sae/extract_features"

# Import extraction utilities
sys.path.insert(0, str(PROJECT_ROOT))
try:
    from scripts.sae.extract_sae_features_cohort import (
        extract_sae_features_for_variant,
        aggregate_patient_sae_features
    )
except ImportError:
    print("⚠️  Could not import extraction utilities")
    print("   Will use API calls directly")


async def extract_sae_for_variant_api(
    client: httpx.AsyncClient,
    variant: Dict[str, Any]
) -> Optional[Dict[str, Any]]:
    """Extract SAE features for a variant via API."""
    try:
        response = await client.post(
            SAE_FEATURES_ENDPOINT,
            json={
                "variant": variant,
                "model_id": "evo2_7b",
                "return_top_k": 64
            },
            timeout=180.0
        )
        
        if response.status_code == 200:
            return response.json()
        else:
            return None
    except Exception as e:
        print(f"      ⚠️  Error extracting SAE for variant: {e}")
        return None


async def expand_brca_cohort() -> Dict[str, Any]:
    """Expand BRCA cohort by extracting SAE features for additional patients."""
    print("=" * 80)
    print("🚀 STEP 2: EXPAND BRCA COHORT")
    print("=" * 80)
    print()
    
    # Step 1: Load current checkpoint
    print("📥 Step 1: Loading current BRCA checkpoint...")
    if not BRCA_CHECKPOINT.exists():
        print(f"❌ BRCA checkpoint not found: {BRCA_CHECKPOINT}")
        sys.exit(1)
    
    with open(BRCA_CHECKPOINT, 'r') as f:
        current_data = json.load(f)
    
    current_patients = current_data.get("data", {})
    print(f"   ✅ Current patients: {len(current_patients)}")
    
    # Step 2: Load recurrence labels
    print()
    print("📥 Step 2: Loading recurrence labels...")
    if not RECURRENCE_LABELS.exists():
        print(f"❌ Recurrence labels not found: {RECURRENCE_LABELS}")
        sys.exit(1)
    
    with open(RECURRENCE_LABELS, 'r') as f:
        labels_data = json.load(f)
    
    all_outcomes = labels_data.get("outcomes", {})
    print(f"   ✅ Total patients with labels: {len(all_outcomes)}")
    
    # Step 3: Identify patients to add
    print()
    print("🔍 Step 3: Identifying patients to add...")
    
    current_patient_ids = set(current_patients.keys())
    all_patient_ids = set(all_outcomes.keys())
    patients_to_add = all_patient_ids - current_patient_ids
    
    print(f"   Current patients: {len(current_patient_ids)}")
    print(f"   Total with labels: {len(all_patient_ids)}")
    print(f"   Patients to add: {len(patients_to_add)}")
    
    if len(patients_to_add) == 0:
        print("   ✅ All patients already extracted!")
        return {"status": "complete", "n_patients": len(current_patients)}
    
    # Step 4: Get mutations for new patients (from cBioPortal)
    print()
    print("📥 Step 4: Fetching mutations for new patients...")
    print("   ⚠️  NOTE: This requires cBioPortal access and mutation extraction")
    print("   ⚠️  For now, we'll create a plan for expansion")
    
    # TODO: Extract mutations for new patients
    # This would require:
    # 1. cBioPortal API access
    # 2. Extract mutations for each new patient
    # 3. Extract SAE features for each variant
    # 4. Aggregate per patient
    
    expansion_plan = {
        "expansion_date": datetime.now().isoformat(),
        "current_patients": len(current_patients),
        "total_patients_with_labels": len(all_patient_ids),
        "patients_to_add": len(patients_to_add),
        "patient_ids_to_add": list(patients_to_add)[:50],  # First 50 for now
        "estimated_variants_per_patient": 30,
        "estimated_total_variants": len(patients_to_add) * 30,
        "estimated_time_hours": len(patients_to_add) * 0.5,  # ~30 min per patient
        "status": "plan_created",
        "next_steps": [
            "1. Extract mutations for new patients from cBioPortal",
            "2. Extract SAE features for each variant (via API)",
            "3. Aggregate features per patient",
            "4. Merge with existing checkpoint",
            "5. Re-run validation with expanded cohort"
        ]
    }
    
    # Save expansion plan
    plan_file = OUTPUT_DIR.parent / "brca_cohort_expansion_plan.json"
    with open(plan_file, 'w') as f:
        json.dump(expansion_plan, f, indent=2)
    
    print(f"   💾 Expansion plan saved to: {plan_file}")
    print()
    print("=" * 80)
    print("📊 EXPANSION PLAN")
    print("=" * 80)
    print(f"Current patients: {len(current_patients)}")
    print(f"Patients to add: {len(patients_to_add)}")
    print(f"Estimated variants: {len(patients_to_add) * 30}")
    print(f"Estimated time: {len(patients_to_add) * 0.5:.1f} hours")
    print()
    print("⚠️  NOTE: Full expansion requires:")
    print("   1. cBioPortal mutation extraction")
    print("   2. SAE feature extraction (slow, ~30 min/patient)")
    print("   3. This is a long-running process (8-12 hours)")
    print()
    print("💡 RECOMMENDATION: Run expansion in batches")
    print("   - Batch 1: 50 patients (~4 hours)")
    print("   - Batch 2: 50 patients (~4 hours)")
    print("   - Batch 3: Remaining patients")
    print("=" * 80)
    
    return expansion_plan


if __name__ == "__main__":
    try:
        result = asyncio.run(expand_brca_cohort())
        print("\n✅ STEP 2 COMPLETE: Expansion plan created")
        print("\n📋 NEXT: Execute expansion (requires mutation extraction + SAE extraction)")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ EXPANSION PLANNING FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
