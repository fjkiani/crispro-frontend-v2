#!/usr/bin/env python3
"""
Extract BRCA Recurrence Labels - Deliverable #1
================================================

Mission: Extract recurrence outcomes from TCGA-BRCA clinical data
Objective: Get DFS_STATUS (Disease-Free Survival) for 200 BRCA patients
Timeline: 4 hours
Status: EXECUTION READY

Based on: SAE_BRCA_NEXT_10_DELIVERABLES.md (Deliverable #1)
Corrections: Use DFS_STATUS as primary (not new_tumor_event), handle multiple field names
"""

import sys
from pathlib import Path
import json
from typing import Dict, List, Any, Optional
from datetime import datetime
import pandas as pd

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import existing extraction utilities
try:
    from oncology_coPilot.oncology_backend_minimal.scripts.benchmark.extract_cbioportal_trial_datasets import (
        extract_clinical_outcomes
    )
except ImportError:
    # Fallback: direct import
    sys.path.insert(0, str(PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal"))
    try:
        from scripts.benchmark.extract_cbioportal_trial_datasets import (
            extract_clinical_outcomes
        )
    except ImportError as e:
        print(f"❌ Error importing extraction utilities: {e}")
        print("   Trying alternative import path...")
        # Try using pyBioPortal directly
        pybioportal_path = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend" / "tests" / "pyBioPortal-master"
        if str(pybioportal_path) not in sys.path:
            sys.path.insert(0, str(pybioportal_path))
        
        try:
            from pybioportal import clinical_data as cd
        except ImportError:
            print("❌ Cannot import pyBioPortal. Please install or check path.")
            sys.exit(1)

# Configuration
STUDY_ID = "brca_tcga_pan_can_atlas_2018"  # TCGA-BRCA PanCan Atlas
OUTPUT_DIR = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "brca_recurrence_labels.json"

# BRCA checkpoint file (to match patient IDs)
BRCA_CHECKPOINT = OUTPUT_DIR / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"


def parse_dfs_status(dfs_status: Optional[str]) -> Optional[bool]:
    """
    Parse DFS_STATUS to binary recurrence outcome.
    
    Formats:
    - "0:DiseaseFree" → recurrence = False
    - "1:Recurred/Progressed" → recurrence = True
    - "1:Recurred" → recurrence = True
    - "0:CENSORED" → recurrence = False (censored = no recurrence observed)
    """
    if not dfs_status or not isinstance(dfs_status, str):
        return None
    
    dfs_lower = dfs_status.lower()
    
    # Recurrence = True
    if dfs_status.startswith("1:") or any(x in dfs_lower for x in ["recurred", "progressed", "recurrence"]):
        return True
    
    # No recurrence = False
    if dfs_status.startswith("0:") or any(x in dfs_lower for x in ["diseasefree", "disease-free", "censored"]):
        return False
    
    return None


def parse_months(value) -> Optional[float]:
    """Parse months field, handling NaN strings and None."""
    if value is None:
        return None
    if isinstance(value, str):
        if value.lower() in ["nan", "na", "", "null", "none"]:
            return None
        try:
            return float(value)
        except ValueError:
            return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def extract_brca_recurrence_labels() -> Dict[str, Any]:
    """
    Extract recurrence outcomes from TCGA-BRCA clinical data.
    
    CORRECTED: Use DFS_STATUS as primary (not new_tumor_event)
    Handles multiple field name variations (defensive)
    """
    print("=" * 80)
    print("🔬 EXTRACTING BRCA RECURRENCE LABELS - Deliverable #1")
    print("=" * 80)
    print(f"Study ID: {STUDY_ID}")
    print(f"Output: {OUTPUT_FILE}")
    print()
    
    # Step 1: Load BRCA checkpoint to get patient IDs
    print("📥 Step 1: Loading BRCA checkpoint to get patient IDs...")
    if not BRCA_CHECKPOINT.exists():
        print(f"⚠️  BRCA checkpoint not found: {BRCA_CHECKPOINT}")
        print("   Will extract for all patients in study")
        target_patient_ids = None
    else:
        with open(BRCA_CHECKPOINT, 'r') as f:
            brca_data = json.load(f)
        
        # Extract patient IDs from checkpoint
        if "data" in brca_data and isinstance(brca_data["data"], dict):
            target_patient_ids = set(brca_data["data"].keys())
            print(f"   ✅ Found {len(target_patient_ids)} patients in BRCA checkpoint")
        else:
            print("   ⚠️  Unexpected checkpoint structure, extracting for all patients")
            target_patient_ids = None
    
    # Step 2: Extract clinical outcomes
    print()
    print("📊 Step 2: Extracting clinical outcomes from cBioPortal...")
    try:
        clinical_df = extract_clinical_outcomes(STUDY_ID)
    except Exception as e:
        print(f"❌ Error extracting clinical outcomes: {e}")
        print(f"   Trying alternative method...")
        
        # Fallback: Use pyBioPortal directly
        try:
            from pybioportal import clinical_data as cd
            clinical_df = cd.get_all_clinical_data_in_study(STUDY_ID)
            if not isinstance(clinical_df, pd.DataFrame):
                raise ValueError("Expected DataFrame from pyBioPortal")
        except Exception as e2:
            print(f"❌ Alternative method also failed: {e2}")
            raise
    
    if clinical_df.empty:
        raise ValueError(f"No clinical data found for {STUDY_ID}")
    
    print(f"   ✅ Retrieved {len(clinical_df)} patients from cBioPortal")
    
    # Step 3: Extract recurrence labels
    print()
    print("🔍 Step 3: Extracting recurrence labels...")
    outcomes = {}
    stats = {
        "total_patients": 0,
        "with_recurrence_label": 0,
        "recurrence_true": 0,
        "recurrence_false": 0,
        "with_dfs_months": 0,
        "with_vital_status": 0,
        "matched_to_checkpoint": 0
    }
    
    for _, row in clinical_df.iterrows():
        patient_id = str(row.get("patientId", ""))
        if not patient_id:
            continue
        
        stats["total_patients"] += 1
        
        # Check if patient is in checkpoint (if we have target list)
        if target_patient_ids and patient_id not in target_patient_ids:
            continue
        
        if target_patient_ids:
            stats["matched_to_checkpoint"] += 1
        
        # Try multiple field names (defensive)
        dfs_status = (
            row.get("DFS_STATUS") or 
            row.get("new_tumor_event_after_initial_treatment") or
            row.get("disease_free_survival_status") or
            row.get("DFS_STATUS_STRING")
        )
        
        dfs_months = (
            parse_months(row.get("DFS_MONTHS")) or
            parse_months(row.get("days_to_new_tumor_event_after_initial_treatment")) or
            parse_months(row.get("disease_free_survival_months"))
        )
        
        vital_status = row.get("OS_STATUS") or row.get("VITAL_STATUS")
        
        # Parse recurrence
        recurrence = parse_dfs_status(dfs_status)
        
        # Build outcome record
        outcome = {
            "patient_id": patient_id,
            "recurrence": recurrence,
            "recurrence_free_survival_months": dfs_months,
            "dfs_status": dfs_status,
            "vital_status": vital_status
        }
        
        # Track stats
        if recurrence is not None:
            stats["with_recurrence_label"] += 1
            if recurrence:
                stats["recurrence_true"] += 1
            else:
                stats["recurrence_false"] += 1
        
        if dfs_months is not None:
            stats["with_dfs_months"] += 1
        
        if vital_status:
            stats["with_vital_status"] += 1
        
        outcomes[patient_id] = outcome
    
    # Step 4: Save results
    print()
    print("💾 Step 4: Saving results...")
    output_data = {
        "study_id": STUDY_ID,
        "extraction_date": datetime.now().isoformat(),
        "n_patients": len(outcomes),
        "stats": stats,
        "outcomes": outcomes,
        "provenance": {
            "source": "cBioPortal",
            "study_id": STUDY_ID,
            "primary_field": "DFS_STATUS",
            "fallback_fields": [
                "new_tumor_event_after_initial_treatment",
                "disease_free_survival_status"
            ],
            "recurrence_definition": "DFS_STATUS == '1:Recurred/Progressed' → recurrence = True",
            "matched_to_checkpoint": target_patient_ids is not None
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"   ✅ Saved to: {OUTPUT_FILE}")
    
    # Step 5: Print summary
    print()
    print("=" * 80)
    print("📊 EXTRACTION SUMMARY")
    print("=" * 80)
    print(f"Total patients in study: {stats['total_patients']}")
    if target_patient_ids:
        print(f"Matched to checkpoint: {stats['matched_to_checkpoint']}")
    print(f"Patients with recurrence labels: {stats['with_recurrence_label']}")
    print(f"  - Recurrence = True: {stats['recurrence_true']}")
    print(f"  - Recurrence = False: {stats['recurrence_false']}")
    print(f"Patients with DFS months: {stats['with_dfs_months']}")
    print(f"Patients with vital status: {stats['with_vital_status']}")
    print()
    
    if stats["with_recurrence_label"] < 100:
        print("⚠️  WARNING: Less than 100 patients with recurrence labels!")
        print("   May need to use alternative field names or data source")
    else:
        print("✅ SUCCESS: Sufficient patients with recurrence labels")
    
    print("=" * 80)
    
    return output_data


if __name__ == "__main__":
    try:
        result = extract_brca_recurrence_labels()
        print("\n✅ EXTRACTION COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ EXTRACTION FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
