#!/usr/bin/env python3
"""
pyBioPortal Treatments Extractor

Uses pyBioPortal library to extract treatment data including platinum response.
"""

import sys
from pathlib import Path
from typing import Dict, List
import pandas as pd

# Add pyBioPortal to path
PYBIOPORTAL_PATH = Path("oncology-coPilot/oncology-backend/tests/pyBioPortal-master")
if PYBIOPORTAL_PATH.exists():
    sys.path.insert(0, str(PYBIOPORTAL_PATH))

def extract_treatments_with_pybioportal(study_id: str = "ov_tcga_pan_can_atlas_2018") -> pd.DataFrame:
    """
    Extract all patient-level treatments using pyBioPortal.
    
    Returns: DataFrame with treatment data including response fields.
    """
    try:
        from pybioportal import treatments as trt
        
        print(f"\n📥 Fetching treatments from pyBioPortal for {study_id}...")
        df = trt.fetch_all_patient_level_treatments(
            study_view_filter={"studyIds": [study_id]}
        )
        
        print(f"   ✅ Fetched {len(df)} treatment records")
        return df
        
    except ImportError:
        print("   ⚠️  pyBioPortal not found - skipping")
        return pd.DataFrame()
    except Exception as e:
        print(f"   ❌ pyBioPortal extraction failed: {e}")
        return pd.DataFrame()

def extract_response_from_treatments(df: pd.DataFrame) -> Dict[str, Dict]:
    """
    Extract platinum response from treatments DataFrame.
    
    Looks for:
    - Response fields in treatment records
    - Platinum drug names (carboplatin, cisplatin, oxaliplatin)
    - Outcome fields
    """
    results = {}
    
    if df.empty:
        return results
    
    print(f"\n📥 Extracting response data from {len(df)} treatment records...")
    
    # Group by patient
    for patient_id, group in df.groupby('patientId'):
        # Look for response fields
        response_value = None
        source_field = None
        
        # Check all columns for response-related data
        for col in df.columns:
            if any(term in col.lower() for term in ["response", "outcome", "status", "result"]):
                values = group[col].dropna().unique()
                if len(values) > 0:
                    # Take first non-null value
                    response_value = str(values[0])
                    source_field = col
                    break
        
        # If no direct response field, check treatment names for platinum
        if not response_value:
            treatment_cols = [c for c in df.columns if "treatment" in c.lower() or "drug" in c.lower()]
            for col in treatment_cols:
                values = group[col].dropna().astype(str).str.lower()
                platinum_drugs = values[values.str.contains("carboplatin|cisplatin|oxaliplatin|platinum", na=False)]
                if len(platinum_drugs) > 0:
                    # Found platinum treatment - check if there's an outcome
                    outcome_cols = [c for c in df.columns if "outcome" in c.lower() or "response" in c.lower()]
                    for outcome_col in outcome_cols:
                        outcomes = group[outcome_col].dropna().unique()
                        if len(outcomes) > 0:
                            response_value = str(outcomes[0])
                            source_field = outcome_col
                            break
        
        if response_value:
            # Import from parent directory
            import sys
            from pathlib import Path as PathLib
            parent_dir = PathLib(__file__).parent.parent.parent
            sys.path.insert(0, str(parent_dir))
            from extract_platinum_response_labels import normalize_response, extract_tcga_patient_id
            normalized = normalize_response(response_value)
            if normalized != "unknown":
                tcga_patient_id = extract_tcga_patient_id(patient_id) or patient_id
                results[tcga_patient_id] = {
                    "platinum_response": normalized,
                    "raw_response_value": response_value,
                    "source_field": source_field or "pybioportal_treatments"
                }
    
    print(f"   ✅ Extracted {len(results)} patients with response labels")
    return results

def extract_platinum_response_from_pybioportal() -> Dict[str, Dict]:
    """
    Main extraction function using pyBioPortal.
    
    Returns: Dict mapping patient_id → {platinum_response, raw_value, source_field}
    """
    print("\n" + "="*80)
    print("PYBIOPORTAL TREATMENTS EXTRACTION")
    print("="*80)
    
    results = {}
    
    # Try both studies
    for study_id in ["ov_tcga_pan_can_atlas_2018", "ov_tcga"]:
        df = extract_treatments_with_pybioportal(study_id)
        if not df.empty:
            study_results = extract_response_from_treatments(df)
            # Merge (don't overwrite existing)
            for patient_id, data in study_results.items():
                if patient_id not in results:
                    results[patient_id] = data
    
    print(f"\n✅ pyBioPortal Extraction: Found {len(results)} patients with response labels")
    return results

if __name__ == "__main__":
    results = extract_platinum_response_from_pybioportal()
    print(f"\n📊 Results: {len(results)} patients")
    for patient_id, data in list(results.items())[:5]:
        print(f"   {patient_id}: {data['platinum_response']} ({data['source_field']})")

