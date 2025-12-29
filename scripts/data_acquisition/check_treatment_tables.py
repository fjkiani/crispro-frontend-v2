#!/usr/bin/env python3
"""
Check Treatment Tables for Platinum Response Data

Many cBioPortal studies have treatment/response data in separate tables.
This script checks if promising studies have treatment endpoints with platinum response.
"""

import requests
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

CBIO_API_BASE = "https://www.cbioportal.org/api"
API_DELAY = 1.0

# Studies to check
STUDIES_TO_CHECK = [
    "msk_spectrum_tme_2022",  # 42 patients, has chemo intent field
    "ovary_cptac_gdc",  # 102 patients
    "lgsoc_mapk_msk_2022",  # 119 patients
    "hgsoc_msk_2021"  # 45 patients
]


def check_treatment_endpoints(study_id: str) -> Dict:
    """Check if study has treatment endpoints."""
    print(f"\n📋 Checking treatment endpoints for {study_id}...")
    
    result = {
        "study_id": study_id,
        "has_treatment_endpoints": False,
        "endpoints": []
    }
    
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/treatment",
            timeout=30
        )
        
        if response.status_code == 200:
            treatments = response.json()
            result["has_treatment_endpoints"] = True
            result["endpoints"] = treatments
            print(f"   ✅ Found treatment endpoints: {len(treatments)}")
            return result
        elif response.status_code == 404:
            print(f"   ❌ No treatment endpoints (404)")
            return result
        else:
            print(f"   ⚠️  API returned {response.status_code}")
            return result
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return result


def check_sample_clinical_data(study_id: str) -> Dict:
    """Check sample-level clinical data (may have treatment response)."""
    print(f"📋 Checking sample-level clinical data for {study_id}...")
    
    result = {
        "study_id": study_id,
        "has_sample_data": False,
        "attributes": []
    }
    
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/clinical-data",
            params={
                "clinicalDataType": "SAMPLE",
                "projection": "DETAILED"
            },
            timeout=60
        )
        
        if response.status_code == 200:
            clinical_data = response.json()
            
            attributes = set()
            for row in clinical_data:
                attr_id = row.get("clinicalAttributeId", "")
                if attr_id:
                    attributes.add(attr_id)
            
            result["has_sample_data"] = True
            result["attributes"] = sorted(list(attributes))
            
            # Check for platinum response fields
            platinum_keywords = ["platinum", "response", "sensitive", "resistant", "refractory", "chemo", "treatment"]
            platinum_fields = [attr for attr in attributes if any(kw in attr.lower() for kw in platinum_keywords)]
            
            result["platinum_fields"] = platinum_fields
            result["has_platinum_response"] = len(platinum_fields) > 0
            
            print(f"   ✅ Found {len(attributes)} sample-level attributes")
            if platinum_fields:
                print(f"   ✅ HAS PLATINUM FIELDS: {', '.join(platinum_fields)}")
            else:
                print(f"   ❌ No platinum response fields")
            
            return result
        else:
            print(f"   ⚠️  API returned {response.status_code}")
            return result
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return result


def main():
    """Check treatment tables for all promising studies."""
    print("="*80)
    print("TREATMENT TABLES CHECK - PLATINUM RESPONSE INVESTIGATION")
    print("="*80)
    
    all_results = []
    
    for study_id in STUDIES_TO_CHECK:
        print(f"\n{'='*80}")
        print(f"Study: {study_id}")
        print('='*80)
        
        # Check treatment endpoints
        treatment_info = check_treatment_endpoints(study_id)
        time.sleep(API_DELAY)
        
        # Check sample-level clinical data
        sample_info = check_sample_clinical_data(study_id)
        time.sleep(API_DELAY)
        
        result = {
            "study_id": study_id,
            "treatment_endpoints": treatment_info,
            "sample_clinical": sample_info
        }
        
        all_results.append(result)
    
    # Save results
    output_file = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw" / "treatment_tables_check.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    studies_with_treatment = [r for r in all_results if r["treatment_endpoints"]["has_treatment_endpoints"]]
    studies_with_sample_platinum = [r for r in all_results if r["sample_clinical"].get("has_platinum_response")]
    
    print(f"\n✅ Studies with treatment endpoints: {len(studies_with_treatment)}")
    for r in studies_with_treatment:
        print(f"   - {r['study_id']}")
    
    print(f"\n✅ Studies with sample-level platinum fields: {len(studies_with_sample_platinum)}")
    for r in studies_with_sample_platinum:
        print(f"   - {r['study_id']}: {', '.join(r['sample_clinical'].get('platinum_fields', []))}")
    
    print(f"\n📁 Results saved to: {output_file}")
    print()


if __name__ == "__main__":
    main()

