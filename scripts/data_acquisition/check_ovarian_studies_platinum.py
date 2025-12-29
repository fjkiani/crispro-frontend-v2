#!/usr/bin/env python3
"""
Check Ovarian Studies for Platinum Response Data

Investigates promising ovarian cancer studies from cBioPortal to find:
- Platinum response labels in clinical data
- Patient counts (need ≥50)
- Mutation data availability
"""

import json
import requests
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

CBIO_API_BASE = "https://www.cbioportal.org/api"
API_DELAY = 1.0

# Promising ovarian studies to investigate
OVARIAN_STUDIES = [
    {
        "study_id": "hgsoc_msk_2021",
        "name": "High-Grade Serous Ovarian Cancer (MSK, NPJ Genome Med 2021)",
        "expected_patients": 45
    },
    {
        "study_id": "msk_spectrum_tme_2022",
        "name": "Ovarian Cancer - MSK SPECTRUM (MSK, Nature 2022)",
        "expected_patients": 42
    },
    {
        "study_id": "lgsoc_mapk_msk_2022",
        "name": "Low-Grade Serous Ovarian Cancer (MSK, Clin Cancer Res 2022)",
        "expected_patients": 119  # Meets minimum!
    },
    {
        "study_id": "ovary_cptac_gdc",
        "name": "Ovarian Cancer (CPTAC GDC, 2025)",
        "expected_patients": "unknown"
    },
    {
        "study_id": "ucec_msk_2024",
        "name": "Endometrial and Ovarian Cancer (MSK, Nature Medicine 2024)",
        "expected_patients": 33
    }
]


def check_study_clinical_attributes(study_id: str) -> Dict:
    """Check what clinical attributes are available for a study."""
    print(f"\n📋 Checking clinical attributes for {study_id}...")
    
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/clinical-data",
            params={
                "clinicalDataType": "PATIENT",
                "projection": "DETAILED"
            },
            timeout=60
        )
        
        if response.status_code == 200:
            clinical_data = response.json()
            
            # Extract unique attribute IDs
            attributes = set()
            for row in clinical_data:
                attr_id = row.get("clinicalAttributeId", "")
                if attr_id:
                    attributes.add(attr_id)
            
            # Check for platinum response related fields
            platinum_keywords = ["platinum", "response", "sensitive", "resistant", "refractory", "chemo", "treatment"]
            platinum_fields = [attr for attr in attributes if any(kw in attr.lower() for kw in platinum_keywords)]
            
            return {
                "study_id": study_id,
                "total_attributes": len(attributes),
                "all_attributes": sorted(list(attributes)),
                "platinum_related_fields": platinum_fields,
                "has_platinum_response": len(platinum_fields) > 0,
                "status": "success"
            }
        else:
            return {
                "study_id": study_id,
                "status": f"error_{response.status_code}",
                "error": response.text[:200]
            }
    
    except Exception as e:
        return {
            "study_id": study_id,
            "status": "error",
            "error": str(e)
        }


def check_study_patient_count(study_id: str) -> Optional[int]:
    """Get patient count for a study."""
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/patients",
            params={"pageSize": 10000},
            timeout=30
        )
        
        if response.status_code == 200:
            patients = response.json()
            return len(patients)
        return None
    except:
        return None


def check_study_has_mutations(study_id: str) -> bool:
    """Check if study has mutation data."""
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/molecular-profiles",
            timeout=30
        )
        
        if response.status_code == 200:
            profiles = response.json()
            for profile in profiles:
                if "mutation" in profile.get("molecularAlterationType", "").lower():
                    return True
        return False
    except:
        return False


def main():
    """Check all promising ovarian studies."""
    print("="*80)
    print("OVARIAN STUDIES - PLATINUM RESPONSE INVESTIGATION")
    print("="*80)
    
    results = []
    
    for study in OVARIAN_STUDIES:
        study_id = study["study_id"]
        print(f"\n{'='*80}")
        print(f"Study: {study['name']}")
        print(f"Study ID: {study_id}")
        print('='*80)
        
        # Check patient count
        print("📊 Checking patient count...")
        patient_count = check_study_patient_count(study_id)
        time.sleep(API_DELAY)
        
        if patient_count:
            print(f"   ✅ Patient count: {patient_count}")
            meets_minimum = patient_count >= 50
            print(f"   {'✅ Meets minimum (≥50)' if meets_minimum else '❌ Below minimum (<50)'}")
        else:
            print(f"   ⚠️  Could not determine patient count")
            meets_minimum = False
        
        # Check for mutations
        print("🧬 Checking for mutation data...")
        has_mutations = check_study_has_mutations(study_id)
        time.sleep(API_DELAY)
        print(f"   {'✅ Has mutation data' if has_mutations else '❌ No mutation data'}")
        
        # Check clinical attributes
        clinical_info = check_study_clinical_attributes(study_id)
        time.sleep(API_DELAY)
        
        if clinical_info.get("status") == "success":
            print(f"   ✅ Found {clinical_info['total_attributes']} clinical attributes")
            
            if clinical_info["has_platinum_response"]:
                print(f"   ✅ HAS PLATINUM-RELATED FIELDS:")
                for field in clinical_info["platinum_related_fields"]:
                    print(f"      - {field}")
            else:
                print(f"   ❌ No platinum response fields found")
                print(f"   Sample attributes: {', '.join(clinical_info['all_attributes'][:10])}")
        
        # Compile result
        result = {
            "study_id": study_id,
            "name": study["name"],
            "patient_count": patient_count,
            "meets_minimum": meets_minimum,
            "has_mutations": has_mutations,
            "clinical_attributes": clinical_info
        }
        
        results.append(result)
    
    # Save results
    output_file = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw" / "ovarian_studies_platinum_check.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    promising = [r for r in results if r["meets_minimum"] and r["has_mutations"] and r["clinical_attributes"].get("has_platinum_response")]
    
    print(f"\n✅ Promising studies (meets minimum + mutations + platinum response): {len(promising)}")
    for r in promising:
        print(f"   - {r['study_id']}: {r['patient_count']} patients")
        print(f"     Platinum fields: {', '.join(r['clinical_attributes'].get('platinum_related_fields', []))}")
    
    print(f"\n⚠️  Studies meeting minimum but no platinum response: {len([r for r in results if r['meets_minimum'] and r['has_mutations'] and not r['clinical_attributes'].get('has_platinum_response')])}")
    
    print(f"\n📁 Results saved to: {output_file}")
    print()


if __name__ == "__main__":
    main()

