#!/usr/bin/env python3
"""
KELIM DATA HUNT: Scan ALL cBioPortal ovarian cancer studies for CA-125 data.
"""

import sys
import os
import json
import requests
from datetime import datetime
from typing import List, Dict, Any

BASE_URL = "https://www.cbioportal.org/api"

def get_all_ovarian_studies() -> List[Dict]:
    """Get all ovarian cancer studies from cBioPortal."""
    print("\n🔍 Fetching all studies from cBioPortal...")
    
    try:
        response = requests.get(f"{BASE_URL}/studies", params={"pageSize": 10000})
        response.raise_for_status()
        all_studies = response.json()
    except Exception as e:
        print(f"❌ Failed to fetch studies: {e}")
        return []
    
    # Filter for ovarian cancer studies
    ovarian_keywords = ["ovarian", "ov_", "_ov", "hgsc", "ovary", "gynecologic"]
    ovarian_studies = []
    
    for study in all_studies:
        study_id = study.get("studyId", "").lower()
        name = study.get("name", "").lower()
        description = study.get("description", "").lower()
        cancer_type = study.get("cancerTypeId", "").lower()
        
        if any(kw in study_id or kw in name or kw in description or kw in cancer_type 
               for kw in ovarian_keywords):
            ovarian_studies.append(study)
    
    print(f"✅ Found {len(ovarian_studies)} ovarian-related studies out of {len(all_studies)} total")
    return ovarian_studies


def scan_study_for_ca125(study_id: str) -> Dict[str, Any]:
    """Scan a single study for CA-125 related clinical attributes."""
    result = {
        "study_id": study_id,
        "ca125_attributes": [],
        "longitudinal_potential": False,
        "sample_count": 0,
        "patient_count": 0,
        "error": None
    }
    
    try:
        response = requests.get(
            f"{BASE_URL}/studies/{study_id}/clinical-attributes",
            params={"pageSize": 10000}
        )
        response.raise_for_status()
        attrs = response.json()
        
        ca125_keywords = ["ca125", "ca-125", "ca_125", "ca 125", "serum_marker", "tumor_marker"]
        
        for attr in attrs:
            attr_id = attr.get("clinicalAttributeId", "").lower()
            display_name = attr.get("displayName", "").lower()
            description = attr.get("description", "").lower()
            
            if any(kw in attr_id or kw in display_name or kw in description for kw in ca125_keywords):
                result["ca125_attributes"].append({
                    "attribute_id": attr.get("clinicalAttributeId"),
                    "display_name": attr.get("displayName"),
                    "description": attr.get("description"),
                    "datatype": attr.get("datatype"),
                    "patient_attribute": attr.get("patientAttribute", False)
                })
        
        time_keywords = ["date", "time", "cycle", "visit", "baseline", "followup", "follow_up", 
                        "day_", "week_", "month_", "pre_", "post_", "c1", "c2", "c3"]
        
        longitudinal_attrs = []
        for attr in attrs:
            attr_id = attr.get("clinicalAttributeId", "").lower()
            if any(kw in attr_id for kw in time_keywords):
                longitudinal_attrs.append(attr.get("clinicalAttributeId"))
        
        if longitudinal_attrs and result["ca125_attributes"]:
            result["longitudinal_potential"] = True
            result["time_related_attributes"] = longitudinal_attrs[:10]
        
        try:
            sample_response = requests.get(f"{BASE_URL}/studies/{study_id}")
            if sample_response.ok:
                study_info = sample_response.json()
                result["sample_count"] = study_info.get("allSampleCount", 0)
                result["patient_count"] = study_info.get("allPatientCount", 0) 
        except:
            pass
            
    except Exception as e:
        result["error"] = str(e)
    
    return result


def fetch_ca125_sample_data(study_id: str, attribute_id: str, limit: int = 50) -> List[Dict]:
    """Fetch actual CA-125 values from a study."""
    for data_type in ["PATIENT", "SAMPLE"]:
        try:
            response = requests.get(
                f"{BASE_URL}/studies/{study_id}/clinical-data",
                params={
                    "clinicalDataType": data_type,
                    "attributeId": attribute_id,
                    "pageSize": limit
                }
            )
            if response.ok and response.json():
                return response.json()
        except:
            pass
    return []


def main():
    print("=" * 80)
    print("🔬 KELIM DATA HUNT: Scanning cBioPortal for CA-125 Longitudinal Data")
    print("=" * 80)
    print(f"Started: {datetime.now().isoformat()}")
    
    ovarian_studies = get_all_ovarian_studies()
    
    if not ovarian_studies:
        print("❌ No ovarian studies found!")
        return
    
    print(f"\n📊 Studies to scan: {len(ovarian_studies)}")
    for s in ovarian_studies:
        print(f"  - {s.get('studyId')}: {s.get('name', 'N/A')[:60]}")
    
    results = {
        "scan_timestamp": datetime.now().isoformat(),
        "total_ovarian_studies": len(ovarian_studies),
        "studies_with_ca125": [],
        "studies_with_longitudinal_potential": [],
        "all_results": []
    }
    
    print("\n" + "=" * 80)
    print("🔍 Scanning for CA-125 attributes...")
    print("=" * 80)
    
    for i, study in enumerate(ovarian_studies):
        study_id = study.get("studyId")
        print(f"\n[{i+1}/{len(ovarian_studies)}] Scanning: {study_id}")
        
        scan_result = scan_study_for_ca125(study_id)
        scan_result["study_name"] = study.get("name")
        scan_result["cancer_type"] = study.get("cancerTypeId")
        
        if scan_result["ca125_attributes"]:
            print(f"  ✅ FOUND CA-125 attributes: {len(scan_result['ca125_attributes'])}")
            for attr in scan_result["ca125_attributes"]:
                print(f"     - {attr['attribute_id']}: {attr['display_name']}")
            results["studies_with_ca125"].append(scan_result)
            
            for attr in scan_result["ca125_attributes"][:1]:
                sample_data = fetch_ca125_sample_data(study_id, attr["attribute_id"])
                if sample_data:
                    print(f"     📈 Sample values: {len(sample_data)} records")
                    scan_result["sample_values"] = sample_data[:5]
        else:
            print(f"  ❌ No CA-125 attributes found")
        
        if scan_result["longitudinal_potential"]:
            print(f"  ⏱️  Longitudinal potential detected!")
            results["studies_with_longitudinal_potential"].append(study_id)
        
        if scan_result["error"]:
            print(f"  ⚠️  Error: {scan_result['error']}")
        
        results["all_results"].append(scan_result)
    
    print("\n" + "=" * 80)
    print("📋 SCAN SUMMARY")
    print("=" * 80)
    print(f"Total ovarian studies scanned: {len(ovarian_studies)}")
    print(f"Studies with CA-125 attributes: {len(results['studies_with_ca125'])}")
    print(f"Studies with longitudinal potential: {len(results['studies_with_longitudinal_potential'])}")
    
    if results["studies_with_ca125"]:
        print("\n🎯 STUDIES WITH CA-125 DATA:")
        for study in results["studies_with_ca125"]:
            print(f"\n  📌 {study['study_id']}")
            print(f"     Name: {study.get('study_name', 'N/A')[:60]}")
            print(f"     Patients: {study.get('patient_count', 'N/A')}")
            print(f"     Samples: {study.get('sample_count', 'N/A')}")
            print(f"     CA-125 attributes: {[a['attribute_id'] for a in study['ca125_attributes']]}")
            print(f"     Longitudinal: {'YES ⏱️' if study['longitudinal_potential'] else 'NO'}")
            if study.get("sample_values"):
                print(f"     Sample values: {[v.get('value') for v in study['sample_values'][:3]]}")
    
    output_dir = "/Users/fahadkiani/Desktop/development/crispr-assistant-main/data/validation/kelim_resurrection"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "cbioportal_ca125_scan_report.json")
    
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Full report saved to: {output_path}")
    print(f"Completed: {datetime.now().isoformat()}")
    
    return results


if __name__ == "__main__":
    main()
