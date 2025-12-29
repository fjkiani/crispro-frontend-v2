#!/usr/bin/env python3
"""
Extract CA-125 data from makeanimpact_ccr_2023 study.
"""

import json
import requests
from datetime import datetime

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "makeanimpact_ccr_2023"

def main():
    print("=" * 80)
    print(f"🔬 Extracting CA-125 data from: {STUDY_ID}")
    print("=" * 80)
    
    # 1. Get study info
    print("\n1️⃣ Getting study info...")
    response = requests.get(f"{BASE_URL}/studies/{STUDY_ID}")
    if response.ok:
        study = response.json()
        print(f"   Study: {study.get('name')}")
        print(f"   Description: {study.get('description', 'N/A')[:200]}")
        print(f"   Patients: {study.get('allPatientCount')}")
        print(f"   Samples: {study.get('allSampleCount')}")
    
    # 2. Get all clinical attributes
    print("\n2️⃣ Getting all clinical attributes...")
    response = requests.get(f"{BASE_URL}/studies/{STUDY_ID}/clinical-attributes")
    if response.ok:
        attrs = response.json()
        print(f"   Total attributes: {len(attrs)}")
        
        # Find CA-125 and related
        relevant = []
        for attr in attrs:
            attr_id = attr.get("clinicalAttributeId", "").lower()
            if any(kw in attr_id for kw in ["ca_125", "ca125", "date", "time", "cycle", 
                                           "platinum", "cancer", "stage", "chemo"]):
                relevant.append(attr)
                print(f"   - {attr.get('clinicalAttributeId')}: {attr.get('displayName')}")
    
    # 3. Fetch CA-125 data
    print("\n3️⃣ Fetching CA-125 values...")
    for data_type in ["PATIENT", "SAMPLE"]:
        response = requests.get(
            f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
            params={
                "clinicalDataType": data_type,
                "attributeId": "CA_125",
                "pageSize": 10000
            }
        )
        if response.ok and response.json():
            data = response.json()
            print(f"   ✅ Found {len(data)} CA-125 records ({data_type} level)")
            
            # Show sample values
            print("\n   Sample CA-125 values:")
            values = []
            for record in data[:20]:
                patient_id = record.get("patientId")
                value = record.get("value")
                print(f"      {patient_id}: {value}")
                try:
                    values.append(float(value))
                except:
                    pass
            
            if values:
                print(f"\n   Stats: min={min(values):.1f}, max={max(values):.1f}, mean={sum(values)/len(values):.1f}")
    
    # 4. Check cancer types in study
    print("\n4️⃣ Checking cancer types...")
    response = requests.get(
        f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
        params={
            "clinicalDataType": "PATIENT",
            "attributeId": "CANCER_TYPE",
            "pageSize": 10000
        }
    )
    if response.ok:
        data = response.json()
        cancer_types = {}
        for record in data:
            ctype = record.get("value", "Unknown")
            cancer_types[ctype] = cancer_types.get(ctype, 0) + 1
        
        print("   Cancer type distribution:")
        for ctype, count in sorted(cancer_types.items(), key=lambda x: -x[1])[:15]:
            print(f"      {ctype}: {count}")
    
    # 5. Get mutation data for patients with CA-125
    print("\n5️⃣ Checking for ovarian cancer patients with CA-125...")
    response = requests.get(
        f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
        params={
            "clinicalDataType": "PATIENT", 
            "pageSize": 100000
        }
    )
    
    if response.ok:
        all_clinical = response.json()
        
        # Build patient profile
        patients = {}
        for record in all_clinical:
            pid = record.get("patientId")
            if pid not in patients:
                patients[pid] = {}
            patients[pid][record.get("clinicalAttributeId")] = record.get("value")
        
        # Filter for ovarian with CA-125
        ovarian_with_ca125 = []
        for pid, data in patients.items():
            cancer = data.get("CANCER_TYPE", "").lower()
            ca125 = data.get("CA_125")
            
            if "ovarian" in cancer and ca125:
                ovarian_with_ca125.append({
                    "patient_id": pid,
                    "cancer_type": data.get("CANCER_TYPE"),
                    "ca125": ca125,
                    "stage": data.get("STAGE"),
                    "os_status": data.get("OS_STATUS"),
                    "os_months": data.get("OS_MONTHS")
                })
        
        print(f"\n   🎯 Found {len(ovarian_with_ca125)} ovarian patients with CA-125!")
        if ovarian_with_ca125:
            print("\n   Sample ovarian patients with CA-125:")
            for p in ovarian_with_ca125[:10]:
                print(f"      {p['patient_id']}: CA-125={p['ca125']}, Stage={p.get('stage')}")
    
    # Save all data
    output = {
        "study_id": STUDY_ID,
        "scan_timestamp": datetime.now().isoformat(),
        "ovarian_patients_with_ca125": ovarian_with_ca125 if 'ovarian_with_ca125' in dir() else [],
        "ca125_patient_count": len(ovarian_with_ca125) if 'ovarian_with_ca125' in dir() else 0
    }
    
    output_path = "/Users/fahadkiani/Desktop/development/crispr-assistant-main/data/validation/kelim_resurrection/impact_ca125_extract.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"\n💾 Saved to: {output_path}")


if __name__ == "__main__":
    main()
