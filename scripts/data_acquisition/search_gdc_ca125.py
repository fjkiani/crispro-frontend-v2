#!/usr/bin/env python3
"""
Search GDC (Genomic Data Commons) for CA-125 clinical data in TCGA-OV.
"""

import json
import requests
from datetime import datetime

GDC_API = "https://api.gdc.cancer.gov"

def main():
    print("=" * 80)
    print("🔬 Searching GDC for TCGA-OV CA-125 Clinical Data")
    print("=" * 80)
    
    # 1. Get TCGA-OV clinical data fields
    print("\n1️⃣ Querying TCGA-OV clinical data structure...")
    
    # Get clinical data fields for TCGA-OV
    query = {
        "filters": {
            "op": "=",
            "content": {
                "field": "project.project_id",
                "value": "TCGA-OV"
            }
        },
        "fields": [
            "case_id",
            "submitter_id", 
            "diagnoses.tumor_stage",
            "diagnoses.vital_status",
            "diagnoses.days_to_death",
            "diagnoses.days_to_last_follow_up",
            "diagnoses.treatments.therapeutic_agents",
            "diagnoses.treatments.treatment_type",
            "follow_ups.days_to_follow_up",
            "follow_ups.progression_or_recurrence",
            "follow_ups.tumor_status",
            "exposures.pack_years_smoked"
        ],
        "format": "JSON",
        "size": 50
    }
    
    response = requests.post(
        f"{GDC_API}/cases",
        json=query
    )
    
    if response.ok:
        data = response.json()
        cases = data.get("data", {}).get("hits", [])
        print(f"   ✅ Retrieved {len(cases)} TCGA-OV cases")
        
        # Check for any biomarker/lab data
        print("\n   Sample case structure:")
        if cases:
            print(json.dumps(cases[0], indent=4)[:2000])
    
    # 2. Check for biospecimen/lab results
    print("\n2️⃣ Searching for lab results / biomarker data...")
    
    query2 = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "project.project_id", "value": "TCGA-OV"}},
            ]
        },
        "fields": [
            "case_id",
            "submitter_id"
        ],
        "facets": "files.data_category",
        "format": "JSON",
        "size": 0
    }
    
    response = requests.post(f"{GDC_API}/cases", json=query2)
    if response.ok:
        facets = response.json().get("data", {}).get("aggregations", {})
        print(f"   Data categories in TCGA-OV:")
        if "files.data_category" in facets:
            for bucket in facets["files.data_category"]["buckets"]:
                print(f"      - {bucket['key']}: {bucket['doc_count']} cases")
    
    # 3. Look for clinical supplement files
    print("\n3️⃣ Searching for clinical supplement files...")
    
    file_query = {
        "filters": {
            "op": "and", 
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-OV"}},
                {"op": "=", "content": {"field": "data_category", "value": "Clinical"}}
            ]
        },
        "fields": [
            "file_id",
            "file_name",
            "data_type",
            "data_format"
        ],
        "format": "JSON",
        "size": 100
    }
    
    response = requests.post(f"{GDC_API}/files", json=file_query)
    if response.ok:
        files = response.json().get("data", {}).get("hits", [])
        print(f"   ✅ Found {len(files)} clinical files")
        
        # Look for interesting file types
        for f in files[:20]:
            print(f"      - {f.get('file_name')}: {f.get('data_type')}")
    
    # 4. Check annotations endpoint
    print("\n4️⃣ Checking for annotations/lab values...")
    
    response = requests.get(
        f"{GDC_API}/projects/TCGA-OV",
        params={"fields": "summary.data_categories,summary.experimental_strategies,name,primary_site"}
    )
    if response.ok:
        project = response.json().get("data", {})
        print(f"   Project: {project.get('name')}")
        print(f"   Primary site: {project.get('primary_site')}")
        
        if "summary" in project:
            print("   Data categories:")
            for cat in project.get("summary", {}).get("data_categories", []):
                print(f"      - {cat.get('data_category')}: {cat.get('file_count')} files")
    
    print("\n" + "=" * 80)
    print("📋 CONCLUSION:")
    print("=" * 80)
    print("""
    GDC TCGA-OV contains clinical XML supplements but these typically include:
    - Diagnosis data (stage, grade)
    - Treatment records (chemotherapy regimens)
    - Follow-up data (vital status, survival)
    - Demographics (age, race)
    
    HOWEVER, GDC does NOT typically include:
    - Serial CA-125 laboratory values
    - Longitudinal biomarker measurements
    - Treatment response assessments with lab values
    
    CA-125 serial measurements are usually only available in:
    - Clinical trial databases (NCTN, Project Data Sphere)
    - Institutional datasets (MSK, MSKCC, Dana-Farber)
    - Collaborative study data (ICON7, GOG-218, AOCS)
    
    For KELIM validation, we need COLLABORATOR DATA or synthetic data.
    """)


if __name__ == "__main__":
    main()
