#!/usr/bin/env python3
"""
Attempt to download PTRC-HGSOC data from TCIA using API

TCIA has a REST API that may allow programmatic access.
This script attempts to use the API to download data.
"""

import requests
import json
import sys
from pathlib import Path
from typing import Dict, Optional

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

OUTPUT_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw" / "ptrc_hgsoc"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# TCIA API endpoints
TCIA_API_BASE = "https://services.cancerimagingarchive.net/services/v4"
COLLECTION_NAME = "PTRC-HGSOC"


def get_collection_info(collection_name: str) -> Optional[Dict]:
    """Get collection information from TCIA API."""
    print(f"📥 Fetching collection info for {collection_name}...")
    
    try:
        # Try to get collection details
        url = f"{TCIA_API_BASE}/TCIA/query/getCollectionValues"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            collections = response.json()
            print(f"   ✅ Found {len(collections)} collections")
            
            # Find our collection
            for collection in collections:
                if collection_name.upper() in collection.get("Collection", "").upper():
                    print(f"   ✅ Found collection: {collection.get('Collection')}")
                    return collection
        
        print(f"   ⚠️  Collection not found in API response")
        return None
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return None


def get_patients(collection_name: str) -> Optional[list]:
    """Get patient list for collection."""
    print(f"📥 Fetching patient list for {collection_name}...")
    
    try:
        url = f"{TCIA_API_BASE}/TCIA/query/getPatient"
        params = {"Collection": collection_name}
        response = requests.get(url, params=params, timeout=30)
        
        if response.status_code == 200:
            patients = response.json()
            print(f"   ✅ Found {len(patients)} patients")
            return patients
        else:
            print(f"   ⚠️  API returned {response.status_code}")
            return None
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return None


def get_clinical_data(collection_name: str) -> Optional[list]:
    """Get clinical data for collection."""
    print(f"📥 Fetching clinical data for {collection_name}...")
    
    try:
        # Try different endpoints
        endpoints = [
            f"{TCIA_API_BASE}/TCIA/query/getClinicalData",
            f"{TCIA_API_BASE}/TCIA/query/getCollectionValues",
        ]
        
        for endpoint in endpoints:
            try:
                params = {"Collection": collection_name}
                response = requests.get(endpoint, params=params, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    if data:
                        print(f"   ✅ Found clinical data via {endpoint}")
                        return data
            except:
                continue
        
        print(f"   ⚠️  Clinical data not available via API")
        return None
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return None


def main():
    """Attempt to download PTRC-HGSOC data via TCIA API."""
    print("="*80)
    print("PTRC-HGSOC TCIA API DOWNLOAD ATTEMPT")
    print("="*80)
    print()
    print("⚠️  NOTE: TCIA may require authentication for downloads.")
    print("   This script attempts to use the public API endpoints.")
    print()
    
    # Get collection info
    collection_info = get_collection_info(COLLECTION_NAME)
    
    if collection_info:
        print(f"\n📋 Collection Info:")
        print(json.dumps(collection_info, indent=2))
    
    # Get patients
    patients = get_patients(COLLECTION_NAME)
    
    if patients:
        print(f"\n💾 Saving patient list...")
        patient_file = OUTPUT_DIR / "patients.json"
        with open(patient_file, "w") as f:
            json.dump(patients, f, indent=2)
        print(f"   ✅ Saved: {patient_file}")
        print(f"   ✅ Total patients: {len(patients)}")
    
    # Get clinical data
    clinical_data = get_clinical_data(COLLECTION_NAME)
    
    if clinical_data:
        print(f"\n💾 Saving clinical data...")
        clinical_file = OUTPUT_DIR / "clinical" / "clinical_data.json"
        clinical_file.parent.mkdir(parents=True, exist_ok=True)
        with open(clinical_file, "w") as f:
            json.dump(clinical_data, f, indent=2)
        print(f"   ✅ Saved: {clinical_file}")
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    if patients:
        print(f"✅ Patients: {len(patients)}")
    else:
        print("❌ Patients: Not available via API")
    
    if clinical_data:
        print("✅ Clinical data: Available")
    else:
        print("❌ Clinical data: Not available via API")
        print("\n⚠️  NOTE: Clinical data may require:")
        print("   1. Manual download from TCIA website")
        print("   2. Authentication/login")
        print("   3. Data Use Agreement acceptance")
    
    print("\n📋 Next Steps:")
    print("1. If API doesn't work, download manually from:")
    print("   https://www.cancerimagingarchive.net/collection/ptrc-hgsoc/")
    print("2. Place downloaded files in:")
    print(f"   {OUTPUT_DIR}/")
    print("3. Run extraction script:")
    print("   python3 scripts/data_acquisition/extract_ptrc_hgsoc.py")
    print()


if __name__ == "__main__":
    main()

