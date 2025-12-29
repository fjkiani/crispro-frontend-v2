#!/usr/bin/env python3
"""
KELIM DATA HUNT V2: Search ALL cBioPortal studies for ANY CA-125 data.
"""

import sys
import os
import json
import requests
from datetime import datetime
from typing import List, Dict, Any

BASE_URL = "https://www.cbioportal.org/api"

def main():
    print("=" * 80)
    print("🔬 KELIM DATA HUNT V2: Searching ALL cBioPortal for CA-125")
    print("=" * 80)
    
    # Get ALL clinical attributes across all studies
    print("\n🔍 Fetching ALL clinical attributes from cBioPortal...")
    
    try:
        response = requests.get(f"{BASE_URL}/clinical-attributes", params={"pageSize": 100000})
        response.raise_for_status()
        all_attrs = response.json()
        print(f"✅ Retrieved {len(all_attrs)} clinical attributes")
    except Exception as e:
        print(f"❌ Failed: {e}")
        return
    
    # Search for CA-125 related attributes
    ca125_keywords = ["ca125", "ca-125", "ca_125", "ca 125", "serum_marker", "tumor_marker", 
                      "mucin", "muc16", "serum marker"]
    
    ca125_attrs = []
    for attr in all_attrs:
        attr_id = attr.get("clinicalAttributeId", "").lower()
        display_name = attr.get("displayName", "").lower()
        description = attr.get("description", "").lower() if attr.get("description") else ""
        
        if any(kw in attr_id or kw in display_name or kw in description for kw in ca125_keywords):
            ca125_attrs.append(attr)
    
    print(f"\n🎯 Found {len(ca125_attrs)} CA-125 related attributes:")
    for attr in ca125_attrs:
        print(f"  - Study: {attr.get('studyId')}")
        print(f"    ID: {attr.get('clinicalAttributeId')}")
        print(f"    Name: {attr.get('displayName')}")
        print(f"    Desc: {attr.get('description', 'N/A')[:80]}")
        print()
    
    # Also search for time-series related attributes
    print("\n" + "=" * 80)
    print("🕒 Searching for time-series/longitudinal attributes...")
    
    time_keywords = ["baseline", "cycle", "day_", "week_", "month_", "pre_treatment", 
                    "post_treatment", "followup", "follow_up", "visit", "timepoint",
                    "serial", "longitudinal", "c1d1", "c2d1", "c3d1"]
    
    time_attrs = []
    for attr in all_attrs:
        attr_id = attr.get("clinicalAttributeId", "").lower()
        if any(kw in attr_id for kw in time_keywords):
            time_attrs.append(attr)
    
    # Group by study
    time_by_study = {}
    for attr in time_attrs:
        study = attr.get("studyId")
        if study not in time_by_study:
            time_by_study[study] = []
        time_by_study[study].append(attr.get("clinicalAttributeId"))
    
    print(f"Found {len(time_by_study)} studies with time-series attributes:")
    for study, attrs in sorted(time_by_study.items(), key=lambda x: len(x[1]), reverse=True)[:20]:
        print(f"  - {study}: {len(attrs)} time attrs → {attrs[:5]}")
    
    # Save results
    results = {
        "scan_timestamp": datetime.now().isoformat(),
        "total_clinical_attributes": len(all_attrs),
        "ca125_related_attributes": ca125_attrs,
        "ca125_count": len(ca125_attrs),
        "time_series_studies": time_by_study
    }
    
    output_dir = "/Users/fahadkiani/Desktop/development/crispr-assistant-main/data/validation/kelim_resurrection"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "cbioportal_ca125_global_scan.json")
    
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Saved to: {output_path}")


if __name__ == "__main__":
    main()
