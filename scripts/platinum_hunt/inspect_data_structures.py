#!/usr/bin/env python3
"""
Inspect actual data structures to find response fields.
"""

import sys
from pathlib import Path
import pandas as pd
import httpx
import json

# Add pyBioPortal
PYBIOPORTAL_PATH = Path("oncology-coPilot/oncology-backend/tests/pyBioPortal-master")
if PYBIOPORTAL_PATH.exists():
    sys.path.insert(0, str(PYBIOPORTAL_PATH))

print("="*80)
print("INSPECTING DATA STRUCTURES FOR RESPONSE FIELDS")
print("="*80)

# 1. Inspect pyBioPortal treatments
print("\n1. PYBIOPORTAL TREATMENTS STRUCTURE")
print("-"*80)
try:
    from pybioportal import treatments as trt
    df = trt.fetch_all_patient_level_treatments(
        study_view_filter={"studyIds": ["ov_tcga_pan_can_atlas_2018"]}
    )
    print(f"   Columns: {list(df.columns)[:20]}...")
    print(f"   Total columns: {len(df.columns)}")
    print(f"\n   Sample row:")
    if len(df) > 0:
        print(df.iloc[0].to_dict())
    
    # Check for response-related columns
    response_cols = [c for c in df.columns if any(term in c.lower() for term in 
                  ["response", "outcome", "status", "result", "effect", "efficacy"])]
    print(f"\n   Response-related columns: {response_cols}")
    
except Exception as e:
    print(f"   ❌ Failed: {e}")

# 2. Inspect Broad Firehose
print("\n2. BROAD FIREHOSE STRUCTURE")
print("-"*80)
try:
    import tarfile
    tar_path = Path("data/validation/broad_firehose_cache/OV.Clinical_Pick_Tier1.tar.gz")
    if tar_path.exists():
        with tarfile.open(tar_path, 'r:gz') as tar:
            members = tar.getmembers()
            print(f"   Files in archive: {[m.name for m in members[:10]]}")
            
            # Extract first clinical file
            clinical_files = [m for m in members if 'clinical' in m.name.lower() and m.name.endswith('.txt')]
            if clinical_files:
                extracted = tar.extractfile(clinical_files[0])
                df = pd.read_csv(extracted, sep='\t', nrows=5, low_memory=False)
                print(f"\n   Columns: {list(df.columns)[:30]}...")
                print(f"   Total columns: {len(df.columns)}")
                
                # Check for response columns
                response_cols = [c for c in df.columns if any(term in c.lower() for term in 
                          ["response", "outcome", "platinum", "therapy", "treatment", "status"])]
                print(f"\n   Response-related columns: {response_cols}")
                
                if response_cols:
                    print(f"\n   Sample values from {response_cols[0]}:")
                    extracted2 = tar.extractfile(clinical_files[0])
                    df2 = pd.read_csv(extracted2, sep='\t', low_memory=False)
                    print(df2[response_cols[0]].value_counts().head(10))
    else:
        print("   ⚠️  Broad Firehose file not downloaded yet")
except Exception as e:
    print(f"   ❌ Failed: {e}")

# 3. Inspect GDC XML structure
print("\n3. GDC XML STRUCTURE")
print("-"*80)
try:
    import xml.etree.ElementTree as ET
    cache_dir = Path("data/validation/gdc_cache")
    xml_files = list(cache_dir.glob("*.xml"))[:3]
    
    if xml_files:
        print(f"   Found {len(xml_files)} cached XML files")
        for xml_file in xml_files[:1]:
            print(f"\n   Parsing {xml_file.name}...")
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            # Print all unique tags
            tags = set()
            for elem in root.iter():
                tags.add(elem.tag)
            
            print(f"   Unique tags: {sorted(list(tags))[:30]}...")
            
            # Look for response-related tags
            response_tags = [t for t in tags if any(term in t.lower() for term in 
                          ["response", "outcome", "platinum", "therapy", "treatment"])]
            print(f"\n   Response-related tags: {response_tags}")
            
            if response_tags:
                print(f"\n   Sample values from {response_tags[0]}:")
                for elem in root.iter(response_tags[0]):
                    if elem.text:
                        print(f"      {elem.text[:100]}")
    else:
        print("   ⚠️  No GDC XML files cached yet")
except Exception as e:
    print(f"   ❌ Failed: {e}")

# 4. Check cBioPortal clinical attributes
print("\n4. CBIOPORTAL CLINICAL ATTRIBUTES")
print("-"*80)
try:
    CBIO_BASE = "https://www.cbioportal.org/api"
    with httpx.Client(timeout=60.0) as client:
        r = client.get(
            f"{CBIO_BASE}/studies/ov_tcga/clinical-attributes",
            params={"projection": "DETAILED"}
        )
        if r.status_code == 200:
            attrs = r.json()
            response_attrs = [a for a in attrs if any(term in a.get("displayName", "").lower() for term in 
                          ["response", "outcome", "platinum", "therapy", "treatment"])]
            print(f"   Response-related attributes: {[a.get('clinicalAttributeId') for a in response_attrs[:10]]}")
            print(f"\n   All attributes with 'outcome' or 'response':")
            for attr in response_attrs[:10]:
                print(f"      {attr.get('clinicalAttributeId')}: {attr.get('displayName')}")
except Exception as e:
    print(f"   ❌ Failed: {e}")

print("\n" + "="*80)
print("INSPECTION COMPLETE")
print("="*80)






