#!/usr/bin/env python3
"""
⚔️ HRD ATTRIBUTE INVESTIGATION SCRIPT ⚔️

Mission: Find HRD-related attributes in cBioPortal for TCGA-OV
Study: ov_tcga_pan_can_atlas_2018
"""

import httpx
import json
import sys
from typing import List, Dict, Any

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga_pan_can_atlas_2018"

def _headers() -> Dict[str, str]:
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    import os
    token = os.getenv("CBIO_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers

def get_clinical_attributes(study_id: str) -> List[Dict]:
    """Get all clinical attributes for a study"""
    print(f"🔍 Fetching clinical attributes for {study_id}...")
    
    with httpx.Client(timeout=60.0, headers=_headers()) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/clinical-attributes")
        r.raise_for_status()
        attrs = r.json() or []
    
    print(f"✅ Found {len(attrs)} clinical attributes")
    return attrs

def search_hrd_attributes(attrs: List[Dict]) -> List[Dict]:
    """Search for HRD-related attributes"""
    hrd_keywords = [
        "hrd", "gis", "homologous recombination", "loh", "lst", "tai",
        "genomic instability", "hr deficiency", "hrd score", "gis score"
    ]
    
    hrd_attrs = []
    
    for attr in attrs:
        attr_id = attr.get("clinicalAttributeId", "").lower()
        display_name = attr.get("displayName", "").lower()
        description = attr.get("description", "").lower()
        
        # Check if any keyword matches
        text = f"{attr_id} {display_name} {description}"
        if any(keyword in text for keyword in hrd_keywords):
            hrd_attrs.append(attr)
    
    return hrd_attrs

def get_molecular_profiles(study_id: str) -> List[Dict]:
    """Get all molecular profiles for a study"""
    print(f"🔍 Fetching molecular profiles for {study_id}...")
    
    with httpx.Client(timeout=60.0, headers=_headers()) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/molecular-profiles")
        r.raise_for_status()
        profiles = r.json() or []
    
    print(f"✅ Found {len(profiles)} molecular profiles")
    return profiles

def search_hrd_profiles(profiles: List[Dict]) -> List[Dict]:
    """Search for HRD-related molecular profiles"""
    hrd_keywords = ["hrd", "gis", "copy number", "genomic instability", "loh", "lst", "tai"]
    
    hrd_profiles = []
    
    for profile in profiles:
        name = profile.get("name", "").lower()
        description = profile.get("description", "").lower()
        
        text = f"{name} {description}"
        if any(keyword in text for keyword in hrd_keywords):
            hrd_profiles.append(profile)
    
    return hrd_profiles

def sample_clinical_data(study_id: str, attribute_id: str, n_samples: int = 10) -> List[Dict]:
    """Sample clinical data for a specific attribute"""
    print(f"🔍 Sampling {n_samples} patients with attribute '{attribute_id}'...")
    
    # First get samples
    with httpx.Client(timeout=60.0, headers=_headers()) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/samples")
        r.raise_for_status()
        samples = r.json() or []
    
    if not samples:
        return []
    
    # Get sample IDs
    sample_ids = [s.get("sampleId") for s in samples[:n_samples] if s.get("sampleId")]
    
    if not sample_ids:
        return []
    
    # Fetch clinical data for these samples
    data = {
        "entityIds": sample_ids,
        "entityType": "SAMPLE",
        "projection": "DETAILED",
        "attributeIds": [attribute_id]
    }
    
    with httpx.Client(timeout=120.0, headers=_headers()) as client:
        r = client.post(f"{CBIO_BASE}/clinical-data/fetch", json=data)
        r.raise_for_status()
        result = r.json() or []
    
    return result if isinstance(result, list) else result.get("result", [])

def main():
    print("⚔️ HRD ATTRIBUTE INVESTIGATION ⚔️\n")
    print(f"Study: {STUDY_ID}\n")
    
    # Step 1: Check clinical attributes
    print("=" * 60)
    print("STEP 1: CLINICAL ATTRIBUTES")
    print("=" * 60)
    
    try:
        attrs = get_clinical_attributes(STUDY_ID)
        hrd_attrs = search_hrd_attributes(attrs)
        
        print(f"\n🎯 Found {len(hrd_attrs)} HRD-related clinical attributes:\n")
        
        for attr in hrd_attrs:
            attr_id = attr.get("clinicalAttributeId", "N/A")
            display_name = attr.get("displayName", "N/A")
            datatype = attr.get("datatype", "N/A")
            description = attr.get("description", "N/A")[:100]
            
            print(f"  • {attr_id}")
            print(f"    Display: {display_name}")
            print(f"    Type: {datatype}")
            print(f"    Description: {description}...")
            print()
        
        # Sample data for first HRD attribute
        if hrd_attrs:
            first_attr = hrd_attrs[0]
            attr_id = first_attr.get("clinicalAttributeId")
            print(f"\n📊 Sampling data for '{attr_id}':")
            sample_data = sample_clinical_data(STUDY_ID, attr_id, n_samples=5)
            
            for row in sample_data[:3]:
                entity_id = row.get("entityId", "N/A")
                value = row.get("value", "N/A")
                print(f"  • {entity_id}: {value}")
            
            if len(sample_data) > 3:
                print(f"  ... and {len(sample_data) - 3} more")
        
    except Exception as e:
        print(f"❌ Error fetching clinical attributes: {e}")
        import traceback
        traceback.print_exc()
    
    # Step 2: Check molecular profiles
    print("\n" + "=" * 60)
    print("STEP 2: MOLECULAR PROFILES")
    print("=" * 60)
    
    try:
        profiles = get_molecular_profiles(STUDY_ID)
        hrd_profiles = search_hrd_profiles(profiles)
        
        print(f"\n🎯 Found {len(hrd_profiles)} HRD-related molecular profiles:\n")
        
        for profile in hrd_profiles:
            name = profile.get("name", "N/A")
            profile_id = profile.get("molecularProfileId", "N/A")
            description = profile.get("description", "N/A")[:100]
            
            print(f"  • {profile_id}")
            print(f"    Name: {name}")
            print(f"    Description: {description}...")
            print()
    
    except Exception as e:
        print(f"❌ Error fetching molecular profiles: {e}")
        import traceback
        traceback.print_exc()
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    print(f"\n✅ Investigation complete!")
    print(f"   Clinical attributes checked: {len(attrs) if 'attrs' in locals() else 0}")
    print(f"   HRD-related attributes found: {len(hrd_attrs) if 'hrd_attrs' in locals() else 0}")
    print(f"   Molecular profiles checked: {len(profiles) if 'profiles' in locals() else 0}")
    print(f"   HRD-related profiles found: {len(hrd_profiles) if 'hrd_profiles' in locals() else 0}")
    
    # Save findings
    findings = {
        "study_id": STUDY_ID,
        "clinical_attributes": {
            "total": len(attrs) if 'attrs' in locals() else 0,
            "hrd_related": [
                {
                    "clinicalAttributeId": a.get("clinicalAttributeId"),
                    "displayName": a.get("displayName"),
                    "datatype": a.get("datatype"),
                    "description": a.get("description")
                }
                for a in (hrd_attrs if 'hrd_attrs' in locals() else [])
            ]
        },
        "molecular_profiles": {
            "total": len(profiles) if 'profiles' in locals() else 0,
            "hrd_related": [
                {
                    "molecularProfileId": p.get("molecularProfileId"),
                    "name": p.get("name"),
                    "description": p.get("description")
                }
                for p in (hrd_profiles if 'hrd_profiles' in locals() else [])
            ]
        }
    }
    
    output_file = ".cursor/ayesha/sae_documentation/HRD_ATTRIBUTE_IDS.json"
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(findings, f, indent=2)
    
    print(f"\n💾 Findings saved to: {output_file}")

if __name__ == "__main__":
    main()








