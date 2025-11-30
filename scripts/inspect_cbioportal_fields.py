#!/usr/bin/env python3
"""
Quick script to inspect what response fields are actually available in cBioPortal.
"""

import httpx
import json
from collections import Counter

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga"
PAN_CANCER_STUDY_ID = "ov_tcga_pan_can_atlas_2018"

def inspect_study(study_id: str, study_name: str):
    print(f"\n{'='*80}")
    print(f"INSPECTING: {study_name} ({study_id})")
    print('='*80)
    
    try:
        with httpx.Client(timeout=120.0) as client:
            # Fetch patient-level clinical
            r = client.get(
                f"{CBIO_BASE}/studies/{study_id}/clinical-data",
                params={
                    "clinicalDataType": "PATIENT",
                    "projection": "DETAILED"
                }
            )
            
            if r.status_code != 200:
                print(f"   ❌ API returned {r.status_code}")
                return
            
            data = r.json()
            
            if not isinstance(data, list):
                print(f"   ⚠️  Unexpected format: {type(data)}")
                return
            
            print(f"   ✅ Fetched {len(data)} clinical data rows")
            
            # Extract all field names
            all_fields = []
            response_candidates = []
            
            for row in data:
                attr_id = row.get("clinicalAttributeId", "")
                value = row.get("value", "")
                
                if attr_id:
                    all_fields.append(attr_id)
                    
                    # Check if it might be response-related
                    attr_upper = attr_id.upper()
                    if any(term in attr_upper for term in ["RESPONSE", "PLATINUM", "OUTCOME", "THERAPY", "TREATMENT"]):
                        response_candidates.append({
                            "field": attr_id,
                            "value": value,
                            "patient": row.get("patientId") or row.get("entityId")
                        })
            
            # Show unique fields
            field_counts = Counter(all_fields)
            print(f"\n   📊 Total unique fields: {len(field_counts)}")
            print(f"\n   🔍 Response-related fields found: {len(response_candidates)}")
            
            if response_candidates:
                print("\n   Response Candidate Fields:")
                seen_fields = set()
                for item in response_candidates[:20]:  # Show first 20
                    if item["field"] not in seen_fields:
                        print(f"      - {item['field']}: '{item['value']}' (patient: {item['patient']})")
                        seen_fields.add(item["field"])
            
            # Show all fields (top 30)
            print(f"\n   All Clinical Fields (top 30 by frequency):")
            for field, count in field_counts.most_common(30):
                print(f"      - {field}: {count} occurrences")
            
    except Exception as e:
        print(f"   ❌ Failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    inspect_study(STUDY_ID, "TCGA-OV Original")
    inspect_study(PAN_CANCER_STUDY_ID, "TCGA PanCancer Atlas")






