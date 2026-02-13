#!/usr/bin/env python3
"""
Check sample-level clinical data for response fields.
"""

import httpx
import json
from collections import Counter

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga"

print("Checking SAMPLE-level clinical data...")

try:
    with httpx.Client(timeout=120.0) as client:
        # Fetch SAMPLE-level clinical
        r = client.get(
            f"{CBIO_BASE}/studies/{STUDY_ID}/clinical-data",
            params={
                "clinicalDataType": "SAMPLE",
                "projection": "DETAILED"
            }
        )
        
        if r.status_code != 200:
            print(f"   ❌ API returned {r.status_code}")
        else:
            data = r.json()
            
            if not isinstance(data, list):
                print(f"   ⚠️  Unexpected format: {type(data)}")
            else:
                print(f"   ✅ Fetched {len(data)} SAMPLE-level clinical data rows")
                
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
                                "sample": row.get("sampleId") or row.get("entityId")
                            })
                
                # Show unique fields
                field_counts = Counter(all_fields)
                print(f"\n   📊 Total unique fields: {len(field_counts)}")
                print(f"\n   🔍 Response-related fields found: {len(response_candidates)}")
                
                if response_candidates:
                    print("\n   Response Candidate Fields (first 20):")
                    seen_fields = set()
                    for item in response_candidates[:20]:
                        if item["field"] not in seen_fields:
                            print(f"      - {item['field']}: '{item['value']}' (sample: {item['sample']})")
                            seen_fields.add(item["field"])
                
                # Check TREATMENT_OUTCOME_FIRST_COURSE specifically
                treatment_outcome = [row for row in data if row.get("clinicalAttributeId") == "TREATMENT_OUTCOME_FIRST_COURSE"]
                print(f"\n   TREATMENT_OUTCOME_FIRST_COURSE: {len(treatment_outcome)} rows")
                if treatment_outcome:
                    values = Counter(row.get("value") for row in treatment_outcome if row.get("value"))
                    print(f"      Values: {dict(values)}")
                
except Exception as e:
    print(f"   ❌ Failed: {e}")
    import traceback
    traceback.print_exc()






