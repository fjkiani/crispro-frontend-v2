#!/usr/bin/env python3
"""
Check all outcome/response fields in TCGA-OV to see which has most coverage.
"""

import httpx
import json
from collections import Counter, defaultdict

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga"

print("Analyzing all outcome/response fields in TCGA-OV...")

try:
    with httpx.Client(timeout=120.0) as client:
        r = client.get(
            f"{CBIO_BASE}/studies/{STUDY_ID}/clinical-data",
            params={
                "clinicalDataType": "PATIENT",
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
                print(f"   ✅ Fetched {len(data)} clinical data rows")
                
                # Group by field
                field_values = defaultdict(list)
                
                for row in data:
                    attr_id = row.get("clinicalAttributeId", "")
                    value = row.get("value", "")
                    
                    if attr_id and value:
                        field_values[attr_id].append(value)
                
                # Find outcome/response-related fields
                outcome_fields = {}
                for field, values in field_values.items():
                    field_upper = field.upper()
                    if any(term in field_upper for term in ["OUTCOME", "RESPONSE", "PLATINUM", "THERAPY", "TREATMENT", "DFS", "OS", "PFS"]):
                        unique_values = set(v for v in values if v and v.strip() and v.upper() not in ["NA", "N/A", "NULL", ""])
                        if unique_values:
                            outcome_fields[field] = {
                                "count": len(unique_values),
                                "coverage": len([v for v in values if v and v.strip() and v.upper() not in ["NA", "N/A", "NULL", ""]]),
                                "total": len(values),
                                "unique_values": list(unique_values)[:10]  # First 10
                            }
                
                print(f"\n   📊 Found {len(outcome_fields)} outcome/response-related fields:")
                print()
                
                # Sort by coverage
                sorted_fields = sorted(outcome_fields.items(), key=lambda x: x[1]["coverage"], reverse=True)
                
                for field, info in sorted_fields[:20]:  # Top 20
                    coverage_pct = (info["coverage"] / info["total"] * 100) if info["total"] > 0 else 0
                    print(f"   {field}:")
                    print(f"      Coverage: {info['coverage']}/{info['total']} ({coverage_pct:.1f}%)")
                    print(f"      Unique values: {info['count']}")
                    print(f"      Sample values: {info['unique_values']}")
                    print()
                
except Exception as e:
    print(f"   ❌ Failed: {e}")
    import traceback
    traceback.print_exc()






