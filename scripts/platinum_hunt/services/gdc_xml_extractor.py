#!/usr/bin/env python3
"""
GDC XML Clinical File Extractor

Downloads and parses GDC XML clinical supplement files to extract platinum response data.
"""

import httpx
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Optional
import tempfile
import zipfile

GDC_BASE = "https://api.gdc.cancer.gov"
GDC_DATA_BASE = "https://api.gdc.cancer.gov/data"

def download_gdc_file(file_id: str, output_path: Path) -> bool:
    """Download a file from GDC by file ID."""
    try:
        with httpx.Client(timeout=300.0, follow_redirects=True) as client:
            r = client.get(
                f"{GDC_DATA_BASE}/{file_id}",
                headers={"Content-Type": "application/json"}
            )
            r.raise_for_status()
            
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(r.content)
            return True
    except Exception as e:
        print(f"   ⚠️  Failed to download {file_id}: {e}")
        return False

def parse_gdc_xml(xml_path: Path) -> Dict[str, str]:
    """Parse GDC XML clinical file and extract response fields."""
    results = {}
    
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        # GDC XML structure: <clinical_study> or <patient> or similar
        # Look for response-related fields
        for elem in root.iter():
            tag_lower = elem.tag.lower()
            text = elem.text.strip() if elem.text else ""
            
            # Check tag names
            if any(term in tag_lower for term in ["response", "outcome", "platinum", "therapy", "treatment"]):
                if text:
                    results[elem.tag] = text
            
            # Check attributes
            for attr, value in elem.attrib.items():
                attr_lower = attr.lower()
                if any(term in attr_lower for term in ["response", "outcome", "platinum"]):
                    if value:
                        results[f"{elem.tag}@{attr}"] = value
        
        # Also check for common GDC fields
        # primary_therapy_outcome_success, response_to_treatment, etc.
        for field_name in ["primary_therapy_outcome_success", "response_to_treatment", 
                          "platinum_status", "treatment_outcome", "therapy_response"]:
            for elem in root.iter(field_name):
                if elem.text:
                    results[field_name] = elem.text.strip()
        
    except Exception as e:
        print(f"   ⚠️  Failed to parse XML {xml_path}: {e}")
    
    return results

def extract_platinum_response_from_gdc_xml() -> Dict[str, Dict]:
    """
    Download and parse GDC XML clinical files for TCGA-OV.
    
    Returns: Dict mapping patient_id → {platinum_response, raw_value, source_field}
    """
    print("\n" + "="*80)
    print("GDC XML CLINICAL FILE EXTRACTION")
    print("="*80)
    
    results = {}
    cache_dir = Path("data/validation/gdc_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        print("\n📥 Step 1: Querying GDC for TCGA-OV clinical supplement files...")
        
        with httpx.Client(timeout=120.0) as client:
            # Query for clinical supplement files
            r = client.post(
                f"{GDC_BASE}/files",
                json={
                    "filters": {
                        "op": "and",
                        "content": [
                            {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-OV"}},
                            {"op": "=", "content": {"field": "files.data_type", "value": "Clinical Supplement"}},
                            {"op": "=", "content": {"field": "files.data_format", "value": "XML"}}
                        ]
                    },
                    "size": 1000,
                    "format": "JSON"
                }
            )
            
            if r.status_code != 200:
                print(f"   ⚠️  GDC API returned {r.status_code}: {r.text[:200]}")
                return results
            
            files_data = r.json()
            files = files_data.get("data", {}).get("hits", [])
            print(f"   ✅ Found {len(files)} clinical supplement XML files")
            
            if len(files) == 0:
                print("   ⚠️  No XML files found - trying TSV format...")
                # Try TSV format
                r2 = client.post(
                    f"{GDC_BASE}/files",
                    json={
                        "filters": {
                            "op": "and",
                            "content": [
                                {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-OV"}},
                                {"op": "=", "content": {"field": "files.data_type", "value": "Clinical Supplement"}}
                            ]
                        },
                        "size": 1000
                    }
                )
                if r2.status_code == 200:
                    files_data2 = r2.json()
                    files = files_data2.get("data", {}).get("hits", [])
                    print(f"   ✅ Found {len(files)} clinical supplement files (any format)")
            
            print(f"\n📥 Step 2: Downloading and parsing {min(len(files), 100)} files...")
            
            for i, file_info in enumerate(files[:100]):  # Limit to 100 for speed
                if i % 20 == 0:
                    print(f"   Progress: {i}/{min(len(files), 100)}")
                
                file_id = file_info.get("id")
                file_name = file_info.get("file_name", "")
                case_id = file_info.get("cases", [{}])[0].get("id") if file_info.get("cases") else None
                submitter_id = file_info.get("cases", [{}])[0].get("submitter_id") if file_info.get("cases") else None
                
                if not file_id or not submitter_id:
                    continue
                
                # Check cache
                cache_file = cache_dir / f"{file_id}.xml"
                if cache_file.exists():
                    xml_path = cache_file
                else:
                    # Download file
                    xml_path = cache_dir / file_name
                    if not download_gdc_file(file_id, xml_path):
                        continue
                
                # Parse XML
                xml_data = parse_gdc_xml(xml_path)
                
                # Extract response
                response_value = None
                source_field = None
                
                for field, value in xml_data.items():
                    if value and isinstance(value, str):
                        field_lower = field.lower()
                        if any(term in field_lower for term in ["response", "outcome", "platinum", "therapy"]):
                            response_value = value
                            source_field = field
                            break
                
                if response_value:
                    # Import from parent directory
                    import sys
                    from pathlib import Path as PathLib
                    parent_dir = PathLib(__file__).parent.parent.parent
                    sys.path.insert(0, str(parent_dir))
                    from extract_platinum_response_labels import normalize_response
                    normalized = normalize_response(response_value)
                    if normalized != "unknown":
                        results[submitter_id] = {
                            "platinum_response": normalized,
                            "raw_response_value": response_value,
                            "source_field": source_field or "gdc_xml",
                            "file_id": file_id
                        }
        
        print(f"\n✅ GDC XML Extraction: Found {len(results)} patients with response labels")
        
    except Exception as e:
        print(f"   ❌ GDC XML extraction failed: {e}")
        import traceback
        traceback.print_exc()
    
    return results

if __name__ == "__main__":
    results = extract_platinum_response_from_gdc_xml()
    print(f"\n📊 Results: {len(results)} patients")
    for patient_id, data in list(results.items())[:5]:
        print(f"   {patient_id}: {data['platinum_response']} ({data['source_field']})")

