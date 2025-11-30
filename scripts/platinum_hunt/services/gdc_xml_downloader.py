#!/usr/bin/env python3
"""
GDC XML File Downloader - Aggressive Extraction

Actually downloads and parses GDC XML clinical supplement files.
"""

import httpx
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional
import tempfile
import zipfile
import asyncio

GDC_BASE = "https://api.gdc.cancer.gov"
GDC_DATA_BASE = "https://api.gdc.cancer.gov/data"

async def query_gdc_clinical_files() -> List[Dict]:
    """Query GDC for TCGA-OV clinical supplement files."""
    print("\n📥 Querying GDC for TCGA-OV clinical supplement files...")
    
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": ["TCGA-OV"]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["Clinical"]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Clinical Supplement"]}}
        ]
    }
    
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,file_size,cases.submitter_id",
        "format": "json",
        "size": 1000
    }
    
    try:
        async with httpx.AsyncClient(timeout=300.0) as client:
            r = await client.get(f"{GDC_BASE}/files", params=params)
            r.raise_for_status()
            data = r.json()
            
            if "data" in data and "hits" in data["data"]:
                files = data["data"]["hits"]
                print(f"   ✅ Found {len(files)} clinical supplement files")
                return files
            return []
    except Exception as e:
        print(f"   ❌ Error querying GDC: {e}")
        return []

async def download_and_parse_gdc_file(file_id: str, file_name: str, cache_dir: Path) -> Dict[str, str]:
    """Download and parse a single GDC XML file."""
    cache_file = cache_dir / f"{file_id}.xml"
    
    # Check cache first
    if cache_file.exists():
        print(f"   📂 Using cached: {file_name}")
        xml_path = cache_file
    else:
        print(f"   ⬇️  Downloading: {file_name} ({file_id})")
        try:
            async with httpx.AsyncClient(timeout=300.0, follow_redirects=True) as client:
                r = await client.get(f"{GDC_DATA_BASE}/{file_id}")
                r.raise_for_status()
                
                cache_dir.mkdir(parents=True, exist_ok=True)
                cache_file.write_bytes(r.content)
                xml_path = cache_file
                print(f"   ✅ Downloaded {len(r.content)} bytes")
        except Exception as e:
            print(f"   ❌ Download failed: {e}")
            return {}
    
    # Parse XML
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        # Extract patient ID from XML (handle namespaces)
        patient_id = None
        
        # Try to find patient ID in various namespace formats
        for elem in root.iter():
            # Check tag name (ignoring namespace)
            tag_name = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
            if tag_name in ["patient_id", "submitter_id", "bcr_patient_barcode"]:
                if elem.text and elem.text.strip():
                    patient_id = elem.text.strip()
                    break
        
        # If not found, try to extract from file name
        if not patient_id:
            import re
            match = re.search(r"TCGA-[\w-]+", file_name)
            if match:
                patient_id = match.group(0)
        
        # Look for response fields (handle namespaces)
        response_value = None
        response_field_names = [
            "primary_therapy_outcome_success",  # Most common in TCGA
            "best_overall_response",
            "response_to_treatment",
            "treatment_outcome",
            "platinum_status",
            "response"
        ]
        
        for field_name in response_field_names:
            for elem in root.iter():
                # Check tag name (ignoring namespace)
                tag_name = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
                if field_name.lower() in tag_name.lower():
                    if elem.text and elem.text.strip():
                        response_value = elem.text.strip()
                        break
            if response_value:
                break
        
        if patient_id and response_value:
            return {
                "patient_id": patient_id,
                "response": response_value,
                "source": "GDC_XML",
                "file_id": file_id
            }
    except Exception as e:
        print(f"   ❌ Parse error for {file_name}: {e}")
    
    return {}

async def extract_platinum_response_from_gdc_xml_aggressive() -> Dict[str, Dict]:
    """Aggressively extract from GDC XML files."""
    print("\n" + "="*80)
    print("GDC XML AGGRESSIVE EXTRACTION")
    print("="*80)
    
    cache_dir = Path("data/validation/gdc_xml_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Query for files
    files = await query_gdc_clinical_files()
    
    if not files:
        print("   ⚠️  No clinical supplement files found")
        return results
    
    # Download and parse ALL files (aggressive extraction)
    print(f"\n   Processing ALL {len(files)} files (aggressive extraction)...")
    tasks = []
    for i, file_info in enumerate(files):
        file_id = file_info.get("file_id")
        file_name = file_info.get("file_name", "unknown")
        tasks.append(download_and_parse_gdc_file(file_id, file_name, cache_dir))
    
    # Process in batches of 10
    batch_size = 10
    for i in range(0, len(tasks), batch_size):
        batch = tasks[i:i+batch_size]
        batch_results = await asyncio.gather(*batch, return_exceptions=True)
        
        for result in batch_results:
            if isinstance(result, dict) and result:
                patient_id = result.get("patient_id")
                if patient_id:
                    results[patient_id] = result
        
        progress = min(i+batch_size, len(tasks))
        print(f"   Progress: {progress}/{len(tasks)} files processed, {len(results)} patients found")
        # REMOVED EARLY EXIT - Process ALL 597 files to maximize overlap with Zo's dataset
        # Need ≥40 patients for statistical validation, so we need maximum coverage
    
    print(f"\n   ✅ Found {len(results)} patients with response data from GDC XML")
    return results

if __name__ == "__main__":
    import asyncio
    results = asyncio.run(extract_platinum_response_from_gdc_xml_aggressive())
    print(f"\nTotal patients found: {len(results)}")

