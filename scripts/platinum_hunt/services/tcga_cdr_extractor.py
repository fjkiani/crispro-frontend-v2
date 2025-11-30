#!/usr/bin/env python3
"""
TCGA-CDR (Pan-Cancer Clinical Data Resource) Extractor

TCGA-CDR has standardized clinical data including response fields.
"""

import httpx
import pandas as pd
from pathlib import Path
from typing import Dict
import json

# TCGA-CDR download URL (from Liu et al. 2018 paper)
TCGA_CDR_URL = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7ecfae82fe72"

def download_tcga_cdr() -> Path:
    """Download TCGA-CDR clinical data file."""
    print("\n📥 Downloading TCGA-CDR clinical data...")
    
    cache_dir = Path("data/validation/tcga_cdr_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "TCGA-CDR-SupplementTableS1.xlsx"
    
    if cache_file.exists():
        print(f"   ✅ Using cached file: {cache_file}")
        return cache_file
    
    try:
        with httpx.Client(timeout=300.0, follow_redirects=True) as client:
            r = client.get(TCGA_CDR_URL)
            r.raise_for_status()
            
            cache_file.write_bytes(r.content)
            print(f"   ✅ Downloaded {len(r.content)} bytes")
            return cache_file
    except Exception as e:
        print(f"   ❌ Download failed: {e}")
        return None

def extract_platinum_response_from_tcga_cdr() -> Dict[str, Dict]:
    """Extract platinum response from TCGA-CDR."""
    print("\n" + "="*80)
    print("TCGA-CDR EXTRACTION")
    print("="*80)
    
    results = {}
    
    # Download file
    file_path = download_tcga_cdr()
    if not file_path:
        return results
    
    try:
        # Read Excel file
        print("   📖 Reading TCGA-CDR Excel file...")
        df = pd.read_excel(file_path, sheet_name=0)  # First sheet
        
        print(f"   ✅ Loaded {len(df)} rows, {len(df.columns)} columns")
        print(f"   Columns: {list(df.columns)[:20]}...")
        
        # Find patient ID column
        patient_id_col = None
        for col in df.columns:
            if "patient" in col.lower() or "bcr_patient_barcode" in col.lower() or "submitter_id" in col.lower():
                patient_id_col = col
                break
        
        if not patient_id_col:
            print("   ⚠️  Could not find patient ID column")
            return results
        
        # Find response columns
        response_cols = [c for c in df.columns if any(term in c.lower() for term in 
                      ["response", "outcome", "platinum", "best_overall", "treatment_outcome"])]
        
        print(f"   ✅ Found {len(response_cols)} response-related columns: {response_cols}")
        
        # Extract data
        for idx, row in df.iterrows():
            patient_id = str(row[patient_id_col]).strip()
            
            # Extract TCGA-XX-XXXX format
            if not patient_id.startswith("TCGA-"):
                continue
            
            # Get response from first available response column
            response_value = None
            response_col = None
            for col in response_cols:
                value = row[col]
                if pd.notna(value) and str(value).strip():
                    response_value = str(value).strip()
                    response_col = col
                    break
            
            if response_value:
                # Normalize response
                import sys
                from pathlib import Path as PathLib
                parent_dir = PathLib(__file__).parent.parent.parent
                sys.path.insert(0, str(parent_dir))
                from extract_platinum_response_labels import normalize_response
                
                normalized = normalize_response(response_value)
                if normalized != "unknown":
                    results[patient_id] = {
                        "response": normalized,
                        "source": "TCGA_CDR",
                        "original_value": response_value,
                        "source_field": response_col
                    }
        
        print(f"\n   ✅ Found {len(results)} patients with response data from TCGA-CDR")
        
    except Exception as e:
        print(f"   ❌ Error parsing TCGA-CDR: {e}")
        import traceback
        traceback.print_exc()
    
    return results

if __name__ == "__main__":
    results = extract_platinum_response_from_tcga_cdr()
    print(f"\nTotal patients found: {len(results)}")






