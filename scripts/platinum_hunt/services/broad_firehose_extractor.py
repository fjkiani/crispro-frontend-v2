#!/usr/bin/env python3
"""
Broad Firehose TCGA Clinical Data Extractor

Downloads and parses Broad Firehose TCGA clinical data files.
"""

import httpx
import pandas as pd
from pathlib import Path
from typing import Dict
import tempfile
import zipfile

BROAD_FIREHOSE_BASE = "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/OV/20160128"

def download_broad_firehose_clinical() -> Path:
    """
    Download Broad Firehose TCGA-OV clinical data.
    
    Returns: Path to downloaded file
    """
    print("\n📥 Downloading Broad Firehose TCGA-OV clinical data...")
    
    # Broad Firehose file structure
    # Clinical file: gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz
    clinical_file_url = f"{BROAD_FIREHOSE_BASE}/gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz"
    
    cache_dir = Path("data/validation/broad_firehose_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    output_file = cache_dir / "OV.Clinical_Pick_Tier1.tar.gz"
    
    if output_file.exists():
        print(f"   ✅ Using cached file: {output_file}")
        return output_file
    
    try:
        with httpx.Client(timeout=300.0, follow_redirects=True) as client:
            r = client.get(clinical_file_url)
            r.raise_for_status()
            
            output_file.write_bytes(r.content)
            print(f"   ✅ Downloaded {len(r.content) / 1024 / 1024:.1f} MB")
            return output_file
            
    except Exception as e:
        print(f"   ❌ Download failed: {e}")
        return None

def parse_broad_firehose_clinical(tar_path: Path) -> Dict[str, Dict]:
    """
    Parse Broad Firehose clinical tar.gz file.
    
    Returns: Dict mapping patient_id → {platinum_response, raw_value, source_field}
    """
    results = {}
    
    if not tar_path or not tar_path.exists():
        return results
    
    print(f"\n📥 Extracting and parsing {tar_path.name}...")
    
    try:
        import tarfile
        
        with tarfile.open(tar_path, 'r:gz') as tar:
            # Find clinical data file
            clinical_files = [m for m in tar.getmembers() if 'clinical' in m.name.lower() and m.name.endswith('.txt')]
            
            if not clinical_files:
                print("   ⚠️  No clinical .txt files found in archive")
                return results
            
            # Extract first clinical file
            clinical_file = clinical_files[0]
            extracted = tar.extractfile(clinical_file)
            
            # Parse as TSV - Broad Firehose is TRANSPOSED
            # Row 0: "Hybridization REF" + patient IDs
            # Row 1: "Composite Element REF" + "value"
            # Row 2+: Field names in first column, values in subsequent columns
            
            df = pd.read_csv(extracted, sep='\t', low_memory=False, header=None)
            print(f"   ✅ Parsed {len(df)} rows, {len(df.columns)} columns")
            
            # First column contains field names (skip row 0 and 1)
            field_names = df.iloc[2:, 0].astype(str).tolist()
            
            # Find response-related fields
            response_field_indices = []
            for idx, field_name in enumerate(field_names):
                if any(term in field_name.lower() for term in 
                      ["response", "outcome", "platinum", "therapy", "treatment", 
                       "recurrence", "progression", "resistant", "sensitive"]):
                    response_field_indices.append(idx + 2)  # +2 because we skipped rows 0-1
            
            print(f"   ✅ Found {len(response_field_indices)} response-related fields")
            if response_field_indices:
                print(f"      Fields: {[field_names[i-2] for i in response_field_indices[:5]]}")
            
            # Patient IDs are in row 0, columns 1+ (skip column 0 which is "Hybridization REF")
            patient_ids = df.iloc[0, 1:].astype(str).tolist()
            print(f"   ✅ Found {len(patient_ids)} patient IDs")
            
            # Extract responses for each patient
            for col_idx, patient_id in enumerate(patient_ids):
                if not patient_id or not patient_id.startswith("TCGA-") or patient_id == "value":
                    continue
                
                # Extract TCGA patient ID (remove sample suffix if present)
                if len(patient_id) > 12:
                    patient_id = patient_id[:12]
                
                # Check response fields for this patient
                response_value = None
                source_field = None
                
                for field_idx in response_field_indices:
                    value = df.iloc[field_idx, col_idx + 1]  # +1 because column 0 is field names
                    if pd.notna(value) and str(value).strip() and str(value).strip().upper() != "NA":
                        response_value = str(value).strip()
                        source_field = field_names[field_idx - 2]  # -2 because we skipped rows 0-1
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
                        results[patient_id] = {
                            "platinum_response": normalized,
                            "raw_response_value": response_value,
                            "source_field": source_field or "broad_firehose"
                        }
            
            print(f"   ✅ Extracted {len(results)} patients with response labels")
            
    except Exception as e:
        print(f"   ❌ Parsing failed: {e}")
        import traceback
        traceback.print_exc()
    
    return results

def extract_platinum_response_from_broad_firehose() -> Dict[str, Dict]:
    """
    Main extraction function from Broad Firehose.
    
    Returns: Dict mapping patient_id → {platinum_response, raw_value, source_field}
    """
    print("\n" + "="*80)
    print("BROAD FIREHOSE EXTRACTION")
    print("="*80)
    
    results = {}
    
    tar_path = download_broad_firehose_clinical()
    if tar_path:
        results = parse_broad_firehose_clinical(tar_path)
    
    print(f"\n✅ Broad Firehose: Found {len(results)} patients with response labels")
    return results

if __name__ == "__main__":
    results = extract_platinum_response_from_broad_firehose()
    print(f"\n📊 Results: {len(results)} patients")
    for patient_id, data in list(results.items())[:5]:
        print(f"   {patient_id}: {data['platinum_response']} ({data['source_field']})")

