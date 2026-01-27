#!/usr/bin/env python3
"""
Download GSE165897 Dataset (11 Paired scRNA-seq Samples)
Single-cell RNA-seq data for treatment-naïve vs post-NACT ovarian cancer
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
import sys

# Try to import GEOquery (R package wrapper) or use direct download
try:
    import GEOparse
    HAS_GEO = True
except ImportError:
    HAS_GEO = False
    print("⚠️  GEOparse not installed. Install with: pip install GEOparse")
    print("   Or download manually from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897")

DATA_DIR = Path("data/serial_sae/gse165897")
DATA_DIR.mkdir(parents=True, exist_ok=True)

GSE_ID = "GSE165897"

def download_via_geoparse():
    """Download GSE165897 using GEOparse."""
    print(f"📥 Downloading {GSE_ID} via GEOparse...")
    
    gse = GEOparse.get_GEO(geo=GSE_ID, destdir=str(DATA_DIR))
    
    print(f"✅ Downloaded {GSE_ID}")
    print(f"  Platforms: {len(gse.platforms)}")
    print(f"  Samples: {len(gse.gsms)}")
    
    # Extract sample metadata
    sample_metadata = []
    for gsm_name, gsm in gse.gsms.items():
        sample_metadata.append({
            "sample_id": gsm_name,
            "title": gsm.metadata.get("title", [""])[0],
            "characteristics": gsm.metadata.get("characteristics_ch1", []),
            "treatment": "treatment-naïve" if "treatment-naïve" in str(gsm.metadata).lower() else "post-NACT",
            "patient_id": None  # Extract from characteristics
        })
    
    # Save metadata
    metadata_df = pd.DataFrame(sample_metadata)
    metadata_df.to_csv(DATA_DIR / "sample_metadata.csv", index=False)
    print(f"💾 Saved sample metadata: {DATA_DIR / 'sample_metadata.csv'}")
    
    # Extract expression data (if available as matrix)
    if hasattr(gse, "table"):
        gse.table.to_csv(DATA_DIR / "expression_matrix.csv")
        print(f"💾 Saved expression matrix: {DATA_DIR / 'expression_matrix.csv'}")
    
    return gse

def download_via_manual_instructions():
    """Provide manual download instructions."""
    print("="*80)
    print("MANUAL DOWNLOAD INSTRUCTIONS FOR GSE165897")
    print("="*80)
    print()
    print("1. Navigate to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897")
    print("2. Click 'Download family' → 'Series Matrix File(s)'")
    print("3. Download: GSE165897_series_matrix.txt.gz")
    print("4. Save to:", DATA_DIR)
    print("5. Extract: gunzip GSE165897_series_matrix.txt.gz")
    print()
    print("OR")
    print()
    print("1. Navigate to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897")
    print("2. Click 'SRA Run Selector'")
    print("3. Download FASTQ files for 11 samples")
    print("4. Process with 10x Genomics Cell Ranger")
    print()
    print("="*80)

def identify_paired_samples(metadata_df: pd.DataFrame) -> dict:
    """Identify 11 paired patients (treatment-naïve + post-NACT)."""
    print("🔍 Identifying paired samples...")
    
    # Group by patient ID (extract from characteristics or title)
    # Expected format: Patient ID in title or characteristics
    paired_patients = {}
    
    # Try to extract patient ID from title or characteristics
    for _, row in metadata_df.iterrows():
        title = str(row.get("title", ""))
        characteristics = str(row.get("characteristics", ""))
        
        # Extract patient ID (heuristic: look for "P" followed by number)
        import re
        patient_match = re.search(r'[Pp](\d+)', title + " " + characteristics)
        if patient_match:
            patient_id = f"P{patient_match.group(1)}"
        else:
            # Fallback: use sample ID as patient ID (if only 2 samples per patient)
            patient_id = row["sample_id"].split("_")[0] if "_" in row["sample_id"] else row["sample_id"]
        
        treatment = row.get("treatment", "unknown")
        
        if patient_id not in paired_patients:
            paired_patients[patient_id] = {
                "treatment_naive": None,
                "post_nact": None
            }
        
        if treatment == "treatment-naïve":
            paired_patients[patient_id]["treatment_naive"] = row["sample_id"]
        elif treatment == "post-NACT":
            paired_patients[patient_id]["post_nact"] = row["sample_id"]
    
    # Filter to patients with both timepoints
    complete_pairs = {
        pid: data for pid, data in paired_patients.items()
        if data["treatment_naive"] and data["post_nact"]
    }
    
    print(f"✅ Identified {len(complete_pairs)} paired patients")
    return complete_pairs

def main():
    print("="*80)
    print("Download GSE165897 Dataset (11 Paired scRNA-seq Samples)")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    if HAS_GEO:
        try:
            gse = download_via_geoparse()
            
            # Load metadata
            metadata_df = pd.read_csv(DATA_DIR / "sample_metadata.csv")
            
            # Identify paired samples
            paired_patients = identify_paired_samples(metadata_df)
            
            # Save paired patients list
            output = {
                "gse_id": GSE_ID,
                "total_samples": len(metadata_df),
                "paired_patients": paired_patients,
                "download_date": datetime.now().isoformat(),
                "data_source": "GEO",
                "data_type": "scRNA-seq",
                "notes": "11 paired patients (treatment-naïve vs post-NACT)"
            }
            
            output_file = DATA_DIR / "paired_patients.json"
            with open(output_file, "w") as f:
                json.dump(output, f, indent=2)
            print(f"💾 Saved paired patients list: {output_file}")
            
            print("\n" + "="*80)
            print("✅ Download complete")
            print("="*80)
            print("\nNext step: Process scRNA-seq data with compute_serial_sae.py")
            
        except Exception as e:
            print(f"❌ Error downloading via GEOparse: {e}")
            print("\nFalling back to manual download instructions...\n")
            download_via_manual_instructions()
    else:
        download_via_manual_instructions()

if __name__ == "__main__":
    main()
