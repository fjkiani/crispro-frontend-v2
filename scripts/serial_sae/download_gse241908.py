#!/usr/bin/env python3
"""
Download GSE241908 (Platinum-resistant ovarian cancer, pre/post bevacizumab)
7 patients, 14 RNA-seq samples (7 paired)

Access: Public GEO
Study: Bevacizumab response in platinum-resistant ovarian cancer
"""

import os
import sys
import requests
import gzip
from pathlib import Path
from typing import Dict, List
import json
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Output directory
OUTPUT_DIR = Path("data/serial_sae/gse241908")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GEO_BASE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241908/suppl/"

def get_geo_metadata(accession: str) -> Dict:
    """Get GEO metadata for accession."""
    url = f"{GEO_BASE_URL}?acc={accession}&targ=self&form=text&view=full"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text

def download_geo_file(filename: str, output_path: Path) -> bool:
    """Download a file from GEO FTP."""
    url = f"{GEO_FTP_BASE}{filename}"
    print(f"  Downloading {filename}...")
    
    try:
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"  ✅ Downloaded: {output_path}")
        return True
    except Exception as e:
        print(f"  ❌ Failed to download {filename}: {e}")
        return False

def main():
    """Download GSE241908 data."""
    print("="*80)
    print("GSE241908 Download: Platinum-Resistant Ovarian Cancer + Bevacizumab")
    print("="*80)
    print("\nCohort: 7 patients, 14 RNA-seq samples (7 paired)")
    print("Treatment: Bevacizumab (anti-VEGF) for platinum-resistant disease")
    print("Access: Public GEO\n")
    
    # Expected files (actual GEO files)
    expected_files = [
        "GSE241908_data1_TPM.csv.gz",  # Expression matrix (TPM values)
        "GSE241908_series_matrix.txt.gz",  # Series matrix (if available)
    ]
    
    downloaded_files = []
    
    for filename in expected_files:
        output_path = OUTPUT_DIR / filename
        if output_path.exists():
            print(f"  ⏭️  Already exists: {filename}")
            downloaded_files.append(str(output_path))
        else:
            if download_geo_file(filename, output_path):
                downloaded_files.append(str(output_path))
    
    # Create summary
    summary = {
        "accession": "GSE241908",
        "cohort": "Platinum-resistant ovarian cancer + bevacizumab",
        "patients": 7,
        "samples": 14,
        "paired": True,
        "treatment": "Bevacizumab (anti-VEGF)",
        "data_type": "Bulk RNA-seq",
        "files_downloaded": downloaded_files,
        "date_downloaded": datetime.now().isoformat(),
        "note": "Bevacizumab-specific cohort, can validate VEGF pathway findings from GSE165897"
    }
    
    summary_path = OUTPUT_DIR / "download_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n💾 Download summary: {summary_path}")
    print(f"\n✅ Downloaded {len(downloaded_files)} files")
    print("\n📋 Next steps:")
    print("  1. Extract expression matrix from GSE241908_series_matrix.txt.gz")
    print("  2. Process into pathway scores (pre/post bevacizumab)")
    print("  3. Compute pathway kinetics (ΔVEGF, ΔDDR, etc.)")
    print("  4. Correlate with bevacizumab response")

if __name__ == "__main__":
    main()
