#!/usr/bin/env python3
"""
Guide for downloading PTRC-HGSOC data

Since TCIA requires authentication and manual download,
this script provides guidance and checks for downloaded files.
"""

import json
import sys
from pathlib import Path
from typing import List

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

OUTPUT_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw" / "ptrc_hgsoc"
CLINICAL_DIR = OUTPUT_DIR / "clinical"
MAF_DIR = OUTPUT_DIR / "maf"


def check_downloaded_files() -> dict:
    """Check what files have been downloaded."""
    status = {
        "clinical_files": [],
        "maf_files": [],
        "other_files": []
    }
    
    # Check clinical directory
    if CLINICAL_DIR.exists():
        status["clinical_files"] = [f.name for f in CLINICAL_DIR.iterdir() if f.is_file()]
    
    # Check MAF directory
    if MAF_DIR.exists():
        status["maf_files"] = [f.name for f in MAF_DIR.iterdir() if f.is_file()]
    
    # Check root directory
    if OUTPUT_DIR.exists():
        status["other_files"] = [
            f.name for f in OUTPUT_DIR.iterdir() 
            if f.is_file() and f.name != "README.md" and f.name != "DOWNLOAD_INSTRUCTIONS.md"
        ]
    
    return status


def print_download_guide():
    """Print download instructions."""
    print("="*80)
    print("PTRC-HGSOC DATA DOWNLOAD GUIDE")
    print("="*80)
    print()
    print("📋 STEP 1: Access TCIA Collection")
    print("   URL: https://www.cancerimagingarchive.net/collection/ptrc-hgsoc/")
    print()
    print("📋 STEP 2: Register/Login (if needed)")
    print("   Registration: https://www.cancerimagingarchive.net/registration/")
    print("   Requires: Email, institution, Data Use Agreement")
    print()
    print("📋 STEP 3: Download Files")
    print()
    print("   A. Clinical Data (CSV) - REQUIRED")
    print("      - Look for: 'Clinical data' → 'Download (45.31kb)'")
    print("      - Format: CSV")
    print("      - Save to:", CLINICAL_DIR)
    print()
    print("   B. Genomic Data (MAF/VCF) - REQUIRED")
    print("      - Look for: Genomic/mutation data download link")
    print("      - May be labeled as: 'MAF', 'VCF', 'Mutation', or 'Genomic'")
    print("      - Format: MAF (.maf) or VCF (.vcf) or TXT")
    print("      - Save to:", MAF_DIR)
    print("      - ⚠️  NOTE: If not visible, check:")
    print("         * Supplementary files section")
    print("         * Associated publication")
    print("         * Contact TCIA support")
    print()
    print("   C. Other Data (Optional)")
    print("      - Histopathology images (SVS) - not needed for our use case")
    print("      - Proteomic data - may be useful but not required")
    print()
    print("📋 STEP 4: Verify Files Downloaded")
    print("   Run this script again to check downloaded files")
    print()
    print("📋 STEP 5: Run Extraction")
    print("   python3 scripts/data_acquisition/extract_ptrc_hgsoc.py")
    print()


def main():
    """Check download status and provide guidance."""
    print("="*80)
    print("PTRC-HGSOC DOWNLOAD STATUS CHECK")
    print("="*80)
    print()
    
    # Check for downloaded files
    status = check_downloaded_files()
    
    print("📁 Current Status:")
    print()
    
    if status["clinical_files"]:
        print("✅ Clinical Files Found:")
        for f in status["clinical_files"]:
            print(f"   - {f}")
    else:
        print("❌ Clinical Files: NOT FOUND")
        print(f"   Expected location: {CLINICAL_DIR}")
    
    print()
    
    if status["maf_files"]:
        print("✅ MAF/Genomic Files Found:")
        for f in status["maf_files"]:
            print(f"   - {f}")
    else:
        print("❌ MAF/Genomic Files: NOT FOUND")
        print(f"   Expected location: {MAF_DIR}")
        print("   ⚠️  These files are REQUIRED for extraction")
    
    print()
    
    if status["other_files"]:
        print("📄 Other Files Found:")
        for f in status["other_files"]:
            print(f"   - {f}")
    
    print()
    print("="*80)
    
    # Check if ready for extraction
    has_clinical = len(status["clinical_files"]) > 0
    has_maf = len(status["maf_files"]) > 0
    
    if has_clinical and has_maf:
        print("✅ READY FOR EXTRACTION!")
        print()
        print("All required files found. You can now run:")
        print("   python3 scripts/data_acquisition/extract_ptrc_hgsoc.py")
    else:
        print("⚠️  NOT READY FOR EXTRACTION")
        print()
        if not has_clinical:
            print("❌ Missing: Clinical data files")
        if not has_maf:
            print("❌ Missing: Genomic/MAF files")
        print()
        print_download_guide()
    
    print()


if __name__ == "__main__":
    main()

