#!/usr/bin/env python3
"""
Download TCGA-OV Files Directly via GDC API
No gdc-client required - direct HTTP download
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import httpx
import hashlib
import time

DATA_DIR = Path("data/serial_sae/tcga_ov")
DOWNLOAD_DIR = DATA_DIR / "downloads"
DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

GDC_API_BASE = "https://api.gdc.cancer.gov"
GDC_DATA_BASE = "https://api.gdc.cancer.gov/data"

def calculate_md5(file_path: Path) -> str:
    """Calculate MD5 hash of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def download_file(
    file_id: str,
    file_name: str,
    output_dir: Path,
    expected_md5: Optional[str] = None,
    chunk_size: int = 8192
) -> bool:
    """Download a single file from GDC API."""
    output_path = output_dir / file_name
    
    # Skip if already downloaded and MD5 matches
    if output_path.exists() and expected_md5:
        actual_md5 = calculate_md5(output_path)
        if actual_md5 == expected_md5:
            print(f"  ✅ Already downloaded: {file_name}")
            return True
    
    try:
        # GDC API download endpoint
        url = f"{GDC_DATA_BASE}/{file_id}"
        
        with httpx.Client(timeout=300.0, follow_redirects=True) as client:
            with client.stream("GET", url) as response:
                response.raise_for_status()
                
                total_size = int(response.headers.get("content-length", 0))
                downloaded = 0
                
                with open(output_path, "wb") as f:
                    for chunk in response.iter_bytes(chunk_size=chunk_size):
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        # Progress indicator for large files
                        if total_size > 0 and downloaded % (10 * 1024 * 1024) == 0:
                            percent = (downloaded / total_size) * 100
                            print(f"    Progress: {percent:.1f}% ({downloaded / (1024*1024):.1f} MB / {total_size / (1024*1024):.1f} MB)")
                
                # Verify MD5 if provided
                if expected_md5:
                    actual_md5 = calculate_md5(output_path)
                    if actual_md5 != expected_md5:
                        print(f"    ⚠️  MD5 mismatch for {file_name}")
                        print(f"       Expected: {expected_md5}")
                        print(f"       Got: {actual_md5}")
                        return False
                
                print(f"  ✅ Downloaded: {file_name} ({downloaded / (1024*1024):.1f} MB)")
                return True
                
    except Exception as e:
        print(f"  ❌ Error downloading {file_name}: {e}")
        if output_path.exists():
            output_path.unlink()  # Remove partial file
        return False

def download_files_from_manifest(manifest_file: Path, output_dir: Path, max_files: Optional[int] = None):
    """Download files listed in GDC manifest."""
    print(f"📋 Loading manifest: {manifest_file}")
    df = pd.read_csv(manifest_file)
    
    total_files = len(df)
    if max_files:
        df = df.head(max_files)
        print(f"📊 Downloading first {len(df)} of {total_files} files")
    else:
        print(f"📊 Downloading {total_files} files")
    
    total_size = df['file_size'].sum() / (1024**3)
    print(f"📦 Total size: {total_size:.2f} GB")
    print()
    
    # Create subdirectories for organization
    rnaseq_dir = output_dir / "rna_seq"
    rnaseq_dir.mkdir(exist_ok=True)
    
    successful = 0
    failed = 0
    
    for idx, row in df.iterrows():
        file_id = row['file_id']
        file_name = row['file_name']
        file_size_mb = row['file_size'] / (1024**2)
        
        print(f"[{idx+1}/{len(df)}] {file_name} ({file_size_mb:.1f} MB)")
        
        # Determine output directory based on file type
        if 'rna_seq' in file_name.lower() or 'expression' in file_name.lower():
            output_subdir = rnaseq_dir
        else:
            output_subdir = output_dir
        
        # Download file
        if download_file(file_id, file_name, output_subdir):
            successful += 1
        else:
            failed += 1
        
        # Rate limiting - be nice to the API
        time.sleep(0.5)
    
    print()
    print("="*80)
    print("Download Summary")
    print("="*80)
    print(f"✅ Successful: {successful}")
    print(f"❌ Failed: {failed}")
    print(f"📁 Output directory: {output_dir}")
    print()

def main():
    print("="*80)
    print("TCGA-OV Direct API Download")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    manifest_file = DATA_DIR / "tcga_ov_rnaseq_manifest.csv"
    
    if not manifest_file.exists():
        print(f"❌ Error: Manifest file not found: {manifest_file}")
        print("Run: python3 scripts/serial_sae/download_tcga_ov_gdc.py first")
        return
    
    # Check if user wants to limit files (for testing)
    import sys
    max_files = None
    if len(sys.argv) > 1:
        try:
            max_files = int(sys.argv[1])
            print(f"⚠️  Limiting to first {max_files} files (for testing)")
        except ValueError:
            print(f"⚠️  Invalid number: {sys.argv[1]}, downloading all files")
    
    # Download files
    download_files_from_manifest(manifest_file, DOWNLOAD_DIR, max_files)
    
    print("="*80)
    print("✅ Download complete")
    print("="*80)
    print()
    print("Next steps:")
    print("1. Process RNA-seq files: scripts/serial_sae/process_tcga_ov_rnaseq.py")
    print("2. Compute SAE pathway scores")
    print("3. Validate against clinical outcomes")

if __name__ == "__main__":
    main()
