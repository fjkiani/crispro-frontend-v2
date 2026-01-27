#!/usr/bin/env python3
"""
Download TCGA-OV Data from GDC Data Portal
Direct download of RNA-seq, mutations, and clinical data
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import subprocess
import sys

DATA_DIR = Path("data/serial_sae/tcga_ov")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# GDC Data Portal API
GDC_API_BASE = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-OV"

def check_gdc_client():
    """Check if GDC Data Transfer Tool is installed."""
    try:
        result = subprocess.run(
            ["gdc-client", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    return False

def download_via_api():
    """Download TCGA-OV data via GDC API."""
    print("="*80)
    print("Download TCGA-OV Data from GDC Data Portal")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    import httpx
    
    # Query for TCGA-OV files
    print(f"🔍 Querying GDC for TCGA-OV files...")
    
    # Get file IDs for RNA-seq
    rna_query = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "cases.project.project_id",
                        "value": [PROJECT_ID]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "files.data_type",
                        "value": ["Gene Expression Quantification"]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "files.experimental_strategy",
                        "value": ["RNA-Seq"]
                    }
                }
            ]
        },
        "format": "JSON",
        "size": "10000"
    }
    
    try:
        with httpx.Client(timeout=120.0) as client:
            # Query files
            r = client.post(
                f"{GDC_API_BASE}/files",
                json=rna_query,
                params={"pretty": "true"}
            )
            r.raise_for_status()
            rna_files = r.json()
            
            print(f"✅ Found {rna_files.get('data', {}).get('pagination', {}).get('count', 0)} RNA-seq files")
            
            # Get mutations
            mut_query = {
                "filters": {
                    "op": "and",
                    "content": [
                        {
                            "op": "in",
                            "content": {
                                "field": "cases.project.project_id",
                                "value": [PROJECT_ID]
                            }
                        },
                        {
                            "op": "in",
                            "content": {
                                "field": "files.data_type",
                                "value": ["Masked Somatic Mutation"]
                            }
                        }
                    ]
                },
                "format": "JSON",
                "size": "10000"
            }
            
            r2 = client.post(
                f"{GDC_API_BASE}/files",
                json=mut_query,
                params={"pretty": "true"}
            )
            r2.raise_for_status()
            mut_files = r2.json()
            
            print(f"✅ Found {mut_files.get('data', {}).get('pagination', {}).get('count', 0)} mutation files")
            
            # Get clinical data
            clinical_query = {
                "filters": {
                    "op": "and",
                    "content": [
                        {
                            "op": "in",
                            "content": {
                                "field": "cases.project.project_id",
                                "value": [PROJECT_ID]
                            }
                        },
                        {
                            "op": "in",
                            "content": {
                                "field": "files.data_category",
                                "value": ["Clinical"]
                            }
                        }
                    ]
                },
                "format": "JSON",
                "size": "10000"
            }
            
            r3 = client.post(
                f"{GDC_API_BASE}/files",
                json=clinical_query,
                params={"pretty": "true"}
            )
            r3.raise_for_status()
            clinical_files = r3.json()
            
            print(f"✅ Found {clinical_files.get('data', {}).get('pagination', {}).get('count', 0)} clinical files")
            
            # Save file manifests
            if rna_files.get("data", {}).get("hits"):
                rna_manifest = []
                for hit in rna_files["data"]["hits"]:
                    rna_manifest.append({
                        "file_id": hit.get("id"),
                        "file_name": hit.get("file_name"),
                        "file_size": hit.get("file_size"),
                        "cases": [c.get("case_id") for c in hit.get("cases", [])]
                    })
                
                manifest_df = pd.DataFrame(rna_manifest)
                manifest_file = DATA_DIR / "tcga_ov_rnaseq_manifest.csv"
                manifest_df.to_csv(manifest_file, index=False)
                print(f"💾 Saved RNA-seq manifest: {manifest_file}")
            
            # Save summary
            summary = {
                "project": PROJECT_ID,
                "download_date": datetime.now().isoformat(),
                "rna_seq_files": rna_files.get("data", {}).get("pagination", {}).get("count", 0),
                "mutation_files": mut_files.get("data", {}).get("pagination", {}).get("count", 0),
                "clinical_files": clinical_files.get("data", {}).get("pagination", {}).get("count", 0),
                "note": "Use GDC Data Transfer Tool (gdc-client) to download actual files",
                "instructions": "Run: gdc-client download -m tcga_ov_rnaseq_manifest.csv"
            }
            
            summary_file = DATA_DIR / "download_summary.json"
            with open(summary_file, "w") as f:
                json.dump(summary, f, indent=2)
            print(f"💾 Saved summary: {summary_file}")
            
            print("\n" + "="*80)
            print("✅ File manifest created")
            print("="*80)
            print("\nNext steps:")
            print("1. Install GDC Data Transfer Tool: pip install gdc-client")
            print("2. Download files: gdc-client download -m data/serial_sae/tcga_ov/tcga_ov_rnaseq_manifest.csv")
            print("3. Process downloaded files for SAE computation")
            
    except Exception as e:
        print(f"❌ Error querying GDC API: {e}")
        print("\nAlternative: Download manually from https://portal.gdc.cancer.gov/")
        print("  - Select project: TCGA-OV")
        print("  - Add files: RNA-Seq, Mutations, Clinical")
        print("  - Download manifest and use gdc-client")

def main():
    print("Checking for GDC Data Transfer Tool...")
    if check_gdc_client():
        print("✅ GDC client found")
    else:
        print("⚠️  GDC client not found")
        print("   Install with: pip install gdc-client")
        print("   Or download manually from: https://portal.gdc.cancer.gov/")
        print()
    
    download_via_api()

if __name__ == "__main__":
    main()
