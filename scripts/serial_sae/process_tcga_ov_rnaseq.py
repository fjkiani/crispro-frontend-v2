#!/usr/bin/env python3
"""
Process TCGA-OV RNA-seq Files
Extract expression matrices and prepare for SAE computation
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import json
import gzip

DATA_DIR = Path("data/serial_sae/tcga_ov")
DOWNLOAD_DIR = DATA_DIR / "downloads" / "rna_seq"
OUTPUT_DIR = DATA_DIR / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_rnaseq_file(file_path: Path) -> Optional[pd.DataFrame]:
    """Parse a TCGA RNA-seq file (TSV format with comment header)."""
    try:
        # TCGA files have comment lines starting with #, then header row
        # Format: gene_id, gene_name, gene_type, unstranded, stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded
        
        # Read file, skipping comment lines
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find first non-comment line (header)
        header_idx = 0
        for i, line in enumerate(lines):
            if not line.strip().startswith('#'):
                header_idx = i
                break
        
        # Read from header onwards
        df = pd.read_csv(file_path, sep="\t", skiprows=header_idx, low_memory=False)
        
        # Skip rows that start with 'N_' (unmapped, multimapping, etc.)
        df = df[~df.iloc[:, 0].astype(str).str.startswith('N_')]
        
        # Use gene_name for matching (Hugo symbols)
        if 'gene_name' not in df.columns:
            print(f"  ⚠️  No gene_name column in {file_path.name}")
            return None
        
        # Use TPM for expression (tpm_unstranded)
        if 'tpm_unstranded' in df.columns:
            expr_col = 'tpm_unstranded'
        elif 'fpkm_unstranded' in df.columns:
            expr_col = 'fpkm_unstranded'
        elif 'fpkm_uq_unstranded' in df.columns:
            expr_col = 'fpkm_uq_unstranded'
        else:
            print(f"  ⚠️  No expression column found in {file_path.name}")
            return None
        
        # Extract gene name and expression
        result = pd.DataFrame({
            'gene_name': df['gene_name'].astype(str),
            'expression': pd.to_numeric(df[expr_col], errors='coerce')
        })
        
        # Remove NaN values and empty gene names
        result = result.dropna()
        result = result[result['gene_name'].str.strip() != '']
        result = result[result['gene_name'] != 'nan']
        
        # Set gene_name as index
        result = result.set_index('gene_name')
        
        return result
        
    except Exception as e:
        print(f"  ❌ Error parsing {file_path.name}: {e}")
        return None

def process_all_rnaseq_files(download_dir: Path) -> pd.DataFrame:
    """Process all RNA-seq files and create combined expression matrix."""
    print("="*80)
    print("Process TCGA-OV RNA-seq Files")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Find all RNA-seq files
    rnaseq_files = list(download_dir.glob("*.tsv")) + list(download_dir.glob("*.txt"))
    
    if not rnaseq_files:
        print(f"❌ No RNA-seq files found in {download_dir}")
        return pd.DataFrame()
    
    print(f"📊 Found {len(rnaseq_files)} RNA-seq files")
    print()
    
    # Process each file
    all_expressions = {}
    sample_ids = []
    
    for i, file_path in enumerate(rnaseq_files):
        print(f"[{i+1}/{len(rnaseq_files)}] Processing {file_path.name}...")
        
        # Extract sample ID from filename (TCGA format)
        sample_id = file_path.stem.split('.')[0] if '.' in file_path.name else file_path.stem
        
        # Parse file
        df = parse_rnaseq_file(file_path)
        
        if df is not None and not df.empty:
            # df already has gene_name as index and expression as column
            all_expressions[sample_id] = df['expression']
            sample_ids.append(sample_id)
            print(f"  ✅ Processed: {len(df)} genes")
        else:
            print(f"  ⚠️  Skipped: No valid data")
    
    if not all_expressions:
        print("❌ No valid expression data extracted")
        return pd.DataFrame()
    
    # Create combined expression matrix
    print()
    print("📊 Creating combined expression matrix...")
    expression_matrix = pd.DataFrame(all_expressions)
    
    # Fill NaN with 0 (or use imputation if preferred)
    expression_matrix = expression_matrix.fillna(0)
    
    print(f"✅ Expression matrix: {len(expression_matrix)} genes × {len(expression_matrix.columns)} samples")
    
    # Save expression matrix
    output_file = OUTPUT_DIR / "tcga_ov_expression_matrix.csv"
    expression_matrix.to_csv(output_file)
    print(f"💾 Saved: {output_file}")
    
    # Save sample metadata
    sample_metadata = pd.DataFrame({
        'sample_id': sample_ids,
        'file_count': len(sample_ids)
    })
    metadata_file = OUTPUT_DIR / "tcga_ov_sample_metadata.csv"
    sample_metadata.to_csv(metadata_file, index=False)
    print(f"💾 Saved: {metadata_file}")
    
    return expression_matrix

def main():
    if not DOWNLOAD_DIR.exists():
        print(f"❌ Error: Download directory not found: {DOWNLOAD_DIR}")
        print("Run: python3 scripts/serial_sae/download_tcga_ov_api.py first")
        return
    
    expression_matrix = process_all_rnaseq_files(DOWNLOAD_DIR)
    
    if not expression_matrix.empty:
        print()
        print("="*80)
        print("✅ Processing complete")
        print("="*80)
        print()
        print("Next steps:")
        print("1. Compute SAE pathway scores: scripts/serial_sae/compute_tcga_ov_sae.py")
        print("2. Match with clinical outcomes")
        print("3. Validate pathway signatures")

if __name__ == "__main__":
    main()
