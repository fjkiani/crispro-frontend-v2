#!/usr/bin/env python3
"""
Process Manually Downloaded cBioPortal Data
Parse TSV files and extract paired sample data
"""

import pandas as pd
import json
from pathlib import Path
from typing import Dict, List

DATA_DIR = Path("data/serial_sae/cbioportal_paired")
OUTPUT_DIR = DATA_DIR

def load_paired_patients() -> Dict:
    """Load paired patients JSON."""
    with open(DATA_DIR / "paired_patients.json") as f:
        data = json.load(f)
    return data.get("paired_patients", {})

def process_mutations_tsv(tsv_path: Path, paired_patients: Dict) -> pd.DataFrame:
    """Process mutations TSV file and extract paired sample data."""
    print(f"📥 Loading mutations from {tsv_path}...")
    
    # Load TSV (could be MAF format or cBioPortal format)
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    print(f"  Loaded {len(df)} mutation rows")
    print(f"  Columns: {list(df.columns)[:10]}...")
    
    # Identify sample ID column (could be "Tumor_Sample_Barcode", "Sample_ID", etc.)
    sample_col = None
    for col in ["Tumor_Sample_Barcode", "Sample_ID", "sampleId", "Sample ID"]:
        if col in df.columns:
            sample_col = col
            break
    
    if not sample_col:
        print("⚠️  Could not identify sample ID column")
        return pd.DataFrame()
    
    # Extract paired sample IDs
    paired_sample_ids = set()
    for patient_data in paired_patients.values():
        for sample in patient_data.get("primary_samples", []):
            paired_sample_ids.add(sample["sample_id"])
        for sample in patient_data.get("recurrent_samples", []):
            paired_sample_ids.add(sample["sample_id"])
    
    # Filter to paired samples only
    df_filtered = df[df[sample_col].isin(paired_sample_ids)].copy()
    print(f"  Filtered to {len(df_filtered)} mutations for paired samples")
    
    return df_filtered

def process_expression_tsv(tsv_path: Path, paired_patients: Dict) -> pd.DataFrame:
    """Process expression TSV file and extract paired sample data."""
    print(f"📥 Loading expression from {tsv_path}...")
    
    # Load TSV (genes as rows, samples as columns)
    df = pd.read_csv(tsv_path, sep="\t", index_col=0, low_memory=False)
    print(f"  Loaded {len(df)} genes, {len(df.columns)} samples")
    
    # Extract paired sample IDs
    paired_sample_ids = set()
    for patient_data in paired_patients.values():
        for sample in patient_data.get("primary_samples", []):
            paired_sample_ids.add(sample["sample_id"])
        for sample in patient_data.get("recurrent_samples", []):
            paired_sample_ids.add(sample["sample_id"])
    
    # Filter columns to paired samples only
    available_samples = [col for col in df.columns if col in paired_sample_ids]
    df_filtered = df[available_samples]
    print(f"  Filtered to {len(df_filtered.columns)} paired samples")
    
    return df_filtered

def main():
    print("="*80)
    print("Process Manually Downloaded cBioPortal Data")
    print("="*80)
    print()
    
    # Load paired patients
    paired_patients = load_paired_patients()
    print(f"📊 Loaded {len(paired_patients)} paired patients")
    
    # Process mutations
    mut_file = DATA_DIR / "msk_spectrum_tme_2022_mutations.tsv"
    if mut_file.exists():
        mutations_df = process_mutations_tsv(mut_file, paired_patients)
        if not mutations_df.empty:
            output_mut = OUTPUT_DIR / "processed_mutations.csv"
            mutations_df.to_csv(output_mut, index=False)
            print(f"💾 Saved processed mutations: {output_mut}")
    else:
        print(f"⚠️  Mutations file not found: {mut_file}")
        print("   Please download from cBioPortal and save to this location")
    
    # Process expression
    expr_file = DATA_DIR / "msk_spectrum_tme_2022_expression.tsv"
    if expr_file.exists():
        expression_df = process_expression_tsv(expr_file, paired_patients)
        if not expression_df.empty:
            output_expr = OUTPUT_DIR / "processed_expression.csv"
            expression_df.to_csv(output_expr)
            print(f"💾 Saved processed expression: {output_expr}")
    else:
        print(f"⚠️  Expression file not found: {expr_file}")
        print("   Please download from cBioPortal and save to this location")
    
    print("\n" + "="*80)
    print("✅ Processing complete")
    print("="*80)
    print("\nNext step: Run compute_serial_sae.py to compute pathway scores")

if __name__ == "__main__":
    main()
