#!/usr/bin/env python3
"""
Process GSE241908: Platinum-resistant ovarian cancer + bevacizumab
7 patients, 14 RNA-seq samples (7 paired pre/post bevacizumab)

Process TPM expression matrix and compute SAE pathway scores
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import json
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import SAE computation utilities
try:
    from api.services.sae_feature_service import get_sae_feature_service
    SAE_SERVICE_AVAILABLE = True
except ImportError:
    SAE_SERVICE_AVAILABLE = False
    print("⚠️  SAE service not available, using simplified pathway computation")

# Data directories
DATA_DIR = Path("data/serial_sae/gse241908")
RESULTS_DIR = DATA_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Pathway gene lists (from SAE service or manual definition)
PATHWAY_GENES = {
    "DDR": ["BRCA1", "BRCA2", "PALB2", "RAD51C", "RAD51D", "BRIP1", "BARD1", "ATM", "ATR", "CHEK1", "CHEK2"],
    "MAPK": ["KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3"],
    "PI3K": ["PIK3CA", "PIK3CB", "PIK3R1", "PTEN", "AKT1", "AKT2", "MTOR"],
    "VEGF": ["VEGFA", "VEGFB", "VEGFC", "VEGFR1", "VEGFR2", "VEGFR3", "PDGFRA", "PDGFRB"],
    "HER2": ["ERBB2", "ERBB3", "ERBB4", "EGFR"],
    "Efflux": ["ABCB1", "ABCC1", "ABCC2", "ABCG2"]
}

def load_expression_matrix(file_path: Path) -> pd.DataFrame:
    """Load TPM expression matrix from gzipped CSV (space-separated with quotes)."""
    print(f"Loading expression matrix: {file_path}")
    
    import csv
    import shlex
    
    with gzip.open(file_path, 'rt') as f:
        # Read first line to get sample names
        first_line = f.readline().strip()
        # Use shlex to properly parse quoted strings
        sample_names = shlex.split(first_line)
        print(f"  Found {len(sample_names)} samples: {sample_names[:5]}...")
        
        # Read data rows
        data_rows = []
        gene_ids = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Parse using shlex to handle quotes
            parts = shlex.split(line)
            if len(parts) < 2:
                continue
            
            gene_id = parts[0]
            values = [float(x) for x in parts[1:]]
            
            # Check if we have the right number of values
            if len(values) == len(sample_names):
                gene_ids.append(gene_id)
                data_rows.append(values)
            else:
                print(f"  ⚠️  Skipping row with {len(values)} values (expected {len(sample_names)}): {gene_id}")
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, index=gene_ids, columns=sample_names)
    
    print(f"  Loaded {df.shape[0]} genes × {df.shape[1]} samples")
    print(f"  Gene IDs (first 3): {list(df.index)[:3]}")
    return df

def identify_timepoints(sample_names: List[str]) -> Dict[str, str]:
    """
    Identify pre/post bevacizumab timepoints.
    Load from sample mapping JSON if available, otherwise infer from names.
    """
    # Try to load sample mapping
    sample_mapping_path = DATA_DIR / "sample_mapping.json"
    if sample_mapping_path.exists():
        with open(sample_mapping_path) as f:
            sample_mapping = json.load(f)
        
        timepoint_map = {}
        for sample in sample_names:
            if sample in sample_mapping:
                timepoint_map[sample] = sample_mapping[sample].get("timepoint", "unknown")
            else:
                timepoint_map[sample] = "unknown"
        return timepoint_map
    
    # Fallback: infer from sample names
    timepoint_map = {}
    for sample in sample_names:
        sample_lower = sample.lower()
        if any(keyword in sample_lower for keyword in ["baseline", "pre", "before", "_t0", "_t1", "pre-bev"]):
            timepoint_map[sample] = "pre"
        elif any(keyword in sample_lower for keyword in ["post", "after", "_t2", "_t3", "post-bev", "bevacizumab"]):
            timepoint_map[sample] = "post"
        else:
            timepoint_map[sample] = "unknown"
    
    return timepoint_map

def load_ensembl_mapping() -> Dict[str, str]:
    """Load Ensembl ID → Gene Symbol mapping."""
    mapping_path = DATA_DIR / "ensembl_to_symbol.json"
    if mapping_path.exists():
        with open(mapping_path) as f:
            return json.load(f)
    return {}

def map_ensembl_to_symbols(ensembl_ids: pd.Index, mapping: Dict[str, str]) -> pd.Index:
    """
    Map Ensembl IDs to gene symbols using provided mapping.
    Returns Index with symbols (or original IDs if not found).
    """
    # Create reverse mapping: symbol → list of Ensembl IDs
    symbol_to_ensembl = {}
    for ens_id, symbol in mapping.items():
        if symbol not in symbol_to_ensembl:
            symbol_to_ensembl[symbol] = []
        symbol_to_ensembl[symbol].append(ens_id)
    
    # Map each Ensembl ID to symbol
    mapped_symbols = []
    for ens_id in ensembl_ids:
        # Remove version number for matching (ENSG00000000003.14 -> ENSG00000000003)
        base_id = ens_id.split('.')[0] if '.' in ens_id else ens_id
        
        if base_id in mapping:
            mapped_symbols.append(mapping[base_id])
        else:
            # Keep original if not found
            mapped_symbols.append(ens_id)
    
    return pd.Index(mapped_symbols)

def compute_pathway_score(expression_df: pd.DataFrame, pathway_name: str, pathway_genes: List[str], ensembl_mapping: Dict[str, str] = None) -> pd.Series:
    """
    Compute pathway score as mean expression of pathway genes.
    
    Args:
        expression_df: Genes × Samples expression matrix (Ensembl IDs as index)
        pathway_name: Name of pathway
        pathway_genes: List of gene symbols in pathway
        ensembl_to_symbol: Optional mapping from Ensembl ID to gene symbol
    
    Returns:
        Series with pathway scores per sample
    """
    # If we have Ensembl IDs, we need to map to symbols
    # For now, use a simplified approach: search for genes by name in a mapping
    # Or use all genes and filter by pathway membership later
    
    # Check if index is Ensembl IDs
    is_ensembl = expression_df.index[0].startswith("ENSG") if len(expression_df.index) > 0 else False
    
    if is_ensembl:
        if not ensembl_mapping:
            print(f"  ⚠️  Ensembl IDs detected but no mapping provided for {pathway_name}")
            return pd.Series(0.0, index=expression_df.columns)
        
        # Create reverse mapping: symbol → Ensembl ID
        symbol_to_ensembl = {v: k for k, v in ensembl_mapping.items()}
        
        # Find Ensembl IDs for pathway genes
        # Expression matrix has Ensembl IDs with version numbers (ENSG00000000003.14)
        # Mapping has base IDs (ENSG00000000003) - need to match by base ID
        pathway_ensembl_ids = []
        expression_base_ids = {idx.split('.')[0]: idx for idx in expression_df.index}
        
        for gene_symbol in pathway_genes:
            if gene_symbol in symbol_to_ensembl:
                mapping_ensembl_id = symbol_to_ensembl[gene_symbol]
                base_id = mapping_ensembl_id.split('.')[0] if '.' in mapping_ensembl_id else mapping_ensembl_id
                
                # Match by base ID
                if base_id in expression_base_ids:
                    full_id = expression_base_ids[base_id]
                    pathway_ensembl_ids.append(full_id)
        
        if len(pathway_ensembl_ids) == 0:
            print(f"  ⚠️  No Ensembl IDs found for {pathway_name} pathway genes")
            return pd.Series(0.0, index=expression_df.columns)
        
        available_genes = pathway_ensembl_ids
    else:
        # Find genes that exist in expression matrix (already symbols)
        available_genes = [g for g in pathway_genes if g in expression_df.index]
    
    if len(available_genes) == 0:
        print(f"  ⚠️  No genes found for {pathway_name}")
        return pd.Series(0.0, index=expression_df.columns)
    
    # Compute mean expression (log2(TPM+1) normalized)
    pathway_expr = expression_df.loc[available_genes]
    
    # Log2 transform
    pathway_expr_log = np.log2(pathway_expr + 1)
    
    # Mean across genes
    pathway_score = pathway_expr_log.mean(axis=0)
    
    print(f"  {pathway_name}: {len(available_genes)}/{len(pathway_genes)} genes, mean score range: [{pathway_score.min():.3f}, {pathway_score.max():.3f}]")
    
    return pathway_score

def compute_pathway_kinetics(pathway_scores_df: pd.DataFrame, timepoint_map: Dict[str, str]) -> pd.DataFrame:
    """
    Compute pathway kinetics (Δ = post - pre) for each patient.
    
    Args:
        pathway_scores_df: Pathways × Samples dataframe
        timepoint_map: Sample → timepoint mapping
    
    Returns:
        DataFrame with patient_id, pathway, Δ value
    """
    kinetics_data = []
    
    # Group samples by patient (assume patient ID in sample name)
    # For GSE241908, we need to infer patient IDs from sample names
    # Common patterns: "Patient1_pre", "P1_baseline", etc.
    
    # Try to extract patient IDs
    patient_samples = {}
    for sample in pathway_scores_df.columns:
        # Try common patterns
        parts = sample.split('_')
        if len(parts) >= 2:
            patient_id = parts[0]  # Assume first part is patient ID
        else:
            patient_id = sample.split('-')[0] if '-' in sample else sample
        
        if patient_id not in patient_samples:
            patient_samples[patient_id] = {"pre": [], "post": []}
        
        timepoint = timepoint_map.get(sample, "unknown")
        if timepoint in ["pre", "post"]:
            patient_samples[patient_id][timepoint].append(sample)
    
    # Compute Δ for each patient
    for patient_id, samples in patient_samples.items():
        if len(samples["pre"]) == 0 or len(samples["post"]) == 0:
            continue
        
        # Use first pre and first post sample (or average if multiple)
        pre_samples = samples["pre"]
        post_samples = samples["post"]
        
        for pathway in pathway_scores_df.index:
            pre_score = pathway_scores_df.loc[pathway, pre_samples].mean()
            post_score = pathway_scores_df.loc[pathway, post_samples].mean()
            
            delta = post_score - pre_score
            
            kinetics_data.append({
                "patient_id": patient_id,
                "pathway": pathway,
                "pre_score": pre_score,
                "post_score": post_score,
                "delta": delta
            })
    
    return pd.DataFrame(kinetics_data)

def main():
    """Process GSE241908 and compute pathway kinetics."""
    print("="*80)
    print("GSE241908 Processing: Pathway Kinetics Analysis")
    print("="*80)
    
    # Load expression matrix
    tpm_file = DATA_DIR / "GSE241908_data1_TPM.csv.gz"
    if not tpm_file.exists():
        print(f"❌ Expression matrix not found: {tpm_file}")
        print("   Run download_gse241908.py first")
        return
    
    expression_df = load_expression_matrix(tpm_file)
    
    # Identify timepoints
    print("\nIdentifying timepoints...")
    timepoint_map = identify_timepoints(list(expression_df.columns))
    
    pre_count = sum(1 for v in timepoint_map.values() if v == "pre")
    post_count = sum(1 for v in timepoint_map.values() if v == "post")
    unknown_count = sum(1 for v in timepoint_map.values() if v == "unknown")
    
    print(f"  Pre-treatment: {pre_count} samples")
    print(f"  Post-treatment: {post_count} samples")
    print(f"  Unknown: {unknown_count} samples")
    
    if unknown_count > 0:
        print(f"\n  ⚠️  Sample timepoint mapping:")
        for sample, tp in timepoint_map.items():
            print(f"    {sample}: {tp}")
    
    # Load Ensembl mapping
    print("\nLoading Ensembl ID → Gene Symbol mapping...")
    ensembl_mapping = load_ensembl_mapping()
    if ensembl_mapping:
        print(f"  ✅ Loaded {len(ensembl_mapping)} gene mappings")
    else:
        print("  ⚠️  No Ensembl mapping found - pathway scores may be incomplete")
    
    # Compute pathway scores
    print("\nComputing pathway scores...")
    pathway_scores = {}
    
    for pathway_name, pathway_genes in PATHWAY_GENES.items():
        pathway_scores[pathway_name] = compute_pathway_score(
            expression_df, pathway_name, pathway_genes, ensembl_mapping
        )
    
    pathway_scores_df = pd.DataFrame(pathway_scores)
    pathway_scores_df = pathway_scores_df.T  # Pathways × Samples
    
    # Save pathway scores
    scores_path = RESULTS_DIR / "pathway_scores.csv"
    pathway_scores_df.to_csv(scores_path)
    print(f"\n💾 Pathway scores saved: {scores_path}")
    
    # Compute pathway kinetics
    print("\nComputing pathway kinetics (Δ = post - pre)...")
    kinetics_df = compute_pathway_kinetics(pathway_scores_df, timepoint_map)
    
    if len(kinetics_df) == 0:
        print("  ⚠️  No paired samples found - check timepoint mapping")
        print("  Saving pathway scores only")
        return
    
    # Save kinetics
    kinetics_path = RESULTS_DIR / "pathway_kinetics.csv"
    kinetics_df.to_csv(kinetics_path, index=False)
    print(f"💾 Pathway kinetics saved: {kinetics_path}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("PATHWAY KINETICS SUMMARY")
    print("="*80)
    
    for pathway in PATHWAY_GENES.keys():
        pathway_deltas = kinetics_df[kinetics_df['pathway'] == pathway]['delta']
        if len(pathway_deltas) > 0:
            mean_delta = pathway_deltas.mean()
            std_delta = pathway_deltas.std()
            print(f"\n{pathway}:")
            print(f"  Mean Δ: {mean_delta:.4f} ± {std_delta:.4f}")
            print(f"  Range: [{pathway_deltas.min():.4f}, {pathway_deltas.max():.4f}]")
            print(f"  Patients: {len(pathway_deltas)}")
    
    # Key finding: VEGF pathway (bevacizumab target)
    vegf_deltas = kinetics_df[kinetics_df['pathway'] == 'VEGF']['delta']
    if len(vegf_deltas) > 0:
        print("\n" + "="*80)
        print("KEY FINDING: VEGF Pathway (Bevacizumab Target)")
        print("="*80)
        print(f"  Mean ΔVEGF: {vegf_deltas.mean():.4f}")
        if vegf_deltas.mean() < 0:
            print("  ✅ VEGF DECREASED post-bevacizumab (expected - anti-angiogenic effect)")
        else:
            print("  ⚠️  VEGF INCREASED post-bevacizumab (unexpected - possible resistance)")
    
    print("\n✅ Processing complete!")

if __name__ == "__main__":
    main()
