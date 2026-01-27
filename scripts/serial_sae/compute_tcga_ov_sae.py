#!/usr/bin/env python3
"""
Compute SAE Pathway Scores for TCGA-OV
Process expression matrix and compute 7D mechanism vectors
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List
import json
import sys

# Add backend to path for SAE service
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

from api.services.sae_feature_service import SAEFeatureService
from api.services.pathway_to_mechanism_vector import convert_pathway_scores_to_mechanism_vector

DATA_DIR = Path("data/serial_sae/tcga_ov")
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def compute_pathway_scores_from_expression(
    expression_df: pd.DataFrame,
    pathway_genes: Dict[str, List[str]]
) -> pd.DataFrame:
    """Compute pathway scores from expression matrix."""
    pathway_scores = {}
    
    for pathway_name, genes in pathway_genes.items():
        # Get expression for pathway genes (case-insensitive match)
        pathway_expr = expression_df[expression_df.index.str.upper().isin([g.upper() for g in genes])]
        
        if not pathway_expr.empty:
            # Mean expression per sample (across pathway genes)
            pathway_scores[pathway_name] = pathway_expr.mean(axis=0)
        else:
            # No genes found - return zeros
            pathway_scores[pathway_name] = pd.Series(0.0, index=expression_df.columns)
    
    return pd.DataFrame(pathway_scores)

def normalize_expression(expression_df: pd.DataFrame, method: str = "log2_tpm") -> pd.DataFrame:
    """Normalize expression data."""
    if method == "log2_tpm":
        # Log2 transform (add pseudocount to avoid log(0))
        return np.log2(expression_df + 1)
    elif method == "zscore":
        # Z-score normalization per gene
        return expression_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    else:
        return expression_df

def main():
    print("="*80)
    print("Compute SAE Pathway Scores for TCGA-OV")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load expression matrix
    expr_file = PROCESSED_DIR / "tcga_ov_expression_matrix.csv"
    if not expr_file.exists():
        print(f"❌ Error: Expression matrix not found: {expr_file}")
        print("Run: python3 scripts/serial_sae/process_tcga_ov_rnaseq.py first")
        return
    
    print(f"📥 Loading expression matrix...")
    expression_df = pd.read_csv(expr_file, index_col=0)
    print(f"✅ Loaded: {len(expression_df)} genes × {len(expression_df.columns)} samples")
    
    # Normalize expression
    print(f"📊 Normalizing expression (log2)...")
    expression_normalized = normalize_expression(expression_df, method="log2_tpm")
    
    # Pathway gene definitions (from codebase)
    pathway_genes = {
        "ddr": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK2", "RAD51", "PALB2", "RAD51C", "RAD51D", "BARD1", "NBN", "TP53", "PARP1", "PARP2", "MBD4", "MLH1", "MSH2", "MSH6", "PMS2"],
        "ras_mapk": ["NF1", "KRAS", "BRAF", "NRAS", "MEK1", "MEK2", "ERK1", "ERK2", "MAPK1", "MAPK3", "MAP2K1", "EGFR"],
        "pi3k": ["PIK3CA", "PIK3CB", "PIK3CD", "PTEN", "AKT1", "AKT2", "AKT3", "MTOR"],
        "vegf": ["VEGFA", "VEGFR1", "VEGFR2", "KDR", "FLT1"],
        "her2": ["ERBB2", "HER2"],
        "io": [],  # IO pathway needs TMB/MSI - set separately
        "efflux": ["ABCB1", "MDR1", "ABCC1", "ABCG2"]
    }
    
    print(f"🧬 Computing pathway scores...")
    pathway_scores_df = compute_pathway_scores_from_expression(expression_normalized, pathway_genes)
    
    print(f"✅ Pathway scores: {len(pathway_scores_df)} pathways × {len(pathway_scores_df.columns)} samples")
    
    # Convert to mechanism vectors
    print(f"🔄 Converting to 7D mechanism vectors...")
    mechanism_vectors = []
    sample_ids = pathway_scores_df.columns.tolist()
    
    for sample_id in sample_ids:
        pathway_scores = pathway_scores_df[sample_id].to_dict()
        mechanism_vector = convert_pathway_scores_to_mechanism_vector(pathway_scores)
        mechanism_vectors.append({
            "sample_id": sample_id,
            "pathway_scores": pathway_scores,
            "mechanism_vector": mechanism_vector
        })
    
    # Save results
    results_df = pd.DataFrame(mechanism_vectors)
    results_file = RESULTS_DIR / "tcga_ov_sae_scores.json"
    with open(results_file, "w") as f:
        json.dump(mechanism_vectors, f, indent=2)
    print(f"💾 Saved: {results_file}")
    
    # Save pathway scores matrix
    pathway_file = RESULTS_DIR / "tcga_ov_pathway_scores.csv"
    pathway_scores_df.to_csv(pathway_file)
    print(f"💾 Saved: {pathway_file}")
    
    print()
    print("="*80)
    print("✅ SAE computation complete")
    print("="*80)
    print(f"\nNext: Correlate with clinical outcomes")

if __name__ == "__main__":
    main()
