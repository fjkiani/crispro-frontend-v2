import pandas as pd
import numpy as np
from scipy import stats
import argparse
from pathlib import Path

# Embedded Hallmark Gene Sets (Minimal viable signatures)
# Source: MSigDB v2023.1.Hs
HALLMARK_DNA_REPAIR = [
    "DDB2", "PCNA", "CHEK1", "CHEK2", "BRCA1", "BRCA2", "FANCA", "FANCC", "RAD51", 
    "XRCC1", "XRCC2", "XRCC3", "ERCC1", "ERCC2", "ATM", "ATR", "MDC1", "TP53", "MLH1", "MSH2"
]
# Broaden validation with E2F (Cell Cycle)
HALLMARK_E2F_TARGETS = [
    "E2F1", "E2F2", "E2F3", "CCNE1", "CCNE2", "CCNA2", "CCNB1", "CCNB2", "CDK1", "CDK2",
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "PCNA", "TYMS", "DHFR", "TK1"
]

def load_expression_data(file_path):
    print(f"Loading expression data from {file_path}...")
    # Support CSV, TSV, GZ
    if str(file_path).endswith('.tsv') or str(file_path).endswith('.tsv.gz'):
        sep = '\t'
    else:
        sep = ','
        
    # Auto-detect orientation: TCGA usually Columns=Patients, Rows=Genes
    df = pd.read_csv(file_path, sep=sep, index_col=0)
    
    # Check if index are genes (strings)
    if df.index.dtype == 'O':
        print(f"Context: {df.shape[0]} genes x {df.shape[1]} samples")
        return df
    else:
        # Maybe transposed?
        print("⚠️ Warning: Index doesn't look like gene symbols. Attempting transpose...")
        return df.T

def calculate_z_score_signature(expression_df, gene_set, name):
    """
    Computes Z-score signature:
    1. Subset genes in gene_set
    2. Z-score each gene across all samples
    3. Average Z-scores for the signature
    """
    valid_genes = [g for g in gene_set if g in expression_df.index]
    print(f"Computing {name}: Found {len(valid_genes)}/{len(gene_set)} genes")
    
    if not valid_genes:
        print(f"❌ No genes found for {name}!")
        return pd.Series(0, index=expression_df.columns)
        
    subset = expression_df.loc[valid_genes]
    
    # Z-score normalization (genes across samples)
    # axis=1 (columns) for stats.zscore usually z-scores the column, so we transpose or be careful
    # We want z-score of a GENE across SAMPLES.
    # If rows=genes, cols=samples:
    z_matrix = subset.apply(stats.zscore, axis=1) # Z-score along rows (samples)
    
    # Average Z-scores
    scores = z_matrix.mean(axis=0)
    return scores

def main():
    parser = argparse.ArgumentParser(description="Score TCGA-OV with Hallmark Pathways")
    parser.add_argument("--input", required=True, help="Path to expression matrix")
    parser.add_argument("--output", default="tcga_ov_pathway_scores.csv", help="Output score file")
    
    args = parser.parse_args()
    
    # 1. Load
    expr = load_expression_data(args.input)
    
    # 2. Score
    ddr_scores = calculate_z_score_signature(expr, HALLMARK_DNA_REPAIR, "Hallmark_DDR")
    e2f_scores = calculate_z_score_signature(expr, HALLMARK_E2F_TARGETS, "Hallmark_E2F")
    
    # 3. Combine
    results = pd.DataFrame({
        "DDR_Score": ddr_scores,
        "E2F_Score": e2f_scores
    })
    
    # 4. Save
    results.to_csv(args.output)
    print(f"✅ Scores saved to {args.output}")

if __name__ == "__main__":
    main()
