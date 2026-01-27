#!/usr/bin/env python3
"""
Pathway Kinetics Validation (GSE165897)
Process scRNA-seq data → pseudo-bulk → pathway scores → kinetics (Δ values)

Mission: Validate SAE pathway kinetics detect mechanism-specific chemotherapy resistance
Dataset: GSE165897 (n=11 HGSOC patients, scRNA-seq, pre-treatment + post-NACT)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import sys

# Set matplotlib to non-interactive backend BEFORE any matplotlib imports
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

# Note: We use pandas directly for TSV files, no scanpy needed

# Data is downloaded in scripts/data_acquisition/sae/
DATA_DIR = Path("scripts/data_acquisition/sae")
RESULTS_DIR = Path("data/serial_sae/gse165897/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Pathway gene lists (from audit)
DDR_GENES = [
    "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", "PARP1"
]

MAPK_GENES = [
    "KRAS", "NRAS", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2"
]

PI3K_GENES = [
    "PIK3CA", "AKT1", "AKT2", "PTEN", "mTOR"
]

VEGF_GENES = [
    "VEGFA", "VEGFR1", "VEGFR2", "HIF1A"
]

PATHWAY_GENES = {
    "ddr": DDR_GENES,
    "mapk": MAPK_GENES,
    "pi3k": PI3K_GENES,
    "vegf": VEGF_GENES
}

def load_gse165897_data(data_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load GSE165897 scRNA-seq data from TSV files.
    
    Args:
        data_dir: Directory containing GSE165897 files
    
    Returns:
        Tuple of (expression_df, cell_info_df)
        - expression_df: genes × cells (UMI counts)
        - cell_info_df: cells × metadata (patient_id, treatment_phase, etc.)
    """
    print("📥 Loading GSE165897 data...")
    
    # Load cell metadata
    cell_info_path = data_dir / "GSE165897_cellInfo_HGSOC.tsv.gz"
    cell_info_df = pd.read_csv(cell_info_path, sep="\t", compression="gzip")
    print(f"  Loaded {len(cell_info_df)} cells with metadata")
    
    # Load expression matrix (genes × cells)
    expression_path = data_dir / "GSE165897_UMIcounts_HGSOC.tsv.gz"
    expression_df = pd.read_csv(expression_path, sep="\t", compression="gzip", index_col=0)
    print(f"  Loaded {len(expression_df)} genes × {len(expression_df.columns)} cells")
    
    # Transpose to cells × genes for easier aggregation
    expression_df = expression_df.T
    
    # Ensure cell IDs match between expression and metadata
    # Cell IDs in expression are column names, in metadata they're in 'cell' column
    expression_df.index.name = "cell"
    expression_df = expression_df.reset_index()
    
    # Merge with cell info to get patient_id and treatment_phase
    merged_df = expression_df.merge(cell_info_df[["cell", "patient_id", "treatment_phase", "sample"]], 
                                    on="cell", how="inner")
    
    print(f"  Matched {len(merged_df)} cells after merging")
    
    return merged_df, cell_info_df

def aggregate_to_pseudobulk(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate scRNA-seq data to pseudo-bulk per patient-timepoint.
    
    Method:
    1. Group by patient_id + treatment_phase
    2. Sum UMI counts per gene per group
    3. Normalize to CPM (counts per million)
    4. Log2 transform: log2(CPM + 1)
    
    Args:
        merged_df: DataFrame with cells × genes + metadata columns
    
    Returns:
        DataFrame: patient_timepoint × genes (log2(CPM + 1))
    """
    print("📊 Aggregating to pseudo-bulk by patient_id + treatment_phase...")
    
    # Create group identifier
    merged_df["patient_timepoint"] = merged_df["patient_id"].astype(str) + "_" + merged_df["treatment_phase"].astype(str)
    
    # Get gene columns (exclude metadata columns)
    metadata_cols = ["cell", "patient_id", "treatment_phase", "sample", "patient_timepoint"]
    gene_cols = [col for col in merged_df.columns if col not in metadata_cols]
    
    # Group by patient_timepoint and sum counts
    pseudobulk_list = []
    for group_name, group_df in merged_df.groupby("patient_timepoint"):
        # Sum UMI counts per gene
        counts_sum = group_df[gene_cols].sum(axis=0)
        
        # Normalize to CPM
        total_counts = counts_sum.sum()
        if total_counts > 0:
            cpm = (counts_sum / total_counts) * 1e6
        else:
            cpm = pd.Series(0.0, index=gene_cols)
        
        # Log2 transform
        log2_cpm = np.log2(cpm + 1)
        
        pseudobulk_list.append({
            "patient_timepoint": group_name,
            "patient_id": group_df["patient_id"].iloc[0],
            "treatment_phase": group_df["treatment_phase"].iloc[0],
            **log2_cpm.to_dict()
        })
    
    pseudobulk_df = pd.DataFrame(pseudobulk_list)
    pseudobulk_df.set_index("patient_timepoint", inplace=True)
    
    print(f"✅ Aggregated to {len(pseudobulk_df)} patient-timepoints")
    return pseudobulk_df

def compute_pathway_score_from_expression(
    expression_vector: pd.Series,
    pathway_genes: List[str]
) -> float:
    """
    Compute pathway burden score from gene expression.
    
    Formula: mean(log2(expression + 1)) for pathway genes
    Normalize to 0-1 scale
    
    Args:
        expression_vector: Series of gene → expression (log2(CPM + 1))
        pathway_genes: List of gene symbols
    
    Returns:
        Pathway score (0.0-1.0)
    """
    # Get expression values for pathway genes (case-insensitive)
    pathway_expressions = []
    for gene in pathway_genes:
        # Try exact match first
        if gene in expression_vector.index:
            pathway_expressions.append(expression_vector[gene])
        else:
            # Try case-insensitive match
            matches = [g for g in expression_vector.index if g.upper() == gene.upper()]
            if matches:
                pathway_expressions.append(expression_vector[matches[0]])
    
    if not pathway_expressions:
        return 0.0
    
    # Mean of log2(expression + 1) - already log2 transformed, so just mean
    mean_expression = np.mean(pathway_expressions)
    
    # Normalize to 0-1 scale (empirical range: 0-15 for log2(CPM+1))
    # Using 0-15 range: 0 → 0.0, 7.5 → 0.5, 15 → 1.0
    normalized = min(1.0, max(0.0, mean_expression / 15.0))
    
    return normalized

def compute_pathway_scores(pseudobulk_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute pathway scores for all groups.
    
    Args:
        pseudobulk_df: DataFrame of groups × genes (log2(CPM + 1))
    
    Returns:
        DataFrame of groups × pathways (scores 0.0-1.0)
    """
    print("🧬 Computing pathway scores...")
    
    pathway_scores = []
    for group in pseudobulk_df.index:
        expression_vector = pseudobulk_df.loc[group]
        
        scores = {
            "group": group,
            "ddr": compute_pathway_score_from_expression(expression_vector, PATHWAY_GENES["ddr"]),
            "mapk": compute_pathway_score_from_expression(expression_vector, PATHWAY_GENES["mapk"]),
            "pi3k": compute_pathway_score_from_expression(expression_vector, PATHWAY_GENES["pi3k"]),
            "vegf": compute_pathway_score_from_expression(expression_vector, PATHWAY_GENES["vegf"])
        }
        pathway_scores.append(scores)
    
    pathway_scores_df = pd.DataFrame(pathway_scores)
    pathway_scores_df.set_index("group", inplace=True)
    
    print(f"✅ Computed pathway scores for {len(pathway_scores_df)} groups")
    return pathway_scores_df

def compute_pathway_kinetics(pathway_scores_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute pathway kinetics (Δ values: post - pre).
    
    Args:
        pathway_scores_df: DataFrame with index = patient_timepoint, columns include patient_id, treatment_phase, pathway scores
    
    Returns:
        DataFrame of patients × pathways (Δ values)
    """
    print("📈 Computing pathway kinetics (Δ values)...")
    
    kinetics = []
    
    # Group by patient_id
    for patient_id in pathway_scores_df["patient_id"].unique():
        patient_data = pathway_scores_df[pathway_scores_df["patient_id"] == patient_id]
        
        # Find pre and post timepoints
        pre_data = patient_data[patient_data["treatment_phase"].str.contains("naive|pre|baseline", case=False, na=False)]
        post_data = patient_data[patient_data["treatment_phase"].str.contains("post|nact", case=False, na=False)]
        
        if len(pre_data) == 0 or len(post_data) == 0:
            continue
        
        # Use first match if multiple (should be unique per patient)
        pre_scores = pre_data.iloc[0]
        post_scores = post_data.iloc[0]
        
        # Compute delta (post - pre)
        kinetics.append({
            "patient_id": patient_id,
            "ddr_delta": post_scores["ddr"] - pre_scores["ddr"],
            "mapk_delta": post_scores["mapk"] - pre_scores["mapk"],
            "pi3k_delta": post_scores["pi3k"] - pre_scores["pi3k"],
            "vegf_delta": post_scores["vegf"] - pre_scores["vegf"],
            "ddr_pre": pre_scores["ddr"],
            "ddr_post": post_scores["ddr"],
            "mapk_pre": pre_scores["mapk"],
            "mapk_post": post_scores["mapk"],
            "pi3k_pre": pre_scores["pi3k"],
            "pi3k_post": post_scores["pi3k"],
            "vegf_pre": pre_scores["vegf"],
            "vegf_post": post_scores["vegf"]
        })
    
    kinetics_df = pd.DataFrame(kinetics)
    
    print(f"✅ Computed kinetics for {len(kinetics_df)} patients")
    return kinetics_df

def classify_kinetic_patterns(kinetics_df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify kinetic patterns (DDR restoration, MAPK activation, etc.).
    
    Args:
        kinetics_df: DataFrame with Δ values
    
    Returns:
        DataFrame with pattern classifications
    """
    print("🔍 Classifying kinetic patterns...")
    
    patterns = []
    for _, row in kinetics_df.iterrows():
        pattern = {
            "patient_id": row["patient_id"],
            "ddr_restoration": row["ddr_delta"] > 0.1,  # DDR increases (restoration)
            "mapk_activation": row["mapk_delta"] > 0.1,  # MAPK increases (activation)
            "pi3k_activation": row["pi3k_delta"] > 0.1,  # PI3K increases (activation)
            "vegf_activation": row["vegf_delta"] > 0.1,  # VEGF increases (activation)
            "ddr_suppression": row["ddr_delta"] < -0.1,  # DDR decreases (suppression)
            "mapk_suppression": row["mapk_delta"] < -0.1,  # MAPK decreases (suppression)
        }
        patterns.append(pattern)
    
    patterns_df = pd.DataFrame(patterns)
    
    print(f"✅ Classified patterns for {len(patterns_df)} patients")
    return patterns_df

def main():
    print("="*80)
    print("Pathway Kinetics Validation (GSE165897)")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Step 1: Load GSE165897 data from TSV files
    if not (DATA_DIR / "GSE165897_UMIcounts_HGSOC.tsv.gz").exists():
        print(f"❌ Expression data not found: {DATA_DIR / 'GSE165897_UMIcounts_HGSOC.tsv.gz'}")
        return
    
    if not (DATA_DIR / "GSE165897_cellInfo_HGSOC.tsv.gz").exists():
        print(f"❌ Cell metadata not found: {DATA_DIR / 'GSE165897_cellInfo_HGSOC.tsv.gz'}")
        return
    
    print("📂 Step 1: Loading GSE165897 data...")
    merged_df, cell_info_df = load_gse165897_data(DATA_DIR)
    print(f"✅ Loaded {len(merged_df)} cells from {len(cell_info_df)} cell metadata records\n")
    
    # Step 2: Aggregate to pseudo-bulk
    print("📊 Step 2: Aggregating to pseudo-bulk...")
    pseudobulk_df = aggregate_to_pseudobulk(merged_df)
    print(f"✅ Aggregated to {len(pseudobulk_df)} pseudo-bulk groups\n")
    
    # Step 3: Compute pathway scores
    print("🧬 Step 3: Computing pathway scores...")
    pathway_scores_df = compute_pathway_scores(pseudobulk_df)
    print(f"✅ Computed pathway scores\n")
    
    # Add patient_id and treatment_phase back for kinetics computation
    # The index is "group" which contains patient_timepoint
    pathway_scores_df = pathway_scores_df.reset_index()
    pathway_scores_df.rename(columns={"group": "patient_timepoint"}, inplace=True)
    pathway_scores_df["patient_id"] = pathway_scores_df["patient_timepoint"].str.split("_").str[0]
    pathway_scores_df["treatment_phase"] = pathway_scores_df["patient_timepoint"].str.split("_", n=1).str[1]
    pathway_scores_df.set_index("patient_timepoint", inplace=True)
    
    pathway_scores_df.to_csv(RESULTS_DIR / "pathway_scores.csv")
    print(f"💾 Saved pathway scores: {RESULTS_DIR / 'pathway_scores.csv'}\n")
    
    # Step 4: Compute kinetics
    print("📈 Step 4: Computing pathway kinetics (Δ values)...")
    kinetics_df = compute_pathway_kinetics(pathway_scores_df)
    print(f"✅ Computed kinetics for {len(kinetics_df)} patients\n")
    kinetics_df.to_csv(RESULTS_DIR / "pathway_kinetics.csv", index=False)
    print(f"💾 Saved pathway kinetics: {RESULTS_DIR / 'pathway_kinetics.csv'}\n")
    
    # Step 6: Classify patterns
    patterns_df = classify_kinetic_patterns(kinetics_df)
    patterns_df.to_csv(RESULTS_DIR / "kinetic_patterns.csv", index=False)
    print(f"💾 Saved kinetic patterns: {RESULTS_DIR / 'kinetic_patterns.csv'}\n")
    
    # Step 7: Summary statistics
    print("="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"\nPathway Kinetics (Δ values):")
    print(kinetics_df[["ddr_delta", "mapk_delta", "pi3k_delta", "vegf_delta"]].describe())
    
    print(f"\nKinetic Patterns:")
    print(f"  DDR Restoration: {patterns_df['ddr_restoration'].sum()} / {len(patterns_df)}")
    print(f"  MAPK Activation: {patterns_df['mapk_activation'].sum()} / {len(patterns_df)}")
    print(f"  PI3K Activation: {patterns_df['pi3k_activation'].sum()} / {len(patterns_df)}")
    print(f"  VEGF Activation: {patterns_df['vegf_activation'].sum()} / {len(patterns_df)}")
    
    # Step 8: Statistical analysis
    print("\n" + "="*80)
    print("STATISTICAL ANALYSIS")
    print("="*80)
    
    # Correlation between pathways
    pathway_deltas = kinetics_df[["ddr_delta", "mapk_delta", "pi3k_delta", "vegf_delta"]]
    correlation_matrix = pathway_deltas.corr()
    print("\n📊 Pathway Kinetics Correlations:")
    print(correlation_matrix.round(3))
    correlation_matrix.to_csv(RESULTS_DIR / "pathway_correlations.csv")
    print(f"💾 Saved correlation matrix: {RESULTS_DIR / 'pathway_correlations.csv'}\n")
    
    # Effect sizes (Cohen's d for each pathway)
    print("📈 Effect Sizes (Cohen's d):")
    for pathway in ["ddr", "mapk", "pi3k", "vegf"]:
        delta_col = f"{pathway}_delta"
        mean_delta = kinetics_df[delta_col].mean()
        std_delta = kinetics_df[delta_col].std()
        cohens_d = mean_delta / std_delta if std_delta > 0 else 0.0
        print(f"  {pathway.upper()}: d = {cohens_d:.3f} (mean Δ = {mean_delta:.4f}, SD = {std_delta:.4f})")
    
    # Generate heatmap
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Heatmap of pathway kinetics
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Left: Pathway kinetics heatmap (patients × pathways)
        kinetics_matrix = kinetics_df.set_index("patient_id")[["ddr_delta", "mapk_delta", "pi3k_delta", "vegf_delta"]]
        sns.heatmap(kinetics_matrix, annot=True, fmt=".3f", cmap="RdBu_r", center=0, 
                   ax=axes[0], cbar_kws={"label": "Δ (post - pre)"})
        axes[0].set_title("Pathway Kinetics (Δ values)\nby Patient", fontsize=12, fontweight="bold")
        axes[0].set_xlabel("Pathway", fontsize=10)
        axes[0].set_ylabel("Patient ID", fontsize=10)
        
        # Right: Correlation matrix
        sns.heatmap(correlation_matrix, annot=True, fmt=".3f", cmap="coolwarm", center=0,
                   ax=axes[1], square=True, cbar_kws={"label": "Correlation"})
        axes[1].set_title("Pathway Kinetics Correlations", fontsize=12, fontweight="bold")
        axes[1].set_xlabel("Pathway", fontsize=10)
        axes[1].set_ylabel("Pathway", fontsize=10)
        
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / "pathway_kinetics_heatmap.png", dpi=300, bbox_inches="tight")
        print(f"💾 Saved heatmap: {RESULTS_DIR / 'pathway_kinetics_heatmap.png'}\n")
        plt.close()
        
    except ImportError:
        print("⚠️  matplotlib/seaborn not available - skipping heatmap generation\n")
    
    # Summary report
    report_path = RESULTS_DIR / "pathway_kinetics_report.txt"
    with open(report_path, "w") as f:
        f.write("="*80 + "\n")
        f.write("PATHWAY KINETICS VALIDATION REPORT (GSE165897)\n")
        f.write("="*80 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Dataset: GSE165897 (n={len(kinetics_df)} patients, paired pre/post-NACT)\n")
        f.write(f"Pathways Analyzed: DDR, MAPK, PI3K, VEGF\n\n")
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 80 + "\n")
        f.write(str(kinetics_df[["ddr_delta", "mapk_delta", "pi3k_delta", "vegf_delta"]].describe()))
        f.write("\n\n")
        f.write("KINETIC PATTERNS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"DDR Restoration (Δ > 0.1): {patterns_df['ddr_restoration'].sum()} / {len(patterns_df)}\n")
        f.write(f"MAPK Activation (Δ > 0.1): {patterns_df['mapk_activation'].sum()} / {len(patterns_df)}\n")
        f.write(f"PI3K Activation (Δ > 0.1): {patterns_df['pi3k_activation'].sum()} / {len(patterns_df)}\n")
        f.write(f"VEGF Activation (Δ > 0.1): {patterns_df['vegf_activation'].sum()} / {len(patterns_df)}\n")
        f.write(f"DDR Suppression (Δ < -0.1): {patterns_df['ddr_suppression'].sum()} / {len(patterns_df)}\n")
        f.write(f"MAPK Suppression (Δ < -0.1): {patterns_df['mapk_suppression'].sum()} / {len(patterns_df)}\n")
        f.write("\n\n")
        f.write("CORRELATIONS:\n")
        f.write("-" * 80 + "\n")
        f.write(str(correlation_matrix.round(3)))
        f.write("\n\n")
        f.write("NOTE: Resistance labels not available in GSE165897 metadata.\n")
        f.write("For resistance correlation, external clinical data or publication\n")
        f.write("supplement would be required.\n")
    
    print(f"💾 Saved report: {report_path}\n")
    
    # Step 9: Resistance-stratified analysis (if labels available)
    print("\n" + "="*80)
    print("RESISTANCE-STRATIFIED ANALYSIS")
    print("="*80)
    
    # Load resistance labels (DECIDER cohort: PFI < 6 months = resistant)
    resistance_labels = load_resistance_labels(kinetics_df['patient_id'].unique())
    
    if resistance_labels:
        # Merge resistance labels
        kinetics_df['resistance_label'] = kinetics_df['patient_id'].map(resistance_labels)
        kinetics_df['resistant_binary'] = (kinetics_df['resistance_label'] == 'resistant').astype(int)
        
        # Filter to patients with labels
        labeled_df = kinetics_df[kinetics_df['resistance_label'].notna()].copy()
        
        if len(labeled_df) > 0:
            print(f"\n📊 Resistance Labels: {len(labeled_df)} patients with labels")
            print(f"   Resistant (PFI < 6 months): {(labeled_df['resistant_binary'] == 1).sum()}")
            print(f"   Sensitive (PFI ≥ 6 months): {(labeled_df['resistant_binary'] == 0).sum()}")
            
            # Statistical tests: Compare pathway kinetics between resistant vs sensitive
            print("\n📈 Pathway Kinetics by Resistance Status:")
            from scipy import stats
            
            results = []
            for pathway in ["ddr", "mapk", "pi3k", "vegf"]:
                delta_col = f"{pathway}_delta"
                
                resistant_deltas = labeled_df[labeled_df['resistant_binary'] == 1][delta_col].values
                sensitive_deltas = labeled_df[labeled_df['resistant_binary'] == 0][delta_col].values
                
                if len(resistant_deltas) > 0 and len(sensitive_deltas) > 0:
                    # Mann-Whitney U test (non-parametric, handles small samples)
                    statistic, p_value = stats.mannwhitneyu(resistant_deltas, sensitive_deltas, alternative='two-sided')
                    
                    # Effect size (Cohen's d)
                    pooled_std = np.sqrt((np.var(resistant_deltas) + np.var(sensitive_deltas)) / 2)
                    cohens_d = (np.mean(resistant_deltas) - np.mean(sensitive_deltas)) / pooled_std if pooled_std > 0 else 0.0
                    
                    results.append({
                        "pathway": pathway.upper(),
                        "resistant_mean": np.mean(resistant_deltas),
                        "sensitive_mean": np.mean(sensitive_deltas),
                        "mean_difference": np.mean(resistant_deltas) - np.mean(sensitive_deltas),
                        "p_value": p_value,
                        "cohens_d": cohens_d,
                        "n_resistant": len(resistant_deltas),
                        "n_sensitive": len(sensitive_deltas)
                    })
                    
                    sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
                    print(f"   {pathway.upper()}:")
                    print(f"      Resistant Δ: {np.mean(resistant_deltas):.4f} (n={len(resistant_deltas)})")
                    print(f"      Sensitive Δ: {np.mean(sensitive_deltas):.4f} (n={len(sensitive_deltas)})")
                    print(f"      Difference: {np.mean(resistant_deltas) - np.mean(sensitive_deltas):.4f} (p={p_value:.4f}{sig}, d={cohens_d:.3f})")
            
            # Save resistance-stratified results
            if results:
                resistance_results_df = pd.DataFrame(results)
                resistance_results_df.to_csv(RESULTS_DIR / "resistance_stratified_analysis.csv", index=False)
                print(f"\n💾 Saved resistance-stratified analysis: {RESULTS_DIR / 'resistance_stratified_analysis.csv'}")
            
            # Update kinetics CSV with resistance labels
            kinetics_df.to_csv(RESULTS_DIR / "pathway_kinetics.csv", index=False)
            
            # Generate resistance-stratified visualization
            try:
                import matplotlib.pyplot as plt
                import seaborn as sns
                
                fig, axes = plt.subplots(2, 2, figsize=(14, 10))
                axes = axes.flatten()
                
                for idx, pathway in enumerate(["ddr", "mapk", "pi3k", "vegf"]):
                    delta_col = f"{pathway}_delta"
                    
                    # Box plot
                    data_for_plot = []
                    labels_for_plot = []
                    
                    for label_type in ['resistant', 'sensitive']:
                        binary_val = 1 if label_type == 'resistant' else 0
                        values = labeled_df[labeled_df['resistant_binary'] == binary_val][delta_col].values
                        data_for_plot.extend(values)
                        labels_for_plot.extend([label_type] * len(values))
                    
                    plot_df = pd.DataFrame({
                        'delta': data_for_plot,
                        'resistance': labels_for_plot
                    })
                    
                    sns.boxplot(data=plot_df, x='resistance', y='delta', ax=axes[idx])
                    sns.stripplot(data=plot_df, x='resistance', y='delta', ax=axes[idx], 
                                 color='black', alpha=0.5, size=4)
                    axes[idx].set_title(f"{pathway.upper()} Kinetics\nby Resistance Status", fontweight="bold")
                    axes[idx].set_xlabel("Resistance Status")
                    axes[idx].set_ylabel("Δ (post - pre)")
                    axes[idx].axhline(0, color='gray', linestyle='--', alpha=0.5)
                
                plt.tight_layout()
                plt.savefig(RESULTS_DIR / "resistance_stratified_kinetics.png", dpi=300, bbox_inches="tight")
                print(f"💾 Saved resistance-stratified visualization: {RESULTS_DIR / 'resistance_stratified_kinetics.png'}")
                plt.close()
            except Exception as e:
                print(f"⚠️  Could not generate resistance visualization: {e}")
        else:
            print("⚠️  No patients with resistance labels found")
    else:
        print("⚠️  Resistance labels not available (create resistance_labels.json or resistance_labels.csv)")
        print("   Expected format: {\"EOC1005\": \"resistant\", \"EOC136\": \"sensitive\", ...}")
        print("   Based on DECIDER table: PFI < 6 months = resistant, PFI ≥ 6 months = sensitive")
    
    print("\n" + "="*80)
    print("✅ Pathway kinetics computation complete")
    print("="*80)
    print(f"\n📁 Results saved to: {RESULTS_DIR}")
    print("   - pathway_scores.csv (pre/post scores)")
    print("   - pathway_kinetics.csv (Δ values)")
    print("   - kinetic_patterns.csv (pattern classifications)")
    print("   - pathway_correlations.csv (correlation matrix)")
    print("   - pathway_kinetics_heatmap.png (visualization)")
    print("   - pathway_kinetics_report.txt (summary report)")
    if resistance_labels:
        print("   - resistance_stratified_analysis.csv (resistance-stratified stats)")
        print("   - resistance_stratified_kinetics.png (resistance-stratified visualization)")

def load_resistance_labels(patient_ids: np.ndarray) -> Optional[Dict[str, str]]:
    """
    Load resistance labels for GSE165897 patients.
    
    Sources (in priority order):
    1. resistance_labels.json (manual mapping)
    2. resistance_labels.csv (manual mapping)
    3. DECIDER cohort table (if patient IDs match)
    
    Returns:
        Dict mapping patient_id -> "resistant" or "sensitive" or None if unavailable
    """
    labels = {}
    
    # Try JSON file first
    json_file = RESULTS_DIR.parent / "resistance_labels.json"
    if json_file.exists():
        try:
            with open(json_file, 'r') as f:
                labels = json.load(f)
            print(f"✅ Loaded resistance labels from {json_file}")
            return labels
        except Exception as e:
            print(f"⚠️  Could not load {json_file}: {e}")
    
    # Try CSV file
    csv_file = RESULTS_DIR.parent / "resistance_labels.csv"
    if csv_file.exists():
        try:
            df = pd.read_csv(csv_file)
            if 'patient_id' in df.columns and 'resistance_label' in df.columns:
                labels = dict(zip(df['patient_id'], df['resistance_label']))
                print(f"✅ Loaded resistance labels from {csv_file}")
                return labels
        except Exception as e:
            print(f"⚠️  Could not load {csv_file}: {e}")
    
    # Note: Individual patient PFI values not in provided table
    # User needs to create mapping file with individual patient data
    print(f"⚠️  No resistance label file found. Create one of:")
    print(f"   - {json_file} (format: {{\"EOC1005\": \"resistant\", \"EOC136\": \"sensitive\", ...}})")
    print(f"   - {csv_file} (columns: patient_id, resistance_label)")
    print(f"   Based on DECIDER table: PFI < 6 months = resistant, PFI ≥ 6 months = sensitive")
    
    return None

if __name__ == "__main__":
    main()
