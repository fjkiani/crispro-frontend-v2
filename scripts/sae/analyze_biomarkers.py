#!/usr/bin/env python3
"""
⚔️ SAE BIOMARKER ANALYSIS - TCGA-OV PLATINUM RESPONSE
=====================================================

Agent: Zo (Lead Commander)
Date: January 15, 2025
Priority: HIGH
Sprint: 1 - Task 1

Objective: Execute biomarker correlation analysis on SAE cohort data.
This script loads the SAE features extracted by extract_sae_features_cohort.py,
runs the biomarker_correlation_service to compute statistical correlations,
and generates visualization plots for top features.

Inputs:
- data/validation/sae_cohort/sae_features_tcga_ov_platinum.json

Outputs:
- data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json
- data/validation/sae_cohort/plots/ (correlation plots, distribution plots, CV stability)

Success Criteria:
- Biomarker analysis completes without errors
- Top 100 features identified with p < 0.01
- Plots generated for top 10 features
- Results saved for downstream interpretation
"""

import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from typing import Dict
from dataclasses import asdict
from loguru import logger

# Add api directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

from api.services.biomarker_correlation_service import BiomarkerCorrelationService

# --- Configuration ---
OUTPUT_DIR = Path("data/validation/sae_cohort")
PLOTS_DIR = OUTPUT_DIR / "plots"
LOG_FILE = OUTPUT_DIR / "analyze_biomarkers_log.log"

# --- Logging Setup ---
logger.remove()
logger.add(sys.stderr, level="INFO")
logger.add(LOG_FILE, rotation="10 MB", level="DEBUG")

# --- Plotting Functions ---
def plot_correlation_distribution(summary: Dict, output_path: Path):
    """Plots distribution of Pearson correlations."""
    all_correlations = [fc["pearson_r"] for fc in summary["top_features"]]
    
    plt.figure(figsize=(10, 6))
    plt.hist(all_correlations, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel("Pearson Correlation (r)")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of SAE Feature Correlations with Platinum Response\n(Top {len(all_correlations)} features)")
    plt.axvline(0, color='red', linestyle='--', alpha=0.5, label='r=0 (no correlation)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    logger.info(f"Saved correlation distribution plot to {output_path}")

def plot_top_features_barplot(summary: Dict, output_path: Path, top_n: int = 10):
    """Plots top N features by absolute correlation."""
    top_features = summary["top_features"][:top_n]
    
    indices = [f"F{fc['feature_index']}" for fc in top_features]
    correlations = [fc["pearson_r"] for fc in top_features]
    ci_lower = [fc["bootstrap_ci_lower"] for fc in top_features]
    ci_upper = [fc["bootstrap_ci_upper"] for fc in top_features]
    
    # Compute error bars
    errors = [
        [correlations[i] - ci_lower[i], ci_upper[i] - correlations[i]]
        for i in range(len(correlations))
    ]
    errors = np.array(errors).T
    
    plt.figure(figsize=(12, 6))
    x_pos = np.arange(len(indices))
    colors = ['green' if r > 0 else 'red' for r in correlations]
    
    plt.bar(x_pos, correlations, yerr=errors, alpha=0.7, color=colors, edgecolor='black')
    plt.xticks(x_pos, indices, rotation=45, ha='right')
    plt.xlabel("SAE Feature Index")
    plt.ylabel("Pearson Correlation (r)")
    plt.title(f"Top {top_n} SAE Features Correlated with Platinum Response")
    plt.axhline(0, color='black', linestyle='-', linewidth=0.5)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    logger.info(f"Saved top features barplot to {output_path}")

def plot_cv_stability(summary: Dict, output_path: Path, top_n: int = 20):
    """Plots CV stability for top features."""
    top_features = summary["top_features"][:top_n]
    
    indices = [f"F{fc['feature_index']}" for fc in top_features]
    stability = [fc["cv_stability"] for fc in top_features]
    
    plt.figure(figsize=(12, 6))
    x_pos = np.arange(len(indices))
    colors = ['green' if s >= 0.6 else 'orange' for s in stability]
    
    plt.bar(x_pos, stability, alpha=0.7, color=colors, edgecolor='black')
    plt.xticks(x_pos, indices, rotation=45, ha='right')
    plt.xlabel("SAE Feature Index")
    plt.ylabel("CV Stability (fraction of folds significant)")
    plt.title(f"Cross-Validation Stability for Top {top_n} SAE Features")
    plt.axhline(0.6, color='red', linestyle='--', alpha=0.5, label='Stability Threshold (0.6)')
    plt.ylim(0, 1.0)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    logger.info(f"Saved CV stability plot to {output_path}")

def plot_effect_sizes(summary: Dict, output_path: Path, top_n: int = 20):
    """Plots Cohen's d effect sizes for top features."""
    top_features = summary["top_features"][:top_n]
    
    indices = [f"F{fc['feature_index']}" for fc in top_features]
    effect_sizes = [fc["cohen_d"] for fc in top_features]
    
    plt.figure(figsize=(12, 6))
    x_pos = np.arange(len(indices))
    colors = ['green' if d >= 0.5 else ('orange' if d >= 0.3 else 'gray') for d in effect_sizes]
    
    plt.bar(x_pos, effect_sizes, alpha=0.7, color=colors, edgecolor='black')
    plt.xticks(x_pos, indices, rotation=45, ha='right')
    plt.xlabel("SAE Feature Index")
    plt.ylabel("Cohen's d (Effect Size)")
    plt.title(f"Effect Sizes for Top {top_n} SAE Features (Sensitive vs Refractory)")
    plt.axhline(0.3, color='orange', linestyle='--', alpha=0.5, label='Small Effect (0.3)')
    plt.axhline(0.5, color='green', linestyle='--', alpha=0.5, label='Medium Effect (0.5)')
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    logger.info(f"Saved effect sizes plot to {output_path}")

def generate_summary_table(summary: Dict, output_path: Path):
    """Generates markdown summary table of top features."""
    top_features = summary["top_features"][:20]  # Top 20 for table
    
    with open(output_path, 'w') as f:
        f.write("# SAE Biomarker Discovery Summary\n\n")
        f.write(f"**Cohort**: TCGA-OV Platinum Response ({summary['cohort_size']} patients)\n\n")
        f.write(f"**Outcome Distribution**: {summary['outcome_distribution']}\n\n")
        f.write(f"**Total Features Analyzed**: {summary['total_features_analyzed']}\n\n")
        f.write(f"**Significant Features (p<{summary['p_value_threshold']}, d>={summary['effect_size_threshold']}, CV>={summary['cv_stability_threshold']})**:  {summary['significant_features_count']}\n\n")
        f.write(f"**Correction Method**: {summary['correction_method']}\n\n")
        f.write("---\n\n")
        f.write("## Top 20 SAE Features\n\n")
        f.write("| Rank | Feature | Pearson r | P-value | Cohen's d | CV Stability | 95% CI |\n")
        f.write("|------|---------|-----------|---------|-----------|--------------|--------|\n")
        
        for fc in top_features:
            f.write(
                f"| {fc['rank']} | {fc['feature_index']} | "
                f"{fc['pearson_r']:.3f} | {fc['pearson_p']:.2e} | "
                f"{fc['cohen_d']:.3f} | {fc['cv_stability']:.2f} | "
                f"[{fc['bootstrap_ci_lower']:.3f}, {fc['bootstrap_ci_upper']:.3f}] |\n"
            )
        
        f.write("\n---\n\n")
        f.write(f"**Generated**: {summary['provenance']['timestamp']}\n\n")
        f.write(f"**Script**: {summary['provenance']['script']}\n\n")
        f.write("⚠️ **RUO DISCLAIMER**: This analysis is for research and validation only.\n")
    
    logger.info(f"Saved summary table to {output_path}")

# --- Main Execution ---
def main():
    """Main analysis execution."""
    import argparse
    
    parser = argparse.ArgumentParser(description="SAE Biomarker Correlation Analysis")
    parser.add_argument("--input", type=str, help="Input SAE cohort file (overrides default)")
    parser.add_argument("--output", type=str, help="Output biomarker file (overrides default)")
    parser.add_argument("--plots-dir", type=str, help="Output plots directory (overrides default)")
    args = parser.parse_args()
    
    # Override paths if provided
    if args.input:
        from api.services.biomarker_correlation_service import SAE_COHORT_FILE
        import api.services.biomarker_correlation_service as bcs
        bcs.SAE_COHORT_FILE = Path(args.input)
        logger.info(f"Using custom input file: {args.input}")
    
    output_file = Path(args.output) if args.output else OUTPUT_DIR / "sae_tcga_ov_platinum_biomarkers.json"
    plots_dir = Path(args.plots_dir) if args.plots_dir else PLOTS_DIR
    
    logger.info("⚔️ STARTING SAE BIOMARKER ANALYSIS")
    logger.info("="*80)
    
    # Create output directories
    plots_dir.mkdir(parents=True, exist_ok=True)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize service
    logger.info("Initializing BiomarkerCorrelationService...")
    service = BiomarkerCorrelationService()
    
    # Run analysis
    logger.info("Running biomarker correlation analysis...")
    summary = service.run_analysis()
    
    # Convert dataclass to plain dict
    summary_dict: Dict = asdict(summary)
    
    # Save results
    with open(output_file, 'w') as f:
        json.dump(summary_dict, f, indent=2)
    logger.info(f"✅ Biomarker analysis saved to {output_file}")
    
    # Generate plots
    logger.info("Generating visualization plots...")
    
    plot_correlation_distribution(
        summary_dict,
        plots_dir / "correlation_distribution.png"
    )
    
    plot_top_features_barplot(
        summary_dict,
        plots_dir / "top_features_barplot.png",
        top_n=10
    )
    
    plot_cv_stability(
        summary_dict,
        plots_dir / "cv_stability.png",
        top_n=20
    )
    
    plot_effect_sizes(
        summary_dict,
        plots_dir / "effect_sizes.png",
        top_n=20
    )
    
    # Generate summary table
    logger.info("Generating summary table...")
    summary_md_path = output_file.parent / "biomarker_summary.md"
    generate_summary_table(
        summary_dict,
        summary_md_path
    )
    
    # Print summary to console
    print("\n" + "="*80)
    print("⚔️ SAE BIOMARKER ANALYSIS COMPLETE")
    print("="*80)
    print(f"Cohort Size: {summary_dict['cohort_size']}")
    print(f"Outcome Distribution: {summary_dict['outcome_distribution']}")
    print(f"Total Features Analyzed: {summary_dict['total_features_analyzed']}")
    print(f"Significant Features: {summary_dict['significant_features_count']}")
    print(f"\nTop 10 Features:")
    print("-"*80)
    
    for i, fc in enumerate(summary_dict['top_features'][:10], 1):
        print(
            f"{i}. Feature {fc['feature_index']}: "
            f"r={fc['pearson_r']:.3f} (p={fc['pearson_p']:.2e}), "
            f"d={fc['cohen_d']:.3f}, "
            f"CV={fc['cv_stability']:.2f}, "
            f"95%CI=[{fc['bootstrap_ci_lower']:.3f}, {fc['bootstrap_ci_upper']:.3f}]"
        )
    
    print("="*80)
    print(f"\n📊 Plots saved to: {plots_dir}")
    print(f"📄 Summary table: {summary_md_path}")
    print(f"📋 Full results: {output_file}")
    print("\n⚠️  RUO DISCLAIMER: Results for validation only. Manager approval required.")
    print("="*80)

if __name__ == "__main__":
    main()

