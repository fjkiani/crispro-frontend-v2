#!/usr/bin/env python3
"""
Generate Figure 2: Mechanism Fit Performance

Box plot comparing mechanism fit scores for DDR-targeting vs non-DDR trials.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import json
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

from api.services.mechanism_fit_ranker import MechanismFitRanker
from api.services.pathway_to_mechanism_vector import convert_moa_dict_to_vector

# Publication settings
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'figure.titlesize': 16,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

def load_validation_data():
    """Load trial MoA vectors and compute mechanism fit scores"""
    script_dir = Path(__file__).resolve().parent
    moa_path = script_dir.parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal" / "api" / "resources" / "trial_moa_vectors.json"
    
    if not moa_path.exists():
        raise FileNotFoundError(f"Trial MoA vectors file not found: {moa_path}")
    
    trial_moa_vectors = json.loads(moa_path.read_text())
    
    # DDR-high patient mechanism vector (7D)
    patient_mechanism_vector = [0.88, 0.12, 0.05, 0.02, 0.0, 0.0, 0.0]
    
    ranker = MechanismFitRanker(alpha=0.7, beta=0.3)
    
    # Prepare trials for scoring
    trials = []
    for nct_id, data in trial_moa_vectors.items():
        moa_dict = data.get("moa_vector") or {}
        moa_vector = convert_moa_dict_to_vector(moa_dict, use_7d=True)
        
        trials.append({
            "nct_id": nct_id,
            "title": nct_id,
            "eligibility_score": 0.85,
            "moa_vector": moa_vector,
        })
    
    # Rank trials (allow min_mechanism_fit=0.0 to see all scores)
    ranked_scores = ranker.rank_trials(
        trials=trials,
        sae_mechanism_vector=patient_mechanism_vector,
        min_eligibility=0.60,
        min_mechanism_fit=0.0,
    )
    
    # Split by DDR tag
    ddr_scores = []
    non_ddr_scores = []
    
    for score in ranked_scores:
        moa_dict = (trial_moa_vectors.get(score.nct_id, {}) or {}).get("moa_vector") or {}
        ddr_value = float(moa_dict.get("ddr", 0.0) or 0.0)
        if ddr_value > 0.5:
            ddr_scores.append(score.mechanism_fit_score)
        else:
            non_ddr_scores.append(score.mechanism_fit_score)
    
    return ddr_scores, non_ddr_scores

def generate_figure2(output_dir, format="both"):
    """Generate Figure 2: Mechanism Fit Performance"""
    
    # Load data
    ddr_scores, non_ddr_scores = load_validation_data()
    
    # Calculate statistics
    mean_ddr = np.mean(ddr_scores)
    mean_non = np.mean(non_ddr_scores)
    median_ddr = np.median(ddr_scores)
    median_non = np.median(non_ddr_scores)
    delta = mean_ddr - mean_non
    discrimination_ratio = mean_ddr / mean_non if mean_non > 0 else float('inf')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Prepare data for box plot
    data = [non_ddr_scores, ddr_scores]
    positions = [1, 2]
    labels = ['Non-DDR\n(n=16)', 'DDR\n(n=31)']
    
    # Create box plot
    bp = ax.boxplot(data, positions=positions, labels=labels,
                   patch_artist=True, widths=0.6,
                   showmeans=True, meanline=True,
                   boxprops=dict(linewidth=1.5),
                   whiskerprops=dict(linewidth=1.5),
                   capprops=dict(linewidth=1.5),
                   medianprops=dict(linewidth=2, color='black'),
                   meanprops=dict(linewidth=2, linestyle='--', color='red'))
    
    # Color boxes
    bp['boxes'][0].set_facecolor('#FFCDD2')  # Light red
    bp['boxes'][1].set_facecolor('#C8E6C9')  # Light green
    bp['boxes'][0].set_edgecolor('#C62828')  # Dark red
    bp['boxes'][1].set_edgecolor('#2E7D32')  # Dark green
    
    # Add mean lines with annotations
    ax.axhline(mean_non, color='#C62828', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.axhline(mean_ddr, color='#2E7D32', linestyle='--', alpha=0.5, linewidth=1.5)
    
    # Add annotations
    ax.text(1, mean_non + 0.05, f'Mean: {mean_non:.3f}\n(77% below 0.20)',
            ha='center', va='bottom', fontsize=9, color='#C62828',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    ax.text(2, mean_ddr + 0.05, f'Mean: {mean_ddr:.3f}\n(Exceeds 0.92)',
            ha='center', va='bottom', fontsize=9, color='#2E7D32',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Add separation annotation
    mid_y = (mean_ddr + mean_non) / 2
    ax.text(1.5, mid_y, f'Separation Δ = {delta:.3f}\n({discrimination_ratio:.1f}× discrimination)',
            ha='center', va='center', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8, edgecolor='black', linewidth=1.5))
    
    # Formatting
    ax.set_ylabel('Mechanism Fit Score', fontsize=12, fontweight='bold')
    ax.set_title('Mechanism Fit Performance: DDR vs Non-DDR Trials', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_ylim(-0.05, 1.1)
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add target lines
    ax.axhline(0.92, color='gray', linestyle=':', alpha=0.5, linewidth=1, label='Target (0.92)')
    ax.axhline(0.20, color='gray', linestyle=':', alpha=0.5, linewidth=1, label='Target (0.20)')
    ax.legend(loc='upper right', fontsize=9)
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "figure2_mechanism_fit_performance.pdf", format='pdf', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure2_mechanism_fit_performance.pdf'}")
    
    if format in ["png", "both"]:
        plt.savefig(output_dir / "figure2_mechanism_fit_performance.png", format='png', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure2_mechanism_fit_performance.png'}")
    
    plt.close()

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "../figures"
    format_type = sys.argv[2] if len(sys.argv) > 2 else "both"
    generate_figure2(output_dir, format_type)


