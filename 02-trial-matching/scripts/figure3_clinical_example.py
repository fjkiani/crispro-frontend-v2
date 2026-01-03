#!/usr/bin/env python3
"""
Generate Figure 3: Clinical Example (MBD4+TP53)

Multi-panel figure showing patient pathway burden, top ranked trials, and mechanism alignment.
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
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

def load_trial_data():
    """Load trial data and compute rankings"""
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
    
    # Rank trials
    ranked_scores = ranker.rank_trials(
        trials=trials,
        sae_mechanism_vector=patient_mechanism_vector,
        min_eligibility=0.60,
        min_mechanism_fit=0.50,
    )
    
    # Get top 5 trials
    top_5 = ranked_scores[:5]
    
    return patient_mechanism_vector, top_5, trial_moa_vectors

def generate_figure3(output_dir, format="both"):
    """Generate Figure 3: Clinical Example (MBD4+TP53)"""
    
    # Load data
    patient_vector, top_5_trials, trial_moa_vectors = load_trial_data()
    
    # Create 3-panel figure
    fig = plt.figure(figsize=(15, 5))
    
    pathways = ['DDR', 'MAPK', 'PI3K', 'VEGF', 'HER2', 'IO', 'Efflux']
    pathway_colors = ['#D32F2F', '#F57C00', '#7B1FA2', '#1976D2', '#C2185B', '#00796B', '#5D4037']
    
    # Panel A: Patient Mechanism Vector
    ax1 = plt.subplot(1, 3, 1)
    bars1 = ax1.barh(pathways, patient_vector, color=pathway_colors, edgecolor='black', linewidth=1.5)
    ax1.set_xlabel('Pathway Burden', fontsize=11, fontweight='bold')
    ax1.set_title('(A) Patient Pathway Burden', fontsize=12, fontweight='bold', pad=10)
    ax1.set_xlim(0, 1.0)
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    
    # Add value labels
    for i, (pathway, value) in enumerate(zip(pathways, patient_vector)):
        if value > 0:
            ax1.text(value + 0.02, i, f'{value:.2f}', va='center', fontsize=9, fontweight='bold')
    
    # Highlight DDR (highest)
    bars1[0].set_edgecolor('black')
    bars1[0].set_linewidth(2.5)
    
    # Panel B: Top 5 Ranked Trials
    ax2 = plt.subplot(1, 3, 2)
    
    trial_ids = [trial.nct_id for trial in top_5_trials]
    mechanism_fits = [trial.mechanism_fit_score for trial in top_5_trials]
    combined_scores = [trial.combined_score for trial in top_5_trials]
    
    x = np.arange(len(trial_ids))
    width = 0.35
    
    bars2a = ax2.barh(x - width/2, mechanism_fits, width, label='Mechanism Fit', 
                      color='#2E7D32', edgecolor='black', linewidth=1.5)
    bars2b = ax2.barh(x + width/2, combined_scores, width, label='Combined Score',
                      color='#1976D2', edgecolor='black', linewidth=1.5)
    
    ax2.set_yticks(x)
    ax2.set_yticklabels([nct[:12] + '...' if len(nct) > 12 else nct for nct in trial_ids], fontsize=9)
    ax2.set_xlabel('Score', fontsize=11, fontweight='bold')
    ax2.set_title('(B) Top 5 Ranked Trials', fontsize=12, fontweight='bold', pad=10)
    ax2.legend(loc='lower right', fontsize=9)
    ax2.set_xlim(0, 1.0)
    ax2.grid(axis='x', alpha=0.3, linestyle='--')
    
    # Add value labels
    for i, (mf, cs) in enumerate(zip(mechanism_fits, combined_scores)):
        ax2.text(mf + 0.02, i - width/2, f'{mf:.3f}', va='center', fontsize=8, fontweight='bold')
        ax2.text(cs + 0.02, i + width/2, f'{cs:.3f}', va='center', fontsize=8, fontweight='bold')
    
    # Panel C: Mechanism Alignment Breakdown (Top Trial)
    ax3 = plt.subplot(1, 3, 3)
    
    # Get top trial's MoA vector
    top_trial_id = top_5_trials[0].nct_id
    top_trial_moa = trial_moa_vectors.get(top_trial_id, {}).get("moa_vector", {})
    
    # Compute alignment (element-wise product, normalized)
    alignment = [p * t for p, t in zip(patient_vector, [
        float(top_trial_moa.get("ddr", 0.0) or 0.0),
        float(top_trial_moa.get("mapk", 0.0) or 0.0),
        float(top_trial_moa.get("pi3k", 0.0) or 0.0),
        float(top_trial_moa.get("vegf", 0.0) or 0.0),
        float(top_trial_moa.get("her2", 0.0) or 0.0),
        float(top_trial_moa.get("io", 0.0) or 0.0),
        float(top_trial_moa.get("efflux", 0.0) or 0.0),
    ])]
    
    bars3 = ax3.barh(pathways, alignment, color=pathway_colors, edgecolor='black', linewidth=1.5)
    ax3.set_xlabel('Alignment Score', fontsize=11, fontweight='bold')
    ax3.set_title(f'(C) Per-Pathway Alignment\n(Top Trial: {top_trial_id[:12]}...)', 
                 fontsize=12, fontweight='bold', pad=10)
    ax3.set_xlim(0, 1.0)
    ax3.grid(axis='x', alpha=0.3, linestyle='--')
    
    # Add value labels
    for i, (pathway, value) in enumerate(zip(pathways, alignment)):
        if value > 0.01:
            ax3.text(value + 0.02, i, f'{value:.2f}', va='center', fontsize=9, fontweight='bold')
    
    # Highlight DDR (highest alignment)
    if alignment[0] > 0:
        bars3[0].set_edgecolor('black')
        bars3[0].set_linewidth(2.5)
    
    plt.tight_layout()
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "figure3_clinical_example.pdf", format='pdf', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure3_clinical_example.pdf'}")
    
    if format in ["png", "both"]:
        plt.savefig(output_dir / "figure3_clinical_example.png", format='png', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure3_clinical_example.png'}")
    
    plt.close()

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "../figures"
    format_type = sys.argv[2] if len(sys.argv) > 2 else "both"
    generate_figure3(output_dir, format_type)


