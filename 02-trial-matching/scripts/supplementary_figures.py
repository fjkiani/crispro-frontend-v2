#!/usr/bin/env python3
"""
Generate Supplementary Figures

All supplementary figures for the publication.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import json
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

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

def generate_supp_figure1(output_dir, format="both"):
    """Supplementary Figure 1: Pathway Aggregation Algorithm"""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Title
    ax.text(5, 7.5, 'Pathway Aggregation Algorithm', 
            ha='center', va='top', fontsize=14, fontweight='bold')
    
    # Gene mutations
    genes_box = plt.Rectangle((1, 5), 1.5, 1, facecolor='#E3F2FD', 
                              edgecolor='black', linewidth=1.5)
    ax.add_patch(genes_box)
    ax.text(1.75, 5.5, 'Gene\nMutations', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(1.75, 5.1, 'MBD4, TP53', ha='center', va='center', fontsize=9)
    
    # Arrow
    ax.arrow(2.5, 5.5, 1, 0, head_width=0.2, head_length=0.2, 
            fc='black', ec='black', linewidth=2)
    ax.text(3, 5.8, 'Gene→Pathway\nMapping', ha='center', va='bottom',
            fontsize=9, style='italic')
    
    # Pathway burden
    pathway_box = plt.Rectangle((3.5, 5), 1.5, 1, facecolor='#FFF3E0',
                               edgecolor='black', linewidth=1.5)
    ax.add_patch(pathway_box)
    ax.text(4.25, 5.5, 'Pathway\nBurden', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(4.25, 5.1, 'DDR: 0.88', ha='center', va='center', fontsize=9)
    
    # Arrow
    ax.arrow(5, 5.5, 1, 0, head_width=0.2, head_length=0.2,
            fc='black', ec='black', linewidth=2)
    ax.text(5.5, 5.8, 'Vector\nConstruction', ha='center', va='bottom',
            fontsize=9, style='italic')
    
    # 7D Vector
    vector_box = plt.Rectangle((6, 5), 1.5, 1, facecolor='#E8F5E9',
                              edgecolor='black', linewidth=1.5)
    ax.add_patch(vector_box)
    ax.text(6.75, 5.5, '7D Vector', ha='center', va='center',
            fontsize=10, fontweight='bold')
    
    # Pathway mapping table
    pathways = ['DDR', 'MAPK', 'PI3K', 'VEGF', 'HER2', 'IO', 'Efflux']
    genes = ['MBD4, TP53', 'KRAS, BRAF', 'PIK3CA, PTEN', 'VEGFA', 'ERBB2', 'TMB≥20', 'ABCB1']
    values = [0.88, 0.12, 0.05, 0.02, 0.0, 0.0, 0.0]
    
    y_start = 3.5
    for i, (pathway, gene, value) in enumerate(zip(pathways, genes, values)):
        y_pos = y_start - i * 0.4
        ax.text(1, y_pos, f'{pathway}:', ha='left', va='center', 
               fontsize=9, fontweight='bold')
        ax.text(2.5, y_pos, gene, ha='left', va='center', fontsize=9)
        ax.text(6, y_pos, f'→ {value:.2f}', ha='left', va='center', 
               fontsize=9, fontweight='bold', color='#2E7D32')
    
    # Save
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "supp_figure1_pathway_aggregation.pdf", format='pdf', bbox_inches='tight')
    if format in ["png", "both"]:
        plt.savefig(output_dir / "supp_figure1_pathway_aggregation.png", format='png', bbox_inches='tight')
    
    plt.close()

def generate_supp_figure2(output_dir, format="both"):
    """Supplementary Figure 2: Mechanism Fit Ranking Algorithm"""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Title
    ax.text(5, 7.5, 'Mechanism Fit Ranking Algorithm', 
            ha='center', va='top', fontsize=14, fontweight='bold')
    
    # Patient vector
    patient_box = plt.Rectangle((1, 5.5), 1.5, 1, facecolor='#E3F2FD',
                                edgecolor='black', linewidth=1.5)
    ax.add_patch(patient_box)
    ax.text(1.75, 6, 'Patient\nVector', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(1.75, 5.6, '[0.88, 0.12, ...]', ha='center', va='center',
            fontsize=8, family='monospace')
    
    # Trial vector
    trial_box = plt.Rectangle((7.5, 5.5), 1.5, 1, facecolor='#E8F5E9',
                             edgecolor='black', linewidth=1.5)
    ax.add_patch(trial_box)
    ax.text(8.25, 6, 'Trial MoA\nVector', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(8.25, 5.6, '[0.95, 0.10, ...]', ha='center', va='center',
            fontsize=8, family='monospace')
    
    # Arrows to cosine
    ax.arrow(2.5, 6, 2, -0.5, head_width=0.15, head_length=0.15,
            fc='black', ec='black', linewidth=2)
    ax.arrow(7.5, 6, -2, -0.5, head_width=0.15, head_length=0.15,
            fc='black', ec='black', linewidth=2)
    
    # Cosine similarity
    cosine_box = plt.Rectangle((4, 3.5), 2, 1, facecolor='#FFF3E0',
                              edgecolor='black', linewidth=2)
    ax.add_patch(cosine_box)
    ax.text(5, 4.2, 'Cosine Similarity', ha='center', va='center',
            fontsize=11, fontweight='bold')
    ax.text(5, 3.8, 'mechanism_fit = 0.983', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#2E7D32')
    
    # Formula
    formula_text = (
        "mechanism_fit = cosine(patient_vector, trial_moa_vector)\n"
        "combined_score = 0.7 × eligibility + 0.3 × mechanism_fit"
    )
    ax.text(5, 2.5, formula_text, ha='center', va='center',
            fontsize=10, family='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
    
    # Save
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "supp_figure2_mechanism_fit_algorithm.pdf", format='pdf', bbox_inches='tight')
    if format in ["png", "both"]:
        plt.savefig(output_dir / "supp_figure2_mechanism_fit_algorithm.png", format='png', bbox_inches='tight')
    
    plt.close()

def generate_supp_figure3(output_dir, format="both"):
    """Supplementary Figure 3: Trial Coverage by Pathway"""
    
    # Load trial data
    script_dir = Path(__file__).resolve().parent
    moa_path = script_dir.parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal" / "api" / "resources" / "trial_moa_vectors.json"
    
    if not moa_path.exists():
        print(f"  ⚠️ Trial MoA vectors file not found, using placeholder data")
        pathway_counts = {'DDR': 31, 'MAPK': 6, 'VEGF': 3, 'HER2': 3, 'IO': 6}
    else:
        trial_moa_vectors = json.loads(moa_path.read_text())
        
        pathway_counts = {'DDR': 0, 'MAPK': 0, 'VEGF': 0, 'HER2': 0, 'IO': 0, 'PI3K': 0, 'Efflux': 0}
        
        for nct_id, data in trial_moa_vectors.items():
            moa_dict = data.get("moa_vector", {})
            if float(moa_dict.get("ddr", 0.0) or 0.0) > 0.5:
                pathway_counts['DDR'] += 1
            if float(moa_dict.get("mapk", 0.0) or 0.0) > 0.5:
                pathway_counts['MAPK'] += 1
            if float(moa_dict.get("vegf", 0.0) or 0.0) > 0.5:
                pathway_counts['VEGF'] += 1
            if float(moa_dict.get("her2", 0.0) or 0.0) > 0.5:
                pathway_counts['HER2'] += 1
            if float(moa_dict.get("io", 0.0) or 0.0) > 0.5:
                pathway_counts['IO'] += 1
    
    pathways = list(pathway_counts.keys())
    counts = list(pathway_counts.values())
    colors = ['#D32F2F', '#F57C00', '#7B1FA2', '#1976D2', '#C2185B', '#00796B', '#5D4037']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    bars = ax.bar(pathways, counts, color=colors[:len(pathways)], 
                  edgecolor='black', linewidth=1.5)
    
    # Add value labels
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
               f'{count}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylabel('Number of Trials', fontsize=12, fontweight='bold')
    ax.set_title('Trial MoA Coverage by Pathway', fontsize=14, fontweight='bold', pad=20)
    ax.set_ylim(0, max(counts) * 1.2)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add total annotation
    total = sum(counts)
    ax.text(0.02, 0.98, f'Total Trials Tagged: {total}',
            transform=ax.transAxes, fontsize=10, fontweight='bold',
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5',
            facecolor='wheat', alpha=0.8))
    
    # Save
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "supp_figure3_trial_coverage.pdf", format='pdf', bbox_inches='tight')
    if format in ["png", "both"]:
        plt.savefig(output_dir / "supp_figure3_trial_coverage.png", format='png', bbox_inches='tight')
    
    plt.close()

def generate_all_supplementary(output_dir, format="both"):
    """Generate all supplementary figures"""
    
    print("  Supplementary Figure 1: Pathway Aggregation Algorithm...")
    generate_supp_figure1(output_dir, format)
    print("    ✅ Complete")
    
    print("  Supplementary Figure 2: Mechanism Fit Ranking Algorithm...")
    generate_supp_figure2(output_dir, format)
    print("    ✅ Complete")
    
    print("  Supplementary Figure 3: Trial Coverage by Pathway...")
    generate_supp_figure3(output_dir, format)
    print("    ✅ Complete")

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "../figures"
    format_type = sys.argv[2] if len(sys.argv) > 2 else "both"
    generate_all_supplementary(output_dir, format_type)


