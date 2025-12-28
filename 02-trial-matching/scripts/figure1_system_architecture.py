#!/usr/bin/env python3
"""
Generate Figure 1: System Architecture

Pathway-based mechanism matching system architecture flow diagram.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, ConnectionPatch
from pathlib import Path
import numpy as np

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

def generate_figure1(output_dir, format="both"):
    """Generate Figure 1: System Architecture"""
    
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Title
    ax.text(5, 9.5, 'Pathway-Based Mechanism Matching System Architecture', 
            ha='center', va='top', fontsize=14, fontweight='bold')
    
    # Color scheme
    patient_color = '#E3F2FD'  # Light blue
    trial_color = '#E8F5E9'    # Light green
    process_color = '#FFF3E0'  # Light orange
    output_color = '#F3E5F5'   # Light purple
    border_color = '#424242'   # Dark gray
    
    # === LEFT SIDE: Patient Data ===
    
    # Patient Mutations Box
    patient_box = FancyBboxPatch((0.5, 6.5), 1.8, 1.2,
                                boxstyle="round,pad=0.05",
                                facecolor=patient_color,
                                edgecolor=border_color,
                                linewidth=1.5)
    ax.add_patch(patient_box)
    ax.text(1.4, 7.4, 'Patient\nMutations', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(1.4, 7.0, 'MBD4 p.R361*\nTP53 p.R175H', ha='center', va='center',
            fontsize=9, style='italic')
    
    # Pathway Aggregation Box
    pathway_box = FancyBboxPatch((0.5, 4.5), 1.8, 1.2,
                                boxstyle="round,pad=0.05",
                                facecolor=process_color,
                                edgecolor=border_color,
                                linewidth=1.5)
    ax.add_patch(pathway_box)
    ax.text(1.4, 5.4, 'Pathway\nAggregation', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(1.4, 5.0, 'Gene → Pathway', ha='center', va='center',
            fontsize=9, style='italic')
    
    # 7D Mechanism Vector (Patient)
    vector_box_patient = FancyBboxPatch((0.5, 2.5), 1.8, 1.2,
                                       boxstyle="round,pad=0.05",
                                       facecolor=patient_color,
                                       edgecolor=border_color,
                                       linewidth=1.5)
    ax.add_patch(vector_box_patient)
    ax.text(1.4, 3.4, '7D Mechanism\nVector (Patient)', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(1.4, 3.0, '[0.88, 0.12, 0.05,\n0.02, 0.0, 0.0, 0.0]', ha='center', va='center',
            fontsize=8, family='monospace')
    
    # === RIGHT SIDE: Trial Data ===
    
    # Trial Interventions Box
    trial_box = FancyBboxPatch((7.7, 6.5), 1.8, 1.2,
                               boxstyle="round,pad=0.05",
                               facecolor=trial_color,
                               edgecolor=border_color,
                               linewidth=1.5)
    ax.add_patch(trial_box)
    ax.text(8.6, 7.4, 'Trial\nInterventions', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(8.6, 7.0, 'PARP inhibitor\nATR inhibitor', ha='center', va='center',
            fontsize=9, style='italic')
    
    # MoA Tagging Box
    moa_box = FancyBboxPatch((7.7, 4.5), 1.8, 1.2,
                            boxstyle="round,pad=0.05",
                            facecolor=process_color,
                            edgecolor=border_color,
                            linewidth=1.5)
    ax.add_patch(moa_box)
    ax.text(8.6, 5.4, 'MoA Tagging', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(8.6, 5.0, 'Gemini API', ha='center', va='center',
            fontsize=9, style='italic')
    
    # 7D MoA Vector (Trial)
    vector_box_trial = FancyBboxPatch((7.7, 2.5), 1.8, 1.2,
                                      boxstyle="round,pad=0.05",
                                      facecolor=trial_color,
                                      edgecolor=border_color,
                                      linewidth=1.5)
    ax.add_patch(vector_box_trial)
    ax.text(8.6, 3.4, '7D MoA Vector\n(Trial)', ha='center', va='center',
            fontsize=10, fontweight='bold')
    ax.text(8.6, 3.0, '[0.95, 0.10, 0.05,\n0.02, 0.0, 0.0, 0.0]', ha='center', va='center',
            fontsize=8, family='monospace')
    
    # === CENTER: Processing ===
    
    # Cosine Similarity Box
    cosine_box = FancyBboxPatch((4.0, 2.0), 2.0, 1.5,
                               boxstyle="round,pad=0.05",
                               facecolor=process_color,
                               edgecolor=border_color,
                               linewidth=2.0)
    ax.add_patch(cosine_box)
    ax.text(5.0, 3.0, 'Cosine Similarity', ha='center', va='center',
            fontsize=11, fontweight='bold')
    ax.text(5.0, 2.6, 'Mechanism Fit', ha='center', va='center',
            fontsize=10, style='italic')
    ax.text(5.0, 2.3, '= 0.983', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#2E7D32')
    
    # Combined Scoring Box
    combined_box = FancyBboxPatch((4.0, 0.5), 2.0, 1.0,
                                 boxstyle="round,pad=0.05",
                                 facecolor=process_color,
                                 edgecolor=border_color,
                                 linewidth=2.0)
    ax.add_patch(combined_box)
    ax.text(5.0, 1.2, 'Combined Scoring', ha='center', va='center',
            fontsize=11, fontweight='bold')
    ax.text(5.0, 0.8, '0.7×eligibility +', ha='center', va='center',
            fontsize=9)
    ax.text(5.0, 0.5, '0.3×mechanism_fit', ha='center', va='center',
            fontsize=9)
    
    # Ranked Trial List Box
    output_box = FancyBboxPatch((4.0, -0.5), 2.0, 0.8,
                                boxstyle="round,pad=0.05",
                                facecolor=output_color,
                                edgecolor=border_color,
                                linewidth=2.0)
    ax.add_patch(output_box)
    ax.text(5.0, 0.1, 'Ranked Trial List', ha='center', va='center',
            fontsize=11, fontweight='bold')
    ax.text(5.0, -0.2, 'Top 5-12 trials', ha='center', va='center',
            fontsize=9, style='italic')
    
    # === ARROWS ===
    
    # Patient side arrows
    arrow1 = FancyArrowPatch((1.4, 6.5), (1.4, 5.7),
                             arrowstyle='->', lw=2, color='#1976D2')
    ax.add_patch(arrow1)
    ax.text(1.7, 6.1, 'Pathway\nAggregation', ha='left', va='center',
            fontsize=8, color='#1976D2')
    
    arrow2 = FancyArrowPatch((1.4, 4.5), (1.4, 3.7),
                             arrowstyle='->', lw=2, color='#1976D2')
    ax.add_patch(arrow2)
    
    # Trial side arrows
    arrow3 = FancyArrowPatch((8.6, 6.5), (8.6, 5.7),
                             arrowstyle='->', lw=2, color='#2E7D32')
    ax.add_patch(arrow3)
    ax.text(8.3, 6.1, 'MoA\nTagging', ha='right', va='center',
            fontsize=8, color='#2E7D32')
    
    arrow4 = FancyArrowPatch((8.6, 4.5), (8.6, 3.7),
                             arrowstyle='->', lw=2, color='#2E7D32')
    ax.add_patch(arrow4)
    
    # Center arrows (from vectors to cosine similarity)
    arrow5 = FancyArrowPatch((2.3, 3.1), (4.0, 2.75),
                             arrowstyle='->', lw=2, color='#1976D2')
    ax.add_patch(arrow5)
    
    arrow6 = FancyArrowPatch((7.7, 3.1), (6.0, 2.75),
                             arrowstyle='->', lw=2, color='#2E7D32')
    ax.add_patch(arrow6)
    
    # Arrow from cosine to combined scoring
    arrow7 = FancyArrowPatch((5.0, 2.0), (5.0, 1.5),
                             arrowstyle='->', lw=2, color='#F57C00')
    ax.add_patch(arrow7)
    
    # Arrow from combined to output
    arrow8 = FancyArrowPatch((5.0, 0.5), (5.0, 0.3),
                             arrowstyle='->', lw=2, color='#7B1FA2')
    ax.add_patch(arrow8)
    
    # === ANNOTATIONS ===
    
    # 7D Vector annotation
    ax.text(5.0, 4.5, '7D Vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]',
            ha='center', va='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "figure1_system_architecture.pdf", format='pdf', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure1_system_architecture.pdf'}")
    
    if format in ["png", "both"]:
        plt.savefig(output_dir / "figure1_system_architecture.png", format='png', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure1_system_architecture.png'}")
    
    plt.close()

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "../figures"
    format_type = sys.argv[2] if len(sys.argv) > 2 else "both"
    generate_figure1(output_dir, format_type)


