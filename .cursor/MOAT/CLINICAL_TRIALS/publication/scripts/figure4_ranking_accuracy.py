#!/usr/bin/env python3
"""
Generate Figure 4: Ranking Accuracy Comparison

Bar chart comparing Top-3 accuracy and MRR vs baseline targets.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

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

def generate_figure4(output_dir, format="both"):
    """Generate Figure 4: Ranking Accuracy Comparison"""
    
    # Data from validation report
    metrics = ['Top-3 Accuracy', 'MRR']
    baseline = [0.70, 0.65]  # MVP targets
    ours = [1.00, 0.75]  # Actual results
    improvements = [42.9, 15.4]  # Percentage improvements
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    x = np.arange(len(metrics))
    width = 0.35
    
    # Create bars
    bars1 = ax.bar(x - width/2, baseline, width, label='Baseline (Target)', 
                   color='#BDBDBD', edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, ours, width, label='Our Method', 
                   color='#2E7D32', edgecolor='black', linewidth=1.5)
    
    # Add value labels on bars
    for i, (base, our, imp) in enumerate(zip(baseline, ours, improvements)):
        ax.text(i - width/2, base + 0.02, f'{base:.2f}', 
               ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.text(i + width/2, our + 0.02, f'{our:.2f}\n(+{imp:.1f}%)', 
               ha='center', va='bottom', fontsize=10, fontweight='bold', color='#2E7D32')
    
    # Formatting
    ax.set_ylabel('Accuracy Metric', fontsize=12, fontweight='bold')
    ax.set_title('Ranking Accuracy: Mechanism-Based vs Baseline Methods', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1.15)
    ax.set_yticks(np.arange(0, 1.2, 0.2))
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.legend(loc='upper left', fontsize=10)
    
    # Add target line
    ax.axhline(0.70, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    ax.axhline(0.65, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "figure4_ranking_accuracy.pdf", format='pdf', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure4_ranking_accuracy.pdf'}")
    
    if format in ["png", "both"]:
        plt.savefig(output_dir / "figure4_ranking_accuracy.png", format='png', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure4_ranking_accuracy.png'}")
    
    plt.close()

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "../figures"
    format_type = sys.argv[2] if len(sys.argv) > 2 else "both"
    generate_figure4(output_dir, format_type)


