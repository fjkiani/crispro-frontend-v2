#!/usr/bin/env python3
"""
Generate Figure 5: Shortlist Compression

Comparison diagram showing generic search vs mechanism-based matching.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

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

def generate_figure5(output_dir, format="both"):
    """Generate Figure 5: Shortlist Compression"""
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis('off')
    
    # Title
    ax.text(5, 5.5, 'Shortlist Compression: Generic vs Mechanism-Based Matching', 
            ha='center', va='top', fontsize=14, fontweight='bold')
    
    # Color scheme
    generic_color = '#FFCDD2'  # Light red
    mechanism_color = '#C8E6C9'  # Light green
    border_color = '#424242'  # Dark gray
    
    # === LEFT: Generic Search ===
    
    generic_box = FancyBboxPatch((0.5, 1.5), 3.5, 3.0,
                                boxstyle="round,pad=0.1",
                                facecolor=generic_color,
                                edgecolor=border_color,
                                linewidth=2.0)
    ax.add_patch(generic_box)
    
    ax.text(2.25, 4.0, 'Generic Search', ha='center', va='center',
            fontsize=12, fontweight='bold')
    ax.text(2.25, 3.5, '50+ Trials', ha='center', va='center',
            fontsize=11, fontweight='bold', color='#C62828')
    ax.text(2.25, 3.0, '• Trial 1', ha='center', va='center', fontsize=9)
    ax.text(2.25, 2.7, '• Trial 2', ha='center', va='center', fontsize=9)
    ax.text(2.25, 2.4, '• Trial 3', ha='center', va='center', fontsize=9)
    ax.text(2.25, 2.1, '• ...', ha='center', va='center', fontsize=9, style='italic')
    ax.text(2.25, 1.8, 'Time: 2-3 hours', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#C62828')
    
    # === RIGHT: Mechanism-Based Matching ===
    
    mechanism_box = FancyBboxPatch((6.0, 1.5), 3.5, 3.0,
                                  boxstyle="round,pad=0.1",
                                  facecolor=mechanism_color,
                                  edgecolor=border_color,
                                  linewidth=2.0)
    ax.add_patch(mechanism_box)
    
    ax.text(7.75, 4.0, 'Mechanism-Based\nMatching', ha='center', va='center',
            fontsize=12, fontweight='bold')
    ax.text(7.75, 3.3, '5-12 Trials', ha='center', va='center',
            fontsize=11, fontweight='bold', color='#2E7D32')
    ax.text(7.75, 2.9, '• Trial 1', ha='center', va='center', fontsize=9)
    ax.text(7.75, 2.7, '  (DDR: 0.989)', ha='center', va='center', 
            fontsize=8, style='italic', color='#2E7D32')
    ax.text(7.75, 2.4, '• Trial 2', ha='center', va='center', fontsize=9)
    ax.text(7.75, 2.2, '  (DDR: 0.989)', ha='center', va='center',
            fontsize=8, style='italic', color='#2E7D32')
    ax.text(7.75, 1.9, '• Trial 3', ha='center', va='center', fontsize=9)
    ax.text(7.75, 1.7, '  (DDR: 0.989)', ha='center', va='center',
            fontsize=8, style='italic', color='#2E7D32')
    ax.text(7.75, 1.5, 'Time: 30-45 min', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#2E7D32')
    
    # === ARROW ===
    
    arrow = FancyArrowPatch((4.0, 3.0), (6.0, 3.0),
                           arrowstyle='->', lw=3, color='#F57C00',
                           mutation_scale=30)
    ax.add_patch(arrow)
    
    # Compression annotation
    ax.text(5.0, 3.5, '60-65% Reduction', ha='center', va='center',
            fontsize=12, fontweight='bold', color='#F57C00',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', 
                     alpha=0.9, edgecolor='black', linewidth=2))
    
    # === STATISTICS ===
    
    stats_text = (
        "Compression Ratio: 10-24%\n"
        "Time Reduction: 60-65%\n"
        "Trial Count: 50+ → 5-12"
    )
    ax.text(5.0, 0.8, stats_text, ha='center', va='center',
            fontsize=10, style='italic',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                     alpha=0.8, edgecolor='gray', linewidth=1))
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if format in ["pdf", "both"]:
        plt.savefig(output_dir / "figure5_shortlist_compression.pdf", format='pdf', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure5_shortlist_compression.pdf'}")
    
    if format in ["png", "both"]:
        plt.savefig(output_dir / "figure5_shortlist_compression.png", format='png', bbox_inches='tight')
        print(f"  Saved: {output_dir / 'figure5_shortlist_compression.png'}")
    
    plt.close()

if __name__ == "__main__":
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "../figures"
    format_type = sys.argv[2] if len(sys.argv) > 2 else "both"
    generate_figure5(output_dir, format_type)


