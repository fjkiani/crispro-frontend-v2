#!/usr/bin/env python3
"""
Generate Figure 1: Metastasis Interception Framework Overview
4-panel figure for Nature Biotechnology submission
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Set up the figure
fig = plt.figure(figsize=(16, 12), facecolor='white')

# Color scheme
COLORS = {
    'primary': '#1a5276',      # Deep blue
    'secondary': '#2e86ab',    # Medium blue
    'accent': '#d35400',       # Orange
    'light': '#ebf5fb',        # Light blue background
    'evo2': '#8e44ad',         # Purple for Evo2
    'enformer': '#27ae60',     # Green for Enformer
    'alphafold': '#e74c3c',    # Red for AlphaFold
    'text': '#2c3e50',         # Dark text
}

# ============================================
# Panel A: 8-Step Metastatic Cascade (top-left)
# ============================================
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.axis('off')
ax1.set_title('A. Metastatic Cascade (8 Steps)', fontsize=14, fontweight='bold', loc='left')

steps = [
    'EMT', 'Local\nInvasion', 'Intravasation', 'Circulation',
    'Extravasation', 'Micrometastasis', 'Dormancy', 'Colonization'
]

# Draw cascade as connected boxes
for i, step in enumerate(steps):
    x = 1 + (i % 4) * 2.2
    y = 7 - (i // 4) * 3.5
    
    box = FancyBboxPatch((x-0.9, y-0.6), 1.8, 1.2,
                         boxstyle="round,pad=0.05,rounding_size=0.2",
                         facecolor=COLORS['light'],
                         edgecolor=COLORS['primary'],
                         linewidth=2)
    ax1.add_patch(box)
    ax1.text(x, y, step, ha='center', va='center', fontsize=9, fontweight='bold',
             color=COLORS['text'])
    ax1.text(x, y-0.4, f'Step {i+1}', ha='center', va='center', fontsize=7,
             color=COLORS['secondary'])

# Add arrows between steps
for i in range(7):
    x1 = 1.9 + (i % 4) * 2.2
    y1 = 7 - (i // 4) * 3.5
    
    if i == 3:  # Move to second row
        ax1.annotate('', xy=(1, 3.5), xytext=(8.1, 7),
                    arrowprops=dict(arrowstyle='->', color=COLORS['accent'], lw=2))
    elif i < 3:
        ax1.annotate('', xy=(x1+0.3, y1), xytext=(x1-0.1, y1),
                    arrowprops=dict(arrowstyle='->', color=COLORS['accent'], lw=2))
    elif i > 3:
        x1 = 1.9 + ((i-4) % 4) * 2.2
        ax1.annotate('', xy=(x1+0.3, 3.5), xytext=(x1-0.1, 3.5),
                    arrowprops=dict(arrowstyle='->', color=COLORS['accent'], lw=2))

# ============================================
# Panel B: Multi-Modal AI Integration (top-right)
# ============================================
ax2 = fig.add_subplot(2, 2, 2)
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis('off')
ax2.set_title('B. Multi-Modal AI Integration', fontsize=14, fontweight='bold', loc='left')

# Three input models
models = [
    ('Evo2', 'Sequence\nDisruption', COLORS['evo2'], 1.5),
    ('Enformer', 'Chromatin\nAccessibility', COLORS['enformer'], 5),
    ('AlphaFold 3', 'Structural\nValidation', COLORS['alphafold'], 8.5),
]

for name, desc, color, x in models:
    # Model box
    box = FancyBboxPatch((x-1.2, 6.5), 2.4, 2.5,
                         boxstyle="round,pad=0.05,rounding_size=0.3",
                         facecolor=color, edgecolor='white',
                         linewidth=2, alpha=0.9)
    ax2.add_patch(box)
    ax2.text(x, 8.2, name, ha='center', va='center', fontsize=11, fontweight='bold',
             color='white')
    ax2.text(x, 7.2, desc, ha='center', va='center', fontsize=8,
             color='white')
    
    # Arrow down
    ax2.annotate('', xy=(x, 5.5), xytext=(x, 6.4),
                arrowprops=dict(arrowstyle='->', color=color, lw=2))

# Integration box
int_box = FancyBboxPatch((2, 3.5), 6, 1.8,
                         boxstyle="round,pad=0.05,rounding_size=0.3",
                         facecolor=COLORS['primary'],
                         edgecolor='white', linewidth=2)
ax2.add_patch(int_box)
ax2.text(5, 4.4, 'Target Lock Score', ha='center', va='center', 
         fontsize=12, fontweight='bold', color='white')

# Output arrow
ax2.annotate('', xy=(5, 1.8), xytext=(5, 3.4),
            arrowprops=dict(arrowstyle='->', color=COLORS['primary'], lw=3))

# Output box
out_box = FancyBboxPatch((2.5, 0.5), 5, 1.2,
                         boxstyle="round,pad=0.05,rounding_size=0.2",
                         facecolor=COLORS['light'],
                         edgecolor=COLORS['primary'], linewidth=2)
ax2.add_patch(out_box)
ax2.text(5, 1.1, 'Gene Prioritization\n(AUROC 0.976)', ha='center', va='center',
         fontsize=10, fontweight='bold', color=COLORS['text'])

# ============================================
# Panel C: Target Lock Formula (bottom-left)
# ============================================
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.axis('off')
ax3.set_title('C. Target Lock Score Computation', fontsize=14, fontweight='bold', loc='left')

# Formula
formula = r'$\mathrm{TL}(g,s) = w_F \cdot F_g + w_E \cdot E_{g,s} + w_C \cdot C_g + w_R \cdot R_g$'
ax3.text(5, 8, formula, ha='center', va='center', fontsize=14, 
         fontweight='bold', color=COLORS['text'])

# Components with weights
components = [
    ('F', 'Functionality\n(Evo2)', 0.25, COLORS['evo2']),
    ('E', 'Essentiality\n(Step-specific)', 0.35, COLORS['accent']),
    ('C', 'Chromatin\n(Enformer)', 0.20, COLORS['enformer']),
    ('R', 'Regulatory\n(Conservation)', 0.20, COLORS['primary']),
]

for i, (symbol, desc, weight, color) in enumerate(components):
    x = 1.25 + i * 2.2
    
    # Component circle
    circle = plt.Circle((x, 5), 0.8, facecolor=color, edgecolor='white', linewidth=2)
    ax3.add_patch(circle)
    ax3.text(x, 5, symbol, ha='center', va='center', fontsize=16, 
             fontweight='bold', color='white')
    
    # Weight
    ax3.text(x, 3.8, f'w={weight}', ha='center', va='center', fontsize=9,
             color=COLORS['text'])
    
    # Description
    ax3.text(x, 2.5, desc, ha='center', va='center', fontsize=8,
             color=COLORS['text'])

# Threshold note
ax3.text(5, 0.8, 'Threshold: TL >= 0.70 for guide design prioritization',
         ha='center', va='center', fontsize=10, style='italic',
         color=COLORS['secondary'])

# ============================================
# Panel D: Assassin Score Pipeline (bottom-right)
# ============================================
ax4 = fig.add_subplot(2, 2, 4)
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 10)
ax4.axis('off')
ax4.set_title('D. Guide Design & Structural Validation', fontsize=14, fontweight='bold', loc='left')

# Pipeline steps
pipeline = [
    ('Target Gene', 'TL >= 0.70', 1.5, 8),
    ('Guide Design', '3 guides/gene', 5, 8),
    ('Efficacy\nPrediction', 'DeepCRISPR', 8.5, 8),
    ('Safety\nScoring', 'Off-target\nanalysis', 1.5, 4.5),
    ('AlphaFold 3\nValidation', 'pLDDT>=50\niPTM>=0.30', 5, 4.5),
    ('Assassin\nScore', 'Final ranking', 8.5, 4.5),
]

for i, (name, desc, x, y) in enumerate(pipeline):
    color = COLORS['alphafold'] if 'AlphaFold' in name else COLORS['primary']
    if 'Assassin' in name:
        color = COLORS['accent']
    
    box = FancyBboxPatch((x-1.2, y-1), 2.4, 2,
                         boxstyle="round,pad=0.05,rounding_size=0.2",
                         facecolor=color if i == 5 else COLORS['light'],
                         edgecolor=color,
                         linewidth=2)
    ax4.add_patch(box)
    
    text_color = 'white' if i == 5 else COLORS['text']
    ax4.text(x, y+0.3, name, ha='center', va='center', fontsize=9, 
             fontweight='bold', color=text_color)
    ax4.text(x, y-0.4, desc, ha='center', va='center', fontsize=7,
             color=text_color if i == 5 else COLORS['secondary'])

# Draw horizontal flow
ax4.annotate('', xy=(3.8, 8), xytext=(2.7, 8),
            arrowprops=dict(arrowstyle='->', color=COLORS['secondary'], lw=2))
ax4.annotate('', xy=(7.3, 8), xytext=(6.2, 8),
            arrowprops=dict(arrowstyle='->', color=COLORS['secondary'], lw=2))

# Down to second row
ax4.annotate('', xy=(9.5, 5.5), xytext=(9.5, 7),
            arrowprops=dict(arrowstyle='->', color=COLORS['secondary'], lw=2))

# Second row flow
ax4.annotate('', xy=(7.3, 4.5), xytext=(6.2, 4.5),
            arrowprops=dict(arrowstyle='<-', color=COLORS['secondary'], lw=2))
ax4.annotate('', xy=(3.8, 4.5), xytext=(2.7, 4.5),
            arrowprops=dict(arrowstyle='<-', color=COLORS['secondary'], lw=2))

# Validation result box
result_box = FancyBboxPatch((3, 1), 4, 1.2,
                            boxstyle="round,pad=0.05,rounding_size=0.2",
                            facecolor=COLORS['light'],
                            edgecolor=COLORS['accent'], linewidth=2)
ax4.add_patch(result_box)
ax4.text(5, 1.6, '15/15 guides validated', ha='center', va='center',
         fontsize=10, fontweight='bold', color=COLORS['text'])
ax4.text(5, 1.2, '100% structural pass rate', ha='center', va='center',
         fontsize=8, color=COLORS['secondary'])

ax4.annotate('', xy=(5, 2.2), xytext=(5, 3.4),
            arrowprops=dict(arrowstyle='<-', color=COLORS['accent'], lw=2))

# Main title
fig.suptitle('Metastasis Interception: Multi-Modal CRISPR Guide Design Framework',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save
plt.savefig('/Users/fahadkiani/Desktop/development/crispr-assistant-main/publication/figures/Kiani_Figure1.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('/Users/fahadkiani/Desktop/development/crispr-assistant-main/publication/figures/Kiani_Figure1.svg',
            bbox_inches='tight', facecolor='white')

print("Figure 1 generated successfully!")
print("   - Kiani_Figure1.png (300 DPI)")
print("   - Kiani_Figure1.svg (vector)")
