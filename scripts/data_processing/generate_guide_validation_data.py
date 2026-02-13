#!/usr/bin/env python3
"""
Generate 80-Guide Validation Dataset for Publication
Strategy: 10 guides per mission step (8 steps × 10 = 80 guides)
Author: Zo (Agent)
Date: October 7, 2025
"""

import asyncio
import httpx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
from typing import List, Dict

# Configuration
API_BASE = "http://127.0.0.1:8000"
TIMEOUT = 180.0

MISSION_STEPS = [
    "primary_growth",
    "local_invasion",
    "intravasation",
    "survival_in_circulation",
    "extravasation",
    "micrometastasis_formation",
    "angiogenesis",
    "metastatic_colonization"
]

# One representative gene per mission step
MISSION_TO_GENE = {
    "primary_growth": ("BRAF", "7", 140753336, "A", "T"),  # V600E
    "local_invasion": ("TWIST1", "7", 19069313, "G", "A"),
    "intravasation": ("MMP2", "16", 55448195, "C", "T"),
    "survival_in_circulation": ("BCL2", "18", 63320128, "G", "A"),  # Placeholder
    "extravasation": ("VCAM1", "1", 100720221, "C", "T"),  # Placeholder
    "micrometastasis_formation": ("CXCR4", "2", 136116763, "A", "G"),
    "angiogenesis": ("VEGFA", "6", 43778335, "G", "A"),
    "metastatic_colonization": ("MET", "7", 116735218, "T", "C")
}


async def generate_guides_for_mission(mission: str, num_guides: int = 10) -> List[Dict]:
    """Generate guide RNAs for a mission step"""
    gene, chrom, pos, ref, alt = MISSION_TO_GENE[mission]
    
    print(f"🎯 {mission} → {gene}")
    
    payload = {
        "mutations": [{
            "gene": gene,
            "hgvs_p": "canonical_target",
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt
        }],
        "mission_step": mission,
        "patient_id": f"GUIDE_VAL_{mission}",
        "options": {
            "model_id": "evo2_1b",
            "profile": "baseline",
            "num_candidates": num_guides
        }
    }
    
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            response = await client.post(
                f"{API_BASE}/api/metastasis/intercept",
                json=payload
            )
            
            if response.status_code == 200:
                data = response.json()
                candidates = data.get("candidates", [])
                
                # Add mission and gene metadata
                for candidate in candidates:
                    candidate["mission_step"] = mission
                    candidate["target_gene"] = gene
                
                print(f"   ✅ Generated {len(candidates)} guides")
                return candidates
            else:
                print(f"   ⚠️  Error {response.status_code}: {response.text[:200]}")
                return []
                
    except Exception as e:
        print(f"   ⚠️  Exception: {e}")
        return []


async def generate_all_guides() -> pd.DataFrame:
    """Generate guides for all 8 mission steps"""
    print("⚔️  Generating 80-Guide Validation Dataset...")
    print(f"📊 Target: 10 guides × 8 mission steps = 80 guides total")
    print()
    
    all_guides = []
    
    for mission in MISSION_STEPS:
        guides = await generate_guides_for_mission(mission, num_guides=10)
        all_guides.extend(guides)
    
    print()
    print(f"✅ Total guides generated: {len(all_guides)}")
    
    df = pd.DataFrame(all_guides)
    
    # Save data
    output_csv = Path("data/guide_validation_dataset.csv")
    output_json = Path("data/guide_validation_dataset.json")
    
    df.to_csv(output_csv, index=False)
    
    with open(output_json, 'w') as f:
        json.dump(all_guides, f, indent=2)
    
    print(f"💾 Saved data to:")
    print(f"   - {output_csv}")
    print(f"   - {output_json}")
    print()
    
    return df


def generate_figures(df: pd.DataFrame):
    """Generate F3-F5 figures for publication"""
    print("📊 Generating publication figures...")
    
    if df.empty:
        print("⚠️  No data to plot!")
        return
    
    # Figure 3: Efficacy Distribution
    print("   Creating F3: Efficacy Distribution...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if 'efficacy_proxy' in df.columns:
        sns.violinplot(data=df, x='mission_step', y='efficacy_proxy', ax=ax)
        ax.set_title('Guide RNA Efficacy Distribution by Mission Step', fontsize=14, weight='bold')
        ax.set_xlabel('Mission Step', fontsize=12)
        ax.set_ylabel('Efficacy Score', fontsize=12)
        ax.set_ylim(0, 1)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('figures/F3_efficacy_distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig('figures/F3_efficacy_distribution.svg', bbox_inches='tight')
        plt.close()
        print("   ✅ F3 saved")
    
    # Figure 4: Safety Distribution
    print("   Creating F4: Safety Distribution...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if 'safety_score' in df.columns:
        sns.violinplot(data=df, x='mission_step', y='safety_score', ax=ax)
        ax.set_title('Guide RNA Safety Distribution by Mission Step', fontsize=14, weight='bold')
        ax.set_xlabel('Mission Step', fontsize=12)
        ax.set_ylabel('Safety Score', fontsize=12)
        ax.set_ylim(0, 1)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('figures/F4_safety_distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig('figures/F4_safety_distribution.svg', bbox_inches='tight')
        plt.close()
        print("   ✅ F4 saved")
    
    # Figure 5: Assassin Score Distribution
    print("   Creating F5: Assassin Score Distribution...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if 'assassin_score' in df.columns:
        sns.boxplot(data=df, x='mission_step', y='assassin_score', ax=ax)
        ax.set_title('Assassin Score Distribution by Mission Step', fontsize=14, weight='bold')
        ax.set_xlabel('Mission Step', fontsize=12)
        ax.set_ylabel('Assassin Score', fontsize=12)
        ax.set_ylim(0, 1)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('figures/F5_assassin_score_distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig('figures/F5_assassin_score_distribution.svg', bbox_inches='tight')
        plt.close()
        print("   ✅ F5 saved")
    
    print("✅ Figures complete")
    print()


def generate_table2(df: pd.DataFrame):
    """Generate Table 2: Performance Metrics"""
    print("📊 Generating Table 2: Performance Metrics...")
    
    if df.empty:
        print("⚠️  No data for table!")
        return
    
    metrics = {}
    
    for col in ['efficacy_proxy', 'safety_score', 'assassin_score']:
        if col in df.columns:
            metrics[col] = {
                'mean': df[col].mean(),
                'std': df[col].std(),
                'min': df[col].min(),
                'max': df[col].max(),
                'median': df[col].median()
            }
    
    # Create table
    table_data = []
    for metric_name, stats in metrics.items():
        table_data.append({
            'Metric': metric_name.replace('_', ' ').title(),
            'Mean ± SD': f"{stats['mean']:.3f} ± {stats['std']:.3f}",
            'Min': f"{stats['min']:.3f}",
            'Max': f"{stats['max']:.3f}",
            'Median': f"{stats['median']:.3f}"
        })
    
    table_df = pd.DataFrame(table_data)
    
    # Save as CSV
    table_df.to_csv('data/table2_performance_metrics.csv', index=False)
    
    # Save as LaTeX
    latex = table_df.to_latex(index=False)
    with open('data/table2_performance_metrics.tex', 'w') as f:
        f.write(latex)
    
    print("✅ Table 2 saved:")
    print("   - data/table2_performance_metrics.csv")
    print("   - data/table2_performance_metrics.tex")
    print()
    print("📈 Summary Statistics:")
    print(table_df.to_string(index=False))
    print()


async def main():
    """Main execution"""
    print("=" * 70)
    print("⚔️  80-GUIDE VALIDATION DATASET GENERATION")
    print("=" * 70)
    print()
    
    # Generate guides
    df = await generate_all_guides()
    
    if not df.empty:
        # Generate figures
        generate_figures(df)
        
        # Generate table
        generate_table2(df)
    
    print("=" * 70)
    print("✅ GUIDE VALIDATION DATASET COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    asyncio.run(main())

