#!/usr/bin/env python3
"""
Generate Target Lock Heatmap Data for Publication Figure 2
Author: Zo (Agent)
Date: October 7, 2025
"""

import asyncio
import httpx
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import json
from typing import List, Dict

# Configuration
API_BASE = "http://127.0.0.1:8000"
TIMEOUT = 120.0

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

TARGET_GENES = ["VEGFA", "MMP2", "TWIST1", "SNAIL1", "CXCR4", "MET", "BRAF"]

# Gene coordinates (GRCh38) - using known hotspots where available
GENE_COORDS = {
    # VEGFA - canonical exon 3 position (angiogenesis target)
    "VEGFA": {"chrom": "6", "pos": 43778335, "ref": "G", "alt": "A"},
    
    # MMP2 - matrix metalloproteinase (invasion/intravasation)
    "MMP2": {"chrom": "16", "pos": 55448195, "ref": "C", "alt": "T"},
    
    # TWIST1 - EMT transcription factor (local invasion)
    "TWIST1": {"chrom": "7", "pos": 19069313, "ref": "G", "alt": "A"},
    
    # SNAI1 (SNAIL1) - EMT transcription factor
    "SNAIL1": {"chrom": "20", "pos": 49985933, "ref": "C", "alt": "T"},
    
    # CXCR4 - chemokine receptor (homing/micrometastasis)
    "CXCR4": {"chrom": "2", "pos": 136116763, "ref": "A", "alt": "G"},
    
    # MET - receptor tyrosine kinase (colonization)
    "MET": {"chrom": "7", "pos": 116735218, "ref": "T", "alt": "C"},
    
    # BRAF V600E - known hotspot (primary growth)
    "BRAF": {"chrom": "7", "pos": 140753336, "ref": "A", "alt": "T"}  # V600E
}


async def fetch_target_lock_score(gene: str, mission: str, session_id: int) -> Dict:
    """Fetch target lock score for gene in mission context"""
    coords = GENE_COORDS[gene]
    
    payload = {
        "mutations": [{
            "gene": gene,
            "hgvs_p": "synthetic_variant",
            "chrom": coords["chrom"],
            "pos": coords["pos"],
            "ref": coords["ref"],
            "alt": coords["alt"]
        }],
        "mission_step": mission,
        "patient_id": f"HEATMAP_{session_id}_{gene}_{mission}",
        "options": {"model_id": "evo2_1b", "profile": "baseline"}
    }
    
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            response = await client.post(
                f"{API_BASE}/api/metastasis/intercept",
                json=payload
            )
            response.raise_for_status()
            data = response.json()
            
            target = data.get("validated_target", {})
            
            return {
                "gene": gene,
                "mission": mission,
                "target_lock_score": target.get("target_lock_score", 0.0),
                "functionality_score": target.get("functionality_score", 0.0),
                "essentiality_score": target.get("essentiality_score", 0.0),
                "chromatin_score": target.get("chromatin_score", 0.0),
                "regulatory_score": target.get("regulatory_score", 0.0),
                "rationale": target.get("rationale", []),
                "provenance": data.get("provenance", {})
            }
    except Exception as e:
        print(f"⚠️  Error fetching {gene} × {mission}: {e}")
        return {
            "gene": gene,
            "mission": mission,
            "target_lock_score": 0.0,
            "functionality_score": 0.0,
            "essentiality_score": 0.0,
            "chromatin_score": 0.0,
            "regulatory_score": 0.0,
            "rationale": [f"Error: {str(e)}"],
            "provenance": {}
        }


async def generate_heatmap_data() -> pd.DataFrame:
    """Generate full heatmap dataset (8 steps × 7 genes = 56 analyses)"""
    print("⚔️  Generating Target Lock Heatmap Data...")
    print(f"📊 Total analyses: {len(MISSION_STEPS)} steps × {len(TARGET_GENES)} genes = {len(MISSION_STEPS) * len(TARGET_GENES)}")
    print()
    
    tasks = []
    session_id = 0
    
    for mission in MISSION_STEPS:
        for gene in TARGET_GENES:
            tasks.append(fetch_target_lock_score(gene, mission, session_id))
            session_id += 1
    
    print(f"🚀 Launching {len(tasks)} parallel API calls...")
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    # Filter out exceptions
    valid_results = [r for r in results if isinstance(r, dict)]
    failed_count = len(results) - len(valid_results)
    
    print(f"✅ Completed {len(valid_results)}/{len(results)} analyses")
    if failed_count > 0:
        print(f"⚠️  {failed_count} analyses failed")
    
    df = pd.DataFrame(valid_results)
    
    # Save raw data
    output_csv = Path("data/target_lock_heatmap_data.csv")
    output_json = Path("data/target_lock_heatmap_data.json")
    
    df.to_csv(output_csv, index=False)
    
    # Save full data including rationale and provenance
    with open(output_json, 'w') as f:
        json.dump(valid_results, f, indent=2)
    
    print(f"💾 Saved data to:")
    print(f"   - {output_csv}")
    print(f"   - {output_json}")
    print()
    
    return df


def plot_heatmap(df: pd.DataFrame):
    """Generate target lock heatmap visualization (Figure 2)"""
    print("📊 Generating heatmap visualization...")
    
    # Pivot data for heatmap
    pivot = df.pivot(index="gene", columns="mission", values="target_lock_score")
    
    # Create figure
    plt.figure(figsize=(16, 8))
    
    # Generate heatmap
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={'label': 'Target Lock Score', 'shrink': 0.8},
        linewidths=0.5,
        linecolor='gray'
    )
    
    # Styling
    plt.title(
        "Target Lock Scores: Metastatic Mission Steps × Target Genes",
        fontsize=18,
        weight='bold',
        pad=20
    )
    plt.xlabel("Metastatic Mission Step", fontsize=14, weight='bold')
    plt.ylabel("Target Gene", fontsize=14, weight='bold')
    
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha='right', fontsize=11)
    plt.yticks(fontsize=12, weight='bold')
    
    plt.tight_layout()
    
    # Save outputs
    output_png = Path("figures/F2_target_lock_heatmap.png")
    output_svg = Path("figures/F2_target_lock_heatmap.svg")
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_svg, format='svg', bbox_inches='tight')
    
    print(f"✅ Heatmap saved to:")
    print(f"   - {output_png}")
    print(f"   - {output_svg}")
    print()
    
    # Generate summary statistics
    print("📈 Summary Statistics:")
    print(f"   Mean target lock score: {df['target_lock_score'].mean():.3f}")
    print(f"   Median target lock score: {df['target_lock_score'].median():.3f}")
    print(f"   Min target lock score: {df['target_lock_score'].min():.3f}")
    print(f"   Max target lock score: {df['target_lock_score'].max():.3f}")
    print(f"   Std dev: {df['target_lock_score'].std():.3f}")
    print()
    
    # Top 5 gene-mission pairs
    print("🎯 Top 5 Gene-Mission Pairs (Highest Target Lock Score):")
    top5 = df.nlargest(5, 'target_lock_score')[['gene', 'mission', 'target_lock_score']]
    for idx, row in top5.iterrows():
        print(f"   {row['gene']} × {row['mission']}: {row['target_lock_score']:.3f}")
    print()


def plot_component_scores(df: pd.DataFrame):
    """Generate component score breakdown visualization (Supplementary)"""
    print("📊 Generating component score breakdown...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    score_cols = [
        ('functionality_score', 'Functionality Score'),
        ('essentiality_score', 'Essentiality Score'),
        ('chromatin_score', 'Chromatin Accessibility Score'),
        ('regulatory_score', 'Regulatory Impact Score')
    ]
    
    for idx, (col, title) in enumerate(score_cols):
        ax = axes[idx // 2, idx % 2]
        pivot = df.pivot(index="gene", columns="mission", values=col)
        
        sns.heatmap(
            pivot,
            annot=True,
            fmt=".3f",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
            cbar_kws={'label': title, 'shrink': 0.8},
            linewidths=0.5,
            ax=ax
        )
        
        ax.set_title(title, fontsize=14, weight='bold')
        ax.set_xlabel("Mission Step", fontsize=11)
        ax.set_ylabel("Gene", fontsize=11)
        ax.tick_params(axis='x', rotation=45, labelsize=9)
        ax.tick_params(axis='y', labelsize=10)
    
    plt.tight_layout()
    
    output_png = Path("figures/F2_supp_component_scores.png")
    output_svg = Path("figures/F2_supp_component_scores.svg")
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_svg, format='svg', bbox_inches='tight')
    
    print(f"✅ Component breakdown saved to:")
    print(f"   - {output_png}")
    print(f"   - {output_svg}")
    print()


async def main():
    """Main execution"""
    print("=" * 60)
    print("⚔️  TARGET LOCK HEATMAP DATA GENERATION")
    print("=" * 60)
    print()
    
    # Generate data
    df = await generate_heatmap_data()
    
    # Generate visualizations
    plot_heatmap(df)
    plot_component_scores(df)
    
    print("=" * 60)
    print("✅ FIGURE 2 GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())

    # Generate visualizations
    plot_heatmap(df)
    plot_component_scores(df)
    
    print("=" * 60)
    print("✅ FIGURE 2 GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())


    # Generate visualizations
    plot_heatmap(df)
    plot_component_scores(df)
    
    print("=" * 60)
    print("✅ FIGURE 2 GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())


    # Generate visualizations
    plot_heatmap(df)
    plot_component_scores(df)
    
    print("=" * 60)
    print("✅ FIGURE 2 GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())


    # Generate visualizations
    plot_heatmap(df)
    plot_component_scores(df)
    
    print("=" * 60)
    print("✅ FIGURE 2 GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())


    # Generate visualizations
    plot_heatmap(df)
    plot_component_scores(df)
    
    print("=" * 60)
    print("✅ FIGURE 2 GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())


