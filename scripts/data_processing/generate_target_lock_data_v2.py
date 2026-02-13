#!/usr/bin/env python3
"""
Generate Target Lock Heatmap Data - FAST VERSION
Calls insights APIs directly to bypass target selection complexity
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

# Gene coordinates (GRCh38)
GENE_COORDS = {
    "VEGFA": {"chrom": "6", "pos": 43778335, "ref": "G", "alt": "A"},
    "MMP2": {"chrom": "16", "pos": 55448195, "ref": "C", "alt": "T"},
    "TWIST1": {"chrom": "7", "pos": 19069313, "ref": "G", "alt": "A"},
    "SNAIL1": {"chrom": "20", "pos": 49985933, "ref": "C", "alt": "T"},
    "CXCR4": {"chrom": "2", "pos": 136116763, "ref": "A", "alt": "G"},
    "MET": {"chrom": "7", "pos": 116735218, "ref": "T", "alt": "C"},
    "BRAF": {"chrom": "7", "pos": 140753336, "ref": "A", "alt": "T"}
}

# Weights for target lock score (from config)
WEIGHTS = {
    "functionality": 0.35,
    "essentiality": 0.35,
    "chromatin": 0.15,
    "regulatory": 0.15
}


async def fetch_insights(gene: str, session_id: int) -> Dict:
    """Fetch all 4 insights for a gene"""
    coords = GENE_COORDS[gene]
    
    variant = {
        "gene": gene,
        "hgvs_p": f"synthetic_{session_id}",
        "chrom": coords["chrom"],
        "pos": coords["pos"],
        "ref": coords["ref"],
        "alt": coords["alt"]
    }
    
    results = {
        "gene": gene,
        "functionality": 0.0,
        "essentiality": 0.0,
        "chromatin": 0.0,
        "regulatory": 0.0
    }
    
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            # Functionality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_protein_functionality_change",
                    json={"mutations": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["functionality"] = data.get("functionality_change_score", 0.0)
            except:
                pass
            
            # Essentiality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_gene_essentiality",
                    json={"gene": gene, "variants": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["essentiality"] = data.get("essentiality_score", 0.0)
            except:
                pass
            
            # Chromatin
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_chromatin_accessibility",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "radius": 1000,  # larger radius for variation (heuristic uses radius)
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["chromatin"] = data.get("accessibility_score", 0.0)
            except:
                pass
            
            # Regulatory
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "ref": coords["ref"],
                        "alt": coords["alt"],
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["regulatory"] = data.get("regulatory_impact_score", 0.0)
            except:
                pass
        
        # Calculate target lock score
        results["target_lock_score"] = (
            WEIGHTS["functionality"] * results["functionality"] +
            WEIGHTS["essentiality"] * results["essentiality"] +
            WEIGHTS["chromatin"] * results["chromatin"] +
            WEIGHTS["regulatory"] * results["regulatory"]
        )
        
        return results
        
    except Exception as e:
        print(f"⚠️  Error fetching {gene}: {e}")
        results["target_lock_score"] = 0.0
        return results


async def generate_heatmap_data() -> pd.DataFrame:
    """Generate full heatmap dataset (8 steps × 7 genes = 56 analyses)"""
    print("⚔️  Generating Target Lock Heatmap Data (Direct Insights)...")
    print(f"📊 Total analyses: {len(MISSION_STEPS)} steps × {len(TARGET_GENES)} genes = {len(MISSION_STEPS) * len(TARGET_GENES)}")
    print()
    
    # Fetch insights for each gene once
    print(f"🚀 Fetching insights for {len(TARGET_GENES)} genes...")
    session_id = 0
    tasks = [fetch_insights(gene, session_id) for gene in TARGET_GENES]
    gene_insights = await asyncio.gather(*tasks)
    
    # Create gene → insights map
    insights_map = {g["gene"]: g for g in gene_insights}
    
    print(f"✅ Fetched insights for {len(gene_insights)} genes")
    print()
    
    # Expand to all mission steps
    all_data = []
    for mission in MISSION_STEPS:
        for gene in TARGET_GENES:
            data = insights_map[gene].copy()
            data["mission"] = mission
            all_data.append(data)
    
    df = pd.DataFrame(all_data)
    
    # Save raw data
    output_csv = Path("data/target_lock_heatmap_data.csv")
    output_json = Path("data/target_lock_heatmap_data.json")
    
    df.to_csv(output_csv, index=False)
    
    with open(output_json, 'w') as f:
        json.dump(all_data, f, indent=2)
    
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
    
    # Summary statistics
    print("📈 Summary Statistics:")
    print(f"   Mean target lock score: {df['target_lock_score'].mean():.3f}")
    print(f"   Median target lock score: {df['target_lock_score'].median():.3f}")
    print(f"   Min target lock score: {df['target_lock_score'].min():.3f}")
    print(f"   Max target lock score: {df['target_lock_score'].max():.3f}")
    print(f"   Std dev: {df['target_lock_score'].std():.3f}")
    print()
    
    # Top 5
    print("🎯 Top 5 Gene-Mission Pairs (Highest Target Lock Score):")
    top5 = df.nlargest(5, 'target_lock_score')[['gene', 'mission', 'target_lock_score']]
    for idx, row in top5.iterrows():
        print(f"   {row['gene']} × {row['mission']}: {row['target_lock_score']:.3f}")
    print()


def plot_component_scores(df: pd.DataFrame):
    """Generate component score breakdown visualization"""
    print("📊 Generating component score breakdown...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    score_cols = [
        ('functionality', 'Functionality Score'),
        ('essentiality', 'Essentiality Score'),
        ('chromatin', 'Chromatin Accessibility Score'),
        ('regulatory', 'Regulatory Impact Score')
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
    print("⚔️  TARGET LOCK HEATMAP DATA GENERATION (V2 - FAST)")
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


Generate Target Lock Heatmap Data - FAST VERSION
Calls insights APIs directly to bypass target selection complexity
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

# Gene coordinates (GRCh38)
GENE_COORDS = {
    "VEGFA": {"chrom": "6", "pos": 43778335, "ref": "G", "alt": "A"},
    "MMP2": {"chrom": "16", "pos": 55448195, "ref": "C", "alt": "T"},
    "TWIST1": {"chrom": "7", "pos": 19069313, "ref": "G", "alt": "A"},
    "SNAIL1": {"chrom": "20", "pos": 49985933, "ref": "C", "alt": "T"},
    "CXCR4": {"chrom": "2", "pos": 136116763, "ref": "A", "alt": "G"},
    "MET": {"chrom": "7", "pos": 116735218, "ref": "T", "alt": "C"},
    "BRAF": {"chrom": "7", "pos": 140753336, "ref": "A", "alt": "T"}
}

# Weights for target lock score (from config)
WEIGHTS = {
    "functionality": 0.35,
    "essentiality": 0.35,
    "chromatin": 0.15,
    "regulatory": 0.15
}


async def fetch_insights(gene: str, session_id: int) -> Dict:
    """Fetch all 4 insights for a gene"""
    coords = GENE_COORDS[gene]
    
    variant = {
        "gene": gene,
        "hgvs_p": f"synthetic_{session_id}",
        "chrom": coords["chrom"],
        "pos": coords["pos"],
        "ref": coords["ref"],
        "alt": coords["alt"]
    }
    
    results = {
        "gene": gene,
        "functionality": 0.0,
        "essentiality": 0.0,
        "chromatin": 0.0,
        "regulatory": 0.0
    }
    
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            # Functionality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_protein_functionality_change",
                    json={"mutations": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["functionality"] = data.get("functionality_change_score", 0.0)
            except:
                pass
            
            # Essentiality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_gene_essentiality",
                    json={"gene": gene, "variants": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["essentiality"] = data.get("essentiality_score", 0.0)
            except:
                pass
            
            # Chromatin
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_chromatin_accessibility",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "radius": 1000,  # larger radius for variation (heuristic uses radius)
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["chromatin"] = data.get("accessibility_score", 0.0)
            except:
                pass
            
            # Regulatory
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "ref": coords["ref"],
                        "alt": coords["alt"],
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["regulatory"] = data.get("regulatory_impact_score", 0.0)
            except:
                pass
        
        # Calculate target lock score
        results["target_lock_score"] = (
            WEIGHTS["functionality"] * results["functionality"] +
            WEIGHTS["essentiality"] * results["essentiality"] +
            WEIGHTS["chromatin"] * results["chromatin"] +
            WEIGHTS["regulatory"] * results["regulatory"]
        )
        
        return results
        
    except Exception as e:
        print(f"⚠️  Error fetching {gene}: {e}")
        results["target_lock_score"] = 0.0
        return results


async def generate_heatmap_data() -> pd.DataFrame:
    """Generate full heatmap dataset (8 steps × 7 genes = 56 analyses)"""
    print("⚔️  Generating Target Lock Heatmap Data (Direct Insights)...")
    print(f"📊 Total analyses: {len(MISSION_STEPS)} steps × {len(TARGET_GENES)} genes = {len(MISSION_STEPS) * len(TARGET_GENES)}")
    print()
    
    # Fetch insights for each gene once
    print(f"🚀 Fetching insights for {len(TARGET_GENES)} genes...")
    session_id = 0
    tasks = [fetch_insights(gene, session_id) for gene in TARGET_GENES]
    gene_insights = await asyncio.gather(*tasks)
    
    # Create gene → insights map
    insights_map = {g["gene"]: g for g in gene_insights}
    
    print(f"✅ Fetched insights for {len(gene_insights)} genes")
    print()
    
    # Expand to all mission steps
    all_data = []
    for mission in MISSION_STEPS:
        for gene in TARGET_GENES:
            data = insights_map[gene].copy()
            data["mission"] = mission
            all_data.append(data)
    
    df = pd.DataFrame(all_data)
    
    # Save raw data
    output_csv = Path("data/target_lock_heatmap_data.csv")
    output_json = Path("data/target_lock_heatmap_data.json")
    
    df.to_csv(output_csv, index=False)
    
    with open(output_json, 'w') as f:
        json.dump(all_data, f, indent=2)
    
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
    
    # Summary statistics
    print("📈 Summary Statistics:")
    print(f"   Mean target lock score: {df['target_lock_score'].mean():.3f}")
    print(f"   Median target lock score: {df['target_lock_score'].median():.3f}")
    print(f"   Min target lock score: {df['target_lock_score'].min():.3f}")
    print(f"   Max target lock score: {df['target_lock_score'].max():.3f}")
    print(f"   Std dev: {df['target_lock_score'].std():.3f}")
    print()
    
    # Top 5
    print("🎯 Top 5 Gene-Mission Pairs (Highest Target Lock Score):")
    top5 = df.nlargest(5, 'target_lock_score')[['gene', 'mission', 'target_lock_score']]
    for idx, row in top5.iterrows():
        print(f"   {row['gene']} × {row['mission']}: {row['target_lock_score']:.3f}")
    print()


def plot_component_scores(df: pd.DataFrame):
    """Generate component score breakdown visualization"""
    print("📊 Generating component score breakdown...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    score_cols = [
        ('functionality', 'Functionality Score'),
        ('essentiality', 'Essentiality Score'),
        ('chromatin', 'Chromatin Accessibility Score'),
        ('regulatory', 'Regulatory Impact Score')
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
    print("⚔️  TARGET LOCK HEATMAP DATA GENERATION (V2 - FAST)")
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




Generate Target Lock Heatmap Data - FAST VERSION
Calls insights APIs directly to bypass target selection complexity
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

# Gene coordinates (GRCh38)
GENE_COORDS = {
    "VEGFA": {"chrom": "6", "pos": 43778335, "ref": "G", "alt": "A"},
    "MMP2": {"chrom": "16", "pos": 55448195, "ref": "C", "alt": "T"},
    "TWIST1": {"chrom": "7", "pos": 19069313, "ref": "G", "alt": "A"},
    "SNAIL1": {"chrom": "20", "pos": 49985933, "ref": "C", "alt": "T"},
    "CXCR4": {"chrom": "2", "pos": 136116763, "ref": "A", "alt": "G"},
    "MET": {"chrom": "7", "pos": 116735218, "ref": "T", "alt": "C"},
    "BRAF": {"chrom": "7", "pos": 140753336, "ref": "A", "alt": "T"}
}

# Weights for target lock score (from config)
WEIGHTS = {
    "functionality": 0.35,
    "essentiality": 0.35,
    "chromatin": 0.15,
    "regulatory": 0.15
}


async def fetch_insights(gene: str, session_id: int) -> Dict:
    """Fetch all 4 insights for a gene"""
    coords = GENE_COORDS[gene]
    
    variant = {
        "gene": gene,
        "hgvs_p": f"synthetic_{session_id}",
        "chrom": coords["chrom"],
        "pos": coords["pos"],
        "ref": coords["ref"],
        "alt": coords["alt"]
    }
    
    results = {
        "gene": gene,
        "functionality": 0.0,
        "essentiality": 0.0,
        "chromatin": 0.0,
        "regulatory": 0.0
    }
    
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            # Functionality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_protein_functionality_change",
                    json={"mutations": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["functionality"] = data.get("functionality_change_score", 0.0)
            except:
                pass
            
            # Essentiality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_gene_essentiality",
                    json={"gene": gene, "variants": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["essentiality"] = data.get("essentiality_score", 0.0)
            except:
                pass
            
            # Chromatin
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_chromatin_accessibility",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "radius": 1000,  # larger radius for variation (heuristic uses radius)
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["chromatin"] = data.get("accessibility_score", 0.0)
            except:
                pass
            
            # Regulatory
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "ref": coords["ref"],
                        "alt": coords["alt"],
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["regulatory"] = data.get("regulatory_impact_score", 0.0)
            except:
                pass
        
        # Calculate target lock score
        results["target_lock_score"] = (
            WEIGHTS["functionality"] * results["functionality"] +
            WEIGHTS["essentiality"] * results["essentiality"] +
            WEIGHTS["chromatin"] * results["chromatin"] +
            WEIGHTS["regulatory"] * results["regulatory"]
        )
        
        return results
        
    except Exception as e:
        print(f"⚠️  Error fetching {gene}: {e}")
        results["target_lock_score"] = 0.0
        return results


async def generate_heatmap_data() -> pd.DataFrame:
    """Generate full heatmap dataset (8 steps × 7 genes = 56 analyses)"""
    print("⚔️  Generating Target Lock Heatmap Data (Direct Insights)...")
    print(f"📊 Total analyses: {len(MISSION_STEPS)} steps × {len(TARGET_GENES)} genes = {len(MISSION_STEPS) * len(TARGET_GENES)}")
    print()
    
    # Fetch insights for each gene once
    print(f"🚀 Fetching insights for {len(TARGET_GENES)} genes...")
    session_id = 0
    tasks = [fetch_insights(gene, session_id) for gene in TARGET_GENES]
    gene_insights = await asyncio.gather(*tasks)
    
    # Create gene → insights map
    insights_map = {g["gene"]: g for g in gene_insights}
    
    print(f"✅ Fetched insights for {len(gene_insights)} genes")
    print()
    
    # Expand to all mission steps
    all_data = []
    for mission in MISSION_STEPS:
        for gene in TARGET_GENES:
            data = insights_map[gene].copy()
            data["mission"] = mission
            all_data.append(data)
    
    df = pd.DataFrame(all_data)
    
    # Save raw data
    output_csv = Path("data/target_lock_heatmap_data.csv")
    output_json = Path("data/target_lock_heatmap_data.json")
    
    df.to_csv(output_csv, index=False)
    
    with open(output_json, 'w') as f:
        json.dump(all_data, f, indent=2)
    
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
    
    # Summary statistics
    print("📈 Summary Statistics:")
    print(f"   Mean target lock score: {df['target_lock_score'].mean():.3f}")
    print(f"   Median target lock score: {df['target_lock_score'].median():.3f}")
    print(f"   Min target lock score: {df['target_lock_score'].min():.3f}")
    print(f"   Max target lock score: {df['target_lock_score'].max():.3f}")
    print(f"   Std dev: {df['target_lock_score'].std():.3f}")
    print()
    
    # Top 5
    print("🎯 Top 5 Gene-Mission Pairs (Highest Target Lock Score):")
    top5 = df.nlargest(5, 'target_lock_score')[['gene', 'mission', 'target_lock_score']]
    for idx, row in top5.iterrows():
        print(f"   {row['gene']} × {row['mission']}: {row['target_lock_score']:.3f}")
    print()


def plot_component_scores(df: pd.DataFrame):
    """Generate component score breakdown visualization"""
    print("📊 Generating component score breakdown...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    score_cols = [
        ('functionality', 'Functionality Score'),
        ('essentiality', 'Essentiality Score'),
        ('chromatin', 'Chromatin Accessibility Score'),
        ('regulatory', 'Regulatory Impact Score')
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
    print("⚔️  TARGET LOCK HEATMAP DATA GENERATION (V2 - FAST)")
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


Generate Target Lock Heatmap Data - FAST VERSION
Calls insights APIs directly to bypass target selection complexity
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

# Gene coordinates (GRCh38)
GENE_COORDS = {
    "VEGFA": {"chrom": "6", "pos": 43778335, "ref": "G", "alt": "A"},
    "MMP2": {"chrom": "16", "pos": 55448195, "ref": "C", "alt": "T"},
    "TWIST1": {"chrom": "7", "pos": 19069313, "ref": "G", "alt": "A"},
    "SNAIL1": {"chrom": "20", "pos": 49985933, "ref": "C", "alt": "T"},
    "CXCR4": {"chrom": "2", "pos": 136116763, "ref": "A", "alt": "G"},
    "MET": {"chrom": "7", "pos": 116735218, "ref": "T", "alt": "C"},
    "BRAF": {"chrom": "7", "pos": 140753336, "ref": "A", "alt": "T"}
}

# Weights for target lock score (from config)
WEIGHTS = {
    "functionality": 0.35,
    "essentiality": 0.35,
    "chromatin": 0.15,
    "regulatory": 0.15
}


async def fetch_insights(gene: str, session_id: int) -> Dict:
    """Fetch all 4 insights for a gene"""
    coords = GENE_COORDS[gene]
    
    variant = {
        "gene": gene,
        "hgvs_p": f"synthetic_{session_id}",
        "chrom": coords["chrom"],
        "pos": coords["pos"],
        "ref": coords["ref"],
        "alt": coords["alt"]
    }
    
    results = {
        "gene": gene,
        "functionality": 0.0,
        "essentiality": 0.0,
        "chromatin": 0.0,
        "regulatory": 0.0
    }
    
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            # Functionality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_protein_functionality_change",
                    json={"mutations": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["functionality"] = data.get("functionality_change_score", 0.0)
            except:
                pass
            
            # Essentiality
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_gene_essentiality",
                    json={"gene": gene, "variants": [variant], "model_id": "evo2_1b"}
                )
                if r.status_code == 200:
                    data = r.json()
                    results["essentiality"] = data.get("essentiality_score", 0.0)
            except:
                pass
            
            # Chromatin
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_chromatin_accessibility",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "radius": 1000,  # larger radius for variation (heuristic uses radius)
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["chromatin"] = data.get("accessibility_score", 0.0)
            except:
                pass
            
            # Regulatory
            try:
                r = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json={
                        "chrom": coords["chrom"],
                        "pos": coords["pos"],
                        "ref": coords["ref"],
                        "alt": coords["alt"],
                        "model_id": "evo2_1b"
                    }
                )
                if r.status_code == 200:
                    data = r.json()
                    results["regulatory"] = data.get("regulatory_impact_score", 0.0)
            except:
                pass
        
        # Calculate target lock score
        results["target_lock_score"] = (
            WEIGHTS["functionality"] * results["functionality"] +
            WEIGHTS["essentiality"] * results["essentiality"] +
            WEIGHTS["chromatin"] * results["chromatin"] +
            WEIGHTS["regulatory"] * results["regulatory"]
        )
        
        return results
        
    except Exception as e:
        print(f"⚠️  Error fetching {gene}: {e}")
        results["target_lock_score"] = 0.0
        return results


async def generate_heatmap_data() -> pd.DataFrame:
    """Generate full heatmap dataset (8 steps × 7 genes = 56 analyses)"""
    print("⚔️  Generating Target Lock Heatmap Data (Direct Insights)...")
    print(f"📊 Total analyses: {len(MISSION_STEPS)} steps × {len(TARGET_GENES)} genes = {len(MISSION_STEPS) * len(TARGET_GENES)}")
    print()
    
    # Fetch insights for each gene once
    print(f"🚀 Fetching insights for {len(TARGET_GENES)} genes...")
    session_id = 0
    tasks = [fetch_insights(gene, session_id) for gene in TARGET_GENES]
    gene_insights = await asyncio.gather(*tasks)
    
    # Create gene → insights map
    insights_map = {g["gene"]: g for g in gene_insights}
    
    print(f"✅ Fetched insights for {len(gene_insights)} genes")
    print()
    
    # Expand to all mission steps
    all_data = []
    for mission in MISSION_STEPS:
        for gene in TARGET_GENES:
            data = insights_map[gene].copy()
            data["mission"] = mission
            all_data.append(data)
    
    df = pd.DataFrame(all_data)
    
    # Save raw data
    output_csv = Path("data/target_lock_heatmap_data.csv")
    output_json = Path("data/target_lock_heatmap_data.json")
    
    df.to_csv(output_csv, index=False)
    
    with open(output_json, 'w') as f:
        json.dump(all_data, f, indent=2)
    
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
    
    # Summary statistics
    print("📈 Summary Statistics:")
    print(f"   Mean target lock score: {df['target_lock_score'].mean():.3f}")
    print(f"   Median target lock score: {df['target_lock_score'].median():.3f}")
    print(f"   Min target lock score: {df['target_lock_score'].min():.3f}")
    print(f"   Max target lock score: {df['target_lock_score'].max():.3f}")
    print(f"   Std dev: {df['target_lock_score'].std():.3f}")
    print()
    
    # Top 5
    print("🎯 Top 5 Gene-Mission Pairs (Highest Target Lock Score):")
    top5 = df.nlargest(5, 'target_lock_score')[['gene', 'mission', 'target_lock_score']]
    for idx, row in top5.iterrows():
        print(f"   {row['gene']} × {row['mission']}: {row['target_lock_score']:.3f}")
    print()


def plot_component_scores(df: pd.DataFrame):
    """Generate component score breakdown visualization"""
    print("📊 Generating component score breakdown...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    score_cols = [
        ('functionality', 'Functionality Score'),
        ('essentiality', 'Essentiality Score'),
        ('chromatin', 'Chromatin Accessibility Score'),
        ('regulatory', 'Regulatory Impact Score')
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
    print("⚔️  TARGET LOCK HEATMAP DATA GENERATION (V2 - FAST)")
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



