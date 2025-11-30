#!/usr/bin/env python3
"""
⚔️ MOCK SAE COHORT GENERATOR (DEVELOPMENT ONLY)
================================================

Agent: Zo (Lead Commander)
Date: January 15, 2025
Purpose: Generate synthetic SAE features for pipeline testing

⚠️ WARNING: This is MOCK DATA for development and testing only.
DO NOT use for clinical decisions or scientific conclusions.

This script creates synthetic SAE features that mimic the structure
of real Evo2+SAE outputs to verify the biomarker correlation pipeline works.

Real biological signals:
- Features 100-102: Synthetic "DDR-like" features (higher in sensitive)
- Features 200-202: Synthetic "MAPK-like" features (higher in resistant)
- Features 300-302: Synthetic "IO-like" features (correlated with outcome)

All other features: Random noise

Output: data/validation/sae_cohort/sae_features_tcga_ov_platinum_MOCK.json
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any

# Configuration
PLATINUM_LABELS_FILE = Path("data/validation/tcga_ov_platinum_response_labels.json")
OUTPUT_FILE = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum_MOCK.json")
SAE_D_HIDDEN = 32768  # 32K SAE features
RANDOM_SEED = 42

# Synthetic "signal" feature ranges
DDR_FEATURES = [100, 101, 102]
MAPK_FEATURES = [200, 201, 202]
IO_FEATURES = [300, 301, 302]

def generate_mock_sae_features(
    outcome: str,
    num_mutations: int = 10
) -> Dict[str, Any]:
    """
    Generates mock SAE features for a patient based on outcome.
    
    Synthetic signals:
    - DDR features (100-102): Higher in sensitive patients
    - MAPK features (200-202): Higher in resistant/refractory patients
    - IO features (300-302): Correlated with outcome spectrum
    
    Args:
        outcome: "sensitive", "resistant", or "refractory"
        num_mutations: Number of mutations for this patient
    
    Returns:
        Dict with mock SAE feature aggregates
    """
    # Outcome encoding for signal strength
    outcome_signal = {
        "sensitive": 1.0,
        "resistant": 0.5,
        "refractory": 0.0
    }.get(outcome, 0.5)
    
    # Generate per-mutation features
    mutation_features = []
    
    for _ in range(num_mutations):
        features = np.random.randn(SAE_D_HIDDEN) * 0.1  # Small random noise
        
        # Inject synthetic signals
        # DDR features: higher in sensitive (positive correlation with outcome)
        for idx in DDR_FEATURES:
            features[idx] = np.random.normal(outcome_signal * 2.0, 0.3)
        
        # MAPK features: higher in resistant/refractory (negative correlation)
        for idx in MAPK_FEATURES:
            features[idx] = np.random.normal((1 - outcome_signal) * 2.0, 0.3)
        
        # IO features: moderate correlation
        for idx in IO_FEATURES:
            features[idx] = np.random.normal(outcome_signal * 1.5, 0.5)
        
        mutation_features.append(features.tolist())
    
    # Aggregate across mutations (mean, max)
    mutation_features_np = np.array(mutation_features)
    mean_features = np.mean(mutation_features_np, axis=0).tolist()
    max_features = np.max(mutation_features_np, axis=0).tolist()
    
    # Top-k features (find indices with highest activation in mean)
    mean_np = np.array(mean_features)
    top_k_indices = np.argsort(np.abs(mean_np))[-64:][::-1]  # Top 64 features (SAE k=64)
    top_k_counts = {int(idx): 1 for idx in top_k_indices}
    
    return {
        "mean_features": mean_features,
        "max_features": max_features,
        "top_k_counts": top_k_counts,
        "num_mutations": num_mutations,
        "provenance": {
            "aggregation_method": "mean_max_topk_counts",
            "script": "scripts/sae/generate_mock_sae_cohort.py (MOCK DATA)",
            "timestamp": datetime.now().isoformat(),
            "mock_signals": {
                "ddr_features": DDR_FEATURES,
                "mapk_features": MAPK_FEATURES,
                "io_features": IO_FEATURES
            }
        }
    }

def generate_mock_mutations(num_mutations: int = 10) -> List[Dict]:
    """Generates mock mutation records for a patient."""
    mutations = []
    genes = ["TP53", "BRCA1", "BRCA2", "KRAS", "PIK3CA", "PTEN", "RB1", "NF1"]
    
    for i in range(num_mutations):
        mutations.append({
            "gene": genes[i % len(genes)],
            "chrom": str(np.random.choice([1, 2, 3, 7, 13, 17, 19])),
            "pos": np.random.randint(1000000, 100000000),
            "ref": np.random.choice(["A", "T", "C", "G"]),
            "alt": np.random.choice(["A", "T", "C", "G"]),
            "hgvs_p": f"p.{np.random.choice(['Arg', 'Lys', 'Glu', 'Asp'])}{np.random.randint(1, 500)}{np.random.choice(['Ter', 'Gly', 'Ala'])}",
            "assembly": "GRCh38",
            "mock": True
        })
    
    return mutations

def main():
    """Main mock data generation."""
    print("="*80)
    print("⚔️ MOCK SAE COHORT GENERATOR")
    print("="*80)
    print("⚠️  WARNING: Generating MOCK DATA for development testing only!")
    print()
    
    # Set random seed
    np.random.seed(RANDOM_SEED)
    
    # Load platinum labels
    if not PLATINUM_LABELS_FILE.exists():
        print(f"❌ Platinum labels file not found: {PLATINUM_LABELS_FILE}")
        print("   Run scripts/extract_platinum_response_labels.py first")
        return
    
    with open(PLATINUM_LABELS_FILE, 'r') as f:
        platinum_data = json.load(f)
    
    patients = platinum_data.get("patients", [])
    print(f"📥 Loaded {len(patients)} patients from platinum labels")
    
    # Filter for fully labeled patients
    labeled_patients = [
        p for p in patients
        if p.get("platinum_response") in ["sensitive", "resistant", "refractory"]
    ]
    print(f"   Filtered to {len(labeled_patients)} fully labeled patients")
    
    # Generate mock SAE features for each patient
    mock_cohort = []
    
    for i, patient in enumerate(labeled_patients):
        if i % 10 == 0:
            print(f"   Generating mock features for patient {i+1}/{len(labeled_patients)}...")
        
        outcome = patient["platinum_response"]
        num_mutations = np.random.randint(5, 15)  # Random 5-15 mutations per patient
        
        # Outcome encoding for signal strength
        outcome_signal = {
            "sensitive": 1.0,
            "resistant": 0.5,
            "refractory": 0.0
        }.get(outcome, 0.5)
        
        # Generate mock mutations
        mutations = generate_mock_mutations(num_mutations)
        
        # Create variants with top_features (matching real SAE extraction format)
        variants = []
        all_variant_features = []
        
        for mut in mutations:
            # Generate top_features for this variant
            variant_features = np.random.randn(SAE_D_HIDDEN) * 0.1
            # Inject synthetic signals
            for idx in DDR_FEATURES:
                variant_features[idx] = np.random.normal(outcome_signal * 2.0, 0.3)
            for idx in MAPK_FEATURES:
                variant_features[idx] = np.random.normal((1 - outcome_signal) * 2.0, 0.3)
            for idx in IO_FEATURES:
                variant_features[idx] = np.random.normal(outcome_signal * 1.5, 0.5)
            
            all_variant_features.append(variant_features)
            
            # Get top 64 features
            top_k_indices = np.argsort(np.abs(variant_features))[-64:][::-1]
            top_features = [
                {"index": int(idx), "value": float(variant_features[idx])}
                for idx in top_k_indices
            ]
            
            variants.append({
                "variant": {
                    "chrom": mut["chrom"],
                    "pos": mut["pos"],
                    "ref": mut["ref"],
                    "alt": mut["alt"],
                    "gene": mut.get("gene"),
                    "hgvs_p": mut.get("hgvs_p")
                },
                "sae_features": variant_features.tolist(),
                "top_features": top_features,
                "stats": {
                    "sparsity": 0.002,  # 64/32768
                    "mean_activation": float(np.mean(variant_features[variant_features != 0])) if np.any(variant_features != 0) else 0.0,
                    "num_active_features": 64
                }
            })
        
        # Aggregate features across variants (matching extract_sae_features_cohort.py format)
        if all_variant_features:
            features_array = np.array(all_variant_features)  # Shape: [num_variants, 32768]
            mean_features = np.mean(features_array, axis=0).tolist()
            max_features = np.max(features_array, axis=0).tolist()
            
            # Aggregate top features across variants
            feature_activation_counts = {}
            for vf in variants:
                for top_feat in vf.get("top_features", []):
                    idx = top_feat.get("index")
                    val = top_feat.get("value", 0.0)
                    if idx is not None:
                        feature_activation_counts[idx] = feature_activation_counts.get(idx, 0.0) + val
            
            top_features_aggregated = [
                {"index": int(idx), "value": float(val)}
                for idx, val in sorted(feature_activation_counts.items(), key=lambda x: x[1], reverse=True)[:64]
            ]
        else:
            mean_features = []
            max_features = []
            top_features_aggregated = []
        
        aggregated = {
            "num_variants": len(variants),
            "mean_features": mean_features,
            "max_features": max_features,
            "top_features_aggregated": top_features_aggregated,
            "sparsity_mean": 0.002,
            "sparsity_std": 0.0
        }
        
        mock_cohort.append({
            "patient_id": patient["tcga_patient_id"],
            "platinum_response": outcome,
            "variants": variants,
            "aggregated_features": aggregated,
            "provenance": {
                "extraction_script": "scripts/sae/generate_mock_sae_cohort.py (MOCK)",
                "timestamp": datetime.now().isoformat(),
                "mock_data": True,
                "mock_signals": {
                    "ddr_features": DDR_FEATURES,
                    "mapk_features": MAPK_FEATURES,
                    "io_features": IO_FEATURES
                }
            }
        })
    
    # Save mock cohort
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(mock_cohort, f, indent=2)
    
    print(f"\n✅ Generated mock SAE features for {len(mock_cohort)} patients")
    print(f"📄 Saved to: {OUTPUT_FILE}")
    print()
    print("="*80)
    print("⚠️  MOCK DATA DISCLAIMER")
    print("="*80)
    print("This is SYNTHETIC data for pipeline testing only.")
    print("Synthetic signals injected:")
    print(f"  - DDR features ({DDR_FEATURES}): Higher in sensitive patients")
    print(f"  - MAPK features ({MAPK_FEATURES}): Higher in resistant/refractory")
    print(f"  - IO features ({IO_FEATURES}): Correlated with outcome")
    print()
    print("DO NOT use for:")
    print("  ❌ Clinical decisions")
    print("  ❌ Scientific publications")
    print("  ❌ Biomarker validation")
    print()
    print("DO use for:")
    print("  ✅ Pipeline testing")
    print("  ✅ Visualization verification")
    print("  ✅ Statistical method validation")
    print("="*80)
    print()
    print("Next step: Run biomarker analysis on mock data:")
    print(f"  python scripts/sae/analyze_biomarkers.py \\")
    print(f"    --input {OUTPUT_FILE} \\")
    print(f"    --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers_MOCK.json")

if __name__ == "__main__":
    main()



