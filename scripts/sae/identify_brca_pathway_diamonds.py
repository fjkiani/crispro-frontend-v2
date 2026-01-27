#!/usr/bin/env python3
"""
Identify BRCA Pathway Diamonds - Deliverable #5
===============================================

Mission: Identify BRCA-specific diamond features for recurrence prediction
Objective: Find proliferation, DDR, and immune pathway features
Timeline: 8 hours (overnight)
Status: EXECUTION READY

Based on: SAE_BRCA_NEXT_10_DELIVERABLES.md (Deliverable #5)
Pathways: DDR, Proliferation (Oncotype DX), Immune
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Set, Tuple
from datetime import datetime
from collections import defaultdict, Counter
import numpy as np
from scipy import stats
from sklearn.model_selection import KFold
import statistics

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "brca_diamond_features.json"

# Gene→Pathway Mappings (expanded for BRCA)
GENE_TO_PATHWAY = {
    # DDR pathway (BRCA-specific)
    "BRCA1": ["ddr"], "BRCA2": ["ddr"], "ATM": ["ddr"], "ATR": ["ddr"],
    "CHEK2": ["ddr"], "PALB2": ["ddr"], "RAD51": ["ddr"], "RAD51C": ["ddr"],
    "RAD51D": ["ddr"], "BRIP1": ["ddr"], "BARD1": ["ddr"], "TP53": ["ddr", "tp53"],
    "MBD4": ["ddr"], "MLH1": ["ddr"], "MSH2": ["ddr"], "MSH6": ["ddr"], "PMS2": ["ddr"],
    
    # Proliferation pathway (Oncotype DX genes)
    "MKI67": ["proliferation"], "AURKA": ["proliferation"], "BIRC5": ["proliferation"],
    "CCNB1": ["proliferation"], "MYBL2": ["proliferation"], "STK15": ["proliferation"],
    "CCND1": ["proliferation"], "CCNE1": ["proliferation"], "CDK4": ["proliferation"],
    "CDK6": ["proliferation"], "RB1": ["proliferation"], "BCL2": ["proliferation"],
    "BAX": ["proliferation"], "CASP3": ["proliferation"], "MYC": ["proliferation"],
    "MYCN": ["proliferation"],
    
    # Immune pathway (BRCA IO response)
    "CD8A": ["immune"], "CD274": ["immune"], "PDCD1": ["immune"], "CTLA4": ["immune"],
    "LAG3": ["immune"], "TIGIT": ["immune"], "IFNG": ["immune"], "GZMB": ["immune"],
    "PRF1": ["immune"], "TNF": ["immune"], "IL2": ["immune"]
}

# Statistical thresholds
P_VALUE_THRESHOLD = 0.05
COHENS_D_THRESHOLD = 0.5  # Medium-large effect
TOP_N_FEATURES = 100


def build_feature_matrix(patients: List[Dict]) -> Tuple[np.ndarray, List[str]]:
    """
    Build feature matrix from BRCA SAE features.
    
    Returns:
        feature_matrix: [n_patients, 32768] array
        patient_ids: List of patient IDs
    """
    n_patients = len(patients)
    n_features = 32768
    
    feature_matrix = np.zeros((n_patients, n_features))
    patient_ids = []
    
    for i, patient in enumerate(patients):
        patient_id = patient.get("patient_id") or patient.get("tcga_patient_id")
        patient_ids.append(patient_id)
        
        # Aggregate features from all variants
        variants = patient.get("variants", [])
        for variant in variants:
            top_features = variant.get("top_features", [])
            for tf in top_features:
                idx = tf.get("index")
                val = tf.get("value", 0.0)
                if idx is not None and 0 <= idx < n_features:
                    # Sum across variants (as per OV extraction)
                    feature_matrix[i, idx] += abs(float(val))
        
        # Normalize by variant count (optional - try both)
        if len(variants) > 0:
            feature_matrix[i, :] = feature_matrix[i, :] / len(variants)
    
    return feature_matrix, patient_ids


def compute_pearson_correlation(feature_matrix: np.ndarray, outcome_vector: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Pearson correlation for each feature."""
    n_features = feature_matrix.shape[1]
    pearson_r = np.zeros(n_features)
    pearson_p = np.ones(n_features)
    
    for i in range(n_features):
        if i % 1000 == 0:
            print(f"   Computing Pearson for feature {i}/{n_features}...")
        
        feature_values = feature_matrix[:, i]
        
        # Skip if constant
        if np.std(feature_values) == 0:
            continue
        
        r, p = stats.pearsonr(feature_values, outcome_vector)
        pearson_r[i] = r
        pearson_p[i] = p
    
    return pearson_r, pearson_p


def compute_cohen_d(feature_matrix: np.ndarray, outcome_labels: List[bool]) -> np.ndarray:
    """Compute Cohen's d effect size (recurrence vs no recurrence)."""
    n_features = feature_matrix.shape[1]
    cohen_d_values = np.zeros(n_features)
    
    # Split into recurrence and no-recurrence groups
    recurrence_mask = np.array(outcome_labels)
    no_recurrence_mask = ~recurrence_mask
    
    for i in range(n_features):
        if i % 1000 == 0:
            print(f"   Computing Cohen's d for feature {i}/{n_features}...")
        
        feature_values = feature_matrix[:, i]
        
        recurrence_values = feature_values[recurrence_mask]
        no_recurrence_values = feature_values[no_recurrence_mask]
        
        if len(recurrence_values) == 0 or len(no_recurrence_values) == 0:
            continue
        
        # Cohen's d = (mean1 - mean2) / pooled_std
        mean_diff = np.mean(recurrence_values) - np.mean(no_recurrence_values)
        pooled_std = np.sqrt(
            (np.var(recurrence_values) + np.var(no_recurrence_values)) / 2
        )
        
        if pooled_std == 0:
            continue
        
        cohen_d_values[i] = abs(mean_diff / pooled_std)
    
    return cohen_d_values


def find_patients_with_feature(patients: List[Dict], feature_index: int, top_n: int = 30) -> List[int]:
    """Find top N patients with highest activation of a feature."""
    patient_scores = []
    
    for i, patient in enumerate(patients):
        score = 0.0
        variants = patient.get("variants", [])
        for variant in variants:
            for tf in variant.get("top_features", []):
                if tf.get("index") == feature_index:
                    score += abs(float(tf.get("value", 0.0)))
        
        if score > 0:
            patient_scores.append((i, score))
    
    # Sort by score, take top N
    patient_scores.sort(key=lambda x: x[1], reverse=True)
    return [i for i, _ in patient_scores[:top_n]]


def extract_genes_from_patients(patients: List[Dict], patient_indices: List[int]) -> Set[str]:
    """Extract genes from high-activation patients."""
    genes = set()
    
    for i in patient_indices:
        if i >= len(patients):
            continue
        
        patient = patients[i]
        variants = patient.get("variants", [])
        for variant in variants:
            # Try multiple gene field names
            gene = (
                variant.get("gene") or
                variant.get("hugoGeneSymbol") or
                variant.get("geneSymbol") or
                (variant.get("variant", {}).get("gene") if isinstance(variant.get("variant"), dict) else None)
            )
            
            if gene:
                genes.add(str(gene).upper())
    
    return genes


def infer_pathway_from_genes(genes: Set[str]) -> Dict[str, float]:
    """Infer pathway weights from genes."""
    pathway_counts = Counter()
    
    for gene in genes:
        pathways = GENE_TO_PATHWAY.get(gene, [])
        for pathway in pathways:
            pathway_counts[pathway] += 1
    
    # Normalize to weights (0-1)
    total = sum(pathway_counts.values())
    if total == 0:
        return {}
    
    pathway_weights = {pathway: count / total for pathway, count in pathway_counts.items()}
    
    # Only return pathways with weight >= 0.1 (at least 10% of genes)
    return {p: w for p, w in pathway_weights.items() if w >= 0.1}


def identify_brca_diamonds() -> Dict[str, Any]:
    """Identify BRCA-specific diamond features for recurrence."""
    print("=" * 80)
    print("🔬 IDENTIFYING BRCA PATHWAY DIAMONDS - Deliverable #5")
    print("=" * 80)
    print()
    
    # Step 1: Load BRCA SAE features
    print("📥 Step 1: Loading BRCA SAE features...")
    if not BRCA_CHECKPOINT.exists():
        print(f"❌ BRCA checkpoint not found: {BRCA_CHECKPOINT}")
        sys.exit(1)
    
    with open(BRCA_CHECKPOINT, 'r') as f:
        brca_data = json.load(f)
    
    # Convert to list format (for compatibility with existing code)
    patients = []
    if "data" in brca_data and isinstance(brca_data["data"], dict):
        for patient_id, patient_data in brca_data["data"].items():
            patient_data["patient_id"] = patient_id
            patients.append(patient_data)
    else:
        patients = brca_data.get("patients", [])
    
    print(f"   ✅ Loaded {len(patients)} patients")
    
    # Step 2: Load recurrence labels
    print()
    print("📥 Step 2: Loading recurrence labels...")
    if not RECURRENCE_LABELS.exists():
        print(f"❌ Recurrence labels not found: {RECURRENCE_LABELS}")
        sys.exit(1)
    
    with open(RECURRENCE_LABELS, 'r') as f:
        labels_data = json.load(f)
    
    outcomes = labels_data.get("outcomes", {})
    print(f"   ✅ Loaded {len(outcomes)} outcome records")
    
    # Step 3: Build feature matrix and outcome vector
    print()
    print("🔍 Step 3: Building feature matrix...")
    feature_matrix, patient_ids = build_feature_matrix(patients)
    print(f"   ✅ Feature matrix: {feature_matrix.shape}")
    
    # Build outcome vector (binary: True=1.0, False=0.0)
    outcome_vector = []
    outcome_labels = []
    valid_patient_indices = []
    
    for i, patient_id in enumerate(patient_ids):
        outcome = outcomes.get(patient_id, {})
        recurrence = outcome.get("recurrence")
        
        if recurrence is None:
            continue  # Skip patients without labels
        
        outcome_vector.append(1.0 if recurrence else 0.0)
        outcome_labels.append(recurrence)
        valid_patient_indices.append(i)
    
    # Filter feature matrix to valid patients
    feature_matrix = feature_matrix[valid_patient_indices, :]
    outcome_vector = np.array(outcome_vector)
    
    print(f"   ✅ Valid patients: {len(outcome_vector)}")
    print(f"   Recurrence: {sum(outcome_labels)} True, {len(outcome_labels) - sum(outcome_labels)} False")
    
    # Step 4: Compute correlations
    print()
    print("📊 Step 4: Computing correlations (this may take a while)...")
    print("   This will analyze all 32,768 features...")
    
    pearson_r, pearson_p = compute_pearson_correlation(feature_matrix, outcome_vector)
    print("   ✅ Pearson correlation complete")
    
    cohen_d = compute_cohen_d(feature_matrix, outcome_labels)
    print("   ✅ Cohen's d complete")
    
    # Step 5: Filter significant features
    print()
    print("🔍 Step 5: Filtering significant features...")
    significant_features = []
    
    for i in range(feature_matrix.shape[1]):
        if pearson_p[i] < P_VALUE_THRESHOLD and cohen_d[i] >= COHENS_D_THRESHOLD:
            significant_features.append({
                "index": i,
                "pearson_r": float(pearson_r[i]),
                "pearson_p": float(pearson_p[i]),
                "cohen_d": float(cohen_d[i])
            })
    
    # Sort by Cohen's d (descending)
    significant_features.sort(key=lambda x: x["cohen_d"], reverse=True)
    
    # Take top N
    top_features = significant_features[:TOP_N_FEATURES]
    
    print(f"   ✅ Found {len(significant_features)} significant features")
    print(f"   ✅ Selected top {len(top_features)} features")
    
    # Step 6: Map features to pathways
    print()
    print("🗺️  Step 6: Mapping features to pathways...")
    pathway_diamonds = {
        "ddr": [],
        "proliferation": [],
        "immune": []
    }
    
    for feature in top_features:
        feature_index = feature["index"]
        
        # Find patients with high activation
        high_activation_patients = find_patients_with_feature(patients, feature_index, top_n=30)
        
        if not high_activation_patients:
            continue
        
        # Extract genes from those patients
        genes = extract_genes_from_patients(patients, high_activation_patients)
        
        if not genes:
            continue
        
        # Infer pathway
        pathway_weights = infer_pathway_from_genes(genes)
        
        if not pathway_weights:
            continue
        
        # Assign to primary pathway (highest weight)
        primary_pathway = max(pathway_weights.items(), key=lambda x: x[1])[0]
        
        if primary_pathway in pathway_diamonds:
            pathway_diamonds[primary_pathway].append({
                "index": feature_index,
                "cohen_d": feature["cohen_d"],
                "pearson_r": feature["pearson_r"],
                "pearson_p": feature["pearson_p"],
                "pathway_weight": pathway_weights[primary_pathway],
                "genes": list(genes),
                "n_high_activation_patients": len(high_activation_patients)
            })
    
    # Sort each pathway by Cohen's d
    for pathway in pathway_diamonds:
        pathway_diamonds[pathway].sort(key=lambda x: x["cohen_d"], reverse=True)
    
    print(f"   ✅ DDR diamonds: {len(pathway_diamonds['ddr'])}")
    print(f"   ✅ Proliferation diamonds: {len(pathway_diamonds['proliferation'])}")
    print(f"   ✅ Immune diamonds: {len(pathway_diamonds['immune'])}")
    
    # Step 7: Save results
    print()
    print("💾 Step 7: Saving results...")
    output_data = {
        "extraction_date": datetime.now().isoformat(),
        "n_patients": len(outcome_vector),
        "n_features_analyzed": feature_matrix.shape[1],
        "n_significant_features": len(significant_features),
        "n_top_features": len(top_features),
        "pathway_diamonds": pathway_diamonds,
        "statistical_thresholds": {
            "p_value": P_VALUE_THRESHOLD,
            "cohen_d": COHENS_D_THRESHOLD,
            "top_n": TOP_N_FEATURES
        },
        "provenance": {
            "brca_checkpoint": str(BRCA_CHECKPOINT),
            "recurrence_labels": str(RECURRENCE_LABELS),
            "outcome": "recurrence (binary: True/False)"
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"   ✅ Saved to: {OUTPUT_FILE}")
    
    # Step 8: Print summary
    print()
    print("=" * 80)
    print("📊 DIAMOND FEATURES SUMMARY")
    print("=" * 80)
    print(f"Total significant features: {len(significant_features)}")
    print(f"Top features analyzed: {len(top_features)}")
    print()
    print("Pathway diamonds:")
    print(f"  DDR: {len(pathway_diamonds['ddr'])} features")
    if pathway_diamonds['ddr']:
        print(f"    Top DDR: Feature {pathway_diamonds['ddr'][0]['index']} (Cohen's d: {pathway_diamonds['ddr'][0]['cohen_d']:.3f})")
    print(f"  Proliferation: {len(pathway_diamonds['proliferation'])} features")
    if pathway_diamonds['proliferation']:
        print(f"    Top Proliferation: Feature {pathway_diamonds['proliferation'][0]['index']} (Cohen's d: {pathway_diamonds['proliferation'][0]['cohen_d']:.3f})")
    print(f"  Immune: {len(pathway_diamonds['immune'])} features")
    if pathway_diamonds['immune']:
        print(f"    Top Immune: Feature {pathway_diamonds['immune'][0]['index']} (Cohen's d: {pathway_diamonds['immune'][0]['cohen_d']:.3f})")
    print("=" * 80)
    
    return output_data


if __name__ == "__main__":
    try:
        result = identify_brca_diamonds()
        print("\n✅ DIAMOND IDENTIFICATION COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ DIAMOND IDENTIFICATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
