#!/usr/bin/env python3
"""
🧬 MM RESISTANCE PREDICTION - Work Item 1: Mine MM Diamond Features
====================================================================

Mission: Identify MM-specific "diamond" features from TRUE SAE
Agent: Plumber (Implementation)
Date: January 29, 2025
Priority: P0 (Critical)

Objective:
- Load MMRF CoMMpass SAE features (from extract_mmrf_sae_features.py)
- Define MM resistance labels (PI-resistant, IMiD-resistant, Universal-resistant)
- For each resistance type:
  - Identify features with significant differential activation (Cohen's d > 0.5, p < 0.05)
  - Map features to MM-specific pathways (Proteasome/UPR, Cereblon, NRF2, MAPK, Plasma Cell Survival, Drug Efflux)
- Save identified diamond features and their pathway mappings

Success Criteria:
- Identify ≥3 MM diamond features for at least one resistance type
- Generate report with feature IDs, effect sizes, p-values, and pathway mappings

Guardrails:
- RUO/validation-only (no production impact)
- Transparent reporting of feature selection criteria
"""

import json
import sys
import os
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from loguru import logger
from collections import defaultdict

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")

# --- Configuration ---
SAE_FEATURES_FILE = Path("data/validation/mm_cohort/mmrf_sae_features.json")
OUTPUT_DIR = Path("data/validation/mm_cohort/diamond_features")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
DIAMOND_FEATURES_FILE = OUTPUT_DIR / "mm_diamond_features.json"
REPORT_FILE = OUTPUT_DIR / "mm_diamond_features_report.md"

# --- Feature Selection Criteria ---
COHEN_D_THRESHOLD = 0.5
P_VALUE_THRESHOLD = 0.05
FDR_ALPHA = 0.10

# --- MM-specific Pathway Genes ---
MM_PATHWAY_GENES = {
    "proteasome_upr": ["PSMB5", "PSMB8", "XBP1", "IRE1", "ATF6", "DDIT3"],
    "cereblon_pathway": ["CRBN", "CUL4A", "DDB1", "IKZF1", "IKZF3"],
    "nrf2_antioxidant": ["NFE2L2", "KEAP1", "NQO1", "HMOX1"],
    "ras_mapk": ["KRAS", "NRAS", "BRAF", "MAP2K1"],
    "plasma_cell_survival": ["BCL2", "MCL1", "BIRC5", "MYC"],
    "drug_efflux": ["ABCB1", "ABCC1", "ABCG2"],
}

def get_mm_resistance_labels(patient_data: Dict[str, Any]) -> Dict[str, bool]:
    """Determine MM resistance labels for a patient."""
    labels = {
        "PI_RESISTANT": False,
        "IMID_RESISTANT": False,
        "UNIVERSAL_RESISTANT": False,
        "PI_SENSITIVE": False,
        "IMID_SENSITIVE": False,
    }

    mutations = patient_data.get("mutations", [])
    has_psmb5 = any(m.get("gene") == "PSMB5" for m in mutations)
    has_crbn = any(m.get("gene") == "CRBN" for m in mutations)
    has_tp53 = any(m.get("gene") == "TP53" for m in mutations)

    if has_psmb5:
        labels["PI_RESISTANT"] = True
    if has_crbn:
        labels["IMID_RESISTANT"] = True
    if has_tp53:
        labels["UNIVERSAL_RESISTANT"] = True
    
    if not labels["PI_RESISTANT"] and not labels["IMID_RESISTANT"]:
        labels["PI_SENSITIVE"] = True
        labels["IMID_SENSITIVE"] = True

    return labels

def cohen_d(x, y):
    """Calculate Cohen's d for two samples."""
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    if dof <= 0:
        return 0.0
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1) ** 2 + (ny-1)*np.std(y, ddof=1) ** 2) / dof)

def map_feature_to_pathways(feature_id: str, patients: List[Dict[str, Any]]) -> List[str]:
    """Map a feature to biological pathways by analyzing associated genes."""
    associated_genes = defaultdict(int)
    for patient in patients:
        feature_value = patient.get("aggregated_features", {}).get("mean_features", {}).get(feature_id, 0.0)
        all_feature_values = [p.get("aggregated_features", {}).get("mean_features", {}).get(feature_id, 0.0) for p in patients]
        if not all_feature_values:
            continue
        
        threshold = np.percentile(all_feature_values, 90)
        if feature_value >= threshold:
            for mut in patient.get("mutations", []):
                gene = mut.get("gene")
                if gene:
                    associated_genes[gene] += 1
    
    mapped_pathways = set()
    for pathway, genes_in_pathway in MM_PATHWAY_GENES.items():
        for gene in genes_in_pathway:
            if gene in associated_genes and associated_genes[gene] > 0:
                mapped_pathways.add(pathway)
                break
    
    return sorted(list(mapped_pathways))

def mine_diamond_features():
    logger.info("=" * 80)
    logger.info("💎 Mining MM Diamond Features from TRUE SAE")
    logger.info("=" * 80)

    if not SAE_FEATURES_FILE.exists():
        logger.error(f"❌ SAE features file not found: {SAE_FEATURES_FILE}. Please run extract_mmrf_sae_features.py first.")
        sys.exit(1)

    with open(SAE_FEATURES_FILE, 'r') as f:
        data = json.load(f)
    
    patients = data.get("patients", [])
    logger.info(f"Loaded {len(patients)} patients with SAE features.")

    all_feature_ids = sorted(list(set(f_id for p in patients for f_id in p.get("aggregated_features", {}).get("mean_features", {}).keys())))
    
    feature_data_by_resistance_type = defaultdict(lambda: defaultdict(list))
    feature_data_by_sensitivity_type = defaultdict(lambda: defaultdict(list))

    for patient in patients:
        labels = get_mm_resistance_labels(patient)
        mean_features = patient.get("aggregated_features", {}).get("mean_features", {})

        for feature_id in all_feature_ids:
            value = mean_features.get(feature_id, 0.0)
            
            if labels["PI_RESISTANT"]:
                feature_data_by_resistance_type["PI_RESISTANT"][feature_id].append(value)
            if labels["IMID_RESISTANT"]:
                feature_data_by_resistance_type["IMID_RESISTANT"][feature_id].append(value)
            
            if labels["PI_SENSITIVE"]:
                feature_data_by_sensitivity_type["PI_SENSITIVE"][feature_id].append(value)
            if labels["IMID_SENSITIVE"]:
                feature_data_by_sensitivity_type["IMID_SENSITIVE"][feature_id].append(value)

    diamond_features_results = {
        "metadata": {
            "run_date": datetime.now().isoformat(),
            "num_patients": len(patients),
            "total_sae_features": len(all_feature_ids),
            "cohen_d_threshold": COHEN_D_THRESHOLD,
            "p_value_threshold": P_VALUE_THRESHOLD,
            "fdr_alpha": FDR_ALPHA
        },
        "resistance_types": {}
    }

    resistance_types_to_mine = {
        "PI_RESISTANT": "PI_SENSITIVE",
        "IMID_RESISTANT": "IMID_SENSITIVE",
    }

    for res_type, sens_type in resistance_types_to_mine.items():
        logger.info(f"\n--- Mining diamonds for {res_type} vs {sens_type} ---")
        
        resistant_feature_data = feature_data_by_resistance_type[res_type]
        sensitive_feature_data = feature_data_by_sensitivity_type[sens_type]

        if not resistant_feature_data or not sensitive_feature_data:
            logger.warning(f"⚠️ Skipping {res_type}: Insufficient data")
            diamond_features_results["resistance_types"][res_type] = {"status": "SKIPPED", "reason": "Insufficient data"}
            continue

        candidate_features = []
        p_values = []
        
        for feature_id in all_feature_ids:
            res_values = resistant_feature_data.get(feature_id, [])
            sens_values = sensitive_feature_data.get(feature_id, [])

            if len(res_values) < 2 or len(sens_values) < 2:
                continue

            d_score = cohen_d(res_values, sens_values)
            
            try:
                u_stat, p_val = mannwhitneyu(res_values, sens_values, alternative='two-sided')
            except ValueError:
                p_val = 1.0

            if abs(d_score) >= COHEN_D_THRESHOLD and p_val < P_VALUE_THRESHOLD:
                candidate_features.append({
                    "feature_id": feature_id,
                    "cohen_d": d_score,
                    "p_value": p_val,
                    "mean_resistant": np.mean(res_values),
                    "mean_sensitive": np.mean(sens_values),
                    "direction": "higher_in_resistant" if d_score > 0 else "higher_in_sensitive",
                    "pathways": map_feature_to_pathways(feature_id, patients)
                })
                p_values.append(p_val)
        
        if p_values:
            reject, pvals_corrected, _, _ = multipletests(p_values, alpha=FDR_ALPHA, method='fdr_bh')
            for i, feature in enumerate(candidate_features):
                feature["fdr_corrected_p_value"] = pvals_corrected[i]
                feature["fdr_significant"] = reject[i]

        diamond_features = [
            f for f in candidate_features
            if f.get("fdr_significant", False) and f["direction"] == "higher_in_resistant"
        ]
        
        diamond_features.sort(key=lambda x: x["cohen_d"], reverse=True)

        logger.info(f"    Found {len(diamond_features)} diamond features for {res_type}.")
        diamond_features_results["resistance_types"][res_type] = {
            "status": "SUCCESS" if len(diamond_features) >= 3 else "BELOW_TARGET",
            "diamond_features": diamond_features
        }

    # Generate Report
    with open(REPORT_FILE, 'w') as f:
        f.write(f"# MM Diamond Features Mining Report\n\n")
        f.write(f"**Run Date:** {diamond_features_results['metadata']['run_date']}\n")
        f.write(f"**Number of Patients:** {diamond_features_results['metadata']['num_patients']}\n")
        f.write(f"**Total SAE Features:** {diamond_features_results['metadata']['total_sae_features']}\n\n")
        f.write("## Diamond Features by Resistance Type\n\n")
        for res_type, type_results in diamond_features_results["resistance_types"].items():
            f.write(f"### {res_type.replace('_', ' ')}\n\n")
            if type_results.get("status") == "SKIPPED":
                f.write(f"**Status:** SKIPPED - {type_results['reason']}\n\n")
                continue
            
            diamond_features = type_results.get("diamond_features", [])
            if not diamond_features:
                f.write("**Status:** BELOW_TARGET - No diamond features found.\n\n")
                continue

            f.write(f"**Status:** {type_results['status']} ({len(diamond_features)} features found)\n\n")
            f.write("| Feature ID | Cohen's d | p-value | Direction | Pathways |\n")
            f.write("|------------|-----------|---------|-----------|----------|\n")
            for feature in diamond_features:
                f.write(f"| {feature['feature_id']} | {feature['cohen_d']:.3f} | {feature['p_value']:.4f} | {feature['direction']} | {', '.join(feature['pathways'])} |\n")
            f.write("\n")

    with open(DIAMOND_FEATURES_FILE, 'w') as f:
        json.dump(diamond_features_results, f, indent=2)
    logger.info(f"✅ Diamond features saved to {DIAMOND_FEATURES_FILE}")

if __name__ == "__main__":
    mine_diamond_features()
