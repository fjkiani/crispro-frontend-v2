#!/usr/bin/env python3
"""
🧬 MM RESISTANCE PREDICTION - Work Item 1: Validate MM TRUE SAE
================================================================

Mission: Validate TRUE SAE features for MM resistance prediction
Agent: Plumber (Implementation)
Date: January 29, 2025
Priority: P0 (Critical)

Objective:
- Load MMRF CoMMpass SAE features (from extract_mmrf_sae_features.py)
- Define resistance labels (PI-resistant, IMiD-resistant, Universal-resistant)
- Train logistic regression classifiers
- Evaluate performance (AUROC, AUPRC, Sensitivity, Specificity)
- Compare against Proxy SAE (gene-level baseline)

Success Criteria:
- MM TRUE SAE AUROC ≥ 0.70 for at least one resistance type
- Generate comprehensive validation report

Guardrails:
- RUO/validation-only (no production impact)
- Transparent reporting of performance metrics
"""

import json
import sys
import os
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, confusion_matrix, roc_curve
from sklearn.preprocessing import StandardScaler
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")

# --- Configuration ---
SAE_FEATURES_FILE = Path("data/validation/mm_cohort/mmrf_sae_features.json")
OUTPUT_DIR = Path("data/validation/mm_cohort/validation_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_FILE = OUTPUT_DIR / "mmrf_true_sae_validation_report.md"
RESULTS_JSON_FILE = OUTPUT_DIR / "mmrf_true_sae_validation_results.json"

# --- MM-specific Resistance Definitions ---
def get_mm_resistance_labels(patient_data: Dict[str, Any]) -> Dict[str, bool]:
    """
    Determine MM resistance labels for a patient.
    Placeholder logic - needs actual clinical data parsing.
    """
    labels = {
        "PI_RESISTANT": False,
        "IMID_RESISTANT": False,
        "UNIVERSAL_RESISTANT": False,
        "PI_SENSITIVE": False,
        "IMID_SENSITIVE": False,
    }

    clinical_data = patient_data.get("clinical_data", {})
    mutations = patient_data.get("mutations", [])
    
    # Example placeholder logic:
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

# --- Proxy SAE (Gene-level baseline) ---
MM_PROXY_GENES = {
    "DIS3": 2.08,
    "TP53": 1.90,
    "PSMB5": 3.5,
    "CRBN": 5.0,
    "IKZF1": 2.5,
    "IKZF3": 2.5,
}

def compute_proxy_sae_score(mutations: List[Dict[str, Any]], target_genes: List[str]) -> float:
    """Compute a simple proxy SAE score based on presence of target mutations."""
    score = 0.0
    for mut in mutations:
        gene = mut.get("gene")
        if gene in target_genes and gene in MM_PROXY_GENES:
            score += MM_PROXY_GENES[gene]
    return min(score / 5.0, 1.0)

def evaluate_classifier(X: np.ndarray, y: np.ndarray, model_name: str) -> Dict[str, Any]:
    """Evaluate a classifier using stratified K-fold cross-validation."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    auroc_scores = []
    auprc_scores = []
    sensitivity_scores = []
    specificity_scores = []
    all_y_true = []
    all_y_pred_proba = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        model = LogisticRegression(solver='liblinear', random_state=42, class_weight='balanced', max_iter=1000)
        model.fit(X_train_scaled, y_train)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        y_pred = model.predict(X_test_scaled)

        auroc_scores.append(roc_auc_score(y_test, y_pred_proba))
        auprc_scores.append(average_precision_score(y_test, y_pred_proba))
        
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        sensitivity_scores.append(tp / (tp + fn) if (tp + fn) > 0 else 0)
        specificity_scores.append(tn / (tn + fp) if (tn + fp) > 0 else 0)

        all_y_true.extend(y_test)
        all_y_pred_proba.extend(y_pred_proba)

    mean_auroc = np.mean(auroc_scores)
    mean_auprc = np.mean(auprc_scores)
    mean_sensitivity = np.mean(sensitivity_scores)
    mean_specificity = np.mean(specificity_scores)

    logger.info(f"    AUROC: {mean_auroc:.3f} ± {np.std(auroc_scores):.3f}")
    logger.info(f"    AUPRC: {mean_auprc:.3f} ± {np.std(auprc_scores):.3f}")
    logger.info(f"    Sensitivity: {mean_sensitivity:.3f}")
    logger.info(f"    Specificity: {mean_specificity:.3f}")

    return {
        "mean_auroc": float(mean_auroc),
        "std_auroc": float(np.std(auroc_scores)),
        "mean_auprc": float(mean_auprc),
        "std_auprc": float(np.std(auprc_scores)),
        "mean_sensitivity": float(mean_sensitivity),
        "mean_specificity": float(mean_specificity),
        "status": "SUCCESS" if mean_auroc >= 0.70 else "BELOW_TARGET"
    }

def run_validation():
    logger.info("=" * 80)
    logger.info("🧬 MMRF CoMMpass TRUE SAE Validation")
    logger.info("=" * 80)

    if not SAE_FEATURES_FILE.exists():
        logger.error(f"❌ SAE features file not found: {SAE_FEATURES_FILE}. Please run extract_mmrf_sae_features.py first.")
        sys.exit(1)

    with open(SAE_FEATURES_FILE, 'r') as f:
        data = json.load(f)
    
    patients = data.get("patients", [])
    logger.info(f"Loaded {len(patients)} patients with SAE features.")

    # Prepare data for classification
    X_true_sae = []
    X_proxy_sae = []
    y_pi_resistant = []
    y_imid_resistant = []
    y_universal_resistant = []
    patient_ids = []

    for patient in patients:
        labels = get_mm_resistance_labels(patient)
        if labels["PI_SENSITIVE"] or labels["PI_RESISTANT"]:
            patient_ids.append(patient["patient_id"])
            
            aggregated_features = patient.get("aggregated_features", {})
            mean_features = aggregated_features.get("mean_features", {})
            
            all_feature_ids = sorted(list(set(f_id for p in patients for f_id in p.get("aggregated_features", {}).get("mean_features", {}).keys())))
            current_true_sae_vector = [mean_features.get(f_id, 0.0) for f_id in all_feature_ids]
            X_true_sae.append(current_true_sae_vector)

            proxy_pi_genes = ["PSMB5", "NFE2L2", "DIS3", "TP53"]
            proxy_imid_genes = ["CRBN", "IKZF1", "IKZF3", "DIS3", "TP53"]
            
            proxy_pi_score = compute_proxy_sae_score(patient.get("mutations", []), proxy_pi_genes)
            proxy_imid_score = compute_proxy_sae_score(patient.get("mutations", []), proxy_imid_genes)
            
            X_proxy_sae.append([proxy_pi_score, proxy_imid_score])

            y_pi_resistant.append(1 if labels["PI_RESISTANT"] else 0)
            y_imid_resistant.append(1 if labels["IMID_RESISTANT"] else 0)
            y_universal_resistant.append(1 if labels["UNIVERSAL_RESISTANT"] else 0)

    if not patient_ids:
        logger.error("❌ No patients with valid labels and SAE features found. Exiting.")
        sys.exit(1)

    X_true_sae = np.array(X_true_sae)
    X_proxy_sae = np.array(X_proxy_sae)
    y_pi_resistant = np.array(y_pi_resistant)
    y_imid_resistant = np.array(y_imid_resistant)
    y_universal_resistant = np.array(y_universal_resistant)

    logger.info(f"Prepared data for {len(patient_ids)} patients.")
    logger.info(f"TRUE SAE feature dimension: {X_true_sae.shape[1]}")
    logger.info(f"Proxy SAE feature dimension: {X_proxy_sae.shape[1]}")

    results = {
        "metadata": {
            "run_date": datetime.now().isoformat(),
            "num_patients": len(patient_ids),
            "true_sae_feature_dim": X_true_sae.shape[1],
            "proxy_sae_feature_dim": X_proxy_sae.shape[1],
        },
        "resistance_types": {}
    }

    # Run validation for each resistance type
    resistance_types = {
        "PI_RESISTANT": y_pi_resistant,
        "IMID_RESISTANT": y_imid_resistant,
        "UNIVERSAL_RESISTANT": y_universal_resistant,
    }

    for res_type, y_labels in resistance_types.items():
        logger.info(f"\n--- Validating for {res_type} ---")
        if np.sum(y_labels) == 0 or np.sum(y_labels) == len(y_labels):
            logger.warning(f"⚠️ Skipping {res_type}: Insufficient class diversity")
            results["resistance_types"][res_type] = {"status": "SKIPPED", "reason": "Insufficient class diversity"}
            continue

        type_results = {}

        logger.info(f"  Running TRUE SAE for {res_type}...")
        true_sae_metrics = evaluate_classifier(X_true_sae, y_labels, f"TRUE_SAE_{res_type}")
        type_results["TRUE_SAE"] = true_sae_metrics

        logger.info(f"  Running Proxy SAE for {res_type}...")
        proxy_sae_metrics = evaluate_classifier(X_proxy_sae, y_labels, f"PROXY_SAE_{res_type}")
        type_results["PROXY_SAE"] = proxy_sae_metrics
        
        results["resistance_types"][res_type] = type_results

    # Generate Report
    with open(REPORT_FILE, 'w') as f:
        f.write(f"# MMRF CoMMpass TRUE SAE Validation Report\n\n")
        f.write(f"**Run Date:** {results['metadata']['run_date']}\n")
        f.write(f"**Number of Patients:** {results['metadata']['num_patients']}\n\n")
        f.write("## Resistance Type Summary\n\n")
        for res_type, type_results in results["resistance_types"].items():
            f.write(f"### {res_type.replace('_', ' ')}\n\n")
            if type_results.get("status") == "SKIPPED":
                f.write(f"**Status:** SKIPPED - {type_results['reason']}\n\n")
                continue

            for model_type, metrics in type_results.items():
                f.write(f"#### {model_type.replace('_', ' ')}\n\n")
                f.write(f"- **Status:** {metrics['status']}\n")
                f.write(f"- **Mean AUROC:** {metrics['mean_auroc']:.3f} (± {metrics['std_auroc']:.3f})\n")
                f.write(f"- **Mean AUPRC:** {metrics['mean_auprc']:.3f} (± {metrics['std_auprc']:.3f})\n\n")

    with open(RESULTS_JSON_FILE, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"✅ Validation results saved to {RESULTS_JSON_FILE}")

if __name__ == "__main__":
    run_validation()
