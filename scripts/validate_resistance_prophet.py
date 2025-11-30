"""
Resistance Prophet Validation Script - Phase 1 Retrospective

Based on Manager Approved Plan (Jan 14, 2025):
- Dataset: Jr's 161 patients (genomics + platinum response)
- Signals: DNA repair + pathway escape ONLY (NO CA-125 for retrospective)
- Outcome: Binary (sensitive vs resistant+refractory)
- Metrics: AUROC >=0.70 (primary), sensitivity, specificity, PPV/NPV (secondary)

Manager Decisions Applied:
Q1: Use ALL 161 patients; subgroup analysis for Stage IIIC+IV
Q2: Binary classification (sensitive vs resistant+refractory)
Q3: Proceed WITHOUT CA-125 (signals = DNA repair + pathway escape)
Q4: Primary AUROC >=0.70; secondary sens/spec/PPV/NPV + calibration
Q13: Go/No-Go: AUROC >=0.70 AND sensitivity >=0.75 AND specificity >=0.70

Author: Zo
Date: January 14, 2025
For: Ayesha's life - validate the Prophet
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import logging
from sklearn.metrics import (
    roc_auc_score, 
    roc_curve,
    precision_recall_curve,
    average_precision_score,
    confusion_matrix,
    classification_report
)
from sklearn.calibration import calibration_curve
import matplotlib.pyplot as plt
import seaborn as sns

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Resistance Prophet thresholds (Manager Q9)
HIGH_RISK_THRESHOLD = 0.70
MEDIUM_RISK_THRESHOLD = 0.50


def load_jr_dataset(file_path: str) -> pd.DataFrame:
    """
    Load Jr's 161 patients with platinum response data.
    
    Expected fields:
    - patient_id
    - mutations (list)
    - platinum_response (sensitive/resistant/refractory)
    - stage
    - OS data
    """
    logger.info(f"Loading Jr's dataset from {file_path}...")
    
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    # Handle different data structures
    if isinstance(data, dict) and 'patients' in data:
        patients = data['patients']
    elif isinstance(data, list):
        patients = data
    else:
        raise ValueError(f"Unknown data structure: {type(data)}")
    
    df = pd.DataFrame(patients)
    logger.info(f"Loaded {len(df)} patients")
    
    # Manager Q1: Report subgroup sizes
    if 'stage' in df.columns:
        stage_counts = df['stage'].value_counts()
        logger.info(f"Stage distribution:\n{stage_counts}")
    
    # Manager Q2: Binary outcome (sensitive vs resistant+refractory)
    if 'platinum_response' in df.columns:
        response_counts = df['platinum_response'].value_counts()
        logger.info(f"Platinum response distribution:\n{response_counts}")
        
        # Create binary label
        df['resistant_binary'] = df['platinum_response'].apply(
            lambda x: 0 if x == 'sensitive' else 1
        )
        logger.info(f"Binary labels: {df['resistant_binary'].value_counts().to_dict()}")
    
    return df


def compute_dna_repair_capacity(mutations: List[Dict]) -> float:
    """
    Compute DNA repair capacity from mutations.
    
    Uses Manager-approved formula (from P0 Fix #1):
    0.6 * pathway_ddr + 0.2 * hrr_essentiality + 0.2 * exon_disruption
    
    Inverted logic (Manager's fix):
    - Has DDR mutation → pathway_ddr = 0.30 (LOW capacity)
    - No DDR mutation → pathway_ddr = 0.90 (HIGH capacity)
    """
    # DDR genes from metastasis rules
    DDR_GENES = {
        "BRCA1", "BRCA2", "RAD51", "RAD51C", "RAD51D", "RAD54L",
        "PALB2", "BARD1", "BRIP1", "ATM", "ATR", "CHEK1", "CHEK2"
    }
    
    # Find DDR mutations
    ddr_mutations = [m for m in mutations if m.get('gene') in DDR_GENES]
    
    # Manager-approved inverted logic (P0 Fix #1)
    pathway_ddr = 0.30 if len(ddr_mutations) > 0 else 0.90
    
    # Placeholder for hrr_essentiality and exon_disruption
    # In real implementation, these come from SAE feature service
    hrr_essentiality = 0.60 if len(ddr_mutations) > 0 else 0.40
    exon_disruption = 0.70 if len(ddr_mutations) > 0 else 0.30
    
    # Manager-approved formula
    dna_repair_capacity = (
        0.6 * pathway_ddr + 
        0.2 * hrr_essentiality + 
        0.2 * exon_disruption
    )
    
    return dna_repair_capacity


def compute_pathway_escape(mutations: List[Dict]) -> float:
    """
    Compute pathway escape signal from mechanism vector shift.
    
    Logic:
    - Detect activation in bypass pathways (MAPK, PI3K)
    - Higher escape score = higher resistance risk
    """
    # Bypass pathway genes
    MAPK_GENES = {"KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3"}
    PI3K_GENES = {"PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "MTOR", "PTEN"}
    
    # Find bypass pathway mutations
    mapk_mutations = [m for m in mutations if m.get('gene') in MAPK_GENES]
    pi3k_mutations = [m for m in mutations if m.get('gene') in PI3K_GENES]
    
    # Compute escape score
    # More bypass pathway mutations = higher escape risk
    escape_score = min((len(mapk_mutations) + len(pi3k_mutations)) / 3.0, 1.0)
    
    return escape_score


def compute_resistance_probability(dna_repair: float, pathway_escape: float) -> float:
    """
    Compute resistance probability from DNA repair + pathway escape signals.
    
    Phase 1: NO CA-125 signal (2-of-2 detection)
    
    Logic:
    - Both signals high → HIGH probability
    - One signal high → MEDIUM probability
    - Both signals low → LOW probability
    """
    # Invert DNA repair capacity to resistance risk
    # LOW capacity (0.30) → HIGH risk (0.70)
    # HIGH capacity (0.90) → LOW risk (0.10)
    dna_repair_risk = 1.0 - dna_repair
    
    # Average the two signals
    probability = (dna_repair_risk + pathway_escape) / 2.0
    
    return min(probability, 1.0)


def stratify_risk(probability: float, signal_count: int) -> str:
    """
    Stratify risk level (Manager Q9).
    
    HIGH: probability >=0.70 AND >=2 signals
    MEDIUM: 0.50-0.69 OR exactly 1 signal
    LOW: <0.50
    """
    if probability >= HIGH_RISK_THRESHOLD and signal_count >= 2:
        return "HIGH"
    elif probability >= MEDIUM_RISK_THRESHOLD or signal_count == 1:
        return "MEDIUM"
    else:
        return "LOW"


def validate_resistance_prophet(df: pd.DataFrame) -> Dict:
    """
    Validate Resistance Prophet on Jr's 161 patients.
    
    Manager Q13: Go/No-Go criteria:
    - AUROC >=0.70 overall
    - Sensitivity >=0.75
    - Specificity >=0.70
    """
    logger.info("=== RESISTANCE PROPHET VALIDATION (Phase 1: NO CA-125) ===")
    
    # Compute predictions for all patients
    predictions = []
    
    for idx, row in df.iterrows():
        mutations = row.get('mutations', [])
        
        # Signal 1: DNA repair capacity
        dna_repair = compute_dna_repair_capacity(mutations)
        
        # Signal 2: Pathway escape
        pathway_escape = compute_pathway_escape(mutations)
        
        # Overall probability
        probability = compute_resistance_probability(dna_repair, pathway_escape)
        
        # Count signals (both DNA repair and pathway escape)
        signal_count = 2  # Both signals always computed in Phase 1
        
        # Risk stratification
        risk_level = stratify_risk(probability, signal_count)
        
        predictions.append({
            'patient_id': row.get('patient_id'),
            'dna_repair_capacity': dna_repair,
            'pathway_escape_score': pathway_escape,
            'resistance_probability': probability,
            'risk_level': risk_level,
            'true_response': row.get('platinum_response'),
            'true_resistant': row.get('resistant_binary')
        })
    
    # Create predictions DataFrame
    pred_df = pd.DataFrame(predictions)
    
    # Compute metrics (Manager Q4)
    y_true = pred_df['true_resistant'].values
    y_pred_proba = pred_df['resistance_probability'].values
    y_pred_binary = (y_pred_proba >= 0.50).astype(int)
    
    # Primary metric: AUROC (Manager Q13: target >=0.70)
    auroc = roc_auc_score(y_true, y_pred_proba)
    
    # Secondary metrics: Sensitivity, Specificity, PPV, NPV
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_binary).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
    
    # Average Precision (AUPRC)
    auprc = average_precision_score(y_true, y_pred_proba)
    
    # Manager Q13: Go/No-Go assessment
    go_no_go = {
        "auroc_pass": auroc >= 0.70,
        "sensitivity_pass": sensitivity >= 0.75,
        "specificity_pass": specificity >= 0.70,
        "overall_pass": auroc >= 0.70 and sensitivity >= 0.75 and specificity >= 0.70
    }
    
    metrics = {
        "auroc": auroc,
        "auprc": auprc,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": ppv,
        "npv": npv,
        "accuracy": (tp + tn) / len(y_true),
        "total_patients": len(pred_df),
        "true_resistant": int(y_true.sum()),
        "true_sensitive": int(len(y_true) - y_true.sum()),
        "predicted_resistant": int(y_pred_binary.sum()),
        "predicted_sensitive": int(len(y_pred_binary) - y_pred_binary.sum()),
        "confusion_matrix": {
            "tp": int(tp), "tn": int(tn), "fp": int(fp), "fn": int(fn)
        },
        "go_no_go": go_no_go
    }
    
    logger.info("=== PRIMARY METRICS (Manager Q13: Go/No-Go) ===")
    logger.info(f"AUROC: {auroc:.3f} (target >=0.70) {'✅ PASS' if go_no_go['auroc_pass'] else '❌ FAIL'}")
    logger.info(f"Sensitivity: {sensitivity:.3f} (target >=0.75) {'✅ PASS' if go_no_go['sensitivity_pass'] else '❌ FAIL'}")
    logger.info(f"Specificity: {specificity:.3f} (target >=0.70) {'✅ PASS' if go_no_go['specificity_pass'] else '❌ FAIL'}")
    logger.info(f"OVERALL: {'✅ GO - DEPLOY' if go_no_go['overall_pass'] else '⚠️ NO-GO - RUO ONLY'}")
    
    return {
        "metrics": metrics,
        "predictions": pred_df,
        "y_true": y_true,
        "y_pred_proba": y_pred_proba
    }


def generate_validation_report(validation_results: Dict, output_dir: Path):
    """
    Generate comprehensive validation report.
    
    Manager Q1: Subgroup analysis for Stage IIIC+IV
    Manager Q4: Include calibration curve
    """
    metrics = validation_results["metrics"]
    pred_df = validation_results["predictions"]
    y_true = validation_results["y_true"]
    y_pred_proba = validation_results["y_pred_proba"]
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Performance metrics report
    report = []
    report.append("=" * 80)
    report.append("RESISTANCE PROPHET VALIDATION REPORT (Phase 1: NO CA-125)")
    report.append("=" * 80)
    report.append(f"Date: {pd.Timestamp.now()}")
    report.append(f"Dataset: Jr's 161 patients (TCGA-OV + GDC platinum response)")
    report.append(f"Signals Used: DNA Repair + Pathway Escape (2-of-2)")
    report.append(f"Manager Policy: MANAGER_ANSWERS_TO_RESISTANCE_PROPHET_QUESTIONS.md")
    report.append("")
    
    report.append("=" * 80)
    report.append("PRIMARY METRICS (Manager Q13: Go/No-Go Criteria)")
    report.append("=" * 80)
    report.append(f"AUROC: {metrics['auroc']:.3f} (target >=0.70) {'✅ PASS' if metrics['go_no_go']['auroc_pass'] else '❌ FAIL'}")
    report.append(f"Sensitivity: {metrics['sensitivity']:.3f} (target >=0.75) {'✅ PASS' if metrics['go_no_go']['sensitivity_pass'] else '❌ FAIL'}")
    report.append(f"Specificity: {metrics['specificity']:.3f} (target >=0.70) {'✅ PASS' if metrics['go_no_go']['specificity_pass'] else '❌ FAIL'}")
    report.append("")
    report.append(f"OVERALL VERDICT: {'✅ GO - PRODUCTION DEPLOY' if metrics['go_no_go']['overall_pass'] else '⚠️ NO-GO - RUO ONLY WITH DISCLAIMERS'}")
    report.append("")
    
    report.append("=" * 80)
    report.append("SECONDARY METRICS (Manager Q4)")
    report.append("=" * 80)
    report.append(f"AUPRC: {metrics['auprc']:.3f}")
    report.append(f"PPV: {metrics['ppv']:.3f}")
    report.append(f"NPV: {metrics['npv']:.3f}")
    report.append(f"Accuracy: {metrics['accuracy']:.3f}")
    report.append("")
    
    report.append("=" * 80)
    report.append("CONFUSION MATRIX")
    report.append("=" * 80)
    cm = metrics['confusion_matrix']
    report.append(f"True Positives (TP): {cm['tp']} (correctly predicted resistant)")
    report.append(f"True Negatives (TN): {cm['tn']} (correctly predicted sensitive)")
    report.append(f"False Positives (FP): {cm['fp']} (predicted resistant but was sensitive)")
    report.append(f"False Negatives (FN): {cm['fn']} (predicted sensitive but was resistant)")
    report.append("")
    
    report.append("=" * 80)
    report.append("DATASET SUMMARY")
    report.append("=" * 80)
    report.append(f"Total Patients: {metrics['total_patients']}")
    report.append(f"True Resistant: {metrics['true_resistant']} ({metrics['true_resistant']/metrics['total_patients']:.1%})")
    report.append(f"True Sensitive: {metrics['true_sensitive']} ({metrics['true_sensitive']/metrics['total_patients']:.1%})")
    report.append("")
    
    # Save report
    report_path = output_dir / "RESISTANCE_PROPHET_VALIDATION_REPORT.txt"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    
    logger.info(f"Report saved: {report_path}")
    
    # 2. Save predictions CSV
    pred_csv_path = output_dir / "resistance_predictions.csv"
    pred_df.to_csv(pred_csv_path, index=False)
    logger.info(f"Predictions saved: {pred_csv_path}")
    
    # 3. Generate figures
    _generate_validation_figures(y_true, y_pred_proba, pred_df, output_dir)
    
    return report_path


def _generate_validation_figures(y_true, y_pred_proba, pred_df, output_dir):
    """Generate validation figures (Manager Q4: calibration curve)"""
    
    # Figure 1: ROC Curve
    fpr, tpr, thresholds = roc_curve(y_true, y_pred_proba)
    auroc = roc_auc_score(y_true, y_pred_proba)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, label=f'Resistance Prophet (AUROC={auroc:.3f})', linewidth=2)
    plt.plot([0, 1], [0, 1], 'k--', label='Random (AUROC=0.50)')
    plt.axhline(y=0.75, color='r', linestyle=':', label='Target Sensitivity (0.75)')
    plt.xlabel('False Positive Rate (1 - Specificity)')
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.title('Resistance Prophet ROC Curve\n(Phase 1: DNA Repair + Pathway Escape)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "roc_curve.png", dpi=300)
    plt.close()
    logger.info("ROC curve saved")
    
    # Figure 2: Precision-Recall Curve
    precision, recall, pr_thresholds = precision_recall_curve(y_true, y_pred_proba)
    auprc = average_precision_score(y_true, y_pred_proba)
    
    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, label=f'Resistance Prophet (AUPRC={auprc:.3f})', linewidth=2)
    plt.axhline(y=y_true.mean(), color='k', linestyle='--', label=f'Baseline Precision ({y_true.mean():.3f})')
    plt.xlabel('Recall (Sensitivity)')
    plt.ylabel('Precision (PPV)')
    plt.title('Resistance Prophet Precision-Recall Curve')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "precision_recall_curve.png", dpi=300)
    plt.close()
    logger.info("Precision-Recall curve saved")
    
    # Figure 3: Calibration Curve (Manager Q4)
    prob_true, prob_pred = calibration_curve(y_true, y_pred_proba, n_bins=10)
    
    plt.figure(figsize=(8, 6))
    plt.plot(prob_pred, prob_true, marker='o', label='Resistance Prophet', linewidth=2)
    plt.plot([0, 1], [0, 1], 'k--', label='Perfect Calibration')
    plt.xlabel('Predicted Probability')
    plt.ylabel('True Fraction of Resistant')
    plt.title('Resistance Prophet Calibration Curve')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "calibration_curve.png", dpi=300)
    plt.close()
    logger.info("Calibration curve saved")
    
    # Figure 4: Probability distribution by outcome
    plt.figure(figsize=(10, 6))
    
    sensitive_probs = pred_df[pred_df['true_resistant'] == 0]['resistance_probability']
    resistant_probs = pred_df[pred_df['true_resistant'] == 1]['resistance_probability']
    
    plt.hist(sensitive_probs, bins=20, alpha=0.5, label=f'Sensitive (n={len(sensitive_probs)})', color='green')
    plt.hist(resistant_probs, bins=20, alpha=0.5, label=f'Resistant (n={len(resistant_probs)})', color='red')
    plt.axvline(x=0.50, color='k', linestyle='--', label='MEDIUM threshold (0.50)')
    plt.axvline(x=0.70, color='r', linestyle='--', label='HIGH threshold (0.70)')
    plt.xlabel('Predicted Resistance Probability')
    plt.ylabel('Count')
    plt.title('Resistance Prophet Probability Distribution')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "probability_distribution.png", dpi=300)
    plt.close()
    logger.info("Probability distribution saved")
    
    logger.info("All figures generated")


def run_subgroup_analysis(df: pd.DataFrame) -> Dict:
    """
    Run subgroup analysis for Stage IIIC+IV (Manager Q1).
    """
    logger.info("Running subgroup analysis (Manager Q1)...")
    
    # Filter to Stage IIIC+IV
    if 'stage' in df.columns:
        advanced_df = df[df['stage'].isin(['IIIC', 'IV', 'IVA', 'IVB'])]
        logger.info(f"Advanced stage (IIIC+IV) patients: {len(advanced_df)}")
        
        if len(advanced_df) >= 40:  # Manager Q1: Ensure N>=40
            logger.info("Sufficient N for subgroup analysis - running validation...")
            subgroup_results = validate_resistance_prophet(advanced_df)
            return subgroup_results
        else:
            logger.warning(f"Insufficient N for subgroup ({len(advanced_df)} < 40) - SKIPPING")
            return {"status": "skipped", "reason": "insufficient_n", "n": len(advanced_df)}
    else:
        logger.warning("No 'stage' column - cannot run subgroup analysis")
        return {"status": "skipped", "reason": "no_stage_data"}


def main():
    """Main validation execution"""
    
    # Input file (Jr's data hunt result)
    input_file = Path("tools/benchmarks/tcga_ov_platinum_response_with_genomics.json")
    
    if not input_file.exists():
        logger.error(f"Dataset not found: {input_file}")
        logger.info("Expected Jr's platinum response dataset with 161 patients")
        return
    
    # Load dataset
    df = load_jr_dataset(str(input_file))
    
    # Validate on ALL 161 patients (Manager Q1)
    logger.info("\n=== VALIDATING ON ALL 161 PATIENTS ===")
    all_results = validate_resistance_prophet(df)
    
    # Generate report
    output_dir = Path("results/resistance_prophet_validation")
    generate_validation_report(all_results, output_dir)
    
    # Subgroup analysis (Manager Q1)
    logger.info("\n=== SUBGROUP ANALYSIS (STAGE IIIC+IV) ===")
    subgroup_results = run_subgroup_analysis(df)
    
    if subgroup_results.get("status") != "skipped":
        subgroup_dir = output_dir / "subgroup_advanced_stage"
        generate_validation_report(subgroup_results, subgroup_dir)
    
    logger.info("\n🎯 VALIDATION COMPLETE - Check results/ directory for full report")
    logger.info(f"Manager Q13 Go/No-Go: {all_results['metrics']['go_no_go']}")


if __name__ == "__main__":
    main()

