#!/usr/bin/env python3
"""
⚔️ SAE VALIDATION SCRIPT - TCGA OVARIAN CANCER COHORT ⚔️

Mission: Validate SAE Phase 2 against REAL TCGA data.
Data Source: tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json (89 real patients)
Ground Truth: HRD scores from TCGA (0-100 scale)

Metrics:
- AUROC (Area Under ROC): Discrimination ability
- AUPRC (Area Under Precision-Recall): Performance on imbalanced classes
- Sensitivity/Specificity at HRD ≥42 threshold
- Per-feature correlation with HRD

Output: JSON results + plots + manager-reviewable metrics

Date: January 14, 2025
Status: VALIDATION SCRIPT V1
"""

import json
import sys
import os
from pathlib import Path
from typing import Dict, List, Any, Optional
import requests
from dataclasses import dataclass
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================================
# CONFIGURATION
# ============================================================================

# Backend URL
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")

# Data paths
TCGA_DATA_PATH = Path(__file__).parent.parent / "tools" / "benchmarks" / "hrd_tcga_ov_labeled_sample_use_evo.json"
OUTPUT_DIR = Path(__file__).parent.parent / ".cursor" / "ayesha" / "sae_documentation" / "validation_results"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# HRD Threshold (per Manager's Policy)
HRD_THRESHOLD = 42

# API Timeout
API_TIMEOUT = 120  # seconds


# ============================================================================
# DATA MODELS
# ============================================================================

@dataclass
class ValidationResult:
    """Results for a single patient"""
    patient_id: str
    true_hrd_score: float
    true_hrd_high: bool  # True if HRD ≥42
    
    # Predicted values
    predicted_hrd_proxy: Optional[float]  # From DNA repair capacity
    predicted_ddr_burden: Optional[float]  # From pathway_burden_ddr
    predicted_her2_burden: Optional[float]  # From pathway_burden_her2
    
    # API response status
    api_success: bool
    api_error: Optional[str]
    
    # Raw SAE features
    sae_features: Optional[Dict[str, Any]]


# ============================================================================
# DATA LOADING
# ============================================================================

def load_tcga_data() -> List[Dict[str, Any]]:
    """Load real TCGA data"""
    print(f"📂 Loading TCGA data from: {TCGA_DATA_PATH}")
    
    if not TCGA_DATA_PATH.exists():
        raise FileNotFoundError(f"TCGA data not found: {TCGA_DATA_PATH}")
    
    with open(TCGA_DATA_PATH, 'r') as f:
        data = json.load(f)
    
    print(f"✅ Loaded {len(data)} TCGA samples")
    return data


# ============================================================================
# API INTERACTION
# ============================================================================

def build_request_payload(tcga_sample: Dict[str, Any]) -> Dict[str, Any]:
    """Convert TCGA sample to API request payload"""
    
    # Extract mutations
    mutations = []
    for mut in tcga_sample.get("somatic_mutations", []):
        mutations.append({
            "gene": mut["gene"],
            "hgvs_p": mut.get("protein_change", ""),
            "variant_type": mut.get("variant_type", "SNV"),
            "chrom": mut.get("chrom", ""),
            "pos": mut.get("pos", 0),
            "ref": mut.get("ref", ""),
            "alt": mut.get("alt", "")
        })
    
    # Build tumor context
    tumor_context = {
        "patient_id": tcga_sample["patient_id"],
        "disease": "ovarian_cancer_hgs",
        "stage": "IV",  # TCGA-OV is mostly advanced
        "somatic_mutations": mutations,
        "hrd_score": tcga_sample.get("hrd_score"),  # Ground truth
        "tmb_mut_per_mb": tcga_sample.get("tmb", 0),
        "platinum_sensitive": tcga_sample.get("platinum_sensitive"),
        "completeness": 0.8  # Level 2 (NGS available)
    }
    
    # Build request
    return {
        "germline_status": "negative",  # TCGA is somatic-focused
        "tumor_context": tumor_context,
        "model_id": "evo2_1b"
    }


def call_backend_api(payload: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Call /api/ayesha/complete_care_v2 endpoint"""
    
    endpoint = f"{BACKEND_URL}/api/ayesha/complete_care_v2"
    
    try:
        response = requests.post(
            endpoint,
            json=payload,
            timeout=API_TIMEOUT,
            headers={"Content-Type": "application/json"}
        )
        
        if response.status_code == 200:
            return response.json()
        else:
            print(f"⚠️ API returned status {response.status_code}: {response.text[:200]}")
            return None
            
    except requests.exceptions.Timeout:
        print(f"⏱️ API timeout after {API_TIMEOUT}s")
        return None
    except Exception as e:
        print(f"❌ API error: {e}")
        return None


# ============================================================================
# VALIDATION LOGIC
# ============================================================================

def validate_patient(tcga_sample: Dict[str, Any]) -> ValidationResult:
    """Validate a single patient"""
    
    patient_id = tcga_sample["patient_id"]
    true_hrd = tcga_sample.get("hrd_score", 0)
    true_hrd_high = true_hrd >= HRD_THRESHOLD
    
    print(f"🧬 Validating {patient_id} (HRD={true_hrd}, HRD-High={true_hrd_high})...")
    
    # Build payload
    payload = build_request_payload(tcga_sample)
    
    # Call API
    response = call_backend_api(payload)
    
    if response is None:
        return ValidationResult(
            patient_id=patient_id,
            true_hrd_score=true_hrd,
            true_hrd_high=true_hrd_high,
            predicted_hrd_proxy=None,
            predicted_ddr_burden=None,
            predicted_her2_burden=None,
            api_success=False,
            api_error="API call failed",
            sae_features=None
        )
    
    # Extract SAE features
    sae_features = response.get("sae_features")
    
    if sae_features is None:
        return ValidationResult(
            patient_id=patient_id,
            true_hrd_score=true_hrd,
            true_hrd_high=true_hrd_high,
            predicted_hrd_proxy=None,
            predicted_ddr_burden=None,
            predicted_her2_burden=None,
            api_success=False,
            api_error="sae_features not in response",
            sae_features=None
        )
    
    # Extract predictions
    dna_repair_capacity = sae_features.get("dna_repair_capacity", 0)
    pathway_burden_ddr = sae_features.get("pathway_burden_ddr", 0)
    pathway_burden_her2 = sae_features.get("pathway_burden_her2", 0)
    
    # DNA repair capacity is INVERSELY correlated with HRD
    # High DNA repair → Low HRD
    # Low DNA repair → High HRD
    # So we'll use (1 - dna_repair_capacity) as HRD proxy
    predicted_hrd_proxy = (1.0 - dna_repair_capacity) * 100  # Scale to 0-100
    
    print(f"  ✅ DNA_repair={dna_repair_capacity:.2f}, HRD_proxy={predicted_hrd_proxy:.1f}, DDR_burden={pathway_burden_ddr:.2f}")
    
    return ValidationResult(
        patient_id=patient_id,
        true_hrd_score=true_hrd,
        true_hrd_high=true_hrd_high,
        predicted_hrd_proxy=predicted_hrd_proxy,
        predicted_ddr_burden=pathway_burden_ddr,
        predicted_her2_burden=pathway_burden_her2,
        api_success=True,
        api_error=None,
        sae_features=sae_features
    )


def run_validation(max_patients: Optional[int] = None) -> List[ValidationResult]:
    """Run validation on TCGA cohort"""
    
    print("⚔️ STARTING SAE VALIDATION ⚔️\n")
    
    # Load data
    tcga_data = load_tcga_data()
    
    if max_patients:
        tcga_data = tcga_data[:max_patients]
        print(f"🔬 Testing on first {max_patients} patients\n")
    
    # Validate each patient
    results = []
    for i, sample in enumerate(tcga_data, 1):
        print(f"[{i}/{len(tcga_data)}] ", end="")
        result = validate_patient(sample)
        results.append(result)
        print()
    
    return results


# ============================================================================
# METRICS COMPUTATION
# ============================================================================

def compute_metrics(results: List[ValidationResult]) -> Dict[str, Any]:
    """Compute validation metrics"""
    
    print("\n📊 COMPUTING METRICS...\n")
    
    # Filter successful API calls
    valid_results = [r for r in results if r.api_success and r.predicted_hrd_proxy is not None]
    
    if len(valid_results) == 0:
        print("❌ NO VALID RESULTS TO COMPUTE METRICS")
        return {"error": "No valid results"}
    
    print(f"✅ {len(valid_results)}/{len(results)} patients with valid predictions\n")
    
    # Extract arrays
    y_true = np.array([r.true_hrd_high for r in valid_results])  # Binary: HRD ≥42
    y_score = np.array([r.predicted_hrd_proxy for r in valid_results])  # Continuous: 0-100
    
    # Normalize scores to 0-1 for AUROC/AUPRC
    y_score_norm = y_score / 100.0
    
    # Check if we have both classes
    if len(np.unique(y_true)) < 2:
        print("⚠️ Only one class present in ground truth - cannot compute AUROC/AUPRC")
        auroc = None
        auprc = None
    else:
        # Compute AUROC
        auroc = roc_auc_score(y_true, y_score_norm)
        print(f"🎯 AUROC: {auroc:.3f}")
        
        # Compute AUPRC
        auprc = average_precision_score(y_true, y_score_norm)
        print(f"🎯 AUPRC: {auprc:.3f}")
    
    # Compute confusion matrix at threshold
    y_pred = y_score >= HRD_THRESHOLD
    
    tp = np.sum((y_true == 1) & (y_pred == 1))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    
    print(f"🎯 Sensitivity: {sensitivity:.3f}")
    print(f"🎯 Specificity: {specificity:.3f}")
    print(f"🎯 PPV: {ppv:.3f}")
    print(f"🎯 NPV: {npv:.3f}")
    
    # Correlation analysis
    true_hrd_scores = np.array([r.true_hrd_score for r in valid_results])
    pred_hrd_scores = np.array([r.predicted_hrd_proxy for r in valid_results])
    
    correlation = np.corrcoef(true_hrd_scores, pred_hrd_scores)[0, 1]
    print(f"🎯 Correlation (HRD_true vs HRD_pred): {correlation:.3f}")
    
    # Compile metrics
    metrics = {
        "n_total": len(results),
        "n_valid": len(valid_results),
        "n_hrd_high": int(np.sum(y_true)),
        "n_hrd_low": int(len(y_true) - np.sum(y_true)),
        "auroc": float(auroc) if auroc is not None else None,
        "auprc": float(auprc) if auprc is not None else None,
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "ppv": float(ppv),
        "npv": float(npv),
        "correlation": float(correlation),
        "confusion_matrix": {
            "tp": int(tp),
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn)
        }
    }
    
    return metrics


# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_results(results: List[ValidationResult], metrics: Dict[str, Any]):
    """Generate validation plots"""
    
    print("\n📈 GENERATING PLOTS...\n")
    
    valid_results = [r for r in results if r.api_success and r.predicted_hrd_proxy is not None]
    
    if len(valid_results) == 0:
        print("⚠️ No valid results to plot")
        return
    
    # Prepare data
    y_true = np.array([r.true_hrd_high for r in valid_results])
    y_score = np.array([r.predicted_hrd_proxy for r in valid_results]) / 100.0
    true_hrd_scores = np.array([r.true_hrd_score for r in valid_results])
    pred_hrd_scores = np.array([r.predicted_hrd_proxy for r in valid_results])
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("SAE Validation Results - TCGA Ovarian Cancer Cohort", fontsize=16, fontweight='bold')
    
    # Plot 1: ROC Curve
    if metrics.get("auroc") is not None:
        fpr, tpr, _ = roc_curve(y_true, y_score)
        axes[0, 0].plot(fpr, tpr, linewidth=2, label=f'AUROC = {metrics["auroc"]:.3f}')
        axes[0, 0].plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random')
        axes[0, 0].set_xlabel('False Positive Rate')
        axes[0, 0].set_ylabel('True Positive Rate')
        axes[0, 0].set_title('ROC Curve')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Precision-Recall Curve
    if metrics.get("auprc") is not None:
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        axes[0, 1].plot(recall, precision, linewidth=2, label=f'AUPRC = {metrics["auprc"]:.3f}')
        axes[0, 1].set_xlabel('Recall')
        axes[0, 1].set_ylabel('Precision')
        axes[0, 1].set_title('Precision-Recall Curve')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Scatter plot (True HRD vs Predicted HRD)
    axes[1, 0].scatter(true_hrd_scores, pred_hrd_scores, alpha=0.6, s=50)
    axes[1, 0].plot([0, 100], [0, 100], 'r--', linewidth=2, label='Perfect prediction')
    axes[1, 0].axhline(y=HRD_THRESHOLD, color='orange', linestyle='--', linewidth=1, label='HRD threshold')
    axes[1, 0].axvline(x=HRD_THRESHOLD, color='orange', linestyle='--', linewidth=1)
    axes[1, 0].set_xlabel('True HRD Score')
    axes[1, 0].set_ylabel('Predicted HRD Proxy')
    axes[1, 0].set_title(f'HRD Prediction (Correlation = {metrics["correlation"]:.3f})')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Confusion Matrix
    cm = np.array([
        [metrics["confusion_matrix"]["tn"], metrics["confusion_matrix"]["fp"]],
        [metrics["confusion_matrix"]["fn"], metrics["confusion_matrix"]["tp"]]
    ])
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1, 1],
                xticklabels=['Pred Neg', 'Pred Pos'],
                yticklabels=['True Neg', 'True Pos'])
    axes[1, 1].set_title('Confusion Matrix (HRD ≥42 threshold)')
    axes[1, 1].set_ylabel('True Label')
    axes[1, 1].set_xlabel('Predicted Label')
    
    plt.tight_layout()
    
    # Save plot
    plot_path = OUTPUT_DIR / "validation_plots.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"✅ Plots saved to: {plot_path}")
    
    plt.close()


# ============================================================================
# RESULTS EXPORT
# ============================================================================

def save_results(results: List[ValidationResult], metrics: Dict[str, Any]):
    """Save validation results to JSON"""
    
    print("\n💾 SAVING RESULTS...\n")
    
    # Convert results to dict
    results_dict = [
        {
            "patient_id": r.patient_id,
            "true_hrd_score": r.true_hrd_score,
            "true_hrd_high": r.true_hrd_high,
            "predicted_hrd_proxy": r.predicted_hrd_proxy,
            "predicted_ddr_burden": r.predicted_ddr_burden,
            "predicted_her2_burden": r.predicted_her2_burden,
            "api_success": r.api_success,
            "api_error": r.api_error,
            "sae_features": r.sae_features
        }
        for r in results
    ]
    
    # Compile full output
    output = {
        "validation_date": "2025-01-14",
        "data_source": str(TCGA_DATA_PATH),
        "backend_url": BACKEND_URL,
        "hrd_threshold": HRD_THRESHOLD,
        "metrics": metrics,
        "results": results_dict
    }
    
    # Save JSON
    output_path = OUTPUT_DIR / "validation_results.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✅ Results saved to: {output_path}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main validation script"""
    
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate SAE against TCGA data")
    parser.add_argument("--max-patients", type=int, default=None,
                       help="Maximum number of patients to test (default: all)")
    parser.add_argument("--backend-url", type=str, default="http://localhost:8000",
                       help="Backend URL (default: http://localhost:8000)")
    
    args = parser.parse_args()
    
    global BACKEND_URL
    BACKEND_URL = args.backend_url
    
    try:
        # Run validation
        results = run_validation(max_patients=args.max_patients)
        
        # Compute metrics
        metrics = compute_metrics(results)
        
        # Generate plots
        if "error" not in metrics:
            plot_results(results, metrics)
        
        # Save results
        save_results(results, metrics)
        
        print("\n⚔️ VALIDATION COMPLETE ⚔️\n")
        
        if metrics.get("auroc"):
            print(f"📊 AUROC: {metrics['auroc']:.3f}")
            print(f"📊 AUPRC: {metrics['auprc']:.3f}")
            print(f"📊 Sensitivity: {metrics['sensitivity']:.3f}")
            print(f"📊 Specificity: {metrics['specificity']:.3f}")
            
            # Interpretation
            if metrics['auroc'] >= 0.80:
                print("\n✅ EXCELLENT: Model shows strong discrimination ability!")
            elif metrics['auroc'] >= 0.70:
                print("\n✅ GOOD: Model shows acceptable discrimination!")
            elif metrics['auroc'] >= 0.60:
                print("\n⚠️ MODERATE: Model shows some discrimination, but needs improvement.")
            else:
                print("\n❌ POOR: Model performance is weak - major improvements needed.")
        
    except Exception as e:
        print(f"\n❌ VALIDATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

