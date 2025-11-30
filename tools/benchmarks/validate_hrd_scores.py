#!/usr/bin/env python3
"""
Validate HRD scores extracted from cBioPortal GISTIC data.

Option 2: Validate HRD Score Quality
- Check score distribution (min/max/mean/median, histogram)
- Check correlation with known HRD markers (BRCA1/2 mutations, platinum response)
- Compare to literature expectations (~50% HRD-high)
"""

import json
import os
import argparse
from typing import Dict, List, Optional
from collections import defaultdict
import statistics

def load_hrd_scores(hrd_file: str) -> Dict[str, Dict]:
    """Load HRD scores from JSON file"""
    with open(hrd_file) as f:
        hrd_data = json.load(f)
    
    # Create lookup: patient_id -> HRD metrics
    hrd_lookup = {}
    for record in hrd_data:
        patient_id = record.get("patient_id", "")
        if patient_id:
            hrd_lookup[patient_id] = {
                "hrd_score": record.get("hrd_score", 0),
                "loh": record.get("loh", 0),
                "lst": record.get("lst", 0),
                "tai": record.get("tai", 0)
            }
    
    return hrd_lookup

def load_patient_data(patient_file: str) -> List[Dict]:
    """Load patient data with mutations and platinum response"""
    with open(patient_file) as f:
        return json.load(f)

def check_brca_mutations(patient: Dict) -> bool:
    """Check if patient has BRCA1 or BRCA2 mutations"""
    mutations = patient.get("somatic_mutations", [])
    for mut in mutations:
        gene = mut.get("gene", "").upper()
        if gene in ["BRCA1", "BRCA2"]:
            return True
    return False

def get_platinum_response(patient: Dict) -> Optional[bool]:
    """Get platinum response status (if available)"""
    # Check various field names
    platinum = patient.get("outcome_platinum")
    if platinum is not None:
        return bool(platinum)
    
    platinum = patient.get("platinum_response")
    if platinum is not None:
        return bool(platinum)
    
    return None

def calculate_statistics(scores: List[float]) -> Dict[str, float]:
    """Calculate summary statistics"""
    if not scores:
        return {}
    
    return {
        "count": len(scores),
        "mean": statistics.mean(scores),
        "median": statistics.median(scores),
        "std": statistics.stdev(scores) if len(scores) > 1 else 0.0,
        "min": min(scores),
        "max": max(scores),
        "q25": statistics.quantiles(scores, n=4)[0] if len(scores) > 1 else scores[0],
        "q75": statistics.quantiles(scores, n=4)[2] if len(scores) > 1 else scores[0]
    }

def validate_hrd_scores(hrd_file: str, patient_file: Optional[str] = None) -> Dict:
    """Main validation function"""
    print("🔍 Loading HRD scores...")
    hrd_lookup = load_hrd_scores(hrd_file)
    
    print(f"✅ Loaded HRD scores for {len(hrd_lookup)} patients")
    
    # Extract all HRD scores
    all_hrd_scores = [r["hrd_score"] for r in hrd_lookup.values()]
    all_loh_scores = [r["loh"] for r in hrd_lookup.values()]
    all_lst_scores = [r["lst"] for r in hrd_lookup.values()]
    all_tai_scores = [r["tai"] for r in hrd_lookup.values()]
    
    # 1. Score Distribution
    print("\n📊 HRD Score Distribution:")
    hrd_stats = calculate_statistics(all_hrd_scores)
    print(f"  Count: {hrd_stats['count']}")
    print(f"  Mean: {hrd_stats['mean']:.2f}")
    print(f"  Median: {hrd_stats['median']:.2f}")
    print(f"  Std Dev: {hrd_stats['std']:.2f}")
    print(f"  Min: {hrd_stats['min']}")
    print(f"  Max: {hrd_stats['max']}")
    print(f"  Q25: {hrd_stats['q25']:.2f}")
    print(f"  Q75: {hrd_stats['q75']:.2f}")
    
    # HRD-High threshold (≥42)
    hrd_high_count = sum(1 for s in all_hrd_scores if s >= 42)
    hrd_high_pct = 100 * hrd_high_count / len(all_hrd_scores) if all_hrd_scores else 0
    print(f"\n  HRD-High (≥42): {hrd_high_count}/{len(all_hrd_scores)} ({hrd_high_pct:.1f}%)")
    print(f"  Expected (literature): ~50%")
    
    # Component scores
    print("\n📊 Component Score Statistics:")
    print(f"  LOH: mean={statistics.mean(all_loh_scores):.2f}, "
          f"min={min(all_loh_scores)}, max={max(all_loh_scores)}")
    print(f"  LST: mean={statistics.mean(all_lst_scores):.2f}, "
          f"min={min(all_lst_scores)}, max={max(all_lst_scores)}")
    print(f"  TAI: mean={statistics.mean(all_tai_scores):.2f}, "
          f"min={min(all_tai_scores)}, max={max(all_tai_scores)}")
    
    # 2. Correlation with BRCA mutations
    brca_correlation = None
    if patient_file and os.path.exists(patient_file):
        print("\n🔬 Checking correlation with BRCA mutations...")
        patient_data = load_patient_data(patient_file)
        
        brca_hrd_scores = []
        non_brca_hrd_scores = []
        
        for patient in patient_data:
            patient_id = patient.get("patient_id", "")
            if patient_id not in hrd_lookup:
                continue
            
            hrd_score = hrd_lookup[patient_id]["hrd_score"]
            has_brca = check_brca_mutations(patient)
            
            if has_brca:
                brca_hrd_scores.append(hrd_score)
            else:
                non_brca_hrd_scores.append(hrd_score)
        
        if brca_hrd_scores and non_brca_hrd_scores:
            brca_mean = statistics.mean(brca_hrd_scores)
            non_brca_mean = statistics.mean(non_brca_hrd_scores)
            
            print(f"  BRCA1/2 mutated: {len(brca_hrd_scores)} patients, "
                  f"mean HRD={brca_mean:.2f}")
            print(f"  Non-BRCA: {len(non_brca_hrd_scores)} patients, "
                  f"mean HRD={non_brca_mean:.2f}")
            print(f"  Difference: {brca_mean - non_brca_mean:.2f} "
                  f"({'✅ Higher' if brca_mean > non_brca_mean else '❌ Lower'})")
            
            brca_correlation = {
                "brca_count": len(brca_hrd_scores),
                "brca_mean": brca_mean,
                "non_brca_count": len(non_brca_hrd_scores),
                "non_brca_mean": non_brca_mean,
                "difference": brca_mean - non_brca_mean
            }
        else:
            print("  ⚠️ Insufficient data for BRCA correlation")
    
    # 3. Correlation with platinum response
    platinum_correlation = None
    if patient_file and os.path.exists(patient_file):
        print("\n💊 Checking correlation with platinum response...")
        patient_data = load_patient_data(patient_file)
        
        responder_hrd_scores = []
        non_responder_hrd_scores = []
        
        for patient in patient_data:
            patient_id = patient.get("patient_id", "")
            if patient_id not in hrd_lookup:
                continue
            
            hrd_score = hrd_lookup[patient_id]["hrd_score"]
            platinum_response = get_platinum_response(patient)
            
            if platinum_response is True:
                responder_hrd_scores.append(hrd_score)
            elif platinum_response is False:
                non_responder_hrd_scores.append(hrd_score)
        
        if responder_hrd_scores and non_responder_hrd_scores:
            responder_mean = statistics.mean(responder_hrd_scores)
            non_responder_mean = statistics.mean(non_responder_hrd_scores)
            
            print(f"  Platinum responders: {len(responder_hrd_scores)} patients, "
                  f"mean HRD={responder_mean:.2f}")
            print(f"  Non-responders: {len(non_responder_hrd_scores)} patients, "
                  f"mean HRD={non_responder_mean:.2f}")
            print(f"  Difference: {responder_mean - non_responder_mean:.2f} "
                  f"({'✅ Higher' if responder_mean > non_responder_mean else '❌ Lower'})")
            
            platinum_correlation = {
                "responder_count": len(responder_hrd_scores),
                "responder_mean": responder_mean,
                "non_responder_count": len(non_responder_hrd_scores),
                "non_responder_mean": non_responder_mean,
                "difference": responder_mean - non_responder_mean
            }
        else:
            print("  ⚠️ Insufficient data for platinum correlation")
    
    # Compile validation report
    validation_report = {
        "hrd_score_distribution": hrd_stats,
        "hrd_high_percentage": hrd_high_pct,
        "hrd_high_count": hrd_high_count,
        "expected_hrd_high": 50.0,
        "component_scores": {
            "loh": calculate_statistics(all_loh_scores),
            "lst": calculate_statistics(all_lst_scores),
            "tai": calculate_statistics(all_tai_scores)
        },
        "brca_correlation": brca_correlation,
        "platinum_correlation": platinum_correlation,
        "validation_status": {
            "distribution_check": "✅ PASS" if 0 <= hrd_stats['min'] <= hrd_stats['max'] <= 200 else "❌ FAIL",
            "hrd_high_check": "✅ PASS" if 30 <= hrd_high_pct <= 70 else "⚠️ WARN",
            "brca_correlation_check": "✅ PASS" if brca_correlation and brca_correlation['difference'] > 0 else "⚠️ WARN",
            "platinum_correlation_check": "✅ PASS" if platinum_correlation and platinum_correlation['difference'] > 0 else "⚠️ WARN"
        }
    }
    
    return validation_report

def main():
    parser = argparse.ArgumentParser(description="Validate HRD scores extracted from cBioPortal")
    parser.add_argument("--hrd-file", default="tools/benchmarks/data/full_hrd_scores.json",
                       help="HRD scores JSON file")
    parser.add_argument("--patient-file", 
                       default="tools/benchmarks/data/tcga_ov_patients_with_hrd.json",
                       help="Patient data file (for correlation analysis)")
    parser.add_argument("--output", 
                       default=".cursor/ayesha/sae_documentation/HRD_SCORE_VALIDATION.md",
                       help="Output validation report")
    
    args = parser.parse_args()
    
    # Run validation
    report = validate_hrd_scores(args.hrd_file, args.patient_file)
    
    # Generate markdown report
    print("\n📝 Generating validation report...")
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    with open(args.output, "w") as f:
        f.write("# HRD Score Validation Report\n\n")
        f.write("**Generated by:** `validate_hrd_scores.py`\n")
        f.write("**Date:** " + str(__import__('datetime').datetime.now()) + "\n\n")
        
        f.write("## 📊 Score Distribution\n\n")
        stats = report["hrd_score_distribution"]
        f.write(f"- **Count:** {stats['count']}\n")
        f.write(f"- **Mean:** {stats['mean']:.2f}\n")
        f.write(f"- **Median:** {stats['median']:.2f}\n")
        f.write(f"- **Std Dev:** {stats['std']:.2f}\n")
        f.write(f"- **Range:** {stats['min']} - {stats['max']}\n")
        f.write(f"- **Q25:** {stats['q25']:.2f}\n")
        f.write(f"- **Q75:** {stats['q75']:.2f}\n\n")
        
        f.write(f"### HRD-High Threshold (≥42)\n\n")
        f.write(f"- **Observed:** {report['hrd_high_count']}/{stats['count']} ({report['hrd_high_percentage']:.1f}%)\n")
        f.write(f"- **Expected (literature):** ~50%\n")
        f.write(f"- **Status:** {report['validation_status']['hrd_high_check']}\n\n")
        
        f.write("## 🔬 Correlation Analysis\n\n")
        
        if report["brca_correlation"]:
            brca = report["brca_correlation"]
            f.write("### BRCA1/2 Mutations\n\n")
            f.write(f"- **BRCA mutated:** {brca['brca_count']} patients, mean HRD={brca['brca_mean']:.2f}\n")
            f.write(f"- **Non-BRCA:** {brca['non_brca_count']} patients, mean HRD={brca['non_brca_mean']:.2f}\n")
            f.write(f"- **Difference:** {brca['difference']:.2f}\n")
            f.write(f"- **Status:** {report['validation_status']['brca_correlation_check']}\n\n")
        else:
            f.write("### BRCA1/2 Mutations\n\n")
            f.write("⚠️ Insufficient data for correlation analysis\n\n")
        
        if report["platinum_correlation"]:
            plat = report["platinum_correlation"]
            f.write("### Platinum Response\n\n")
            f.write(f"- **Responders:** {plat['responder_count']} patients, mean HRD={plat['responder_mean']:.2f}\n")
            f.write(f"- **Non-responders:** {plat['non_responder_count']} patients, mean HRD={plat['non_responder_mean']:.2f}\n")
            f.write(f"- **Difference:** {plat['difference']:.2f}\n")
            f.write(f"- **Status:** {report['validation_status']['platinum_correlation_check']}\n\n")
        else:
            f.write("### Platinum Response\n\n")
            f.write("⚠️ Insufficient data for correlation analysis\n\n")
        
        f.write("## ✅ Validation Summary\n\n")
        f.write("| Check | Status |\n")
        f.write("|-------|--------|\n")
        f.write(f"| Distribution (0-200 range) | {report['validation_status']['distribution_check']} |\n")
        f.write(f"| HRD-High (~50%) | {report['validation_status']['hrd_high_check']} |\n")
        f.write(f"| BRCA correlation | {report['validation_status']['brca_correlation_check']} |\n")
        f.write(f"| Platinum correlation | {report['validation_status']['platinum_correlation_check']} |\n\n")
        
        f.write("## 📝 Notes\n\n")
        f.write("- HRD scores calculated from GISTIC discrete copy-number data (gene-level proxy)\n")
        f.write("- LOH, LST, TAI components calculated using simplified gene-level approximations\n")
        f.write("- Literature expectation: ~50% of TCGA-OV patients are HRD-high (≥42)\n")
        f.write("- BRCA1/2 mutations should correlate with higher HRD scores\n")
        f.write("- Platinum responders should have higher HRD scores on average\n")
    
    print(f"✅ Validation report saved to: {args.output}")
    
    # Print summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    for check, status in report["validation_status"].items():
        print(f"  {check}: {status}")

if __name__ == "__main__":
    main()







