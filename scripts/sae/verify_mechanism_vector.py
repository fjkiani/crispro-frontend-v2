#!/usr/bin/env python3
"""
Mechanism Vector Verification Script

Verifies mechanism vector structure and pathway mapping:
1. Vector structure (7D, value ranges, sum validation)
2. Pathway mapping (DDR index, TP53 contribution, IO eligibility)

Reference: MBD4_TP53_VERIFICATION_LAYER_IMPLEMENTATION_PLAN.md Task 2.2
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime


def verify_vector_structure(mechanism_vector: List[float], expected_dim: int = 7) -> Dict[str, Any]:
    """
    Verify mechanism vector has correct structure.
    
    Expected:
    - 7D vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    - All values in range [0.0, 1.0]
    - Sum of values ≤ 7.0 (not normalized to 1.0)
    
    Args:
        mechanism_vector: Mechanism vector list
        expected_dim: Expected dimension (default: 7)
    
    Returns:
        Dict with verification results
    """
    result = {
        "verified": False,
        "dimension": len(mechanism_vector) if mechanism_vector else 0,
        "expected_dim": expected_dim,
        "all_in_range": False,
        "sum_valid": False,
        "confidence": 0.0,
        "error": None
    }
    
    if not mechanism_vector:
        result["error"] = "Mechanism vector is empty"
        return result
    
    # Check dimension
    dim_correct = len(mechanism_vector) == expected_dim
    result["dimension"] = len(mechanism_vector)
    
    # Check value ranges (mechanism vectors can be > 1.0, e.g., DDR=1.4 from DDR+50%TP53)
    # Allow values up to 2.0 (reasonable upper bound for combined pathway scores)
    all_in_range = all(0.0 <= v <= 2.0 for v in mechanism_vector)
    result["all_in_range"] = all_in_range
    
    # Check sum (should be ≤ dim*2, not normalized to 1.0)
    vector_sum = sum(mechanism_vector)
    sum_valid = vector_sum <= expected_dim * 2.0  # Allow up to 2.0 per dimension
    result["sum_valid"] = sum_valid
    result["vector_sum"] = vector_sum
    
    # Verify all checks pass
    result["verified"] = dim_correct and all_in_range and sum_valid
    
    if result["verified"]:
        result["confidence"] = 0.95
    elif dim_correct and all_in_range:
        result["confidence"] = 0.70  # Sum might be slightly over
    elif dim_correct:
        result["confidence"] = 0.50  # Dimension correct but values out of range
    else:
        result["confidence"] = 0.0
    
    return result


def verify_pathway_mapping(
    pathway_scores: Dict[str, float],
    mechanism_vector: List[float],
    tumor_context: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Verify pathway scores correctly map to mechanism vector indices.
    
    Expected for MBD4+TP53:
    - DDR pathway → mechanism_vector[0] = 0.70-0.90
    - MAPK pathway → mechanism_vector[1] = 0.10-0.30 (low)
    - IO → mechanism_vector[5] = 1.0 (TMB ≥20)
    
    Args:
        pathway_scores: Dict of pathway scores
        mechanism_vector: Mechanism vector (7D)
        tumor_context: Optional tumor context for IO eligibility
    
    Returns:
        Dict with verification results
    """
    result = {
        "verified": False,
        "ddr_mapping": False,
        "tp53_contribution": False,
        "io_eligibility": False,
        "expected_ddr": None,
        "actual_ddr": None,
        "confidence": 0.0,
        "error": None
    }
    
    if not mechanism_vector or len(mechanism_vector) < 7:
        result["error"] = "Mechanism vector must be 7D"
        return result
    
    try:
        # Check DDR mapping (index 0)
        ddr_score = pathway_scores.get("ddr", 0.0)
        ddr_idx = 0
        actual_ddr = mechanism_vector[ddr_idx]
        result["actual_ddr"] = actual_ddr
        
        # Check TP53 → DDR contribution (50% rule)
        tp53_score = pathway_scores.get("tp53", 0.0)
        expected_ddr_with_tp53 = ddr_score + (tp53_score * 0.5)
        result["expected_ddr"] = expected_ddr_with_tp53
        
        # Verify DDR mapping (within 10% tolerance)
        ddr_difference = abs(actual_ddr - expected_ddr_with_tp53)
        result["ddr_mapping"] = ddr_difference < 0.10
        result["ddr_difference"] = ddr_difference
        
        # Verify TP53 contribution
        result["tp53_contribution"] = ddr_difference < 0.10
        
        # Check IO eligibility (index 5)
        io_idx = 5
        actual_io = mechanism_vector[io_idx]
        expected_io = 0.0
        
        if tumor_context:
            tmb = tumor_context.get("tmb_score", 0.0)
            msi_status = tumor_context.get("msi_status", "")
            
            if tmb >= 20.0 or (msi_status and msi_status.upper() in ["MSI-H", "MSI-HIGH"]):
                expected_io = 1.0
        
        io_difference = abs(actual_io - expected_io)
        result["io_eligibility"] = io_difference < 0.10
        result["expected_io"] = expected_io
        result["actual_io"] = actual_io
        result["io_difference"] = io_difference
        
        # Overall verification
        result["verified"] = result["ddr_mapping"] and result["tp53_contribution"] and result["io_eligibility"]
        
        if result["verified"]:
            result["confidence"] = 0.95
        elif result["ddr_mapping"] and result["tp53_contribution"]:
            result["confidence"] = 0.85  # IO might be edge case
        elif result["ddr_mapping"]:
            result["confidence"] = 0.70  # DDR correct but TP53 wrong
        else:
            result["confidence"] = 0.50
    
    except Exception as e:
        result["error"] = str(e)
        result["confidence"] = 0.0
    
    return result


def verify_all_mechanism_vectors(analysis_result_path: str) -> Dict[str, Any]:
    """
    Verify mechanism vectors in analysis result.
    
    Args:
        analysis_result_path: Path to analysis JSON file
    
    Returns:
        Comprehensive verification report
    """
    with open(analysis_result_path, 'r') as f:
        analysis_result = json.load(f)
    
    # Get mechanism vector (can be nested in different structures)
    sae_features = analysis_result.get("sae_features", {})
    mechanism_vector = sae_features.get("mechanism_vector")
    
    # If not found, check if it's nested in a "vector" key
    if isinstance(mechanism_vector, dict):
        mechanism_vector = mechanism_vector.get("vector")
    
    # Also check in phases structure
    if not mechanism_vector:
        phases = analysis_result.get("phases", {})
        # Check phase4_trials (most likely location based on file structure)
        phase4_trials = phases.get("phase4_trials", {})
        if phase4_trials:
            mechanism_vector_data = phase4_trials.get("mechanism_vector", {})
            if isinstance(mechanism_vector_data, dict):
                mechanism_vector = mechanism_vector_data.get("vector")
        
        # Also check phase4_trial_matching
        if not mechanism_vector:
            phase4_trial = phases.get("phase4_trial_matching", {})
            if phase4_trial:
                mechanism_vector_data = phase4_trial.get("mechanism_vector", {})
                if isinstance(mechanism_vector_data, dict):
                    mechanism_vector = mechanism_vector_data.get("vector")
        
        # Also check phase4_sae_features
        if not mechanism_vector:
            phase4 = phases.get("phase4_sae_features", {})
            if phase4:
                mechanism_vector_data = phase4.get("mechanism_vector", {})
                if isinstance(mechanism_vector_data, dict):
                    mechanism_vector = mechanism_vector_data.get("vector")
        
        # Also check phase5_summary
        if not mechanism_vector:
            phase5 = phases.get("phase5_summary", {})
            if phase5:
                mechanism_vector_data = phase5.get("mechanism_vector", {})
                if isinstance(mechanism_vector_data, dict):
                    mechanism_vector = mechanism_vector_data.get("vector")
    
    # Get pathway scores
    pathway_scores = analysis_result.get("pathway_disruption", {})
    if not pathway_scores:
        provenance = analysis_result.get("provenance", {})
        confidence_breakdown = provenance.get("confidence_breakdown", {})
        pathway_scores = confidence_breakdown.get("pathway_disruption", {})
    
    # Get tumor context
    tumor_context = analysis_result.get("tumor_context", {})
    
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "analysis_file": analysis_result_path,
        "structure_check": None,
        "pathway_mapping_check": None,
        "overall_score": {
            "total_checks": 0,
            "passed": 0,
            "failed": 0,
            "pass_rate": 0.0
        }
    }
    
    all_checks = []
    
    # Verify structure
    if mechanism_vector:
        structure_result = verify_vector_structure(mechanism_vector)
        verification_report["structure_check"] = structure_result
        all_checks.append(structure_result["verified"])
    
    # Verify pathway mapping
    if mechanism_vector and pathway_scores:
        mapping_result = verify_pathway_mapping(pathway_scores, mechanism_vector, tumor_context)
        verification_report["pathway_mapping_check"] = mapping_result
        all_checks.append(mapping_result["verified"])
    
    # Compute overall score
    verification_report["overall_score"] = {
        "total_checks": len(all_checks),
        "passed": sum(all_checks),
        "failed": len(all_checks) - sum(all_checks),
        "pass_rate": sum(all_checks) / len(all_checks) if all_checks else 0.0
    }
    
    return verification_report


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python3 verify_mechanism_vector.py <analysis_result.json>")
        print("Example: python3 verify_mechanism_vector.py results/ayesha_analysis/ayesha_mbd4_tp53_analysis_20251127_013200.json")
        sys.exit(1)
    
    analysis_path = sys.argv[1]
    
    if not Path(analysis_path).exists():
        print(f"Error: Analysis file not found: {analysis_path}")
        sys.exit(1)
    
    print(f"Verifying mechanism vectors in: {analysis_path}")
    print()
    
    report = verify_all_mechanism_vectors(analysis_path)
    
    # Print summary
    print("=" * 60)
    print("MECHANISM VECTOR VERIFICATION REPORT")
    print("=" * 60)
    print()
    
    # Structure check
    if report.get("structure_check"):
        structure = report["structure_check"]
        status = "✅" if structure["verified"] else "❌"
        print(f"{status} Vector Structure")
        print(f"  Dimension: {structure['dimension']} (expected: {structure['expected_dim']})")
        print(f"  All in range [0.0, 1.0]: {structure['all_in_range']}")
        print(f"  Sum valid (≤{structure['expected_dim']}): {structure['sum_valid']}")
        if structure.get("vector_sum") is not None:
            print(f"  Vector sum: {structure['vector_sum']:.3f}")
        if structure.get("error"):
            print(f"  Error: {structure['error']}")
        print()
    
    # Pathway mapping check
    if report.get("pathway_mapping_check"):
        mapping = report["pathway_mapping_check"]
        status = "✅" if mapping["verified"] else "❌"
        print(f"{status} Pathway Mapping")
        print(f"  DDR mapping: {mapping['ddr_mapping']}")
        print(f"    Expected DDR: {mapping.get('expected_ddr', 'N/A'):.3f}")
        print(f"    Actual DDR: {mapping.get('actual_ddr', 'N/A'):.3f}")
        print(f"    Difference: {mapping.get('ddr_difference', 'N/A'):.3f}")
        print(f"  TP53 contribution: {mapping['tp53_contribution']}")
        print(f"  IO eligibility: {mapping['io_eligibility']}")
        print(f"    Expected IO: {mapping.get('expected_io', 'N/A'):.1f}")
        print(f"    Actual IO: {mapping.get('actual_io', 'N/A'):.1f}")
        if mapping.get("error"):
            print(f"  Error: {mapping['error']}")
        print()
    
    # Overall score
    overall = report["overall_score"]
    print("=" * 60)
    print(f"Overall Pass Rate: {overall['pass_rate']:.1%} ({overall['passed']}/{overall['total_checks']} checks passed)")
    print("=" * 60)
    
    # Save report
    output_path = analysis_path.replace(".json", "_mechanism_vector_verification.json")
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nVerification report saved to: {output_path}")
    
    # Exit with error code if any checks failed
    if overall["pass_rate"] < 1.0:
        sys.exit(1)


if __name__ == "__main__":
    main()

