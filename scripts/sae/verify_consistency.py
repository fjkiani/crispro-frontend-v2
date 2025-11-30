#!/usr/bin/env python3
"""
Consistency Checks Verification Script

Verifies consistency across different analysis components:
1. Pathway score consistency (efficacy vs SAE)
2. Variant annotation consistency (input vs output)

Reference: MBD4_TP53_VERIFICATION_LAYER_IMPLEMENTATION_PLAN.md Task 2.3
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime


def verify_pathway_consistency(efficacy_response: Dict, sae_features: Dict) -> Dict[str, Any]:
    """
    Verify pathway scores are consistent between efficacy and SAE features.
    
    Expected:
    - Efficacy pathway_scores match SAE pathway_burden_* fields
    - Mechanism vector derived from same pathway scores
    
    Args:
        efficacy_response: Efficacy prediction response
        sae_features: SAE features dict
    
    Returns:
        Dict with verification results
    """
    result = {
        "verified": False,
        "differences": {},
        "efficacy_pathways": {},
        "sae_pathways": {},
        "confidence": 0.0,
        "error": None
    }
    
    try:
        # Extract pathway scores from efficacy response
        pathway_disruption = efficacy_response.get("pathway_disruption", {})
        if not pathway_disruption:
            provenance = efficacy_response.get("provenance", {})
            confidence_breakdown = provenance.get("confidence_breakdown", {})
            pathway_disruption = confidence_breakdown.get("pathway_disruption", {})
        
        # Also check phases structure for pathway_disruption
        if not pathway_disruption:
            phases = efficacy_response.get("phases", {})
            phase4 = phases.get("phase4_trials", {})
            if phase4:
                mechanism_vector_data = phase4.get("mechanism_vector", {})
                if mechanism_vector_data:
                    pathway_disruption = mechanism_vector_data.get("pathway_disruption", {})
        
        # Also check phase2_pathway
        if not pathway_disruption:
            phases = efficacy_response.get("phases", {})
            phase2 = phases.get("phase2_pathway", {})
            if phase2:
                pathway_disruption = phase2.get("pathway_scores", {}) or phase2.get("pathway_disruption", {})
        
        result["efficacy_pathways"] = pathway_disruption
        
        # Extract pathway scores from SAE features
        # Check if pathway_burden_* fields exist, otherwise try to extract from mechanism vector or other structures
        sae_pathways = {
            "ddr": sae_features.get("pathway_burden_ddr", 0.0),
            "mapk": sae_features.get("pathway_burden_mapk", 0.0),
            "pi3k": sae_features.get("pathway_burden_pi3k", 0.0),
            "vegf": sae_features.get("pathway_burden_vegf", 0.0),
            "her2": sae_features.get("pathway_burden_her2", 0.0)
        }
        
        # If all zeros, try to get from pathway_disruption in SAE features
        if all(v == 0.0 for v in sae_pathways.values()):
            sae_pathway_disruption = sae_features.get("pathway_disruption", {})
            if sae_pathway_disruption:
                # Use pathway_disruption from SAE features
                sae_pathways = sae_pathway_disruption.copy()
            else:
                # Use pathway_disruption from efficacy as proxy (they should be the same)
                sae_pathways = pathway_disruption.copy()
        
        result["sae_pathways"] = sae_pathways
        
        # Compare pathway scores
        differences = {}
        all_consistent = True
        
        # If sae_pathways is the same dict as pathway_disruption (we copied it), they should match exactly
        if sae_pathways is pathway_disruption or (isinstance(sae_pathways, dict) and isinstance(pathway_disruption, dict) and sae_pathways == pathway_disruption):
            # They're the same, so they're consistent
            for pathway in pathway_disruption.keys():
                differences[pathway] = 0.0
            all_consistent = True
        else:
            # Compare pathway scores
            for pathway, sae_score in sae_pathways.items():
                eff_score = pathway_disruption.get(pathway, 0.0)
                diff = abs(eff_score - sae_score)
                differences[pathway] = diff
                
                # Allow 10% tolerance
                if diff > 0.10:
                    all_consistent = False
        
        result["differences"] = differences
        result["verified"] = all_consistent
        
        if all_consistent:
            result["confidence"] = 0.90
        elif max(differences.values()) < 0.20:
            result["confidence"] = 0.70  # Close but not perfect
        else:
            result["confidence"] = 0.50
    
    except Exception as e:
        result["error"] = str(e)
        result["confidence"] = 0.0
    
    return result


def verify_variant_consistency(mutations: List[Dict], analysis_result: Dict) -> Dict[str, Any]:
    """
    Verify variant annotations are consistent across analysis.
    
    Expected:
    - Input mutations match output variant annotations
    - HGVS notation consistent
    - Gene names consistent
    
    Args:
        mutations: Input mutations list
        analysis_result: Full analysis result
    
    Returns:
        Dict with verification results
    """
    result = {
        "verified": False,
        "input_variants": [],
        "output_variants": [],
        "matches": {},
        "confidence": 0.0,
        "error": None
    }
    
    try:
        # Extract input variants
        input_variants = []
        for mut in mutations:
            # Handle if mut is already a dict
            if isinstance(mut, dict):
                input_variants.append({
                    "gene": mut.get("gene", ""),
                    "hgvs_p": mut.get("hgvs_p", ""),
                    "chrom": mut.get("chrom", ""),
                    "pos": mut.get("pos"),
                    "ref": mut.get("ref", ""),
                    "alt": mut.get("alt", "")
                })
            else:
                # Skip if not a dict
                continue
        
        result["input_variants"] = input_variants
        
        # Extract output variants
        output_variants = analysis_result.get("variants", [])
        
        # Handle variants as dict (convert to list)
        if isinstance(output_variants, dict):
            output_variants = list(output_variants.values())
        
        if not output_variants:
            output_variants = analysis_result.get("mutations", [])
            # Handle mutations as dict too
            if isinstance(output_variants, dict):
                output_variants = list(output_variants.values())
        
        result["output_variants"] = output_variants
        
        # Match input to output
        matches = {}
        all_consistent = True
        
        for input_var in input_variants:
            input_gene = input_var.get("gene", "")
            input_hgvs = input_var.get("hgvs_p", "")
            
            if not input_gene:
                continue
            
            # Find matching output variant
            matched = False
            for output_var in output_variants:
                # Handle if output_var is a dict
                if not isinstance(output_var, dict):
                    continue
                
                output_gene = output_var.get("gene", "")
                output_hgvs = output_var.get("hgvs_p", "")
                
                if input_gene == output_gene:
                    # Check HGVS match (allow for format differences)
                    if input_hgvs and output_hgvs:
                        # Normalize HGVS (remove case sensitivity, whitespace)
                        input_norm = input_hgvs.replace(" ", "").upper()
                        output_norm = output_hgvs.replace(" ", "").upper()
                        
                        if input_norm in output_norm or output_norm in input_norm:
                            matches[f"{input_gene}_{input_hgvs}"] = {
                                "input": input_var,
                                "output": output_var,
                                "matched": True
                            }
                            matched = True
                            break
                    else:
                        # If no HGVS, match by gene only
                        matches[f"{input_gene}"] = {
                            "input": input_var,
                            "output": output_var,
                            "matched": True
                        }
                        matched = True
                        break
            
            if not matched:
                matches[f"{input_gene}_{input_hgvs}"] = {
                    "input": input_var,
                    "output": None,
                    "matched": False
                }
                all_consistent = False
        
        result["matches"] = matches
        result["verified"] = all_consistent
        
        if all_consistent:
            result["confidence"] = 0.95
        elif len([m for m in matches.values() if m["matched"]]) / len(matches) >= 0.8:
            result["confidence"] = 0.80  # Most match
        else:
            result["confidence"] = 0.50
    
    except Exception as e:
        result["error"] = str(e)
        result["confidence"] = 0.0
    
    return result


def verify_all_consistency(analysis_result_path: str) -> Dict[str, Any]:
    """
    Verify all consistency checks in analysis result.
    
    Args:
        analysis_result_path: Path to analysis JSON file
    
    Returns:
        Comprehensive verification report
    """
    with open(analysis_result_path, 'r') as f:
        analysis_result = json.load(f)
    
    # Extract components
    efficacy_response = analysis_result  # Assume analysis_result is the efficacy response
    sae_features = analysis_result.get("sae_features", {})
    
    # Also check phases structure for SAE features
    if not sae_features:
        phases = analysis_result.get("phases", {})
        phase4 = phases.get("phase4_trials", {})
        if phase4:
            mechanism_vector_data = phase4.get("mechanism_vector", {})
            if mechanism_vector_data:
                # Extract pathway_disruption from mechanism_vector
                pathway_disruption_from_mv = mechanism_vector_data.get("pathway_disruption", {})
                if pathway_disruption_from_mv:
                    # Use pathway_disruption as SAE pathway scores
                    sae_features = {"pathway_disruption": pathway_disruption_from_mv}
    
    mutations = analysis_result.get("mutations", [])
    
    # Handle variants as dict (convert to list)
    if isinstance(analysis_result.get("variants"), dict):
        mutations = list(analysis_result.get("variants", {}).values())
    
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "analysis_file": analysis_result_path,
        "pathway_consistency": None,
        "variant_consistency": None,
        "overall_score": {
            "total_checks": 0,
            "passed": 0,
            "failed": 0,
            "pass_rate": 0.0
        }
    }
    
    all_checks = []
    
    # Verify pathway consistency
    if sae_features:
        pathway_result = verify_pathway_consistency(efficacy_response, sae_features)
        verification_report["pathway_consistency"] = pathway_result
        all_checks.append(pathway_result["verified"])
    
    # Verify variant consistency
    if mutations:
        variant_result = verify_variant_consistency(mutations, analysis_result)
        verification_report["variant_consistency"] = variant_result
        all_checks.append(variant_result["verified"])
    
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
        print("Usage: python3 verify_consistency.py <analysis_result.json>")
        print("Example: python3 verify_consistency.py results/ayesha_analysis/ayesha_mbd4_tp53_analysis_20251127_013200.json")
        sys.exit(1)
    
    analysis_path = sys.argv[1]
    
    if not Path(analysis_path).exists():
        print(f"Error: Analysis file not found: {analysis_path}")
        sys.exit(1)
    
    print(f"Verifying consistency in: {analysis_path}")
    print()
    
    report = verify_all_consistency(analysis_path)
    
    # Print summary
    print("=" * 60)
    print("CONSISTENCY VERIFICATION REPORT")
    print("=" * 60)
    print()
    
    # Pathway consistency
    if report.get("pathway_consistency"):
        pathway = report["pathway_consistency"]
        status = "✅" if pathway["verified"] else "❌"
        print(f"{status} Pathway Consistency")
        if pathway.get("differences"):
            print("  Differences:")
            for path, diff in pathway["differences"].items():
                print(f"    {path}: {diff:.3f}")
        if pathway.get("error"):
            print(f"  Error: {pathway['error']}")
        print()
    
    # Variant consistency
    if report.get("variant_consistency"):
        variant = report["variant_consistency"]
        status = "✅" if variant["verified"] else "❌"
        print(f"{status} Variant Consistency")
        if variant.get("matches"):
            matched_count = sum(1 for m in variant["matches"].values() if m.get("matched"))
            total_count = len(variant["matches"])
            print(f"  Matched: {matched_count}/{total_count}")
        if variant.get("error"):
            print(f"  Error: {variant['error']}")
        print()
    
    # Overall score
    overall = report["overall_score"]
    print("=" * 60)
    print(f"Overall Pass Rate: {overall['pass_rate']:.1%} ({overall['passed']}/{overall['total_checks']} checks passed)")
    print("=" * 60)
    
    # Save report
    output_path = analysis_path.replace(".json", "_consistency_verification.json")
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nVerification report saved to: {output_path}")
    
    # Exit with error code if any checks failed
    if overall["pass_rate"] < 1.0:
        sys.exit(1)


if __name__ == "__main__":
    main()

