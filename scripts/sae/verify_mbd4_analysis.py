#!/usr/bin/env python3
"""
Unified Verification Script for MBD4+TP53 Analysis

Runs all verification checks and generates comprehensive report.
Integrates all verification modules:
- Variant classification (ClinVar, COSMIC, Evo2)
- Pathway mapping (KEGG, Reactome, formula, TCGA)
- Mechanism vector (structure, pathway mapping)
- Consistency checks

Reference: MBD4_TP53_VERIFICATION_LAYER_IMPLEMENTATION_PLAN.md Task 4.1
"""

import json
import sys
import subprocess
from pathlib import Path
from typing import Dict, Any
from datetime import datetime


def run_verification_script(script_path: str, analysis_path: str, api_base: str = None) -> Dict[str, Any]:
    """
    Run a verification script and return its results.
    
    Args:
        script_path: Path to verification script
        analysis_path: Path to analysis JSON file
        api_base: Optional API base URL
    
    Returns:
        Dict with verification results or error
    """
    try:
        cmd = ["python3", script_path, analysis_path]
        if api_base:
            cmd.append(api_base)
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode == 0:
            # Try to load JSON output if script generates it
            output_path = analysis_path.replace(".json", f"_{Path(script_path).stem}.json")
            if Path(output_path).exists():
                with open(output_path, 'r') as f:
                    return json.load(f)
            return {"status": "success", "stdout": result.stdout}
        else:
            return {
                "status": "error",
                "returncode": result.returncode,
                "stderr": result.stderr,
                "stdout": result.stdout
            }
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "error": "Script execution timed out"}
    except Exception as e:
        return {"status": "error", "error": str(e)}


def run_all_verifications(analysis_result_path: str, api_base: str = "http://127.0.0.1:8000") -> Dict[str, Any]:
    """
    Run all verification checks on analysis results.
    
    Args:
        analysis_result_path: Path to analysis JSON output
        api_base: Base URL for API (default: localhost)
    
    Returns:
        Comprehensive verification report
    """
    script_dir = Path(__file__).parent
    
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "analysis_file": analysis_result_path,
        "checks": {},
        "overall_score": {
            "total_checks": 0,
            "passed": 0,
            "failed": 0,
            "pass_rate": 0.0
        }
    }
    
    # Phase 1: Deterministic Verification
    print("Running Phase 1: Deterministic Verification...")
    
    # Task 1.1: Variant Classification
    variant_script = script_dir / "verify_variant_classification.py"
    if variant_script.exists():
        print("  - Variant Classification Verification...")
        variant_result = run_verification_script(str(variant_script), analysis_result_path, api_base)
        verification_report["checks"]["variant_classification"] = variant_result
    
    # Task 1.2: Pathway Mapping
    pathway_script = script_dir / "verify_pathway_mapping.py"
    if pathway_script.exists():
        print("  - Pathway Mapping Verification...")
        pathway_result = run_verification_script(str(pathway_script), analysis_result_path)
        verification_report["checks"]["pathway_mapping"] = pathway_result
    
    # Task 1.3: Functional Annotation
    functional_script = script_dir / "verify_functional_annotation.py"
    if functional_script.exists():
        print("  - Functional Annotation Verification...")
        functional_result = run_verification_script(str(functional_script), analysis_result_path)
        verification_report["checks"]["functional_annotation"] = functional_result
    
    # Task 1.4: Eligibility & IO
    eligibility_script = script_dir / "verify_eligibility_io.py"
    if eligibility_script.exists():
        print("  - Eligibility & IO Verification...")
        eligibility_result = run_verification_script(str(eligibility_script), analysis_result_path)
        verification_report["checks"]["eligibility_io"] = eligibility_result
    
    # Phase 2: Formula & Consistency Verification
    print("Running Phase 2: Formula & Consistency Verification...")
    
    # Task 2.2: Mechanism Vector
    mechanism_script = script_dir / "verify_mechanism_vector.py"
    if mechanism_script.exists():
        print("  - Mechanism Vector Verification...")
        mechanism_result = run_verification_script(str(mechanism_script), analysis_result_path)
        verification_report["checks"]["mechanism_vector"] = mechanism_result
    
    # Task 2.3: Consistency Checks
    consistency_script = script_dir / "verify_consistency.py"
    if consistency_script.exists():
        print("  - Consistency Checks...")
        consistency_result = run_verification_script(str(consistency_script), analysis_result_path)
        verification_report["checks"]["consistency"] = consistency_result
    
    # Compute overall verification score
    all_checks = []
    
    # Extract verified flags from all check results
    for category, check_result in verification_report["checks"].items():
        if isinstance(check_result, dict):
            # Try to load the actual JSON report file if it exists
            # Try different naming patterns
            possible_paths = [
                analysis_result_path.replace(".json", f"_{category}_verification.json"),
                analysis_result_path.replace(".json", f"_{Path(script_dir / f'verify_{category}.py').stem}.json"),
            ]
            
            # Load actual verification report if available
            loaded_report = None
            for path_str in possible_paths:
                path = Path(path_str)
                if path.exists():
                    try:
                        with open(path, 'r') as f:
                            loaded_report = json.load(f)
                            break
                    except:
                        pass
            
            # Use loaded report if available, otherwise use check_result
            report_to_use = loaded_report if loaded_report else check_result
            
            # Check if it's a verification report with overall_score
            if "overall_score" in report_to_use:
                overall = report_to_use["overall_score"]
                total = overall.get("total_checks", 0)
                passed = overall.get("passed", 0)
                if total > 0:
                    all_checks.extend([True] * passed + [False] * (total - passed))
            # Check if it's a verification report with variants_verified
            elif "variants_verified" in report_to_use:
                for variant_check in report_to_use["variants_verified"]:
                    for check_name, check_data in variant_check.get("checks", {}).items():
                        if isinstance(check_data, dict) and "verified" in check_data:
                            all_checks.append(check_data["verified"])
            # Check if it's a verification report with pathway_mappings
            elif "pathway_mappings" in report_to_use:
                for pathway_check in report_to_use["pathway_mappings"]:
                    for check_name, check_data in pathway_check.get("checks", {}).items():
                        if isinstance(check_data, dict) and "verified" in check_data:
                            all_checks.append(check_data["verified"])
                # Also check formula_check and tcga_validation
                if "formula_check" in report_to_use:
                    formula = report_to_use["formula_check"]
                    if isinstance(formula, dict) and "verified" in formula:
                        all_checks.append(formula["verified"])
                if "tcga_validation" in report_to_use:
                    tcga = report_to_use["tcga_validation"]
                    if isinstance(tcga, dict) and "verified" in tcga:
                        all_checks.append(tcga["verified"])
            # Check if it's a verification report with functional_annotations
            elif "functional_annotations" in report_to_use:
                for func_check in report_to_use["functional_annotations"]:
                    for check_name, check_data in func_check.get("checks", {}).items():
                        if isinstance(check_data, dict) and "verified" in check_data:
                            all_checks.append(check_data["verified"])
            # Check if it's a verification report with structure_check and pathway_mapping_check
            elif "structure_check" in report_to_use or "pathway_mapping_check" in report_to_use:
                if "structure_check" in report_to_use:
                    structure = report_to_use["structure_check"]
                    if isinstance(structure, dict) and "verified" in structure:
                        all_checks.append(structure["verified"])
                if "pathway_mapping_check" in report_to_use:
                    mapping = report_to_use["pathway_mapping_check"]
                    if isinstance(mapping, dict) and "verified" in mapping:
                        all_checks.append(mapping["verified"])
            # Check if it's a verification report with pathway_consistency and variant_consistency
            elif "pathway_consistency" in report_to_use or "variant_consistency" in report_to_use:
                if "pathway_consistency" in report_to_use:
                    pathway = report_to_use["pathway_consistency"]
                    if isinstance(pathway, dict) and "verified" in pathway:
                        all_checks.append(pathway["verified"])
                if "variant_consistency" in report_to_use:
                    variant = report_to_use["variant_consistency"]
                    if isinstance(variant, dict) and "verified" in variant:
                        all_checks.append(variant["verified"])
            # Check if it's a verification report with io_eligibility and drug_recommendations
            elif "io_eligibility" in report_to_use or "drug_recommendations" in report_to_use:
                if "io_eligibility" in report_to_use:
                    io = report_to_use["io_eligibility"]
                    if isinstance(io, dict) and "verified" in io:
                        all_checks.append(io["verified"])
                if "drug_recommendations" in report_to_use:
                    for drug_check in report_to_use["drug_recommendations"]:
                        verification = drug_check.get("verification", {})
                        if isinstance(verification, dict) and "verified" in verification:
                            all_checks.append(verification["verified"])
    
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
        print("Usage: python3 verify_mbd4_analysis.py <analysis_result.json> [api_base]")
        print("Example: python3 verify_mbd4_analysis.py results/ayesha_analysis/ayesha_mbd4_tp53_analysis_20251127_013200.json")
        sys.exit(1)
    
    analysis_path = sys.argv[1]
    api_base = sys.argv[2] if len(sys.argv) > 2 else "http://127.0.0.1:8000"
    
    if not Path(analysis_path).exists():
        print(f"Error: Analysis file not found: {analysis_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("MBD4+TP53 ANALYSIS VERIFICATION")
    print("=" * 60)
    print(f"Analysis file: {analysis_path}")
    print(f"API base: {api_base}")
    print()
    
    report = run_all_verifications(analysis_path, api_base)
    
    # Print summary
    print()
    print("=" * 60)
    print("VERIFICATION SUMMARY")
    print("=" * 60)
    
    overall = report["overall_score"]
    print(f"Overall Pass Rate: {overall['pass_rate']:.1%}")
    print(f"Total Checks: {overall['total_checks']}")
    print(f"Passed: {overall['passed']}")
    print(f"Failed: {overall['failed']}")
    print()
    
    # Print category summaries
    script_dir = Path(__file__).parent
    for category, check_result in report["checks"].items():
        if isinstance(check_result, dict):
            # Try to load the actual JSON report file to get accurate status
            possible_paths = [
                analysis_path.replace(".json", f"_{category}_verification.json"),
                analysis_path.replace(".json", f"_{Path(script_dir / f'verify_{category}.py').stem}.json"),
            ]
            
            loaded_report = None
            for path_str in possible_paths:
                path = Path(path_str)
                if path.exists():
                    try:
                        with open(path, 'r') as f:
                            loaded_report = json.load(f)
                            break
                    except:
                        pass
            
            # Use loaded report if available
            if loaded_report and "overall_score" in loaded_report:
                cat_overall = loaded_report["overall_score"]
                pass_rate = cat_overall.get("pass_rate", 0)
                status = "✅" if pass_rate >= 0.8 else "⚠️" if pass_rate >= 0.5 else "❌"
                print(f"{status} {category.replace('_', ' ').title()}: {pass_rate:.1%} ({cat_overall.get('passed', 0)}/{cat_overall.get('total_checks', 0)})")
            elif "overall_score" in check_result:
                cat_overall = check_result["overall_score"]
                pass_rate = cat_overall.get("pass_rate", 0)
                status = "✅" if pass_rate >= 0.8 else "⚠️" if pass_rate >= 0.5 else "❌"
                print(f"{status} {category.replace('_', ' ').title()}: {pass_rate:.1%} ({cat_overall.get('passed', 0)}/{cat_overall.get('total_checks', 0)})")
            elif "status" in check_result:
                # Check if stdout indicates success
                stdout = check_result.get("stdout", "")
                if "Overall Pass Rate: 100.0%" in stdout or "Overall Pass Rate: 100%" in stdout:
                    status = "✅"
                elif check_result["status"] == "success":
                    status = "✅"
                else:
                    status = "❌"
                print(f"{status} {category.replace('_', ' ').title()}: {check_result.get('status', 'unknown')}")
    
    print("=" * 60)
    
    # Save verification report
    output_path = analysis_path.replace(".json", "_verification.json")
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nVerification report saved to: {output_path}")
    
    # Exit with error code if any checks failed
    if overall["pass_rate"] < 1.0:
        print(f"\n⚠️  Some verification checks failed. Review report for details.")
        sys.exit(1)
    else:
        print(f"\n✅ All verification checks passed!")


if __name__ == "__main__":
    main()

