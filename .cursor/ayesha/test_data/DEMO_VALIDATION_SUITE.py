#!/usr/bin/env python3
"""
⚔️ AYESHA DEMO VALIDATION SUITE ⚔️

Automated testing suite to validate complete demo workflow end-to-end.

Usage:
    cd crispr-assistant-main
    python .cursor/ayesha/test_data/DEMO_VALIDATION_SUITE.py

Requirements:
- Backend running on http://127.0.0.1:8000
- Frontend running on http://localhost:5173
"""

import requests
import json
import time
from pathlib import Path
from typing import Dict, Any, List
import sys

# API Base URL
API_BASE = "http://127.0.0.1:8000"

# Test data directory
TEST_DATA_DIR = Path(__file__).parent

# ANSI color codes for output
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

def print_header(text: str):
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}\n")

def print_success(text: str):
    print(f"{Colors.GREEN}✅ {text}{Colors.END}")

def print_error(text: str):
    print(f"{Colors.RED}❌ {text}{Colors.END}")

def print_info(text: str):
    print(f"{Colors.YELLOW}ℹ️  {text}{Colors.END}")

def print_result(test_name: str, passed: bool, details: str = ""):
    status = "PASS" if passed else "FAIL"
    color = Colors.GREEN if passed else Colors.RED
    print(f"{color}{status}{Colors.END} - {test_name}")
    if details:
        print(f"      {details}")

# ============================================================================
# TEST 1: BACKEND HEALTH CHECK
# ============================================================================

def test_backend_health() -> bool:
    """Test 1: Verify backend is running and healthy"""
    print_header("TEST 1: Backend Health Check")
    
    try:
        response = requests.get(f"{API_BASE}/healthz", timeout=5)
        
        if response.status_code == 200:
            data = response.json()
            print_success(f"Backend healthy: {data}")
            return True
        else:
            print_error(f"Backend unhealthy: {response.status_code}")
            return False
            
    except requests.exceptions.ConnectionError:
        print_error("Backend not reachable - is it running?")
        print_info("Start backend: cd oncology-backend-minimal && venv/bin/python -m uvicorn api.main:app --reload")
        return False
    except Exception as e:
        print_error(f"Health check failed: {e}")
        return False

# ============================================================================
# TEST 2: QUICK INTAKE (Level 0)
# ============================================================================

def test_quick_intake() -> Dict[str, Any]:
    """Test 2: Quick Intake generates Level 0 TumorContext"""
    print_header("TEST 2: Quick Intake (Level 0)")
    
    # Load test data
    with open(TEST_DATA_DIR / "ayesha_level0_intake.json", 'r') as f:
        intake_data = json.load(f)
    
    print_info(f"Sending Quick Intake request for {intake_data['tumor_type']}")
    
    try:
        response = requests.post(
            f"{API_BASE}/api/tumor/quick_intake",
            json=intake_data,
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            tumor_context = result.get("tumor_context", {})
            
            # Validate Level 0 outputs
            checks = []
            
            # Check TMB estimate
            tmb = tumor_context.get("tmb")
            checks.append(("TMB estimated", tmb is not None and tmb > 0))
            print_info(f"TMB: {tmb}")
            
            # Check HRD estimate
            hrd = tumor_context.get("hrd_score")
            checks.append(("HRD estimated (platinum proxy)", hrd is not None and hrd > 0))
            print_info(f"HRD Score: {hrd}")
            
            # Check MSI null (no inference)
            msi = tumor_context.get("msi_status")
            checks.append(("MSI null (no inference)", msi is None or msi == "null"))
            print_info(f"MSI Status: {msi}")
            
            # Check completeness score
            completeness = tumor_context.get("completeness_score", 0.0)
            checks.append(("Completeness < 0.5 (Level 0)", completeness < 0.5))
            print_info(f"Completeness: {completeness}")
            
            # Check provenance
            provenance = result.get("provenance", {})
            checks.append(("No report mode flag", provenance.get("no_report_mode") == True))
            checks.append(("Disease priors used", provenance.get("disease_priors_used") == True))
            
            all_passed = all(check[1] for check in checks)
            
            for check_name, passed in checks:
                print_result(check_name, passed)
            
            if all_passed:
                print_success("Quick Intake: ALL CHECKS PASSED")
            else:
                print_error("Quick Intake: SOME CHECKS FAILED")
            
            return result if all_passed else None
            
        else:
            print_error(f"Quick Intake failed: {response.status_code}")
            print_error(f"Response: {response.text}")
            return None
            
    except Exception as e:
        print_error(f"Quick Intake error: {e}")
        return None

# ============================================================================
# TEST 3: EFFICACY PREDICTION (Level 0 with PARP Penalty)
# ============================================================================

def test_efficacy_level0(tumor_context: Dict[str, Any]) -> Dict[str, Any]:
    """Test 3: Efficacy prediction with Level 0 data (PARP penalty expected)"""
    print_header("TEST 3: Efficacy Prediction (Level 0 - PARP Penalty)")
    
    # Build efficacy request
    efficacy_request = {
        "mutations": [
            {
                "gene": "TP53",
                "hgvs_p": "R248W",
                "chrom": "17",
                "pos": 7676154,
                "ref": "C",
                "alt": "T"
            }
        ],
        "model_id": "evo2_1b",
        "germline_status": "negative",
        "tumor_context": tumor_context,
        "options": {
            "adaptive": True,
            "ensemble": True
        }
    }
    
    print_info("Running efficacy prediction with germline_status='negative' and Level 0 tumor context")
    
    try:
        response = requests.post(
            f"{API_BASE}/api/efficacy/predict",
            json=efficacy_request,
            timeout=120
        )
        
        if response.status_code == 200:
            result = response.json()
            drugs = result.get("drugs", [])
            
            # Find Olaparib (PARP inhibitor)
            olaparib = next((d for d in drugs if "olaparib" in d.get("name", "").lower()), None)
            
            checks = []
            
            if olaparib:
                # Check PARP penalty applied
                provenance = olaparib.get("sporadic_gates_provenance", {})
                gates_applied = provenance.get("gates_applied", [])
                
                checks.append(("Olaparib found in results", True))
                checks.append(("PARP penalty gate applied", "PARP_HRD_LOW" in gates_applied or "PARP_UNKNOWN_HRD" in gates_applied))
                
                efficacy = olaparib.get("efficacy_score", 0)
                confidence = olaparib.get("confidence", 0)
                
                checks.append(("Efficacy reduced (<0.5)", efficacy < 0.5))
                checks.append(("Confidence capped (≤0.4)", confidence <= 0.4))
                
                print_info(f"Olaparib Efficacy: {efficacy:.3f}")
                print_info(f"Olaparib Confidence: {confidence:.3f}")
                print_info(f"Gates Applied: {gates_applied}")
                
                # Check rationale
                rationale = provenance.get("rationale", [])
                parp_rationale = next((r for r in rationale if "PARP" in r.get("gate", "")), None)
                if parp_rationale:
                    print_info(f"PARP Gate Reason: {parp_rationale.get('reason', 'N/A')}")
                
            else:
                print_error("Olaparib not found in results!")
                checks.append(("Olaparib found in results", False))
            
            all_passed = all(check[1] for check in checks)
            
            for check_name, passed in checks:
                print_result(check_name, passed)
            
            if all_passed:
                print_success("Efficacy (L0): ALL CHECKS PASSED - PARP PENALTY WORKING")
            else:
                print_error("Efficacy (L0): SOME CHECKS FAILED")
            
            return result if all_passed else None
            
        else:
            print_error(f"Efficacy prediction failed: {response.status_code}")
            print_error(f"Response: {response.text}")
            return None
            
    except Exception as e:
        print_error(f"Efficacy prediction error: {e}")
        return None

# ============================================================================
# TEST 4: NGS INGESTION (Level 2)
# ============================================================================

def test_ngs_ingestion() -> Dict[str, Any]:
    """Test 4: Ingest tumor NGS report and generate Level 2 TumorContext"""
    print_header("TEST 4: NGS Report Ingestion (Level 2)")
    
    # Load mock NGS report
    with open(TEST_DATA_DIR / "ayesha_tumor_ngs.json", 'r') as f:
        ngs_data = json.load(f)
    
    print_info("Ingesting Foundation Medicine CDx report")
    
    try:
        response = requests.post(
            f"{API_BASE}/api/tumor/ingest_ngs",
            json={
                "report_json": ngs_data,
                "report_source": "Foundation Medicine"
            },
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            tumor_context = result.get("tumor_context", {})
            
            checks = []
            
            # Validate measured values
            tmb = tumor_context.get("tmb")
            checks.append(("TMB measured (6.8)", tmb == 6.8))
            print_info(f"TMB: {tmb}")
            
            hrd = tumor_context.get("hrd_score")
            checks.append(("HRD measured (58)", hrd == 58.0))
            checks.append(("HRD ≥42 (HRD-HIGH)", hrd >= 42))
            print_info(f"HRD Score: {hrd} (HRD-HIGH ✅)")
            
            msi = tumor_context.get("msi_status")
            checks.append(("MSI measured (MSS)", msi == "MSS"))
            print_info(f"MSI Status: {msi}")
            
            # Check BRCA1 biallelic loss
            mutations = tumor_context.get("somatic_mutations", [])
            brca1 = next((m for m in mutations if m.get("gene") == "BRCA1"), None)
            if brca1:
                checks.append(("BRCA1 mutation detected", True))
                checks.append(("BRCA1 LOH detected", brca1.get("loh") == True))
                checks.append(("BRCA1 biallelic loss flagged", brca1.get("biallelic_loss") == True))
                print_info(f"BRCA1: {brca1.get('hgvs_p')} (frameshift + LOH = biallelic loss)")
            
            # Check completeness
            completeness = tumor_context.get("completeness_score", 0.0)
            checks.append(("Completeness ≥0.7 (Level 2)", completeness >= 0.7))
            print_info(f"Completeness: {completeness}")
            
            # Check provenance
            provenance = result.get("provenance", {})
            checks.append(("Source: Foundation Medicine", provenance.get("source") == "Foundation Medicine"))
            checks.append(("Report hash present", "report_hash" in provenance))
            
            all_passed = all(check[1] for check in checks)
            
            for check_name, passed in checks:
                print_result(check_name, passed)
            
            if all_passed:
                print_success("NGS Ingestion: ALL CHECKS PASSED")
            else:
                print_error("NGS Ingestion: SOME CHECKS FAILED")
            
            return result if all_passed else None
            
        else:
            print_error(f"NGS ingestion failed: {response.status_code}")
            print_error(f"Response: {response.text}")
            return None
            
    except Exception as e:
        print_error(f"NGS ingestion error: {e}")
        return None

# ============================================================================
# TEST 5: EFFICACY PREDICTION (Level 2 with PARP Rescue)
# ============================================================================

def test_efficacy_level2(tumor_context: Dict[str, Any]) -> Dict[str, Any]:
    """Test 5: Efficacy prediction with Level 2 data (PARP rescue expected)"""
    print_header("TEST 5: Efficacy Prediction (Level 2 - PARP RESCUE)")
    
    # Build efficacy request with Level 2 data
    efficacy_request = {
        "mutations": [
            {
                "gene": "TP53",
                "hgvs_p": "R248W",
                "chrom": "17",
                "pos": 7676154,
                "ref": "C",
                "alt": "T"
            },
            {
                "gene": "BRCA1",
                "hgvs_p": "Q1756fs",
                "chrom": "17",
                "pos": 43063373,
                "ref": "G",
                "alt": "GC"
            }
        ],
        "model_id": "evo2_1b",
        "germline_status": "negative",
        "tumor_context": tumor_context,
        "options": {
            "adaptive": True,
            "ensemble": True
        }
    }
    
    print_info("Running efficacy prediction with HRD=58 and BRCA1 biallelic loss")
    
    try:
        response = requests.post(
            f"{API_BASE}/api/efficacy/predict",
            json=efficacy_request,
            timeout=120
        )
        
        if response.status_code == 200:
            result = response.json()
            drugs = result.get("drugs", [])
            
            # Find Olaparib
            olaparib = next((d for d in drugs if "olaparib" in d.get("name", "").lower()), None)
            
            checks = []
            
            if olaparib:
                provenance = olaparib.get("sporadic_gates_provenance", {})
                gates_applied = provenance.get("gates_applied", [])
                
                checks.append(("Olaparib found in results", True))
                checks.append(("PARP rescue gate applied", "PARP_HRD_RESCUE" in gates_applied))
                
                efficacy = olaparib.get("efficacy_score", 0)
                confidence = olaparib.get("confidence", 0)
                
                checks.append(("Efficacy RESCUED (≥0.7)", efficacy >= 0.7))
                checks.append(("Confidence HIGH (≥0.7)", confidence >= 0.7))
                checks.append(("No confidence cap (Level 2)", confidence > 0.4))
                
                print_info(f"Olaparib Efficacy: {efficacy:.3f} (vs ~0.32 in Level 0)")
                print_info(f"Olaparib Confidence: {confidence:.3f} (vs 0.4 cap in Level 0)")
                print_info(f"Gates Applied: {gates_applied}")
                
                # Calculate improvement
                improvement = ((efficacy - 0.32) / 0.32) * 100
                print_info(f"Efficacy Improvement: +{improvement:.1f}% (Level 0 → Level 2)")
                
                # Check rationale
                rationale = provenance.get("rationale", [])
                parp_rationale = next((r for r in rationale if "PARP" in r.get("gate", "")), None)
                if parp_rationale:
                    print_info(f"PARP Gate Reason: {parp_rationale.get('reason', 'N/A')}")
                    checks.append(("PARP rescue reason documented", "rescued" in parp_rationale.get('reason', '').lower()))
                
            else:
                print_error("Olaparib not found in results!")
                checks.append(("Olaparib found in results", False))
            
            all_passed = all(check[1] for check in checks)
            
            for check_name, passed in checks:
                print_result(check_name, passed)
            
            if all_passed:
                print_success("Efficacy (L2): ALL CHECKS PASSED - PARP RESCUE WORKING!")
            else:
                print_error("Efficacy (L2): SOME CHECKS FAILED")
            
            return result if all_passed else None
            
        else:
            print_error(f"Efficacy prediction failed: {response.status_code}")
            print_error(f"Response: {response.text}")
            return None
            
    except Exception as e:
        print_error(f"Efficacy prediction error: {e}")
        return None

# ============================================================================
# TEST 6: TMB-HIGH IO BOOST (Bonus Test)
# ============================================================================

def test_io_boost() -> bool:
    """Test 6: Verify immunotherapy boost for TMB-high"""
    print_header("TEST 6: Immunotherapy Boost (TMB-High)")
    
    # Create tumor context with TMB ≥20
    tumor_context = {
        "tmb": 22.0,
        "msi_status": "MSS",
        "hrd_score": 30.0,
        "completeness_score": 0.8
    }
    
    efficacy_request = {
        "mutations": [
            {
                "gene": "TP53",
                "hgvs_p": "R248W",
                "chrom": "17",
                "pos": 7676154,
                "ref": "C",
                "alt": "T"
            }
        ],
        "model_id": "evo2_1b",
        "germline_status": "negative",
        "tumor_context": tumor_context,
        "options": {"adaptive": True}
    }
    
    print_info("Running efficacy with TMB=22 (TMB-HIGH)")
    
    try:
        response = requests.post(
            f"{API_BASE}/api/efficacy/predict",
            json=efficacy_request,
            timeout=120
        )
        
        if response.status_code == 200:
            result = response.json()
            drugs = result.get("drugs", [])
            
            # Find checkpoint inhibitor (Pembrolizumab, Nivolumab, etc.)
            io_drug = next((d for d in drugs 
                           if any(name in d.get("name", "").lower() 
                                 for name in ["pembrolizumab", "nivolumab", "checkpoint"])), 
                          None)
            
            checks = []
            
            if io_drug:
                provenance = io_drug.get("sporadic_gates_provenance", {})
                gates_applied = provenance.get("gates_applied", [])
                
                checks.append(("Checkpoint inhibitor found", True))
                checks.append(("IO boost gate applied", "IO_TMB_HIGH_BOOST" in gates_applied))
                
                print_info(f"Drug: {io_drug.get('name')}")
                print_info(f"Efficacy: {io_drug.get('efficacy_score'):.3f}")
                print_info(f"Gates Applied: {gates_applied}")
                
                # Check boost was applied
                rationale = provenance.get("rationale", [])
                io_rationale = next((r for r in rationale if "IO" in r.get("gate", "")), None)
                if io_rationale:
                    boost = io_rationale.get("boost", 1.0)
                    checks.append(("Boost factor ≥1.3", boost >= 1.3))
                    print_info(f"IO Boost Factor: {boost}x")
                    print_info(f"IO Boost Reason: {io_rationale.get('reason', 'N/A')}")
            else:
                print_error("Checkpoint inhibitor not found in results!")
                checks.append(("Checkpoint inhibitor found", False))
            
            all_passed = all(check[1] for check in checks)
            
            for check_name, passed in checks:
                print_result(check_name, passed)
            
            if all_passed:
                print_success("IO Boost: ALL CHECKS PASSED")
            else:
                print_error("IO Boost: SOME CHECKS FAILED")
            
            return all_passed
            
        else:
            print_error(f"Efficacy prediction failed: {response.status_code}")
            return False
            
    except Exception as e:
        print_error(f"Efficacy prediction error: {e}")
        return False

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print(f"""
{Colors.BOLD}{Colors.BLUE}
⚔️ ═══════════════════════════════════════════════════════════════ ⚔️
   AYESHA DEMO VALIDATION SUITE
   Complete E2E Testing for Sporadic Cancer Workflow
⚔️ ═══════════════════════════════════════════════════════════════ ⚔️
{Colors.END}
    """)
    
    results = {
        "test_1_health": False,
        "test_2_quick_intake": False,
        "test_3_efficacy_l0": False,
        "test_4_ngs_ingestion": False,
        "test_5_efficacy_l2": False,
        "test_6_io_boost": False
    }
    
    # Test 1: Health check
    results["test_1_health"] = test_backend_health()
    if not results["test_1_health"]:
        print_error("\n❌ ABORTING: Backend not healthy!\n")
        sys.exit(1)
    
    time.sleep(1)
    
    # Test 2: Quick Intake (Level 0)
    quick_intake_result = test_quick_intake()
    results["test_2_quick_intake"] = quick_intake_result is not None
    
    if not results["test_2_quick_intake"]:
        print_error("\n❌ Quick Intake failed - cannot continue!\n")
        sys.exit(1)
    
    time.sleep(1)
    
    # Test 3: Efficacy with Level 0 (PARP penalty)
    tumor_context_l0 = quick_intake_result.get("tumor_context", {})
    efficacy_l0_result = test_efficacy_level0(tumor_context_l0)
    results["test_3_efficacy_l0"] = efficacy_l0_result is not None
    
    time.sleep(1)
    
    # Test 4: NGS Ingestion (Level 2)
    ngs_result = test_ngs_ingestion()
    results["test_4_ngs_ingestion"] = ngs_result is not None
    
    if not results["test_4_ngs_ingestion"]:
        print_error("\n❌ NGS Ingestion failed - cannot test Level 2!\n")
    else:
        time.sleep(1)
        
        # Test 5: Efficacy with Level 2 (PARP rescue)
        tumor_context_l2 = ngs_result.get("tumor_context", {})
        efficacy_l2_result = test_efficacy_level2(tumor_context_l2)
        results["test_5_efficacy_l2"] = efficacy_l2_result is not None
    
    time.sleep(1)
    
    # Test 6: IO Boost (TMB-high)
    results["test_6_io_boost"] = test_io_boost()
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    
    print_header("VALIDATION SUMMARY")
    
    total_tests = len(results)
    passed_tests = sum(1 for passed in results.values() if passed)
    pass_rate = (passed_tests / total_tests) * 100
    
    print(f"\nTests Run: {total_tests}")
    print(f"Tests Passed: {passed_tests}")
    print(f"Tests Failed: {total_tests - passed_tests}")
    print(f"Pass Rate: {pass_rate:.1f}%\n")
    
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        color = Colors.GREEN if passed else Colors.RED
        print(f"{color}{status}{Colors.END} - {test_name}")
    
    if pass_rate == 100:
        print(f"\n{Colors.BOLD}{Colors.GREEN}")
        print("⚔️ ═══════════════════════════════════════════════════════════════ ⚔️")
        print("   🎯 ALL TESTS PASSED - DEMO READY FOR AYESHA! 🎯")
        print("⚔️ ═══════════════════════════════════════════════════════════════ ⚔️")
        print(f"{Colors.END}\n")
        return 0
    else:
        print(f"\n{Colors.BOLD}{Colors.RED}")
        print("⚔️ ═══════════════════════════════════════════════════════════════ ⚔️")
        print("   ⚠️  SOME TESTS FAILED - FIX BEFORE DEMO! ⚠️")
        print("⚔️ ═══════════════════════════════════════════════════════════════ ⚔️")
        print(f"{Colors.END}\n")
        return 1

if __name__ == "__main__":
    sys.exit(main())



