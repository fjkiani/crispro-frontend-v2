#!/usr/bin/env python3
"""
⚔️ MM RESISTANCE PREDICTION - Work Item 4: Validation Framework
================================================================

Mission: Implement comprehensive validation framework for MM resistance prediction
Agent: Plumber (Implementation)
Date: January 29, 2025
Priority: P0 (Critical)

Objective:
- Define 5 key validation tests for MM resistance prediction
- Simulate patient scenarios with specific mutations/cytogenetics
- Call the /api/resistance/predict endpoint
- Verify expected outcomes (risk level, detected signals, alternative therapies)
- Generate markdown report summarizing test results

Tests:
1. PSMB5 mutation → PI Resistance (RR ≥ 2.0)
2. CRBN mutation → IMiD Resistance (RR ≥ 2.5)
3. del(17p) cytogenetic abnormality → Universal Resistance (HR ≥ 2.0)
4. DIS3 + TP53 Co-occurrence → Ultra-high resistance risk
5. TRUE SAE vs PROXY Comparison (TRUE SAE AUROC > PROXY SAE AUROC) - Placeholder

Success Criteria:
- validate_mm_resistance.py created and executable
- All 5 validation tests implemented (some may be placeholders/skipped if data not ready)
- Validation report generated (markdown output)
- At least 3/5 tests pass (with acceptable reasons for failure, e.g., INSUFFICIENT_DATA)

Guardrails:
- RUO/validation-only (no production impact)
- Transparent reporting of test results and assumptions
"""

import asyncio
import httpx
import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
from dataclasses import dataclass
from loguru import logger

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")

# --- Configuration ---
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")
RESISTANCE_PREDICT_ENDPOINT = f"{BACKEND_URL}/api/resistance/predict"
OUTPUT_DIR = Path("data/validation/mm_cohort/validation_reports")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_FILE = OUTPUT_DIR / "mm_resistance_validation_report.md"
RESULTS_JSON_FILE = OUTPUT_DIR / "mm_resistance_validation_results.json"

# --- Test Definitions ---
@dataclass
class TestCase:
    name: str
    description: str
    mutations: List[Dict[str, Any]]
    cytogenetics: Optional[Dict[str, bool]]
    drug_class: Optional[str]
    treatment_line: int
    prior_therapies: Optional[List[str]]
    expected_risk_level: str
    expected_signals: List[str]
    expected_min_probability: float
    expected_max_probability: float
    expected_min_alternatives: int
    skip_reason: Optional[str] = None

TEST_CASES = [
    TestCase(
        name="Test 1: PSMB5 -> PI Resistance",
        description="Detects high resistance to Proteasome Inhibitors due to PSMB5 mutation.",
        mutations=[{"gene": "PSMB5", "hgvs_p": "p.Ala49Thr", "consequence": "missense_variant"}],
        cytogenetics=None,
        drug_class="proteasome_inhibitor",
        treatment_line=1,
        prior_therapies=None,
        expected_risk_level="HIGH",
        expected_signals=["MM_DRUG_CLASS_RESISTANCE"],
        expected_min_probability=0.70,
        expected_max_probability=0.90,
        expected_min_alternatives=1,
    ),
    TestCase(
        name="Test 2: CRBN -> IMiD Resistance",
        description="Detects very high resistance to IMiDs due to CRBN loss-of-function mutation.",
        mutations=[{"gene": "CRBN", "hgvs_p": "p.Trp400*", "consequence": "stop_gained"}],
        cytogenetics=None,
        drug_class="imid",
        treatment_line=2,
        prior_therapies=["imid"],
        expected_risk_level="HIGH",
        expected_signals=["MM_DRUG_CLASS_RESISTANCE"],
        expected_min_probability=0.80,
        expected_max_probability=0.95,
        expected_min_alternatives=1,
    ),
    TestCase(
        name="Test 3: del(17p) -> Universal Resistance",
        description="Detects high universal resistance due to del(17p) cytogenetic abnormality.",
        mutations=[{"gene": "TP53", "hgvs_p": "p.Arg175His", "consequence": "missense_variant"}],
        cytogenetics={"del_17p": True, "t_4_14": False, "1q_gain": False},
        drug_class="anti_cd38",
        treatment_line=3,
        prior_therapies=["proteasome_inhibitor", "imid"],
        expected_risk_level="HIGH",
        expected_signals=["MM_CYTOGENETICS"],
        expected_min_probability=0.75,
        expected_max_probability=0.95,
        expected_min_alternatives=2,
    ),
    TestCase(
        name="Test 4: DIS3 + TP53 Co-occurrence -> Ultra-high Risk",
        description="Detects ultra-high risk due to co-occurrence of DIS3 and TP53 mutations.",
        mutations=[
            {"gene": "DIS3", "hgvs_p": "p.Lys75del", "consequence": "frameshift_variant"},
            {"gene": "TP53", "hgvs_p": "p.Arg248Gln", "consequence": "missense_variant"}
        ],
        cytogenetics=None,
        drug_class="proteasome_inhibitor",
        treatment_line=1,
        prior_therapies=None,
        expected_risk_level="HIGH",
        expected_signals=["MM_HIGH_RISK_GENE"],
        expected_min_probability=0.70,
        expected_max_probability=0.90,
        expected_min_alternatives=1,
    ),
    TestCase(
        name="Test 5: TRUE SAE vs PROXY Comparison (Placeholder)",
        description="Compares TRUE SAE performance against Proxy SAE baseline (requires WI1 data).",
        mutations=[],
        cytogenetics=None,
        drug_class=None,
        treatment_line=1,
        prior_therapies=None,
        expected_risk_level="LOW",
        expected_signals=[],
        expected_min_probability=0.0,
        expected_max_probability=1.0,
        expected_min_alternatives=0,
        skip_reason="Requires data from Work Item 1 (MMRF TRUE SAE extraction and validation) to be meaningful."
    ),
]

async def run_test_case(test_case: TestCase) -> Dict[str, Any]:
    logger.info(f"\n--- Running {test_case.name} ---")
    test_result = {
        "name": test_case.name,
        "description": test_case.description,
        "status": "FAIL",
        "details": {},
        "expected": {
            "risk_level": test_case.expected_risk_level,
            "signals": test_case.expected_signals,
            "probability_range": [test_case.expected_min_probability, test_case.expected_max_probability],
            "min_alternatives": test_case.expected_min_alternatives
        }
    }

    if test_case.skip_reason:
        test_result["status"] = "SKIPPED"
        test_result["details"]["skip_reason"] = test_case.skip_reason
        logger.warning(f"SKIPPED: {test_case.name} - {test_case.skip_reason}")
        return test_result

    try:
        async with httpx.AsyncClient(timeout=120.0) as client:
            payload = {
                "disease": "myeloma",
                "mutations": test_case.mutations,
                "current_drug_class": test_case.drug_class,
                "treatment_line": test_case.treatment_line,
                "prior_therapies": test_case.prior_therapies,
                "cytogenetics": test_case.cytogenetics
            }
            logger.debug(f"Sending payload: {json.dumps(payload, indent=2)}")
            response = await client.post(RESISTANCE_PREDICT_ENDPOINT, json=payload)
            response.raise_for_status()
            prediction = response.json()
            logger.debug(f"Received prediction: {json.dumps(prediction, indent=2)}")

            test_result["details"]["prediction"] = prediction

            # Assertions
            if prediction["risk_level"] == test_case.expected_risk_level:
                test_result["details"]["risk_level_match"] = True
            else:
                test_result["details"]["risk_level_match"] = False
                test_result["details"]["risk_level_actual"] = prediction["risk_level"]
                logger.error(f"  ❌ Risk level mismatch: Expected {test_case.expected_risk_level}, got {prediction['risk_level']}")

            if test_case.expected_min_probability <= prediction["probability"] <= test_case.expected_max_probability:
                test_result["details"]["probability_range_match"] = True
            else:
                test_result["details"]["probability_range_match"] = False
                test_result["details"]["probability_actual"] = prediction["probability"]
                logger.error(f"  ❌ Probability out of range: Expected {test_case.expected_min_probability}-{test_case.expected_max_probability}, got {prediction['probability']:.2f}")

            actual_signals = {s["signal_type"] for s in prediction["signals_detected"] if s["detected"]}
            expected_signals_set = set(test_case.expected_signals)
            
            if expected_signals_set.issubset(actual_signals):
                test_result["details"]["signals_match"] = True
            else:
                test_result["details"]["signals_match"] = False
                test_result["details"]["signals_expected"] = list(expected_signals_set)
                test_result["details"]["signals_actual"] = list(actual_signals)
                logger.error(f"  ❌ Signals mismatch: Expected {expected_signals_set}, got {actual_signals}")

            actual_alternatives_count = len(prediction.get("next_line_options", []))
            if actual_alternatives_count >= test_case.expected_min_alternatives:
                test_result["details"]["min_alternatives_met"] = True
            else:
                test_result["details"]["min_alternatives_met"] = False
                test_result["details"]["min_alternatives_actual"] = actual_alternatives_count
                logger.error(f"  ❌ Min alternatives not met: Expected {test_case.expected_min_alternatives}, got {actual_alternatives_count}")

            if (test_result["details"].get("risk_level_match", False) and
                test_result["details"].get("probability_range_match", False) and
                test_result["details"].get("signals_match", False) and
                test_result["details"].get("min_alternatives_met", False)):
                test_result["status"] = "PASS"
                logger.info(f"  ✅ PASS: {test_case.name}")
            else:
                test_result["status"] = "FAIL"
                logger.error(f"  ❌ FAIL: {test_case.name}")

    except httpx.HTTPStatusError as e:
        test_result["details"]["api_error"] = f"HTTP Error: {e.response.status_code} - {e.response.text}"
        logger.error(f"  ❌ API HTTP Error for {test_case.name}: {e.response.status_code} - {e.response.text}")
    except httpx.RequestError as e:
        test_result["details"]["api_error"] = f"Request Error: {e}"
        logger.error(f"  ❌ API Request Error for {test_case.name}: {e}")
    except Exception as e:
        test_result["details"]["exception"] = str(e)
        logger.error(f"  ❌ Unexpected Error for {test_case.name}: {e}", exc_info=True)
    
    return test_result

async def main():
    logger.info("=" * 80)
    logger.info("⚔️ MM Resistance Prediction Validation Framework")
    logger.info("=" * 80)

    all_test_results = []
    for test_case in TEST_CASES:
        result = await run_test_case(test_case)
        all_test_results.append(result)

    # Generate Report
    with open(REPORT_FILE, 'w') as f:
        f.write(f"# MM Resistance Prediction Validation Report\n\n")
        f.write(f"**Run Date:** {datetime.now().isoformat()}\n")
        f.write(f"**Total Tests:** {len(all_test_results)}\n")
        f.write(f"**Passed:** {sum(1 for r in all_test_results if r['status'] == 'PASS')}\n")
        f.write(f"**Skipped:** {sum(1 for r in all_test_results if r['status'] == 'SKIPPED')}\n")
        f.write(f"**Failed:** {sum(1 for r in all_test_results if r['status'] == 'FAIL')}\n\n")

        f.write("## Test Results Summary\n\n")
        for result in all_test_results:
            status_icon = "✅" if result["status"] == "PASS" else "❌" if result["status"] == "FAIL" else "⚠️"
            f.write(f"### {status_icon} {result['name']} ({result['status']})\n\n")
            f.write(f"**Description:** {result['description']}\n\n")
            
            if result["status"] == "SKIPPED":
                f.write(f"**Reason:** {result['details'].get('skip_reason', 'N/A')}\n\n")
            elif result["status"] == "FAIL":
                f.write(f"**Failure Details:**\n")
                for key, value in result["details"].items():
                    if key not in ["prediction", "risk_level_match", "probability_range_match", "signals_match", "min_alternatives_met"]:
                        f.write(f"- {key.replace('_', ' ').title()}: {value}\n")
                f.write("\n")
            
            f.write(f"**Expected:**\n")
            f.write(f"- Risk Level: {result['expected']['risk_level']}\n")
            f.write(f"- Probability Range: {result['expected']['probability_range'][0]:.2f}-{result['expected']['probability_range'][1]:.2f}\n")
            f.write(f"- Signals: {', '.join(result['expected']['signals'])}\n")
            f.write(f"- Min Alternatives: {result['expected']['min_alternatives']}\n\n")

            if "prediction" in result["details"]:
                pred = result["details"]["prediction"]
                f.write(f"**Actual Prediction:**\n")
                f.write(f"- Risk Level: {pred['risk_level']} ({'✅' if result['details'].get('risk_level_match') else '❌'})\n")
                f.write(f"- Probability: {pred['probability']:.2f} ({'✅' if result['details'].get('probability_range_match') else '❌'})\n")
                f.write(f"- Signals Detected ({pred['signal_count']}): {', '.join([s['signal_type'] for s in pred['signals_detected'] if s['detected']])} ({'✅' if result['details'].get('signals_match') else '❌'})\n")
                f.write(f"- Alternatives Count: {len(pred.get('next_line_options', []))} ({'✅' if result['details'].get('min_alternatives_met') else '❌'})\n\n")
            
            f.write("---\n\n")

    final_output = {
        "metadata": {
            "run_date": datetime.now().isoformat(),
            "total_tests": len(TEST_CASES),
            "passed_tests": sum(1 for r in all_test_results if r["status"] == "PASS"),
            "skipped_tests": sum(1 for r in all_test_results if r["status"] == "SKIPPED"),
            "failed_tests": sum(1 for r in all_test_results if r["status"] == "FAIL"),
        },
        "test_results": all_test_results
    }

    with open(RESULTS_JSON_FILE, 'w') as f:
        json.dump(final_output, f, indent=2)
    logger.info(f"✅ Validation results saved to {RESULTS_JSON_FILE}")
    logger.info(f"✅ Validation report generated: {REPORT_FILE}")

    logger.info("\nValidation framework run complete. Check the report for details.")
    
    if final_output["metadata"]["failed_tests"] > 0:
        sys.exit(1)

if __name__ == "__main__":
    import os
    asyncio.run(main())
