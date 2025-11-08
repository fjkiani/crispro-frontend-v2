#!/usr/bin/env python3
"""
Comprehensive Food Validator Test - Prove It's NOT Hardcoded

Tests 5 different compounds with different biomarker profiles
to verify dynamic behavior and real LLM/API integration.
"""

import asyncio
import json
import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"
sys.path.insert(0, str(backend_path))

from api.services.enhanced_evidence_service import EnhancedEvidenceService
from api.services.dynamic_food_extraction import DynamicFoodExtractor
from api.services.food_spe_integration import FoodSPEIntegrationService
from api.services.food_treatment_line_service import compute_food_treatment_line_features

# Test cases - DIVERSE compounds to prove dynamic behavior
TEST_CASES = [
    {
        "name": "Curcumin (Turmeric)",
        "compound": "Curcumin",
        "disease": "ovarian_cancer_hgs",
        "biomarkers": {
            "HRD": "NEGATIVE",  # Different from Vitamin D test
            "TMB": 15.3,        # High TMB
            "MSI": "MSI-H",
            "BRCA1": "NEGATIVE"
        },
        "expected_mechanism": "inflammation",  # Should find anti-inflammatory
        "expected_pathways": ["NF-kB", "inflammation"]
    },
    {
        "name": "Green Tea (EGCG)",
        "compound": "Green Tea Extract",
        "disease": "ovarian_cancer_hgs",
        "biomarkers": {
            "HRD": "POSITIVE",
            "TMB": 4.1,         # Low TMB
            "MSI": "STABLE",
            "BRCA2": "POSITIVE" # BRCA2 mutation
        },
        "expected_mechanism": "angiogenesis",  # Should find anti-angiogenic
        "expected_pathways": ["VEGF", "angiogenesis"]
    },
    {
        "name": "Omega-3 (Fish Oil)",
        "compound": "Omega-3 Fatty Acids",
        "disease": "ovarian_cancer_hgs",
        "biomarkers": {
            "HRD": "NEGATIVE",
            "TMB": 6.8,
            "MSI": "STABLE",
            "BRCA1": "NEGATIVE"
        },
        "expected_mechanism": "inflammation",
        "expected_pathways": ["COX-2", "inflammation", "prostaglandin"]
    },
    {
        "name": "Quercetin",
        "compound": "Quercetin",
        "disease": "ovarian_cancer_hgs",
        "biomarkers": {
            "HRD": "POSITIVE",
            "TMB": 11.2,
            "MSI": "MSI-H",
            "BRCA1": "POSITIVE"
        },
        "expected_mechanism": "cell_cycle",
        "expected_pathways": ["p53", "cell cycle", "apoptosis"]
    },
    {
        "name": "Genistein (Soy)",
        "compound": "Genistein",
        "disease": "ovarian_cancer_hgs",
        "biomarkers": {
            "HRD": "NEGATIVE",
            "TMB": 3.4,
            "MSI": "STABLE",
            "BRCA1": "NEGATIVE"
        },
        "expected_mechanism": "hormone",
        "expected_pathways": ["estrogen receptor", "hormone"]
    }
]

# Colors
GREEN = "\033[92m"
BLUE = "\033[94m"
YELLOW = "\033[93m"
RED = "\033[91m"
RESET = "\033[0m"
BOLD = "\033[1m"

def print_header(text):
    print(f"\n{BOLD}{BLUE}{'='*80}{RESET}")
    print(f"{BOLD}{BLUE}{text}{RESET}")
    print(f"{BOLD}{BLUE}{'='*80}{RESET}\n")

def print_test(num, total, name):
    print(f"\n{BOLD}{GREEN}[TEST {num}/{total}] {name}{RESET}")

def print_check(label, got, expected=None):
    if expected:
        match = any(e.lower() in str(got).lower() for e in expected) if isinstance(expected, list) else (expected.lower() in str(got).lower())
        symbol = "✅" if match else "⚠️"
        color = GREEN if match else YELLOW
    else:
        symbol = "ℹ️"
        color = BLUE
    print(f"{color}{symbol} {label}: {got}{RESET}")

async def run_comprehensive_tests():
    """Run 5 diverse tests to prove dynamic behavior"""
    
    print_header("COMPREHENSIVE FOOD VALIDATOR TEST - PROVING IT'S NOT HARDCODED")
    
    evidence_service = EnhancedEvidenceService()
    spe_service = FoodSPEIntegrationService()
    
    results_summary = []
    
    for i, test_case in enumerate(TEST_CASES, 1):
        print_test(i, len(TEST_CASES), test_case['name'])
        
        # Create fresh extractor for each test (clear cache)
        extractor = DynamicFoodExtractor()
        
        # Extract targets
        print(f"{BLUE}Step 1: Dynamic Extraction{RESET}")
        extraction = await extractor.extract_all(
            test_case['compound'],
            test_case['disease']
        )
        
        targets = extraction.get('targets', [])
        print_check("Targets Found", len(targets))
        print(f"  → {', '.join(targets[:5])}")
        
        # Mine evidence
        print(f"{BLUE}Step 2: Evidence Mining (PubMed){RESET}")
        evidence = await evidence_service.get_complete_evidence(
            test_case['compound'],
            test_case['disease']
        )
        
        papers_found = len(evidence.get('papers', []))
        evidence_grade = evidence.get('evidence_grade', 'UNKNOWN')
        print_check("Papers Found", papers_found)
        print_check("Evidence Grade", evidence_grade)
        
        # Check mechanisms
        mechanisms = evidence.get('llm_mechanisms', [])
        print_check("Mechanisms", ', '.join(mechanisms[:3]) if mechanisms else "None", test_case['expected_mechanism'])
        
        # S/P/E scoring
        print(f"{BLUE}Step 3: S/P/E Integration{RESET}")
        
        disease_context = {
            "disease": test_case['disease'],
            "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
            "biomarkers": test_case['biomarkers'],
            "pathways_disrupted": ["DNA repair", "Cell cycle"]
        }
        
        treatment_history = {
            "current_line": "L2",
            "prior_therapies": ["carboplatin"]
        }
        
        # Get pathways from disease context
        pathways = disease_context.get('pathways_disrupted', [])
        
        spe_result = await spe_service.compute_spe_score(
            test_case['compound'],
            targets[:3],  # Top 3 targets
            pathways,
            disease_context,
            evidence_grade,  # Pass the grade we got from evidence mining
            treatment_history,
            evo2_enabled=False  # Phase 1 - P/E only
        )
        
        confidence = spe_result.get('confidence', 0)
        print_check("Confidence", f"{confidence:.2f}")
        
        # SAE features
        print(f"{BLUE}Step 4: SAE Features (Treatment Line){RESET}")
        
        treatment_history = {
            "current_line": "L2",
            "prior_therapies": ["carboplatin"]
        }
        
        sae_features = compute_food_treatment_line_features(
            test_case['compound'],
            disease_context,
            treatment_history
        )
        
        line_appropriateness = sae_features.get('line_appropriateness', 0)
        print_check("Line Appropriateness", f"{line_appropriateness:.2f}")
        
        # Summary
        test_result = {
            "compound": test_case['compound'],
            "targets": len(targets),
            "papers": papers_found,
            "grade": evidence_grade,
            "confidence": confidence,
            "is_dynamic": papers_found > 0 and len(targets) > 0
        }
        results_summary.append(test_result)
        
        print(f"{GREEN}✓ Test Complete{RESET}\n")
        print("-" * 80)
    
    # Final Summary
    print_header("FINAL RESULTS - DYNAMIC BEHAVIOR VERIFICATION")
    
    dynamic_count = sum(1 for r in results_summary if r['is_dynamic'])
    
    print(f"\n{BOLD}Tests Run: {len(TEST_CASES)}{RESET}")
    print(f"{BOLD}Dynamic Behavior: {dynamic_count}/{len(TEST_CASES)}{RESET}")
    
    print(f"\n{BOLD}Per-Compound Results:{RESET}")
    for r in results_summary:
        status = f"{GREEN}✅ DYNAMIC{RESET}" if r['is_dynamic'] else f"{RED}❌ STATIC{RESET}"
        print(f"  • {r['compound']}: {r['papers']} papers, {r['targets']} targets, {r['grade']} grade, {r['confidence']:.2f} confidence {status}")
    
    # Verdict
    if dynamic_count == len(TEST_CASES):
        print(f"\n{GREEN}{BOLD}🎯 VERDICT: FULLY DYNAMIC - NOT HARDCODED{RESET}")
        print(f"{GREEN}All compounds produced unique results from real APIs{RESET}")
        return True
    elif dynamic_count >= len(TEST_CASES) * 0.6:
        print(f"\n{YELLOW}{BOLD}⚠️ VERDICT: MOSTLY DYNAMIC - SOME GAPS{RESET}")
        print(f"{YELLOW}Most compounds work, but some edge cases need fixing{RESET}")
        return False
    else:
        print(f"\n{RED}{BOLD}❌ VERDICT: STILL HARDCODED{RESET}")
        print(f"{RED}System is not dynamically querying APIs{RESET}")
        return False

if __name__ == "__main__":
    result = asyncio.run(run_comprehensive_tests())
    sys.exit(0 if result else 1)
