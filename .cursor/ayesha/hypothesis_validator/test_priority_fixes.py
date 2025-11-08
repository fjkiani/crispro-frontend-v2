#!/usr/bin/env python3
"""
Test script for Priority Fixes 1-4

Tests:
1. Evidence Synthesis (Fix 1) - Heuristic grading works
2. Dosage Extraction (Fix 2) - Regex patterns extract doses
3. SAE Coverage (Fix 3) - 20+ compounds have rules
4. Timing Recommendations (Fix 4) - Evidence-based fallback works
"""

import sys
from pathlib import Path

# Get absolute path to backend
script_dir = Path(__file__).parent.absolute()
workspace_root = script_dir.parent.parent.parent.absolute()
backend_path = workspace_root / "oncology-coPilot/oncology-backend-minimal"

if not backend_path.exists():
    print(f"❌ Backend path not found: {backend_path}")
    sys.exit(1)

# Add backend to path
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

import json

try:
    from api.services.enhanced_evidence_service import EnhancedEvidenceService
    from api.services.dietician_recommendations import DieticianRecommendationsService
    from api.services.food_treatment_line_service import compute_food_treatment_line_features
except ImportError as e:
    print(f"❌ Import error: {e}")
    print(f"Backend path: {backend_path}")
    print(f"Python path: {sys.path[:3]}")
    sys.exit(1)

def test_fix1_evidence_synthesis():
    """Test Fix 1: Evidence synthesis with heuristic grading"""
    print("\n" + "="*60)
    print("TEST 1: Evidence Synthesis (Fix 1)")
    print("="*60)
    
    service = EnhancedEvidenceService()
    
    test_cases = [
        {
            "name": "0 papers → INSUFFICIENT",
            "papers": [],
            "expected": "INSUFFICIENT"
        },
        {
            "name": "2 papers → WEAK",
            "papers": [
                {"pmid": "1", "title": "Study 1", "abstract": "..."},
                {"pmid": "2", "title": "Study 2", "abstract": "..."}
            ],
            "expected": "WEAK"
        },
        {
            "name": "5 papers → MODERATE",
            "papers": [
                {"pmid": str(i), "title": f"Study {i}", "abstract": "..."}
                for i in range(5)
            ],
            "expected": "MODERATE"
        },
        {
            "name": "3 RCTs → STRONG",
            "papers": [
                {"pmid": "1", "title": "Randomized controlled trial", "abstract": "randomized..."},
                {"pmid": "2", "title": "RCT Study", "abstract": "..."},
                {"pmid": "3", "title": "Randomized Study", "abstract": "..."},
                {"pmid": "4", "title": "Study 4", "abstract": "..."},
                {"pmid": "5", "title": "Study 5", "abstract": "..."}
            ],
            "expected": "STRONG"
        }
    ]
    
    passed = 0
    failed = 0
    
    for test in test_cases:
        grade = service._heuristic_grade(test["papers"])
        if grade == test["expected"]:
            print(f"✅ PASS: {test['name']} → {grade}")
            passed += 1
        else:
            print(f"❌ FAIL: {test['name']} → Expected {test['expected']}, got {grade}")
            failed += 1
    
    print(f"\nFix 1 Results: {passed}/{len(test_cases)} passed")
    return failed == 0


def test_fix2_dosage_extraction():
    """Test Fix 2: Dosage extraction with regex patterns"""
    print("\n" + "="*60)
    print("TEST 2: Dosage Extraction (Fix 2)")
    print("="*60)
    
    service = DieticianRecommendationsService()
    
    test_cases = [
        {
            "name": "Range pattern: 2000-4000 IU",
            "papers": [
                {
                    "pmid": "123",
                    "title": "Vitamin D Study",
                    "abstract": "Patients received 2000-4000 IU daily for 12 months."
                }
            ],
            "compound": "Vitamin D",
            "should_have_dose": True
        },
        {
            "name": "Single dose: 500 mg",
            "papers": [
                {
                    "pmid": "456",
                    "title": "NAC Study",
                    "abstract": "Dosage was 500 mg daily per day with food."
                }
            ],
            "compound": "NAC",
            "should_have_dose": True
        },
        {
            "name": "Decimal: 2.5 mg",
            "papers": [
                {
                    "pmid": "789",
                    "title": "Folate Study",
                    "abstract": "Administered 2.5 mg methylfolate daily."
                }
            ],
            "compound": "Folate",
            "should_have_dose": True
        },
        {
            "name": "No dose in papers",
            "papers": [
                {
                    "pmid": "999",
                    "title": "General Study",
                    "abstract": "This study examined the effects of the compound."
                }
            ],
            "compound": "Unknown",
            "should_have_dose": False
        }
    ]
    
    passed = 0
    failed = 0
    
    for test in test_cases:
        result = service.extract_dosage_from_evidence(test["papers"], test["compound"])
        has_dose = bool(result.get("recommended_dose"))
        
        if has_dose == test["should_have_dose"]:
            if has_dose:
                print(f"✅ PASS: {test['name']} → Found: {result['recommended_dose']}")
            else:
                print(f"✅ PASS: {test['name']} → No dose (expected)")
            passed += 1
        else:
            print(f"❌ FAIL: {test['name']} → Expected dose={test['should_have_dose']}, got {has_dose}")
            failed += 1
    
    print(f"\nFix 2 Results: {passed}/{len(test_cases)} passed")
    return failed == 0


def test_fix3_sae_coverage():
    """Test Fix 3: SAE coverage for 20+ compounds"""
    print("\n" + "="*60)
    print("TEST 3: SAE Coverage (Fix 3)")
    print("="*60)
    
    # Test compounds
    test_compounds = [
        {"name": "NAC", "expected_score": 1.0, "has_rule": True},
        {"name": "Vitamin D", "expected_score": 0.9, "has_rule": True},
        {"name": "Omega-3", "expected_score": 0.85, "has_rule": True},
        {"name": "Curcumin", "expected_score": 0.7, "has_rule": True},
        {"name": "Resveratrol", "expected_score": 0.7, "has_rule": True},  # New compound
        {"name": "Quercetin", "expected_score": 0.7, "has_rule": True},  # New compound
        {"name": "CoQ10", "expected_score": 0.75, "has_rule": True},  # New compound
        {"name": "UnknownCompoundXYZ", "expected_score": 0.6, "has_rule": False}  # Should use default
    ]
    
    disease_context = {
        "biomarkers": {"HRD": "POSITIVE"},
        "pathways_disrupted": ["DNA repair"]
    }
    
    treatment_history = {
        "prior_therapies": ["carboplatin", "paclitaxel"]
    }
    
    passed = 0
    failed = 0
    
    for test in test_compounds:
        result = compute_food_treatment_line_features(
            compound=test["name"],
            disease_context=disease_context,
            treatment_history=treatment_history
        )
        
        line_app = result.get("line_appropriateness", 0)
        expected_min = test["expected_score"]
        expected_max = min(expected_min + 0.15, 1.0)  # Allow for biomarker boosts
        
        if expected_min <= line_app <= expected_max:
            print(f"✅ PASS: {test['name']} → line_appropriateness={line_app:.2f} (expected ~{expected_min:.2f})")
            passed += 1
        else:
            print(f"❌ FAIL: {test['name']} → Expected ~{expected_min:.2f}, got {line_app:.2f}")
            failed += 1
    
    # Also test that supplement_treatment_rules.json exists and has 20+ compounds
    rules_file = Path(__file__).parent / "data/supplement_treatment_rules.json"
    if rules_file.exists():
        with open(rules_file) as f:
            rules_data = json.load(f)
            compound_count = len(rules_data.get("supplement_rules", {}))
            if compound_count >= 20:
                print(f"✅ PASS: supplement_treatment_rules.json has {compound_count} compounds")
                passed += 1
            else:
                print(f"❌ FAIL: Expected ≥20 compounds, got {compound_count}")
                failed += 1
    else:
        print(f"❌ FAIL: supplement_treatment_rules.json not found")
        failed += 1
    
    print(f"\nFix 3 Results: {passed}/{len(test_compounds) + 1} passed")
    return failed == 0


def test_fix4_timing_recommendations():
    """Test Fix 4: Timing recommendations with evidence fallback"""
    print("\n" + "="*60)
    print("TEST 4: Timing Recommendations (Fix 4)")
    print("="*60)
    
    service = DieticianRecommendationsService()
    
    test_cases = [
        {
            "name": "Vitamin D (hardcoded) → Morning with breakfast",
            "compound": "Vitamin D",
            "evidence": {},
            "expected_time": "Morning"
        },
        {
            "name": "Resveratrol (unknown) with morning in papers → Morning",
            "compound": "Resveratrol",
            "evidence": {
                "papers": [
                    {
                        "pmid": "111",
                        "abstract": "Patients took resveratrol in the morning with breakfast for optimal absorption."
                    }
                ]
            },
            "expected_time": "Morning"
        },
        {
            "name": "Unknown compound with evening in papers → Evening",
            "compound": "UnknownSupplement",
            "evidence": {
                "papers": [
                    {
                        "pmid": "222",
                        "abstract": "Administration at bedtime showed better results."
                    }
                ]
            },
            "expected_time": "Evening"
        },
        {
            "name": "Unknown compound with 'with food' → With meals",
            "compound": "AnotherUnknown",
            "evidence": {
                "papers": [
                    {
                        "pmid": "333",
                        "abstract": "Take with meals to reduce gastrointestinal side effects."
                    }
                ]
            },
            "expected_time": "With meals"
        },
        {
            "name": "Unknown compound, no evidence → As directed",
            "compound": "NoEvidenceCompound",
            "evidence": {},
            "expected_time": "As directed"
        }
    ]
    
    passed = 0
    failed = 0
    
    for test in test_cases:
        result = service.generate_timing_recommendations(test["compound"], test["evidence"])
        best_time = result.get("best_time", "")
        
        # Check if expected time is in result (for hardcoded, might have more detail)
        if test["expected_time"].lower() in best_time.lower():
            print(f"✅ PASS: {test['name']} → {best_time}")
            passed += 1
        else:
            print(f"❌ FAIL: {test['name']} → Expected '{test['expected_time']}', got '{best_time}'")
            failed += 1
    
    print(f"\nFix 4 Results: {passed}/{len(test_cases)} passed")
    return failed == 0


def main():
    """Run all tests"""
    print("="*60)
    print("PRIORITY FIXES TEST SUITE")
    print("="*60)
    
    results = []
    
    # Test Fix 1
    results.append(("Fix 1: Evidence Synthesis", test_fix1_evidence_synthesis()))
    
    # Test Fix 2
    results.append(("Fix 2: Dosage Extraction", test_fix2_dosage_extraction()))
    
    # Test Fix 3
    results.append(("Fix 3: SAE Coverage", test_fix3_sae_coverage()))
    
    # Test Fix 4
    results.append(("Fix 4: Timing Recommendations", test_fix4_timing_recommendations()))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    passed_count = sum(1 for _, result in results if result)
    total_count = len(results)
    
    for name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {name}")
    
    print(f"\nOverall: {passed_count}/{total_count} test suites passed")
    
    if passed_count == total_count:
        print("\n🎉 ALL TESTS PASSED - PRIORITY FIXES VALIDATED")
        return 0
    else:
        print(f"\n⚠️ {total_count - passed_count} test suite(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

