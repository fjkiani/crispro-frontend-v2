#!/usr/bin/env python3
"""
Smoke Test for Ayesha Trial Filtering Engine

Tests core logic without requiring full backend stack.
"""
import sys
import os

# Add backend to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'oncology-coPilot', 'oncology-backend-minimal'))

def test_ca125_intelligence():
    """Test CA-125 intelligence service."""
    print("🧪 Testing CA-125 Intelligence Service...")
    
    try:
        from api.services.ca125_intelligence import CA125IntelligenceService
        
        service = CA125IntelligenceService()
        
        # Test Ayesha's value (2842)
        result = service.analyze(2842.0)
        
        assert result['ca125_value'] == 2842.0, "CA-125 value mismatch"
        assert result['disease_burden'] == 'EXTENSIVE', f"Expected EXTENSIVE, got {result['disease_burden']}"
        assert result['marker_utility'] == 'HIGHLY_TRACKABLE', "Marker utility mismatch"
        assert 'expected_response' in result, "Missing expected_response"
        assert 'trial_preferences' in result, "Missing trial_preferences"
        
        print("  ✅ CA-125 Intelligence: PASS")
        print(f"     Burden: {result['disease_burden']}")
        print(f"     Cycle 3 Expected: {result['expected_response'].get('cycle3_expected_drop', 'N/A')}")
        print(f"     Cycle 6 Expected: {result['expected_response'].get('cycle6_expected_drop', 'N/A')}")
        
        return True
    except Exception as e:
        print(f"  ❌ CA-125 Intelligence: FAIL - {e}")
        return False

def test_schemas():
    """Test Pydantic schemas."""
    print("\n🧪 Testing Schemas...")
    
    try:
        from api.schemas.ayesha_trials import (
            AyeshaTrialProfile,
            TrialMatchReasoning,
            AyeshaTrialMatch,
            AyeshaTrialSearchResponse
        )
        
        # Test default profile
        profile = AyeshaTrialProfile()
        assert profile.ca125 == 2842.0, "Default CA-125 mismatch"
        assert profile.germline_status == "NEGATIVE", "Default germline status mismatch"
        assert profile.stage == "IVB", "Default stage mismatch"
        
        # Test reasoning
        reasoning = TrialMatchReasoning(
            match_score=0.85,
            why_eligible=["✅ Test reason"],
            why_good_fit=["🎯 Test fit"],
            evidence_tier="STANDARD",
            enrollment_likelihood="HIGH"
        )
        assert reasoning.match_score == 0.85, "Match score mismatch"
        
        print("  ✅ Schemas: PASS")
        print(f"     Default profile: CA-125={profile.ca125}, Stage={profile.stage}")
        
        return True
    except Exception as e:
        print(f"  ❌ Schemas: FAIL - {e}")
        import traceback
        traceback.print_exc()
        return False

def test_scoring_engine():
    """Test scoring engine (without HybridTrialSearchService dependency)."""
    print("\n🧪 Testing Scoring Engine Logic...")
    
    try:
        from api.services.ayesha_trial_matching.scoring_engine import ScoringEngine
        from api.schemas.ayesha_trials import AyeshaTrialProfile
        from api.services.ca125_intelligence import CA125IntelligenceService
        
        scoring_engine = ScoringEngine()
        profile = AyeshaTrialProfile()
        ca125_intel = CA125IntelligenceService().analyze(profile.ca125)
        
        # Mock trial
        mock_trial = {
            "title": "First-line Ovarian Cancer Trial",
            "description": "Stage IV first-line treatment",
            "phase": "Phase 3",
            "interventions": ["Carboplatin", "Paclitaxel", "Bevacizumab"],
            "locations": [{"state": "NY", "city": "New York", "facility": "Mount Sinai"}],
            "eligibility_text": "All comers, germline BRCA not required"
        }
        
        score = scoring_engine.calculate_match_score(mock_trial, profile, ca125_intel)
        
        assert 0.0 <= score <= 1.0, f"Score out of range: {score}"
        assert score > 0.5, f"Score too low: {score} (should be boosted)"
        
        print("  ✅ Scoring Engine: PASS")
        print(f"     Match Score: {score:.2f} ({score*100:.0f}%)")
        
        return True
    except Exception as e:
        print(f"  ❌ Scoring Engine: FAIL - {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all smoke tests."""
    print("⚔️ AYESHA TRIAL FILTERING ENGINE - SMOKE TESTS ⚔️\n")
    
    results = []
    
    # Test 1: CA-125 Intelligence
    results.append(test_ca125_intelligence())
    
    # Test 2: Schemas
    results.append(test_schemas())
    
    # Test 3: Scoring Engine
    results.append(test_scoring_engine())
    
    # Summary
    print("\n" + "="*60)
    passed = sum(results)
    total = len(results)
    
    if passed == total:
        print(f"✅ ALL TESTS PASSED ({passed}/{total})")
        return 0
    else:
        print(f"⚠️  SOME TESTS FAILED ({passed}/{total} passed)")
        return 1

if __name__ == "__main__":
    sys.exit(main())


