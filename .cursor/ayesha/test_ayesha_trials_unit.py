#!/usr/bin/env python3
"""
Unit Tests for Ayesha Trial Filtering Engine

Tests core logic without requiring database connections.
"""
import sys
import os

# Add backend to path
backend_path = os.path.join(os.path.dirname(__file__), '..', 'oncology-coPilot', 'oncology-backend-minimal')
sys.path.insert(0, backend_path)

def test_ca125_intelligence():
    """Test CA-125 intelligence service."""
    print("🧪 Test 1: CA-125 Intelligence Service")
    
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
        assert 'boost_keywords' in result, "Missing boost_keywords"
        
        print("  ✅ PASS")
        print(f"     Burden: {result['disease_burden']}")
        print(f"     Cycle 3 Expected: {result['expected_response']['cycle3_expected_drop']}")
        print(f"     Cycle 6 Expected: {result['expected_response']['cycle6_expected_drop']}")
        print(f"     Boost Keywords: {len(result['boost_keywords'])} keywords")
        
        return True
    except Exception as e:
        print(f"  ❌ FAIL - {e}")
        import traceback
        traceback.print_exc()
        return False

def test_schemas():
    """Test Pydantic schemas."""
    print("\n🧪 Test 2: Pydantic Schemas")
    
    try:
        from api.schemas.ayesha_trials import (
            AyeshaTrialProfile,
            TrialMatchReasoning,
            AyeshaTrialMatch,
            AyeshaTrialSearchResponse
        )
        
        # Test default profile
        profile = AyeshaTrialProfile()
        assert profile.ca125 == 2842.0, f"Default CA-125 mismatch: {profile.ca125}"
        assert profile.germline_status == "NEGATIVE", f"Default germline status mismatch: {profile.germline_status}"
        assert profile.stage == "IVB", f"Default stage mismatch: {profile.stage}"
        assert profile.treatment_line == 0, "Default treatment line mismatch"
        
        # Test reasoning
        reasoning = TrialMatchReasoning(
            match_score=0.85,
            why_eligible=["✅ Test reason 1", "✅ Test reason 2"],
            why_good_fit=["🎯 Test fit 1"],
            conditional_requirements=["⚠️ Test conditional"],
            red_flags=[],
            evidence_tier="STANDARD",
            enrollment_likelihood="HIGH"
        )
        assert reasoning.match_score == 0.85, "Match score mismatch"
        assert len(reasoning.why_eligible) == 2, "Why eligible count mismatch"
        assert reasoning.evidence_tier == "STANDARD", "Evidence tier mismatch"
        
        print("  ✅ PASS")
        print(f"     Default profile: CA-125={profile.ca125}, Stage={profile.stage}, Germline={profile.germline_status}")
        print(f"     Reasoning: Score={reasoning.match_score}, Tier={reasoning.evidence_tier}")
        
        return True
    except Exception as e:
        print(f"  ❌ FAIL - {e}")
        import traceback
        traceback.print_exc()
        return False

def test_scoring_engine_logic():
    """Test scoring engine logic (without database dependencies)."""
    print("\n🧪 Test 3: Scoring Engine Logic")
    
    try:
        from api.services.ayesha_trial_matching.scoring_engine import ScoringEngine
        from api.schemas.ayesha_trials import AyeshaTrialProfile
        from api.services.ca125_intelligence import CA125IntelligenceService
        
        scoring_engine = ScoringEngine()
        profile = AyeshaTrialProfile()
        ca125_intel = CA125IntelligenceService().analyze(profile.ca125)
        
        # Mock trial with all boosts
        mock_trial = {
            "title": "First-line Ovarian Cancer Stage IV Trial",
            "description": "Stage IV first-line treatment with carboplatin and paclitaxel",
            "phase": "Phase 3",
            "interventions": ["Carboplatin", "Paclitaxel", "Bevacizumab"],
            "locations": [{"state": "NY", "city": "New York", "facility": "Mount Sinai"}],
            "eligibility_text": "All comers, germline BRCA not required, first-line therapy"
        }
        
        score = scoring_engine.calculate_match_score(mock_trial, profile, ca125_intel)
        
        assert 0.0 <= score <= 1.0, f"Score out of range: {score}"
        assert score > 0.5, f"Score too low: {score} (should be boosted by multiple factors)"
        
        # Test with minimal trial (should have lower score)
        minimal_trial = {
            "title": "Ovarian Cancer Trial",
            "description": "Treatment trial",
            "phase": "Phase 1",
            "interventions": ["Unknown Drug"],
            "locations": [{"state": "CA", "city": "Los Angeles"}],
            "eligibility_text": "Ovarian cancer"
        }
        
        minimal_score = scoring_engine.calculate_match_score(minimal_trial, profile, ca125_intel)
        assert minimal_score < score, "Minimal trial should score lower than boosted trial"
        
        print("  ✅ PASS")
        print(f"     Boosted Trial Score: {score:.2f} ({score*100:.0f}%)")
        print(f"     Minimal Trial Score: {minimal_score:.2f} ({minimal_score*100:.0f}%)")
        print(f"     Score Difference: {score - minimal_score:.2f}")
        
        return True
    except Exception as e:
        print(f"  ❌ FAIL - {e}")
        import traceback
        traceback.print_exc()
        return False

def test_reasoning_generator():
    """Test reasoning generator."""
    print("\n🧪 Test 4: Reasoning Generator")
    
    try:
        from api.services.ayesha_trial_matching.reasoning_generator import ReasoningGenerator
        from api.schemas.ayesha_trials import AyeshaTrialProfile
        from api.services.ca125_intelligence import CA125IntelligenceService
        
        generator = ReasoningGenerator()
        profile = AyeshaTrialProfile()
        ca125_intel = CA125IntelligenceService().analyze(profile.ca125)
        
        mock_trial = {
            "title": "First-line Ovarian Cancer Stage IV Trial",
            "description": "Stage IV first-line treatment",
            "phase": "Phase 3",
            "interventions": ["Carboplatin", "Paclitaxel", "Bevacizumab"],
            "locations": [{"state": "NY", "city": "New York", "facility": "Mount Sinai"}],
            "eligibility_text": "All comers, germline BRCA not required"
        }
        
        reasoning = generator.generate_complete_reasoning(
            mock_trial, profile, ca125_intel, 0.85
        )
        
        assert reasoning.match_score == 0.85, "Match score mismatch"
        assert len(reasoning.why_eligible) > 0, "Should have eligibility reasons"
        assert len(reasoning.why_good_fit) > 0, "Should have good fit reasons"
        assert reasoning.evidence_tier in ["STANDARD", "SUPPORTED", "INVESTIGATIONAL"], "Invalid evidence tier"
        assert reasoning.enrollment_likelihood in ["HIGH", "MEDIUM", "LOW"], "Invalid enrollment likelihood"
        
        print("  ✅ PASS")
        print(f"     Why Eligible: {len(reasoning.why_eligible)} reasons")
        print(f"     Why Good Fit: {len(reasoning.why_good_fit)} reasons")
        print(f"     Evidence Tier: {reasoning.evidence_tier}")
        print(f"     Enrollment Likelihood: {reasoning.enrollment_likelihood}")
        
        return True
    except Exception as e:
        print(f"  ❌ FAIL - {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all unit tests."""
    print("⚔️ AYESHA TRIAL FILTERING ENGINE - UNIT TESTS ⚔️\n")
    
    results = []
    
    # Test 1: CA-125 Intelligence
    results.append(test_ca125_intelligence())
    
    # Test 2: Schemas
    results.append(test_schemas())
    
    # Test 3: Scoring Engine
    results.append(test_scoring_engine_logic())
    
    # Test 4: Reasoning Generator
    results.append(test_reasoning_generator())
    
    # Summary
    print("\n" + "="*60)
    passed = sum(results)
    total = len(results)
    
    if passed == total:
        print(f"✅ ALL TESTS PASSED ({passed}/{total})")
        print("\n🎯 Core logic verified! Ready for integration testing.")
        return 0
    else:
        print(f"⚠️  SOME TESTS FAILED ({passed}/{total} passed)")
        return 1

if __name__ == "__main__":
    sys.exit(main())


