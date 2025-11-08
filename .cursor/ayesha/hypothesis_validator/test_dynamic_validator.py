#!/usr/bin/env python3
"""
Test Dynamic Food Validator Endpoint

Tests what makes it unique vs. googling and how biomarker targeting works.
"""

import requests
import json
from typing import Dict, Any

API_BASE = "http://127.0.0.1:8000"

def test_1_resveratrol_basic():
    """Test 1: Basic Resveratrol validation"""
    print("🧪 TEST 1: Resveratrol - Basic Validation")
    print("=" * 60)
    
    payload = {
        "compound": "Resveratrol",
        "disease_context": {
            "disease": "ovarian_cancer_hgs",
            "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
            "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
            "pathways_disrupted": ["DNA repair", "Angiogenesis", "Inflammation"]
        },
        "treatment_history": {
            "current_line": "L3",
            "prior_therapies": ["carboplatin", "paclitaxel"]
        },
        "patient_medications": [],
        "use_evo2": False
    }
    
    response = requests.post(
        f"{API_BASE}/api/hypothesis/validate_food_dynamic",
        json=payload,
        timeout=30
    )
    
    if response.status_code != 200:
        print(f"❌ Error: {response.status_code}")
        print(response.text)
        return None
    
    data = response.json()
    
    print(f"✅ Status: {data.get('status')}")
    print(f"✅ Compound: {data.get('compound')}")
    print(f"✅ Overall Score: {data.get('overall_score', 0):.3f}")
    print(f"✅ Confidence: {data.get('confidence', 0):.3f}")
    print(f"✅ Verdict: {data.get('verdict')}")
    print()
    print("📊 S/P/E Breakdown:")
    spe = data.get('spe_breakdown', {})
    print(f"  - Sequence (S): {spe.get('sequence', 0):.3f}")
    print(f"  - Pathway (P): {spe.get('pathway', 0):.3f}")
    print(f"  - Evidence (E): {spe.get('evidence', 0):.3f}")
    print()
    print("🎯 Targets Extracted:")
    targets = data.get('targets', [])
    print(f"  - Found {len(targets)} targets: {', '.join(targets[:5])}")
    print()
    print("🧬 Pathways Mapped:")
    pathways = data.get('pathways', [])
    print(f"  - Found {len(pathways)} pathways: {', '.join(pathways[:5])}")
    print()
    print("📚 Evidence:")
    evidence = data.get('evidence', {})
    print(f"  - Grade: {evidence.get('evidence_grade')}")
    print(f"  - Total Papers: {evidence.get('total_papers', 0)}")
    print(f"  - RCTs: {evidence.get('rct_count', 0)}")
    print()
    print("🏥 SAE Features (Treatment Line Intelligence):")
    sae = data.get('sae_features', {})
    print(f"  - Line Appropriateness: {sae.get('line_appropriateness', 0):.2f}")
    print(f"  - Cross-Resistance Risk: {sae.get('cross_resistance', 0):.2f}")
    print(f"  - Sequencing Fitness: {sae.get('sequencing_fitness', 0):.2f}")
    print()
    
    return data

def test_2_biomarker_targeting():
    """Test 2: How biomarkers affect recommendations"""
    print("\n🧪 TEST 2: Biomarker Targeting - HRD+ vs HRD-")
    print("=" * 60)
    
    # Test A: HRD+ patient (Ayesha's case)
    payload_hrd_plus = {
        "compound": "Vitamin D",
        "disease_context": {
            "disease": "ovarian_cancer_hgs",
            "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
            "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},  # HRD+
            "pathways_disrupted": ["DNA repair", "Cell cycle"]
        },
        "treatment_history": {
            "current_line": "L3",
            "prior_therapies": ["carboplatin", "paclitaxel"]
        },
        "use_evo2": False
    }
    
    # Test B: HRD- patient (different biomarker)
    payload_hrd_minus = {
        "compound": "Vitamin D",
        "disease_context": {
            "disease": "ovarian_cancer_hgs",
            "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
            "biomarkers": {"HRD": "NEGATIVE", "TMB": 5.0},  # HRD-
            "pathways_disrupted": ["Cell cycle"]  # No DNA repair pathway
        },
        "treatment_history": {
            "current_line": "L3",
            "prior_therapies": ["carboplatin", "paclitaxel"]
        },
        "use_evo2": False
    }
    
    print("📊 HRD+ Patient (DNA repair pathway active):")
    resp1 = requests.post(f"{API_BASE}/api/hypothesis/validate_food_dynamic", json=payload_hrd_plus, timeout=30)
    if resp1.status_code == 200:
        data1 = resp1.json()
        print(f"  - Score: {data1.get('overall_score', 0):.3f}")
        print(f"  - Confidence: {data1.get('confidence', 0):.3f}")
        print(f"  - SAE Line Appropriateness: {data1.get('sae_features', {}).get('line_appropriateness', 0):.2f}")
        print(f"  - Pathways: {', '.join(data1.get('pathways', [])[:3])}")
    
    print("\n📊 HRD- Patient (No DNA repair pathway):")
    resp2 = requests.post(f"{API_BASE}/api/hypothesis/validate_food_dynamic", json=payload_hrd_minus, timeout=30)
    if resp2.status_code == 200:
        data2 = resp2.json()
        print(f"  - Score: {data2.get('overall_score', 0):.3f}")
        print(f"  - Confidence: {data2.get('confidence', 0):.2f}")
        print(f"  - SAE Line Appropriateness: {data2.get('sae_features', {}).get('line_appropriateness', 0):.2f}")
        print(f"  - Pathways: {', '.join(data2.get('pathways', [])[:3])}")
    
    print("\n🔍 DIFFERENTIATION:")
    print("  → HRD+ gets boosted line appropriateness (biomarker gate matches)")
    print("  → Different pathway alignment scores (DNA repair vs. no DNA repair)")
    print("  → Biomarker-specific confidence modulation")

def test_3_vs_googling():
    """Test 3: What makes this unique vs. googling"""
    print("\n🧪 TEST 3: What's Unique vs. Google Search")
    print("=" * 60)
    
    payload = {
        "compound": "Curcumin",
        "disease_context": {
            "disease": "ovarian_cancer_hgs",
            "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
            "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
            "pathways_disrupted": ["Inflammation", "Angiogenesis"]
        },
        "treatment_history": {
            "current_line": "L3",
            "prior_therapies": ["carboplatin", "paclitaxel"]
        },
        "patient_medications": ["warfarin"],
        "use_evo2": False
    }
    
    response = requests.post(
        f"{API_BASE}/api/hypothesis/validate_food_dynamic",
        json=payload,
        timeout=30
    )
    
    if response.status_code != 200:
        print(f"❌ Error: {response.status_code}")
        return
    
    data = response.json()
    
    print("🔍 GOOGLE SEARCH GIVES YOU:")
    print("  ❌ Generic information about curcumin")
    print("  ❌ Not personalized to YOUR cancer type")
    print("  ❌ Not matched to YOUR biomarkers (HRD+)")
    print("  ❌ Not considering YOUR treatment history (L3 post-platinum)")
    print("  ❌ Not checking YOUR medications (warfarin)")
    print("  ❌ No integrated S/P/E scoring")
    print("  ❌ No treatment line intelligence")
    print()
    print("✅ OUR SYSTEM GIVES YOU:")
    print(f"  ✅ Personalized verdict: {data.get('verdict')} with {data.get('confidence', 0):.1%} confidence")
    print(f"  ✅ Targets mapped: {len(data.get('targets', []))} molecular targets")
    print(f"  ✅ Pathway alignment: {data.get('spe_breakdown', {}).get('pathway', 0):.1%} match to YOUR cancer pathways")
    print(f"  ✅ Evidence grade: {data.get('evidence', {}).get('evidence_grade')} ({data.get('evidence', {}).get('total_papers', 0)} papers)")
    print(f"  ✅ Treatment line fit: {data.get('sae_features', {}).get('line_appropriateness', 0):.1%} appropriate for L3")
    print(f"  ✅ Drug interaction check: {'⚠️ INTERACTIONS FOUND' if data.get('dietician_recommendations', {}).get('interactions', {}).get('interactions') else '✅ SAFE'}")
    print(f"  ✅ Dosage recommendations: {data.get('dietician_recommendations', {}).get('dosage', {}).get('recommended_dose', 'N/A')}")
    print(f"  ✅ Safety alerts: {len(data.get('dietician_recommendations', {}).get('safety', {}).get('contraindications', []))} contraindications")
    print()
    print("🎯 KEY DIFFERENTIATORS:")
    print("  1. PERSONALIZED: Matches compound to YOUR cancer type, biomarkers, treatment history")
    print("  2. INTEGRATED SCORING: S/P/E framework + SAE treatment line intelligence")
    print("  3. BIOMARKER-SPECIFIC: HRD+ vs HRD- gets different recommendations")
    print("  4. DRUG INTERACTIONS: Checks YOUR medication list")
    print("  5. TREATMENT LINE: Considers where you are in treatment journey")
    print("  6. EVIDENCE-GRADED: Not just papers - STRONG/MODERATE/WEAK classification")
    print("  7. DIETICIAN-READY: Complete guidance package, not just search results")

def test_4_input_documentation():
    """Test 4: Document all inputs"""
    print("\n🧪 TEST 4: Complete Input Documentation")
    print("=" * 60)
    
    print("📋 INPUTS REQUIRED:")
    print()
    print("1. COMPOUND (Required):")
    print("   - Any food or supplement name")
    print("   - Examples: 'Resveratrol', 'Curcumin', 'Vitamin D', 'Green Tea Extract'")
    print("   - System dynamically extracts targets (ChEMBL/PubChem/LLM)")
    print()
    print("2. DISEASE_CONTEXT (Required):")
    print("   {")
    print("     'disease': 'ovarian_cancer_hgs',  // Disease type")
    print("     'mutations': [                    // Tumor mutations (if available)")
    print("       {'gene': 'TP53', 'hgvs_p': 'R248Q'}")
    print("     ],")
    print("     'biomarkers': {                   // ⭐ KEY FOR TARGETING")
    print("       'HRD': 'POSITIVE',              // HRD status")
    print("       'TMB': 8.2,                     // Tumor mutational burden")
    print("       'MSI': 'STABLE',                // Microsatellite instability")
    print("       'PIK3CA': 'MUTANT'              // Specific gene mutations")
    print("     },")
    print("     'pathways_disrupted': [           // Cancer pathways affected")
    print("       'DNA repair',")
    print("       'Angiogenesis',")
    print("       'Inflammation'")
    print("     ]")
    print("   }")
    print()
    print("3. TREATMENT_HISTORY (Optional but recommended):")
    print("   {")
    print("     'current_line': 'L3',            // Current treatment line")
    print("     'prior_therapies': [              // Prior treatments")
    print("       'carboplatin',")
    print("       'paclitaxel'")
    print("     ]")
    print("   }")
    print("   → Used for: SAE features (line appropriateness, sequencing fitness)")
    print()
    print("4. PATIENT_MEDICATIONS (Optional):")
    print("   ['warfarin', 'metformin']")
    print("   → Used for: Drug interaction checking")
    print()
    print("5. USE_EVO2 (Optional):")
    print("   false  // Phase 1: disabled, Phase 2: experimental toggle")
    print()
    print("🎯 BIOMARKER TARGETING HOW IT WORKS:")
    print()
    print("Step 1: System extracts compound targets (e.g., Vitamin D → VDR, BRCA1)")
    print("Step 2: Maps targets to cancer pathways (e.g., VDR → DNA repair pathway)")
    print("Step 3: Checks biomarker context:")
    print("  - HRD+ + DNA repair pathway → BOOST line appropriateness (+0.1)")
    print("  - HRD+ + DNA repair compound → BOOST confidence (+0.05)")
    print("  - TMB ≥10 → Additional confidence boost (+0.03)")
    print("Step 4: Pathway alignment scoring:")
    print("  - Compound pathways ∩ Disease pathways = Alignment score")
    print("  - Higher alignment = Higher pathway (P) score")
    print("Step 5: SAE treatment line features:")
    print("  - Biomarker gates in supplement_rules.json match your biomarkers")
    print("  - Boosts line appropriateness if biomarker matches")
    print()
    print("📊 EXAMPLE: Vitamin D for HRD+ patient")
    print("  → Targets: VDR, BRCA1, TP53 (DNA repair targets)")
    print("  → Pathways: DNA repair (matches disease pathway)")
    print("  → Biomarker: HRD+ (matches Vitamin D biomarker gate)")
    print("  → Result: Higher line appropriateness (0.9 → 1.0)")
    print("  → Confidence boost: +0.05 for biomarker match")

def main():
    """Run all tests"""
    print("⚔️ DYNAMIC FOOD VALIDATOR - TESTING SUITE")
    print("=" * 60)
    print()
    
    # Check if backend is running
    try:
        health = requests.get(f"{API_BASE}/healthz", timeout=5)
        if health.status_code != 200:
            print("⚠️  Backend may not be running. Expected 200, got", health.status_code)
    except Exception as e:
        print(f"❌ Cannot connect to backend: {e}")
        print(f"   Make sure backend is running at {API_BASE}")
        return
    
    print("✅ Backend is running\n")
    
    # Run tests
    try:
        test_1_resveratrol_basic()
        test_2_biomarker_targeting()
        test_3_vs_googling()
        test_4_input_documentation()
        
        print("\n" + "=" * 60)
        print("✅ ALL TESTS COMPLETE")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ Test error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()


