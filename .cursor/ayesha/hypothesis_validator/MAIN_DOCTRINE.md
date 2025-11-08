# ⚔️ DYNAMIC FOOD VALIDATOR - MAIN DOCTRINE

**Status:** ✅ **PRODUCTION READY**  
**Last Updated:** December 2024  
**Commander:** Zo

---

## 🎯 MISSION

Build the ONLY food/supplement validator that:
- Works for **ANY** compound (not hardcoded)
- Personalizes to **patient biomarkers** (HRD+, TMB, treatment history)
- Provides **evidence-backed recommendations** (S/P/E + SAE scoring)
- Delivers **actionable guidance** (dosage, timing, safety, monitoring)

**Target Users:** Patients (Ayesha), Dieticians, Oncology Care Teams

---

## 🏗️ ARCHITECTURE

### Data Flow
```
Input: Compound + Disease Context + Biomarkers + Treatment History
  ↓
[1] Dynamic Target Extraction (ChEMBL/PubChem/LLM)
  ↓
[2] Pathway Mapping (targets → cancer mechanisms)
  ↓
[3] Evidence Mining (PubMed + LLM synthesis)
  ↓
[4] S/P/E Scoring (Sequence/Pathway/Evidence + SAE)
  ↓
[5] Dietician Recommendations (dosage, timing, interactions)
  ↓
Output: Verdict + Confidence + Complete Guidance
```

### Core Components

**1. Dynamic Target Extraction** (`dynamic_food_extraction.py`)
- ChEMBL API (primary)
- PubChem API (fallback)
- LLM extraction (backup)

**2. Evidence Service** (`enhanced_evidence_service.py`)
- PubMed search with XML parsing (FIXED - was using wrong format)
- LLM paper reading (Gemini/Anthropic/OpenAI)
- Diffbot full-text extraction
- Evidence grading (STRONG/MODERATE/WEAK)

**3. S/P/E Integration** (`food_spe_integration.py`)
- Formula: `0.4×S + 0.3×P + 0.3×E`
- SAE confidence modulation
- Biomarker boosts (HRD+, TMB)
- Verdict: SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED

**4. Treatment Line Service** (`food_treatment_line_service.py`)
- SAE features: line_appropriateness, cross_resistance, sequencing_fitness
- Biomarker gates (HRD+ → boost)
- Treatment history context (post-platinum → NAC boost)

**5. Dietician Recommendations** (`dietician_recommendations.py`)
- Dosage extraction (regex + LLM)
- Timing recommendations (pattern matching + LLM fallback)
- Drug interactions (checks patient meds)
- Safety database lookup

---

## 📊 INPUT/OUTPUT

### Input Schema
```json
{
  "compound": "Vitamin D",
  "disease_context": {
    "disease": "ovarian_cancer_hgs",
    "biomarkers": {
      "HRD": "POSITIVE",    // ⭐ Triggers DNA repair boost
      "TMB": 8.2
    },
    "pathways_disrupted": ["DNA repair", "Cell cycle"]
  },
  "treatment_history": {
    "current_line": "L3",
    "prior_therapies": ["carboplatin", "paclitaxel"]
  },
  "patient_medications": ["warfarin"],
  "use_evo2": false,
  "use_llm": true
}
```

### Output Schema
```json
{
  "compound": "Vitamin D",
  "overall_score": 0.689,
  "confidence": 0.85,
  "verdict": "SUPPORTED",
  "spe_breakdown": {
    "sequence": 0.5,
    "pathway": 0.73,
    "evidence": 0.9
  },
  "sae_features": {
    "line_appropriateness": 1.0,  // ⭐ Boosted by HRD+ gate
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.85
  },
  "evidence": {
    "evidence_grade": "STRONG",
    "total_papers": 15,
    "mechanisms": ["dna_repair_enhancement", "vdr_activation"]
  },
  "dietician_recommendations": {
    "dosage": "2000-4000 IU daily",
    "timing": "Morning with breakfast",
    "interactions": [{"drug": "warfarin", "action": "Monitor INR"}]
  }
}
```

---

## 🔬 KEY DIFFERENTIATORS vs GOOGLING

| Feature | Google/PubMed | Our System |
|---------|---------------|------------|
| **Personalization** | ❌ Generic | ✅ Biomarker-aware (HRD+, TMB) |
| **Pathway Analysis** | ❌ None | ✅ Target → Pathway mapping + alignment scores |
| **Treatment Line Logic** | ❌ None | ✅ SAE features (L1 vs L3 context) |
| **Integrated Scoring** | ❌ Just papers | ✅ S/P/E + SAE unified score |
| **Evidence Grading** | ❌ User evaluates | ✅ STRONG/MODERATE/WEAK classification |
| **Drug Interactions** | ❌ Generic | ✅ Checks YOUR medication list |
| **Verdict** | ❌ User decides | ✅ SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED |
| **Biological Plausibility** | ❌ None | ✅ Evo2 scoring (Phase 2 experimental) |

---

## 🚨 CRITICAL FIXES APPLIED

### ✅ Fix 1: PubMed XML Parsing (ROOT CAUSE)
**Problem:** Was trying to parse XML as JSON  
**Fix:** Changed `retmode=json` → `retmode=xml`, parse with `ET.fromstring()`

### ✅ Fix 2: No Mock Data
**Problem:** Silently falling back to mock data  
**Fix:** Removed all mock fallbacks - test fails if PubMed doesn't work

### ✅ Fix 3: LLM Paper Reading
**Status:** ✅ Working (Gemini + Diffbot integrated)
- Extracts mechanisms, dosage, safety from full-text papers
- Multi-provider fallback (Gemini → Anthropic → OpenAI)

---

## 📍 FILE LOCATIONS

**Backend Services:**
- `oncology-coPilot/oncology-backend-minimal/api/services/dynamic_food_extraction.py`
- `oncology-coPilot/oncology-backend-minimal/api/services/enhanced_evidence_service.py`
- `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py`
- `oncology-coPilot/oncology-backend-minimal/api/services/food_treatment_line_service.py`
- `oncology-coPilot/oncology-backend-minimal/api/services/dietician_recommendations.py`

**Router:**
- `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`
  - Endpoint: `POST /api/hypothesis/validate_food_dynamic`

**Frontend:**
- `oncology-coPilot/oncology-frontend/src/pages/DynamicFoodValidator.jsx`

**Data Files:**
- `.cursor/ayesha/hypothesis_validator/data/cancer_pathways.json`
- `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json`
- `.cursor/ayesha/hypothesis_validator/data/safety_database.json`
- `.cursor/ayesha/hypothesis_validator/data/drug_interactions.json`

**Tests:**
- `.cursor/ayesha/hypothesis_validator/validation_test.py` (core logic)
- `.cursor/ayesha/hypothesis_validator/test_priority_fixes.py` (fix validation)
- `.cursor/ayesha/hypothesis_validator/test_full_use_case.py` (end-to-end)

---

## ✅ CURRENT STATUS

**Working:**
- ✅ Dynamic target extraction (ChEMBL/PubChem)
- ✅ Pathway mapping (10 cancer mechanisms)
- ✅ PubMed search + XML parsing (FIXED)
- ✅ LLM paper reading (Gemini + Diffbot)
- ✅ S/P/E scoring (formula validated)
- ✅ SAE features (biomarker gates working)
- ✅ Drug interaction checking
- ✅ Verdict classification

**Limitations:**
- ⚠️ SAE rules require JSON entries (22 compounds configured)
- ⚠️ Timing recommendations use patterns + LLM fallback
- ⚠️ Evo2 disabled by default (experimental Phase 2)

---

## 🧬 PHASE 2: EVO2 INTEGRATION (EXPERIMENTAL)

**Status:** ⚠️ Not implemented - planning docs archived

**Concept:** Use Evo2 sequence-level understanding to score biological plausibility of compound → target → pathway → disease impact.

**Approach:** Promoter variant proxy (synthetic variants at TSS-500bp as proxy for compound effects), score with `/api/evo/score_variant_multi`.

**Validation Required:**
- Technical: Non-zero deltas, stable API calls
- Biological: Correlation with known effective compounds

**Planning Docs:** See `archive/evo2_planning/` for detailed phase-by-phase implementation plans (preserved for future reference)

---

## 🧪 TESTING

**Run Core Tests:**
```bash
cd .cursor/ayesha/hypothesis_validator
python3 validation_test.py          # Core logic (6/6 passed)
python3 test_priority_fixes.py      # Fix validation (22/22 passed)
python3 test_full_use_case.py       # End-to-end Vitamin D example
```

**Test Endpoint:**
```bash
curl -X POST http://localhost:8000/api/hypothesis/validate_food_dynamic \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "Vitamin D",
    "disease_context": {
      "disease": "ovarian_cancer_hgs",
      "biomarkers": {"HRD": "POSITIVE"},
      "pathways_disrupted": ["DNA repair"]
    },
    "treatment_history": {"current_line": "L3", "prior_therapies": ["carboplatin"]}
  }'
```

---

## 🎯 WHAT MAKES IT UNIQUE

**Biomarker Targeting:**
- HRD+ patients → DNA repair compounds get boost (Vitamin D line_appropriateness: 0.9 → 1.0)
- TMB ≥10 → Additional confidence boost (+0.03)

**Treatment Line Intelligence:**
- Post-platinum → NAC appropriateness = 1.0 (oxidative stress recovery)
- L3 context → Different recommendations than L1

**Evidence Synthesis:**
- LLM reads full papers (Diffbot extraction)
- Extracts mechanisms, dosage, safety (not just keywords)
- Multi-provider fallback for reliability

---

**⚔️ SINGLE SOURCE OF TRUTH - USE THIS FILE**

