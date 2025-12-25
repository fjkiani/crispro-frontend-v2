# 🥗 MODULE 06: TOXICITY-AWARE NUTRITION AGENT - PHASE 3 LLM ENHANCEMENT

**Status:** ✅ **PHASE 3 COMPLETE**  
**Date:** January 28, 2025  
**Agent:** Agent Jr  
**Connection:** Implements Phase 3 from `TOXICITY_MOAT_IMPLEMENTATION_TASKS.md`

---

## 📋 WHAT WAS BUILT

### Phase 3: LLM Enhancement for Toxicity Mitigation

**Goal:** Add personalized, context-aware LLM explanations to toxicity mitigation recommendations.

**Files Created/Modified:**

1. **`oncology-coPilot/oncology-backend-minimal/api/services/llm_toxicity_service.py`** ✅ NEW
   - `generate_toxicity_rationale()` - Personalized rationale generation
   - `generate_mitigation_dossier()` - Complete dossier with LLM summaries
   - Graceful fallback to static rationales if LLM unavailable

2. **`oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`** ✅ MODIFIED
   - Added LLM enhancement integration in `validate_food_dynamic` endpoint
   - Optional `enable_llm_enhancement` flag in request
   - Adds `llm_rationale`, `patient_summary`, `llm_enhanced` fields to `toxicity_mitigation`

3. **`oncology-coPilot/oncology-backend-minimal/test_llm_toxicity.py`** ✅ NEW
   - Comprehensive test suite for LLM enhancement
   - Edge case testing (empty genes, unknown drugs, special chars)
   - Integration tests

---

## 🎯 CAPABILITIES ADDED

### 1. Personalized Rationale Generation

**Before (Phase 2):**
```json
{
  "toxicity_mitigation": {
    "mitigates": true,
    "target_drug": "carboplatin",
    "mechanism": "Glutathione precursor, supports DNA repair enzymes"
  }
}
```

**After (Phase 3 with LLM):**
```json
{
  "toxicity_mitigation": {
    "mitigates": true,
    "target_drug": "carboplatin",
    "mechanism": "Glutathione precursor, supports DNA repair enzymes",
    "llm_rationale": "NAC supports DNA repair by boosting glutathione levels, which are depleted by carboplatin's platinum-based mechanism. For BRCA1 carriers, this is particularly important as...",
    "patient_summary": "Carboplatin works by damaging cancer cell DNA, but it can also stress your body's DNA repair systems. NAC is a supplement that may help protect your cells during treatment by supporting these repair mechanisms.",
    "llm_enhanced": true,
    "llm_confidence": 0.75
  }
}
```

### 2. Patient-Friendly Explanations

- **Clinician View**: Technical rationale with pathway connections
- **Patient View**: Simple 2-sentence explanation (8th grade reading level)
- **Context-Aware**: Incorporates cancer type, treatment phase, germline variants

### 3. Complete Mitigation Dossiers

Generates structured documents with:
- Executive summary (LLM-generated)
- Food recommendations with personalized rationales
- Timing protocols
- Monitoring recommendations

---

## 🔗 INTEGRATION WITH MODULE 06

### Original Module 06 Spec (Basic)

The original `06_NUTRITION_AGENT.mdc` specified:
- Drug toxicity mapping
- Protective supplement recommendations
- Food-drug interaction checking
- Timing rule generation
- Germline-aware adjustments

### Phase 3 Enhancement

**What Phase 3 Adds:**
- ✅ **LLM-Powered Personalization**: Rationales tailored to patient context
- ✅ **Multi-Audience Explanations**: Clinician + Patient + Researcher views
- ✅ **Evidence Synthesis**: LLM summarizes multiple foods into coherent plan
- ✅ **Graceful Degradation**: Falls back to static rationales if LLM unavailable

**Connection Points:**
- Uses existing `get_mitigating_foods()` from Phase 1+2
- Enhances `toxicity_mitigation` field already in food validation response
- Integrates with existing `validate_food_dynamic` endpoint
- No breaking changes - LLM enhancement is optional

---

## 📊 ARCHITECTURE

```
┌─────────────────────────────────────────────────────────────┐
│  validate_food_dynamic (hypothesis_validator.py)           │
│                                                              │
│  [Phase 2] Toxicity Mitigation Check                        │
│  ├── compute_pathway_overlap()                              │
│  ├── get_mitigating_foods()                                  │
│  └── toxicity_mitigation = {...}                            │
│                                                              │
│  [Phase 3] LLM Enhancement (if enabled)                     │
│  ├── get_llm_toxicity_service()                             │
│  ├── generate_toxicity_rationale()                          │
│  │   ├── query_llm() [Gemini API]                           │
│  │   └── Returns: rationale + patient_summary               │
│  └── toxicity_mitigation["llm_rationale"] = ...              │
└─────────────────────────────────────────────────────────────┘
```

---

## 🧪 TESTING

### Test File: `test_llm_toxicity.py`

**Test Coverage:**
- ✅ Single rationale generation (NAC + carboplatin + BRCA1)
- ✅ Cardiotoxicity scenario (CoQ10 + doxorubicin)
- ✅ Full dossier generation
- ✅ Edge cases (empty genes, unknown drugs, special chars)

**Run Tests:**
```bash
cd oncology-coPilot/oncology-backend-minimal
python3 test_llm_toxicity.py
```

**Pre-requisite:** `GEMINI_API_KEY` must be set in `.env` (falls back to static if missing)

---

## 📝 USAGE

### API Request

```json
POST /api/hypothesis/validate_food_dynamic
{
  "compound": "NAC",
  "disease_context": {
    "disease": "ovarian_cancer",
    "mutations": [{"gene": "BRCA1"}]
  },
  "patient_medications": ["carboplatin"],
  "enable_llm_enhancement": true  // NEW FLAG
}
```

### Response (Enhanced)

```json
{
  "toxicity_mitigation": {
    "mitigates": true,
    "target_drug": "carboplatin",
    "target_moa": "platinum_agent",
    "pathway": "dna_repair",
    "mechanism": "Glutathione precursor, supports DNA repair enzymes",
    "timing": "post-chemo (not during infusion)",
    "dose": "600mg twice daily",
    "llm_rationale": "NAC supports DNA repair by...",  // NEW
    "patient_summary": "Carboplatin works by...",      // NEW
    "llm_enhanced": true,                               // NEW
    "llm_confidence": 0.75                              // NEW
  }
}
```

---

## ✅ ACCEPTANCE CRITERIA (MET)

1. ✅ LLM service created with proper path resolution
2. ✅ Integration into `validate_food_dynamic` endpoint
3. ✅ Graceful fallback if LLM unavailable
4. ✅ Patient-friendly summaries generated
5. ✅ Edge cases handled (empty genes, unknown drugs)
6. ✅ Test suite created and passing
7. ✅ No breaking changes to existing API

---

## 🔄 NEXT STEPS (Future Phases)

### Phase 4: Dossier Generation (Planned)
- Generate complete markdown dossiers
- Export functionality
- Integration with frontend

### Phase 5: Frontend Integration (Planned)
- Display LLM-enhanced rationales in FoodRankingPanel
- Toggle between clinician/patient views
- Export dossier functionality

---

## 📚 REFERENCES

- **Implementation Plan**: `.cursor/plans/TOXICITY_MOAT_IMPLEMENTATION_TASKS.md`
- **Module 06 Spec**: `.cursor/MOAT/orchestration/06_NUTRITION_AGENT.mdc`
- **Master Index**: `.cursor/MOAT/orchestration/00_MASTER_INDEX.mdc`
- **Advanced Care Plan**: `.cursor/MOAT/ADVANCED_CARE_PLAN_EXPLAINED.md`

---

## 👤 OWNERSHIP

**Built By:** Agent Jr  
**Date:** January 28, 2025  
**Status:** ✅ Phase 3 Complete - Ready for Phase 4  
**Files Owned:**
- `api/services/llm_toxicity_service.py`
- `test_llm_toxicity.py`
- Integration in `api/routers/hypothesis_validator.py` (lines ~949-990)

---

**Last Updated:** January 28, 2025


