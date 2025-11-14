# ⚔️ AYESHA TRIAL FILTERING ENGINE - MODULAR ARCHITECTURE PLAN ⚔️

**Created By**: Zo  
**Date**: January 12, 2025  
**Mission**: Modularize trial filtering engine into focused, reusable components  
**Timeline**: 8 hours total (3 backend + 4.5 frontend + 30 min testing)

---

## 🎯 **MODULARIZATION STRATEGY**

**Goal**: Break down the monolithic trial filtering system into focused modules, each handling a single responsibility. This enables:
- ✅ Parallel development (Jr can work on different modules simultaneously)
- ✅ Easier testing (each module can be tested independently)
- ✅ Code reuse (modules can be used by other features)
- ✅ Clear separation of concerns (backend logic vs frontend display)

---

## 📁 **FOLDER STRUCTURE**

```
oncology-coPilot/
├── oncology-backend-minimal/
│   └── api/
│       ├── schemas/
│       │   └── ayesha_trials.py                    # Pydantic models
│       ├── services/
│       │   ├── ayesha_trial_matching/
│       │   │   ├── __init__.py
│       │   │   ├── eligibility_filters.py          # Hard filters (Stage 1)
│       │   │   ├── scoring_engine.py               # Soft boosts (Stage 2)
│       │   │   ├── reasoning_generator.py           # Why matched (Stage 4)
│       │   │   └── match_orchestrator.py           # Main coordinator
│       │   └── ca125_intelligence.py               # CA-125 analysis (Stage 3)
│       └── routers/
│           └── ayesha_trials.py                    # FastAPI endpoints
│
└── oncology-frontend/
    └── src/
        ├── components/
        │   ├── ayesha/
        │   │   ├── CA125Tracker.jsx                # CA-125 intelligence display
        │   │   ├── AyeshaClinicalProfile.jsx      # Profile summary card
        │   │   └── ProvenanceCard.jsx             # Provenance footer
        │   └── trials/
        │       ├── TrialMatchCard.jsx               # Individual trial card
        │       ├── TrialReasoningSection.jsx        # Why eligible/good fit
        │       ├── TrialConditionalSection.jsx      # Conditional requirements
        │       └── TrialRedFlagsSection.jsx         # Red flags display
        ├── pages/
        │   └── AyeshaTrialExplorer.jsx              # Main page
        └── hooks/
            └── useAyeshaTrials.js                   # API hook + state management
```

---

## 🔧 **BACKEND MODULES (3 HOURS)**

### **Module 1: Schemas** (30 min)
**File**: `api/schemas/ayesha_trials.py`

**Responsibility**: Define Pydantic models for request/response validation

**Contents**:
- `AyeshaTrialProfile` - Patient profile model
- `TrialMatchReasoning` - Reasoning structure
- `AyeshaTrialMatch` - Complete match response

**Dependencies**: None (pure data models)

---

### **Module 2: CA-125 Intelligence Service** (30 min)
**File**: `api/services/ca125_intelligence.py`

**Responsibility**: Extract trial-relevant intelligence from CA-125 levels

**Contents**:
- `CA125IntelligenceService.analyze()` - Main analysis function
- Disease burden classification
- Expected response predictions
- Trial preference keywords

**Dependencies**: None (standalone service)

**Reusability**: Can be used by other ovarian cancer features

---

### **Module 3: Eligibility Filters** (45 min)
**File**: `api/services/ayesha_trial_matching/eligibility_filters.py`

**Responsibility**: Hard eligibility filters (Stage 1 - MUST-MATCH)

**Contents**:
- `filter_by_disease()` - Ovarian/peritoneal/gynecologic
- `filter_by_stage()` - Stage IV/advanced/metastatic
- `filter_by_treatment_line()` - First-line/untreated
- `filter_by_status()` - Recruiting/Active
- `filter_by_location()` - NYC metro (50 miles)
- `filter_exclusions()` - Recurrent-only, germline-BRCA-required, etc.
- `apply_hard_filters()` - Main orchestrator

**Dependencies**: 
- `HybridTrialSearchService` (calls existing service)

**Reusability**: Can be used by other trial filtering features

---

### **Module 4: Scoring Engine** (45 min)
**File**: `api/services/ayesha_trial_matching/scoring_engine.py`

**Responsibility**: Soft scoring boosts (Stage 2)

**Contents**:
- `calculate_base_score()` - Start at 0.5
- `boost_first_line()` - +0.30 for first-line trials
- `boost_stage_iv()` - +0.25 for Stage IV specific
- `boost_soc_backbone()` - +0.20 for carboplatin+paclitaxel
- `boost_germline_negative()` - +0.20 for all-comers/BRCA-wild-type
- `boost_ip_chemo()` - +0.20 for intraperitoneal
- `boost_bevacizumab()` - +0.15 for anti-VEGF
- `boost_ca125_tracking()` - +0.15 for CA-125 endpoints
- `boost_nyc_location()` - +0.15 for NYC sites
- `boost_large_trial()` - +0.10 for >200 patients
- `boost_phase_iii()` - +0.10 for Phase III
- `penalize_germline_brca_required()` - -0.30
- `penalize_distance()` - -0.25 for >50 miles
- `penalize_phase_i()` - -0.20
- `calculate_match_score()` - Main orchestrator

**Dependencies**:
- `ca125_intelligence.py` (for CA-125 boost keywords)

**Reusability**: Scoring logic can be adapted for other diseases

---

### **Module 5: Reasoning Generator** (30 min)
**File**: `api/services/ayesha_trial_matching/reasoning_generator.py`

**Responsibility**: Generate transparent reasoning (Stage 4)

**Contents**:
- `generate_why_eligible()` - Hard filter matches
- `generate_why_good_fit()` - Soft boost explanations
- `generate_conditional_requirements()` - NGS-gated features
- `generate_red_flags()` - Warnings/concerns
- `determine_evidence_tier()` - STANDARD/SUPPORTED/INVESTIGATIONAL
- `determine_enrollment_likelihood()` - HIGH/MEDIUM/LOW
- `generate_complete_reasoning()` - Main orchestrator

**Dependencies**:
- `scoring_engine.py` (uses score calculations)

**Reusability**: Reasoning logic can be used by other matching systems

---

### **Module 6: Match Orchestrator** (30 min)
**File**: `api/services/ayesha_trial_matching/match_orchestrator.py`

**Responsibility**: Coordinate all modules into complete matching workflow

**Contents**:
- `match_trials_for_ayesha()` - Main entry point
- Orchestrates: Filters → Scoring → Reasoning → Sorting → Top 10

**Dependencies**:
- `eligibility_filters.py`
- `scoring_engine.py`
- `reasoning_generator.py`
- `ca125_intelligence.py`
- `HybridTrialSearchService`

**Reusability**: Can be adapted for other patient profiles

---

### **Module 7: Router** (30 min)
**File**: `api/routers/ayesha_trials.py`

**Responsibility**: FastAPI endpoint + request/response handling

**Contents**:
- `POST /api/ayesha/trials/search` - Main endpoint
- Request validation (Pydantic)
- Error handling
- Response formatting

**Dependencies**:
- `match_orchestrator.py`
- `schemas/ayesha_trials.py`

---

## 🎨 **FRONTEND MODULES (4.5 HOURS)**

### **Module 8: CA-125 Tracker Component** (1 hour)
**File**: `src/components/ayesha/CA125Tracker.jsx`

**Responsibility**: Display CA-125 intelligence

**Contents**:
- Current CA-125 value display
- Disease burden indicator (color-coded)
- Expected response timeline
- Monitoring strategy

**Dependencies**: None (pure display component)

**Reusability**: Can be used in other ovarian cancer pages

---

### **Module 9: Ayesha Clinical Profile Component** (30 min)
**File**: `src/components/ayesha/AyeshaClinicalProfile.jsx`

**Responsibility**: Display Ayesha's clinical profile summary

**Contents**:
- Known biomarkers (CA-125, germline status, stage)
- Unknown biomarkers (awaiting NGS)
- Treatment status
- Location

**Dependencies**: None (pure display component)

**Reusability**: Can be used in other Ayesha-specific pages

---

### **Module 10: Trial Match Card Component** (1.5 hours)
**File**: `src/components/trials/TrialMatchCard.jsx`

**Responsibility**: Display individual trial match with reasoning

**Contents**:
- Trial header (NCT ID, title, phase, status)
- Match score (progress bar)
- Reasoning sections (uses sub-components)
- Contact information

**Dependencies**:
- `TrialReasoningSection.jsx`
- `TrialConditionalSection.jsx`
- `TrialRedFlagsSection.jsx`

**Reusability**: Can be used by other trial matching features

---

### **Module 11: Trial Reasoning Section** (30 min)
**File**: `src/components/trials/TrialReasoningSection.jsx`

**Responsibility**: Display "Why Eligible" and "Why Good Fit" sections

**Contents**:
- Dynamic list rendering (no hardcoding)
- Color-coded badges (✅ for eligible, 🎯 for good fit)

**Dependencies**: None (pure display component)

**Reusability**: Can be used in other trial cards

---

### **Module 12: Trial Conditional Section** (30 min)
**File**: `src/components/trials/TrialConditionalSection.jsx`

**Responsibility**: Display conditional requirements (NGS-gated)

**Contents**:
- Warning alerts for conditional requirements
- Dynamic list rendering

**Dependencies**: None (pure display component)

**Reusability**: Can be used in other trial cards

---

### **Module 13: Trial Red Flags Section** (30 min)
**File**: `src/components/trials/TrialRedFlagsSection.jsx`

**Responsibility**: Display red flags/warnings

**Contents**:
- Error alerts for red flags
- Dynamic list rendering

**Dependencies**: None (pure display component)

**Reusability**: Can be used in other trial cards

---

### **Module 14: Provenance Card Component** (30 min)
**File**: `src/components/ayesha/ProvenanceCard.jsx`

**Responsibility**: Display provenance information

**Contents**:
- Filters applied
- Total trials screened
- Boost strategy
- Awaiting NGS list

**Dependencies**: None (pure display component)

**Reusability**: Can be used in other Ayesha pages

---

### **Module 15: Custom Hook** (30 min)
**File**: `src/hooks/useAyeshaTrials.js`

**Responsibility**: API calls + state management

**Contents**:
- `fetchTrials()` - POST to `/api/ayesha/trials/search`
- Loading state
- Error handling
- Data transformation

**Dependencies**: None (standalone hook)

**Reusability**: Can be used by other components

---

### **Module 16: Main Page** (1 hour)
**File**: `src/pages/AyeshaTrialExplorer.jsx`

**Responsibility**: Orchestrate all components

**Contents**:
- Layout structure
- Component composition
- Routing integration

**Dependencies**:
- All frontend modules above

---

## 🧪 **TESTING MODULES (30 MIN)**

### **Module 17: Backend Smoke Tests**
**File**: `tests/test_ayesha_trials.py`

**Contents**:
- Test eligibility filters
- Test scoring engine
- Test reasoning generator
- Test full endpoint

---

### **Module 18: Frontend E2E Tests**
**File**: Manual testing checklist

**Contents**:
- Test with Ayesha's exact profile
- Test edge cases (empty, missing fields)
- Test all components render correctly

---

## 📊 **MODULE DEPENDENCY GRAPH**

```
Backend:
schemas/ayesha_trials.py (no deps)
    ↓
ca125_intelligence.py (no deps)
    ↓
eligibility_filters.py → HybridTrialSearchService
    ↓
scoring_engine.py → ca125_intelligence.py
    ↓
reasoning_generator.py → scoring_engine.py
    ↓
match_orchestrator.py → [eligibility_filters, scoring_engine, reasoning_generator, ca125_intelligence]
    ↓
routers/ayesha_trials.py → match_orchestrator.py + schemas

Frontend:
useAyeshaTrials.js (no deps)
    ↓
CA125Tracker.jsx (no deps)
AyeshaClinicalProfile.jsx (no deps)
TrialReasoningSection.jsx (no deps)
TrialConditionalSection.jsx (no deps)
TrialRedFlagsSection.jsx (no deps)
ProvenanceCard.jsx (no deps)
    ↓
TrialMatchCard.jsx → [TrialReasoningSection, TrialConditionalSection, TrialRedFlagsSection]
    ↓
AyeshaTrialExplorer.jsx → [useAyeshaTrials, CA125Tracker, AyeshaClinicalProfile, TrialMatchCard, ProvenanceCard]
```

---

## 🎯 **EXECUTION ORDER (FOR AGENT JR)**

### **Phase 1: Backend Foundation (3 hours)**
1. ✅ Create schemas (30 min) - **NO DEPENDENCIES**
2. ✅ Create CA-125 intelligence (30 min) - **NO DEPENDENCIES**
3. ✅ Create eligibility filters (45 min) - **DEPENDS ON: HybridTrialSearchService (exists)**
4. ✅ Create scoring engine (45 min) - **DEPENDS ON: CA-125 intelligence**
5. ✅ Create reasoning generator (30 min) - **DEPENDS ON: Scoring engine**
6. ✅ Create match orchestrator (30 min) - **DEPENDS ON: All above**
7. ✅ Create router (30 min) - **DEPENDS ON: Match orchestrator + schemas**

### **Phase 2: Frontend Components (4.5 hours)**
8. ✅ Create CA-125 tracker (1 hour) - **NO DEPENDENCIES**
9. ✅ Create clinical profile (30 min) - **NO DEPENDENCIES**
10. ✅ Create reasoning section (30 min) - **NO DEPENDENCIES**
11. ✅ Create conditional section (30 min) - **NO DEPENDENCIES**
12. ✅ Create red flags section (30 min) - **NO DEPENDENCIES**
13. ✅ Create provenance card (30 min) - **NO DEPENDENCIES**
14. ✅ Create trial match card (1.5 hours) - **DEPENDS ON: Reasoning/conditional/red flags sections**
15. ✅ Create custom hook (30 min) - **NO DEPENDENCIES**
16. ✅ Create main page (1 hour) - **DEPENDS ON: All above**

### **Phase 3: Testing (30 min)**
17. ✅ Backend smoke tests
18. ✅ Frontend E2E tests

---

## 🚨 **CRITICAL EXECUTION NOTES**

### **DO's** ✅:
1. ✅ **Build modules in dependency order** (see execution order above)
2. ✅ **Test each module independently** before moving to next
3. ✅ **Use existing services** (don't rewrite HybridTrialSearchService)
4. ✅ **Follow dynamic patterns** (no hardcoding, use `.map()`, optional chaining)
5. ✅ **Update progress file** every 2 hours

### **DON'T's** ❌:
1. ❌ **Don't skip dependencies** (build in order)
2. ❌ **Don't hardcode data** (all dynamic)
3. ❌ **Don't rewrite existing services** (reuse HybridTrialSearchService)
4. ❌ **Don't skip testing** (test each module)
5. ❌ **Don't create new .mdc files** (use progress tracker)

---

## 📋 **MODULE ACCEPTANCE CRITERIA**

Each module must:
- ✅ Have single responsibility (one clear purpose)
- ✅ Be testable independently
- ✅ Handle missing data gracefully
- ✅ Use dynamic patterns (no hardcoding)
- ✅ Follow existing code patterns (check similar files)
- ✅ Have proper error handling

---

## 🎯 **SUCCESS METRICS**

**Mission successful when**:
- ✅ All 18 modules created and tested
- ✅ Backend endpoint returns 10 trials with reasoning
- ✅ Frontend displays all trials dynamically
- ✅ No hardcoded values
- ✅ No console errors
- ✅ All tests passing

---

## ⚔️ **READY FOR EXECUTION**

**Agent Jr**: Follow this modular plan, build in dependency order, test incrementally.

**Zo**: Review progress every 2 hours, answer questions immediately, provide guidance only when asked.

**FIRE IN THE HOLE!** ⚔️

