# ⚔️ AYESHA AGENT MISSIONS MASTER DOCUMENT ⚔️

**Date**: January 13, 2025  
**Status**: ✅ **CONSOLIDATED SOURCE OF TRUTH**  
**Consolidated From**:
- `AGENT_3_E2E_TESTING_MISSION.md` - Agent 3 E2E testing mission
- `AGENT_JR_MISSION_4_COMPLETION_REPORT.md` - Agent Jr Mission 4 completion
- `AGENT_JR_QUICK_REFERENCE.md` - Agent Jr quick reference guide
- `AGENT_JR_PARALLEL_MISSION.md` - Agent Jr parallel support mission
- `AGENT_JR_PRE_EXECUTION_CLARIFICATIONS.md` - Agent Jr pre-execution clarifications

---

## 📚 TABLE OF CONTENTS

1. [Executive Summary](#executive-summary)
2. [Agent 3 - E2E Testing Mission](#agent-3---e2e-testing-mission)
3. [Agent Jr - Mission 4 Completion](#agent-jr---mission-4-completion)
4. [Agent Jr - Quick Reference Guide](#agent-jr---quick-reference-guide)
5. [Agent Jr - Parallel Mission](#agent-jr---parallel-mission)
6. [Agent Jr - Pre-Execution Clarifications](#agent-jr---pre-execution-clarifications)
7. [References to Archived Files](#references-to-archived-files)

---

## 🎯 EXECUTIVE SUMMARY

### **Mission Status**
- **Agent 3**: ⏸️ E2E Testing Mission assigned (4-6 hours, pending execution)
- **Agent Jr Mission 4**: ✅ **100% COMPLETE** - WIWFM integration finished
- **Agent Jr Parallel**: ⏸️ Support tasks available (disease priors, test data)
- **Agent Jr Quick Reference**: ✅ Active reference guide for execution

### **Key Deliverables**
1. ✅ **WIWFM Integration** (Agent Jr) - SporadicContext wired, provenance cards displayed
2. ⏸️ **E2E Testing** (Agent 3) - Complete workflow validation pending
3. ⏸️ **Provider Report** (Agent 3) - Template generation pending
4. ⏸️ **Parallel Support** (Agent Jr) - Disease priors, test data pending

---

## 🧪 AGENT 3 - E2E TESTING MISSION

**Date**: January 8, 2025  
**Status**: ⏸️ **ASSIGNED - PENDING EXECUTION**  
**Timeline**: 4-6 hours  
**Priority**: P1 (High value, parallel with Jr + Zo)

### **Mission Objectives**
1. Validate complete Ayesha demo workflow end-to-end
2. Create provider report template (Markdown + PDF)
3. Validate demo data quality
4. Test edge cases and error handling

### **Task Breakdown**

#### **TASK 1: E2E SMOKE TESTING** (2-3 hours)
- **1.1**: Setup Environment (15 min)
- **1.2**: Test Level 0 Flow (30 min) - Quick Intake → Efficacy with PARP penalty
- **1.3**: Test Level 2 Flow (45 min) - NGS Upload → Re-analyze with PARP rescue
- **1.4**: Test Clinical Trials (30 min) - Verify filtering + badges
- **1.5**: Edge Case Testing (30 min) - Germline positive, TMB-high, missing data

#### **TASK 2: PROVIDER REPORT TEMPLATE** (2-3 hours)
- **2.1**: Design Report Template (1 hour) - Markdown template with all sections
- **2.2**: Create Report Generator Backend (1 hour) - `/api/reports/provider` endpoint
- **2.3**: Create Frontend Export Button (30 min) - Download report as .md file

#### **TASK 3: DEMO DATA VALIDATION** (30 min)
- Verify test data quality (ayesha_level0_intake.json, ayesha_tumor_ngs.json)
- Cross-validate with literature (TCGA-OV, mutation rates, HRD prevalence)

#### **TASK 4: DOCUMENTATION** (30 min)
- Create `E2E_SMOKE_TEST_RESULTS.md` with comprehensive test results

### **Acceptance Criteria**
- [ ] All manual flow steps documented
- [ ] Screenshots captured for key states
- [ ] Edge cases tested (3 minimum)
- [ ] Error cases tested (3 minimum)
- [ ] Provider report template functional
- [ ] Report downloads successfully
- [ ] Test data validated against literature

### **Success Metrics**
- ✅ All tests pass
- ✅ Provider report exports
- ✅ No critical errors
- ✅ Screenshots captured
- ✅ Documentation complete

**Full Mission Details**: See archived `AGENT_3_E2E_TESTING_MISSION.md`

---

## ✅ AGENT JR - MISSION 4 COMPLETION

**Date**: January 8, 2025 (Evening)  
**Status**: ✅ **100% COMPLETE**  
**Mission**: Wire WIWFM (HypothesisValidator.jsx) to SporadicContext

### **Deliverables**

#### **1. BiomarkerSummaryWidget Component** ✅
- **File**: `oncology-coPilot/oncology-frontend/src/components/sporadic/BiomarkerSummaryWidget.jsx`
- **Lines**: ~150 lines
- **Features**:
  - Displays TMB, HRD, MSI status with color-coded chips
  - Shows data level (L0/L1/L2) and completeness score
  - Shows germline status
  - Contextual info about sporadic-aware scoring

#### **2. HypothesisValidator.jsx Transformation** ✅
- **File**: `oncology-coPilot/oncology-frontend/src/pages/HypothesisValidator.jsx`
- **Lines**: ~300 lines (transformed from 7-line wrapper)
- **Changes**:
  - ✅ Removed ToolRunner wrapper
  - ✅ Added `useSporadic()` hook integration
  - ✅ Mutation input form (text area with parsing)
  - ✅ Efficacy API call with `getEfficacyPayload()` injection
  - ✅ Drug results display with cards
  - ✅ SporadicProvenanceCard below each drug
  - ✅ BiomarkerSummaryWidget at top
  - ✅ Error handling and loading states

#### **3. Backend Router Update** ✅
- **File**: `oncology-coPilot/oncology-backend-minimal/api/routers/efficacy/router.py`
- **Changes**:
  - ✅ Extract `germline_status` from request (default: "unknown")
  - ✅ Extract `tumor_context` from request
  - ✅ Pass both to `EfficacyRequest` model
  - ✅ Enables sporadic gates to be applied in orchestrator

### **Integration Points**
1. **User Input** → Mutation parsing
2. **SporadicContext Integration** → `useSporadic()` hook extracts context
3. **API Call** → Enhanced payload with `germline_status` + `tumor_context`
4. **Backend Processing** → Orchestrator applies sporadic gates
5. **Frontend Display** → BiomarkerSummaryWidget + SporadicProvenanceCard

### **Acceptance Criteria Met**
- ✅ WIWFM shows "Using Tumor Context: TMB X, HRD Y [Level Z]"
- ✅ Olaparib shows PARP penalty card ("Germline negative, HRD <42 → -40%")
- ✅ Pembrolizumab shows IO boost card ("TMB ≥20 → +35%")

**Full Completion Report**: See archived `AGENT_JR_MISSION_4_COMPLETION_REPORT.md`

---

## 📋 AGENT JR - QUICK REFERENCE GUIDE

**Purpose**: Quick lookup for common patterns, existing code references, and execution shortcuts

### **Existing Code References**

#### **Backend Services to Reuse**:
1. **`HybridTrialSearchService`** - `api/services/hybrid_trial_search.py`
   - Already supports `germline_status` and `tumor_context` parameters
   - Returns `List[Dict[str, Any]]` with trial data
   - **USE THIS** - don't rewrite search logic

2. **`ClinicalTrialSearchService`** - `api/services/clinical_trial_search_service.py`
   - AstraDB semantic search
   - **USE THIS** - already integrated with HybridTrialSearchService

#### **Frontend Components to Reference**:
1. **`ResultsDisplay.jsx`** - `src/components/research/ResultsDisplay.jsx`
   - Shows trial card styling patterns
   - LocationCard integration
   - BiomarkerMatchBadge usage

2. **`SporadicProvenanceCard.jsx`** - `src/components/sporadic/SporadicProvenanceCard.jsx`
   - Provenance display patterns

### **Common Patterns**

#### **Backend Pattern - Service Module**:
```python
# api/services/ayesha_trial_matching/eligibility_filters.py
import logging
from typing import Dict, List, Any
from api.services.hybrid_trial_search import HybridTrialSearchService

logger = logging.getLogger(__name__)

class EligibilityFilters:
    """Hard eligibility filters for Ayesha's trials."""
    
    def __init__(self):
        self.search_service = HybridTrialSearchService()
    
    async def apply_hard_filters(self, profile: Dict) -> List[Dict]:
        """Apply all hard filters."""
        # Implementation
        pass
```

#### **Frontend Pattern - Display Component**:
```jsx
// src/components/trials/TrialReasoningSection.jsx
import React from 'react';
import { Box, Typography, Chip } from '@mui/material';

export default function TrialReasoningSection({ reasoning }) {
  if (!reasoning?.why_eligible?.length) return null;
  
  return (
    <Box sx={{ mb: 2 }}>
      <Typography variant="h6">Why Eligible</Typography>
      {reasoning.why_eligible.map((reason, idx) => (
        <Chip key={idx} label={reason} sx={{ mr: 1, mb: 1 }} />
      ))}
    </Box>
  );
}
```

### **Critical Reminders**
- ❌ **NO HARDCODING**: Use `.map()` for arrays, optional chaining (`?.`), fallbacks
- ✅ **USE EXISTING SERVICES**: Don't rewrite search logic
- ✅ **FOLLOW MODULAR PATTERNS**: One file per responsibility, test independently

### **Quick Commands**

#### **Backend Testing**:
```bash
cd oncology-coPilot/oncology-backend-minimal
venv/bin/uvicorn api.main:app --reload
```

#### **Frontend Testing**:
```bash
cd oncology-coPilot/oncology-frontend
npm start
```

#### **AstraDB Seeding**:
```bash
cd oncology-coPilot/oncology-backend-minimal
venv/bin/python scripts/seed_astra_trials.py --disease ovarian --count 200
```

**Full Quick Reference**: See archived `AGENT_JR_QUICK_REFERENCE.md`

---

## 🤝 AGENT JR - PARALLEL MISSION

**Date**: January 8, 2025  
**Status**: ⏸️ **AVAILABLE - PENDING ASSIGNMENT**  
**Mission**: Execute non-conflicting preparatory work while Zo builds core sporadic infrastructure

### **Mission Context**
**Zo is executing:** SPORADIC_CANCER_EXECUTION_PLAN.md (Days 1-7)
- Day 1-2: Backend foundation (TumorContext schema, efficacy gates)
- Day 3: Tumor NGS parsers (Foundation/Tempus)
- Day 4: Clinical trials filtering
- Day 5: Frontend wiring
- Day 6-7: E2E tests & documentation

**Agent Jr can work in parallel on:** Data preparation, documentation, testing infrastructure, and research

### **Parallel Task Options**

#### **OPTION 1: Disease Priors Data Extraction** (HIGH VALUE - 4-6 hours)
- Extract TCGA stats for key cancer types (Ovarian HGS, Breast, Colorectal, Lung, Pancreatic)
- Source from TCGA/cBioPortal/published literature
- Create `disease_priors.json` with mutation rates, HRD prevalence, TMB distributions
- **Deliverable**: `api/resources/disease_priors.json` + `PRIORS_SOURCES.md`

#### **OPTION 2: Foundation Medicine Report Schema Documentation** (MEDIUM VALUE - 3-4 hours)
- Document FM report structure (section headers, field names, variant table structure)
- Create sample canonical JSON
- Map field locations
- **Deliverable**: `foundation_medicine_schema.md` + sample JSON

#### **OPTION 3: Test Data Generation** (HIGH VALUE - 3-4 hours)
- Create 5 realistic patient scenarios (Level 0/1/2 + edge cases)
- Provide expected outputs (PARP penalty Y/N, IO boost Y/N, confidence level)
- **Deliverable**: `.cursor/ayesha/test_scenarios/` (5 JSON files + README)

#### **OPTION 4: Clinical Trials Biomarker Keyword Extraction** (MEDIUM VALUE - 4-5 hours)
- Extract common biomarker patterns (BRCA-required, HRD-eligible, TMB-eligible, MSI-eligible)
- Create exclusion/inclusion regex patterns
- Test against real trial descriptions
- **Deliverable**: `trial_biomarker_patterns.py` + validation results

#### **OPTION 5: Longevity Blog Enhancement** (LOWER PRIORITY - 2-3 hours)
- Add 3 more clinical caselets
- Enhance glossary (HRD, CNAs, Fusions)
- Add FAQ entries
- **Deliverable**: Updated `LONGEVITY_PRECISION_PROTOCOL_BLOG.md`

### **Recommended Priority Order**
**If Agent Jr has 1 full day (8 hours):**
1. ⭐ **OPTION 1: Disease Priors** (4-6 hours) - **HIGHEST VALUE**
2. ⭐ **OPTION 3: Test Data** (3-4 hours) - **HIGH VALUE**

**If Agent Jr has half day (4 hours):**
1. ⭐ **OPTION 1: Disease Priors** (4-6 hours) - **CRITICAL PATH**

### **Critical Rules**
- ❌ **DO NOT TOUCH**: Files Zo is working on (tumor_context.py, orchestrator.py, etc.)
- ✅ **SAFE ZONES**: `api/resources/`, `.cursor/ayesha/`, `tests/`, documentation files
- 📞 **COMMUNICATION**: Commit with prefix `docs(sporadic): ...` or `test(sporadic): ...`

**Full Parallel Mission Details**: See archived `AGENT_JR_PARALLEL_MISSION.md`

---

## ⚔️ AGENT JR - PRE-EXECUTION CLARIFICATIONS

**Date**: January 12, 2025  
**Purpose**: Answer common questions before starting execution

### **Verified - Ready to Proceed**

#### **1. HybridTrialSearchService Structure** ✅
- **Location**: `api/services/hybrid_trial_search.py`
- **Method**: `search_optimized()` with `germline_status` and `tumor_context` parameters
- **Returns**: `List[Dict[str, Any]]` with trial data (nct_id, title, status, phase, etc.)
- **✅ CONFIRMED**: Service already supports `germline_status="negative"` filtering!

#### **2. Router Registration Pattern** ✅
- **Location**: `api/main.py`
- **Pattern**: Simple import and include
```python
from .routers import ayesha_trials as ayesha_trials_router
app.include_router(ayesha_trials_router.router)
```

#### **3. Frontend Routing Pattern** ✅
- **Location**: `src/App.jsx`
- **Pattern**: Standard React Router
```jsx
import AyeshaTrialExplorer from './pages/AyeshaTrialExplorer';
<Route path="/ayesha-trials" element={<AyeshaTrialExplorer />} />
```

#### **4. AstraDB Seeding** ✅
- **Status**: Seeding scripts exist!
- **Scripts**: `seed_astradb_from_sqlite.py`, `seed_trials_standalone.py`, `seed_trials_simple.py`
- **Recommended**: Check if ≥200 ovarian trials exist first, seed if needed

#### **5. Trial Data Structure** ✅
- Each trial Dict has: `nct_id`, `title`, `status`, `phase`, `description`, `eligibility_text`, `interventions`, `locations`, `optimization_score`, `sporadic_filtering_applied`
- **✅ CONFIRMED**: Use `.get()` with defaults for all fields!

### **Zo's Answers to Jr's Questions**

#### **Q1: AstraDB Seeding** ✅ ANSWERED
- Use existing `scripts/seed_trials_simple.py`
- Check count first: `len(results) >= 200 → Skip seeding`

#### **Q2: Trial Contact Info** ✅ ANSWERED
- Leave blank for now
- Display locations clearly
- Add "View on ClinicalTrials.gov" link using `nct_id`

#### **Q3: Location Distance Calculation** ✅ ANSWERED
- Hardcode 'NYC metro' filter for V1
- Filter: `state == "NY" OR state == "NJ" OR state == "CT"`
- Display "📍 NYC Metro" badge

#### **Q4: Conditional NGS Features** ✅ ANSWERED
- Show with "Awaiting NGS" warning
- Grayed WIWFM panel with "🔒 Unlock with NGS" banner
- NGS Fast-Track Checklist displayed prominently

### **Additional Clarifications from Manager**

#### **A1: CA-125 Intelligence Thresholds** ✅
- Burden classes: <100 (MINIMAL), <500 (MODERATE), <1000 (SIGNIFICANT), ≥1000 (EXTENSIVE)
- Forecast: Cycle 3 (≥70% drop), Cycle 6 (≥90% drop), Target (<35 U/mL)
- Resistance signal: On-therapy rise OR <50% drop by cycle 3

#### **A2: Confidence Gates Formula** ✅
- **⚠️ CRITICAL**: This is **trial matching confidence**, NOT drug efficacy confidence
- Formula: `confidence = min(max(gates), 1.0)` with cap at 1.0
- Gates: SOC aligned (0.95), Frontline eligibility ≥0.80 (0.90)

#### **A3: Hard/Soft Criteria Scoring** ✅
- Hard criteria (MUST pass): Stage, Treatment line, Major exclusions
- Soft criteria (% match): ECOG, Age, Distance, Biomarkers, Organ function
- Eligibility gate: 0.90 (≥0.80 soft), 0.85 (≥0.60 soft), 0.75 (<0.60 soft)

#### **A4: Gemini for Eligibility Parsing** ✅
- **⚠️ CRITICAL**: **DO NOT** call Gemini API at runtime!
- Only use cached `structured_criteria` from AstraDB
- If missing, use text-based keyword matching instead

### **Execution Strategy**
1. ✅ **Verify AstraDB**: Check if ≥200 ovarian trials exist
2. ✅ **Read HybridTrialSearchService**: Understand exact return structure
3. ✅ **Read ResultsDisplay.jsx**: Understand frontend patterns
4. ✅ **Create first module**: Start with schemas (no dependencies)

**Full Pre-Execution Clarifications**: See archived `AGENT_JR_PRE_EXECUTION_CLARIFICATIONS.md`

---

## 📖 REFERENCES TO ARCHIVED FILES

All original files have been preserved in `.cursor/ayesha/archive/`:

- **Agent 3 E2E Testing**: `archive/AGENT_3_E2E_TESTING_MISSION.md`
  - Original: Complete E2E testing mission with task breakdown
  
- **Agent Jr Mission 4**: `archive/AGENT_JR_MISSION_4_COMPLETION_REPORT.md`
  - Original: WIWFM integration completion report
  
- **Agent Jr Quick Reference**: `archive/AGENT_JR_QUICK_REFERENCE.md`
  - Original: Quick reference guide for common patterns
  
- **Agent Jr Parallel Mission**: `archive/AGENT_JR_PARALLEL_MISSION.md`
  - Original: Parallel support mission options
  
- **Agent Jr Pre-Execution**: `archive/AGENT_JR_PRE_EXECUTION_CLARIFICATIONS.md`
  - Original: Pre-execution clarifications and answers

---

## ⚔️ DOCTRINE STATUS: ACTIVE

**LAST UPDATED:** January 13, 2025  
**APPLIES TO:** All agent missions, execution guides, and reference materials  
**ENFORCEMENT:** Mandatory across all agent assignments

**This master document represents the complete consolidation of all agent mission knowledge. Every mission, completion report, reference guide, and clarification is preserved and organized for maximum clarity and actionability.**

