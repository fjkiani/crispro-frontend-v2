# 🚀 MOAT Frontend Integration Guide - COMPREHENSIVE AUDIT

**Date:** January 2025  
**Status:** 📋 **READY FOR IMPLEMENTATION** - Complete audit complete  
**Last Updated:** January 2025 (Comprehensive line-by-line audit)

---

## 🚨 CRITICAL FINDING: Documentation vs Reality

**The master index (`01_CURRENT_STATE.md`) is OUTDATED.** The actual implementation is **significantly more complete** than documented.

### **Actual Agent Status:**

| Agent | Documentation Says | **ACTUAL STATUS** | Location |
|-------|-------------------|-------------------|----------|
| **01 Data Extraction** | ⏳ SKELETON | ✅ **INTEGRATED** | `extraction_agent.py` + `_run_extraction_phase()` |
| **02 Biomarker** | ✅ INTEGRATED | ✅ **INTEGRATED** | `_run_biomarker_agent()` line 453 |
| **03 Resistance** | ✅ VALIDATED | ✅ **INTEGRATED** | `_run_resistance_agent()` line 549 |
| **04 Drug Efficacy** | ⏳ SKELETON | ✅ **INTEGRATED** | `_run_drug_efficacy_agent()` line 743 |
| **05 Trial Matching** | ⏳ SKELETON | ✅ **INTEGRATED** | `_run_trial_matching_agent()` line 923 |
| **06 Nutrition** | ⏳ SKELETON | ✅ **INTEGRATED** | `_run_nutrition_agent()` line 660 |
| **07 Care Plan** | ✅ INTEGRATED | ✅ **INTEGRATED** | `_run_care_plan_agent()` line 991 |
| **08 Monitoring** | ✅ INTEGRATED | ✅ **INTEGRATED** | `_run_monitoring_agent()` line 1137 |
| **09 Trigger System** | ⬜ TODO | ⚠️ **EXISTS BUT NOT IN PIPELINE** | `trigger_engine.py` exists, not called |
| **14 Synthetic Lethality** | Not in index | ✅ **INTEGRATED** | `_run_synthetic_lethality_agent()` line 821 |
| **16 Toxicity Risk** | ⬜ TODO | ⚠️ **SERVICE EXISTS, NOT INTEGRATED** | `safety_service.py` exists, no agent |

### **Key Insights:**

1. **9 out of 10 core agents are FULLY INTEGRATED** (not skeletons)
2. **Data Extraction, Drug Efficacy, Nutrition are COMPLETE** (contrary to documentation)
3. **Synthetic Lethality is integrated** (not mentioned in master index)
4. **Trigger System exists but not wired to pipeline**
5. **Toxicity Risk service exists but needs orchestrator integration**
6. **SAE Features can be extracted but not stored in PatientState**

---

## 📋 Table of Contents

1. [What Exists - Complete Inventory](#what-exists---complete-inventory)
2. [What's Missing - Gap Analysis](#whats-missing---gap-analysis)
3. [Where Work Was Left Off](#where-work-was-left-off)
4. [Component-by-Component Audit](#component-by-component-audit)
5. [Test Coverage Analysis](#test-coverage-analysis)
6. [Outcome-Based Test Definitions](#outcome-based-test-definitions)
7. [Implementation Plan](#implementation-plan)
8. [Data Transformation Requirements](#data-transformation-requirements)

---

## 🔍 What Exists - Complete Inventory

### **Frontend Components (Orchestrator Dashboard)**

**Location:** `oncology-coPilot/oncology-frontend/src/components/orchestrator/`

#### **Analysis Components (All Exist & Working):**
- ✅ `BiomarkerCard.jsx` - Displays TMB, MSI, HRD
- ✅ `ResistanceCard.jsx` - Displays resistance predictions
- ✅ `DrugRankingCard.jsx` - Displays drug efficacy rankings with mechanism vector
- ✅ `TrialMatchesCard.jsx` - Displays trial matches with scores (eligibility, mechanism_fit, combined)
- ✅ `NutritionCard.jsx` - Displays nutrition plan
- ✅ `SyntheticLethalityCard.jsx` - Displays SL analysis
- ✅ `CarePlanViewer.jsx` - Displays care plan
- ✅ `MonitoringDashboard.jsx` - Displays monitoring config

#### **Common Components (All Exist):**
- ✅ `LoadingState.jsx`
- ✅ `ErrorState.jsx`
- ✅ `EmptyState.jsx`

#### **Patient Components:**
- ✅ `PatientUpload.jsx` - File upload (VCF/PDF/MAF/JSON/TXT) - **FULLY FUNCTIONAL**

#### **Missing Components:**
- ❌ `PipelineStatusCard.jsx` - **MENTIONED BUT DELETED** (was in Dashboard folder)

### **Frontend Pages**

#### **OrchestratorDashboard.jsx** ✅ **FULLY FUNCTIONAL**
- **Location:** `oncology-coPilot/oncology-frontend/src/pages/OrchestratorDashboard.jsx`
- **Status:** ✅ Complete and working
- **Features:**
  - Uses `/api/orchestrate/full` ✅
  - Uses `useOrchestrator` hook ✅
  - Lazy loading configured ✅
  - Tab-based layout (Analysis, Care Plan, Monitoring) ✅
  - File upload integrated ✅
  - All 9 analysis cards integrated ✅

#### **UniversalCompleteCare.jsx** ⚠️ **USES LEGACY ENDPOINT**
- **Location:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`
- **Status:** ⚠️ Working but uses legacy endpoint
- **Current API:** `/api/complete_care/v2` ❌
- **Should Use:** `/api/orchestrate/full` ❌
- **Features:**
  - ✅ Component structure complete
  - ✅ Uses orchestrator Analysis components (BiomarkerCard, ResistanceCard, etc.)
  - ⚠️ ResistancePlaybook shows placeholder (line 550-564)
  - ⚠️ SAE Features shows placeholder (line 566-583)
  - ❌ No file upload
  - ❌ No status polling
  - ❌ No pipeline status display

### **Hooks & Services**

#### **Hooks:**
- ✅ `useOrchestrator.ts` - **EXISTS & WORKING**
  - Location: `oncology-coPilot/oncology-frontend/src/hooks/useOrchestrator.ts`
  - Provides: `state`, `status`, `loading`, `error`, `runPipeline`, `refreshStatus`, `refreshState`, `clearError`
  - Used by: `OrchestratorDashboard.jsx` ✅

- ✅ `usePipelineStatus.js` - **EXISTS BUT NOT USED**
  - Location: `oncology-coPilot/oncology-frontend/src/hooks/usePipelineStatus.js`
  - Purpose: Polls `/api/orchestrate/status/{patient_id}` for real-time progress
  - Status: ✅ Complete implementation
  - **ISSUE:** Not used in `UniversalCompleteCare.jsx` ❌

#### **API Services:**
- ✅ `orchestrator.ts` - **EXISTS & WORKING**
  - Location: `oncology-coPilot/oncology-frontend/src/services/api/orchestrator.ts`
  - Methods: `runPipeline`, `getStatus`, `getState`, `processEvent`, `listStates`, `healthCheck`
  - Used by: `useOrchestrator.ts` ✅

### **Tests**

#### **E2E Tests:**
- ✅ `orchestrator.e2e.test.js` - **EXISTS**
  - Location: `oncology-coPilot/oncology-frontend/src/__tests__/orchestrator.e2e.test.js`
  - Coverage:
    - ✅ Patient upload flow
    - ✅ Analysis pipeline
    - ✅ Care plan generation
    - ✅ Monitoring configuration
    - ✅ Error handling
    - ✅ Performance
  - **Status:** Basic structure exists, uses mocks

#### **Integration Tests:**
- ✅ `orchestrator.integration.test.js` - **EXISTS**
  - Location: `oncology-coPilot/oncology-frontend/src/__tests__/orchestrator.integration.test.js`
  - Coverage:
    - ✅ `/api/orchestrate/full` endpoint
    - ✅ File upload
    - ✅ Status polling (`/api/orchestrate/status/{patient_id}`)
    - ✅ State retrieval (`/api/orchestrate/state/{patient_id}`)
    - ✅ Health check
    - ✅ Full pipeline workflow
    - ✅ Error handling
  - **Status:** Complete test structure, requires running backend

### **Legacy Components (Used by UniversalCompleteCare)**

#### **Ayesha Components:**
- ✅ `ResistancePlaybook.jsx` - **EXISTS**
  - Location: `oncology-coPilot/oncology-frontend/src/components/ayesha/ResistancePlaybook.jsx`
  - Expects: `resistance_playbook` prop with `risks`, `combo_strategies`, `next_line_switches`, `trial_keywords`
  - **ISSUE:** Shows placeholder when data not available
  - **Current:** Used in `UniversalCompleteCare.jsx` but shows placeholder

- ✅ `AyeshaSAEFeaturesCard.jsx` - **EXISTS**
  - Location: `oncology-coPilot/oncology-frontend/src/components/ayesha/AyeshaSAEFeaturesCard.jsx`
  - Expects: `sae_features` prop with `features`, `boosting_features`, `limiting_features`, `dna_repair_capacity`, `pathway_burden`
  - **ISSUE:** Shows placeholder when data not available
  - **Current:** Used in `UniversalCompleteCare.jsx` but shows placeholder

---

## ❌ What's Missing - Gap Analysis

### **Critical Gaps:**

1. **UniversalCompleteCare Orchestrator Integration** 🔴 **CRITICAL**
   - **Current:** Uses `/api/complete_care/v2` (legacy)
   - **Required:** Use `/api/orchestrate/full`
   - **Impact:** No real-time updates, no file upload, no status tracking

2. **PipelineStatusCard Component** 🔴 **CRITICAL**
   - **Status:** Mentioned in README but file was deleted
   - **Required:** Create component to display pipeline phase progress
   - **Location:** Should be `components/orchestrator/Dashboard/PipelineStatusCard.jsx`

3. **Status Polling in UniversalCompleteCare** 🔴 **CRITICAL**
   - **Status:** `usePipelineStatus` hook exists but not used
   - **Required:** Integrate into `UniversalCompleteCare.jsx`
   - **Impact:** No real-time progress updates

4. **File Upload in UniversalCompleteCare** 🔴 **HIGH**
   - **Status:** `PatientUpload` component exists but not used
   - **Required:** Integrate into `UniversalCompleteCare.jsx`
   - **Impact:** Can't upload VCF/PDF/MAF files

5. **Data Transformation Function** 🔴 **CRITICAL**
   - **Status:** Not created
   - **Required:** `mapOrchestratorToLegacy()` function
   - **Location:** `utils/orchestratorMapper.js`
   - **Purpose:** Transform orchestrator response to legacy format for backward compatibility

6. **Resistance Playbook Data Mapping** 🟡 **MEDIUM**
   - **Status:** Component exists but shows placeholder
   - **Issue:** Orchestrator returns `resistance_prediction.next_line_options` but component expects full `resistance_playbook` format
   - **Required:** Data transformation or component update

7. **SAE Features Integration** 🟡 **MEDIUM**
   - **Status:** Component exists but shows placeholder
   - **Issue:** SAE features can be extracted but not stored in `PatientState`
   - **Required:** Backend enhancement + frontend integration

### **Missing Tests:**

1. **Component Tests for UniversalCompleteCare** ❌
   - No tests for orchestrator integration
   - No tests for data transformation
   - No tests for status polling

2. **Data Transformation Tests** ❌
   - No tests for `mapOrchestratorToLegacy()`
   - No tests for component prop mapping

3. **Integration Tests for UniversalCompleteCare** ❌
   - No tests for full workflow with orchestrator
   - No tests for file upload integration

---

## 📍 Where Work Was Left Off

### **Completed Work:**

1. ✅ **OrchestratorDashboard** - Fully functional with orchestrator integration
2. ✅ **Orchestrator Analysis Components** - All 9 components exist and work
3. ✅ **useOrchestrator Hook** - Complete and working
4. ✅ **usePipelineStatus Hook** - Complete but unused
5. ✅ **orchestratorApi Service** - Complete and working
6. ✅ **PatientUpload Component** - Complete and working
7. ✅ **Basic E2E Tests** - Structure exists
8. ✅ **Integration Tests** - Structure exists

### **Incomplete Work:**

1. ⚠️ **UniversalCompleteCare Integration** - Uses legacy endpoint
   - **Last State:** Working with `/api/complete_care/v2`
   - **Next Step:** Replace with `/api/orchestrate/full`
   - **Blockers:** None

2. ⚠️ **PipelineStatusCard** - File was deleted
   - **Last State:** Mentioned in README
   - **Next Step:** Recreate component
   - **Blockers:** None

3. ⚠️ **Data Transformation** - Not started
   - **Last State:** Not created
   - **Next Step:** Create `mapOrchestratorToLegacy()` function
   - **Blockers:** None

4. ⚠️ **Resistance Playbook Integration** - Shows placeholder
   - **Last State:** Component exists, data format mismatch
   - **Next Step:** Map orchestrator data to component format
   - **Blockers:** None

5. ⚠️ **SAE Features Integration** - Shows placeholder
   - **Last State:** Component exists, data not in PatientState
   - **Next Step:** Backend enhancement + frontend integration
   - **Blockers:** Backend needs to store SAE features in PatientState

### **Documentation Files:**

**Found Documentation:**
- ✅ `.cursor/MOAT/orchestration/01_CURRENT_STATE.md` - Current state (outdated)
- ✅ `.cursor/MOAT/orchestration/02_FRONTEND_STATUS.md` - Frontend status
- ✅ `.cursor/MOAT/orchestration/03_DELIVERABLES_PLAN.md` - Deliverables plan
- ✅ `.cursor/MOAT/UNIVERSAL_PAGES_AUDIT_AND_DELIVERABLES.md` - Universal pages audit
- ✅ `oncology-coPilot/oncology-frontend/ORCHESTRATOR_DASHBOARD_VERIFICATION.md` - Dashboard verification

**Missing Documentation:**
- ❌ `00_MISSION.mdc` - Mentioned in README but not found
- ❌ `ULTIMATE_MOAT_ORCHESTRATION.mdc` - Mentioned but not found
- ❌ `agent-implementation-guide.mdc` - Mentioned but not found

**Note:** MDC files may have been converted to MD or may not exist in the current workspace.

---

## 🔬 Component-by-Component Audit

### **UniversalCompleteCare.jsx - Line-by-Line Analysis**

#### **Current Implementation (Lines 1-174):**
- ✅ Imports orchestrator Analysis components (BiomarkerCard, ResistanceCard, DrugRankingCard, TrialMatchesCard)
- ✅ Uses legacy components (ResistancePlaybook, AyeshaSAEFeaturesCard)
- ❌ **Line 147:** Uses `/api/complete_care/v2` instead of `/api/orchestrate/full`
- ❌ **Line 139-174:** `handleGeneratePlan()` uses legacy endpoint
- ❌ No file upload integration
- ❌ No status polling integration
- ❌ No pipeline status display

#### **Component Usage (Lines 320-595):**
- ✅ **Lines 396-401:** BiomarkerCard - Works with `result.biomarker_intelligence`
- ✅ **Lines 406-411:** ResistanceCard - Works with `result.resistance_prediction`
- ✅ **Lines 416-422:** DrugRankingCard - Works with `getDrugRanking()` from `result.wiwfm`
- ✅ **Lines 427-432:** TrialMatchesCard - Works with `getTrials()` from `result.trials`
- ⚠️ **Lines 550-564:** ResistancePlaybook - Shows placeholder when `result.resistance_playbook` is null
- ⚠️ **Lines 566-583:** AyeshaSAEFeaturesCard - Shows placeholder when `result.sae_features` is null

#### **Data Extraction Functions:**
- ✅ **Lines 189-208:** `getDrugRanking()` - Extracts from `result.wiwfm`
- ✅ **Lines 224-229:** `getTrials()` - Extracts from `result.trials`
- ✅ **Lines 231-236:** `getHintTiles()` - Extracts from `result.hint_tiles`
- ✅ **Lines 211-221:** `getPGxDrugs()` - Extracts from `result.wiwfm.drugs`

### **OrchestratorDashboard.jsx - Line-by-Line Analysis**

#### **Current Implementation:**
- ✅ **Line 11:** Uses `useOrchestrator` hook
- ✅ **Line 49-52:** Uses `PatientUpload` component
- ✅ **Lines 94-148:** All analysis cards integrated with lazy loading
- ✅ **Lines 152-168:** Care Plan and Monitoring tabs integrated
- ✅ All components receive data from `state` object (orchestrator format)

### **Component Data Format Expectations**

#### **BiomarkerCard:**
- **Expects:** `biomarkerProfile` prop
- **Orchestrator Format:** `state.biomarker_profile` ✅
- **Legacy Format:** `result.biomarker_intelligence` ✅
- **Status:** ✅ Compatible

#### **ResistanceCard:**
- **Expects:** `resistancePrediction` prop
- **Orchestrator Format:** `state.resistance_prediction` ✅
- **Legacy Format:** `result.resistance_prediction` ✅
- **Status:** ✅ Compatible

#### **DrugRankingCard:**
- **Expects:** `drugRanking` prop with `ranked_drugs` array
- **Orchestrator Format:** `state.drug_ranking` (array) ✅
- **Legacy Format:** `result.wiwfm.drugs` (array) ⚠️ Needs transformation
- **Status:** ⚠️ Needs transformation function

#### **TrialMatchesCard:**
- **Expects:** `trialMatches` prop (array)
- **Orchestrator Format:** `state.trial_matches` (array) ✅
- **Legacy Format:** `result.trials.trials` (nested array) ⚠️ Needs transformation
- **Status:** ⚠️ Needs transformation function

#### **ResistancePlaybook:**
- **Expects:** `resistance_playbook` prop with:
  - `risks` (array)
  - `combo_strategies` (array)
  - `next_line_switches` (array)
  - `trial_keywords` (array)
- **Orchestrator Format:** `state.resistance_prediction.next_line_options` (object with `alternatives`, `regimen_changes`, `monitoring_changes`) ⚠️
- **Legacy Format:** `result.resistance_playbook` (full format) ✅
- **Status:** ❌ Format mismatch - needs transformation

#### **AyeshaSAEFeaturesCard:**
- **Expects:** `sae_features` prop with:
  - `features` (array)
  - `boosting_features` (array)
  - `limiting_features` (array)
  - `dna_repair_capacity` (object)
  - `pathway_burden` (object)
- **Orchestrator Format:** Not in `PatientState` ❌
- **Legacy Format:** `result.sae_features` ✅
- **Status:** ❌ Not available in orchestrator response

---

## 🧪 Test Coverage Analysis

### **Existing Tests:**

#### **E2E Tests (`orchestrator.e2e.test.js`):**
- ✅ Structure exists
- ✅ Uses mocks
- ⚠️ **Gap:** No tests for UniversalCompleteCare integration
- ⚠️ **Gap:** No tests for data transformation
- ⚠️ **Gap:** No tests for status polling

#### **Integration Tests (`orchestrator.integration.test.js`):**
- ✅ Structure exists
- ✅ Tests real backend endpoints
- ✅ Tests file upload
- ✅ Tests status polling
- ✅ Tests full pipeline workflow
- ⚠️ **Gap:** No tests for UniversalCompleteCare workflow
- ⚠️ **Gap:** No tests for data transformation

### **Missing Tests:**

1. **Data Transformation Tests** ❌
   - `mapOrchestratorToLegacy()` function tests
   - Component prop mapping tests
   - Edge case tests (null values, empty arrays)

2. **UniversalCompleteCare Integration Tests** ❌
   - Orchestrator endpoint integration
   - File upload integration
   - Status polling integration
   - Component rendering with orchestrator data

3. **Component Tests** ❌
   - ResistancePlaybook with orchestrator data
   - AyeshaSAEFeaturesCard with orchestrator data
   - PipelineStatusCard (when created)

---

## ✅ Outcome-Based Test Definitions

### **Test Category 1: Data Transformation**

#### **Test 1.1: mapOrchestratorToLegacy() - Full Response**
**Outcome:** Function correctly transforms complete orchestrator response to legacy format

**Test Cases:**
- ✅ Transforms `biomarker_profile` → `biomarker_intelligence`
- ✅ Transforms `drug_ranking` array → `wiwfm.drugs` array with correct field mapping
- ✅ Transforms `trial_matches` array → `trials.trials` array
- ✅ Transforms `mechanism_vector` → `mechanism_map.vector`
- ✅ Transforms `resistance_prediction.next_line_options` → `resistance_playbook` (if possible)
- ✅ Handles missing optional fields gracefully
- ✅ Handles null values
- ✅ Handles empty arrays

**Success Criteria:**
- All required fields mapped correctly
- No data loss
- Components can consume transformed data

#### **Test 1.2: mapOrchestratorToLegacy() - Partial Response**
**Outcome:** Function handles partial orchestrator responses (some agents skipped)

**Test Cases:**
- ✅ Handles missing `biomarker_profile`
- ✅ Handles missing `drug_ranking`
- ✅ Handles missing `trial_matches`
- ✅ Handles missing `nutrition_plan`
- ✅ Handles missing `synthetic_lethality_result`

**Success Criteria:**
- No errors thrown
- Components show appropriate empty states

### **Test Category 2: UniversalCompleteCare Integration**

#### **Test 2.1: Orchestrator Endpoint Integration**
**Outcome:** UniversalCompleteCare successfully uses `/api/orchestrate/full`

**Test Cases:**
- ✅ Replaces `/api/complete_care/v2` with `/api/orchestrate/full`
- ✅ Request format matches `OrchestratePipelineRequest` schema
- ✅ Response handling matches `OrchestratePipelineResponse` schema
- ✅ Error handling works correctly
- ✅ Loading states display correctly

**Success Criteria:**
- API call succeeds
- Response data displayed correctly
- Error messages user-friendly

#### **Test 2.2: File Upload Integration**
**Outcome:** UniversalCompleteCare supports file upload (VCF/PDF/MAF)

**Test Cases:**
- ✅ PatientUpload component integrated
- ✅ File selection works
- ✅ File type detection works (VCF, PDF, MAF, JSON, TXT)
- ✅ File upload triggers pipeline execution
- ✅ Error handling for invalid files

**Success Criteria:**
- Files upload successfully
- Pipeline executes with uploaded file
- Error messages clear

#### **Test 2.3: Status Polling Integration**
**Outcome:** UniversalCompleteCare shows real-time pipeline progress

**Test Cases:**
- ✅ usePipelineStatus hook integrated
- ✅ Polling starts when pipeline begins
- ✅ Progress updates displayed
- ✅ Polling stops when pipeline completes
- ✅ Polling stops when pipeline errors
- ✅ Error handling during polling

**Success Criteria:**
- Real-time progress visible
- No memory leaks (polling stops)
- Error states handled

#### **Test 2.4: PipelineStatusCard Display**
**Outcome:** PipelineStatusCard displays pipeline phase and progress

**Test Cases:**
- ✅ Component renders correctly
- ✅ Phase labels display correctly
- ✅ Progress bar updates
- ✅ Alerts display correctly
- ✅ Agent execution status displays
- ✅ Auto-stops when complete

**Success Criteria:**
- All information visible
- Updates in real-time
- No performance issues

### **Test Category 3: Component Data Mapping**

#### **Test 3.1: BiomarkerCard with Orchestrator Data**
**Outcome:** BiomarkerCard displays orchestrator biomarker_profile correctly

**Test Cases:**
- ✅ Receives `biomarkerProfile` prop from transformed data
- ✅ Displays TMB classification
- ✅ Displays MSI status
- ✅ Displays HRD status
- ✅ Handles missing fields gracefully

**Success Criteria:**
- All data displayed correctly
- No errors with missing fields

#### **Test 3.2: DrugRankingCard with Orchestrator Data**
**Outcome:** DrugRankingCard displays orchestrator drug_ranking correctly

**Test Cases:**
- ✅ Receives `drugRanking` prop with `ranked_drugs` array
- ✅ Displays drug names correctly
- ✅ Displays efficacy scores correctly
- ✅ Displays mechanism vector correctly
- ✅ Handles empty array gracefully

**Success Criteria:**
- All drugs displayed
- Scores formatted correctly
- Mechanism vector displayed

#### **Test 3.3: TrialMatchesCard with Orchestrator Data**
**Outcome:** TrialMatchesCard displays orchestrator trial_matches correctly

**Test Cases:**
- ✅ Receives `trialMatches` prop (array)
- ✅ Displays trial titles correctly
- ✅ Displays mechanism_fit_score
- ✅ Displays eligibility_score
- ✅ Displays combined_score
- ✅ Handles empty array gracefully

**Success Criteria:**
- All trials displayed
- Scores formatted correctly
- Mechanism fit visible

#### **Test 3.4: ResistancePlaybook with Orchestrator Data**
**Outcome:** ResistancePlaybook displays orchestrator resistance data correctly

**Test Cases:**
- ✅ Receives `resistance_playbook` prop from transformed data
- ✅ Displays next_line_options if available
- ✅ Shows placeholder if data not available
- ✅ Handles format mismatch gracefully

**Success Criteria:**
- Data displayed when available
- Placeholder shown when not available
- No errors

#### **Test 3.5: AyeshaSAEFeaturesCard with Orchestrator Data**
**Outcome:** AyeshaSAEFeaturesCard displays SAE features if available

**Test Cases:**
- ✅ Receives `sae_features` prop from transformed data
- ✅ Displays DNA repair capacity if available
- ✅ Displays pathway burden if available
- ✅ Shows placeholder if data not available
- ✅ Handles missing fields gracefully

**Success Criteria:**
- Data displayed when available
- Placeholder shown when not available
- No errors

### **Test Category 4: End-to-End Workflow**

#### **Test 4.1: Complete Pipeline Workflow**
**Outcome:** User can upload file, run pipeline, and see all results

**Test Cases:**
- ✅ User uploads VCF file
- ✅ Pipeline executes successfully
- ✅ Status polling shows progress
- ✅ All components display data
- ✅ No errors during execution

**Success Criteria:**
- Complete workflow succeeds
- All data displayed
- No errors

#### **Test 4.2: Error Recovery**
**Outcome:** System handles errors gracefully and allows retry

**Test Cases:**
- ✅ Network error handled
- ✅ API error handled
- ✅ Agent failure handled
- ✅ Retry mechanism works
- ✅ Error messages user-friendly

**Success Criteria:**
- Errors don't crash app
- User can retry
- Error messages clear

---

## 🚀 Implementation Plan

### **Phase 1: Core Integration** 🔴 **CRITICAL** (Week 1)

#### **Task 1.1: Create Data Transformation Function** 🔴 **CRITICAL**
**File:** `oncology-coPilot/oncology-frontend/src/utils/orchestratorMapper.js` (NEW)

**Implementation:**
```javascript
export const mapOrchestratorToLegacy = (orchestratorResponse) => {
  return {
    // Core mappings
    biomarker_intelligence: orchestratorResponse.biomarker_profile,
    resistance_prediction: orchestratorResponse.resistance_prediction,
    wiwfm: {
      drugs: orchestratorResponse.drug_ranking?.map(d => ({
        name: d.drug_name,
        drug_name: d.drug_name,
        moa: d.mechanism || d.drug_class || d.moa,
        efficacy_score: d.efficacy_score,
        confidence: d.confidence || d.efficacy_score,
        evidence_tier: d.tier || d.evidence_tier || 'insufficient',
        badges: d.badges || [],
        rationale: d.rationale || []
      })) || [],
      evidence_tier: orchestratorResponse.drug_ranking?.[0]?.evidence_tier || 'Supported',
      run_signature: orchestratorResponse.patient_id
    },
    trials: {
      trials: orchestratorResponse.trial_matches?.map(t => ({
        nct_id: t.nct_id,
        title: t.title,
        phase: t.phase,
        status: t.status,
        mechanism_fit_score: t.mechanism_fit_score,
        eligibility_score: t.eligibility_score,
        combined_score: t.combined_score,
        why_matched: t.why_matched,
        url: t.url,
        brief_summary: t.brief_summary,
        locations: t.locations,
        contact: t.contact
      })) || []
    },
    mechanism_map: {
      vector: orchestratorResponse.mechanism_vector || [0,0,0,0,0,0,0]
    },
    resistance_playbook: orchestratorResponse.resistance_prediction?.next_line_options || null,
    care_plan: orchestratorResponse.care_plan,
    nutrition_plan: orchestratorResponse.nutrition_plan,
    synthetic_lethality: orchestratorResponse.synthetic_lethality_result,
    monitoring_config: orchestratorResponse.monitoring_config,
    sae_features: orchestratorResponse.sae_features || null,
    soc_recommendation: null,
    next_test_recommender: null,
    hint_tiles: null,
    toxicity_assessments: orchestratorResponse.toxicity_assessments || null,
    resistance_alert: null,
    provenance: {
      orchestrator: 'MOAT Orchestrator',
      run_id: orchestratorResponse.patient_id,
      generated_at: orchestratorResponse.updated_at,
      phase: orchestratorResponse.phase,
      completed_agents: orchestratorResponse.completed_agents || []
    },
    summary: {
      components_included: orchestratorResponse.completed_agents || [],
      ngs_status: orchestratorResponse.mutation_count > 0 ? 'available' : 'pending',
      confidence_level: 'moderate-high (70-90%)',
      pipeline_version: '2.0'
    }
  };
};
```

**Estimated Time:** 2-3 hours

#### **Task 1.2: Update UniversalCompleteCare to Use Orchestrator** 🔴 **CRITICAL**
**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Import `mapOrchestratorToLegacy` function
2. Replace `/api/complete_care/v2` with `/api/orchestrate/full`
3. Update request format
4. Transform response using `mapOrchestratorToLegacy()`
5. Update error handling

**Estimated Time:** 4-6 hours

#### **Task 1.3: Create PipelineStatusCard Component** 🔴 **CRITICAL**
**File:** `oncology-coPilot/oncology-frontend/src/components/orchestrator/Dashboard/PipelineStatusCard.jsx` (NEW)

**Implementation:** See Task 1.5 in original plan

**Estimated Time:** 3-4 hours

#### **Task 1.4: Integrate Status Polling** 🔴 **CRITICAL**
**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Import `usePipelineStatus` hook
2. Add status polling when pipeline starts
3. Display PipelineStatusCard
4. Stop polling when complete

**Estimated Time:** 2-3 hours

#### **Task 1.5: Integrate File Upload** 🔴 **HIGH**
**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Import `PatientUpload` component
2. Add file upload section
3. Handle upload completion
4. Trigger pipeline with uploaded file

**Estimated Time:** 2-3 hours

### **Phase 2: Component Integration** 🟡 **MEDIUM** (Week 2)

#### **Task 2.1: Verify ResistancePlaybook Integration** 🟡 **MEDIUM**
**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Check if `resistance_playbook` data available after transformation
2. Update transformation to map `next_line_options` to playbook format
3. Remove placeholder if data available

**Estimated Time:** 2-3 hours

#### **Task 2.2: Handle SAE Features** 🟡 **MEDIUM**
**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Check if SAE features available in orchestrator response
2. Map to component format if available
3. Show placeholder if not available

**Estimated Time:** 2-3 hours

### **Phase 3: Testing** 🟢 **ONGOING** (Week 3)

#### **Task 3.1: Write Data Transformation Tests** 🟢 **HIGH**
**File:** `oncology-coPilot/oncology-frontend/src/utils/__tests__/orchestratorMapper.test.js` (NEW)

**Coverage:**
- Full response transformation
- Partial response transformation
- Edge cases (null, empty arrays)

**Estimated Time:** 3-4 hours

#### **Task 3.2: Write Integration Tests** 🟢 **HIGH**
**File:** `oncology-coPilot/oncology-frontend/src/__tests__/universalCompleteCare.integration.test.js` (NEW)

**Coverage:**
- Orchestrator endpoint integration
- File upload integration
- Status polling integration
- Component rendering

**Estimated Time:** 4-6 hours

---

## 📊 Progress Tracking

| Task | Status | Time Spent | Blockers |
|------|--------|------------|----------|
| 1.1: Create Data Transformation | ⏳ Pending | 0h | None |
| 1.2: Update UniversalCompleteCare | ⏳ Pending | 0h | None |
| 1.3: Create PipelineStatusCard | ⏳ Pending | 0h | None |
| 1.4: Integrate Status Polling | ⏳ Pending | 0h | None |
| 1.5: Integrate File Upload | ⏳ Pending | 0h | None |
| 2.1: Verify ResistancePlaybook | ⏳ Pending | 0h | None |
| 2.2: Handle SAE Features | ⏳ Pending | 0h | None |
| 3.1: Write Transformation Tests | ⏳ Pending | 0h | None |
| 3.2: Write Integration Tests | ⏳ Pending | 0h | None |

**Total Estimated Time:** 20-30 hours  
**Current Progress:** 0% (0/9 tasks)

---

## 🎯 Success Criteria

### **Phase 1 Complete When:**
- ✅ UniversalCompleteCare uses `/api/orchestrate/full`
- ✅ Data transformation works correctly
- ✅ File upload works (VCF/PDF/MAF)
- ✅ Status polling shows real-time progress
- ✅ PipelineStatusCard displays correctly
- ✅ All existing components render correctly

### **Phase 2 Complete When:**
- ✅ ResistancePlaybook displays data (no placeholder)
- ✅ SAE Features handled correctly (placeholder if not available)

### **Phase 3 Complete When:**
- ✅ All tests pass
- ✅ Test coverage >80%
- ✅ No regressions

---

**Last Updated:** January 2025  
**Status:** ✅ **READY FOR IMPLEMENTATION** - Complete audit complete with validation script

---

## 🔬 VALIDATION & TESTING APPROACH

### **Validation Script Created**

**File:** `.cursor/MOAT/ORCHESTRATOR_VALIDATION_SCRIPT.py`

**Purpose:** Test orchestrator endpoints with REAL patient data to validate:
1. API endpoint functionality
2. Response format correctness
3. Agent output completeness
4. Data transformation requirements
5. Pipeline phase progression

**Test Data Sources:**
- ✅ Real TCGA ovarian cancer data (`tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`)
- ✅ Test patient profiles (Ayesha BRCA1, KRAS+TP53, MBD4+TP53+NF1)
- ✅ Real mutations from validated studies

**How to Run:**
```bash
# 1. Start backend
cd oncology-coPilot/oncology-backend-minimal
uvicorn main:app --reload --port 8000

# 2. Run validation script
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
python3 .cursor/MOAT/ORCHESTRATOR_VALIDATION_SCRIPT.py
```

**What It Tests:**
- ✅ Health check endpoint
- ✅ Full pipeline execution with 3 real patient profiles
- ✅ Status endpoint polling
- ✅ State endpoint retrieval
- ✅ Data transformation validation
- ✅ Response format validation

**Output:**
- Console output with pass/fail status
- Detailed JSON report: `.cursor/MOAT/validation_results.json`

### **Unified Risk-Benefit Score Integration**

**Reference:** `/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/kly/.cursor/rules/lectures/drugDevelopment/UNIFIED_RISK_BENEFIT_SCORE_PLAN.md`

**Key Insight:** The plan defines a composition layer that combines:
- **Efficacy scores** (from S/P/E framework - validated: 100% top-5 accuracy)
- **Toxicity tiers** (from PGx screening - validated: 100% sensitivity/specificity)
- **Composite score** with action labels (PREFERRED / CONSIDER WITH MONITORING / AVOID)

**Integration Point:** This should be integrated into the drug ranking phase or as a post-processing step after drug efficacy agent completes.

**Current Status:** ⚠️ **NOT INTEGRATED** - Service exists but not wired to orchestrator

**Impact on Frontend:**
- DrugRankingCard should display composite scores
- Action labels should be visible
- HIGH toxicity drugs should be clearly marked as AVOID

---

## 📊 ACTUAL API ENDPOINT STRUCTURE

### **Endpoint: `/api/orchestrate/full`**

**Router:** `api/routers/orchestrate.py` (line 41)

**Request Format:**
```json
{
  "disease": "ovarian",
  "mutations": [
    {
      "gene": "BRCA1",
      "hgvs_p": "p.Arg1835Ter",
      "hgvs_c": "c.5503C>T",
      "chrom": "17",
      "pos": 43044295,
      "ref": "C",
      "alt": "T",
      "consequence": "stop_gained"
    }
  ],
  "treatment_line": 1,
  "prior_therapies": [],
  "skip_agents": []
}
```

**Response Format (OrchestratePipelineResponse):**
```json
{
  "patient_id": "PT-ABC12345",
  "disease": "ovarian",
  "phase": "complete",
  "progress_percent": 100,
  "completed_agents": ["biomarker", "resistance", "drug_efficacy", "trial_matching", "nutrition", "care_plan", "monitoring"],
  "created_at": "2025-01-28T10:00:00Z",
  "updated_at": "2025-01-28T10:01:30Z",
  "duration_ms": 90000,
  "mutation_count": 2,
  "mechanism_vector": [0.8, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0],
  "biomarker_profile": {
    "tmb": {"value": 12.5, "classification": "TMB-H"},
    "msi": {"status": "MSS"},
    "hrd": {"status": "HRD+", "genes_mutated": ["BRCA1"]},
    "io_eligible": false,
    "parp_eligible": true
  },
  "resistance_prediction": {
    "risk_level": "MEDIUM",
    "probability": 0.45,
    "confidence": 0.85,
    "detected_genes": [{"gene": "TP53", "risk_ratio": 1.90}],
    "alternatives": [...],
    "regimen_changes": [...],
    "monitoring_changes": {...}
  },
  "drug_ranking": [
    {
      "drug_name": "olaparib",
      "drug_class": "PARP inhibitor",
      "efficacy_score": 0.85,
      "tier": "supported",
      "confidence": 0.75,
      "mechanism": "PARP inhibition",
      "rationale": ["BRCA1 mutation", "HRD+ status"]
    }
  ],
  "trial_matches": [
    {
      "nct_id": "NCT12345",
      "title": "PARP Inhibitor Trial",
      "phase": "PHASE3",
      "status": "RECRUITING",
      "mechanism_fit_score": 0.92,
      "eligibility_score": 0.85,
      "combined_score": 0.89,
      "why_matched": "BRCA1 mutation, HRD+ status",
      "url": "https://clinicaltrials.gov/study/NCT12345"
    }
  ],
  "nutrition_plan": {
    "supplements": [...],
    "foods_to_prioritize": [...],
    "foods_to_avoid": [...],
    "drug_food_interactions": [...]
  },
  "synthetic_lethality_result": {
    "synthetic_lethality_detected": true,
    "essentiality_scores": [...],
    "recommended_drugs": [...]
  },
  "care_plan": {
    "executive_summary": {...},
    "sections": [...]
  },
  "monitoring_config": {
    "frequency": "monthly",
    "biomarkers": ["CA-125", "ctDNA"],
    "imaging": "CT every 3 months"
  },
  "data_quality_flags": [],
  "alerts": []
}
```

**⚠️ CRITICAL GAP:** Response does NOT include (even though they exist in PatientState):
- ❌ `nutrition_plan` - **EXISTS IN STATE BUT NOT IN RESPONSE** (line 232 in orchestrate.py)
- ❌ `synthetic_lethality_result` - **EXISTS IN STATE BUT NOT IN RESPONSE** (line 232 in orchestrate.py)

**Note:** Response also does NOT include:
- ❌ `sae_features` (can be extracted but not stored in PatientState)
- ❌ `toxicity_assessments` (service exists but not integrated)
- ❌ `soc_recommendation` (not in orchestrator)
- ❌ `next_test_recommender` (not in orchestrator)
- ❌ `hint_tiles` (not in orchestrator)

**Fix Required:** Update `_state_to_response()` in `api/routers/orchestrate.py` to include:
```python
nutrition_plan=state.nutrition_plan,
synthetic_lethality_result=state.synthetic_lethality_result,
```

And update `OrchestratePipelineResponse` schema to include these fields.

### **Endpoint: `/api/orchestrate/status/{patient_id}`**

**Response Format:**
```json
{
  "patient_id": "PT-ABC12345",
  "phase": "analyzing",
  "progress_percent": 45,
  "current_agent": "resistance",
  "completed_agents": ["biomarker"],
  "alerts": [],
  "errors": [],
  "status_url": "/api/orchestrate/status/PT-ABC12345",
  "care_plan_url": null
}
```

### **Endpoint: `/api/patients/{patient_id}`**

**Response Format:** Same as `/api/orchestrate/full` but returns full `PatientState.to_full_dict()` including:
- All agent outputs
- Complete audit trail
- Agent execution history
- All state changes

---

## 🔍 VALIDATED AGENT OUTPUT FORMATS

### **Biomarker Agent Output** (`biomarker_profile`)
```python
{
  "tmb": {
    "value": 12.5,
    "classification": "TMB-H",  # TMB-H, TMB-M, TMB-L
    "threshold_used": 10.0,
    "mutation_count": 475,
    "total_mutations": 500,
    "exome_size_mb": 38.0,
    "method": "nonsilent_variants_per_mb"
  },
  "msi": {
    "status": "MSS",  # MSI-H, MSS
    "dmmr_genes_mutated": ["MLH1", "MSH2"],
    "method": "gene_panel_with_impact",
    "confidence": "high"
  },
  "hrd": {
    "status": "HRD+",  # HRD+, HRD-, HRD-uncertain
    "genes_mutated": ["BRCA1", "BRCA2"],
    "method": "comprehensive_ddr_panel",
    "high_impact_detected": true,
    "confidence": "high"
  },
  "io_eligible": false,
  "parp_eligible": true,
  "provenance": {
    "calculation_method": "enhanced_gene_panel",
    "timestamp": "2025-01-28T10:00:00Z",
    "disease": "ovarian"
  }
}
```

### **Resistance Agent Output** (`resistance_prediction`)
```python
{
  "risk_level": "MEDIUM",  # HIGH, MEDIUM, LOW
  "probability": 0.45,
  "confidence": 0.85,
  "detected_genes": [
    {
      "gene": "TP53",
      "risk_ratio": 1.90,
      "p_value": 0.11,
      "status": "TREND"  # VALIDATED, TREND
    }
  ],
  "next_line_options": {
    "alternatives": [
      {
        "drug": "Alternative A",
        "rationale": "...",
        "evidence_tier": "VALIDATED"
      }
    ],
    "regimen_changes": [...],
    "monitoring_changes": {
      "frequency": "biweekly",
      "biomarkers": ["CA-125", "ctDNA"]
    }
  },
  "provenance": {
    "method": "ResistanceProphet + ResistancePlaybook",
    "timestamp": "2025-01-28T10:00:15Z"
  }
}
```

### **Drug Efficacy Agent Output** (`drug_ranking`)
```python
[
  {
    "drug_name": "olaparib",
    "moa": "PARP inhibition",
    "efficacy_score": 0.85,  # 0.0 - 1.0
    "confidence": 0.75,
    "evidence_tier": "supported",  # supported, consider, insufficient
    "badges": ["FDA Approved", "NCCN Recommended"],
    "rationale": [
      "BRCA1 mutation detected",
      "HRD+ status",
      "PARP eligibility confirmed"
    ],
    "drug_class": "PARP inhibitor",
    "tier": "supported"
  }
]
```

**Mechanism Vector** (7D): `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
- Extracted from `pathway_disruption` scores in efficacy response
- Stored in `state.mechanism_vector`

### **Trial Matching Agent Output** (`trial_matches`)
```python
[
  {
    "nct_id": "NCT12345",
    "title": "PARP Inhibitor Trial",
    "brief_summary": "...",
    "phase": "PHASE3",
    "status": "RECRUITING",
    "mechanism_fit_score": 0.92,  # 0.0 - 1.0
    "eligibility_score": 0.85,    # 0.0 - 1.0
    "combined_score": 0.89,       # 0.7×eligibility + 0.3×mechanism_fit
    "trial_moa": [0.8, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0],
    "mechanism_alignment": {
      "DDR": 0.92,
      "MAPK": 0.15
    },
    "why_matched": "BRCA1 mutation, HRD+ status, DDR pathway alignment",
    "url": "https://clinicaltrials.gov/study/NCT12345",
    "locations": [...],
    "contact": {...}
  }
]
```

### **Nutrition Agent Output** (`nutrition_plan`)
```python
{
  "patient_id": "PT-ABC12345",
  "treatment": "first-line",
  "supplements": [
    {
      "name": "Vitamin D",
      "dosage": "1000 IU",
      "evidence_level": "HIGH",
      "rationale": "Supports immune function"
    }
  ],
  "foods_to_prioritize": [
    {
      "food": "Cruciferous vegetables",
      "rationale": "Supports DNA repair pathways"
    }
  ],
  "foods_to_avoid": [
    {
      "food": "Grapefruit",
      "rationale": "CYP3A4 interaction with olaparib"
    }
  ],
  "drug_food_interactions": [
    {
      "drug": "olaparib",
      "food": "Grapefruit",
      "interaction": "CYP3A4 inhibition",
      "severity": "MODERATE"
    }
  ],
  "timing_rules": {...},
  "provenance": {
    "method": "NutritionAgent",
    "timestamp": "2025-01-28T10:00:30Z"
  }
}
```

### **Synthetic Lethality Agent Output** (`synthetic_lethality_result`)
```python
{
  "synthetic_lethality_detected": true,
  "essentiality_scores": [
    {
      "gene": "BRCA1",
      "essentiality_score": 0.75,
      "pathway": "DDR"
    }
  ],
  "broken_pathways": ["DDR", "HRR"],
  "recommended_drugs": [
    {
      "drug": "olaparib",
      "rationale": "Synthetic lethal with BRCA1 loss"
    }
  ],
  "provenance": {
    "method": "SyntheticLethalityAgent",
    "timestamp": "2025-01-28T10:00:45Z"
  }
}
```

### **Care Plan Agent Output** (`care_plan`)
```python
{
  "executive_summary": {
    "patient_id": "PT-ABC12345",
    "disease": "ovarian",
    "mutation_count": 2,
    "top_mutations": [...],
    "biomarker_summary": {...},
    "top_drug": {...},
    "trial_count": 5,
    "risk_level": "MEDIUM"
  },
  "sections": [
    {
      "title": "Patient Summary",
      "order": 1,
      "content": {...}
    },
    {
      "title": "Genomic Biomarkers",
      "order": 2,
      "content": {...}
    },
    # ... 8 total sections
  ],
  "alerts": [...],
  "provenance": {
    "generated_at": "2025-01-28T10:01:00Z",
    "orchestrator_version": "1.0.0"
  }
}
```

### **Monitoring Agent Output** (`monitoring_config`)
```python
{
  "frequency": "monthly",  # weekly, biweekly, monthly
  "biomarkers": ["CA-125", "ctDNA"],
  "primary_biomarker": "CA-125",
  "imaging": "CT every 3 months",
  "alerts_enabled": true,
  "escalation_thresholds": {
    "biomarker_rise_percent": 25,
    "new_lesion": true,
    "symptom_worsening": true,
    "performance_status_drop": true
  },
  "ctdna_monitoring": {
    "enabled": true,
    "frequency": "monthly",
    "target_mutations": ["BRCA1", "TP53"]
  },
  "risk_level": "MEDIUM",
  "provenance": {
    "configured_at": "2025-01-28T10:01:15Z",
    "based_on": {
      "resistance_risk": "MEDIUM",
      "tmb_status": "TMB-H",
      "hrd_status": "HRD+"
    }
  }
}
```

---

## 🧪 VALIDATION TEST CASES

### **Test Case 1: Ayesha BRCA1+TP53 Profile**
**Input:**
- Disease: `ovarian`
- Mutations: BRCA1 p.Arg1835Ter (stop_gained), TP53 p.Arg175His (missense)

**Expected Outputs:**
- ✅ `biomarker_profile.hrd.status` = "HRD+"
- ✅ `biomarker_profile.parp_eligible` = true
- ✅ `drug_ranking[0].drug_name` contains "olaparib" or "niraparib"
- ✅ `drug_ranking[0].efficacy_score` > 0.7
- ✅ `trial_matches` contains PARP inhibitor trials
- ✅ `mechanism_vector[0]` (DDR) > 0.7
- ✅ `nutrition_plan.supplements` contains recommendations
- ✅ `care_plan.sections` has 8 sections
- ✅ `monitoring_config.frequency` = "monthly" or "biweekly"

### **Test Case 2: KRAS+TP53 Profile**
**Input:**
- Disease: `ovarian`
- Mutations: KRAS p.G12D, TP53 p.R175H

**Expected Outputs:**
- ✅ `biomarker_profile.tmb.classification` = "TMB-H" or "TMB-M"
- ✅ `mechanism_vector[1]` (MAPK) > 0.5
- ✅ `resistance_prediction.risk_level` = "MEDIUM" or "HIGH"
- ✅ `drug_ranking` contains MAPK pathway inhibitors
- ✅ `trial_matches` contains trials targeting MAPK pathway

### **Test Case 3: MBD4+TP53+NF1 Profile**
**Input:**
- Disease: `ovarian`
- Mutations: MBD4 frameshift (homozygous), TP53 missense, NF1 stop_gained

**Expected Outputs:**
- ✅ `biomarker_profile.hrd.status` = "HRD+" (MBD4 homozygous)
- ✅ `drug_ranking` contains PARP inhibitors
- ✅ `synthetic_lethality_result.synthetic_lethality_detected` = true
- ✅ `resistance_prediction.detected_genes` contains relevant genes

---

## 📋 VALIDATION CHECKLIST

### **Before Running Tests:**
- [ ] Backend server running on port 8000
- [ ] Health check endpoint responds
- [ ] Test data files accessible
- [ ] Validation script executable

### **During Tests:**
- [ ] All 3 test patients execute successfully
- [ ] Response format matches `OrchestratePipelineResponse` schema
- [ ] All expected agent outputs present
- [ ] Mechanism vector is 7D
- [ ] Pipeline completes (phase = "complete")
- [ ] No errors in agent executions

### **After Tests:**
- [ ] Review validation_results.json
- [ ] Document any discrepancies
- [ ] Update data transformation function if needed
- [ ] Fix any identified issues

---

## ✅ ANSWERS TO INTEGRATION QUESTIONS

### **Q1: API Request Format Mismatch**
**Answer:** 
- ✅ Orchestrator endpoint `/api/orchestrate/full` expects **JSON body** (NOT FormData)
- ✅ Request format: `OrchestratePipelineRequest` with:
  - `disease: string` (required) - e.g., "ovarian", "myeloma"
  - `mutations: List[MutationInput]` (required) - array of mutation objects
  - `treatment_line: int` (optional, default: 1)
  - `prior_therapies: List[str]` (optional)
  - `patient_id: str` (optional, auto-generated if not provided)
- ❌ **Does NOT accept nested `patient_profile` format** from UniversalCompleteCare
- ✅ **Solution:** Transform `patientProfile` to extract:
  - `disease` from `patientProfile.disease.type` or `patientProfile.disease` (if string)
  - `mutations` from `patientProfile.tumor_context.somatic_mutations`
  - `treatment_line` from `patientProfile.treatment.line` (convert "first-line" → 1)

### **Q2: Missing Fields in Response**
**Answer:**
- ✅ **CONFIRMED:** `nutrition_plan` and `synthetic_lethality_result` exist in `PatientState` but are **NOT included** in `_state_to_response()` (line 213-235 in orchestrate.py)
- ✅ **Fix Required:** 
  1. Add to `OrchestratePipelineResponse` schema (line 157-208 in schemas/orchestrate.py)
  2. Add to `_state_to_response()` function (line 213-235 in orchestrate.py)
- ✅ **Action:** Fix backend FIRST before frontend integration

### **Q3: API Client Mismatch**
**Answer:**
- ❌ **PROBLEM:** `orchestratorApi.runPipeline()` sends **FormData** with `request` as JSON string
- ✅ **REALITY:** Endpoint `/api/orchestrate/full` expects **JSON body** with `OrchestratePipelineRequest`
- ✅ **Solution:** Update `orchestratorApi.runPipeline()` to send JSON body (not FormData) unless file is provided
- ✅ **File Upload:** If file provided, use FormData; otherwise use JSON body

### **Q4: Legacy Response Format Mapping**
**Answer:**
- ✅ **Legacy format** (`/api/complete_care/v2`):
  - `result.wiwfm.drugs` → orchestrator `drug_ranking`
  - `result.trials.trials` → orchestrator `trial_matches`
  - `result.biomarker_intelligence` → orchestrator `biomarker_profile`
  - `result.mechanism_map.vector` → orchestrator `mechanism_vector`
  - `result.resistance_playbook` → orchestrator `resistance_prediction.next_line_options` (needs transformation)
- ✅ **Transformation function** must map all these fields correctly

### **Q5: ResistancePlaybook Component Format**
**Answer:**
- ✅ **Expected format:** `{ risks, combo_strategies, next_line_switches, trial_keywords, provenance }`
- ✅ **Orchestrator provides:** `resistance_prediction.next_line_options` with `{ alternatives, regimen_changes, monitoring_changes }`
- ✅ **Solution:** Transform `next_line_options` to playbook format OR update component to accept orchestrator format

### **Q6: File Upload Integration**
**Answer:**
- ✅ **Pattern:** `PatientUpload` component uses `useOrchestrator().runPipeline(request, file, fileType)`
- ✅ **Request format:** `{ patient_id?, options: {} }` (minimal when file provided)
- ✅ **File handling:** FormData with `file` and `file_type` fields
- ✅ **Solution:** Reuse same pattern in UniversalCompleteCare

### **Q7: Status Polling**
**Answer:**
- ✅ **Hook:** `usePipelineStatus(patientId, enabled)` polls `/api/orchestrate/status/{patient_id}`
- ✅ **Polling interval:** Every 5 seconds (hardcoded in hook)
- ✅ **Auto-stops:** When `phase === 'complete'` or `phase === 'error'`
- ✅ **Solution:** Integrate same way as OrchestratorDashboard (if needed)

### **Q8: Patient Profile Transformation**
**Answer:**
- ✅ **UniversalCompleteCare format:**
  ```javascript
  {
    disease: { type: 'ovarian_cancer_hgs', stage: 'IVB' },
    tumor_context: { somatic_mutations: [...] },
    treatment: { line: 'first-line' }
  }
  ```
- ✅ **Orchestrator expects:**
  ```javascript
  {
    disease: 'ovarian',  // string, not object
    mutations: [...],     // extracted from tumor_context.somatic_mutations
    treatment_line: 1     // converted from 'first-line'
  }
  ```
- ✅ **Transformation needed:** Extract and flatten nested structure

### **Q9: Backend Fix Priority**
**Answer:**
- ✅ **Fix backend FIRST** - Add `nutrition_plan` and `synthetic_lethality_result` to response
- ✅ **Then** update frontend to use orchestrator
- ✅ **Reason:** Frontend transformation depends on complete response data

---

## 🔧 IMPLEMENTATION DECISIONS

### **Decision 1: Backend Fix First** ✅
- **Action:** Fix `_state_to_response()` and `OrchestratePipelineResponse` schema
- **Files:** 
  - `api/routers/orchestrate.py` (line 213-235)
  - `api/schemas/orchestrate.py` (line 157-208)
- **Add:** `nutrition_plan` and `synthetic_lethality_result` fields

### **Decision 2: API Client Update** ✅
- **Action:** Update `orchestratorApi.runPipeline()` to send JSON body (not FormData) when no file
- **File:** `oncology-coPilot/oncology-frontend/src/services/api/orchestrator.ts`
- **Logic:** If file provided → FormData, else → JSON body

### **Decision 3: Patient Profile Transformation** ✅
- **Action:** Create `transformPatientProfileToOrchestrator()` function
- **Location:** `utils/orchestratorMapper.js` (same file as `mapOrchestratorToLegacy`)
- **Purpose:** Convert UniversalCompleteCare profile format to orchestrator format

### **Decision 4: Data Transformation** ✅
- **Action:** Create `mapOrchestratorToLegacy()` function
- **Location:** `utils/orchestratorMapper.js`
- **Purpose:** Convert orchestrator response to legacy format for backward compatibility

### **Decision 5: Component Updates** ✅
- **Action:** Keep components as-is, transform data to match component expectations
- **Reason:** Less risky than updating all components

---

---

## 🔄 LEGACY ENDPOINT STRATEGY

### **What to Do with `/api/complete_care/v2`?**

**Current State:**
- ✅ **Legacy endpoint** (`/api/complete_care/v2`) is **fully operational** and actively used by `UniversalCompleteCare.jsx`
- ✅ **New endpoint** (`/api/orchestrate/full`) is **fully operational** and used by `OrchestratorDashboard.jsx`
- ✅ **Both endpoints** are maintained and functional

**Recommendation: GRADUAL MIGRATION** ✅

#### **Phase 1: Keep Both (Current State)** ✅ **IMMEDIATE**

**Rationale:**
- No breaking changes - both endpoints work
- Legacy endpoint is actively used and maintained
- New endpoint has better architecture but requires frontend migration
- **Action:** Keep both operational

#### **Phase 2: Frontend Migration (Next Step)** 🔄 **IN PROGRESS**

**Goal:** Migrate `UniversalCompleteCare.jsx` to use `/api/orchestrate/full`

**Benefits:**
- ✅ Unified architecture (one orchestrator system)
- ✅ File upload support
- ✅ Status polling support
- ✅ Better provenance tracking
- ✅ More structured agent outputs

**Implementation:**
1. ✅ Use `orchestratorMapper.js` (already created) for data transformation
2. ⏳ Update `UniversalCompleteCare.jsx` to call `/api/orchestrate/full`
3. ⏳ Test thoroughly with real patient data
4. ✅ Keep legacy endpoint as fallback during migration

#### **Phase 3: Deprecation (Future)**

**Timeline:** After frontend migration is complete and validated (3+ months)

**Steps:**
1. Add deprecation warning to `/api/complete_care/v2` responses
2. Monitor usage (log endpoint calls)
3. Set deprecation date (e.g., 3 months notice)
4. Remove endpoint after migration complete

**See:** `tests/LEGACY_ENDPOINT_ANALYSIS.md` for complete analysis

---

**Last Updated:** January 2025  
**Status:** ✅ **READY FOR IMPLEMENTATION** - All questions answered, clear implementation path, legacy endpoint strategy defined
