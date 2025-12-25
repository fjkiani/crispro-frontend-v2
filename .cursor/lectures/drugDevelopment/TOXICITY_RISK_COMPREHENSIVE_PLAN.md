# ⚠️ Toxicity Risk Integration - Comprehensive Implementation Plan

**Purpose:** Complete plan for integrating toxicity risk into orchestrator pipeline + frontend  
**Date:** January 28, 2025  
**Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`  
**Reference:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md`  
**Implementation Guide:** `.cursor/MOAT/orchestration/agent-implementation-guide.mdc`

---

## 🎯 Executive Summary

### **Current State Analysis:**

**Backend:**
- ✅ **Safety Service:** 100% complete (`api/services/safety_service.py`)
- ✅ **Pathway Mappings:** 100% complete (`api/services/toxicity_pathway_mappings.py`)
- ✅ **Mitigating Foods:** 100% complete (THE MOAT implemented)
- ✅ **API Endpoint:** `/api/safety/toxicity_risk` fully operational
- ❌ **Orchestrator Integration:** NOT integrated (no toxicity agent in pipeline)

**Frontend:**
- ✅ **ToxicityRiskCard:** Component exists and wired to API
- ✅ **useToxicity Hook:** Complete and working
- ✅ **CoPilot Integration:** 3 quick actions implemented
- ⚠️ **ToxicityRiskCard:** Missing mitigating foods display
- ❌ **Standalone Page:** Not created
- ❌ **Care Plan Integration:** Not displayed in care plan

**Orchestrator:**
- ✅ **Pattern Established:** Biomarker, Resistance, Nutrition agents in Phase 2 (parallel)
- ✅ **State Management:** PatientState has fields for agent outputs
- ❌ **Toxicity Agent:** Not in pipeline
- ❌ **State Field:** `toxicity_assessments` not in PatientState

---

## 📊 Architecture: Where Toxicity Fits

### **Current Pipeline Structure (From 00_MASTER_INDEX.mdc):**

```
Phase 1: EXTRACTING
└── [01] Data Extraction Agent
    └── Extract mutations, germline variants

Phase 2: ANALYZING (Parallel Execution) ← ADD TOXICITY HERE
├── [02] Biomarker Agent (TMB, MSI, HRD)
├── [03] Resistance Agent (MAPK, DDR pathways)
└── [06] Nutrition Agent (Food validation, timing rules)

Phase 3: RANKING
├── [04] Drug Efficacy Agent (S/P/E framework)
└── [14] Synthetic Lethality Agent

Phase 4: MATCHING
└── [05] Trial Matching Agent

Phase 5: PLANNING
└── [07] Care Plan Agent (Aggregates all outputs)

Phase 6: MONITORING
└── [08] Monitoring Agent
```

### **Enhanced Pipeline (With Toxicity):**

```
Phase 2: ANALYZING (Parallel Execution) ← ADD TOXICITY HERE
├── [02] Biomarker Agent ✅
├── [03] Resistance Agent ✅
├── [06] Nutrition Agent ✅
└── [16] Toxicity Risk Agent ⚠️ NEW ← Module 16
    └── Assess toxicity for recommended drugs
    └── Return risk scores, factors, mitigating foods
```

**Key Insight:** Toxicity risk should run in Phase 2 (parallel with biomarker/resistance/nutrition) because:
1. It only needs germline variants (from Phase 1)
2. It can run independently of other agents
3. Results are needed by Care Plan Agent (Phase 5)
4. It can assess drugs from patient profile OR use top drugs from ranking if available (flexible)

---

## 🏗️ Phase 1: Backend Orchestrator Integration (8-10 hours)

### **1.1 Create Toxicity Agent Module**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/toxicity_agent.py` (NEW)

**Pattern:** Follow `_run_biomarker_agent` or `_run_resistance_agent` pattern from agent-implementation-guide.mdc

**Key Requirements:**
- ✅ Use direct service imports (NO HTTP calls)
- ✅ Import `get_safety_service()` from `api.services.safety_service`
- ✅ Use `ToxicityRiskRequest`, `PatientContext`, `GermlineVariant`, `TherapeuticCandidate`, `ClinicalContext` from `api.schemas.safety`
- ✅ Use `get_drug_moa()` from `api.services.toxicity_pathway_mappings`
- ✅ Return dict matching orchestrator output schema
- ✅ Handle errors gracefully (log and continue)

**See:** `.cursor/MOAT/orchestration/16_TOXICITY_RISK_AGENT_INTEGRATION_PLAN.md` Section 1.1 for full implementation

**Time Estimate:** 2-3 hours

---

### **1.2 Wire to Orchestrator**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/orchestrator.py`

**Location 1: Add to `_run_analysis_phase` method (around line 256-305)**

**Pattern (From agent-implementation-guide.mdc lines 315-348):**
- Add `toxicity_risk` to parallel execution tasks
- Store result in `state.toxicity_assessments`

**Location 2: Add `_run_toxicity_risk_agent` method (after `_run_nutrition_agent`)**

**Pattern (From agent-implementation-guide.mdc lines 617-664):**
- Start execution tracking
- Import agent class
- Extract data from PatientState
- Call agent method
- Handle errors gracefully

**See:** `.cursor/MOAT/orchestration/16_TOXICITY_RISK_AGENT_INTEGRATION_PLAN.md` Section 1.2 for full implementation

**Time Estimate:** 1-2 hours

---

### **1.3 Update PatientState**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/state.py`

**Add field:**
```python
toxicity_assessments: Optional[Dict] = None  # From 16_TOXICITY_RISK
```

**Also update:**
- `to_full_dict()` method
- `StateStore._serialize()` method

**Time Estimate:** 30 minutes

---

### **1.4 Update Care Plan Agent**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/orchestrator.py`

**Auto-consume toxicity from state:**
```python
if state.toxicity_assessments:
    sections.append({
        'title': 'Toxicity Risk Assessment',
        'content': state.toxicity_assessments
    })
```

**Time Estimate:** 30 minutes

---

## 🎨 Phase 2: Frontend Standalone Page (8-12 hours)

### **2.1 Create Standalone Page**

**File:** `oncology-coPilot/oncology-frontend/src/pages/ToxicityRiskAssessment.jsx` (NEW)

**Reuse Components:**
- ✅ `useToxicity` hook (already exists)
- ✅ `ToxicityRiskCard` component (will enhance)
- ✅ Input patterns from other pages

**See:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md` lines 310-361 for full specification

**Time Estimate:** 4-6 hours

---

### **2.2 Enhance ToxicityRiskCard**

**File:** `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx`

**Add:**
1. Mitigating Foods Section (THE MOAT) - See TOXICITY_RISK_FRONTEND_AUDIT.md lines 369-401
2. Prominent Pharmacogene Warnings - See TOXICITY_RISK_FRONTEND_AUDIT.md lines 403-417
3. Export functionality (future)

**Time Estimate:** 1-2 hours

---

## 🔗 Phase 3: Complete Care Plan Integration (3-4 hours)

### **3.1 Backend: Add to Complete Care Plan Universal**

**File:** `oncology-coPilot/oncology-backend-minimal/api/routers/complete_care_universal.py`

**Add `_assess_toxicity_risks()` function and integrate into `get_complete_care_v2()`**

**See:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md` lines 423-513

**Time Estimate:** 2-3 hours

---

### **3.2 Frontend: Display in Care Plan**

**File:** `oncology-coPilot/oncology-frontend/src/pages/AyeshaCompleteCare.jsx` (or UniversalCarePlan page)

**Add toxicity risk section to care plan display**

**See:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md` lines 517-543

**Time Estimate:** 1 hour

---

## 📋 Implementation Checklist

### **Backend (Orchestrator Integration):**
- [ ] Create `toxicity_agent.py` (2-3 hours)
- [ ] Add `_run_toxicity_risk_agent()` to orchestrator (1-2 hours)
- [ ] Update `PatientState` with `toxicity_assessments` field (30 min)
- [ ] Update `StateStore._serialize()` to include toxicity (15 min)
- [ ] Update `PatientState.to_full_dict()` to include toxicity (15 min)
- [ ] Update care plan agent to consume toxicity (30 min)
- [ ] Add to `_run_analysis_phase` parallel execution (30 min)
- [ ] Unit tests for toxicity agent (2 hours)
- [ ] Integration test with orchestrator (1 hour)

**Total Backend Time:** 8-10 hours

### **Frontend (Standalone Page + Enhancements):**
- [ ] Create `ToxicityRiskAssessment.jsx` page (4-6 hours)
- [ ] Add route to `App.jsx` (15 min)
- [ ] Enhance `ToxicityRiskCard` with mitigating foods (1-2 hours)
- [ ] Add prominent pharmacogene warnings (1 hour)
- [ ] Add export functionality (PDF, JSON) (2-3 hours)

**Total Frontend Time:** 8-12 hours

### **Care Plan Integration:**
- [ ] Backend: Add `_assess_toxicity_risks()` to complete care (2-3 hours)
- [ ] Frontend: Display toxicity risks in care plan (1 hour)

**Total Care Plan Time:** 3-4 hours

---

## 🎯 Quick Win Strategy

### **Phase 1: Backend Orchestrator (8-10 hours)**
**Why Quick Win:**
- ✅ Backend service 100% complete (reuse existing)
- ✅ Follow established patterns (biomarker, resistance)
- ✅ No new dependencies needed
- ✅ Immediate pipeline integration

**Deliverables:**
- Toxicity risk assessed for all recommended drugs
- Available in orchestrator state
- Auto-consumed by care plan

### **Phase 2: Frontend Standalone Page (4-6 hours)**
**Why Quick Win:**
- ✅ Reuse existing `useToxicity` hook
- ✅ Reuse existing `ToxicityRiskCard` component
- ✅ Simple input form (drug selection)
- ✅ Immediate user value

**Deliverables:**
- Standalone `/toxicity-risk` page
- Single/multi-drug assessment
- Results display

### **Phase 3: Enhancements (8-12 hours)**
**Why Quick Win:**
- ✅ Incremental improvements
- ✅ Each enhancement independent
- ✅ Can ship after Phase 1+2

**Deliverables:**
- Mitigating foods display
- Pharmacogene warnings
- Export functionality
- Care plan integration

---

## ✅ Success Criteria (From TOXICITY_RISK_FRONTEND_AUDIT.md)

### **Standalone Page:**
- [x] User can input germline variants (VCF or manual)
- [x] User can select single or multiple drugs
- [x] Real-time toxicity assessment
- [x] Risk level chips (HIGH/MODERATE/LOW) with color coding
- [x] Contributing factors displayed
- [x] Mitigating foods displayed with timing guidance
- [x] Export functionality (PDF, JSON)
- [x] Shareable link generation

### **Care Plan Integration:**
- [x] Complete Care Plan calls toxicity risk assessment
- [x] Toxicity risks displayed for all recommended drugs
- [x] Mitigating foods shown in care plan summary
- [x] High-risk drugs flagged prominently
- [x] Link to detailed toxicity assessment page

### **Enhanced ToxicityRiskCard:**
- [x] Displays mitigating foods section
- [x] Shows timing guidance ("post-chemo, not during")
- [x] Prominent warnings for high-impact pharmacogenes
- [x] Export functionality
- [x] Link to food validation for detailed recommendations

---

## 🔗 References

- **Module Spec:** `.cursor/MOAT/orchestration/16_TOXICITY_RISK_AGENT.mdc`
- **Integration Plan:** `.cursor/MOAT/orchestration/16_TOXICITY_RISK_AGENT_INTEGRATION_PLAN.md`
- **Frontend Audit:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md`
- **Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`
- **Agent Implementation Guide:** `.cursor/MOAT/orchestration/agent-implementation-guide.mdc`
- **Master Index:** `.cursor/MOAT/orchestration/00_MASTER_INDEX.mdc`
- **Backend Service:** `api/services/safety_service.py`
- **Frontend Hook:** `components/ClinicalGenomicsCommandCenter/hooks/useToxicity.js`

---

## 📊 What to Reuse

### **Backend:**
- ✅ `api/services/safety_service.py` - 100% complete
- ✅ `api/services/toxicity_pathway_mappings.py` - Complete
- ✅ `api/schemas/safety.py` - Complete
- ✅ Orchestrator patterns (biomarker, resistance, nutrition)

### **Frontend:**
- ✅ `useToxicity` hook - Complete
- ✅ `ToxicityRiskCard` component - Needs enhancement
- ✅ Input patterns from other pages
- ✅ Route structure from App.jsx

### **Architecture:**
- ✅ Parallel execution pattern (Phase 2)
- ✅ State management pattern (PatientState)
- ✅ Care plan consumption pattern (auto-reads from state)

---

**Last Updated:** January 28, 2025  
**Status:** ✅ Comprehensive Plan Complete - Ready for Implementation  
**Total Estimated Time:** 19-26 hours (can be done incrementally)


