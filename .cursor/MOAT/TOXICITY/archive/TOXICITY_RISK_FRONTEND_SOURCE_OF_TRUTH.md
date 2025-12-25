# ⚠️ TOXICITY RISK ASSESSMENT - FRONTEND IMPLEMENTATION SOURCE OF TRUTH

**Purpose:** Single source of truth for toxicity risk frontend implementation  
**Date:** January 28, 2025  
**Status:** 🎯 **85% COMPLETE** - Backend 100%, Orchestrator 100%, Frontend 85%  
**Priority:** P0 (Critical - Blocks Product Launch)  
**Last Updated:** January 28, 2025 (Merged with TOXICITY_RISK_FRONTEND_AUDIT.md)

---

## 📋 EXECUTIVE SUMMARY

**What We Have:**
- ✅ Backend: 100% complete (`/api/safety/toxicity_risk`, three-factor model, mitigating foods)
- ✅ **Orchestrator Integration: 100% complete** (Jan 28, 2025)
  - ✅ Toxicity agent created and wired to analysis phase
  - ✅ PatientState updated with `toxicity_assessments` field
  - ✅ Care plan agent auto-consumes toxicity data
  - ✅ Integration tests created (7 test cases)
- ✅ Frontend Components: **85% complete** (enhanced with LLM, mitigating foods)
  - ✅ ToxicityRiskCard: Enhanced with mitigating foods + LLM explanations
  - ✅ Standalone Page: Created (`ToxicityRiskAssessment.jsx`)
  - ✅ Routes: Added to App.jsx (`/toxicity-risk`, `/toxicity-risk/:patientId`)
  - ✅ useToxicity Hook: Complete
  - ✅ useToxicityLLM Hook: **NEW** - AI-powered explanations
  - ⚠️ ToxicityChip: Wired to API (completed)
  - ⚠️ UniversalCompleteCare: Toxicity section exists but needs verification

**What We Need:**
1. ⚠️ Verify UniversalCompleteCare displays toxicity section correctly
2. ⚠️ Complete Care Plan Frontend Display (backend ready, needs UI verification)
3. ⚠️ Export functionality (PDF, JSON) - Not implemented
4. ⚠️ Multi-drug comparison view enhancements

**Estimated Total Time Remaining:** 4-6 hours (frontend polish + verification)

---

## ✅ CURRENT STATE (VERIFIED - Jan 28, 2025)

### **Backend (100% Complete)** ✅

| Component | Status | Location | Verified |
|-----------|--------|----------|----------|
| API Endpoint | ✅ | `api/routers/safety.py` | `/api/safety/toxicity_risk` |
| Safety Service | ✅ | `api/services/safety_service.py` | Three-factor model |
| Pathway Mappings | ✅ | `api/services/toxicity_pathway_mappings.py` | 30+ pharmacogenes, 11 MoA |
| Mitigating Foods | ✅ | `toxicity_pathway_mappings.py` | `get_mitigating_foods()` |
| Schemas | ✅ | `api/schemas/safety.py` | Request/Response models |
| Orchestrator Integration | ✅ | `api/services/orchestrator/orchestrator.py` | `_run_toxicity_risk_agent()` wired to analysis phase |
| PatientState | ✅ | `api/services/orchestrator/state.py` | `toxicity_assessments` field added |
| Care Plan Integration | ✅ | `orchestrator.py` | Auto-consumes toxicity in care plan generation |

**Backend Capabilities:**
- ✅ Risk score (0-1)
- ✅ Risk level (HIGH/MODERATE/LOW)
- ✅ Contributing factors
- ✅ Confidence adjustment
- ✅ Complete provenance
- ✅ Mitigating foods (THE MOAT)

---

### **Orchestrator Integration (100% Complete)** ✅ **NEW - Jan 28, 2025**

| Component | Status | Location | Notes |
|-----------|--------|----------|-------|
| Toxicity Agent | ✅ | `orchestrator.py:743` | `_run_toxicity_risk_agent()` method |
| Analysis Phase Wiring | ✅ | `orchestrator.py:256-305` | Runs in parallel with biomarker/resistance/nutrition |
| PatientState Field | ✅ | `state.py:158` | `toxicity_assessments: Optional[Dict]` |
| State Serialization | ✅ | `state.py:to_full_dict()` | Includes toxicity_assessments |
| Care Plan Integration | ✅ | `orchestrator.py:1199-1350` | Section 5: Toxicity Risk Assessment |
| Integration Tests | ✅ | `tests/test_toxicity_orchestrator_integration.py` | 7 test cases covering full pipeline |

**Orchestrator Capabilities:**
- ✅ Automatically runs toxicity assessment in analysis phase
- ✅ Extracts germline variants from PatientState
- ✅ Gets drugs from patient profile or drug ranking
- ✅ Assesses each drug for toxicity risk
- ✅ Stores results in `state.toxicity_assessments`
- ✅ Includes toxicity section in care plan (Section 5)
- ✅ Shows mitigating foods in care plan
- ✅ Handles errors gracefully (doesn't break pipeline)

**What This Means:**
- ✅ **Backend orchestrator integration is COMPLETE**
- ✅ Toxicity risk assessment runs automatically when patient has germline variants + drugs
- ✅ Care plan includes toxicity section with high-risk drugs flagged
- ⚠️ **Frontend display needs verification** (see Phase 3 below)

---

### **Frontend (85% Complete)** ✅ **ENHANCED - Jan 28, 2025**

| Component | Status | Location | Notes |
|-----------|--------|----------|-------|
| ToxicityRiskCard | ✅ **ENHANCED** | `components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx` | **NEW:** Mitigating foods display + LLM explanations |
| useToxicity Hook | ✅ Working | `components/ClinicalGenomicsCommandCenter/hooks/useToxicity.js` | Calls API correctly |
| useToxicityLLM Hook | ✅ **NEW** | `components/ClinicalGenomicsCommandCenter/hooks/useToxicityLLM.js` | **NEW:** AI-powered explanations (clinician/patient/researcher) |
| ToxicityChip | ✅ **WIRED** | `components/vus/ToxicityChip.jsx` | **UPDATED:** Now calls API, shows dynamic risk |
| Standalone Page | ✅ **CREATED** | `pages/ToxicityRiskAssessment.jsx` | **NEW:** Full standalone page with multi-drug support |
| Routes | ✅ **ADDED** | `App.jsx` | Routes: `/toxicity-risk`, `/toxicity-risk/:patientId` |
| CoPilot Integration | ✅ Complete | `integrations/ClinicalGenomicsCoPilotIntegration.jsx` | 3 quick actions |
| UniversalCompleteCare | ⚠️ **PARTIAL** | `pages/UniversalCompleteCare.jsx` | Toxicity section exists, needs verification |

**Current Usage:**
- ✅ Used in `MechanisticEvidenceTab.jsx` (ClinicalGenomicsCommandCenter)
- ✅ Used in `AnalysisResults.jsx` (VUS analysis - now wired)
- ✅ Used in `ToxicityRiskAssessment.jsx` (standalone page)
- ⚠️ Used in `UniversalCompleteCare.jsx` (needs verification)

**Recent Enhancements (Jan 28, 2025):**
1. ✅ **ToxicityRiskCard Enhanced:**
   - Added mitigating foods display (THE MOAT)
   - Added LLM-powered explanations (3 audience types)
   - Added audience selector (clinician/patient/researcher)
   - Added collapsible explanation section

2. ✅ **Standalone Page Created:**
   - Full-page layout with patient input form
   - Multi-drug assessment support
   - Comparison table view
   - Real-time assessment

3. ✅ **ToxicityChip Wired:**
   - Replaced placeholder with actual API call
   - Dynamic risk level display
   - Detailed tooltip with factors and mitigating foods

---

## 🎯 WHAT WE'RE BUILDING (REMAINING WORK)

### **Vision (From ADVANCED_CARE_PLAN_TOXCITY.md)**

> **"What should I eat to protect myself from THIS specific drug's side effects?"**

**The Answer:**
- Your carboplatin + BRCA1 = DNA repair stress
- NAC helps - it boosts glutathione which supports DNA repair
- Take 600mg twice daily AFTER infusion, not during
- Here's why this matters for YOU

**The MOAT:** No competitor answers this question. We connect:
- Toxicity Detection → Knows your drug damages DNA repair pathways
- Food Validation → Knows NAC supports DNA repair
- The Bridge → Personalized recommendations with timing guidance

---

## 📋 IMPLEMENTATION PLAN (UPDATED STATUS)

### **PHASE 1: ToxicityRiskCard Enhancement** ✅ **COMPLETE** (Jan 28, 2025)

**Status:** ✅ **COMPLETE**

**What Was Built:**
1. ✅ **Mitigating Foods Display** (2-3 hours) - **DONE**
   - Displays `result.mitigating_foods` array
   - Shows: compound, dose, timing, mechanism
   - Format: List with clear labels

2. ✅ **LLM-Powered Explanations** (2-3 hours) - **NEW - DONE**
   - Added `useToxicityLLM` hook
   - Audience selector (clinician/patient/researcher)
   - Collapsible explanation section
   - Error handling and loading states

3. ⚠️ **Prominent Pharmacogene Warnings** (1-2 hours) - **PARTIAL**
   - High-impact pharmacogenes detected
   - ⚠️ Needs: Red alert styling, dose adjustment recommendations

4. ❌ **Export Functionality** (1 hour) - **NOT DONE**
   - PDF export
   - JSON export
   - Shareable link

**Success Criteria:**
- [x] Mitigating foods displayed when present ✅
- [x] Timing guidance visible ("post-chemo, not during") ✅
- [x] LLM explanations working ✅
- [ ] High-impact pharmacogenes show red alert ⚠️
- [ ] Export buttons functional ❌

**Files Modified:**
- ✅ `components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx` - Enhanced
- ✅ `components/ClinicalGenomicsCommandCenter/hooks/useToxicityLLM.js` - Created

---

### **PHASE 2: Standalone Toxicity Risk Page** ✅ **COMPLETE** (Jan 28, 2025)

**Status:** ✅ **COMPLETE**

**What Was Built:**
1. ✅ **Page Component** (4-6 hours) - **DONE**
   - File: `pages/ToxicityRiskAssessment.jsx`
   - Patient input form:
     - Germline variants (manual entry/VCF upload)
     - Drug selection (single/multi-select with MoA mapping)
     - Disease context dropdown
     - Treatment line (optional)
   - Real-time assessment on form submit
   - Results display using enhanced ToxicityRiskCard

2. ✅ **Route Addition** (30 minutes) - **DONE**
   - File: `App.jsx`
   - Routes:
     - `/toxicity-risk` (standalone)
     - `/toxicity-risk/:patientId` (with patient context)

3. ✅ **Multi-Drug Comparison View** (3-4 hours) - **DONE**
   - Comparison table (if multiple drugs selected)
   - Columns: Drug Name, Risk Score, Risk Level, Key Factors, Mitigating Foods
   - Risk ranking (lowest to highest)
   - Expandable rows for detailed view

**Success Criteria:**
- [x] Page accessible at `/toxicity-risk` ✅
- [x] User can input germline variants (manual/VCF) ✅
- [x] User can select single or multiple drugs ✅
- [x] Real-time assessment works ✅
- [x] Results display correctly (single and multi-drug) ✅
- [ ] Export functionality works ❌

**Files Created:**
- ✅ `pages/ToxicityRiskAssessment.jsx` - Created

**Files Modified:**
- ✅ `App.jsx` - Routes added

---

### **PHASE 3: Complete Care Plan Integration** ⚠️ **BACKEND COMPLETE, FRONTEND NEEDS VERIFICATION**

**Status:** ✅ **Backend integration complete** (Jan 28, 2025), ⚠️ **Frontend display needs verification**

**What Exists:**
- ✅ Orchestrator automatically assesses toxicity in analysis phase
- ✅ Care plan agent includes toxicity section (Section 5)
- ✅ `state.toxicity_assessments` populated automatically
- ⚠️ **Frontend display needs verification**

**Tasks:**
1. ✅ **Backend Integration** ✅ **COMPLETE**
   - ✅ File: `api/services/orchestrator/orchestrator.py`
   - ✅ `_run_toxicity_risk_agent()` method created
   - ✅ Wired to analysis phase (runs in parallel)
   - ✅ Care plan includes toxicity section automatically
   - ✅ `state.toxicity_assessments` field added to PatientState

2. ⚠️ **Frontend Display** (3-4 hours) ⚠️ **NEEDS VERIFICATION**
   - File: `pages/UniversalCompleteCare.jsx`
   - ⚠️ Toxicity section exists (imports ToxicityRiskCard)
   - ⚠️ Needs verification: Does it display correctly?
   - ⚠️ Needs verification: Are mitigating foods shown?
   - ⚠️ Needs verification: Are high-risk drugs flagged?

**Success Criteria:**
- [x] Complete Care Plan calls toxicity risk assessment ✅ (backend)
- [ ] Toxicity risks displayed for all recommended drugs ⚠️ (needs verification)
- [ ] Mitigating foods shown in care plan summary ⚠️ (needs verification)
- [ ] High-risk drugs flagged prominently ⚠️ (needs verification)
- [ ] Link to detailed toxicity assessment works ⚠️ (needs verification)

**Files to Verify:**
- ⚠️ `pages/UniversalCompleteCare.jsx` - Toxicity section exists, needs verification

---

### **PHASE 4: ToxicityChip Wiring** ✅ **COMPLETE** (Jan 28, 2025)

**Status:** ✅ **COMPLETE**

**What Was Built:**
1. ✅ **Wire ToxicityChip to API** (2-3 hours) - **DONE**
   - File: `components/vus/ToxicityChip.jsx`
   - Replaced placeholder with actual API call using `useToxicity` hook
   - Shows risk level chip (HIGH/MODERATE/LOW) with color coding
   - Tooltip with details (risk score, key factors, mitigating foods)

**Success Criteria:**
- [x] ToxicityChip calls API when germline variants present ✅
- [x] Risk level chip displays correctly ✅
- [x] Tooltip shows details ✅
- [x] No errors when no variants present ✅

**Files Modified:**
- ✅ `components/vus/ToxicityChip.jsx` - Wired to API

---

## 🎨 UI/UX SPECIFICATIONS (IMPLEMENTED)

### **ToxicityRiskCard Enhancement** ✅ **IMPLEMENTED**

**Mitigating Foods Section** ✅ **DONE:**
```jsx
{result.mitigating_foods && result.mitigating_foods.length > 0 && (
  <Box sx={{ mt: 2 }}>
    <Typography variant="subtitle2" gutterBottom>
      Mitigating Foods/Supplements:
    </Typography>
    <List dense>
      {result.mitigating_foods.map((food, idx) => (
        <ListItem key={idx}>
          <ListItemText
            primary={food.compound}
            secondary={`${food.dose} - ${food.timing} | ${food.mechanism}`}
          />
        </ListItem>
      ))}
    </List>
  </Box>
)}
```

**LLM Explanation Section** ✅ **DONE:**
- Audience selector (clinician/patient/researcher)
- "Generate AI Explanation" button
- Collapsible explanation display
- Error handling

**High-Impact Pharmacogene Warnings** ⚠️ **PARTIAL:**
- ⚠️ Needs: Red Alert styling for high-impact pharmacogenes (DPYD, TPMT)
- ⚠️ Needs: Dose adjustment recommendations
- ⚠️ Needs: Alternative drug suggestions

---

### **Standalone Page Layout** ✅ **IMPLEMENTED**

**Page Structure:**
```
/toxicity-risk
├── Header
│   ├── Title: "Toxicity Risk Assessment (RUO)"
│   ├── Subtitle: "Germline-based toxicity prediction for precision safety"
│   └── RUO Disclaimer
│
├── Input Section ✅
│   ├── Patient Selection (if patientId in URL) ✅
│   ├── Germline Variants Input ✅
│   │   ├── Manual Entry ✅
│   │   └── Load from Patient Profile ✅
│   ├── Drug Selection ✅
│   │   ├── Single Drug (dropdown) ✅
│   │   ├── Multiple Drugs (multi-select) ✅
│   │   └── MoA Auto-Detection ✅
│   └── Clinical Context ✅
│       ├── Disease Selection ✅
│       └── Treatment Line (optional) ✅
│
├── Assessment Results ✅
│   ├── Single Drug View ✅
│   │   └── ToxicityRiskCard (enhanced) ✅
│   │       ├── Risk Score Visualization ✅
│   │       ├── Risk Level Chip ✅
│   │       ├── Confidence Chip ✅
│   │       ├── Contributing Factors ✅
│   │       ├── Mitigating Foods ✅
│   │       └── LLM Explanations ✅
│   │
│   └── Multi-Drug Comparison ✅
│       ├── Comparison Table ✅
│       └── Risk Ranking ✅
│
└── Actions
    ├── Export PDF ❌ (Not implemented)
    ├── Export JSON ❌ (Not implemented)
    └── Share Link ❌ (Not implemented)
```

---

## 🔗 INTEGRATION SPECIFICATIONS

### **Backend: Complete Care Plan Integration** ✅ **COMPLETE**

**File:** `api/services/orchestrator/orchestrator.py`

**Implementation:**
- ✅ `_run_toxicity_risk_agent()` method created
- ✅ Wired to `_run_analysis_phase()` (runs in parallel)
- ✅ Stores results in `state.toxicity_assessments`
- ✅ Care plan agent auto-consumes toxicity data
- ✅ Section 5: Toxicity Risk Assessment included in care plan

**What This Means:**
- ✅ Backend automatically assesses toxicity for all drugs
- ✅ Results stored in PatientState
- ✅ Care plan includes toxicity section automatically
- ⚠️ Frontend needs to display this data (see below)

---

### **Frontend: Complete Care Plan Display** ⚠️ **NEEDS VERIFICATION**

**File:** `pages/UniversalCompleteCare.jsx`

**Current Status:**
- ✅ Imports `ToxicityRiskCard` component
- ⚠️ Needs verification: Does it display toxicity section?
- ⚠️ Needs verification: Are mitigating foods shown?
- ⚠️ Needs verification: Are high-risk drugs flagged?

**Expected Implementation:**
```jsx
{result.toxicity_assessments && result.toxicity_assessments.toxicity_assessments?.length > 0 && (
  <Box sx={{ mt: 4 }}>
    <Typography variant="h5" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
      <WarningIcon color="warning" />
      Toxicity Risk Assessment
    </Typography>
    
    <Grid container spacing={2}>
      {result.toxicity_assessments.toxicity_assessments.map((risk, idx) => (
        <Grid item xs={12} key={idx}>
          <ToxicityRiskCard
            result={{
              risk_score: risk.risk_score,
              confidence: risk.confidence,
              reason: risk.reason,
              factors: risk.factors,
              mitigating_foods: risk.mitigating_foods
            }}
          />
        </Grid>
      ))}
    </Grid>
  </Box>
)}
```

---

## 📊 CAPABILITY MATRIX (UPDATED - Jan 28, 2025)

| Capability | Backend | Frontend | Standalone Page | Care Plan Integration |
|------------|---------|----------|-----------------|----------------------|
| **Risk Score Calculation** | ✅ | ✅ | ✅ | ✅ |
| **Risk Level Classification** | ✅ | ✅ | ✅ | ⚠️ |
| **Contributing Factors** | ✅ | ✅ | ✅ | ⚠️ |
| **Mitigating Foods** | ✅ | ✅ | ✅ | ⚠️ |
| **LLM Explanations** | N/A | ✅ | ✅ | ❌ |
| **Multi-Drug Assessment** | ✅ | ✅ | ✅ | ✅ |
| **Patient Input Form** | N/A | ✅ | ✅ | N/A |
| **Export Functionality** | N/A | ❌ | ❌ | ❌ |
| **Pharmacogene Warnings** | ✅ | ⚠️ Partial | ⚠️ Partial | ⚠️ |
| **Complete Care Plan Integration** | ✅ | ⚠️ | N/A | ⚠️ |

**Legend:**
- ✅ Complete
- ⚠️ Partial (needs enhancement/verification)
- ❌ Missing
- N/A Not applicable

---

## 🎯 SUCCESS CRITERIA (UPDATED)

### **Standalone Page:**
- [x] User can input germline variants (manual) ✅
- [x] User can select single or multiple drugs ✅
- [x] Real-time toxicity assessment ✅
- [x] Risk level chips (HIGH/MODERATE/LOW) with color coding ✅
- [x] Contributing factors displayed ✅
- [x] Mitigating foods displayed with timing guidance ✅
- [x] LLM explanations available ✅
- [ ] Export functionality (PDF, JSON) ❌
- [ ] Shareable link generation ❌

### **Care Plan Integration:**
- [x] Complete Care Plan calls toxicity risk assessment ✅ (backend)
- [ ] Toxicity risks displayed for all recommended drugs ⚠️ (needs verification)
- [ ] Mitigating foods shown in care plan summary ⚠️ (needs verification)
- [ ] High-risk drugs flagged prominently ⚠️ (needs verification)
- [ ] Link to detailed toxicity assessment works ⚠️ (needs verification)

### **Enhanced ToxicityRiskCard:**
- [x] Displays mitigating foods section ✅
- [x] Shows timing guidance ("post-chemo, not during") ✅
- [x] LLM explanations available ✅
- [ ] Prominent warnings for high-impact pharmacogenes ⚠️ (needs red alert styling)
- [ ] Export functionality ❌
- [x] Link to food validation (via mitigating foods) ✅

---

## 📝 IMPLEMENTATION PRIORITY (UPDATED)

### **P0 (Critical - Blocks Product Launch):**
1. ✅ Backend implementation (DONE)
2. ✅ ToxicityRiskCard enhancement (DONE)
3. ✅ Standalone page creation (DONE)
4. ⚠️ Complete Care Plan frontend verification (PENDING)

### **P1 (Important - Product Enhancement):**
5. ⚠️ Prominent pharmacogene warnings (red alert styling)
6. ❌ Export functionality (PDF, JSON)
7. ⚠️ UniversalCompleteCare toxicity display verification

### **P2 (Nice to Have):**
8. Advanced filtering (by risk level, pharmacogene type)
9. Historical tracking (risk scores over time)
10. Patient-specific recommendations based on toxicity risk

---

## 🔗 REFERENCES

- **Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`
- **Contribution Document:** `.cursor/lectures/drugDevelopment/toxicity_risk_contribution.mdc`
- **Concept Document:** `.cursor/rules/research/toxicity_risk_concept.mdc`
- **Backend API:** `api/routers/safety.py` - `/api/safety/toxicity_risk`
- **Frontend Components:** 
  - `components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx` ✅ Enhanced
  - `components/vus/ToxicityChip.jsx` ✅ Wired
  - `components/ClinicalGenomicsCommandCenter/hooks/useToxicity.js` ✅ Complete
  - `components/ClinicalGenomicsCommandCenter/hooks/useToxicityLLM.js` ✅ **NEW**
  - `pages/ToxicityRiskAssessment.jsx` ✅ **NEW**

---

## ✅ IMPLEMENTATION CHECKLIST (UPDATED - Jan 28, 2025)

### **✅ COMPLETED (Jan 28, 2025)**
- [x] ✅ Backend Orchestrator Integration (Deliverables 1-5)
  - [x] Toxicity agent created
  - [x] Wired to analysis phase
  - [x] PatientState updated
  - [x] Care plan integration complete
  - [x] Integration tests created
- [x] ✅ ToxicityRiskCard Enhancement
  - [x] Mitigating foods display
  - [x] LLM explanations
  - [x] Audience selector
- [x] ✅ Standalone Page Creation
  - [x] Page component created
  - [x] Routes added
  - [x] Multi-drug support
- [x] ✅ ToxicityChip Wiring
  - [x] API integration
  - [x] Dynamic risk display
  - [x] Tooltip with details

### **⚠️ PENDING VERIFICATION**
- [ ] UniversalCompleteCare toxicity display
- [ ] Care plan toxicity section rendering
- [ ] Mitigating foods in care plan
- [ ] High-risk drug flagging

### **❌ NOT STARTED**
- [ ] Export functionality (PDF, JSON)
- [ ] Prominent pharmacogene warnings (red alert)
- [ ] Shareable link generation

**Total Completed:** ~85% (Backend 100%, Orchestrator 100%, Frontend 85%)  
**Remaining Work:** 4-6 hours (verification + polish)

---

**Last Updated:** January 28, 2025  
**Status:** 🎯 **85% COMPLETE** - Ready for Verification & Polish  
**Next Action:** Verify UniversalCompleteCare toxicity display, add export functionality


