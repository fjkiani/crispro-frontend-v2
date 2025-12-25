# 🎯 THERAPY FIT - PRODUCTION READY IMPLEMENTATION PLAN

**Date:** January 2025  
**Status:** Production-Ready Plan  
**Inspired By:** Advanced Care Plan Universalization MOAT  
**Last Updated:** January 2025 *(Initial comprehensive plan)*

---

## 🏆 EXECUTIVE SUMMARY

> **The question:** "Which drugs are most likely to work for this patient's specific mutations?"  
> **Answer:** **Ranked drug list with transparent S/P/E scoring, confidence, evidence tiers, and badges.**

| Before | After |
|--------|-------|
| "Try these drugs (no ranking)" | "Ranked by S/P/E: 1. Drug A (0.75 confidence), 2. Drug B (0.68 confidence)" |
| Embedded in MyelomaDigitalTwin only | Standalone `/therapy-fit` page for any cancer type |
| "Why this drug?" → "Clinical judgment" | "Why this drug?" → "S: 0.85, P: 0.70, E: 0.65 (RCT evidence)" |

**Goal:** Transform Therapy Fit from a component embedded in `MyelomaDigitalTwin` into a **first-class, universal, production-ready feature** accessible to any cancer patient.

---

## ✅ CURRENT STATUS: WHAT'S BUILT

### **Backend: 100% Complete** ✅

| Component | Status | Location | Validation |
|-----------|--------|----------|------------|
| **S/P/E Orchestrator** | ✅ Done | `api/services/efficacy_orchestrator/orchestrator.py` | ✅ Tested |
| **Sequence Scorer (Evo2)** | ✅ Done | `api/services/efficacy_orchestrator/sequence_processor.py` | ✅ Tested |
| **Pathway Aggregator** | ✅ Done | `api/services/pathway/aggregation.py` | ✅ 100% MM accuracy |
| **Evidence Scorer** | ✅ Done | `api/services/evidence/literature_client.py` | ✅ Tested |
| **Drug Scorer** | ✅ Done | `api/services/efficacy_orchestrator/drug_scorer.py` | ✅ Tested |
| **API Endpoint** | ✅ Done | `api/routers/efficacy/router.py` | ✅ Live |
| **Insights Integration** | ✅ Done | 4 chips with confidence lifts | ✅ Tested |
| **Provenance Tracking** | ✅ Done | Complete audit trail | ✅ Tested |

**Endpoint:** `POST /api/efficacy/predict`

**Tested Use Cases:**
- ✅ BRAF V600E (melanoma) → BRAF inhibitors rank #1
- ✅ KRAS G12D (myeloma) → MEK > BRAF ranking (correct)
- ✅ BRCA1/BRCA2 (ovarian) → Platinum agents rank #1
- ✅ TP53 R248Q (ovarian) → DNA repair pathway drugs

### **Frontend: 80% Complete** ⚠️

| Component | Status | Location | Notes |
|-----------|--------|----------|-------|
| **EfficacyPanel Component** | ✅ Done | `src/features/efficacy/components/EfficacyPanel.jsx` | Reusable |
| **EfficacyCard Component** | ✅ Done | `src/features/efficacy/components/EfficacyCard.jsx` | Reusable |
| **useEfficacy Hook** | ✅ Done | `src/features/efficacy/useEfficacy.js` | Reusable |
| **MyelomaDigitalTwin Integration** | ✅ Done | `src/pages/MyelomaDigitalTwin.jsx` | Working |
| **Standalone Therapy Fit Page** | ❌ Missing | N/A | **NEEDED** |
| **Disease Selector** | ⚠️ Partial | May need creation/enhancement | **NEEDED** |
| **S/P/E Explanation Section** | ❌ Missing | N/A | **NEEDED** |

**Current Usage:**
- ✅ Used in `MyelomaDigitalTwin` (disease-specific)
- ✅ Used in `OrchestratorDashboard` (orchestrator context)
- ❌ **No standalone page** for direct access

---

## 🎯 WHAT WE NEED: PRODUCTION-READY STANDALONE PAGE

### **The Gap**

**Current State:**
- S/P/E framework works perfectly ✅
- Components are reusable ✅
- But: No standalone "Therapy Fit" page to showcase this capability

**What's Missing:**
1. **Standalone Route**: `/therapy-fit` page (not embedded in other pages)
2. **Disease Agnostic**: Works for any cancer type (not just myeloma)
3. **S/P/E Explanation**: Visual breakdown of how S/P/E works
4. **Better UX**: Optimized for drug ranking use case

### **The Solution**

Create a production-ready standalone `TherapyFitPage` that:
- ✅ Reuses existing `EfficacyPanel` component (no duplication)
- ✅ Works for any disease (universal)
- ✅ Explains S/P/E framework clearly
- ✅ Provides optimal UX for drug ranking

---

## 📋 PHASE 1: CORE PAGE IMPLEMENTATION (4-6 hours)

### **1.1 Create TherapyFitPage Component**

**File:** `oncology-coPilot/oncology-frontend/src/pages/TherapyFitPage.jsx`

**Requirements:**
- Standalone page component (not embedded)
- Reuses existing `EfficacyPanel` component
- Disease selector (universal)
- Mutation input (reuse `VariantInputList` or create universal)
- Model selector (reuse existing)
- S/P/E framework explanation section
- Provenance bar display

**Key Implementation Notes:**
- Follow MUI design patterns (consistent with existing pages)
- Use responsive layout (mobile/tablet/desktop)
- Handle loading and error states
- Pass disease context to `EfficacyPanel` if needed

**Success Criteria:**
- ✅ Page renders without errors
- ✅ Disease selector works
- ✅ Mutation input works
- ✅ EfficacyPanel displays results
- ✅ S/P/E explanation visible

---

### **1.2 Create S/P/E Framework Explanation Component**

**File:** `oncology-coPilot/oncology-frontend/src/components/therapy-fit/SPEFrameworkExplanation.jsx`

**Requirements:**
- Visual breakdown of S/P/E formula
- Explanation of each component (30% S, 40% P, 30% E)
- Evidence tiers explanation (Supported/Consider/Insufficient)
- Badges explanation (RCT, Guideline, ClinVar-Strong, PathwayAligned)
- Insights chips explanation (4 chips with confidence lifts)

**Design Notes:**
- Use icons for each component (Activity, Map, BookOpen)
- Color-coded sections (primary, secondary, success)
- Formula displayed in monospace font
- Chips for tiers and badges

**Success Criteria:**
- ✅ S/P/E breakdown clearly visible
- ✅ Formula displayed
- ✅ Evidence tiers explained
- ✅ Badges explained

---

### **1.3 Create/Enhance Universal Disease Selector**

**File:** `oncology-coPilot/oncology-frontend/src/components/common/DiseaseSelector.jsx`

**Requirements:**
- Dropdown/select for disease type
- Support all cancer types (myeloma, ovarian, melanoma, pancreatic, breast, colorectal, etc.)
- Maps to backend disease codes
- Clean, accessible UI

**Disease List:**
- Multiple Myeloma
- Ovarian Cancer (HGS)
- Melanoma
- Pancreatic Cancer
- Breast Cancer
- Colorectal Cancer
- Lung Cancer
- Prostate Cancer
- (Add more as needed)

**Success Criteria:**
- ✅ All diseases selectable
- ✅ Maps correctly to backend
- ✅ Accessible (keyboard navigation, screen readers)

---

### **1.4 Add Route to App.jsx**

**File:** `oncology-coPilot/oncology-frontend/src/App.jsx`

**Requirements:**
- Add import for `TherapyFitPage`
- Add route: `<Route path="/therapy-fit" element={<TherapyFitPage />} />`
- Consider adding to navigation menu (optional)

**Success Criteria:**
- ✅ Route accessible at `/therapy-fit`
- ✅ Page loads without errors
- ✅ Navigation works (if added to menu)

---

## 📋 PHASE 2: ENHANCEMENTS (2-3 hours)

### **2.1 Add Provenance Bar Component**

**File:** `oncology-coPilot/oncology-frontend/src/components/therapy-fit/ProvenanceBar.jsx`

**Requirements:**
- Display run_id, profile (SP/SPE), fusion status
- Show methods used
- Link to detailed provenance (optional)
- Match style of existing provenance displays

**Success Criteria:**
- ✅ Provenance visible after results load
- ✅ All key metadata displayed

---

### **2.2 Add Export Functionality**

**File:** `oncology-coPilot/oncology-frontend/src/components/therapy-fit/ExportButton.jsx`  
**File:** `oncology-coPilot/oncology-frontend/src/utils/exportUtils.js`

**Requirements:**
- Export to CSV (drug rankings with S/P/E breakdown)
- Include all metadata (run_id, profile, provenance)
- Optional: Export to PDF (can defer to Phase 3)

**CSV Format:**
- Columns: Rank, Drug Name, Efficacy Score, Confidence, Tier, S Score, P Score, E Score, Badges, Citations
- Include header row with metadata

**Success Criteria:**
- ✅ CSV export works
- ✅ Data formatted correctly
- ✅ All metadata included

---

### **2.3 Add Drug Comparison View** (Optional - Can Defer)

**File:** `oncology-coPilot/oncology-frontend/src/components/therapy-fit/DrugComparisonTable.jsx`

**Requirements:**
- Side-by-side comparison of top 3-5 drugs
- S/P/E breakdown for each
- Badges and citations visible
- Expandable rows for details

**Success Criteria:**
- ✅ Comparison table displays
- ✅ All key metrics visible

---

## 📋 PHASE 3: TESTING & VALIDATION (2-3 hours)

### **Test Cases:**

1. **BRAF V600E (Melanoma)**
   - Mutation: `BRAF V600E` (chr7:140453136 T>A)
   - Disease: `melanoma`
   - Expected: BRAF inhibitors rank #1
   - Validation: Check top drug is BRAF inhibitor, confidence >0.6

2. **KRAS G12D (Myeloma)**
   - Mutation: `KRAS G12D`
   - Disease: `multiple_myeloma`
   - Expected: MEK inhibitors rank higher than BRAF inhibitors
   - Validation: MEK confidence > BRAF confidence

3. **BRCA1 C61G (Ovarian)**
   - Mutation: `BRCA1 C61G` (chr17:43104911 T>G)
   - Disease: `ovarian_cancer`
   - Expected: Platinum agents rank #1
   - Validation: Top drug is platinum-based, tier is "supported"

4. **TP53 R248Q (Ovarian)**
   - Mutation: `TP53 R248Q` (chr17:7577120 G>A)
   - Disease: `ovarian_cancer_hgs`
   - Expected: DNA repair pathway drugs
   - Validation: DNA repair drugs in top 3

5. **Multiple Mutations**
   - Mutations: `BRCA1 C61G` + `TP53 R248Q`
   - Disease: `ovarian_cancer`
   - Expected: Combined scores, pathway aggregation
   - Validation: Scores reflect both mutations

### **Testing Checklist:**

- [ ] Page loads without errors
- [ ] Disease selector works for all diseases
- [ ] Mutation input accepts various formats
- [ ] API calls succeed
- [ ] Results display correctly
- [ ] S/P/E breakdown visible
- [ ] Evidence tiers correct
- [ ] Badges display correctly
- [ ] Provenance bar shows metadata
- [ ] Export works (if implemented)
- [ ] Responsive design (mobile/tablet/desktop)
- [ ] Error handling (invalid mutations, API errors)

---

## 🎯 SUCCESS CRITERIA

### **Functional Requirements:**
- ✅ Standalone `/therapy-fit` page accessible
- ✅ Works for any cancer type (universal)
- ✅ S/P/E framework explanation visible
- ✅ Drug rankings display correctly
- ✅ Evidence tiers work correctly
- ✅ Badges display correctly
- ✅ Provenance tracking works

### **UX Requirements:**
- ✅ Page loads in <2 seconds
- ✅ Results display in <5 seconds (API call)
- ✅ Responsive design (mobile/tablet/desktop)
- ✅ Accessible (WCAG 2.1 AA)

### **Business Requirements:**
- ✅ Page discoverable (route accessible)
- ✅ Clear value proposition
- ✅ Documentation accessible

---

## 📊 IMPLEMENTATION TIMELINE

| Phase | Hours | Tasks | Status |
|-------|-------|-------|--------|
| **Phase 1: Core Page** | 4-6h | TherapyFitPage, SPEFrameworkExplanation, DiseaseSelector, Route | ⏳ Pending |
| **Phase 2: Enhancements** | 2-3h | ProvenanceBar, ExportButton, DrugComparisonTable | ⏳ Pending |
| **Phase 3: Testing** | 2-3h | Test cases, validation, error handling | ⏳ Pending |
| **Total** | **9-14h** | **All phases** | **⏳ Pending** |

---

## 🏗️ ARCHITECTURE

### **Component Structure**

```
TherapyFitPage (New Standalone Page)
├── Header Section
│   ├── Title: "Therapy Fit: Drug Efficacy Ranking"
│   ├── Subtitle: "S/P/E Framework Explanation"
│   └── SPEFrameworkExplanation Component
├── Input Section
│   ├── DiseaseSelector
│   ├── VariantInputList (Reused)
│   └── ModelSelector (Reused)
├── Results Section
│   ├── EfficacyPanel (Reused - Main Component)
│   ├── ProvenanceBar (New - Phase 2)
│   └── DrugComparisonTable (New - Phase 2)
└── Actions Section
    └── ExportButton (New - Phase 2)
```

### **Data Flow**

```
User Input (Disease + Mutations)
  ↓
TherapyFitPage
  ↓
EfficacyPanel (Reused Component)
  ↓
useEfficacy Hook (Reused)
  ↓
POST /api/efficacy/predict
  ↓
EfficacyOrchestrator (Backend)
  ↓
Response (Drugs + S/P/E + Provenance)
  ↓
EfficacyPanel → Display Results
```

**Key Principle:** Reuse existing components, don't duplicate.

---

## 📝 FILE STRUCTURE

### **Files to Create**

```
oncology-coPilot/oncology-frontend/src/
├── pages/
│   └── TherapyFitPage.jsx                    # NEW - Main page
├── components/
│   ├── common/
│   │   └── DiseaseSelector.jsx               # NEW or ENHANCE
│   └── therapy-fit/                          # NEW - Directory
│       ├── SPEFrameworkExplanation.jsx       # NEW
│       ├── ProvenanceBar.jsx                 # NEW - Phase 2
│       ├── ExportButton.jsx                  # NEW - Phase 2
│       └── DrugComparisonTable.jsx           # NEW - Phase 2
└── utils/
    └── exportUtils.js                        # NEW - Phase 2
```

### **Files to Modify**

```
oncology-coPilot/oncology-frontend/src/
└── App.jsx                                    # MODIFY - Add route
```

### **Files to Reuse (No Changes)**

```
oncology-coPilot/oncology-frontend/src/
├── features/
│   └── efficacy/
│       ├── components/
│       │   ├── EfficacyPanel.jsx             # REUSE
│       │   └── EfficacyCard.jsx              # REUSE
│       ├── useEfficacy.js                    # REUSE
│       └── api.js                            # REUSE
└── components/
    └── myeloma/
        ├── VariantInputList.jsx              # REUSE
        └── ModelSelector.jsx                 # REUSE
```

---

## 🚀 DEPLOYMENT CHECKLIST

### **Pre-Deployment**
- [ ] All Phase 1 tasks complete
- [ ] All tests passing
- [ ] Code reviewed
- [ ] Documentation updated
- [ ] Route added to navigation (if applicable)

### **Deployment**
- [ ] Deploy to staging
- [ ] Test on staging
- [ ] Deploy to production
- [ ] Monitor for errors

### **Post-Deployment**
- [ ] Verify page accessible
- [ ] Test with real use cases
- [ ] Monitor analytics
- [ ] Gather user feedback

---

## 📚 REFERENCES

- **Backend Implementation**: `api/services/efficacy_orchestrator/orchestrator.py`
- **Frontend Component**: `src/features/efficacy/components/EfficacyPanel.jsx`
- **Concept Documentation**: `.cursor/rules/research/therapy_fit_concept.mdc`
- **Advanced Care Plan Pattern**: `.cursor/MOAT/ADVANCED_CARE_PLAN_UNIVERSAL.md`

---

## ✅ SUMMARY

**Current State:**
- ✅ Backend: 100% complete and tested
- ⚠️ Frontend: 80% complete (missing standalone page)

**What We Need:**
- ❌ Standalone `/therapy-fit` page
- ❌ Disease-agnostic interface
- ❌ S/P/E framework explanation

**Solution:**
- Create `TherapyFitPage` that reuses existing `EfficacyPanel`
- Add disease selector and mutation input
- Add S/P/E explanation section
- Enhance with visualization (Phase 2) and export (Phase 2)

**Effort:** 4-6 hours for Phase 1 (core page), 2-3 hours for Phase 2 (enhancements), 2-3 hours for Phase 3 (testing)

**Impact:** Makes Therapy Fit a first-class feature, improves discoverability, provides better UX for drug ranking use case.

**Status:** Ready for implementation  
**Priority:** High  
**Effort:** Low (reuses existing components)  
**Risk:** Low (well-tested backend, reusable frontend components)

---

**Document History:**
- January 2025: Initial comprehensive production-ready plan created
