# 🎯 AYESHA'S NEXT 10 HIGH-IMPACT DELIVERABLES

**Date**: January 29, 2025  
**Status**: 📋 **PRIORITIZED FOR IMMEDIATE IMPLEMENTATION**  
**Patient**: AK (Ayesha) - Stage IVB HGSOC  
**Profile Status**: L1 Completeness (0.55) - Has IHC + germline, missing NGS/CA-125

---

## 🚨 EXECUTIVE SUMMARY

**Current State**: Ayesha has 12 capabilities operational, but **critical gaps** remain.

**Her Unique Profile:**
- ✅ **Germline**: MBD4 homozygous pathogenic (rare DDR pathway mutation)
- ✅ **IHC**: TP53 mutant (IHC evidence), PD-L1 positive (CPS=10), ER weakly positive
- ⚠️ **Missing**: CA-125 value, NGS somatic mutations, HRD score, TMB
- ⚠️ **Clinical Context**: Progression from baseline (2/2024) to widespread metastases (10/2025)

**Impact Priority**: These 10 deliverables will transform her experience from "capabilities exist" to "personalized, actionable care plan."

---

## 📊 HER CURRENT PROFILE ANALYSIS

### **What She Has (L1 Completeness)**
- ✅ Germline testing: MBD4 homozygous pathogenic, PDGFRA VUS
- ✅ IHC biomarkers: TP53 mutant, PD-L1 positive, ER weakly positive, MMR preserved
- ✅ Imaging: Baseline (2/2024) + Progression CT (10/2025) + PET (11/2025)
- ✅ Pathology: Multiple biopsies confirming Mullerian origin

### **What She's Missing**
- ❌ CA-125 value (null) - Needed for monitoring
- ❌ NGS somatic mutations - Only IHC evidence of TP53
- ❌ HRD score (null) - Critical for PARPi eligibility
- ❌ TMB (null) - Needed for IO decisions
- ❌ Full genomic coordinates - Only gene-level evidence

### **What She Needs (High Impact)**
1. **DDR Status** - MBD4 is DDR pathway, needs classification
2. **CA-125 Entry** - Can't monitor without value
3. **Next Test Recommendations** - Needs NGS, HRD, ctDNA
4. **MBD4 Education** - Rare mutation, needs patient-friendly explanation
5. **Progression Tracking** - Visual timeline of disease progression
6. **PARP Eligibility** - Based on DDR status
7. **Patient "What Changes Now" Report** - Personalized action plan
8. **Clinical Dossier** - Comprehensive report export
9. **Treatment Timeline** - Visual treatment history
10. **Synthetic Lethality Display** - MBD4+TP53 combination analysis

---

## 🎯 PRIORITY 1: DDR_BIN STATUS CALCULATION & DISPLAY (HIGHEST IMPACT)

### **Why Critical for Ayesha**
- She has **MBD4 homozygous pathogenic** mutation (DDR pathway gene)
- DDR_bin will classify her as `DDR_defective` → **PARPi eligible**
- **Clinical Impact**: Determines if she can receive PARPi maintenance therapy

### **Current State**
- ❌ Backend engine ready but **API endpoint not created**
- ❌ Frontend has **no DDR status components**
- ❌ Not displayed in care plan

### **Deliverable**
**Create DDR_bin status calculation and display integrated with her profile:**

**Backend:**
1. Create `POST /api/resistance/ddr-status` endpoint
2. Auto-calculate from her profile:
   - `disease_site: "ovary"`
   - `mutations: [{gene: "MBD4", variant_classification: "pathogenic"}]`
   - `tumor_subtype: "HGSOC"`

**Frontend:**
1. Create `DDRStatusCard` component - Show primary status (DDR_defective/DDR_proficient/unknown)
2. Auto-populate form from `patientProfile.germline.mutations`
3. Display in `GenomicFoundationSection` of care plan
4. Show PARPi eligibility badge

**Data Source:**
```javascript
// From ayesha_11_17_25.js
{
  disease_site: "ovary",
  tumor_subtype: "HGSOC",
  mutations: [
    {
      gene_symbol: "MBD4",
      variant_classification: "pathogenic", // From germline.mutations[0]
      variant_type: "indel" // c.1293delA
    }
  ],
  // HRD score: null (will be unknown until NGS)
}
```

**Expected Result:**
- `DDR_bin_status: "DDR_defective"` (MBD4 pathogenic = extended DDR)
- `BRCA_pathogenic: false`
- `extended_DDR_pathogenic: true`
- **PARPi Eligibility**: ✅ **ELIGIBLE** (DDR_defective)

**Files to Create/Modify:**
- `src/components/ddr/DDRStatusCard.jsx` (new)
- `src/components/complete-care/sections/GenomicFoundationSection.jsx` (add DDR status)
- `src/hooks/useDDRStatus.js` (new)
- `api/routers/resistance.py` (add endpoint)

**Estimated Effort**: 12-16 hours  
**Patient Impact**: **HIGH** - Determines PARPi eligibility

---

## 🎯 PRIORITY 2: CA-125 VALUE ENTRY & MONITORING (CRITICAL DATA)

### **Why Critical for Ayesha**
- Her profile has `ca125_value: null` - **Cannot monitor burden/response**
- CA-125 intelligence requires a value
- **Current gap**: She can't see CA-125 tracking without entering value

### **Current State**
- ⚠️ CA-125 intelligence backend exists
- ❌ **No UI for entering/updating CA-125 value**
- ❌ Profile has `ca125_value: null` - blocks monitoring

### **Deliverable**
**Create CA-125 value entry component and integrate with profile:**

**Component:**
1. `CA125ValueEntry.jsx` - Input form for CA-125 value + date
2. Update `patientProfile.labs.ca125_value` via PatientContext
3. Auto-trigger CA-125 intelligence after value entered
4. Show in patient profile summary card

**Integration:**
- Add to patient profile page
- Add quick entry button to care plan page
- Show "CA-125: N/A" → "CA-125: [Enter Value]" prompt

**Data Flow:**
```
User enters CA-125 value
    ↓
Update patientProfile.labs.ca125_value
    ↓
Auto-trigger CA-125 intelligence calculation
    ↓
Display CA125Tracker with burden/forecast/resistance
```

**Files to Create/Modify:**
- `src/components/patient/CA125ValueEntry.jsx` (new)
- `src/pages/PatientProfile.jsx` (add entry form)
- `src/components/complete-care/CompleteCareHeader.jsx` (add quick entry)
- `src/context/PatientContext.jsx` (add updateCA125 function)

**Estimated Effort**: 4-6 hours  
**Patient Impact**: **HIGH** - Unlocks CA-125 monitoring capability

---

## 🎯 PRIORITY 3: NEXT TEST RECOMMENDATIONS DISPLAY (ACTIONABLE GUIDANCE)

### **Why Critical for Ayesha**
- She's **L1 completeness** (missing NGS, HRD, CA-125)
- Needs clear guidance on what tests to order next
- Backend has recommendations but **not prominently displayed**

### **Current State**
- ✅ Backend `next_test_recommender` exists
- ⚠️ Displayed in care plan but **not prominent**
- ❌ **No actionable "Order Test" buttons** or lab integration

### **Deliverable**
**Create prominent Next Test Recommendations component with actionable buttons:**

**Component:**
1. `NextTestRecommendationsCard.jsx` - Large, prominent card at top of care plan
2. Show prioritized tests:
   - **Priority 1**: NGS (FoundationOne, Tempus, Guardant360)
   - **Priority 2**: HRD assay (Myriad, Leuven)
   - **Priority 3**: ctDNA (Guardant360, FoundationOne Liquid)
   - **Priority 4**: SLFN11 IHC (if available)
3. For each test: Show **rationale** (why needed) + **"Order Test" button** (links to lab portal or generates order form)

**Display Logic:**
- If `hrd_score === null` → Show "HRD Assay Recommended" with rationale
- If `tumor_context.somatic_mutations.length < 2` → Show "NGS Recommended" with rationale
- If `ca125_value === null` → Show "CA-125 Lab Order" with rationale

**Ayesha's Specific Recommendations:**
1. ✅ **HRD Assay** - Priority 1 (needed for PARPi eligibility, current: null)
2. ✅ **NGS (FoundationOne/Tempus)** - Priority 1 (only IHC evidence, current: TP53 IHC only)
3. ✅ **CA-125 Lab** - Priority 2 (monitoring, current: null)
4. ✅ **ctDNA (Guardant360)** - Priority 3 (baseline for monitoring)

**Files to Create/Modify:**
- `src/components/patient/NextTestRecommendationsCard.jsx` (new)
- `src/components/complete-care/sections/GenomicFoundationSection.jsx` (add at top)
- `src/pages/AyeshaCompleteCare.jsx` (display prominently)

**Estimated Effort**: 6-8 hours  
**Patient Impact**: **HIGH** - Provides actionable next steps

---

## 🎯 PRIORITY 4: MBD4 PATIENT EDUCATION SECTION (RARE MUTATION)

### **Why Critical for Ayesha**
- She has **MBD4 homozygous pathogenic** - extremely rare mutation
- **MBD4-associated neoplasia syndrome (MANS)** - needs explanation
- Patients need to understand: "What does this mean for me?"

### **Current State**
- ✅ Backend has mutation explanation capabilities
- ❌ **No patient-friendly MBD4 education component**
- ❌ **No MANS syndrome explanation** in UI

### **Deliverable**
**Create patient-friendly MBD4 education component:**

**Component:**
1. `MBD4PatientEducationCard.jsx` - Patient-friendly explanation
2. Content sections:
   - **"What is MBD4?"** - Simple explanation (DNA repair gene)
   - **"What does your mutation mean?"** - MBD4 homozygous = both copies broken
   - **"MANS Syndrome"** - What it means for her cancer risk (AML, CRC)
   - **"What this means for your treatment"** - DDR pathway disruption → PARPi may help
   - **"Questions to ask your doctor"** - Patient empowerment

**Data Source:**
```javascript
// From ayesha_11_17_25.js
{
  gene: "MBD4",
  variant: "c.1293delA",
  protein_change: "p.K431Nfs*54",
  zygosity: "homozygous",
  classification: "pathogenic",
  syndrome: "MBD4-associated neoplasia syndrome (MANS)",
  risk_increases: ["Acute myelogenous leukemia (AML)", "Colorectal cancer (CRC)"]
}
```

**Display:**
- Show in `GenomicFoundationSection` or dedicated "Your Genetics" section
- Collapsible/expandable for detailed info
- Link to clinical guidelines or research papers

**Files to Create/Modify:**
- `src/components/patient/MBD4PatientEducationCard.jsx` (new)
- `src/components/complete-care/sections/GenomicFoundationSection.jsx` (add education card)

**Estimated Effort**: 4-6 hours  
**Patient Impact**: **MEDIUM-HIGH** - Patient empowerment, understanding

---

## 🎯 PRIORITY 5: DISEASE PROGRESSION TIMELINE VISUALIZATION (CLINICAL CONTEXT)

### **Why Critical for Ayesha**
- She **progressed from baseline (2/2024) to widespread mets (10/2025)**
- Timeline shows: baseline → progression → diagnosis → treatment
- **Visual timeline** helps her understand disease course

### **Current State**
- ✅ Profile has `diagnostic_timeline` array
- ❌ **No visual timeline component** - just text list
- ❌ **No progression indicators** or visual cues

### **Deliverable**
**Create visual disease progression timeline component:**

**Component:**
1. `DiseaseProgressionTimeline.jsx` - Interactive timeline visualization
2. Show key events:
   - **2/1/2024**: Baseline CT - Ovarian cysts, NO carcinomatosis
   - **10/28/2025**: Progression CT - **CARCINOMATOSIS DETECTED** ⚠️
   - **11/11/2025**: PET Scan - **WIDESPREAD METASTASES** 🔴
   - **11/17/2025**: Pathology - Diagnosis confirmed
   - **11/20/2025**: Surgical biopsies - Tissue confirmed
   - **11/24/2025**: Germline testing - MBD4 detected

**Visual Elements:**
- Timeline bar with dates
- Color-coded events (baseline=green, progression=red, diagnosis=orange)
- Hover tooltips with details
- Zoom in/out for different time scales

**Data Source:**
```javascript
// From ayesha_11_17_25.js
diagnostic_timeline: [
  { date: "2024-02-01", report_type: "CT_SCAN", finding: "Baseline...", stage: "baseline" },
  { date: "2025-10-28", report_type: "CT_SCAN", finding: "PROGRESSION DETECTED...", stage: "metastatic_detected" },
  // ... etc
]
```

**Display:**
- Show in patient dashboard or dedicated "Disease Timeline" page
- Link to imaging reports and pathology reports

**Files to Create/Modify:**
- `src/components/patient/DiseaseProgressionTimeline.jsx` (new)
- `src/pages/PatientDashboard.jsx` (add timeline section)

**Estimated Effort**: 6-8 hours  
**Patient Impact**: **MEDIUM** - Clinical context, understanding disease course

---

## 🎯 PRIORITY 6: PARP ELIGIBILITY BADGE & RECOMMENDATIONS (TREATMENT DECISION)

### **Why Critical for Ayesha**
- DDR_bin status determines **PARPi eligibility**
- She likely **DDR_defective** (MBD4 pathogenic) → **PARPi eligible**
- Needs clear **"You are eligible for PARPi"** messaging

### **Current State**
- ⚠️ DDR_bin can determine eligibility but **not displayed**
- ❌ **No PARPi eligibility badge** in care plan
- ❌ **No PARPi-specific recommendations** panel

### **Deliverable**
**Create PARP eligibility badge and recommendations component:**

**Component:**
1. `PARPEligibilityBadge.jsx` - Large, prominent eligibility indicator
2. Display logic:
   - If `DDR_bin_status === "DDR_defective"` → ✅ **"PARP INHIBITOR ELIGIBLE"** (green badge)
   - If `DDR_bin_status === "DDR_proficient"` → ❌ **"PARP INHIBITOR NOT RECOMMENDED"** (red badge)
   - If `DDR_bin_status === "unknown"` → ⚠️ **"INSUFFICIENT DATA - ADDITIONAL TESTING NEEDED"** (yellow badge)

3. `PARPRecommendationsPanel.jsx` - Detailed recommendations:
   - **If eligible**: Show PARPi options (olaparib, niraparib, rucaparib)
   - **Timing**: Frontline maintenance vs recurrent
   - **Evidence**: Link to SOLO-1, PRIMA, NOVA trials
   - **Monitoring**: What to watch for (resistance, side effects)

**Integration:**
- Show in `TherapeuticIntelligenceSection` after DDR status calculated
- Link to relevant clinical trials
- Link to resistance monitoring (Resistance Prophet)

**Files to Create/Modify:**
- `src/components/treatment/PARPEligibilityBadge.jsx` (new)
- `src/components/treatment/PARPRecommendationsPanel.jsx` (new)
- `src/components/complete-care/sections/TherapeuticIntelligenceSection.jsx` (add PARP section)

**Estimated Effort**: 4-6 hours  
**Patient Impact**: **HIGH** - Treatment decision support

---

## 🎯 PRIORITY 7: PATIENT "WHAT CHANGES NOW" PERSONALIZED REPORT (ACTION PLAN)

### **Why Critical for Ayesha**
- She needs **clear, personalized action plan**
- Not just "here are capabilities" but **"here's what YOU should do next"**
- Combines all recommendations into single personalized report

### **Current State**
- ✅ All capabilities exist but **scattered across different sections**
- ❌ **No unified "What Changes Now" report**
- ❌ **No patient-specific action plan**

### **Deliverable**
**Create personalized "What Changes Now" report generator:**

**Component:**
1. `WhatChangesNowReport.jsx` - Personalized action plan
2. Sections:
   - **"Your DDR Status"** - DDR_defective → PARPi eligible
   - **"Tests You Need"** - HRD, NGS, CA-125 (prioritized)
   - **"Treatment Options"** - SOC + PARPi maintenance eligibility
   - **"What to Ask Your Doctor"** - Personalized questions
   - **"Monitoring Plan"** - CA-125, imaging schedule

**Generation Logic:**
- Read from `patientProfile` and `carePlan.result`
- Generate personalized recommendations based on:
  - Her DDR status
  - Her missing tests (NGS, HRD, CA-125)
  - Her mutation (MBD4) → specific guidance
  - Her progression status → monitoring frequency

**Display:**
- Show as first section in care plan or dedicated page
- Exportable as PDF (for sharing with doctor)
- Update automatically when profile/care plan changes

**Files to Create/Modify:**
- `src/components/patient/WhatChangesNowReport.jsx` (new)
- `src/pages/AyeshaCompleteCare.jsx` (add as first section)
- `src/utils/patientReportGenerator.js` (new - report generation logic)

**Estimated Effort**: 8-10 hours  
**Patient Impact**: **HIGH** - Clear, actionable patient guidance

---

## 🎯 PRIORITY 8: CLINICAL DOSSIER EXPORT (SHARING WITH DOCTOR)

### **Why Critical for Ayesha**
- She needs to **share comprehensive report with her oncologist**
- Clinical dossier = all data + recommendations in single PDF/document
- **Export capability** for doctor visits

### **Current State**
- ⚠️ Export logic exists but **not user-friendly**
- ❌ **No "Export Clinical Dossier" button** prominently displayed
- ❌ **No PDF generation** with formatting

### **Deliverable**
**Create clinical dossier export functionality:**

**Component:**
1. `ClinicalDossierExporter.js` - PDF generation utility
2. Content sections:
   - Patient demographics + profile summary
   - Disease timeline
   - All imaging reports (CT, PET)
   - All pathology reports
   - Germline mutations (MBD4, PDGFRA VUS)
   - IHC biomarkers summary
   - DDR status + PARPi eligibility
   - Next test recommendations
   - Treatment recommendations (SOC + PARPi)
   - Clinical trials (top 10)
   - Resistance playbook
   - SAE features

3. Export button in care plan header
4. Format: PDF with table of contents, page numbers, doctor-friendly formatting

**Integration:**
- Add to `CompleteCareActions` component
- Generate from `patientProfile` + `carePlan.result`
- Download as `Ayesha_Clinical_Dossier_[Date].pdf`

**Files to Create/Modify:**
- `src/utils/clinicalDossierExporter.js` (new - PDF generation)
- `src/components/complete-care/CompleteCareActions.jsx` (add export button)
- `src/hooks/useClinicalDossier.js` (new - export logic)

**Estimated Effort**: 6-8 hours  
**Patient Impact**: **MEDIUM** - Doctor communication, documentation

---

## 🎯 PRIORITY 9: SYNTHETIC LETHALITY DISPLAY (MBD4+TP53 COMBINATION)

### **Why Critical for Ayesha**
- She has **MBD4 homozygous + TP53 mutant** combination
- This creates **synthetic lethality** opportunity (BER loss + checkpoint bypass)
- **Backend has analysis** but **not prominently displayed**

### **Current State**
- ✅ Backend synthetic lethality service exists (`/api/synthetic-lethality/analyze`)
- ⚠️ Component exists but **not in care plan** (orphaned)
- ❌ **Not prominently displayed** for her combination

### **Deliverable**
**Create prominent synthetic lethality display for MBD4+TP53:**

**Component:**
1. `SyntheticLethalityCard.jsx` - Enhanced display (already exists, needs integration)
2. Show for MBD4+TP53:
   - **"Synthetic Lethality Detected"** - MBD4 (BER loss) + TP53 (checkpoint bypass)
   - **Mechanism**: "BER pathway disrupted + DNA damage checkpoint bypassed"
   - **Therapeutic Opportunity**: DDR-targeting drugs (PARPi, ATR inhibitors)
   - **Evidence**: Link to research on MBD4+TP53 combinations

**Display:**
- Show in `PathwayMechanismSection` of care plan
- Prominent badge: "SYNTHETIC LETHALITY DETECTED"
- Link to mechanism map visualization

**Data Source:**
```javascript
// From synthetic lethality analysis
{
  gene_pair: ["MBD4", "TP53"],
  mechanism: "BER_loss + checkpoint_bypass",
  therapeutic_opportunity: "DDR_targeting",
  confidence: 0.85
}
```

**Files to Create/Modify:**
- `src/components/complete-care/sections/PathwayMechanismSection.jsx` (ensure SL card displayed)
- Verify `SyntheticLethalityCard.jsx` receives MBD4+TP53 data

**Estimated Effort**: 2-4 hours (integration, not new component)  
**Patient Impact**: **MEDIUM** - Treatment opportunity identification

---

## 🎯 PRIORITY 10: COMPLETENESS SCORE DISPLAY & IMPROVEMENT ROADMAP (DATA GAPS)

### **Why Critical for Ayesha**
- She has **L1 completeness (0.55)** - missing NGS/CA-125
- Needs **clear visualization** of what's missing
- Needs **roadmap** to improve completeness → better recommendations

### **Current State**
- ✅ Profile has `completeness_score: 0.55` (L1)
- ❌ **No visual display** of completeness score
- ❌ **No "how to improve" roadmap**

### **Deliverable**
**Create completeness score display and improvement roadmap:**

**Component:**
1. `CompletenessScoreCard.jsx` - Visual score display
   - Progress bar: 55% complete (L1)
   - Breakdown: "You have: IHC + Germline. Missing: NGS + CA-125"
   - Target: L2 completeness (0.75) requires NGS, L3 (1.0) requires all

2. `CompletenessRoadmap.jsx` - Improvement steps
   - **Step 1**: Order NGS (FoundationOne/Tempus) → +0.20 → L2 (0.75)
   - **Step 2**: Order HRD assay → Unlocks HRD status
   - **Step 3**: Enter CA-125 value → +0.05 → L2+ (0.80)
   - **Step 4**: Order ctDNA baseline → Monitoring capability

**Display:**
- Show in patient profile summary or care plan header
- Link to Next Test Recommendations (Priority 3 above)
- Show "Your recommendations will improve with more data" messaging

**Files to Create/Modify:**
- `src/components/patient/CompletenessScoreCard.jsx` (new)
- `src/components/patient/CompletenessRoadmap.jsx` (new)
- `src/components/complete-care/CompleteCareHeader.jsx` (add score display)

**Estimated Effort**: 4-6 hours  
**Patient Impact**: **MEDIUM** - Data gap awareness, improvement guidance

---

## 📋 IMPLEMENTATION ROADMAP (Priority Order)

### **Week 1: Critical Treatment Decisions**
- [ ] **Priority 1**: DDR_BIN Status Calculation (12-16 hours)
- [ ] **Priority 2**: CA-125 Value Entry (4-6 hours)
- [ ] **Priority 6**: PARP Eligibility Badge (4-6 hours)

**Total Week 1**: 20-28 hours  
**Patient Impact**: Unlocks PARPi eligibility determination

### **Week 2: Actionable Guidance**
- [ ] **Priority 3**: Next Test Recommendations (6-8 hours)
- [ ] **Priority 7**: "What Changes Now" Report (8-10 hours)
- [ ] **Priority 4**: MBD4 Patient Education (4-6 hours)

**Total Week 2**: 18-24 hours  
**Patient Impact**: Clear action plan for patient

### **Week 3: Clinical Context & Polish**
- [ ] **Priority 5**: Disease Progression Timeline (6-8 hours)
- [ ] **Priority 8**: Clinical Dossier Export (6-8 hours)
- [ ] **Priority 9**: Synthetic Lethality Display (2-4 hours)
- [ ] **Priority 10**: Completeness Score Display (4-6 hours)

**Total Week 3**: 18-26 hours  
**Patient Impact**: Clinical context, documentation, data awareness

---

## 🎯 SUCCESS METRICS

### **Patient Experience**
- ✅ DDR status calculated automatically from profile
- ✅ CA-125 value can be entered/updated easily
- ✅ Clear "what to do next" action plan
- ✅ PARPi eligibility prominently displayed
- ✅ Patient-friendly MBD4 education available

### **Clinical Impact**
- ✅ PARPi eligibility determined (DDR_defective)
- ✅ Next test recommendations prioritized
- ✅ Treatment options clear (SOC + PARPi)
- ✅ Resistance monitoring plan established

### **Data Completeness**
- ✅ Completeness score visible
- ✅ Roadmap to improve completeness
- ✅ Clear guidance on missing tests

---

**Last Updated**: January 29, 2025  
**Status**: 📋 **READY FOR IMPLEMENTATION**  
**Total Estimated Effort**: 56-78 hours (3 weeks)  
**Expected Patient Impact**: **TRANSFORMATIVE** - From "capabilities exist" to "personalized, actionable care plan"
