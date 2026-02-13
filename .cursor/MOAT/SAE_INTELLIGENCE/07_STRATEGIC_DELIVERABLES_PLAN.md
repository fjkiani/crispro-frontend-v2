# Strategic Deliverables Plan: Next 3, 5, 10 Deliverables

**Date:** January 29, 2026 (Audit Update)  
**Status:** ✅ **PIVOTED** - True SAE Deprecated. MFAP4 & Proxy Logic Prioritized.  
**Location:** `.cursor/MOAT/SAE_INTELLIGENCE/07_STRATEGIC_DELIVERABLES_PLAN.md`

---

## 🎯 EXECUTIVE SUMMARY

**Strategic Focus (Jan 2026):**
The "Diamond SAE" track (0.783 AUROC) has been **shut down** following the audit.
Resources are re-allocated to:
1.  **MFAP4 Deployment**: Visualizing the Gold-standard biomarker (AUROC 0.763).
2.  **Proxy Scaling**: Expanding the "Gene -> Pathway" logic to new cancers.

---

## 🎯 NEXT 3 DELIVERABLES (IMMEDIATE)

### **Deliverable 1: MFAP4 Frontend Integration** 🔴 **HIGH PRIORITY**

**Status:** 🟡 In Progress  
**Goal:** Display "EMT Risk Score" based on MFAP4 expression.  
**Why:** It is our *only* verified predictive biomarker for platinum resistance.

**What to Build:**
*   **Backend**: `get_transcriptomic_risk(patient_id)` function.
*   **UI**: `EMTRiskGauge.jsx` (Red/Green zone based on Z-score > 1.5).
*   **Logic**: High MFAP4 = High Risk.

### **Deliverable 2: Proxy Steerability UI** 🟡 **MEDIUM PRIORITY**

**Status:** 🟡 In Progress  
**Goal:** "What-If" sliders for Genes (not Features).  
**Why:** Clinicians trust genes. They don't trust "Latent Feature 27607".

**What to Build:**
*   **UI**: `GeneToggle.jsx` (e.g., "Toggle NF1 status").
*   **Logic**: Re-compute Mechanism Vector -> Re-rank Trials.

### **Deliverable 3: Manuscript Submission** ✅ **COMPLETE**

**Status:** ✅ Ready  
**Goal:** Submit MFAP4 paper to bioRxiv.  
**Action:** Detailed in `MANUSCRIPT_DRAFT.md`.

---

## 🛑 DEPRECATED DELIVERABLES (DO NOT BUILD)

### ~~Deliverable 0: TRUE SAE Diamonds Mapping~~
*   **Status**: ❌ **KILLED**
*   **Reason**: Validation Failure (AUROC 0.555).

### ~~Deliverable 1.5: TRUE SAE Frontend Integration~~
*   **Status**: ❌ **KILLED**
*   **Reason**: Feature 27607 is noise. Do not display it.

---

## 🎯 LONG TERM (3+ MONTHS)

### **Deliverable 4: Multi-Cancer Proxy Expansion**
*   **Goal**: Lung (KRAS/EGFR) and Breast (BRCA/PI3K).
*   **Method**: Replicate the "Proxy Logic" (Mutations -> Pathways) for these types.

### **Deliverable 5: Longitudinal Kinetics**
*   **Goal**: CA-125 Velocity tracking.
*   **Method**: Time-series analysis of EMR data.
