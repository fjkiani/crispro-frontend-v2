# SAE Resistance Pipeline Roadmap (2026)

**Strategy**: "Biology First"
**Focus**: Deployment of MFAP4 & Proxy Logic.

---

## 🟢 PHASE 1: DEPLOYMENT (Current)

### **1. MFAP4 Injection**
*   **Goal**: Integrate the AUROC 0.763 biomarker into the clinical UI.
*   **Status**: ✅ Code Complete (`detect_transcriptomic_risk`).
*   **Action**: Ensure UI displays "EMT Risk" prominently.

### **2. Proxy Logic Scaling**
*   **Goal**: Expand "Proxy SAE" (Mutation -> Pathway) map to Lung and Breast.
*   **Status**: 🚧 In Progress.
*   **Action**: Map KRAS/EGFR for Lung; PIK3CA/BRCA for Breast.

### **3. Diamond Purge**
*   **Goal**: Remove all "True SAE" code.
*   **Status**: ✅ Complete.
*   **Action**: Monitor logs to ensure no feature-27607 lookups occur.

---

## 🟡 PHASE 2: EXPANSION (Q2 2026)

### **1. Integration with Timing Engine**
*   **Concept**: Combine `Risk Score` (Resistance) with `Time-to-Event` (Prognosis).
*   **Goal**: "High Risk + Short PFI = High Urgency".

### **2. Multi-Omic Dashboard**
*   **Concept**: Visualizing Layer 1 (RNA) + Layer 2 (DNA) + Layer 3 (Time).
*   **Goal**: A single "Mission Control" view for the clinician.

---

## 🔴 PHASE 3: RESEARCH (Future)

### **1. True SAE Resurrection?**
*   **Condition**: Only revisited if N > 500 patients available with full WES+RNA.
*   **Lesson**: Never trust small-N feature extraction again.

---
