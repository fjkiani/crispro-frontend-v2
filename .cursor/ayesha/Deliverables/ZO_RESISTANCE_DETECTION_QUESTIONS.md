# ⚔️ RESISTANCE DETECTION 2-OF-3 RULE - QUESTIONS FOR MANAGER

**Author**: Zo (Lead AI Agent)  
**Date**: January 14, 2025  
**Purpose**: Questions about clinical rationale and validation of 2-of-3 trigger rule

---

## 🎯 CONTEXT

**Current State**: Resistance detection uses 2-of-3 trigger rule (Manager's Policy C7)

**Implementation**: `api/services/resistance_detection_service.py:82-266`

**Triggers**:
1. HRD drop ≥15 points
2. DNA repair capacity drop ≥0.20
3. CA-125 inadequate response (on-therapy rise OR <50% drop by cycle 3)

**Logic**: Alert if **2 of 3** conditions met

**Gap**: I understand the implementation and thresholds, but not the clinical rationale for 2-of-3 vs 1-of-3 vs 3-of-3.

---

## 📋 CURRENT IMPLEMENTATION

### **2-of-3 Trigger Rule** (`resistance_detection_service.py:189-193`)

```python
trigger_count = len(triggers_met)
resistance_detected = trigger_count >= 2
```

**Thresholds**:
- HRD drop: ≥15 points (e.g., 58 → 43)
- DNA repair drop: ≥0.20 (e.g., 0.75 → 0.50)
- CA-125 inadequate: On-therapy rise OR <50% drop by cycle 3

**HR Restoration Pattern** (R2):
- Special case: HRD drop + DNA repair drop (coherent signal)
- Immediate alert (don't wait for radiology)
- Recommended: Switch to ATR/CHK1 inhibitors

---

## ❓ QUESTIONS FOR MANAGER

### **Question 1: Clinical Rationale for 2-of-3**

**What I See**: 
- Logic: Alert if 2 of 3 triggers met
- Not 1-of-3 (too sensitive?) or 3-of-3 (too specific?)

**My Question**:
1. **What is the clinical rationale for 2-of-3 vs 1-of-3 vs 3-of-3?**
   - Why not 1-of-3? (Would catch more cases, but more false positives?)
   - Why not 3-of-3? (Would be more specific, but miss early signals?)
2. **What sensitivity/specificity trade-off does 2-of-3 achieve?**
   - What's the target sensitivity? (≥80%?)
   - What's the acceptable false positive rate? (<20%?)
3. **Is 2-of-3 based on:**
   - Clinical experience (oncologist consensus)?
   - Literature evidence (resistance detection studies)?
   - Empirical validation (retrospective cohort analysis)?
   - First principles (signal fusion theory)?

**Why This Matters**:
- Need to justify rule in clinical documentation
- Need to know if rule is evidence-based or heuristic
- Need to understand if rule can be validated/refined

---

### **Question 2: Validation Status**

**What I See**: 
- Rule is implemented (Manager's Policy C7)
- No validation mentioned in code/docs

**My Question**:
1. **Has the 2-of-3 rule been validated against real patient outcomes?**
   - Do patients with 2-of-3 triggers actually develop resistance?
   - Do patients without 2-of-3 triggers remain sensitive?
2. **If validated, what was the methodology?**
   - Retrospective cohort analysis?
   - Clinical trial data?
   - Literature meta-analysis?
3. **If validated, what were the results?**
   - Sensitivity? (≥75%?)
   - Specificity? (≥80%?)
   - Lead time? (3-6 weeks earlier than imaging?)
4. **If not validated, should we validate it?**
   - What would validation look like?
   - What metrics would we use? (sensitivity, specificity, lead time, PPV, NPV?)

**Why This Matters**:
- Clinical credibility requires validation
- Need to know if rule is evidence-based or placeholder
- Need to plan validation studies if missing

---

### **Question 3: Threshold Selection**

**What I See**: 
- HRD drop: ≥15 points
- DNA repair drop: ≥0.20
- CA-125: <50% drop by cycle 3

**My Question**:
1. **How were these thresholds determined?**
   - Clinical experience?
   - Literature evidence?
   - Data-driven (ROC curve analysis)?
   - First principles (biological significance)?
2. **Are thresholds disease-specific?**
   - Do ovarian cancer thresholds differ from breast cancer?
   - Should we adjust thresholds based on disease type?
3. **Are thresholds treatment-specific?**
   - Do PARP inhibitor thresholds differ from chemotherapy?
   - Should we adjust thresholds based on treatment?

**Why This Matters**:
- Need to justify thresholds in clinical documentation
- Need to know if thresholds are universal or disease-specific
- Need to understand if thresholds can be refined

---

### **Question 4: False Positive Management**

**What I See**: 
- 2-of-3 rule triggers alert
- Immediate action recommended (order tests, switch therapy)

**My Question**:
1. **What's the acceptable false positive rate?**
   - How many false alarms are acceptable?
   - What's the cost of false positive? (unnecessary tests, therapy switches)
2. **How do we manage false positives?**
   - Confirm with additional tests before switching?
   - Wait for imaging confirmation?
   - Use confidence scores to gate alerts?
3. **What's the false positive rate we're seeing?**
   - Have we tracked false positives?
   - What's the actual rate? (<20%?)

**Why This Matters**:
- Need to balance sensitivity vs specificity
- Need to avoid unnecessary interventions
- Need to manage clinician trust (too many false alarms = ignored alerts)

---

### **Question 5: Disease-Specific Adjustments**

**What I See**: 
- Rule is implemented generically (not disease-specific)

**My Question**:
1. **Should thresholds be disease-specific?**
   - Ovarian cancer: Different thresholds than breast cancer?
   - Different CA-125 thresholds for different cancers?
2. **Should trigger combinations be disease-specific?**
   - Ovarian: HRD + CA-125 more important?
   - Breast: DNA repair + imaging more important?
3. **Should we have disease-specific rules?**
   - Ovarian-specific resistance detection?
   - Breast-specific resistance detection?

**Why This Matters**:
- Ayesha has ovarian cancer
- Need to know if generic rule applies or needs adjustment
- Need to know if we should create disease-specific rules

---

### **Question 6: Clinical Validation Plan**

**What I See**: 
- Rule is implemented and used
- Validation status unclear

**My Question**:
1. **Should we validate the 2-of-3 rule prospectively?**
   - Track patients over time?
   - Compare 2-of-3 triggers to actual resistance?
   - Measure lead time vs imaging?
2. **What would validation look like?**
   - Retrospective cohort analysis?
   - Prospective observational study?
   - Clinical trial integration?
3. **What metrics would we track?**
   - Sensitivity (≥75%?)
   - Specificity (≥80%?)
   - Lead time (3-6 weeks?)
   - PPV (≥70%?)
   - NPV (≥85%?)

**Why This Matters**:
- Need to prove clinical value
- Need to refine rule based on outcomes
- Need to publish validation results

---

## ✅ PROPOSED ANSWERS FOR AGENT

### Clinical rationale (2‑of‑3 vs 1‑of‑3 vs 3‑of‑3)
- Single signals (e.g., CA‑125 alone) are noisy → 1‑of‑3 over‑alerts.
- 3‑of‑3 is too strict and alerts late (often after radiographic progression).
- 2‑of‑3 fuses orthogonal modalities (biomarker trend + genomic capacity change) to balance sensitivity and specificity while preserving lead time.

### Target operating characteristics (policy)
- Sensitivity: 75–85%
- Specificity: 70–85%
- False‑positive rate: <20%
- Lead time vs imaging: 3–6 weeks (median)

### Validation status and plan
- If not completed: run longitudinal retrospective cohort.
  - Derive trigger timestamps; compare to radiographic progression.
  - Metrics: lead time distribution, sensitivity/specificity, PPV/NPV, KM separation post‑alert.

### Threshold selection (policy basis)
- HRD drop ≥15: exceeds test–retest noise; clinically meaningful.
- DNA repair capacity drop ≥0.20 (0–1): large effect size consistent with pathway restoration/escape.
- CA‑125 inadequate (on‑therapy rise or <50% drop by cycle 3): aligns with ovarian monitoring guidance.
- Disease‑specific tuning: CA‑125 criterion is ovarian‑specific; substitute relevant markers (e.g., CEA/CA19‑9) elsewhere or lean more on genomic/ctDNA dynamics.

### False‑positive management (tiered actions)
- Watch alert: first 2‑of‑3 → repeat labs (1–2 weeks), order confirmatory test (HRD/ctDNA).
- Action alert: persistent 2‑of‑3 or HR restoration pattern (concordant HRD↓ + DNA repair↓ on PARP) → discuss switch/ATR‑ or CHK1‑trial; include rationale and required confirmations.

### Disease‑specific adjustments
- Ovarian: keep CA‑125 rule; thresholds above apply.
- Other cancers: replace CA‑125 with disease‑relevant markers or elevate weight of genomic trajectory; re‑validate thresholds per disease.

### Prospective validation (recommended)
- Track alerts and clinician responses prospectively; measure outcome impact and refine thresholds under change control.

---

## 📊 CURRENT STATE DOCUMENTATION

**File**: `RESISTANCE_DETECTION_RATIONALE.md` (to be created)

**Will Document**:
- Current 2-of-3 logic
- Thresholds and their sources
- Clinical use cases
- Integration points (resistance playbook, resistance prophet)

---

## 🎯 EXPECTED OUTCOMES

**After Manager Answers**:
1. ✅ Understand clinical rationale for 2-of-3 vs 1-of-3 vs 3-of-3
2. ✅ Know if rule is validated or needs validation
3. ✅ Understand threshold selection process
4. ✅ Know acceptable false positive rate
5. ✅ Know if rule is disease-specific or universal
6. ✅ Plan validation study if needed

---

**Status**: ✅ **ANSWERS PROVIDED**  
**Last Updated**: January 14, 2025  
**By**: Zo (Lead AI Agent)  
**Answers Added**: January 14, 2025

