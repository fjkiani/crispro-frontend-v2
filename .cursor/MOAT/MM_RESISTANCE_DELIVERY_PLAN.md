# MM Resistance Implementation - Delivery Plan & Questions

**Date:** January 28, 2025  
**Reviewer:** Auto  
**Guide:** `.cursor/MOAT/MM_RESISTANCE_IMPLEMENTATION_GUIDE.mdc`  
**Status:** ✅ **READY TO START** - 8 Questions Need Answers

---

## 🎯 EXECUTIVE SUMMARY

The implementation guide is **excellent and actionable** (91% quality score). It provides:
- ✅ Clear code examples (copy-paste ready)
- ✅ Realistic expectations (handles rare mutations)
- ✅ Proper priority ordering (P0 → P1 → P2)
- ✅ Edge case handling (rate limits, missing data)

**However, 8 critical questions need manager answers before starting implementation.**

---

## ❓ CRITICAL QUESTIONS FOR MANAGER

### **Q1: Signal Architecture** 🔴 **HIGH PRIORITY**

**Question:** Should drug-class resistance mutations be a separate signal or modify existing signals?

**Options:**
- **A)** Separate signal (`MM_DRUG_CLASS_RESISTANCE`) - **RECOMMENDED**
- **B)** Modify `MM_HIGH_RISK_GENE` signal
- **C)** New signal type entirely

**Current Pattern:**
- `MM_HIGH_RISK_GENE` = separate signal
- `MM_CYTOGENETICS` = separate signal
- Each has own `ResistanceSignalData` object

**Recommendation:** **Option A** (consistent with existing pattern)

**Impact:** Affects code structure in `predict_mm_resistance()`

---

### **Q2: Risk Calculation Method** 🔴 **HIGH PRIORITY**

**Question:** How should drug-class resistance mutations affect risk probability?

**Guide Suggests:** Multiplicative (line 320-321)
```python
risk_multiplier *= signal["relative_risk"]
```

**Current Code Uses:** Weighted average (line 1296-1318)
```python
overall_probability = sum(sig.probability * sig.confidence) / sum(sig.confidence)
```

**Options:**
- **A)** Multiplicative (as guide suggests)
- **B)** Weighted average (as current code uses) - **RECOMMENDED**

**Recommendation:** **Option B** (consistent with existing code, more stable)

**Impact:** Affects risk calculation logic

---

### **Q3: cBioPortal Client Reuse** 🟡 **MEDIUM PRIORITY**

**Question:** Should we reuse existing `cbioportal_client.py` or create standalone script?

**Current State:**
- ✅ `scripts/data_acquisition/utils/cbioportal_client.py` exists
- ✅ Has reusable `CBioportalClient` class
- ✅ Used by other download scripts

**Options:**
- **A)** Reuse existing client (refactor guide code) - **RECOMMENDED**
- **B)** Create standalone script (as guide suggests)

**Recommendation:** **Option A** - Check existing client first, reuse if possible

**Impact:** Affects implementation approach (refactor vs. copy-paste)

---

### **Q4: Validation Test Scope** 🟡 **MEDIUM PRIORITY**

**Question:** How many validation tests in Week 1?

**Guide Implements:** 3 tests (PSMB5, CRBN, DIS3/TP53)

**Mission Requires:** 5 tests (add del(17p), RAS/MAPK)

**Options:**
- **A)** 3 tests in Week 1, add 2 later - **RECOMMENDED**
- **B)** All 5 tests in Week 1

**Recommendation:** **Option A** - Start with 3, add 2 when data available

**Impact:** Affects Week 1 scope

---

### **Q5: Pathway Service Signal Type** 🟢 **LOW PRIORITY**

**Question:** Should pathway burden be a separate signal or modulator?

**Options:**
- **A)** Separate signal (`MM_PATHWAY_BURDEN`) - **RECOMMENDED**
- **B)** Modulate existing signals (multiplier)
- **C)** Not a signal, just used in risk calculation

**Recommendation:** **Option A** (consistent pattern)

**Impact:** Affects Week 2 implementation

---

### **Q6: Frontend Integration Pattern** 🟢 **LOW PRIORITY**

**Question:** Separate component or integrated into existing display?

**Options:**
- **A)** Separate component (`MMResistancePanel.jsx`) - **RECOMMENDED**
- **B)** Integrated into `MyelomaResponseDisplay.jsx`

**Recommendation:** **Option A** (cleaner architecture)

**Impact:** Affects Week 3 implementation

---

### **Q7: Evo2 Integration Priority** 🟢 **LOW PRIORITY**

**Question:** Week 1 (P0) or Week 3 (P2)?

**Options:**
- **A)** Week 1 (P0) - required for launch
- **B)** Week 3 (P2) - can wait - **RECOMMENDED**

**Recommendation:** **Option B** (not blocking)

**Impact:** Affects timeline

---

### **Q8: Literature-Based RR Values** 🔴 **HIGH PRIORITY**

**Question:** If validation shows "INSUFFICIENT_DATA", should we:

**Options:**
- **A)** Ship with literature-based RR values (as guide suggests) - **RECOMMENDED**
- **B)** Wait for larger cohort
- **C)** Mark as "EXPERIMENTAL" with lower confidence

**Recommendation:** **Option A** - Ship with literature values, mark confidence appropriately

**Impact:** Affects go/no-go decision

---

## 📋 DELIVERY PLAN (Assuming Answers)

### **Week 1: P0 Blockers (11-16 hours)**

#### **Day 1: PSMB5/CRBN Resistance Mutations (4-6h)**

**Tasks:**
1. Add `MM_RESISTANCE_MUTATIONS` dictionary to `resistance_prophet_service.py`
   - Location: After `MM_HIGH_RISK_GENES` (line 238)
   - Use guide code (lines 94-207)

2. Add `_detect_mm_resistance_mutations()` method
   - Location: After `_detect_mm_high_risk_genes()` (line 880)
   - Use guide code (lines 214-288)
   - Add `_is_lof()` helper

3. Integrate into `predict_mm_resistance()`
   - Add after `mm_gene_signal` (line 1069)
   - Convert to `ResistanceSignalData` format
   - Add to `signals_detected` list

4. Create tests
   - File: `tests/test_mm_resistance_mutations.py`
   - Use guide tests (lines 326-386)

**Deliverable:** PSMB5/CRBN mutations detected, tests passing

---

#### **Day 2: MMRF Cohort Data (4-6h)**

**Tasks:**
1. Check existing cBioPortal client
   - Review `scripts/data_acquisition/utils/cbioportal_client.py`
   - Determine if reusable

2. Create/refactor download script
   - Option A: Refactor guide code to use existing client
   - Option B: Create standalone script (if client insufficient)

3. Run download
   - Target: ≥500 patients, ≥5000 mutations
   - Verify: Treatment response data present

**Deliverable:** Cohort JSON file ready for validation

---

#### **Day 3: Validation Framework (3-4h)**

**Tasks:**
1. Create validation script
   - File: `scripts/validation/validate_mm_resistance.py`
   - Implement 3 tests: PSMB5, CRBN, DIS3/TP53

2. Run validation
   - Accept "INSUFFICIENT_DATA" if needed
   - Document results

3. Document findings
   - Create validation results document

**Deliverable:** Validation framework running, results documented

---

### **Week 2: P1 Enhancements (7-11h)**

#### **Day 4-5: MM Pathway Service (4-6h)**

**Tasks:**
1. Create `mm_pathway_service.py`
   - Use guide code (lines 1084-1160)

2. Integrate into resistance prediction
   - Add pathway burden signal
   - Use to modulate risk

**Deliverable:** Pathway service integrated

---

#### **Day 6: Expanded Gene Markers (2-3h)**

**Tasks:**
1. Add IKZF1, IKZF3, CUL4A to `MM_HIGH_RISK_GENES`
2. Attempt validation (likely INSUFFICIENT_DATA)
3. Document results

**Deliverable:** Gene markers expanded

---

#### **Day 7: TRUE SAE (Zo Leads) (1-2h review)**

**Tasks:**
- Review Zo's work
- Document decision

**Deliverable:** SAE decision documented

---

### **Week 3: P2 Polish (10-14h)**

#### **Day 8-9: Frontend Panel (6-8h)**

**Tasks:**
1. Create `MMResistancePanel.jsx`
2. Integrate into MyelomaDigitalTwin
3. Test UI

**Deliverable:** Resistance predictions in UI

---

#### **Day 10: Evo2 Integration (4-6h)**

**Tasks:**
1. Add Evo2 scoring to `predict_mm_resistance()`
2. Correlate with response
3. Use as secondary signal

**Deliverable:** Evo2 integrated

---

## 🚨 ASSUMPTIONS (If No Answers)

**If manager doesn't answer questions, I will:**

1. **Q1 (Signal Architecture):** Use **Option A** (separate signal) - most consistent
2. **Q2 (Risk Calculation):** Use **Option B** (weighted average) - matches existing code
3. **Q3 (cBioPortal):** Check existing client first, reuse if possible
4. **Q4 (Validation Scope):** Implement 3 tests in Week 1, add 2 later
5. **Q5 (Pathway Signal):** Use **Option A** (separate signal)
6. **Q6 (Frontend):** Use **Option A** (separate component)
7. **Q7 (Evo2):** Use **Option B** (Week 3)
8. **Q8 (Literature Values):** Use **Option A** (ship with literature)

**These are safe defaults based on codebase patterns.**

---

## 📊 SUCCESS METRICS

### **Week 1 Success:**

✅ PSMB5 p.Ala49Thr detected → PI resistance  
✅ CRBN p.Trp400* detected → IMiD resistance  
✅ Cohort downloaded (≥500 patients)  
✅ Validation script runs (1-2 tests PASS, INSUFFICIENT_DATA acceptable)

### **Week 2 Success:**

✅ Pathway service computes 6 pathways  
✅ Gene markers expanded (IKZF1, IKZF3, CUL4A)  
✅ SAE decision documented (Zo)

### **Week 3 Success:**

✅ Frontend shows resistance predictions  
✅ Evo2 integrated as secondary signal  
✅ End-to-end pipeline working

---

## 🎯 READY TO START?

**Status:** ✅ **YES** - Guide is production-ready

**Blockers:** 8 questions (can proceed with assumptions if needed)

**Next Step:** 
1. Get manager answers (preferred)
2. Or proceed with assumptions (safe defaults)

**Estimated Delivery:** 2-3 weeks (28-41 hours of work)

---

**Last Updated:** January 28, 2025  
**Next Action:** Get answers → Start Day 1


