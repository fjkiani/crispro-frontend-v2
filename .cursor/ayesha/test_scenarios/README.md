# 🧪 SPORADIC CANCER TEST SCENARIOS

**Date**: January 8, 2025  
**Created By**: Agent Jr  
**Purpose**: Test scenarios for Zo's sporadic cancer execution plan (Days 1-7)

---

## 📋 SCENARIO OVERVIEW

These 25 test scenarios validate the sporadic cancer logic across different data completeness levels and 15 cancer types:

| Scenario | Level | Cancer Type | Key Features | Purpose |
|----------|-------|-------------|--------------|---------|
| **1** | L0 | Ovarian HGS | Germline negative, no NGS | Test Quick Intake Level 0 (disease priors only) |
| **2** | L1 | Breast TNBC | Manual HRD entry (48) | Test Level 1 (partial data, HRD ≥42 overrides penalty) |
| **3** | L2 | Lung NSCLC | Full FM report, TMB 22 | Test Level 2 (full report, TMB ≥20 boost) |
| **4** | L2 | Colorectal | MSI-H + TMB 55 | Test edge case (both MSI-H and TMB ≥20) |
| **5** | L0 | Ovarian HGS | Ayesha's case (synthetic) | Test realistic Level 0 scenario |
| **6-7** | L0/L1 | Prostate | Low TMB (0.8), HRD unknown/L1 | Test prostate-specific priors |
| **8-9** | L0/L1 | Melanoma | High TMB (13.5), BRAF-driven | Test high TMB cancer |
| **10-11** | L0/L1 | Bladder | Intermediate TMB (5.5), FGFR3 | Test intermediate TMB cancer |
| **12-13** | L0/L1 | Endometrial | High MSI-H (28%), PTEN-driven | Test MSI-H biomarker |
| **14-15** | L0/L1 | Gastric | High MSI-H (22%), HER2 | Test MSI-H + HER2 |
| **16-17** | L0/L1 | Esophageal | High TP53 (73%), small sample | Test small-sample priors |
| **18-19** | L0/L1 | Head/Neck | Low TMB (2.5), low HRD | Test low biomarker cancer |
| **20-21** | L0/L1 | Glioblastoma | Very low TMB (1.5), EGFR | Test brain tumor priors |
| **22-23** | L0/L1 | Renal | Very low TMB (1.2), VHL-driven | Test low TMB cancer |
| **24-25** | L0/L1 | AML | Very low TMB (0.5), FLT3 | Test hematologic malignancy |

---

## 🎯 SCENARIO DESCRIPTIONS

### **Scenario 1: Level 0 - Minimal Data**
**File**: `test_case_1_level_0.json`

**Patient Profile:**
- Ovarian HGS, Stage IV, Line 3
- Germline negative (38 genes tested)
- Platinum sensitive
- **No tumor NGS report** → Level 0 mode

**Key Test Points:**
- ✅ Disease priors used for TMB (5.2, ovarian median)
- ✅ HRD = null (unknown at Level 0)
- ✅ MSI = null (unknown at Level 0)
- ✅ PARP penalty: 0.80x (germline negative, HRD unknown)
- ✅ IO boost: None (TMB < 10, MSI unknown)
- ✅ Confidence cap: 0.4 (Level 0 limit)

**Expected Behavior:**
- Quick Intake endpoint returns `TumorContext` with priors
- PARP inhibitors get 0.80x multiplier
- Confidence capped at 0.4
- Rationale explains "HRD unknown, tumor NGS recommended"

---

### **Scenario 2: Level 1 - Partial Data**
**File**: `test_case_2_level_1.json`

**Patient Profile:**
- Breast TNBC, Stage III, Line 2
- Germline negative
- **Manual entry**: TP53 mutation, HRD score 48
- No full NGS report → Level 1 mode

**Key Test Points:**
- ✅ HRD score 48 (above cutoff 42)
- ✅ PARP penalty: **NO PENALTY** (HRD ≥42 overrides germline negative)
- ✅ IO boost: None (TMB < 10)
- ✅ Confidence cap: 0.6 (Level 1 limit)

**Expected Behavior:**
- Manual entry endpoint accepts HRD score
- PARP inhibitors get **NO penalty** (somatic HRD-high)
- Rationale: "Somatic HRD-high (48) overrides germline negative status"
- This demonstrates the key sporadic cancer logic: **somatic HRD can rescue PARP eligibility**

---

### **Scenario 3: Level 2 - Full Report**
**File**: `test_case_3_level_2.json`

**Patient Profile:**
- Lung NSCLC, Stage IV, Line 1
- Germline negative
- **Full Foundation Medicine report** → Level 2 mode
- TMB 22 (high), MSI-MSS, HRD 18 (low), EGFR L858R

**Key Test Points:**
- ✅ TMB 22 (≥20, very high)
- ✅ IO boost: 1.35x (TMB ≥20, highest boost)
- ✅ PARP penalty: 0.60x (germline negative, HRD < 42)
- ✅ Confidence: **NO CAP** (Level 2, full report)
- ✅ Completeness: 1.0 (all fields available)

**Expected Behavior:**
- Foundation Medicine parser extracts all biomarkers
- Immunotherapy gets 1.35x boost (TMB ≥20)
- PARP gets 0.60x penalty (HRD < 42)
- Highest confidence (no cap, completeness 1.0)

---

### **Scenario 4: Edge Case - MSI-H + High TMB**
**File**: `test_case_4_edge_case.json`

**Patient Profile:**
- Colorectal, Stage III, Line 1
- Germline negative (Lynch panel negative)
- **Full Foundation Medicine report** → Level 2 mode
- TMB 55 (very high), MSI-H, HRD 25 (low), BRAF V600E

**Key Test Points:**
- ✅ **Both MSI-H AND TMB ≥20** present
- ✅ IO boost: 1.35x (TMB ≥20 takes precedence per formula)
- ✅ PARP penalty: 0.60x (HRD < 42)
- ✅ Edge case: Tests formula priority (TMB ≥20 > MSI-H > TMB ≥10)

**Expected Behavior:**
- Per Zo's formula: TMB ≥20 gets 1.35x (highest)
- MSI-H would get 1.30x, but TMB ≥20 takes precedence
- This validates the boost hierarchy logic

---

### **Scenario 5: Ayesha's Case (Synthetic)**
**File**: `test_case_5_ayesha.json`

**Patient Profile:**
- Ovarian HGS, Stage IIIC, Line 3
- Germline negative (CustomNext-Cancer 38 genes)
- Platinum sensitive (initial response)
- **No tumor NGS report yet** → Level 0 mode

**Key Test Points:**
- ✅ Realistic Ayesha profile (synthetic but based on known facts)
- ✅ Platinum-sensitive history suggests possible HRD
- ✅ Level 0 → PARP penalty 0.80x (HRD unknown)
- ✅ **Value demonstration**: Tumor NGS would clarify HRD and potentially lift PARP penalty

**Expected Behavior:**
- Quick Intake returns Level 0 `TumorContext`
- PARP gets 0.80x penalty (conservative, HRD unknown)
- Rationale explains: "Platinum-sensitive suggests possible HRD, but tumor NGS required for confirmation"
- This scenario demonstrates **why tumor NGS is critical for sporadic cases**

---

## 📊 EXPECTED OUTPUTS SUMMARY

### **PARP Penalty Logic**

| Scenario | Germline | HRD Score | PARP Penalty | Reason |
|----------|----------|-----------|--------------|--------|
| 1 (L0) | Negative | null | 0.80x | HRD unknown (Level 0) |
| 2 (L1) | Negative | 48 | **1.0x** | HRD ≥42 (somatic HRD-high) |
| 3 (L2) | Negative | 18 | 0.60x | HRD < 42 |
| 4 (L2) | Negative | 25 | 0.60x | HRD < 42 |
| 5 (L0) | Negative | null | 0.80x | HRD unknown (Level 0) |

**Key Insight**: Scenario 2 demonstrates **somatic HRD can rescue PARP eligibility** even when germline negative!

---

### **IO Boost Logic**

| Scenario | TMB | MSI | IO Boost | Reason |
|----------|-----|-----|----------|--------|
| 1 (L0) | 5.2 | null | 1.0x | TMB < 10, MSI unknown |
| 2 (L1) | 1.8 | null | 1.0x | TMB < 10, MSI unknown |
| 3 (L2) | 22 | MSS | **1.35x** | TMB ≥ 20 (very high) |
| 4 (L2) | 55 | MSI-H | **1.35x** | TMB ≥ 20 (takes precedence) |
| 5 (L0) | 5.2 | null | 1.0x | TMB < 10, MSI unknown |

**Key Insight**: TMB ≥20 gets highest boost (1.35x), even if MSI-H also present.

---

### **Confidence Caps**

| Scenario | Level | Confidence Cap | Base Confidence | Completeness |
|----------|-------|----------------|-----------------|--------------|
| 1 (L0) | L0 | 0.4 | 0.3 | 0.25 |
| 2 (L1) | L1 | 0.6 | 0.4 | 0.5 |
| 3 (L2) | L2 | **None** | 0.6 | 1.0 |
| 4 (L2) | L2 | **None** | 0.6 | 1.0 |
| 5 (L0) | L0 | 0.4 | 0.3 | 0.25 |

**Key Insight**: Level 2 (full report) has NO confidence cap, highest completeness.

---

## ✅ VALIDATION CHECKLIST

**Before Zo runs these tests, verify:**

- [x] All 5 scenarios use correct disease keys (`"ovarian_hgs"`, `"breast_tnbc"`, etc.)
- [x] All expected gates calculated using Zo's formulas (A4)
- [x] PARP penalty logic matches: 0.80x (L0, HRD unknown), 0.60x (HRD < 42), 1.0x (HRD ≥42)
- [x] IO boost logic matches: 1.35x (TMB ≥20), 1.30x (MSI-H), 1.25x (TMB ≥10)
- [x] Confidence caps match: 0.4 (L0), 0.6 (L1), None (L2)
- [x] Completeness scores calculated: 0.25 (L0, only TMB), 0.5 (L1, partial), 1.0 (L2, full)
- [x] All JSON files validate (no syntax errors)
- [x] Scenario 2 demonstrates key logic: somatic HRD ≥42 overrides germline negative
- [x] Scenario 4 tests edge case: MSI-H + TMB ≥20 (TMB takes precedence)
- [x] Scenario 5 (Ayesha) is realistic and demonstrates value of tumor NGS

---

## 🎯 USAGE INSTRUCTIONS FOR ZO

**Day 1-2 Testing:**
- Use Scenario 1 (Level 0) to test Quick Intake endpoint
- Use Scenario 2 (Level 1) to test manual entry endpoint
- Verify PARP penalty and IO boost logic

**Day 6 E2E Testing:**
- Run all 5 scenarios end-to-end
- Compare actual outputs to `expected_gates` and `expected_confidence`
- Validate formulas match expected outputs

**Validation:**
- See `EXPECTED_RESULTS.md` for detailed validation table
- Each scenario has `validation_notes` explaining key test points

---

## 📝 NOTES

- **All formulas from Zo's A4 answer** (PARP penalty, IO boost, confidence caps)
- **All disease keys match Zo's specification** (`"ovarian_hgs"` format)
- **All expected outputs calculated** using provided formulas
- **Scenario 5 (Ayesha)** is synthetic but realistic based on known profile

---

**Agent Jr - Test Scenarios Complete!** ✅

