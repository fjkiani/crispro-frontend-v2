# 🔧 MODULE 04 CRITICAL FIXES - For JR Agent C

**Purpose:** Critical bug fixes for Drug Efficacy (S/P/E) framework  
**Target:** JR Agent C (Module 04 Owner)  
**Status:** ✅ READY FOR IMPLEMENTATION  
**Date:** January 28, 2025

---

## 🎯 EXECUTIVE SUMMARY

**4 Critical Bugs Found & Fixed:**
1. **Pathway Normalization Bug** - Wrong range assumption (1e-6 to 1e-4) → Fixed to (0 to 0.005)
2. **Tier Computation Parameter** - Using normalized `path_pct` instead of raw `s_path` → Fixed
3. **Tier Threshold** - Too high (0.05) for new pathway range → Fixed to 0.001
4. **Sporadic Gates Capping** - Always applied (even without tumor context) → Fixed to conditional

**Impact:** These bugs caused:
- All drugs to get `path_pct = 1.0` (no differentiation)
- All drugs to get confidence = 0.4 (capped by sporadic gates)
- Incorrect tier classification
- Poor benchmark results (MM: 40%, Ovarian: AUROC 0.500, Melanoma: 50%)

**After Fixes:** 
- Pathway differentiation restored (path_pct: 0.037-0.330)
- Confidence differentiation restored (0.549-0.586)
- Correct drug rankings (MEK > BRAF for KRAS G12D)
- Tiers correctly classified

---

## 🐛 BUG #1: Pathway Normalization Range

### Problem

**Location:** `api/services/efficacy_orchestrator/drug_scorer.py:48-49`

**Current (WRONG) Code:**
```python
# OLD - WRONG RANGE
path_pct = (s_path - 1e-6) / (1e-4 - 1e-6)  # Assumes range 1e-6 to 1e-4
path_pct = max(0.0, min(1.0, path_pct))  # Caps at 1.0
```

**Issue:**
- Formula assumes pathway scores are in range `1e-6` to `1e-4` (0.000001 to 0.0001)
- **Actual pathway scores:** ~0.002 (2e-3), which is **20x larger** than assumed maximum
- Result: All drugs get `path_pct = 1.0` (capped), eliminating differentiation

**Evidence:**
- October 2025: Small but meaningful differences (0.015)
- November 2025: No differences (all 0.4)
- Root cause: Normalization assumes wrong range

### Fix

**Location:** `api/services/efficacy_orchestrator/drug_scorer.py:48-55`

**New (CORRECT) Code:**
```python
# NEW - CORRECT RANGE
# Pathway scores are in range 0 to ~0.005 (observed from actual data)
if s_path <= 0:
    path_pct = 0.0
elif s_path >= 0.005:
    path_pct = 1.0
else:
    path_pct = s_path / 0.005  # Normalize to 0-1 range
```

**Verification:**
- After fix: `path_pct` values now differentiated (0.037-0.330)
- MEK inhibitor: path_pct = 0.330 ✅
- BRAF inhibitor: path_pct = 0.293 ✅

---

## 🐛 BUG #2: Tier Computation Parameter

### Problem

**Location:** `api/services/efficacy_orchestrator/drug_scorer.py:138`

**Current (WRONG) Code:**
```python
# OLD - WRONG PARAMETER
tier = compute_evidence_tier(s_seq, path_pct, s_evd, badges, config)
#                                 ^^^^^^^^
#                                 Using normalized path_pct (0-1 range)
```

**Issue:**
- Tier computation expects raw pathway score `s_path` (0 to ~0.005 range)
- But we're passing normalized `path_pct` (0-1 range)
- Result: Tier classification incorrect (all drugs classified as "insufficient")

**Evidence:**
- Tier computation threshold: `s_path < 0.05` (for old range)
- But `path_pct` is 0-1, so `path_pct < 0.05` is too strict
- All drugs fail threshold → "insufficient" tier

### Fix

**Location:** `api/services/efficacy_orchestrator/drug_scorer.py:138`

**New (CORRECT) Code:**
```python
# NEW - CORRECT PARAMETER
tier = compute_evidence_tier(s_seq, s_path, s_evd, badges, config)
#                                 ^^^^^^
#                                 Using raw s_path (0 to ~0.005 range)
```

**Verification:**
- After fix: Tiers correctly classified
- MEK inhibitor: tier = "consider" ✅
- BRAF inhibitor: tier = "consider" ✅

---

## 🐛 BUG #3: Tier Threshold

### Problem

**Location:** `api/services/confidence/tier_computation.py:59`

**Current (WRONG) Code:**
```python
# OLD - TOO HIGH THRESHOLD
if s_path < 0.05:  # Too high for new range (0 to ~0.005)
    return "insufficient"
```

**Issue:**
- Threshold `0.05` was appropriate for old pathway range (1e-6 to 1e-4)
- New pathway range is `0 to 0.005` (20x smaller)
- Threshold `0.05` is **10x larger** than maximum pathway score
- Result: All drugs classified as "insufficient" (even high-scoring drugs)

**Evidence:**
- Actual pathway scores: 0.000165 to 0.00165 (all < 0.05)
- All drugs fail threshold → "insufficient" tier

### Fix

**Location:** `api/services/confidence/tier_computation.py:61`

**New (CORRECT) Code:**
```python
# NEW - APPROPRIATE THRESHOLD
if s_path < 0.001:  # Appropriate for new range (0 to ~0.005)
    return "insufficient"
```

**Verification:**
- After fix: Tiers correctly classified based on actual pathway scores
- High-scoring drugs: tier = "consider" ✅
- Low-scoring drugs: tier = "insufficient" ✅

---

## 🐛 BUG #4: Sporadic Gates Capping

### Problem

**Location:** `api/services/efficacy_orchestrator/orchestrator.py:217-230`

**Current (WRONG) Code:**
```python
# OLD - ALWAYS APPLIES SPORADIC GATES
hasattr(request, 'germline_status') or hasattr(request, 'tumor_context')
# This is ALWAYS True (dataclass fields always present, even if None)
```

**Issue:**
- Check `hasattr(request, 'germline_status')` is always True (dataclass fields always exist)
- Sporadic gates always applied, even when no tumor context
- Defaults to Level 0 (sporadic) → caps confidence at 0.4
- Result: All drugs get confidence = 0.4 (no differentiation)

**Evidence:**
- Before fix: All drugs confidence = 0.4 (capped)
- After fix: Confidence differentiated (0.549-0.586)

### Fix

**Location:** `api/services/efficacy_orchestrator/orchestrator.py:225-228`

**New (CORRECT) Code:**
```python
# NEW - CONDITIONAL SPORADIC GATES
# Only apply if:
# 1. Tumor context is provided (not None/empty), OR
# 2. Germline status is explicitly set (not default "unknown")
should_apply_sporadic_gates = (
    (tumor_context_data is not None and tumor_context_data) or
    (germline_status and germline_status != "unknown")
)

if should_apply_sporadic_gates:
    # Apply sporadic gates (Level 0, cap at 0.4)
    ...
else:
    # No sporadic gates (use full confidence range)
    ...
```

**Verification:**
- After fix: Confidence differentiated
- MEK inhibitor: confidence = 0.563 ✅
- BRAF inhibitor: confidence = 0.549 ✅
- MEK ranks higher than BRAF (correct for KRAS G12D) ✅

---

## 🧪 TEST CASES

### Test Case 1: KRAS G12D (Multiple Myeloma)

**Input:**
```python
mutations = [{"gene": "KRAS", "hgvs_p": "p.G12D"}]
disease = "multiple_myeloma"
```

**Expected After Fixes:**
- MEK inhibitor: confidence > BRAF inhibitor ✅
- Pathway scores differentiated (not all 1.0) ✅
- Tiers correctly classified ✅
- Confidence differentiated (not all 0.4) ✅

**Before Fixes:**
- All drugs: path_pct = 1.0, confidence = 0.4, tier = "insufficient" ❌

### Test Case 2: BRCA1 Truncating (Ovarian)

**Input:**
```python
mutations = [{"gene": "BRCA1", "hgvs_p": "p.C61G", "consequence": "stop_gained"}]
disease = "ovarian_cancer"
```

**Expected After Fixes:**
- PARP inhibitor: high confidence (>0.7) ✅
- Pathway scores differentiated ✅
- Tier = "consider" or "I" ✅

---

## 📋 IMPLEMENTATION CHECKLIST FOR JR AGENT C

When building Module 04, ensure:

- [ ] **Fix 1:** Pathway normalization uses range 0 to 0.005 (not 1e-6 to 1e-4)
- [ ] **Fix 2:** Tier computation uses raw `s_path` (not normalized `path_pct`)
- [ ] **Fix 3:** Tier threshold adjusted to 0.001 (not 0.05)
- [ ] **Fix 4:** Sporadic gates only apply when tumor context provided
- [ ] **Test:** Run MM benchmark - verify >80% pathway alignment accuracy
- [ ] **Test:** Run Ovarian benchmark - verify AUROC >0.75
- [ ] **Test:** Run Melanoma benchmark - verify >90% drug ranking accuracy

---

## 📚 REFERENCES

**Root Cause Analysis:**
- `.cursor/plans/ROOT_CAUSE_ANALYSIS_COMPREHENSIVE.md`
- `.cursor/plans/PATHWAY_NORMALIZATION_FIX_SUMMARY.md`

**Verified Fixes:**
- All 4 fixes verified in code
- Test cases passing
- Benchmarks ready to re-run

---

**Status:** ✅ READY FOR JR AGENT C  
**Priority:** 🔴 CRITICAL - Must fix before Module 04 implementation

