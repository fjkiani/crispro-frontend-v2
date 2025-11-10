# 🐛 SPORADIC GATES VALIDATION - BUG REPORT

**Date**: January 8, 2025 (Evening)  
**Reporter**: Zo (Mission 3 Validation)  
**Module**: `api/services/efficacy_orchestrator/sporadic_gates.py`  
**Test Scenarios**: Agent Jr (Mission 2)

---

## 📊 SUMMARY

**Total Failures**: 5 IO boost failures out of 25 scenarios  
**Failure Rate**: 20% (5/25 scenarios)  
**Critical**: No - All PARP and Confidence tests pass (50/50)

---

## 🐛 BUG #1: TMB BOOST FACTOR MISMATCH (2 FAILURES)

### **Issue**
Test scenarios expect 1.35x boost for TMB ≥20, but `sporadic_gates.py` uses 1.3x.

### **Affected Scenarios**
- **Test Case 3**: Lung NSCLC, TMB 22, Expected 0.81, Got 0.78
- **Test Case 4**: Colorectal, TMB 55, Expected 0.81, Got 0.78

### **Root Cause**
**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sporadic_gates.py`  
**Line**: 153**

```python
if tmb is not None and tmb >= 20:
    tmb_boost = 1.3  # ❌ Should be 1.35 per Zo's A4 answer
```

**Expected (from Zo's A4 answer)**: TMB ≥20 → 1.35x boost  
**Actual (implementation)**: TMB ≥20 → 1.3x boost

### **Impact**
- 2 scenarios fail IO boost validation
- Efficacy scores are 0.03 lower than expected (0.78 vs 0.81)
- This is a **moderate impact** - scores are close but not exact

### **Proposed Fix**

**Option A (Recommended)**: Update implementation to match test scenarios
```python
if tmb is not None and tmb >= 20:
    tmb_boost = 1.35  # ✅ Match Zo's A4 answer
    rationale.append({
        "gate": "IO_TMB_BOOST",
        "verdict": "BOOSTED",
        "boost": 1.35,  # ✅ Updated
        "tmb": tmb,
        "reason": f"TMB-high (≥20): {tmb:.1f} mut/Mb → Checkpoint inhibitor boost 1.35x"
    })
```

**Option B**: Update test scenarios to match implementation (1.3x)
- **Not recommended** - Zo's A4 answer specifies 1.35x

### **Validation**
After fix, re-run:
- Test Case 3: Expected 0.81, Should get 0.81 ✅
- Test Case 4: Expected 0.81, Should get 0.81 ✅

---

## 🐛 BUG #2: MISSING TMB ≥10 TIER (1 FAILURE)

### **Issue**
Test scenarios expect 1.25x boost for TMB ≥10 but <20, but `sporadic_gates.py` doesn't have this tier.

### **Affected Scenarios**
- **Test Case 9**: Melanoma, TMB 13.5, Expected 0.75 (0.60 * 1.25), Got 0.60 (no boost)

### **Root Cause**
**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sporadic_gates.py`  
**Lines**: 150-160

Current logic only checks:
- TMB ≥20 → 1.3x boost
- TMB <20 → No boost

**Missing**: TMB ≥10 but <20 → 1.25x boost (per Zo's A4 answer)

### **Impact**
- 1 scenario fails IO boost validation
- Efficacy score is 0.15 lower than expected (0.60 vs 0.75)
- This is a **moderate impact** - significant score difference

### **Proposed Fix**

Add TMB ≥10 tier to `sporadic_gates.py`:

```python
# Check TMB
tmb = tumor_context.get("tmb")
if tmb is not None:
    if tmb >= 20:
        tmb_boost = 1.35  # ✅ Updated from 1.3
        rationale.append({
            "gate": "IO_TMB_BOOST",
            "verdict": "BOOSTED",
            "boost": 1.35,
            "tmb": tmb,
            "reason": f"TMB-high (≥20): {tmb:.1f} mut/Mb → Checkpoint inhibitor boost 1.35x"
        })
    elif tmb >= 10:  # ✅ NEW TIER
        tmb_boost = 1.25
        rationale.append({
            "gate": "IO_TMB_BOOST",
            "verdict": "BOOSTED",
            "boost": 1.25,
            "tmb": tmb,
            "reason": f"TMB-intermediate (≥10): {tmb:.1f} mut/Mb → Checkpoint inhibitor boost 1.25x"
        })
```

### **Validation**
After fix, re-run:
- Test Case 9: Expected 0.75, Should get 0.75 ✅

---

## 🐛 BUG #3: MSI STATUS STRING MISMATCH (2 FAILURES)

### **Issue**
Test scenarios use `"MSI-H"` but `sporadic_gates.py` checks for `"MSI-High"` (exact match).

### **Affected Scenarios**
- **Test Case 13**: Endometrial, MSI-H, Expected 0.78, Got 0.60
- **Test Case 15**: Gastric, MSI-H, Expected 0.78, Got 0.60

### **Root Cause**
**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sporadic_gates.py`  
**Line**: 164

```python
if msi_status == "MSI-High":  # ❌ Exact match only, doesn't accept "MSI-H"
```

**Test Scenarios Use**: `"MSI-H"`  
**Implementation Expects**: `"MSI-High"`  
**Mismatch**: String comparison fails, no boost applied

### **Impact**
- 2 scenarios fail IO boost validation
- Efficacy scores are 0.18 lower than expected (0.60 vs 0.78)
- This is a **high impact** - significant score difference, MSI-H patients not getting proper boost

### **Proposed Fix**

Update MSI status check to accept both formats (case-insensitive):

```python
# Check MSI
msi_status = tumor_context.get("msi_status")
if msi_status:
    msi_upper = str(msi_status).upper()
    if msi_upper in ["MSI-H", "MSI-HIGH"]:  # ✅ Accept both formats
        msi_boost = 1.3
        rationale.append({
            "gate": "IO_MSI_BOOST",
            "verdict": "BOOSTED",
            "boost": 1.3,
            "msi_status": msi_status,
            "reason": f"MSI-High ({msi_status}) → Checkpoint inhibitor boost 1.3x"
        })
```

### **Validation**
After fix, re-run:
- Test Case 13: Expected 0.78, Should get 0.78 ✅
- Test Case 15: Expected 0.78, Should get 0.78 ✅

---

## 📋 FIX PRIORITY

| Bug | Priority | Impact | Effort | Scenarios Affected |
|-----|----------|--------|--------|-------------------|
| **#3: MSI String Mismatch** | **P0** | High | Low (5 min) | 2 scenarios |
| **#1: TMB Boost Factor** | **P1** | Moderate | Low (2 min) | 2 scenarios |
| **#2: Missing TMB ≥10 Tier** | **P1** | Moderate | Low (5 min) | 1 scenario |

**Total Fix Time**: ~15 minutes  
**Expected Pass Rate After Fix**: 100% (25/25 scenarios)

---

## 🔧 IMPLEMENTATION CHECKLIST

### **Fix #3: MSI String Mismatch** (P0 - Do First)
- [ ] Update line 164 in `sporadic_gates.py` to accept both `"MSI-H"` and `"MSI-High"`
- [ ] Use case-insensitive matching
- [ ] Test with Test Case 13 and 15
- [ ] Verify both scenarios pass

### **Fix #1: TMB Boost Factor** (P1)
- [ ] Update line 153 in `sporadic_gates.py` to use 1.35x instead of 1.3x
- [ ] Update rationale message to reflect 1.35x
- [ ] Test with Test Case 3 and 4
- [ ] Verify both scenarios pass

### **Fix #2: Missing TMB ≥10 Tier** (P1)
- [ ] Add TMB ≥10 but <20 tier to `sporadic_gates.py`
- [ ] Set boost to 1.25x for this tier
- [ ] Add rationale message
- [ ] Test with Test Case 9
- [ ] Verify scenario passes

### **After All Fixes**
- [ ] Re-run full validation suite
- [ ] Verify 100% pass rate (25/25 scenarios)
- [ ] Update `VALIDATION_TEST_RESULTS.md` with final results
- [ ] Mark Mission 3 as complete

---

## 📊 EXPECTED RESULTS AFTER FIXES

| Scenario | Current Status | After Fix | Expected Efficacy | Expected Confidence |
|----------|----------------|-----------|-------------------|---------------------|
| Test 3 | ❌ IO FAIL | ✅ PASS | 0.81 (TMB 22 → 1.35x) | 0.8 |
| Test 4 | ❌ IO FAIL | ✅ PASS | 0.81 (TMB 55 → 1.35x) | 0.8 |
| Test 9 | ❌ IO FAIL | ✅ PASS | 0.75 (TMB 13.5 → 1.25x) | 0.6 |
| Test 13 | ❌ IO FAIL | ✅ PASS | 0.78 (MSI-H → 1.3x) | 0.6 |
| Test 15 | ❌ IO FAIL | ✅ PASS | 0.78 (MSI-H → 1.3x) | 0.6 |

**Final Pass Rate**: 100% (25/25 scenarios) ✅

---

## 🎯 RECOMMENDATIONS FOR ZO

1. **Fix MSI String Mismatch First** (P0) - Quick win, high impact
2. **Align TMB Boost Factors** (P1) - Match Zo's A4 answer (1.35x for TMB ≥20)
3. **Add TMB ≥10 Tier** (P1) - Complete the IO boost logic per Zo's A4 answer

**All fixes are straightforward and should take ~15 minutes total.**

---

**BUG REPORT STATUS: ⚔️ 5 FAILURES DOCUMENTED, 3 FIXES PROPOSED, READY FOR ZO TO IMPLEMENT** ⚔️

