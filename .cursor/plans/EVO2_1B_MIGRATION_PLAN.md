# Evo2 1B Model Migration Plan

**Date**: January 25, 2025  
**Objective**: Switch default Evo2 model from 7B to 1B for faster inference and lower cost  
**Status**: 🚀 **IN PROGRESS**

---

## Executive Summary

**Current State**: System defaults to `evo2_7b` model  
**Target State**: System defaults to `evo2_1b` model  
**Rationale**: 1B model is faster, cheaper, and sufficient for most use cases (per Gap 2 in FAIL_NOW_VS_LATER_ASSESSMENT.md)

---

## Configuration Files to Update

### 1. Backend Config (Primary Default)

**File**: `oncology-coPilot/oncology-backend-minimal/api/config.py`  
**Line**: 158  
**Current**: `DEFAULT_EVO_MODEL = os.getenv("DEFAULT_EVO_MODEL", "evo2_7b")`  
**Change To**: `DEFAULT_EVO_MODEL = os.getenv("DEFAULT_EVO_MODEL", "evo2_1b")`

**Impact**: This is the system-wide default. All endpoints that don't specify `model_id` will use 1B.

### 2. Modal Service Config (Verify)

**File**: `src/services/evo_service/main.py`  
**Line**: 68  
**Current**: `_os.environ['EVO_MODEL_ID'] = 'evo2_1b_base'`  
**Status**: ✅ Already set to 1B! (uses `evo2_1b_base` which is the correct model name)

**Note**: The Modal service uses `evo2_1b_base` (the actual model name), while the backend uses `evo2_1b` (URL routing shorthand). Both refer to the same 1B model.

---

## Model Name Mapping

**Backend URL Routing** (config.py):
- `"evo2_1b"` → `EVO_URL_1B` (Modal service endpoint)
- `"evo2_7b"` → `EVO_URL_7B` (Modal service endpoint)
- `"evo2_40b"` → `EVO_URL_40B` (Modal service endpoint)

**Modal Service Model Names** (utils.py):
- `"evo2_1b_base"` → Uses `configs/evo2-1b-8k.yml` config
- `"evo2_7b_base"` → Uses `configs/evo2-7b-8k.yml` config
- `"evo2_40b_base"` → Uses `configs/evo2-40b-8k.yml` config

**Mapping**: Backend `"evo2_1b"` routes to Modal service that loads `"evo2_1b_base"` model.

---

## Test Plan

### Test 1: Config Default Verification

**Script**: `scripts/test_evo2_1b_config.py` (NEW)

**Test**:
1. Import `DEFAULT_EVO_MODEL` from `api.config`
2. Assert `DEFAULT_EVO_MODEL == "evo2_1b"`
3. Verify environment variable override works

### Test 2: Endpoint Default Behavior

**Script**: `scripts/test_evo2_1b_endpoint.py` (NEW)

**Test**:
1. Call `/api/evo/score_variant` without `model_id` parameter
2. Verify response uses 1B model (check response metadata)
3. Verify scoring works correctly

### Test 3: MBD4+TP53 Variant Scoring

**Script**: `scripts/test_evo2_1b_mbd4_tp53.py` (NEW)

**Test**:
1. Score MBD4 frameshift variant with 1B model
2. Score TP53 R175H hotspot with 1B model
3. Verify scores are reasonable (MBD4 ≥0.8, TP53 ≥0.7)
4. Compare to 7B model scores (if available) for validation

### Test 4: Performance Comparison

**Script**: `scripts/test_evo2_1b_performance.py` (NEW)

**Test**:
1. Measure latency for 1B vs 7B (if both available)
2. Verify 1B is faster
3. Verify 1B scores are within acceptable range of 7B

---

## Implementation Steps

### Step 1: Update Config File ✅

**File**: `oncology-coPilot/oncology-backend-minimal/api/config.py`  
**Change**: Line 158

```python
# BEFORE:
DEFAULT_EVO_MODEL = os.getenv("DEFAULT_EVO_MODEL", "evo2_7b")

# AFTER:
DEFAULT_EVO_MODEL = os.getenv("DEFAULT_EVO_MODEL", "evo2_1b")
```

### Step 2: Create Test Scripts

**Files to Create**:
1. `scripts/test_evo2_1b_config.py` - Config verification
2. `scripts/test_evo2_1b_endpoint.py` - Endpoint default behavior
3. `scripts/test_evo2_1b_mbd4_tp53.py` - MBD4+TP53 validation
4. `scripts/test_evo2_1b_performance.py` - Performance comparison

### Step 3: Run Tests

**Commands**:
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/test_evo2_1b_config.py
python scripts/test_evo2_1b_endpoint.py
python scripts/test_evo2_1b_mbd4_tp53.py
python scripts/test_evo2_1b_performance.py
```

### Step 4: Verify Modal Service

**Check**: Verify Modal service `EVO_URL_1B` is accessible and returns correct model.

**Command**:
```bash
curl -X POST https://crispro--evo-service-evoservice1b-api-1b.modal.run/score_variant \
  -H "Content-Type: application/json" \
  -d '{"chrom": "3", "pos": 129430456, "ref": "A", "alt": ""}'
```

---

## Rollback Plan

If issues arise, revert config change:

```python
DEFAULT_EVO_MODEL = os.getenv("DEFAULT_EVO_MODEL", "evo2_7b")
```

Or set environment variable:
```bash
export DEFAULT_EVO_MODEL=evo2_7b
```

---

## Success Criteria

- [ ] Config file updated to `evo2_1b`
- [ ] All test scripts pass
- [ ] MBD4+TP53 variants score correctly with 1B
- [ ] Endpoints default to 1B when `model_id` not specified
- [ ] Performance is acceptable (1B faster than 7B)
- [ ] No regressions in existing functionality

---

## Notes

- **1B Model Benefits**: Faster inference, lower cost, sufficient accuracy for most variants
- **7B Model**: Still available via explicit `model_id="evo2_7b"` parameter
- **40B Model**: Available for high-accuracy requirements via `model_id="evo2_40b"`
- **Fallback Logic**: System will fallback to 7B or 40B if 1B unavailable (existing logic)

---

**Status**: ✅ Plan created, ready to execute



