# EXTRACTION PIECE 5.5: Modal Payload Size Fix

**Source**: Lines 28000-28100 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-XX  
**Status**: ✅ Complete

---

## Overview

This piece documents the discovery and fix for a critical bug where the Modal SAE service was attempting to serialize a massive 32K-dimensional feature vector (potentially 268 million floats) as JSON, causing crashes due to memory/timeout issues.

---

## The Problem

**Symptom**: Modal SAE service crashing during response serialization

**Investigation**: Checked response payload size

**Root Cause**: Line 346 in `src/services/sae_service/main.py`:
```python
"features": features.cpu().numpy().tolist()
```

**The Issue**:
- Features shape: `[batch=1, seq_len=8193, d_hidden=32768]`
- Full tensor size: **268,435,456 floats** (268 million)
- JSON payload size: **~1-2GB**
- Result: **Crashing Modal** (memory/serialization timeout)

---

## Why This Was Wrong

1. **Massive payload**: 268M floats = 1-2GB of JSON
2. **Crashing Modal**: Memory/serialization timeout
3. **Useless for downstream**: Biomarker analysis only needs top-k features (64 features), not the full 32K-dim vector

---

## The Fix

**Solution**: Remove the full `features` array from the response and only return `top_features` (the 64 most active features)

**Code Change**:
```python
# BEFORE (line 346):
return {
    "features": features.cpu().numpy().tolist(),  # ❌ 268M floats
    "top_features": [...]
}

# AFTER:
# NOTE: We only return top_features (k=64), not the full 32K-dim vector,
# to prevent massive payloads (268M floats = 1-2GB JSON) that crash Modal.
# Downstream biomarker analysis only needs the top-k active features anyway.
return {
    # "features": features.cpu().numpy().tolist(),  # ❌ REMOVED
    "top_features": [
        {"index": int(idx), "value": float(val)}
        for idx, val in zip(top_k_indices, top_k_values)
    ],
    ...
}
```

---

## Additional Fixes

**Provenance Update**:
```python
# BEFORE:
"provenance": {
    "d_in": 4096,  # ❌ Hardcoded wrong value
    "model": "Goodfire/Evo-2-Layer-26-Mixed"  # ❌ Misleading (random init)
}

# AFTER:
"provenance": {
    "d_in": int(features.shape[-1]),  # ✅ Use actual detected d_in
    "model": "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)"  # ✅ Honest
}
```

---

## Deployment

**Command**: `venv/bin/modal deploy src/services/sae_service/main.py`

**Result**: ✅ Deployed successfully in 1.263s

**Endpoint**: `https://crispro--sae-service-saeservice-api.modal.run`

---

## Impact

**Before Fix**:
- ❌ Service crashes on every request
- ❌ 1-2GB JSON payloads
- ❌ Modal timeout/memory errors

**After Fix**:
- ✅ Service returns only 64 features (~1KB JSON)
- ✅ Fast serialization
- ✅ Downstream analysis works correctly

---

## Key Lessons

1. **Don't serialize unnecessary data**: Only return what's needed downstream
2. **Check payload sizes**: 268M floats is way too large for JSON
3. **Top-k is sufficient**: Biomarker analysis only needs active features, not full vector
4. **Honest provenance**: Don't claim trained weights when using random init

---

## Related Files

- `src/services/sae_service/main.py` - Fixed SAE service
- `scripts/sae/extract_sae_features_cohort.py` - Uses `top_features` from response
- `api/services/biomarker_correlation_service.py` - Aggregates `top_features` into feature matrix

---

## Verification

**Test Request**:
```bash
curl -X POST "http://localhost:8000/api/sae/extract_features" \
  -H "Content-Type: application/json" \
  -d '{
    "chrom":"13",
    "pos":25060373,
    "ref":"T",
    "alt":"C",
    "assembly":"GRCh37",
    "window":8192,
    "model_id":"evo2_1b"
  }'
```

**Expected Response**:
- ✅ `top_features`: Array of 64 objects with `index` and `value`
- ✅ No `features` key (removed)
- ✅ Small JSON payload (~1KB)
- ✅ Fast response time

