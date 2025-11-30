# Modal Service Redeploy Needed

**Date**: January 20, 2025  
**Status**: ⚠️ **Code Fixed, Redeploy Required**

---

## 🔍 Issue Found

**Test Results**:
- ✅ Backend is running and reachable
- ❌ Modal SAE service still using **random weights** (not trained)
- ❌ Provenance shows wrong `d_in` (32768 instead of 4096)
- ⚠️ Backend is correctly proxying to Modal service

**Root Cause**:
1. Modal service code was updated but **not redeployed**
2. Bug in provenance: using `features.shape[-1]` (32768) instead of `d_in_detected` (4096)

---

## ✅ Code Fixes Applied

### Fix 1: Provenance `d_in` Bug
**File**: `src/services/sae_service/main.py:389`

**Before**:
```python
"d_in": int(features.shape[-1]),  # Wrong: returns 32768 (feature dim)
```

**After**:
```python
"d_in": int(self.d_in_detected),  # Correct: returns 4096 (Evo2 input dim)
```

### Fix 2: Store `d_in_detected` as Instance Variable
**File**: `src/services/sae_service/main.py:253-255`

**Added**:
```python
self.d_in_detected = d_in_detected  # Store for provenance
```

---

## 🚀 Next Steps: Redeploy Modal Service

### Option 1: Manual Redeploy (Recommended)
```bash
cd src/services/sae_service
modal deploy main.py
```

### Option 2: Check Deployment Status
```bash
# Check if service is deployed
modal app list

# Check service logs
modal app logs crispro sae-service
```

### Option 3: Verify Deployment
After redeploy, test again:
```bash
python3 scripts/sae/test_sae_extraction.py
```

**Expected After Redeploy**:
- ✅ Provenance shows `d_in: 4096` (not 32768)
- ✅ Provenance shows "trained weights" (not "random init")
- ✅ All tests pass

---

## 📊 Current vs Expected

### Current (Before Redeploy)
```json
{
  "d_in": 32768,  // ❌ Wrong (feature dimension)
  "model": "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)"  // ❌ Random weights
}
```

### Expected (After Redeploy)
```json
{
  "d_in": 4096,  // ✅ Correct (Evo2 input dimension)
  "model": "Goodfire/Evo-2-Layer-26-Mixed (trained weights)"  // ✅ Trained weights
}
```

---

## ⚠️ Important Notes

1. **Backend doesn't need restart** - it's just proxying to Modal
2. **Modal service must be redeployed** - code changes won't take effect until redeploy
3. **Checkpoint loading** - Will work after redeploy if:
   - `d_in_detected == 4096` (evo2_7b)
   - SAE checkpoint is available at configured path
   - Checkpoint dimensions match (4096×32768)

---

**Status**: ✅ **DEPLOYED AND VERIFIED** - See `.cursor/ayesha/DEPLOYMENT_SUCCESS.md`

