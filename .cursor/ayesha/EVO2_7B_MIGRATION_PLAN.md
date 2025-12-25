# evo2_7b Migration Plan - Switch from Random to Trained SAE Weights

**Date**: January 20, 2025  
**Purpose**: Migrate from evo2_1b (random SAE weights) to evo2_7b (trained SAE weights)  
**Status**: Ready to implement

---

## 🎯 OBJECTIVE

**Current Problem**:
- Using evo2_1b (1920-dim activations)
- SAE checkpoint is 4096×32768 (for evo2_7b/40b)
- Dimension mismatch → random SAE weights used
- Features may be noise, not biologically meaningful

**Solution**:
- Switch to evo2_7b (4096-dim activations)
- SAE checkpoint matches (4096×32768)
- Load trained SAE weights from Hugging Face
- Features will be biologically meaningful

---

## 📝 CODE CHANGES REQUIRED

### 1. Modal SAE Service (`src/services/sae_service/main.py`)

**Change 1: Default Model ID** (Line 155)
```python
# BEFORE:
model_id = os.getenv("EVO_MODEL_ID", "evo2_1b_base")

# AFTER:
model_id = os.getenv("EVO_MODEL_ID", "evo2_7b_base")
```

**Change 2: Fallback Dimension** (Line 213)
```python
# BEFORE:
d_in_detected = 1920  # Typical for evo2_1b_base

# AFTER:
d_in_detected = 4096  # Typical for evo2_7b_base
```

**Change 3: Enable Checkpoint Loading** (Lines 224-230)
```python
# BEFORE (DISABLED):
# CHECKPOINT LOADING COMPLETELY DISABLED FOR RUO SAFETY
# The Goodfire checkpoint is 4096×32768 but evo2_1b_base emits 1920-dim activations.
# For RUO exploratory analysis, we use a randomly initialized SAE that matches the detected dim.
logger.warning(f"⚠️  Using randomly initialized SAE ({d_in_detected}×32768) for RUO exploratory analysis.")
logger.warning(f"⚠️  Goodfire pre-trained checkpoint (4096×32768) intentionally NOT loaded due to dimension mismatch.")
# Explicitly set sae_weights_path to None to ensure no accidental loading
sae_weights_path = None

# AFTER (ENABLED):
# Load trained SAE weights if dimension matches
if sae_weights_path and d_in_detected == 4096:
    try:
        logger.info(f"Loading SAE checkpoint from: {sae_weights_path}")
        checkpoint = torch.load(sae_weights_path, map_location=device, weights_only=True)
        
        # Strip prefix if needed (matches notebook pattern)
        new_dict = {}
        for key, item in checkpoint.items():
            new_key = key.replace("_orig_mod.", "").replace("module.", "")
            new_dict[new_key] = item
        
        # Load state dict (strict=False to handle minor mismatches)
        self.sae_model.load_state_dict(new_dict, strict=False)
        logger.info("✅ Loaded trained SAE weights (4096×32768) from Goodfire/Evo-2-Layer-26-Mixed")
    except Exception as e:
        logger.warning(f"Failed to load checkpoint: {e}. Using random initialization.")
        logger.warning(f"⚠️  Using randomly initialized SAE ({d_in_detected}×32768)")
else:
    if d_in_detected != 4096:
        logger.warning(f"⚠️  Dimension mismatch: detected {d_in_detected}, checkpoint expects 4096. Using random initialization.")
    else:
        logger.warning(f"⚠️  SAE checkpoint not available. Using random initialization ({d_in_detected}×32768)")
```

**Change 4: Update Provenance** (Line 369)
```python
# BEFORE:
"model": "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)"

# AFTER:
if sae_weights_path and d_in_detected == 4096:
    model_desc = "Goodfire/Evo-2-Layer-26-Mixed (trained weights)"
else:
    model_desc = f"Goodfire/Evo-2-Layer-26-Mixed (random init, d_in={d_in_detected})"
```

### 2. Backend Router (`api/routers/sae.py`)

**Change: Default Model ID** (Line 51)
```python
# BEFORE:
model_id: str = Field("evo2_1b", description="Evo2 model to use")

# AFTER:
model_id: str = Field("evo2_7b", description="Evo2 model to use")
```

### 3. Cohort Extraction Script (`scripts/sae/extract_sae_features_cohort.py`)

**Check if model_id is hardcoded**:
- If hardcoded to "evo2_1b", change to "evo2_7b"
- Or allow model_id to be passed as argument

---

## 🧪 TESTING PLAN

### Step 1: Test SAE Service Initialization

**Command**:
```bash
# Test with evo2_7b
curl -X POST $SAE_SERVICE_URL/extract_features \
  -H "Content-Type: application/json" \
  -d '{
    "chrom": "17",
    "pos": 43044295,
    "ref": "T",
    "alt": "G",
    "assembly": "GRCh38",
    "window": 8192,
    "model_id": "evo2_7b"
  }'
```

**Expected Response**:
```json
{
  "top_features": [...],
  "provenance": {
    "method": "batch_topk_tied_sae",
    "d_in": 4096,
    "d_hidden": 32768,
    "k": 64,
    "model": "Goodfire/Evo-2-Layer-26-Mixed (trained weights)",
    "sae_version": "v1",
    "source": "modal_sae_service"
  }
}
```

**Validation**:
- ✅ `d_in` should be 4096 (not 1920)
- ✅ `model` should say "trained weights" (not "random init")
- ✅ Features should be different from random weights

### Step 2: Test Small Batch (1-2 Patients)

**Command**:
```bash
# Extract features for 1-2 patients only (cost control)
python3 scripts/sae/extract_sae_features_cohort.py \
  --input data/validation/tcga_ov_platinum_response_labels.json \
  --output data/validation/sae_cohort/sae_features_tcga_ov_platinum_test.json \
  --max-patients 2 \
  --model-id evo2_7b
```

**Validation**:
- ✅ Features extracted successfully
- ✅ Feature distributions different from random weights
- ✅ Provenance shows "trained weights"
- ✅ No errors or crashes

### Step 3: Verify Feature Distributions

**Command**:
```bash
python3 -c "
import json
import numpy as np

with open('data/validation/sae_cohort/sae_features_tcga_ov_platinum_test.json', 'r') as f:
    data = json.load(f)

# Check feature values
values = [f['value'] for p in data['patients'] for v in p['variants'] for f in v.get('top_features', [])]
print(f'Feature value range: {min(values):.2e} to {max(values):.2e}')
print(f'Mean: {np.mean(values):.2e}, Std: {np.std(values):.2e}')

# Compare with old data (if available)
# Old (random): Mean=1.62e+15, Std=6.32e+14
# New (trained): Should be different
"
```

**Expected**:
- Feature values should be different from random weights
- May be normalized or in different range
- Should show biological patterns (not pure noise)

### Step 4: Full Cohort Re-Extraction (After Testing)

**Command**:
```bash
# Full extraction with evo2_7b
python3 scripts/sae/extract_sae_features_cohort.py \
  --input data/validation/tcga_ov_platinum_response_labels.json \
  --output data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --model-id evo2_7b
```

**Cost Control**:
- ✅ Use feature flags (`ENABLE_TRUE_SAE`) to gate extraction
- ✅ Circuit breaker implemented (max errors per 100 requests)
- ✅ Monitor Modal usage/costs
- ✅ Test with small batch first

---

## ⚠️ COST CONSIDERATIONS

**User Requirement**: "keep it tight where we always test to make sure we dont burn through credits - we have flags implemented as well"

**Mitigation**:
1. ✅ **Feature Flags**: `ENABLE_TRUE_SAE` gates SAE extraction
2. ✅ **Circuit Breaker**: Max 50 errors per 100 requests (lines 238-257)
3. ✅ **Test First**: Small batch (1-2 patients) before full extraction
4. ✅ **Monitor Usage**: Check Modal dashboard for costs
5. ✅ **Checkpointing**: Script supports checkpointing (resume if interrupted)

**Recommendations**:
- Test with 1-2 patients first
- Verify features are different (not noise)
- Monitor Modal costs during test
- Proceed with full extraction only if test passes

---

## 📊 EXPECTED CHANGES

### Feature Distributions

**Before (Random Weights)**:
- Mean: 1.62e+15
- Std: 6.32e+14
- Range: 1.25e+14 to 5.58e+15
- May be noise (not biologically meaningful)

**After (Trained Weights)**:
- Values may be normalized (different range)
- Should show biological patterns
- Features should correlate with known biology
- More interpretable

### Provenance

**Before**:
```json
{
  "model": "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)",
  "d_in": 1920
}
```

**After**:
```json
{
  "model": "Goodfire/Evo-2-Layer-26-Mixed (trained weights)",
  "d_in": 4096
}
```

---

## ✅ SUCCESS CRITERIA

**Migration Successful If**:
1. ✅ SAE service loads evo2_7b model
2. ✅ Dimension detected as 4096 (not 1920)
3. ✅ Checkpoint loads successfully (no errors)
4. ✅ Provenance shows "trained weights"
5. ✅ Features extracted successfully
6. ✅ Feature distributions different from random weights
7. ✅ No crashes or errors

**Biomarker Analysis Successful If** (After Re-Extraction):
1. ✅ ≥10 significant features found (p < 0.05, FDR-corrected)
2. ✅ Correlations |r| ≥ 0.3 (moderate strength)
3. ✅ Effect sizes Cohen's d ≥ 0.5
4. ✅ Features show biological patterns (not noise)

---

## 📋 IMPLEMENTATION CHECKLIST

- [ ] Update `src/services/sae_service/main.py`:
  - [ ] Change default model to `"evo2_7b_base"`
  - [ ] Update fallback dimension to 4096
  - [ ] Enable checkpoint loading
  - [ ] Update provenance message
- [ ] Update `api/routers/sae.py`:
  - [ ] Change default model_id to `"evo2_7b"`
- [ ] Check `scripts/sae/extract_sae_features_cohort.py`:
  - [ ] Update model_id if hardcoded
- [ ] Test SAE service initialization
- [ ] Test small batch (1-2 patients)
- [ ] Verify feature distributions
- [ ] Monitor costs
- [ ] Proceed with full extraction (if test passes)
- [ ] Re-run biomarker analysis
- [ ] Document findings

---

**Status**: ✅ **READY TO IMPLEMENT**

**Next Step**: Make code changes, test with small batch, verify features are different from random weights.

