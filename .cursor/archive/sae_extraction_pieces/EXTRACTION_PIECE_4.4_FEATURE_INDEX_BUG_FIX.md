# 📖 EXTRACTION PIECE 4.4: Feature Index Bug Fix
**Source**: Lines 30900-31000 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents the critical bug discovery where SAE feature indices were incorrectly flattened from a 3D tensor, the fix implementation, verification, and re-extraction decision.

---

## 🔍 KEY FINDINGS

### **Bug Discovery**

**Problem**: Existing 69-patient dataset had **wrong feature indices**

**Root Cause**: 
- SAE features were flattened from entire 3D tensor `[1, 8193, 32768]`
- Indices were flattened positions (e.g., 183M, 254M) instead of SAE feature indices (0-32767)
- This made all feature matrices zero-filled (indices out of range)

**Impact**:
- All biomarker analysis results were **invalid**
- Feature matrices were all zeros
- No meaningful correlations could be computed

---

### **Fix Implementation**

**Location**: `src/services/sae_service/main.py`

**Change**:
```python
# OLD (WRONG):
features_flat = features.flatten()  # [1 * 8193 * 32768] = flattened indices
top_k_values, top_k_indices = torch.topk(features_flat, k=64)

# NEW (CORRECT):
# Aggregate across sequence positions (mean pooling) FIRST
features_aggregated = features.mean(dim=1).squeeze(0)  # [32768]
# Then get top-k SAE features (indices in 0-32767 range)
top_k_values, top_k_indices = torch.topk(features_aggregated, k=64)
```

**Key Fix**:
1. **Aggregate first**: Mean pool across sequence dimension → single 32K-dim vector
2. **Then top-k**: Apply `torch.topk` to aggregated vector → indices in 0-32767 range

---

### **Verification**

**Test Variant**: `13:25060373 T>C`

**Results**:
```
Number of top features: 64
First 5 indices:
  Index: 17250, Value: 2148608124125184.0
  Index: 27857, Value: 1721234182111232.0
  Index: 29831, Value: 1441678653128704.0
  Index: 23396, Value: 1269208201560064.0
  Index: 20867, Value: 1168599666393088.0
Max index: 32573
Min index: 416
Indices in valid range (0-32767)? True ✅
```

**Verification**: ✅ All indices now in correct range (0-32767)

---

### **Re-Extraction Decision**

**Options Considered:**

**Option A**: Re-extract all 69 patients with fixed Modal service (2-3 hours)
- ✅ **Pros**: Clean, correct data
- ⚠️ **Cons**: Time-consuming

**Option B**: Fix indices in existing data file (quick hack)
- ✅ **Pros**: Fast
- ❌ **Cons**: Fragile, may miss edge cases

**Option C**: Test fix first, then decide on full re-extraction
- ✅ **Pros**: Validates fix before committing time
- ✅ **Cons**: None

**Decision**: **Option C** → Test fix first, then re-extract

**Action Taken**:
1. ✅ Redeployed Modal service with fix
2. ✅ Tested on single variant → verified indices correct
3. ✅ Deleted old incorrect data
4. ✅ Started fresh extraction with corrected service

---

### **Re-Extraction Process**

**Command**:
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main && \
rm -f data/validation/sae_cohort/sae_features_extraction_checkpoint.json && \
rm -f data/validation/sae_cohort/sae_features_tcga_ov_platinum.json && \
ENABLE_EVO2_SAE=1 ENABLE_TRUE_SAE=1 ENABLE_SAE_COHORT_RUN=1 \
MAX_PATIENTS=100 MAX_TOTAL_VARIANTS=5000 \
python3 scripts/sae/extract_sae_features_cohort.py
```

**Status**: Extraction started with corrected service

**Expected**: 2-3 hours for 100 patients

---

## 📊 KEY INSIGHTS

### **Bug Root Cause**

1. **3D Tensor Confusion**: Features shape `[batch=1, seq_len=8193, d_hidden=32768]`
2. **Flattening Error**: Flattened entire tensor instead of aggregating first
3. **Index Mismatch**: Flattened indices (0 to 268M) vs SAE feature indices (0-32767)
4. **Silent Failure**: All-zero feature matrices, no error raised

### **Fix Strategy**

1. **Aggregate First**: Mean pool across sequence dimension
2. **Then Index**: Apply top-k to aggregated vector
3. **Verify**: Test on single variant before full re-extraction
4. **Clean Slate**: Delete old data, start fresh

### **Lessons Learned**

1. **Always Verify**: Test indices are in expected range
2. **Understand Dimensions**: Know tensor shapes before operations
3. **Test Before Scale**: Verify fix on small sample before full extraction
4. **Clean Data**: Delete incorrect data rather than trying to fix in place

---

## 🔗 CONTEXT & CONNECTIONS

- **Blocks**: All biomarker analysis until fixed
- **Related to**: Bug discovery (Piece 3.3), Circuit breaker (Piece 3.4)
- **Enables**: Valid biomarker analysis after re-extraction
- **Key Insight**: Critical bug that invalidated all previous analysis

---

## 📝 NOTES

- Bug discovered before biomarker analysis ran
- Fix verified on single variant before full re-extraction
- Old data deleted to ensure clean extraction
- Re-extraction started with corrected service

---

## 🎯 QUESTIONS RESOLVED

- ✅ What was the bug? → Flattened 3D tensor indices instead of SAE feature indices
- ✅ How was it fixed? → Aggregate across sequence first, then top-k
- ✅ Was fix verified? → Yes, tested on single variant, indices correct
- ✅ What's next? → Re-extract all patients with corrected service

