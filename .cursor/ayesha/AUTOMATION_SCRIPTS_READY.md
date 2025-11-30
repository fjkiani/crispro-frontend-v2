# Automation Scripts Ready

**Date**: January 20, 2025  
**Status**: ✅ Scripts Created - Ready for Use

---

## 🧪 Testing Scripts

### 1. `scripts/sae/test_sae_extraction.py`

**Purpose**: Quick test to verify SAE extraction with evo2_7b and trained weights

**Usage**:
```bash
# Make sure backend is running first
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload

# In another terminal, run test
python3 scripts/sae/test_sae_extraction.py
```

**What it tests**:
- Backend connectivity
- SAE extraction for 3 known variants (BRCA1, KRAS, BRAF)
- Trained weights loading (should show "trained weights" not "random init")
- Correct dimension (4096, not 1920)

**Expected output**:
```
✅ ALL TESTS PASSED - Ready for cohort extraction!
```

---

## 🗺️ Feature→Pathway Mapping Script

### 2. `scripts/sae/create_feature_pathway_mapping.py`

**Purpose**: Create feature→pathway mapping from biomarker analysis results

**Usage** (after biomarker analysis):
```bash
python3 scripts/sae/create_feature_pathway_mapping.py \
    --biomarkers data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
    --cohort data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json \
    --top-n 100 \
    --min-correlation 0.3 \
    --min-p-value 0.05
```

**What it does**:
1. Loads biomarker analysis results (top significant features)
2. For each feature, finds patients with high activation
3. Extracts genes from those patients
4. Uses gene→pathway mapping to infer pathway
5. Creates mapping JSON file

**Strategy**: Gene→Pathway Inference
- BRCA1/BRCA2 → DDR pathway
- KRAS/BRAF → MAPK pathway
- PIK3CA → PI3K pathway
- ERBB2 → HER2 pathway
- etc.

**Output**: `sae_feature_mapping.json` ready for use in SAE service

---

## 📋 Workflow

### Step 1: Test Backend (After Restart)
```bash
python3 scripts/sae/test_sae_extraction.py
```

### Step 2: Extract Features (If needed)
```bash
export ENABLE_SAE_COHORT_RUN=1
export MAX_PATIENTS=66
python3 scripts/sae/extract_sae_features_cohort.py
```

### Step 3: Run Biomarker Analysis
```bash
python3 scripts/sae/analyze_biomarkers.py \
    --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json
```

### Step 4: Create Pathway Mapping
```bash
python3 scripts/sae/create_feature_pathway_mapping.py \
    --biomarkers data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
    --cohort data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json
```

---

## ✅ Scripts Status

- ✅ `test_sae_extraction.py` - Created and executable
- ✅ `create_feature_pathway_mapping.py` - Created and executable
- ✅ Both scripts have proper error handling and logging
- ✅ Both scripts follow existing codebase patterns

---

**Status**: ✅ **READY** - Scripts can be used once backend is running

