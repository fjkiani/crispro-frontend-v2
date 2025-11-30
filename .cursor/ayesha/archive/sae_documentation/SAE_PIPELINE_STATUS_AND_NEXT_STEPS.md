# SAE Pipeline Status & Next Steps

**Date**: January 20, 2025  
**Status**: ✅ Code Complete, ⏸️ Awaiting Backend Restart & Testing

---

## ✅ COMPLETED

### 1. evo2_7b Migration (Code Changes)
- ✅ **SAE Service** (`src/services/sae_service/main.py`):
  - Default model: `evo2_1b_base` → `evo2_7b_base`
  - Fallback dimension: `1920` → `4096`
  - Checkpoint loading: **ENABLED** (trained weights will load)
  - Provenance tracking: Added `sae_weights_loaded` flag

- ✅ **Backend Router** (`api/routers/sae.py`):
  - Default model_id: `"evo2_1b"` → `"evo2_7b"`

- ✅ **Cohort Script** (`scripts/sae/extract_sae_features_cohort.py`):
  - Default model_id: `"evo2_1b"` → `"evo2_7b"` (4 locations)

### 2. Environment Configuration
- ✅ **Backend `.env` file updated**:
  - `ENABLE_TRUE_SAE=1`
  - `ENABLE_EVO2_SAE=1`
  - `SAE_SERVICE_URL=https://crispro--sae-service-saeservice-api.modal.run`
  - `EVO_URL_7B=https://crispro--evo-service-evoservice7b-api-7b.modal.run`

### 3. Documentation
- ✅ Modal service URLs documented
- ✅ Testing plan created
- ✅ Next steps documented

---

## ⏸️ REQUIRED: Backend Restart

**The backend needs to be restarted to pick up the new `.env` variables.**

**Location**: `oncology-coPilot/oncology-backend-minimal/.env`

**After restart, verify with**:
```bash
curl -X POST http://localhost:8000/api/sae/extract_features \
  -H "Content-Type: application/json" \
  -d '{
    "chrom": "17",
    "pos": 43044295,
    "ref": "T",
    "alt": "G",
    "assembly": "GRCh38",
    "window": 8192,
    "model_id": "evo2_7b"
  }' | jq '.provenance.model'
```

**Expected**: `"Goodfire/Evo-2-Layer-26-Mixed (trained weights)"` (not "random init")

---

## 🧪 TESTING PLAN (After Backend Restart)

### Step 1: Single Variant Test ✅
**Purpose**: Verify trained weights load correctly

**Command**:
```bash
python3 -c "
import asyncio
import httpx

async def test():
    async with httpx.AsyncClient(timeout=180.0) as client:
        response = await client.post(
            'http://localhost:8000/api/sae/extract_features',
            json={
                'chrom': '17',
                'pos': 43044295,
                'ref': 'T',
                'alt': 'G',
                'assembly': 'GRCh38',
                'window': 8192,
                'model_id': 'evo2_7b'
            }
        )
        if response.status_code == 200:
            result = response.json()
            prov = result.get('provenance', {})
            print(f'Model: {prov.get(\"model\")}')
            print(f'd_in: {prov.get(\"d_in\")}')
            if 'trained weights' in prov.get('model', '').lower():
                print('✅ TRAINED WEIGHTS LOADED!')
            else:
                print('⚠️  Still using random weights')
        else:
            print(f'Error: {response.text}')

asyncio.run(test())
"
```

**Success Criteria**:
- ✅ Status 200
- ✅ Provenance shows "trained weights"
- ✅ `d_in` is 4096 (not 1920)

### Step 2: Small Batch Test (1 patient)
**Purpose**: Verify full pipeline works with trained weights

**Command**:
```bash
export ENABLE_SAE_COHORT_RUN=1
export MAX_PATIENTS=1
export MAX_TOTAL_VARIANTS=10

python3 scripts/sae/extract_sae_features_cohort.py
```

**Success Criteria**:
- ✅ Features extracted successfully
- ✅ Provenance shows "trained weights"
- ✅ Feature distributions different from random weights

### Step 3: Full Cohort Re-Extraction (66 patients)
**Purpose**: Re-extract all features with trained weights

**Command**:
```bash
export MAX_PATIENTS=66  # Or remove limit
export MAX_TOTAL_VARIANTS=3000

python3 scripts/sae/extract_sae_features_cohort.py
```

**Note**: This will overwrite existing cohort file. Consider backing up first.

---

## 📊 BIOMARKER ANALYSIS (After Re-Extraction)

### Step 4: Re-Run Biomarker Analysis
**Purpose**: Find significant SAE features correlated with platinum response

**Command**:
```bash
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json
```

**Expected Output**:
- Top 100 features with p < 0.01
- Correlation statistics (Pearson r, Cohen's d)
- CV stability scores
- Visualization plots

**Success Criteria**:
- ✅ ≥10 significant features found (p < 0.05, FDR-corrected)
- ✅ Correlations |r| ≥ 0.3 (moderate strength)
- ✅ Effect sizes Cohen's d ≥ 0.5

---

## 🗺️ FEATURE→PATHWAY MAPPING (After Biomarker Analysis)

### Step 5: Create Feature→Pathway Mapping
**Purpose**: Map significant SAE features to biological pathways

**Strategy**: Gene→Pathway Inference

1. **For each significant feature**:
   - Find which patients have it activated (top 10% feature values)
   - Extract genes/mutations from those patients
   - Use gene→pathway mapping to infer pathway

2. **Gene→Pathway Mappings** (from existing code):
   - **DDR**: BRCA1, BRCA2, RAD51, ATM, ATR, CHEK2, PALB2, FANCA, FANCD2
   - **MAPK**: KRAS, NRAS, BRAF, MAP2K1, MAP2K2, MAPK1, MAPK3
   - **PI3K**: PIK3CA, PIK3CB, AKT1, AKT2, MTOR, PTEN, TSC1, TSC2
   - **VEGF**: VEGFA, VEGFR1, VEGFR2, HIF1A, FLT1, KDR
   - **HER2**: ERBB2, HER2
   - **TP53**: TP53, MDM2, MDM4, CDKN1A, CDKN2A, RB1, E2F1

3. **Create Mapping File**:
   - Template: `api/resources/sae_feature_mapping_template.json`
   - Output: `api/resources/sae_feature_mapping.json`
   - Populate `feature_indices` from biomarker results
   - Add correlation statistics
   - Set confidence levels

**Validation**:
- Test cases: BRCA1 → DDR (high), KRAS → MAPK (high), HER2 → HER2 (high)
- Success: ≥80% of known cases map correctly

---

## 📋 EXPECTED RESULTS

### Before (Random Weights)
- Model: "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)"
- d_in: 1920
- Features: May be noise, not biologically meaningful

### After (Trained Weights)
- Model: "Goodfire/Evo-2-Layer-26-Mixed (trained weights)"
- d_in: 4096
- Features: Biologically meaningful, interpretable

---

## 🎯 SUCCESS METRICS

1. **Migration Success**:
   - ✅ Trained weights load correctly
   - ✅ Features extracted successfully
   - ✅ Feature distributions different from random

2. **Biomarker Analysis Success**:
   - ✅ ≥10 significant features found
   - ✅ Strong correlations (|r| ≥ 0.3)
   - ✅ Large effect sizes (Cohen's d ≥ 0.5)

3. **Pathway Mapping Success**:
   - ✅ Mapping file created
   - ✅ ≥80% validation accuracy on known cases
   - ✅ Pathway scores match expected biology

---

## 📚 REFERENCE DOCUMENTS

- **Migration Plan**: `.cursor/ayesha/EVO2_7B_MIGRATION_COMPLETE.md`
- **Testing Guide**: `.cursor/ayesha/READY_FOR_TESTING.md`
- **Next Steps**: `.cursor/ayesha/NEXT_STEPS_AFTER_BACKEND_RESTART.md`
- **Modal URLs**: `.cursor/ayesha/MODAL_SERVICE_URLS.md`
- **Template**: `oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping_template.json`

---

**Status**: ✅ **DEPLOYED AND VERIFIED** - Modal service deployed, trained weights confirmed, small batch test passed

**Latest**: 
- ✅ Modal service redeployed with fixed code
- ✅ Trained weights loading confirmed (`d_in: 4096`, "trained weights")
- ✅ Small batch test passed (1 patient, 47/50 variants)
- ✅ End-to-end pipeline verified
- ⏸️ Ready for full cohort extraction (awaiting approval)

**New**: Quick start guide and automation scripts available:
- `scripts/sae/test_sae_extraction.py` - Quick verification test
- `scripts/sae/create_feature_pathway_mapping.py` - Pathway mapping automation
- See `.cursor/ayesha/QUICK_START_GUIDE.md` for step-by-step instructions

