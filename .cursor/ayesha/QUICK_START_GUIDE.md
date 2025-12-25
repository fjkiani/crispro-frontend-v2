# Quick Start Guide - SAE Pipeline Testing

**Date**: January 20, 2025  
**Status**: ✅ Ready to Test

---

## 🚀 Quick Start (5 Steps)

### Step 1: Start Backend
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload
```

**Verify**: Backend should start and show "Application startup complete"

---

### Step 2: Test SAE Extraction
```bash
# In another terminal
python3 scripts/sae/test_sae_extraction.py
```

**Expected**: 
- ✅ Backend is reachable
- ✅ All 3 tests pass
- ✅ Trained weights loaded
- ✅ Correct dimension (4096)

**If fails**: Check backend `.env` file has:
- `ENABLE_TRUE_SAE=1`
- `ENABLE_EVO2_SAE=1`
- `SAE_SERVICE_URL=https://crispro--sae-service-saeservice-api.modal.run`

---

### Step 3: Small Batch Test (Optional)
```bash
export ENABLE_SAE_COHORT_RUN=1
export MAX_PATIENTS=1
export MAX_TOTAL_VARIANTS=10

python3 scripts/sae/extract_sae_features_cohort.py
```

**Purpose**: Verify full pipeline works with 1 patient

---

### Step 4: Full Cohort Extraction (When Ready)
```bash
export ENABLE_SAE_COHORT_RUN=1
export MAX_PATIENTS=66  # Or remove limit
export MAX_TOTAL_VARIANTS=3000

python3 scripts/sae/extract_sae_features_cohort.py
```

**Note**: This will take time and use Modal credits. Monitor costs.

---

### Step 5: Run Biomarker Analysis
```bash
python3 scripts/sae/analyze_biomarkers.py \
    --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json
```

**Expected**: Top 100 features with correlations and p-values

---

### Step 6: Create Pathway Mapping
```bash
python3 scripts/sae/create_feature_pathway_mapping.py \
    --biomarkers data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
    --cohort data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json \
    --top-n 100
```

**Expected**: `sae_feature_mapping.json` file created with feature→pathway mappings

---

## ✅ Verification Checklist

- [ ] Backend starts successfully
- [ ] `test_sae_extraction.py` passes all tests
- [ ] Provenance shows "trained weights" (not "random init")
- [ ] Dimension is 4096 (not 1920)
- [ ] Small batch test completes successfully
- [ ] Full cohort extraction completes (if running)
- [ ] Biomarker analysis finds significant features
- [ ] Pathway mapping file created

---

## 🐛 Troubleshooting

### Backend Not Reachable
- Check if backend is running: `curl http://localhost:8000/api/ayesha/complete_care_v2/health`
- Check if port 8000 is in use: `lsof -i :8000`

### Feature Flags Not Working
- Check `.env` file in `oncology-coPilot/oncology-backend-minimal/`
- Restart backend after changing `.env`
- Verify flags: `ENABLE_TRUE_SAE=1` and `ENABLE_EVO2_SAE=1`

### Modal Service Errors
- Check Modal service URLs are correct
- Verify services are deployed and running
- Check Modal dashboard for service status

### Dimension Mismatch
- Should be 4096 for evo2_7b
- If 1920, check model_id is "evo2_7b" not "evo2_1b"
- Verify SAE service is using evo2_7b_base model

---

## 📊 Expected Results

### Before (Random Weights)
- Model: "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)"
- d_in: 1920
- Features: May be noise

### After (Trained Weights)
- Model: "Goodfire/Evo-2-Layer-26-Mixed (trained weights)"
- d_in: 4096
- Features: Biologically meaningful

---

**Status**: ✅ **READY** - Start with Step 1!

