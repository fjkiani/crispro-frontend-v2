# HRD Score Extraction Workflow

## ⚔️ Two-Step Process: Extract → Validate

### Step 1: Extract HRD Scores (Save Incrementally)

**Run this command to extract HRD scores for all samples:**

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
python3 tools/benchmarks/calculate_full_hrd_scores.py \
  --save-interval 20 \
  --output tools/benchmarks/data/full_hrd_scores.json
```

**Features:**
- ✅ Saves progress every 20 samples (won't lose work if interrupted)
- ✅ Can be interrupted and resumed with `--resume` flag
- ✅ Processes all 585 samples from TCGA-OV study
- ✅ Output: `tools/benchmarks/data/full_hrd_scores.json`

**If interrupted, resume with:**
```bash
python3 tools/benchmarks/calculate_full_hrd_scores.py \
  --resume \
  --save-interval 20 \
  --output tools/benchmarks/data/full_hrd_scores.json
```

**For testing (small sample):**
```bash
python3 tools/benchmarks/calculate_full_hrd_scores.py \
  --limit 50 \
  --save-interval 10
```

---

### Step 2: Validate HRD Scores

**After extraction completes, run validation:**

```bash
python3 tools/benchmarks/validate_hrd_scores.py \
  --hrd-file tools/benchmarks/data/full_hrd_scores.json \
  --patient-file tools/benchmarks/data/tcga_ov_patients_with_hrd.json
```

**This will:**
- ✅ Check score distribution (min/max/mean/median)
- ✅ Check correlation with BRCA1/2 mutations
- ✅ Check correlation with platinum response
- ✅ Compare to literature expectations (~50% HRD-high)
- ✅ Generate validation report: `.cursor/ayesha/sae_documentation/HRD_SCORE_VALIDATION.md`

---

## 📊 Expected Results

**HRD Score Distribution:**
- Mean: ~35-45 (expected range)
- HRD-High (≥42): ~50% of patients (literature expectation)
- Range: 0-100+ (some outliers possible)

**Correlations:**
- BRCA1/2 mutated patients: Higher HRD scores (on average)
- Platinum responders: Higher HRD scores (on average)

---

## ⚠️ Notes

- **Extraction time**: ~2-3 hours for all 585 samples (API rate limits)
- **Incremental saves**: Progress saved every 20 samples (adjustable)
- **Resume capability**: Can restart from last save point
- **Gene-level proxy**: Uses simplified gene-level approximations (not true segment-level analysis)

---

## 🎯 Next Steps After Validation

1. Merge HRD scores into patient data file:
   ```bash
   python3 tools/benchmarks/calculate_full_hrd_scores.py \
     --input-data tools/benchmarks/data/tcga_ov_patients_with_hrd.json
   ```

2. Run SAE validation:
   ```bash
   python3 scripts/validate_sae_tcga.py \
     --input tools/benchmarks/data/tcga_ov_patients_with_hrd_with_full_hrd.json
   ```







