# Full Cohort Extraction - IN PROGRESS

**Date**: January 20, 2025  
**Status**: 🔄 **EXTRACTION RUNNING**

---

## 🚀 Extraction Started

**Command**:
```bash
export ENABLE_SAE_COHORT_RUN=1
export MAX_PATIENTS=66
export MAX_TOTAL_VARIANTS=3000
python3 scripts/sae/extract_sae_features_cohort.py
```

**Target**:
- 66 patients
- ~3000 variants (max)
- Trained SAE weights (evo2_7b)

---

## 📊 Expected Progress

**Estimated Time**: 10-20 minutes for full extraction

**Checkpoints**:
- Checkpoint saved every 10 patients
- Progress logged to console
- Output file: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

---

## ✅ Verification After Completion

**Check**:
1. Output file exists and was recently modified
2. Provenance shows trained weights for all patients
3. `d_in: 4096` for all variants
4. Patient count matches target (66 patients)

**Command to verify**:
```bash
python3 -c "
import json
from pathlib import Path
data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json'))
print(f'Patients: {data.get(\"num_patients\")}')
prov = data['patients'][0]['variants'][0].get('provenance', {})
print(f'd_in: {prov.get(\"d_in\")}')
print(f'Model: {prov.get(\"model\")}')
"
```

---

## 📝 Next Steps (After Extraction)

1. Verify extraction quality
2. Run biomarker analysis
3. Create pathway mapping
4. Validate results

---

**Status**: 🔄 **RUNNING** - Monitor progress and verify upon completion



