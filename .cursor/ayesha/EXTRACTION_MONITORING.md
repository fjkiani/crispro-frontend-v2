# Full Cohort Extraction - Monitoring

**Started**: January 20, 2025, 9:16 PM  
**Status**: 🔄 **RUNNING**

---

## 📊 Progress Tracking

**Target**: 66 patients, ~3000 variants  
**Process ID**: Check with `ps aux | grep extract_sae_features_cohort`

**Check Progress**:
```bash
# Check checkpoint
python3 -c "
import json
checkpoint = json.load(open('data/validation/sae_cohort/sae_features_extraction_checkpoint.json'))
print(f'Completed: {len(checkpoint.get(\"completed_patients\", []))}/66')
print(f'Failed: {len(checkpoint.get(\"failed_patients\", []))}')
"

# Check if process still running
ps aux | grep extract_sae_features_cohort | grep -v grep
```

---

## ⏱️ Estimated Timeline

- **Per patient**: ~1-2 minutes (depending on variant count)
- **Total time**: ~66-132 minutes (1-2 hours)
- **Checkpoints**: Saved every 10 patients

---

## ✅ Completion Verification

**When extraction completes, verify**:
1. Output file exists: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`
2. Patient count: 66 patients
3. Provenance: All show `d_in: 4096` and "trained weights"
4. No errors in final summary

**Verification command**:
```bash
python3 -c "
import json
data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json'))
print(f'Patients: {data.get(\"num_patients\")}')
print(f'Failed: {data.get(\"num_failed\", 0)}')
prov = data['patients'][0]['variants'][0].get('provenance', {})
print(f'Sample d_in: {prov.get(\"d_in\")}')
print(f'Sample model: {prov.get(\"model\")}')
"
```

---

## 🚀 Next Steps (After Completion)

1. **Verify extraction quality**
2. **Run biomarker analysis**: `python3 scripts/sae/analyze_biomarkers.py`
3. **Create pathway mapping**: `python3 scripts/sae/create_feature_pathway_mapping.py`
4. **Validate results**

---

**Status**: 🔄 **MONITORING** - Extraction in progress



