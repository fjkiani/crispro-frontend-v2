# Next Phase Ready - Biomarker Analysis & Pathway Mapping

**Date**: January 20, 2025  
**Status**: ⏸️ **AWAITING EXTRACTION COMPLETION**

---

## 🎯 Ready to Execute (After Extraction Completes)

### Phase 2: Biomarker Discovery Analysis

**Script**: `scripts/sae/analyze_biomarkers.py`

**Command**:
```bash
python3 scripts/sae/analyze_biomarkers.py \
    --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json
```

**Expected Output**:
- Top 100 features with correlations
- Statistical significance (p-values, FDR-corrected)
- Effect sizes (Cohen's d)
- Visualization plots

**Success Criteria**:
- ≥10 significant features (p < 0.05, FDR-corrected)
- Correlations |r| ≥ 0.3
- Effect sizes Cohen's d ≥ 0.5

---

### Phase 3: Feature→Pathway Mapping

**Script**: `scripts/sae/create_feature_pathway_mapping.py`

**Command**:
```bash
python3 scripts/sae/create_feature_pathway_mapping.py \
    --biomarkers data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
    --cohort data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
    --output oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json \
    --top-n 100 \
    --min-correlation 0.3 \
    --min-p-value 0.05
```

**Expected Output**:
- `sae_feature_mapping.json` file
- Feature indices mapped to pathways (DDR, MAPK, PI3K, etc.)
- Correlation statistics per pathway
- Confidence levels

**Validation**:
- Test on known cases (BRCA1→DDR, KRAS→MAPK)
- Verify ≥80% accuracy on known cases

---

## 📋 Execution Sequence

1. **Wait for extraction completion** (66 patients)
2. **Verify extraction quality** (provenance, dimensions)
3. **Run biomarker analysis** (discover predictive features)
4. **Create pathway mapping** (map features to biology)
5. **Validate mapping** (test on known cases)
6. **Document results** (update status, create summary)

---

## ⚠️ Potential Issues & Mitigations

### If Biomarker Analysis Finds 0 Features
- **Mitigation**: Relax thresholds (p < 0.1, |r| ≥ 0.2)
- **Check**: Sample size, outcome distribution
- **Alternative**: Stratified analysis (sensitive vs resistant)

### If Mapping Validation Fails
- **Mitigation**: Review gene→pathway assignments
- **Check**: Feature activation patterns
- **Iterate**: Refine mapping strategy

---

**Status**: ✅ **SCRIPTS READY** - Execute after extraction completes



