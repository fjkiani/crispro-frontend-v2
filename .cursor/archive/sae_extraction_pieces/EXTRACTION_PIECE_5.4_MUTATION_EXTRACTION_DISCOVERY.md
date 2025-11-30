# EXTRACTION PIECE 5.4: Mutation Extraction Discovery

**Source**: Lines 12000-12200 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-XX  
**Status**: ✅ Complete

---

## Overview

This piece documents the critical discovery that the TCGA-OV platinum labels file contained **only patient IDs and response labels**, but **ZERO mutation data**. This led to the realization that mutation extraction from pyBioPortal was a required prerequisite step.

---

## The Discovery

**User Question**: "did we even extract the dataset that we needed or you just skipped it and hallucinated"

**Investigation**: Checked `data/validation/sae_cohort/` directory

**Files Found**:
- `sae_features_tcga_ov_platinum_MOCK.json` (mock data only)
- `sae_tcga_ov_platinum_biomarkers_MOCK.json` (mock analysis only)
- Logs and plots for mock run

**Critical Finding**: **NO real cohort file** (`sae_features_tcga_ov_platinum.json`)

---

## The Problem

**Platinum Labels File Structure** (`data/validation/tcga_ov_platinum_response_labels.json`):

```json
{
  "tcga_patient_id": "TCGA-23-2078",
  "tcga_sample_id": "TCGA-23-2078-01",
  "platinum_response": "sensitive",
  "raw_response_value": "Complete Remission/Response",
  "source_field": "2df4e5a5-f831-4f0a-81b7-35f516ef42ec"
}
```

**Keys Available**: `['tcga_patient_id', 'tcga_sample_id', 'platinum_response', 'raw_response_value', 'source_field']`

**Missing**: **ZERO mutation data** (no `chrom`, `pos`, `ref`, `alt`)

---

## The Impact

**Original Assumption**: `extract_sae_features_cohort.py` would read mutations from `patient.get("mutations", [])`

**Reality**: That key doesn't exist. The script would process **zero mutations** and write an **empty cohort**.

**What Was Actually Done**:
- ✅ Built + verified the **mock** pipeline end-to-end (stats + plots + RUO endpoint)
- ✅ Planned and partially wired the **real** path (using platinum labels + pyBioPortal mutation extraction + Evo2/SAE services)
- ❌ **NOT executed** the real extraction to completion

---

## The Solution

**Required Steps**:

1. **Fetch somatic mutations** for each of the 469 TCGA-OV patients using their `tcga_sample_id` (e.g., `TCGA-23-2078-01`) via cBioPortal/pyBioPortal

2. **Merge** those mutations into the platinum labels file (or keep them separate and join at extraction time)

3. **Then** run `extract_sae_features_cohort.py` with real mutation data to call Evo2→SAE and produce the real cohort SAE file

---

## pyBioPortal Pattern

**Existing Script**: `oncology-coPilot/oncology-backend-minimal/scripts/tcga_extraction/extract_mutation_frequencies.py`

**Proven Pattern**:
```python
df_muts = mut.get_muts_in_mol_prof_by_sample_list_id(
    profile_id,    # e.g., "ov_tcga_pan_can_atlas_2018_mutations"
    sample_list_id, # e.g., "ov_tcga_pan_can_atlas_2018_all"
    projection="DETAILED",
    pageSize=10000
)
```

**Returns**: DataFrame with columns:
- `sampleId`
- `gene`
- `chr`
- `startPosition`
- `referenceAllele`
- `variantAllele`

---

## Targeted Mutation Extractor Plan

**Script to Build**: Mutation extractor that:

1. Reads the 469 labeled patients from `tcga_ov_platinum_response_labels.json`

2. Uses pyBioPortal to fetch **per-sample mutations** for `ov_tcga_pan_can_atlas_2018` (filtering to our 469 sample IDs)

3. Merges mutations back into a new file: `data/validation/sae_cohort/tcga_ov_platinum_with_mutations.json`

4. Then `extract_sae_features_cohort.py` can read that file and call Evo2→SAE for each variant

---

## Honest Assessment

**What Was Done**:
- ✅ Mock pipeline verified end-to-end
- ✅ Real path planned and partially wired
- ❌ Real extraction **NOT executed** to completion

**What Was NOT Done**:
- ❌ Real cohort SAE extraction
- ❌ Real biomarker analysis on real data
- ❌ Mutation extraction from pyBioPortal

**Next Concrete Steps**:
```bash
# Step 1: Extract mutations from pyBioPortal
python3 scripts/sae/extract_tcga_mutations.py

# Step 2: Extract SAE features for cohort
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
ENABLE_EVO2_SAE=1 ENABLE_TRUE_SAE=1 \
python3 scripts/sae/extract_sae_features_cohort.py

# Step 3: Analyze biomarkers
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots_real
```

---

## Key Lessons

1. **Always verify data structure**: Don't assume keys exist without checking
2. **Mock data is useful**: Proves pipeline works before expensive real extraction
3. **Honest communication**: Admit when steps haven't been executed
4. **Clear next steps**: Provide concrete commands for what needs to happen next

---

## Related Documents

- `scripts/sae/extract_sae_features_cohort.py` - Cohort extraction script (needs mutation data)
- `oncology-coPilot/oncology-backend-minimal/scripts/tcga_extraction/extract_mutation_frequencies.py` - pyBioPortal pattern
- `data/validation/tcga_ov_platinum_response_labels.json` - Labels file (no mutations)

