# 📖 EXTRACTION PIECE 3.2: Real Data Extraction
**Source**: Lines 11000-11400 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents the planning and approach for extracting real SAE features from TCGA-OV cohort, including pyBioPortal integration for mutation extraction and wiring into the cohort extraction pipeline.

---

## 🔍 KEY FINDINGS

### **Real Cohort File (RUO) and Data Status**

**Target File:**
- `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`
- Non-mock version meant to hold per-patient SAE features for TCGA-OV platinum response cohort

**What We Already Have:**

1. **Labels/Source Cohort:**
   - `data/validation/tcga_ov_platinum_response_labels.json` (built by platinum_hunt scripts)
   - Contains: TCGA patient/sample IDs, platinum response class (`sensitive` / `resistant` / `refractory`)

2. **Variant Sources:**
   - cBioPortal / GDC extractors that can pull somatic mutations per TCGA case
   - Code + wiring ready, used elsewhere in project

3. **Extraction & Analysis Code:**
   - `scripts/sae/extract_sae_features_cohort.py`: Orchestrates Evo2 activations + SAE feature extraction
   - `api/services/biomarker_correlation_service.py` + `scripts/sae/analyze_biomarkers.py`: Consume cohort file and run correlation analysis

**What's Missing:**
- Real SAE cohort file hasn't been materialized yet (only MOCK version exists)
- Need to pull real somatic mutations for each labeled TCGA-OV patient
- Need to wire mutation extractors into `extract_sae_features_cohort.py`

---

### **How We've Already Used pyBioPortal**

#### **1. For Platinum Labels** (`scripts/platinum_hunt/services/pybioportal_treatments_extractor.py`)

**Approach:**
- Adds `pyBioPortal-master` to `sys.path`, imports `pybioportal.treatments`
- Calls `treatments.fetch_all_patient_level_treatments(study_view_filter={"studyIds": ["ov_tcga_pan_can_atlas_2018", "ov_tcga"]})`
- Parses treatment rows per `patientId`, looks for response/outcome columns and platinum drug names
- Normalizes text responses via `normalize_response()` and `extract_tcga_patient_id()` into `{"platinum_response": "sensitive|resistant|refractory", ...}`
- Feeds into `scripts/platinum_hunt/orchestrator.py`, which merges GDC XML, TCGA-CDR, pyBioPortal, Broad Firehose → `data/validation/tcga_ov_platinum_response_labels.json`

#### **2. For Mutation Frequencies** (`oncology-backend-minimal/scripts/tcga_extraction/extract_mutation_frequencies.py`)

**Approach:**
- Adds `pyBioPortal-master` to `sys.path`, imports `pybioportal.molecular_profiles`, `sample_lists`, `mutations`
- For each study (including `"ov_tcga_pan_can_atlas_2018"`):
  - Finds the `*_mutations` molecular profile and an `*_all` sample list
  - Calls `mut.get_muts_in_mol_prof_by_sample_list_id(profile_id, sample_list_id, projection="DETAILED", pageSize=10000)`
  - Filters to a list of target genes, then computes **per-pathway mutation frequencies** (how many samples have any mutation in that gene set)
- Output: `mutation_frequencies.json` for universal disease database (pathway weights), **not** per-patient mutation lists

---

### **Can We Reuse This for Real Cohort?**

**Yes—everything we need is in place:**

**What We Know How To Do:**
- Talk to the same cBioPortal instance via `pybioportal` (see both scripts)
- Resolve the **TCGA OV study ID** (`"ov_tcga_pan_can_atlas_2018"`)
- Grab the correct mutation profile + sample list
- Parse the returned mutation dataframe

**What's Missing:**
- Different aggregation: Instead of collapsing to pathway frequencies, we need to keep the **per-sample/per-patient mutation records** and feed them into `extract_sae_features_cohort.py`

---

### **Concrete Plan: How to Get the Real Cohort File**

#### **Step 1: Add Mutation Extractor for TCGA-OV**

**Inside `extract_sae_features_cohort.py` or a sibling script:**

**Imports (same as `extract_mutation_frequencies.py`):**
```python
from pathlib import Path
import sys
import pandas as pd
pybioportal_parent = Path("oncology-coPilot/oncology-backend/tests/pyBioPortal-master")
if str(pybioportal_parent) not in sys.path:
    sys.path.insert(0, str(pybioportal_parent))
from pybioportal import molecular_profiles as mp, sample_lists as sl, mutations as mut
```

**Profile/Sample-List Logic:**
```python
study_id = "ov_tcga_pan_can_atlas_2018"
df_prof = mp.get_all_molecular_profiles_in_study(study_id)
profile_id = next(pid for pid in df_prof["molecularProfileId"].astype(str) if pid.endswith("_mutations"))

df_lists = sl.get_all_sample_lists_in_study(study_id)
sample_list_id = next(sid for sid in df_lists["sampleListId"].astype(str) if sid.endswith("_all"))
df_muts = mut.get_muts_in_mol_prof_by_sample_list_id(
    profile_id, sample_list_id, projection="DETAILED", pageSize=10000
)
```

**Build Mutation Map:**
- Build dict: `sampleId -> [mutation rows]`
- Build dict: `tcga_patient_id -> sampleIds` using `extract_tcga_patient_id` helper
- For each mutation row, map to expected fields:
  ```python
  {
    "gene": row[gene_col],                 # Hugo_Symbol / hugoGeneSymbol
    "chrom": row["chromosome"],
    "pos": int(row["startPosition"]),
    "ref": row["referenceAllele"],
    "alt": row["tumorSeqAllele2"],        # or equivalent
    "hgvs_p": row.get("proteinChange")    # if present
  }
  ```

#### **Step 2: Wire Into `extract_sae_features_cohort.py`**

**Instead of reading mutations from labels JSON (currently missing):**
- Load `tcga_ov_platinum_response_labels.json` (for `tcga_patient_id`, `platinum_response`, `tcga_sample_id`)
- For each patient, look up their sample ID in the mutation map built above
- Attach the list of mutation dicts as `patient["mutations"]`
- Proceed exactly as script already does:
  - For each mutation: call `/api/evo/score_variant_with_activations`
  - Feed `layer_26_alt_activations` into `/api/sae/extract_features`
  - Aggregate SAE features per patient
  - Write to: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

#### **Step 3: Run Extraction (RUO)**

**From repo root, with Evo2 + SAE endpoints enabled:**

```bash
export ENABLE_EVO2_SAE=true
export ENABLE_TRUE_SAE=true
python3 scripts/sae/extract_sae_features_cohort.py
```

**Output:** Real, non-mock cohort file: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

#### **Step 4: Run Analysis on Real Data**

**Once file exists:**

```bash
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots_real
```

**Same `BiomarkerCorrelationService` will give real, RUO SAE↔platinum-response biomarkers.**

---

## 📊 KEY INSIGHTS

### **Reusability**

- **pyBioPortal patterns are reusable**: Same library, same study ID, same mutation profile logic
- **Only difference**: Aggregation level (per-patient mutations vs pathway frequencies)
- **Integration**: Can be added directly to existing extraction script

### **Data Flow**

1. **Labels** → `tcga_ov_platinum_response_labels.json` (already exists)
2. **Mutations** → pyBioPortal extraction → per-patient mutation lists
3. **SAE Features** → Evo2 activations + SAE decoding → per-patient SAE features
4. **Analysis** → Biomarker correlation → top features

### **Dependencies**

- **pyBioPortal**: Already integrated, just need to use mutations endpoint
- **Evo2 Service**: Must be deployed with activations endpoint
- **SAE Service**: Must be deployed and accessible
- **Backend**: Must have flags enabled (`ENABLE_EVO2_SAE`, `ENABLE_TRUE_SAE`)

---

## 🔗 CONTEXT & CONNECTIONS

- **Builds on**: Mock data testing (Piece 3.1), Deployment instructions (Piece 2.5)
- **Enables**: Real biomarker discovery with actual patient data
- **Depends on**: Modal services deployment, pyBioPortal mutation extraction
- **Key Insight**: Reusing existing pyBioPortal patterns makes integration straightforward

---

## 📝 NOTES

- Real cohort file is the non-mock version of the mock file
- pyBioPortal integration follows existing patterns
- Mutation extraction is the missing piece
- Once mutations are extracted, existing pipeline handles SAE extraction
- Analysis pipeline already tested with mock data

---

## 🎯 QUESTIONS RESOLVED

- ✅ Do we have data sources? → Yes, labels exist, mutations can be extracted via pyBioPortal
- ✅ Can we reuse existing code? → Yes, pyBioPortal patterns are reusable
- ✅ What's missing? → Mutation extraction integration into cohort script
- ✅ How to proceed? → Add mutation extractor, wire into extraction script, run extraction

