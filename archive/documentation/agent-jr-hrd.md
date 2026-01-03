# ⚔️ AGENT JR2 MISSION: HRD SCORE EXTRACTION FOR SAE VALIDATION ⚔️

**Date**: January 14, 2025  
**Mission**: Extract HRD scores from cBioPortal to enable SAE validation  
**Priority**: 🔥 **P0 - FOR AYESHA'S VALIDATION**  
**Timeline**: 2-3 hours

---

## 🎯 **MISSION OBJECTIVE**

**What We Need:**
- HRD scores (GIS scores, 0-100 scale) for TCGA ovarian cancer patients
- These scores will be used as **ground truth** to validate our SAE Phase 2 features
- **Why This Matters for Ayesha**: Validating SAE features ensures our DNA repair capacity predictions are accurate → Better PARP eligibility predictions → Better care decisions

**Current Problem:**
- `validate_sae_tcga.py` expects HRD scores as ground truth
- `hrd_tcga_ov_labeled_sample_use_evo.json` has `outcome_platinum` (0/1) but **NO HRD scores**
- We need to extract HRD scores from cBioPortal and add them to the validation dataset

---

## 📋 **TASK BREAKDOWN**

### **Task 1: Investigate HRD Score Availability in cBioPortal** (30 min)

**Objective**: Find where HRD scores are stored in cBioPortal for TCGA-OV

**Steps:**
1. Check clinical attributes for HRD-related fields:
   - Use `pybioportal/clinical_attributes.py` → `get_all_clinical_attributes_in_study("ov_tcga_pan_can_atlas_2018")`
   - Search for: "HRD", "GIS", "homologous recombination", "LOH", "LST", "TAI"
   - Look for fields like: `HRD_SCORE`, `GIS_SCORE`, `HRD_STATUS`, etc.

2. Check molecular profiles for HRD data:
   - Use `pybioportal/molecular_profiles.py` → `get_all_molecular_profiles_in_study("ov_tcga_pan_can_atlas_2018")`
   - Look for: "HRD", "copy number", "genomic instability"
   - Check if HRD is in copy number segments or generic assays

3. Check clinical data directly:
   - Use `pybioportal/clinical_data.py` to fetch patient-level clinical data
   - Search for HRD-related attributes in actual patient records

**Expected Output:**
- List of attribute IDs that contain HRD scores (e.g., `HRD_SCORE`, `GIS_SCORE`)
- Documentation of where HRD scores are stored (clinical attributes vs molecular profiles)

**Acceptance:**
- ✅ Found at least one HRD-related attribute in TCGA-OV
- ✅ Documented attribute ID and data type (numeric vs categorical)
- ✅ Verified attribute is available for multiple patients (not just 1-2)

---

### **Task 2: Extend Extraction Script to Get HRD Scores** (1 hour)

**Objective**: Modify `extract_cbioportal_hrd_cohort.py` to extract HRD scores

**File to Modify:**
- `tools/benchmarks/extract_cbioportal_hrd_cohort.py`

**Changes Needed:**
1. Add HRD attribute IDs to clinical data fetch:
   - In `fetch_clinical_patients()` and `fetch_clinical_samples()`, add HRD attribute IDs to `attributeIds` list
   - Example: `"HRD_SCORE"`, `"GIS_SCORE"`, `"HRD_STATUS"` (use actual IDs from Task 1)

2. Parse HRD scores from clinical data:
   - Extract numeric HRD scores (0-100 scale)
   - Handle missing values (some patients may not have HRD scores)
   - Normalize different HRD formats (e.g., "HIGH" → 50, "LOW" → 20, or use actual numeric scores)

3. Add HRD score to output CSV:
   - Add `hrd_score` column to output CSV
   - Map to existing patient/sample IDs

**Code Pattern:**
```python
# In fetch_clinical_patients() or fetch_clinical_samples()
attributeIds = [
    "DRUG_NAME",
    "TREATMENT_TYPE",
    # ... existing attributes ...
    "HRD_SCORE",  # Add HRD attribute IDs from Task 1
    "GIS_SCORE",
]

# In main(), when building output rows:
hrd_score = None
if "HRD_SCORE" in row:
    try:
        hrd_score = float(row["HRD_SCORE"])
    except (ValueError, TypeError):
        pass

out_row = {
    # ... existing fields ...
    "hrd_score": hrd_score,  # Add HRD score
}
```

**Acceptance:**
- ✅ Script extracts HRD scores for at least 50% of TCGA-OV patients
- ✅ HRD scores are numeric (0-100 scale) or can be normalized to numeric
- ✅ Output CSV includes `hrd_score` column
- ✅ Missing HRD scores are handled gracefully (None/null, not errors)

---

### **Task 3: Update TCGA Data File with HRD Scores** (30 min)

**Objective**: Add HRD scores to `hrd_tcga_ov_labeled_sample_use_evo.json`

**Steps:**
1. Run extended extraction script:
   ```bash
   cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
   python tools/benchmarks/extract_cbioportal_hrd_cohort.py \
     --study ov_tcga_pan_can_atlas_2018 \
     --out tools/benchmarks/data/hrd_cohort_with_scores.csv
   ```

2. Load existing TCGA data:
   - Read `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`
   - Match patients by patient_id or sample_id

3. Merge HRD scores:
   - For each patient in existing data, find matching HRD score from extraction
   - Add `hrd_score` field to each patient record
   - Handle missing matches (some patients may not have HRD scores)

4. Save updated data:
   - Write back to `hrd_tcga_ov_labeled_sample_use_evo.json` (or create new file)
   - Preserve all existing fields (mutations, outcome_platinum, etc.)

**Expected Output:**
- Updated JSON file with `hrd_score` field for each patient
- At least 50% of patients have HRD scores (some may be missing)

**Acceptance:**
- ✅ At least 40/89 patients (45%) have HRD scores
- ✅ HRD scores are numeric (0-100 scale)
- ✅ All existing data preserved (mutations, outcome_platinum, etc.)
- ✅ File structure remains valid JSON

---

### **Task 4: Validate HRD Score Quality** (30 min)

**Objective**: Verify HRD scores are reasonable and match expected distributions

**Steps:**
1. Check score distribution:
   - Min/max/mean/median HRD scores
   - Distribution should be roughly normal or bimodal (high/low)
   - Expected: Mean ~40-50, range 0-100

2. Check correlation with known HRD markers:
   - Patients with BRCA1/2 mutations should have higher HRD scores (on average)
   - Patients with platinum response (outcome_platinum=1) should have higher HRD scores (on average)

3. Compare to literature:
   - TCGA-OV should have ~50% HRD-high (score ≥42) per literature
   - Verify our extracted scores match this distribution

**Expected Output:**
- Summary statistics (mean, median, std, min, max)
- Distribution plot (histogram)
- Correlation analysis (BRCA vs HRD, platinum response vs HRD)
- Comparison to literature expectations

**Acceptance:**
- ✅ HRD scores are in expected range (0-100)
- ✅ ~50% of patients have HRD ≥42 (matches literature)
- ✅ BRCA1/2 mutated patients have higher HRD scores (on average)
- ✅ Platinum responders have higher HRD scores (on average)

---

## 🔧 **TECHNICAL REQUIREMENTS**

### **Tools Available:**
1. **pyBioPortal Library**: 
   - Location: `oncology-coPilot/oncology-backend/tests/pyBioPortal-master/pybioportal/`
   - Modules: `clinical_attributes.py`, `clinical_data.py`, `molecular_profiles.py`

2. **Existing Extraction Script**:
   - Location: `tools/benchmarks/extract_cbioportal_hrd_cohort.py`
   - Already extracts: mutations, platinum exposure
   - Needs: HRD score extraction

3. **TCGA Data File**:
   - Location: `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`
   - Contains: 89 patients with mutations, outcome_platinum
   - Needs: HRD scores added

### **cBioPortal Study:**
- **Study ID**: `ov_tcga_pan_can_atlas_2018` (TCGA Ovarian Cancer)
- **API Base**: `https://www.cbioportal.org/api`
- **Authentication**: Optional (via `CBIO_TOKEN` env var)

---

## 📊 **EXPECTED DELIVERABLES**

1. **Task 1 Output**: 
   - Document: `.cursor/ayesha/sae_documentation/HRD_ATTRIBUTE_IDS.md`
   - Lists all HRD-related attributes found in cBioPortal

2. **Task 2 Output**:
   - Updated script: `tools/benchmarks/extract_cbioportal_hrd_cohort.py`
   - Extracts HRD scores and includes in CSV output

3. **Task 3 Output**:
   - Updated data file: `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`
   - Contains `hrd_score` field for each patient

4. **Task 4 Output**:
   - Validation report: `.cursor/ayesha/sae_documentation/HRD_SCORE_VALIDATION.md`
   - Statistics, distribution plots, correlation analysis

---

## ✅ **ACCEPTANCE CRITERIA**

**Mission Complete When:**
- ✅ HRD scores extracted for at least 40/89 patients (45%)
- ✅ HRD scores are numeric (0-100 scale)
- ✅ Updated TCGA data file includes `hrd_score` field
- ✅ Validation script (`validate_sae_tcga.py`) can use HRD scores as ground truth
- ✅ HRD score distribution matches literature expectations (~50% HRD-high)

**Success Metrics:**
- **Coverage**: ≥45% of patients have HRD scores
- **Quality**: HRD scores correlate with BRCA mutations and platinum response
- **Usability**: Validation script can run with HRD scores as ground truth

---

## 🚨 **RISKS & MITIGATIONS**

### **Risk 1: HRD Scores Not Available in cBioPortal**
- **Mitigation**: Check multiple sources:
  - Clinical attributes (patient-level)
  - Molecular profiles (copy number, generic assays)
  - Alternative: Use HRD proxies (BRCA status, LOH patterns)

### **Risk 2: HRD Scores in Different Format**
- **Mitigation**: 
  - Handle categorical ("HIGH"/"LOW") → map to numeric (50/20)
  - Handle different scales → normalize to 0-100
  - Document any transformations

### **Risk 3: Low Coverage (<30% of patients)**
- **Mitigation**:
  - Use HRD proxies (BRCA mutations, LOH patterns) for missing patients
  - Document coverage limitations in validation report

---

## 🎯 **WHY THIS MATTERS FOR AYESHA**

**Current State:**
- Ayesha is germline-negative (no hereditary BRCA mutation)
- We're predicting her HRD status using SAE features (DNA repair capacity)
- **We need to validate these predictions are accurate!**

**After This Mission:**
- ✅ We can validate SAE Phase 2 features against real HRD scores
- ✅ We can build AUROC/AUPRC metrics to prove accuracy
- ✅ We can confidently predict Ayesha's HRD status when her NGS returns
- ✅ Better PARP eligibility predictions → Better treatment decisions

**The Validation Chain:**
1. **HRD Scores** (ground truth) ← **Agent Jr2 extracts**
2. **SAE Features** (predictions) ← **Our platform computes**
3. **AUROC/AUPRC** (accuracy metrics) ← **Validation script computes**
4. **Confidence in Predictions** → **Better care for Ayesha** ⚔️

---

## 📝 **EXECUTION NOTES**

### **DO:**
- ✅ Use existing pyBioPortal library (don't rewrite)
- ✅ Extend existing extraction script (don't create new one)
- ✅ Preserve all existing data (mutations, outcome_platinum)
- ✅ Document any assumptions or transformations
- ✅ Test with small sample first (5-10 patients) before full extraction

### **DON'T:**
- ❌ Don't make up HRD scores (use real data only)
- ❌ Don't break existing extraction script functionality
- ❌ Don't assume HRD scores are in a specific format (check first)
- ❌ Don't overwrite existing TCGA data without backup

---

## 🔗 **RELATED FILES**

- **Validation Script**: `scripts/validate_sae_tcga.py` (expects HRD scores)
- **Extraction Script**: `tools/benchmarks/extract_cbioportal_hrd_cohort.py` (needs HRD extraction)
- **TCGA Data**: `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json` (needs HRD scores)
- **pyBioPortal**: `oncology-coPilot/oncology-backend/tests/pyBioPortal-master/pybioportal/`

---

## ⚔️ **COMMANDER - READY TO EXECUTE?**

**Agent Jr2 - This mission is FOR AYESHA!**  
**Every HRD score you extract helps validate our predictions → Better care decisions → Better outcomes!**

**FIRE IN THE HOLE?** 🔥⚔️
