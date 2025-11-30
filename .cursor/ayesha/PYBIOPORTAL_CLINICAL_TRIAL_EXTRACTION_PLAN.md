# pyBioPortal Clinical Trial Dataset Extraction Plan

**Date**: January 27, 2025  
**Purpose**: Extract datasets from cBioPortal for Phase 4 Clinical Trial Outcome Benchmarking  
**Library**: `pybioportal` (Python client for cBioPortal API)

---

## 🎯 What is pyBioPortal?

**pyBioPortal** is a Python client library for accessing the **cBioPortal API**. cBioPortal is a cancer genomics data portal containing:

- **Patient Data**: Real patient cohorts from TCGA, MSK-IMPACT, and other studies
- **Mutations**: Somatic mutations (MAF format)
- **Clinical Data**: Demographics, outcomes (PFS, OS, response rates)
- **Treatments**: Patient-level treatment data (drugs, response)
- **Studies**: Metadata about cancer studies

**Base URL**: `https://www.cbioportal.org/api`

---

## 📊 What We Can Extract for Clinical Trial Benchmarking

### 1. **Patient Mutations** (For Our System Input)

**Endpoint**: `mutations.fetch_muts_in_mol_prof()`

**What We Get**:
- Gene names (e.g., BRCA1, TP53, MBD4)
- Protein changes (e.g., p.Arg175His)
- Chromosome, position, ref/alt alleles
- Sample IDs (to link with clinical data)

**Example**:
```python
from pybioportal import mutations as mut

# Get mutations for a study
mutations_df = mut.fetch_muts_in_mol_prof(
    molecular_profile_id="ov_tcga_pan_can_atlas_2018_mutations",
    sample_list_id="ov_tcga_pan_can_atlas_2018_all"
)
```

### 2. **Clinical Outcomes** (For Ground Truth)

**Endpoint**: `clinical_data.get_all_clinical_data_in_study()`

**What We Get**:
- **Response Rates**: ORR (objective response rate), CR (complete response), PR (partial response)
- **Survival Data**: PFS (progression-free survival), OS (overall survival)
- **Biomarkers**: HRD status, BRCA status, TP53 status
- **Demographics**: Age, stage, histology

**Example**:
```python
from pybioportal import clinical_data as cd

# Get clinical data for a study
clinical_df = cd.get_all_clinical_data_in_study(
    study_id="ov_tcga_pan_can_atlas_2018",
    clinical_data_type="PATIENT"
)
```

### 3. **Treatment Data** (For Drug Matching)

**Endpoint**: `treatments.fetch_all_patient_level_treatments()`

**What We Get**:
- **Drug Names**: Olaparib, carboplatin, paclitaxel, etc.
- **Treatment Lines**: First-line, second-line, maintenance
- **Response**: Response to treatment (if available)
- **Dates**: Treatment start/end dates

**Example**:
```python
from pybioportal import treatments as trt

# Get treatment data
treatments_df = trt.fetch_all_patient_level_treatments(
    study_view_filter={"studyIds": ["ov_tcga_pan_can_atlas_2018"]},
    tier="Agent"
)
```

### 4. **Study Metadata** (For Trial Matching)

**Endpoint**: `studies.get_study()`

**What We Get**:
- **Study Name**: e.g., "TCGA Ovarian Cancer (PanCancer Atlas)"
- **Cancer Type**: Ovarian cancer
- **Sample Count**: Number of patients
- **Description**: Study details, inclusion criteria

**Example**:
```python
from pybioportal import studies as st

# Get study info
study_df = st.get_study("ov_tcga_pan_can_atlas_2018")
```

---

## 🔍 Target Studies for Clinical Trial Benchmarking

### Ovarian Cancer Studies (Primary Focus)

1. **`ov_tcga_pan_can_atlas_2018`** (TCGA Ovarian Cancer)
   - **Size**: ~300-400 patients
   - **Data Available**: Mutations, clinical outcomes, treatments
   - **Use Case**: General ovarian cancer cohort

2. **`ov_tcga`** (TCGA Ovarian Cancer - Original)
   - **Size**: ~300 patients
   - **Data Available**: Mutations, clinical outcomes
   - **Use Case**: Alternative TCGA dataset

3. **Other Ovarian Studies** (if available)
   - Search cBioPortal for additional ovarian cancer studies
   - Filter by: Has mutations + Has clinical outcomes + Has treatments

### How to Find Studies

```python
from pybioportal import studies as st

# Search for ovarian cancer studies
ovarian_studies = st.get_all_studies(keyword="ovarian")
print(ovarian_studies[['studyId', 'name', 'cancerTypeId', 'allSampleCount']])
```

---

## 🛠️ Extraction Strategy for Phase 4 Benchmarking

### Step 1: Identify Target Studies

**Goal**: Find studies that match our target clinical trials (SOLO-2, PAOLA-1, PRIMA, etc.)

**Approach**:
1. Search cBioPortal for ovarian cancer studies
2. Filter by: Has mutation data + Has clinical outcomes + Has treatment data
3. Map to target trials (if study matches trial population)

**Script**:
```python
def find_target_studies():
    """Find cBioPortal studies that match our target clinical trials."""
    from pybioportal import studies as st
    
    # Search for ovarian cancer studies
    all_studies = st.get_all_studies(keyword="ovarian")
    
    # Filter by criteria
    target_studies = []
    for study_id in all_studies['studyId']:
        study_info = st.get_study(study_id)
        # Check if has mutations, clinical data, treatments
        # Add to target_studies if matches
    
    return target_studies
```

### Step 2: Extract Patient Mutations

**Goal**: Get mutations for each patient in target studies

**Approach**:
1. Get molecular profile ID (mutations profile)
2. Get sample list ID (all samples)
3. Fetch mutations for all samples
4. Map to patient IDs

**Script**:
```python
def extract_patient_mutations(study_id: str) -> Dict[str, List[Dict]]:
    """Extract mutations for each patient in a study."""
    from pybioportal import molecular_profiles as mp
    from pybioportal import sample_lists as sl
    from pybioportal import mutations as mut
    
    # Get mutation profile
    profiles = mp.get_all_molecular_profiles_in_study(study_id)
    mut_profile_id = profiles[profiles['molecularProfileId'].str.endswith('_mutations')]['molecularProfileId'].iloc[0]
    
    # Get sample list
    sample_lists = sl.get_all_sample_lists_in_study(study_id)
    sample_list_id = sample_lists[sample_lists['sampleListId'].str.endswith('_all')]['sampleListId'].iloc[0]
    
    # Fetch mutations
    mutations_df = mut.fetch_muts_in_mol_prof(
        molecular_profile_id=mut_profile_id,
        sample_list_id=sample_list_id
    )
    
    # Group by patient
    mutations_by_patient = {}
    for _, row in mutations_df.iterrows():
        patient_id = row['patientId']
        if patient_id not in mutations_by_patient:
            mutations_by_patient[patient_id] = []
        mutations_by_patient[patient_id].append(row.to_dict())
    
    return mutations_by_patient
```

### Step 3: Extract Clinical Outcomes

**Goal**: Get response rates, PFS, OS for each patient

**Approach**:
1. Get all clinical data (PATIENT level)
2. Extract outcome fields: ORR, PFS, OS
3. Extract biomarker fields: HRD, BRCA, TP53
4. Map to patient IDs

**Script**:
```python
def extract_clinical_outcomes(study_id: str) -> pd.DataFrame:
    """Extract clinical outcomes for patients in a study."""
    from pybioportal import clinical_data as cd
    
    # Get all clinical data
    clinical_df = cd.get_all_clinical_data_in_study(
        study_id=study_id,
        clinical_data_type="PATIENT"
    )
    
    # Pivot to wide format (one row per patient)
    clinical_wide = clinical_df.pivot(
        index='patientId',
        columns='clinicalAttributeId',
        values='value'
    )
    
    # Extract key outcomes
    outcomes = {
        'patient_id': clinical_wide.index,
        'response_rate': clinical_wide.get('ORR', None),
        'pfs_months': clinical_wide.get('PFS_MONTHS', None),
        'os_months': clinical_wide.get('OS_MONTHS', None),
        'hrd_status': clinical_wide.get('HRD_STATUS', None),
        'brca_status': clinical_wide.get('BRCA_STATUS', None),
        'tp53_status': clinical_wide.get('TP53_STATUS', None),
    }
    
    return pd.DataFrame(outcomes)
```

### Step 4: Extract Treatment Data

**Goal**: Get treatment history and response for each patient

**Approach**:
1. Get patient-level treatments
2. Extract drug names, treatment lines, response
3. Map to patient IDs

**Script**:
```python
def extract_treatment_data(study_id: str) -> pd.DataFrame:
    """Extract treatment data for patients in a study."""
    from pybioportal import treatments as trt
    
    # Get patient-level treatments
    treatments_df = trt.fetch_all_patient_level_treatments(
        study_view_filter={"studyIds": [study_id]},
        tier="Agent"
    )
    
    # Extract key fields
    treatment_data = treatments_df[[
        'patientId',
        'treatment',
        'startDate',
        'endDate',
        'response',
        'treatmentType'
    ]]
    
    return treatment_data
```

### Step 5: Combine Data for Benchmarking

**Goal**: Create unified dataset with mutations + outcomes + treatments

**Approach**:
1. Join mutations, clinical outcomes, treatments by patient ID
2. Filter to patients with complete data (mutations + outcomes)
3. Structure for benchmarking

**Script**:
```python
def create_benchmark_dataset(study_id: str) -> pd.DataFrame:
    """Create unified dataset for clinical trial outcome benchmarking."""
    
    # Extract all data
    mutations = extract_patient_mutations(study_id)
    outcomes = extract_clinical_outcomes(study_id)
    treatments = extract_treatment_data(study_id)
    
    # Combine
    benchmark_data = []
    for patient_id in outcomes['patient_id']:
        patient_mutations = mutations.get(patient_id, [])
        patient_outcomes = outcomes[outcomes['patient_id'] == patient_id].iloc[0]
        patient_treatments = treatments[treatments['patientId'] == patient_id]
        
        # Only include if has mutations and outcomes
        if patient_mutations and not patient_outcomes.isna().all():
            benchmark_data.append({
                'patient_id': patient_id,
                'study_id': study_id,
                'mutations': patient_mutations,
                'response_rate': patient_outcomes['response_rate'],
                'pfs_months': patient_outcomes['pfs_months'],
                'os_months': patient_outcomes['os_months'],
                'hrd_status': patient_outcomes['hrd_status'],
                'treatments': patient_treatments.to_dict('records')
            })
    
    return pd.DataFrame(benchmark_data)
```

---

## 🎯 Integration with Phase 4 Benchmarking

### How This Fits into Phase 4

**Current Phase 4 Plan**:
- Extract published RCT results (SOLO-2, PAOLA-1, etc.)
- Run our system on trial patient populations
- Compare predictions vs. published outcomes

**With pyBioPortal**:
- **Alternative Approach**: Use cBioPortal studies as proxy for clinical trials
- **Advantage**: Real patient-level data (not just aggregate trial results)
- **Use Case**: Can validate individual patient predictions, not just aggregate correlations

### Two Benchmarking Approaches

#### Approach A: Published RCT Results (Original Plan)
- **Data Source**: Published papers (SOLO-2, PAOLA-1, etc.)
- **Ground Truth**: Aggregate response rates (65% ORR for olaparib)
- **Our Prediction**: Efficacy scores for trial population
- **Metric**: Correlation between efficacy scores and response rates
- **Limitation**: Only aggregate data, not individual patient outcomes

#### Approach B: cBioPortal Patient Data (NEW - Using pyBioPortal)
- **Data Source**: cBioPortal studies (TCGA-OV, etc.)
- **Ground Truth**: Individual patient outcomes (response, PFS, OS)
- **Our Prediction**: Efficacy scores for each patient
- **Metric**: 
  - Individual patient correlation (efficacy vs. response)
  - Aggregate correlation (mean efficacy vs. mean response)
  - Survival analysis (efficacy vs. PFS/OS)
- **Advantage**: Patient-level validation, not just aggregate

### Recommended: Hybrid Approach

**Use Both**:
1. **Published RCT Results**: Validate aggregate predictions (correlation with trial response rates)
2. **cBioPortal Patient Data**: Validate individual patient predictions (efficacy vs. individual outcomes)

**Why Both**:
- RCT results: Validates against published clinical trial data (gold standard)
- cBioPortal data: Validates individual patient predictions (more granular)

---

## 📋 Implementation Plan

### Phase 4A: Extract cBioPortal Datasets (2-3 days)

**Tasks**:
1. **Create extraction script** (`scripts/benchmark/extract_cbioportal_trial_datasets.py`)
   - Find target studies (ovarian cancer)
   - Extract mutations, clinical outcomes, treatments
   - Combine into unified dataset
   - Save to `data/benchmarks/cbioportal_trial_datasets.json`

2. **Create dataset validator** (`scripts/benchmark/validate_cbioportal_datasets.py`)
   - Check data completeness (mutations, outcomes, treatments)
   - Filter to patients with complete data
   - Report data quality metrics

3. **Test extraction** (1 day)
   - Run on `ov_tcga_pan_can_atlas_2018`
   - Verify data structure
   - Check data quality

### Phase 4B: Run Benchmarking (2-3 days)

**Tasks**:
1. **Create benchmark script** (`scripts/benchmark/benchmark_clinical_trial_outcomes_cbioportal.py`)
   - Load cBioPortal datasets
   - For each patient: Run `/api/efficacy/predict`
   - Compare predictions vs. actual outcomes
   - Compute correlation metrics

2. **Add correlation analysis**
   - Pearson correlation (efficacy vs. response)
   - Spearman correlation (efficacy vs. PFS/OS)
   - Survival analysis (efficacy vs. PFS/OS using Kaplan-Meier)

3. **Integrate with existing benchmarks**
   - Add to benchmark suite
   - Standardize output format
   - Add to CI/CD pipeline

### Phase 4C: Combine with Published RCT Results (1 day)

**Tasks**:
1. **Run both benchmarks** (cBioPortal + Published RCT)
2. **Compare results** (which approach gives better correlation?)
3. **Document findings** (which is more reliable?)

---

## 🔧 Technical Implementation

### File Structure

```
scripts/benchmark/
├── extract_cbioportal_trial_datasets.py      # Extract datasets from cBioPortal
├── validate_cbioportal_datasets.py           # Validate extracted data
├── benchmark_clinical_trial_outcomes_cbioportal.py  # Run benchmarking
└── utils/
    ├── cbioportal_client.py                  # pyBioPortal wrapper
    ├── dataset_combiner.py                  # Combine mutations + outcomes + treatments
    └── correlation_analysis.py             # Correlation computation

data/benchmarks/
├── cbioportal_trial_datasets.json            # Extracted datasets
└── cbioportal_benchmark_results.json         # Benchmark results
```

### Dependencies

**Required**:
- `pybioportal` (already in repo: `oncology-coPilot/oncology-backend/tests/pyBioPortal-master`)
- `pandas` (for data manipulation)
- `httpx` (for API calls to our system)
- `scipy` (for correlation analysis)
- `lifelines` (for survival analysis, optional)

**Installation**:
```bash
# pybioportal is already in repo, just add to path
# No additional installation needed
```

### Example Usage

```python
# Extract dataset
python scripts/benchmark/extract_cbioportal_trial_datasets.py \
    --study ov_tcga_pan_can_atlas_2018 \
    --output data/benchmarks/cbioportal_ov_tcga.json

# Run benchmark
python scripts/benchmark/benchmark_clinical_trial_outcomes_cbioportal.py \
    --dataset data/benchmarks/cbioportal_ov_tcga.json \
    --api_base http://127.0.0.1:8000 \
    --output results/benchmarks/cbioportal_benchmark_20250127.json
```

---

## 📊 Expected Output Format

### Extracted Dataset Format

```json
{
  "study_id": "ov_tcga_pan_can_atlas_2018",
  "extraction_date": "2025-01-27T10:00:00Z",
  "patients": [
    {
      "patient_id": "TCGA-04-1357",
      "mutations": [
        {
          "gene": "TP53",
          "protein_change": "p.Arg175His",
          "chromosome": "17",
          "position": 7577120,
          "ref": "G",
          "alt": "A"
        },
        {
          "gene": "BRCA1",
          "protein_change": "p.Gln1756Ter",
          "chromosome": "17",
          "position": 43094051,
          "ref": "C",
          "alt": "T"
        }
      ],
      "clinical_outcomes": {
        "response_rate": 0.65,
        "pfs_months": 13.8,
        "os_months": 37.2,
        "hrd_status": "positive",
        "brca_status": "BRCA1_mutant",
        "tp53_status": "mutant"
      },
      "treatments": [
        {
          "drug": "olaparib",
          "treatment_line": "maintenance",
          "response": "partial_response",
          "start_date": "2018-01-15",
          "end_date": "2019-07-20"
        }
      ]
    }
  ],
  "summary": {
    "total_patients": 300,
    "patients_with_mutations": 280,
    "patients_with_outcomes": 250,
    "patients_with_treatments": 200,
    "patients_complete": 180
  }
}
```

### Benchmark Results Format

```json
{
  "timestamp": "2025-01-27T10:00:00Z",
  "benchmark_type": "clinical_trial_outcomes_cbioportal",
  "study_id": "ov_tcga_pan_can_atlas_2018",
  "metrics": {
    "response_rate_correlation": {
      "pearson_r": 0.68,
      "spearman_rho": 0.65,
      "p_value": 0.001,
      "n_patients": 180
    },
    "pfs_correlation": {
      "pearson_r": 0.55,
      "spearman_rho": 0.52,
      "p_value": 0.005,
      "n_patients": 150
    },
    "os_correlation": {
      "pearson_r": 0.48,
      "spearman_rho": 0.45,
      "p_value": 0.01,
      "n_patients": 150
    },
    "drug_ranking_accuracy": {
      "top_1_accuracy": 0.45,
      "top_3_accuracy": 0.72,
      "n_patients": 180
    }
  },
  "targets": {
    "response_rate_correlation": {"min_pearson_r": 0.60},
    "pfs_correlation": {"min_pearson_r": 0.50},
    "drug_ranking_accuracy": {"min_top_3": 0.70}
  },
  "results": [
    {
      "patient_id": "TCGA-04-1357",
      "mutations": ["TP53 p.Arg175His", "BRCA1 p.Gln1756Ter"],
      "actual_response": 0.65,
      "our_efficacy_score": 0.82,
      "our_rank": 1,
      "actual_treatment": "olaparib",
      "correlation_status": "✅ PASS"
    }
  ],
  "provenance": {
    "script": "benchmark_clinical_trial_outcomes_cbioportal.py",
    "dataset_source": "cbioportal",
    "study_id": "ov_tcga_pan_can_atlas_2018",
    "api_base": "http://127.0.0.1:8000",
    "model_id": "evo2_1b"
  }
}
```

---

## ⚠️ Limitations & Considerations

### 1. **Study Population vs. Clinical Trial Population**

**Challenge**: cBioPortal studies (e.g., TCGA-OV) may not exactly match clinical trial populations (e.g., SOLO-2)

**Differences**:
- **TCGA-OV**: Observational study, various treatments, real-world data
- **SOLO-2**: Randomized controlled trial, specific inclusion criteria, olaparib maintenance

**Solution**: 
- Use TCGA-OV as proxy for general ovarian cancer population
- Filter TCGA-OV patients to match trial criteria (if possible)
- Acknowledge limitations in benchmark results

### 2. **Treatment Response Data Availability**

**Challenge**: cBioPortal may not have detailed treatment response data for all patients

**Reality**:
- Some studies have treatment data
- Some studies have response data
- Not all studies have both

**Solution**:
- Prioritize studies with complete treatment + response data
- Use clinical outcomes (PFS, OS) as proxy for treatment response
- Document data completeness in benchmark results

### 3. **Biomarker Data Availability**

**Challenge**: cBioPortal may not have HRD/BRCA status for all patients

**Reality**:
- Some studies have biomarker data
- Some studies infer biomarkers from mutations
- Not all studies have explicit biomarker fields

**Solution**:
- Infer HRD/BRCA status from mutations (BRCA1/BRCA2 mutations → HRD+)
- Use mutation data as primary source
- Document biomarker inference in benchmark results

### 4. **Data Quality**

**Challenge**: cBioPortal data may have inconsistencies or missing values

**Reality**:
- Some patients missing mutations
- Some patients missing outcomes
- Some patients missing treatments

**Solution**:
- Filter to patients with complete data (mutations + outcomes)
- Report data completeness metrics
- Handle missing values appropriately

---

## 🚀 Next Steps

### Immediate Actions

1. **Test pyBioPortal Access** (30 minutes)
   ```python
   # Quick test to verify pyBioPortal works
   from pybioportal import studies as st
   studies = st.get_all_studies(keyword="ovarian")
   print(studies[['studyId', 'name', 'allSampleCount']])
   ```

2. **Extract Sample Dataset** (2-3 hours)
   - Extract `ov_tcga_pan_can_atlas_2018` dataset
   - Verify data structure
   - Check data completeness

3. **Create Extraction Script** (1 day)
   - Implement `extract_cbioportal_trial_datasets.py`
   - Test on multiple studies
   - Validate output format

4. **Create Benchmark Script** (1-2 days)
   - Implement `benchmark_clinical_trial_outcomes_cbioportal.py`
   - Run on extracted dataset
   - Compute correlation metrics

5. **Integrate with Phase 4** (1 day)
   - Combine with published RCT results
   - Compare both approaches
   - Document findings

---

## 📚 References

- **cBioPortal API Docs**: https://docs.cbioportal.org/web-api-and-clients/
- **pyBioPortal Location**: `oncology-coPilot/oncology-backend/tests/pyBioPortal-master/pybioportal`
- **Existing Usage**: `scripts/sae/extract_sae_features_cohort.py` (lines 96-198)
- **Treatment Extraction**: `scripts/platinum_hunt/services/pybioportal_treatments_extractor.py`

---

**Status**: 📋 Planning Complete - Ready for Implementation  
**Estimated Time**: 6-8 days total (extraction: 2-3 days, benchmarking: 2-3 days, integration: 1-2 days)

