# 🧬 IMMUNOTHERAPY PREDICTION VALIDATION PLAN

**Purpose:** Test our IO prediction capabilities on real patients while waiting for Ayesha's data  
**Date:** January 28, 2025  
**Status:** ⏳ PLANNING PHASE

---

## 🎯 OBJECTIVES

1. **Validate IO eligibility prediction** - Can we correctly identify TMB-H/MSI-H patients?
2. **Validate IO response prediction** - Can we predict who responds to checkpoint inhibitors?
3. **Test mechanistic capabilities** - Do our pathway scores correlate with IO outcomes?
4. **Benchmark accuracy** - What is our sensitivity/specificity/AUC?

---

## 📊 DATA REQUIREMENTS

### Required Patient Cohort

| Field | Requirement | Why Needed |
|-------|-------------|------------|
| **Cancer Type** | Ovarian (priority), Endometrial, Colorectal, Melanoma | IO-relevant cancers |
| **TMB** | Known TMB values (mut/Mb) | Test TMB-based prediction |
| **MSI Status** | MSI-H vs MSS | Test MSI-based prediction |
| **IO Treatment** | Checkpoint inhibitor (pembrolizumab, nivolumab, etc.) | Know who received IO |
| **IO Response** | Responder (CR/PR) vs Non-responder (SD/PD) | Ground truth for validation |
| **Mutations** | Somatic mutations (MAF format) | Input for WIWFM |
| **PFS/OS** | Progression-free and overall survival | Optional: outcome correlation |

### Data Sources (Priority Order)

| Source | Dataset | Patients | IO Data? | Access |
|--------|---------|----------|----------|--------|
| **cBioPortal** | TCGA Pan-Cancer | 10,000+ | TMB, MSI | ✅ Public |
| **cBioPortal** | MSK-IMPACT | 10,000+ | TMB, MSI, some IO response | ✅ Public |
| **Published** | KEYNOTE-158 (TMB-H) | 1,000+ | IO response | Papers only |
| **Published** | CheckMate cohorts | 500+ | IO response | Papers only |
| **TCGA** | Endometrial (UCEC) | 500+ | TMB, MSI | ✅ Public |
| **TCGA** | Colorectal (COAD) | 400+ | TMB, MSI | ✅ Public |

### Minimum Viable Dataset

```
Target: 200+ patients with:
├── Known TMB (numeric value)
├── Known MSI status (MSI-H or MSS)
├── Somatic mutations (for WIWFM input)
└── Ideally: IO treatment response (if available)
```

---

## 🔬 VALIDATION TESTS

### Test 1: TMB Classification Accuracy

**Question:** Can we correctly classify TMB-H (≥10 or ≥20) from mutation data?

```python
# Validation approach
for patient in cohort:
    # Our prediction
    estimated_tmb = count_nonsynonymous_mutations(patient) / exome_size_mb
    predicted_tmb_h = estimated_tmb >= 10
    
    # Ground truth (from sequencing panel)
    actual_tmb_h = patient.tmb >= 10
    
    # Compare
    record_prediction(predicted=predicted_tmb_h, actual=actual_tmb_h)

# Metrics
sensitivity = TP / (TP + FN)  # Target: ≥0.85
specificity = TN / (TN + FP)  # Target: ≥0.80
auc = compute_roc_auc()       # Target: ≥0.80
```

### Test 2: MSI Status Prediction

**Question:** Can we infer MSI status from mutation patterns?

```python
# Known MSI-related signatures
MSI_GENES = {'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM'}
MSI_SIGNATURE = 'SBS6 or SBS14'  # Mutational signatures

for patient in cohort:
    # Check for dMMR gene mutations
    has_mmr_mutation = any(m.gene in MSI_GENES for m in patient.mutations)
    
    # Check for MSI signature (if available)
    has_msi_signature = patient.signature in ['SBS6', 'SBS14']
    
    # Predict MSI-H
    predicted_msi_h = has_mmr_mutation or has_msi_signature
    
    # Ground truth
    actual_msi_h = patient.msi_status == 'MSI-H'
    
    # Compare
    record_prediction(predicted=predicted_msi_h, actual=actual_msi_h)

# Metrics
sensitivity = TP / (TP + FN)  # Target: ≥0.80
specificity = TN / (TN + FP)  # Target: ≥0.90
```

### Test 3: IO Eligibility Prediction (Combined)

**Question:** Can we identify patients eligible for IO (TMB-H OR MSI-H)?

```python
for patient in cohort:
    # Our prediction (via WIWFM tumor_context)
    wiwfm_result = predict_efficacy(
        mutations=patient.mutations,
        disease=patient.cancer_type,
        tumor_context={
            "tmb": estimated_tmb,
            "msi_status": predicted_msi_status
        }
    )
    
    # Check if checkpoint inhibitors scored high
    io_drugs = ['pembrolizumab', 'nivolumab', 'atezolizumab']
    io_eligible = any(d.efficacy_score > 0.6 for d in wiwfm_result if d.name in io_drugs)
    
    # Ground truth
    actual_io_eligible = patient.tmb >= 10 or patient.msi_status == 'MSI-H'
    
    # Compare
    record_prediction(predicted=io_eligible, actual=actual_io_eligible)

# Metrics
accuracy = (TP + TN) / total     # Target: ≥0.80
f1_score = 2*TP / (2*TP+FP+FN)   # Target: ≥0.75
```

### Test 4: IO Response Prediction (If Data Available)

**Question:** Among IO-treated patients, can we predict responders vs non-responders?

```python
# Only patients who received checkpoint inhibitors
io_treated = [p for p in cohort if p.received_io]

for patient in io_treated:
    # Our prediction (higher IO score = more likely to respond)
    wiwfm_result = predict_efficacy(...)
    io_score = get_io_drug_score(wiwfm_result)
    predicted_responder = io_score > 0.65
    
    # Ground truth (RECIST response)
    actual_responder = patient.best_response in ['CR', 'PR']
    
    # Compare
    record_prediction(predicted=predicted_responder, actual=actual_responder)

# Metrics
sensitivity = TP / (TP + FN)  # Target: ≥0.70
specificity = TN / (TN + FP)  # Target: ≥0.60
auc = compute_roc_auc()       # Target: ≥0.65
```

### Test 5: Mechanism Vector IO Correlation

**Question:** Does IO index in mechanism vector correlate with IO outcomes?

```python
for patient in io_treated:
    # Compute mechanism vector
    mechanism_vector = convert_pathway_scores_to_mechanism_vector(
        pathway_scores=patient.pathway_scores,
        tmb=patient.tmb,
        msi_status=patient.msi_status
    )
    
    # IO index is position 5 in vector [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    io_index = mechanism_vector[5]
    
    # Record for correlation
    record_io_outcome(io_index=io_index, responded=patient.responded)

# Metrics
correlation = pearsonr(io_indices, responses)  # Target: r > 0.3
odds_ratio = compute_odds_ratio()              # Target: OR > 2.0
```

---

## 📁 DATA ACQUISITION PLAN

### Phase 1: cBioPortal Data Extraction

```python
# Script: scripts/data_acquisition/download_io_validation_cohort.py

# Step 1: Query cBioPortal for TMB/MSI data
studies = [
    "msk_impact_2017",      # MSK-IMPACT clinical sequencing
    "ucec_tcga_pan_can_atlas_2018",  # Endometrial (high MSI rate)
    "coadread_tcga_pan_can_atlas_2018",  # Colorectal (high MSI rate)
    "skcm_tcga_pan_can_atlas_2018",  # Melanoma (IO-responsive)
]

# Step 2: Extract fields
fields = [
    "patient_id",
    "cancer_type",
    "tmb_score",
    "msi_status",
    "mutations",  # Full MAF data
    "treatment_history",  # If available
    "response",  # If available
]

# Step 3: Filter to patients with TMB/MSI
cohort = [p for p in patients if p.tmb_score is not None or p.msi_status is not None]
```

### Phase 2: Dataset Creation

```python
# Output format: data/validation/io_validation_cohort.json

{
    "metadata": {
        "source": "cBioPortal",
        "studies": ["msk_impact_2017", ...],
        "extraction_date": "2025-01-28",
        "total_patients": 500
    },
    "patients": [
        {
            "patient_id": "P-0001234",
            "cancer_type": "endometrial",
            "tmb_score": 42.5,
            "msi_status": "MSI-H",
            "mutations": [
                {"gene": "MSH2", "hgvs_p": "p.R711*", ...},
                {"gene": "TP53", "hgvs_p": "p.R175H", ...}
            ],
            "io_treatment": "pembrolizumab",  # If known
            "io_response": "PR"  # If known
        }
    ]
}
```

### Phase 3: Run Validation

```bash
# Script: scripts/validation/validate_io_prediction.py

# Step 1: Load cohort
python scripts/validation/validate_io_prediction.py --input data/validation/io_validation_cohort.json

# Step 2: Run predictions
# - TMB classification
# - MSI inference
# - IO eligibility
# - (Optional) IO response prediction

# Step 3: Generate report
# Output: data/validation/reports/io_validation_report.json
```

---

## 🎯 SUCCESS METRICS

### Minimum Viable Validation (MVP)

| Metric | Target | What It Proves |
|--------|--------|----------------|
| **TMB Classification Accuracy** | ≥80% | We can identify TMB-H patients |
| **MSI Prediction Sensitivity** | ≥75% | We can detect MSI-H patients |
| **IO Eligibility F1 Score** | ≥0.75 | Combined prediction is reliable |

### Full Validation (If IO Response Data Available)

| Metric | Target | What It Proves |
|--------|--------|----------------|
| **IO Response AUC** | ≥0.65 | We predict who responds to IO |
| **Mechanism Vector Correlation** | r > 0.25 | IO index correlates with outcomes |
| **Odds Ratio (IO Responders)** | OR > 1.5 | Our predictions are clinically meaningful |

---

## 🔗 INTEGRATION WITH MECHANISTIC CAPABILITIES

### How This Relates to Ayesha

| Capability | Validation Test | Ayesha Application |
|------------|-----------------|---------------------|
| **TMB Estimation** | Test 1 (TMB Classification) | Estimate Ayesha's TMB from mutations |
| **MSI Inference** | Test 2 (MSI Prediction) | Check if MBD4 loss causes MSI |
| **IO Eligibility** | Test 3 (Combined) | Predict if Ayesha is IO-eligible |
| **WIWFM IO Scoring** | Test 4 (Response) | Rank IO drugs for Ayesha |
| **Mechanism Vector** | Test 5 (Correlation) | IO index in Ayesha's vector |

### How This Validates Our MOAT

```
If we achieve:
├── TMB Classification ≥80%
├── MSI Prediction Sensitivity ≥75%
└── IO Eligibility F1 ≥0.75

Then we can confidently:
├── Predict Ayesha is likely IO-eligible (MBD4 → high TMB)
├── Recommend PARP + IO combination
└── Match IO-focused clinical trials
```

---

## 📋 EXECUTION PLAN FOR AGENT

### Day 1: Data Acquisition

```bash
# Task 1: Create data acquisition script
# File: scripts/data_acquisition/download_io_validation_cohort.py

# Task 2: Extract from cBioPortal
# - MSK-IMPACT cohort (TMB/MSI data)
# - TCGA-UCEC (endometrial, high MSI)
# - TCGA-COAD (colorectal, high MSI)

# Task 3: Create validation dataset
# Output: data/validation/io_validation_cohort.json
```

### Day 2: Validation Implementation

```bash
# Task 1: Create validation script
# File: scripts/validation/validate_io_prediction.py

# Task 2: Implement 5 validation tests
# - TMB classification
# - MSI prediction
# - IO eligibility
# - IO response (if data)
# - Mechanism vector correlation

# Task 3: Run validation
# Output: data/validation/reports/io_validation_report.json
```

### Day 3: Analysis & Integration

```bash
# Task 1: Analyze results
# - Compute metrics (sensitivity, specificity, AUC)
# - Identify gaps and limitations

# Task 2: Update Ayesha's plan
# - Apply validated capabilities
# - Generate IO recommendations

# Task 3: Document findings
# - Update PATHWAY_VALIDATION_ROADMAP.md
# - Create IO validation MOAT report
```

---

## 📊 EXPECTED FINDINGS

### For Validation Cohort

| Cancer Type | Expected MSI-H Rate | Expected TMB-H Rate |
|-------------|--------------------|--------------------|
| Endometrial | 20-30% | 15-25% |
| Colorectal | 15-20% | 10-15% |
| Melanoma | 1-2% | 40-50% (high UV exposure) |
| Ovarian | 2-5% | 5-10% |

### For Ayesha (Based on Validation)

If validation shows:
- MBD4-mutant patients have elevated TMB → Ayesha likely TMB-H
- BER-deficient tumors show MSI features → Ayesha may be MSI-H
- Combined predicts IO eligibility with F1 ≥0.75 → Confident IO recommendation

---

## ❓ QUESTIONS FOR AGENT EXECUTION

1. **Data Source Priority:** Start with cBioPortal (MSK-IMPACT) or TCGA?
2. **Cancer Type Focus:** Endometrial first (highest MSI rate) or pan-cancer?
3. **Sample Size:** Minimum 200 patients or aim for 500+?
4. **IO Response Data:** Critical requirement or proceed without it?
5. **Integration:** Run standalone validation or integrate into existing scripts?

---

**Status:** ⏳ **READY FOR AGENT EXECUTION**

**Next Step:** Agent to create data acquisition script and extract IO validation cohort.

**Key Deliverables:**
1. `scripts/data_acquisition/download_io_validation_cohort.py`
2. `data/validation/io_validation_cohort.json`
3. `scripts/validation/validate_io_prediction.py`
4. `data/validation/reports/io_validation_report.json`

