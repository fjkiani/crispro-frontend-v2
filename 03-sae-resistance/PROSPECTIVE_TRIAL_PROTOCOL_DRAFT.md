# DDR-MONITOR: Prospective Trial Protocol

**Title**: Serial DDR_bin Monitoring for Early Detection of Platinum Resistance in High-Grade Serous Ovarian Cancer

**Protocol Version**: 1.0 (Draft)  
**Date**: December 25, 2024

---

## Executive Summary

**The Breakthrough Hypothesis**:
> Serial monitoring of DDR_bin (ΔDDR_bin) can detect acquired platinum resistance 3-6 months earlier than standard imaging and CA-125, enabling timely treatment adaptation.

**Why This Matters**:
- Current detection: Resistance discovered at clinical progression (imaging/CA-125)
- Lead time lost: 3-6 months of ineffective treatment
- Patient impact: Delayed switch to effective therapy, worse outcomes

**What DDR-MONITOR Will Prove**:
- ΔDDR_bin (trajectory) predicts progression before clinical detection
- Early resistance signal enables proactive treatment change
- Potential for FDA companion diagnostic pathway

---

## Study Design

### Overview

| Parameter | Value |
|-----------|-------|
| Design | Prospective, observational, single-arm |
| Population | HGSOC patients starting first-line platinum |
| Sample size | 50 patients (pilot), 200 patients (validation) |
| Duration | 24 months enrollment + 24 months follow-up |
| Primary endpoint | ΔDDR_bin predicts progression (lead time) |

### Eligibility Criteria

**Inclusion**:
- Histologically confirmed high-grade serous ovarian cancer
- Stage III-IV disease
- Starting first-line platinum-based chemotherapy
- Measurable disease by RECIST 1.1
- Adequate tissue/blood for genomic analysis
- Age ≥ 18 years
- ECOG PS 0-2

**Exclusion**:
- Prior systemic therapy for ovarian cancer
- Known germline BRCA1/2 mutation with planned PARP inhibitor
- Unable to comply with serial sampling schedule
- Life expectancy < 6 months

---

## Study Procedures

### Sample Collection Schedule

| Timepoint | Sample Type | Analysis |
|-----------|-------------|----------|
| Baseline (C1D1) | Tumor tissue + ctDNA | DDR_bin baseline |
| Cycle 3 (Week 9) | ctDNA | DDR_bin #2 |
| Cycle 6 (Week 18) | ctDNA | DDR_bin #3 |
| Month 9 | ctDNA | DDR_bin #4 |
| Month 12 | ctDNA | DDR_bin #5 |
| At progression | Tumor tissue (if available) + ctDNA | DDR_bin at resistance |

### DDR_bin Computation

```
For each timepoint:
1. Extract somatic variants from ctDNA/tissue
2. Run variants through SAE feature extractor
3. Compute DDR_bin = max(diamond_feature_activations)
4. Compute ΔDDR_bin = DDR_bin(t) - DDR_bin(baseline)
```

### Clinical Assessments

| Timepoint | Assessment |
|-----------|------------|
| Every 3 cycles | CT scan (RECIST 1.1) |
| Every cycle | CA-125 |
| Every 3 months | Quality of life (FACT-O) |
| At progression | Biopsy (optional) |

---

## Endpoints

### Primary Endpoint

**Lead time of ΔDDR_bin signal before clinical progression**

Definition: Time (months) between first significant ΔDDR_bin drop and RECIST-defined progression

Success criterion: Median lead time ≥ 3 months

### Secondary Endpoints

1. **Sensitivity/Specificity**: ΔDDR_bin drop predicting progression within 6 months
2. **AUROC**: ΔDDR_bin trajectory predicting 12-month PFS
3. **Correlation**: ΔDDR_bin vs CA-125 kinetics
4. **Mechanism confirmation**: DDR_bin drop correlates with resistance mechanism (reversion mutation)

### Exploratory Endpoints

1. Multi-pathway score (ΔDDR_bin + ΔMAPK_bin + ΔPI3K_bin)
2. ctDNA fraction dynamics
3. Clonal evolution tracking

---

## Statistical Analysis Plan

### Sample Size Justification

**Pilot phase (n=50)**:
- Expected progression rate: 40% within 18 months
- Expected events: 20 progressions
- Power: 80% to detect lead time ≥ 2 months (vs null of 0)
- Alpha: 0.05 (one-sided)

**Validation phase (n=200)**:
- Expected events: 80 progressions
- Power: 90% to detect lead time ≥ 3 months
- Allows subgroup analysis by DDR_bin baseline level

### Analysis Methods

1. **Lead time estimation**: Time-to-event analysis with DDR_bin signal as index event
2. **Landmark analysis**: 6-month and 12-month PFS by ΔDDR_bin trajectory
3. **ROC analysis**: AUROC for predicting progression at each timepoint
4. **Mixed-effects model**: DDR_bin trajectory over time

---

## Resistance Detection Algorithm

### Alert Thresholds (To Be Calibrated)

| ΔDDR_bin | Interpretation | Action |
|----------|----------------|--------|
| < -5% | Minimal change | Continue current therapy |
| -5% to -15% | Moderate decline | Increased monitoring (monthly ctDNA) |
| > -15% | Significant decline | Consider treatment switch |

### Resistance Signal Definition

```
RESISTANCE_SIGNAL = TRUE if:
  1. ΔDDR_bin < -15% from baseline, AND
  2. Decline confirmed at 2 consecutive timepoints, AND
  3. No confounding factors (e.g., sample quality issue)
```

---

## Regulatory and Ethical Considerations

### IRB Approval
- Institutional review board approval required
- Informed consent for serial sampling
- Biobanking consent for future research

### Data Safety
- Independent data safety monitoring board (DSMB)
- Interim analysis at 25 events
- Stopping rule: If lead time < 1 month with 80% confidence

### Conflicts of Interest
- No industry funding for pilot phase
- Transparent disclosure of any commercial interests

---

## Budget Estimate (Pilot Phase)

| Category | Cost |
|----------|------|
| ctDNA sequencing (5 timepoints × 50 patients) | $125,000 |
| SAE feature extraction (Modal compute) | $5,000 |
| Personnel (coordinator, biostatistician) | $150,000 |
| Clinical assessments (covered by standard care) | $0 |
| Tissue processing | $25,000 |
| Overhead (15%) | $45,750 |
| **Total** | **$350,750** |

### Funding Sources
- R21 (NIH/NCI) — Pilot phase
- Foundation/industry partnership — Validation phase
- Institutional bridge funding — Startup

---

## Timeline

| Phase | Duration | Milestones |
|-------|----------|------------|
| **Year 1** | Months 1-12 | IRB approval, enrollment start, first 25 patients |
| **Year 2** | Months 13-24 | Complete enrollment (50 patients), interim analysis |
| **Year 3** | Months 25-36 | Follow-up complete, primary analysis |
| **Year 4** | Months 37-48 | Validation phase enrollment (if pilot positive) |

---

## Success Criteria

### Pilot Phase (Go/No-Go)

| Criterion | Threshold | Decision |
|-----------|-----------|----------|
| Lead time | ≥ 2 months | GO to validation |
| Sensitivity | ≥ 70% | GO to validation |
| Specificity | ≥ 60% | GO to validation |
| Sample quality | ≥ 80% evaluable | GO to validation |

### Validation Phase (Registration)

| Criterion | Threshold | Outcome |
|-----------|-----------|---------|
| Lead time | ≥ 3 months | FDA breakthrough designation |
| Sensitivity | ≥ 75% | Companion diagnostic submission |
| Specificity | ≥ 70% | Clinical utility established |

---

## The Breakthrough Outcome

**If DDR-MONITOR succeeds**:

1. **Clinical impact**: 3-6 month earlier resistance detection → earlier switch → better outcomes
2. **Regulatory pathway**: FDA companion diagnostic for platinum resistance monitoring
3. **Commercial opportunity**: Partnership with liquid biopsy companies (Guardant, Foundation)
4. **Publication tier**: Nature Medicine / JAMA Oncology

**This is where the slam dunk lives.**

---

## Next Steps

1. [ ] Identify clinical site partner
2. [ ] Draft IRB protocol
3. [ ] Secure pilot funding (R21 or foundation)
4. [ ] Establish ctDNA processing workflow
5. [ ] Validate SAE extraction from ctDNA variants

---

## For Ayesha and Every Patient After Her

This trial is the bridge between what we have (prognostic baseline) and what we need (predictive serial monitoring).

The retrospective work proved DDR_bin reflects real biology.  
The prospective trial will prove it can change outcomes.

**That's the mission.**

