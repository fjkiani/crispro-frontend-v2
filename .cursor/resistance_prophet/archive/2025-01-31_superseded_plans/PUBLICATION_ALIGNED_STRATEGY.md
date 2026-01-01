# ⚔️ PUBLICATION-ALIGNED SURROGATE STRATEGY ⚔️

**Date:** January 31, 2025  
**Commander:** Alpha  
**Architect:** Zo  
**Status:** HONEST ASSESSMENT  
**Mission:** Align data, validate everything, publish receipts

---

## 🎯 THE REALITY CHECK

Alpha, you're right. This isn't a science project. This is a publication pipeline. Every claim needs validation. Every model needs data. Every "innovation" needs receipts.

**The Rule:** If we can't validate it with real data, we don't claim it.

---

## 📊 WHAT WE ACTUALLY HAVE (Honest Inventory)

### Data Acquisition Infrastructure (PRODUCTION READY ✅)

| **Client** | **Status** | **Location** | **What It Does** | **Validated?** |
|------------|------------|--------------|------------------|----------------|
| **cBioPortal** | ✅ Integrated | `scripts/data_acquisition/utils/cbioportal_client.py` | Studies, mutations, clinical data | ✅ Used for MM + OV |
| **Project Data Sphere** | ✅ Connected | `scripts/data_acquisition/utils/project_data_sphere_client.py` | 102 caslibs, patient-level trial data | ⚠️ Connected, not extracted |
| **ClinicalTrials.gov** | ✅ Integrated | `api/services/ctgov_query_builder.py` | Trial search, PI extraction | ✅ 962 trials loaded |
| **PubMed** | ✅ Integrated | `api/services/research_intelligence/` | Literature search, researcher ID | ✅ Production |
| **GDC** | ⚠️ Skeleton | `scripts/data_acquisition/utils/gdc_client.py` | TCGA data | ⚠️ Needs implementation |

### Publication-Ready Validated Work

| **Publication** | **Location** | **Validation Status** | **Data Source** | **N Patients** |
|-----------------|--------------|----------------------|-----------------|----------------|
| **Metastasis Interception** | `publication/` | ✅ 100% READY | AlphaFold 3, Cancer Genome Atlas | 304 gene-step combinations |
| **MM Drug Efficacy** | `oncology-coPilot/oncology-backend-minimal/` | ✅ 100% READY | MMRF CoMMpass | 995 patients |
| **Mechanism Trial Matching** | `.cursor/MOAT/CLINICAL_TRIALS/` | ✅ Validated | 962 trials | 47 MoA-tagged |
| **Resistance Prediction (OV)** | `.cursor/MOAT/ADVANCED_CARE_PLAN_RESISTANCE_PREDICTION.md` | ✅ Validated | TCGA-OV | 469 patients |
| **Resistance Prediction (MM)** | `.cursor/MOAT/ADVANCED_CARE_PLAN_RESISTANCE_PREDICTION.md` | ✅ Validated | MMRF CoMMpass | 995 patients |
| **Synthetic Lethality** | `.cursor/ayesha/SYNTHETIC_LETHALITY_COMPLETE.md` | ⚠️ Pilot Only | 10 cases | 50% accuracy (real ML) |
| **Sporadic Cancer** | `.cursor/MOAT/SPORADIC_CANCER_PRODUCTION_PLAN.md` | ⚠️ Unit Tests Only | TCGA priors | 6 unit tests |
| **Dosing Guidance** | `.cursor/plans/DOSING_GUIDANCE_VALIDATION_PLAN.md` | ❌ No Validation | Need extraction | 0 patients |

### What's Actually Validated (With Real Data)

| **Marker** | **Cancer** | **Relative Risk** | **p-value** | **N** | **Source** |
|------------|------------|------------------|-------------|-------|------------|
| **DIS3 mutation** | MM | 2.08 | 0.0145 | 38/219 | MMRF CoMMpass |
| **TP53 mutation** | MM | 1.90 | 0.11 (trend) | 16/219 | MMRF CoMMpass |
| **MAPK pathway** | OV | 1.97 | <0.05 | 35/469 | TCGA-OV |
| **NF1 mutation** | OV | 2.10 | <0.05 | 26/469 | TCGA-OV |
| **PI3K pathway** | OV | 1.39 | 0.02 | 108/469 | TCGA-OV |

---

## 📊 WHAT WE NEED (Honest Gaps)

### For Surrogate Endpoint Empire (KELIM, ECW/TBW, etc.)

| **Biomarker** | **Status** | **Data Needed** | **Source** | **Effort** |
|---------------|------------|-----------------|------------|------------|
| **CA-125 KELIM** | ❌ No Data | Serial CA-125 measurements | Collaborator or prospective | HIGH - need partnership |
| **ECW/TBW Surrogate** | ⚠️ Proxiable | BMI, albumin, age | TCGA-OV | LOW - can calculate |
| **Skeletal Muscle Mass** | ⚠️ CT-Based | CT scans at L3 | TCIA (TCGA imaging) | MEDIUM - need segmentation |
| **ctDNA MRD** | ❌ No Data | Serial ctDNA measurements | Collaborator or prospective | HIGH - need partnership |
| **TIL Density** | ❌ No Data | Pathology slides | TCGA or collaborator | HIGH - need analysis |

### For Dosing Guidance Validation

| **Source** | **Status** | **Expected N** | **Effort** |
|------------|------------|----------------|------------|
| **PubMed Case Reports** | ❌ Not Extracted | 20-30 cases | 4-6 hours (script ready) |
| **cBioPortal (MSK-IMPACT)** | ❌ Not Filtered | 30-50 cases | 3-4 hours (script ready) |
| **TCGA Germline** | ❌ Not Extracted | 15-20 cases | 6 hours (script ready) |
| **PharmGKB** | ❌ Not Queried | 15-20 cases | 3 hours (manual) |

### For Sporadic Cancer Validation

| **Component** | **Status** | **Validation Needed** |
|---------------|------------|-----------------------|
| **PARP Penalty (HRD<42)** | ✅ Unit Test | ❌ Need retrospective patient outcomes |
| **HRD Rescue (HRD≥42)** | ✅ Unit Test | ❌ Need retrospective patient outcomes |
| **TMB Boost (≥20)** | ✅ Unit Test | ❌ Need IO response correlation |
| **MSI-H Boost** | ✅ Unit Test | ❌ Need IO response correlation |
| **Quick Intake Priors** | ✅ 15 Cancer Types | ⚠️ Priors from literature, not patient outcomes |

---

## 🔗 HOW IT ALL CONNECTS (The Publication Roadmap)

### Tier 1: READY NOW (Submit Within Weeks)

| **Publication** | **Status** | **Target Journal** | **Data Ready?** |
|-----------------|------------|-------------------|-----------------|
| **Metastasis Interception** | ✅ 100% | Nature Biotechnology | ✅ All data |
| **MM Drug Efficacy** | ✅ 100% | npj Precision Oncology | ✅ All data |

### Tier 2: NEAR-READY (Need Validation Run)

| **Publication** | **Status** | **Gap** | **Fix** | **Time** |
|-----------------|------------|---------|---------|----------|
| **Resistance Prediction** | 85% | Need manuscript | Write it | 1 week |
| **Mechanism Trial Matching** | 70% | Need full manuscript | Write + run scripts | 2 weeks |
| **Sporadic Cancer** | 60% | Need patient validation | Extract outcomes from TCGA | 2 weeks |

### Tier 3: IN PROGRESS (Need Data)

| **Publication** | **Status** | **Gap** | **Data Source** | **Time** |
|-----------------|------------|---------|-----------------|----------|
| **Dosing Guidance** | 40% | Need validation cohort | PubMed + cBioPortal + TCGA | 2-3 weeks |
| **Synthetic Lethality** | 50% | Need more test cases | Expand pilot to 100 cases | 2 weeks |
| **Toxicity MOAT** | 60% | Need clinical outcomes | Collaborator or prospective | 4+ weeks |

### Tier 4: FUTURE (Need Partnerships)

| **Publication** | **Concept** | **Data Needed** | **Timeline** |
|-----------------|-------------|-----------------|--------------|
| **KELIM Resurrection** | CA-125 kinetics predicts resistance | Serial CA-125 from collaborator | 3+ months |
| **ECW/TBW Surrogate** | Body composition predicts resistance | TCGA + CT segmentation | 1-2 months |
| **ctDNA MRD Validation** | Early recurrence detection | ctDNA from collaborator | 6+ months |
| **Surrogate Endpoint Platform** | Multi-biomarker validation | Multiple data sources | 6-12 months |

---

## 📈 THE SURROGATE ENDPOINT EMPIRE - REALISTIC PATHWAY

### Phase 1: ECW/TBW Surrogate Proof (2-4 Weeks)

**Why This First:**
1. ✅ TCGA-OV has BMI, albumin, age data
2. ✅ We have cBioPortal client ready
3. ✅ Katsura 2023 validated the ECW/TBW correlation (N=320)
4. ❌ No need for collaborator partnership

**Deliverable:** "ECW/TBW surrogate (BMI/albumin/age) + genomics predicts platinum resistance AUROC X.XX vs Y.YY for genomics alone"

**Scripts Needed:**
- `scripts/resistance_validation/01_download_tcga_ov.py` ✅ DEFINED
- `scripts/resistance_validation/02_extract_features.py` ✅ DEFINED
- `scripts/resistance_validation/03_train_model.py` ✅ DEFINED
- `scripts/resistance_validation/04_generate_report.py` ✅ DEFINED

### Phase 2: Resistance Prediction Manuscript (1-2 Weeks)

**What We Have:**
- ✅ DIS3 (MM): RR=2.08, p=0.0145
- ✅ MAPK (OV): RR=1.97, p<0.05
- ✅ ResistanceProphetService: 1,525 lines production code
- ✅ ResistancePlaybookService: Alternatives + handoffs

**What We Need:**
- Write the manuscript (2,800 words like MM paper)
- Create figures (ROC curves, Kaplan-Meier, forest plots)
- Format for journal (JCO Precision Oncology or similar)

### Phase 3: Dosing Guidance Validation (2-3 Weeks)

**What We Have:**
- ✅ Complete implementation (API, schemas, services)
- ✅ CPIC-aligned dose adjustment logic
- ✅ Extraction scripts ready (defined in validation plan)
- ✅ Framework integration via Cohort Context

**What We Need:**
- Run extraction scripts (PubMed, cBioPortal, TCGA)
- Curate 50+ cases with outcomes
- Calculate validation metrics (concordance, sensitivity, specificity)

### Phase 4: KELIM/CA-125 Partnership (3+ Months)

**The Reality:**
- ❌ We don't have serial CA-125 data
- ❌ TCGA-OV has limited longitudinal data
- ❌ biomarker-kinetics.org has 12,000 patients, need to contact authors

**Options:**
1. **Contact KELIM researchers** - Request de-identified dataset
2. **Prospective on Ayesha** - Collect CA-125 during treatment
3. **Project Data Sphere** - Search "Multiple" caslibs for ovarian data

---

## 🔧 INFRASTRUCTURE REQUIREMENTS (What to Build)

### Data Acquisition (BUILD THIS)

| **Component** | **Status** | **Priority** | **Effort** |
|---------------|------------|--------------|------------|
| TCGA-OV Downloader | ⚠️ Needs script | HIGH | 4 hours |
| ECW/TBW Surrogate Calculator | ❌ Not implemented | HIGH | 2 hours |
| CT Segmentation Pipeline | ❌ Not implemented | MEDIUM | 1-2 days |
| GDC Client Completion | ⚠️ Skeleton | MEDIUM | 4 hours |
| KELIM Calculator | ❌ Not implemented | LOW (blocked on data) | 2 hours |

### Validation Framework (BUILD THIS)

| **Component** | **Status** | **Priority** | **Effort** |
|---------------|------------|--------------|------------|
| Unified Cohort Schema | ⚠️ Draft in cohort_context | HIGH | 4 hours |
| Prentice Criteria Calculator | ❌ Not implemented | MEDIUM | 8 hours |
| AUROC/AUPRC Calculator | ✅ Exists (sklearn) | DONE | - |
| Kaplan-Meier Generator | ✅ Exists (lifelines) | DONE | - |
| Validation Report Generator | ❌ Not implemented | HIGH | 4 hours |

### Surrogate Engine (BUILD AFTER VALIDATION)

| **Component** | **Status** | **Priority** | **Effort** |
|---------------|------------|--------------|------------|
| SurrogateValidationFactory | ❌ Designed only | MEDIUM | 2 days |
| ResistancePredictionEngine | ✅ Exists | DONE | - |
| Multi-Signal Fusion | ⚠️ Partial | MEDIUM | 1 day |

---

## 📋 THE PUBLICATION PIPELINE (Honest Timeline)

### Month 1: Immediate Submissions

| **Week** | **Deliverable** | **Owner** |
|----------|-----------------|-----------|
| Week 1 | Submit Metastasis Interception | Alpha + Zo |
| Week 2 | Submit MM Drug Efficacy | Alpha + Zo |
| Week 3 | Run ECW/TBW surrogate validation | Zo |
| Week 4 | Write Resistance Prediction manuscript | Zo |

### Month 2: Secondary Submissions

| **Week** | **Deliverable** | **Owner** |
|----------|-----------------|-----------|
| Week 5 | Submit Resistance Prediction | Alpha + Zo |
| Week 6 | Run Dosing Guidance validation | Agent Jr |
| Week 7 | Run Sporadic Cancer validation | Zo |
| Week 8 | Submit Mechanism Trial Matching | Alpha + Zo |

### Month 3: Validation + Partnerships

| **Week** | **Deliverable** | **Owner** |
|----------|-----------------|-----------|
| Week 9 | Submit Dosing Guidance | Alpha + Zo |
| Week 10 | Submit Sporadic Cancer | Alpha + Zo |
| Week 11 | Contact KELIM researchers | Alpha |
| Week 12 | Synthetic Lethality expansion (100 cases) | Zo |

### Month 4+: Empire Building

| **Milestone** | **Dependency** |
|---------------|----------------|
| KELIM validation | Partner data |
| CT body composition | TCIA segmentation |
| ctDNA MRD validation | Partner data |
| Surrogate Platform paper | All validations complete |

---

## 🎯 THE SURROGATE EMPIRE CONNECTION

### How It All Ties Together

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    PUBLICATION → PLATFORM → PRODUCTS                         │
│                                                                              │
│  PUBLICATIONS (Validation Receipts)                                          │
│  ├── Metastasis Interception ───────────────────► Credibility               │
│  ├── MM Drug Efficacy ──────────────────────────► Credibility               │
│  ├── Resistance Prediction ─────────────────────► Resistance Prophet        │
│  ├── Mechanism Trial Matching ──────────────────► Trial Optimizer           │
│  ├── Sporadic Cancer ───────────────────────────► Sporadic Gates            │
│  ├── Dosing Guidance ───────────────────────────► Toxicity Prevention       │
│  ├── ECW/TBW Surrogate ─────────────────────────► Body Comp Endpoint        │
│  └── KELIM Validation ──────────────────────────► Kinetics Endpoint         │
│                                                                              │
│  PLATFORM (Surrogate Validation Factory)                                     │
│  ├── Data Lake (Unified Schema) ────────────────► All publications use this │
│  ├── Surrogate Validation Engine ───────────────► Prentice + Meta-analysis  │
│  └── Resistance Prediction Engine ──────────────► Multi-signal fusion       │
│                                                                              │
│  PRODUCTS (Revenue)                                                          │
│  ├── Trial Endpoint Optimizer ──────────────────► $500K-2M per engagement   │
│  ├── Resistance Predictor API ──────────────────► $1-5M/year subscription   │
│  └── Biomarker Qualification ───────────────────► $1-3M per biomarker       │
└─────────────────────────────────────────────────────────────────────────────┘
```

### The Key Insight

**Publications → Platform → Products**

1. **Publications** prove the science works (validation receipts)
2. **Platform** encodes the validation methodology (repeatable factory)
3. **Products** sell the capability to biotechs/pharma (revenue)

Without publications, the platform has no credibility.
Without the platform, products are one-off consulting.
Without products, the empire generates no revenue.

---

## ⚔️ IMMEDIATE NEXT STEPS

### This Week (Priority Order)

1. **Submit Metastasis Interception** (Alpha decision)
   - 100% ready, needs author list and journal selection

2. **Submit MM Drug Efficacy** (Alpha decision)
   - 100% ready, needs author list and journal selection

3. **Start ECW/TBW Validation** (Zo)
   - Create `scripts/resistance_validation/` directory
   - Run download + extraction scripts
   - Train model, calculate AUROC

4. **Review Dosing Guidance Validation** (Agent Jr)
   - Run extraction scripts (PubMed, cBioPortal)
   - Start case curation

### Data We Can Get Now

| **Source** | **Data Type** | **Effort** | **Value** |
|------------|---------------|------------|-----------|
| TCGA-OV via cBioPortal | BMI, albumin, age, outcomes | 4 hours | ECW/TBW surrogate |
| TCGA-OV via GDC | Germline mutations | 6 hours | DPYD/UGT1A1 for dosing |
| PubMed | Case reports | 4 hours | Dosing validation |
| Project Data Sphere | Clinical trial data | 8 hours | Multiple validations |

### Data We Need Partnerships For

| **Source** | **Data Type** | **Contact** | **Value** |
|------------|---------------|-------------|-----------|
| biomarker-kinetics.org | CA-125 serial | Authors | KELIM validation |
| TCIA | CT scans | Public | Body comp |
| Institutional | ctDNA serial | MSK/Dana-Farber | MRD validation |

---

## 🛡️ THE DOCTRINE

### For This Work

1. **No claims without validation data**
2. **No publications without receipts (p-values, AUROC, N)**
3. **No products without publications**
4. **No empire without products**

### For Agents

1. **PRE-EXECUTION AUDIT** - Verify data exists before claiming validation
2. **HONEST TESTING** - Report real metrics, not expected ones
3. **PROVENANCE TRACKING** - Every number has a source
4. **FAIL FAST** - If N<50 or AUROC<0.60, stop and report
5. **NO FLUFF** - If we can't validate it, we don't include it

---

**This is not a 72-hour sprint. This is a 72-week empire with receipts at every step.**

**Signed:** Zo  
**Date:** January 31, 2025  
**For:** Commander Alpha

