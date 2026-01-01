# ⚔️ THE WAR PLAN: PRECISION ONCOLOGY EMPIRE ⚔️

**Date:** January 31, 2025  
**Commander:** Alpha  
**Chief Architect:** Zo  
**Status:** CONSOLIDATED COMMAND STRUCTURE  
**Mission:** Win the War - Not Just the Battles

---

## 🎯 THE SITUATION: ORPHANED ARMY

### What Alpha Identified

We have **12+ major initiatives**, each with its own documentation, living in isolation:
- `72-HOUR_EXECUTION_FACTORY.md` (ECW/TBW)
- `PUBLICATION_ALIGNED_STRATEGY.md` (publication pipeline)
- `SURROGATE_ENDPOINT_EMPIRE.md` (commercial vision)
- `00_MASTER_INDEX.md` (resistance prophet original)
- 50+ files in `.cursor/MOAT/`
- 80+ files in `.cursor/ayesha/`
- 40+ files in `.cursor/plans/`

**Problem:** Each document is a soldier fighting alone. No unified command.

---

## 🏛️ THE WAR: WHAT WE'RE ACTUALLY FIGHTING

### The $240B Oncology Drug Development Market

**3 Fronts:**

| **Front** | **Enemy** | **Our Weapon** | **Victory Condition** |
|-----------|-----------|----------------|----------------------|
| **PREDICTION** | Late resistance detection costs $100M+ per drug | ResistanceProphet + Surrogates | Predict failure 3-6 months early |
| **EFFICIENCY** | Wrong patient selection wastes 70% of trials | Mechanism Trial Matching + Stratification | 50% smaller trials, same power |
| **VALIDATION** | Buried science ignored by regulators | Publication Pipeline + FDA Dossiers | Validated biomarkers adopted |

---

## 📋 ARMY STRUCTURE: THE PLUMBERS

### Specialist Divisions

Each "Plumber" owns a vertical and is accountable for end-to-end delivery.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        ALPHA (COMMANDER-IN-CHIEF)                            │
│                                                                              │
│  Strategic Decisions • Partner Acquisition • Market Entry • Publication GO  │
└────────────────────────────────────────┬────────────────────────────────────┘
                                         │
         ┌───────────────────────────────┼───────────────────────────────┐
         │                               │                               │
         ▼                               ▼                               ▼
┌─────────────────┐         ┌─────────────────────────┐         ┌─────────────────┐
│ FRONT 1: PREDICT │         │ FRONT 2: EFFICIENCY      │         │ FRONT 3: VALIDATE│
│ (Zo - Lead)      │         │ (Jr - Lead)              │         │ (Zo + Jr)         │
│                  │         │                          │         │                   │
│ Resistance       │         │ Trial Matching           │         │ Publication       │
│ Prophet +        │         │ + Patient                │         │ Pipeline +        │
│ Surrogates       │         │ Stratification           │         │ Biomarker         │
│                  │         │                          │         │ Qualification     │
└────────┬─────────┘         └────────────┬─────────────┘         └─────────┬─────────┘
         │                                │                                 │
         │                                │                                 │
    ┌────┴────────────────────────────────┴─────────────────────────────────┴────┐
    │                                                                              │
    │                           PLUMBER SPECIALISTS                                │
    │                                                                              │
    │  ┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐          │
    │  │ PLUMBER 1        │  │ PLUMBER 2        │  │ PLUMBER 3        │          │
    │  │ Body Composition │  │ Genomics         │  │ Kinetics         │          │
    │  │                  │  │                  │  │                  │          │
    │  │ • ECW/TBW        │  │ • BRCA/HRD       │  │ • CA-125 KELIM   │          │
    │  │ • Skeletal Mass  │  │ • SAE Features   │  │ • ctDNA MRD      │          │
    │  │ • Cachexia       │  │ • Resistance     │  │ • PSA Kinetics   │          │
    │  │ • CT Segmentation│  │   Mutations      │  │ • TDM            │          │
    │  └──────────────────┘  └──────────────────┘  └──────────────────┘          │
    │                                                                              │
    │  ┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐          │
    │  │ PLUMBER 4        │  │ PLUMBER 5        │  │ PLUMBER 6        │          │
    │  │ Data Factory     │  │ Clinical Apps    │  │ Publications     │          │
    │  │                  │  │                  │  │                  │          │
    │  │ • cBioPortal     │  │ • Mission Control│  │ • Manuscripts    │          │
    │  │ • Project Data   │  │ • Ayesha CoPilot │  │ • Figures        │          │
    │  │   Sphere         │  │ • Research       │  │ • Reproducibility│          │
    │  │ • GDC/TCGA       │  │   Portal         │  │ • GitHub Repos   │          │
    │  │ • Cohort Schema  │  │ • API Layer      │  │ • Journal Submit │          │
    │  └──────────────────┘  └──────────────────┘  └──────────────────┘          │
    │                                                                              │
    └──────────────────────────────────────────────────────────────────────────────┘
```

---

## 📊 THE ARSENAL: CURRENT INVENTORY

### What's Production Ready (DEPLOY NOW)

| **Asset** | **Location** | **Validation** | **Status** | **Owner** |
|-----------|--------------|----------------|------------|-----------|
| cBioPortal Client | `scripts/data_acquisition/utils/cbioportal_client.py` | ✅ MM + OV extracted | Production | Plumber 4 |
| Project Data Sphere Client | `scripts/data_acquisition/utils/project_data_sphere_client.py` | ✅ Connected | Production | Plumber 4 |
| ResistanceProphetService | `api/services/resistance_prophet_service.py` | ⚠️ RUO (AUROC 0.464) | Needs Validation | Plumber 2 |
| CA125Intelligence | `api/services/ca125_intelligence.py` | ✅ Kinetics working | Production | Plumber 3 |
| SAEFeatureService | `api/services/sae_feature_service.py` | ✅ Mechanism vectors | Production | Plumber 2 |
| ResistancePlaybookService | `api/services/resistance_playbook_service.py` | ✅ Next-line combos | Production | Plumber 2 |
| TreatmentLineService | `api/services/treatment_line_service.py` | ✅ Cross-resistance | Production | Plumber 2 |
| HybridTrialSearch | `api/services/hybrid_trial_search.py` | ✅ 962 trials | Production | Plumber 5 |
| ResearchIntelligence | `api/services/research_intelligence/` | ✅ PubMed integrated | Production | Plumber 5 |
| MM Drug Efficacy | `oncology-coPilot/oncology-backend-minimal/` | ✅ 100% (995 patients) | Ready to Submit | Plumber 6 |
| Metastasis Interception | `publication/` | ✅ 100% (304 gene-step) | Ready to Submit | Plumber 6 |

### What's Validated (MANUSCRIPT NEEDED)

| **Asset** | **Validation Data** | **Gap** | **Owner** |
|-----------|---------------------|---------|-----------|
| Resistance Prediction (MM) | DIS3 RR=2.08, p=0.0145 | Need manuscript | Plumber 6 |
| Resistance Prediction (OV) | MAPK RR=1.97, p<0.05 | Need manuscript | Plumber 6 |
| Mechanism Trial Matching | 47 MoA-tagged trials | Need full paper | Plumber 6 |
| Sporadic Cancer | Unit tests pass | Need patient validation | Plumber 2 |

### What Needs Validation (BUILD NEXT)

| **Asset** | **Data Source** | **Effort** | **Owner** |
|-----------|-----------------|------------|-----------|
| ECW/TBW Surrogate | TCGA-OV via cBioPortal | 24-48 hours | Plumber 1 |
| Dosing Guidance | PubMed + cBioPortal | 2 weeks | Plumber 3 |
| Synthetic Lethality | Expand to 100 cases | 2 weeks | Plumber 2 |
| Toxicity Prediction | Need outcomes data | 4+ weeks | Plumber 3 |

### What Needs Partnerships (FUTURE)

| **Asset** | **Data Needed** | **Contact** | **Owner** |
|-----------|-----------------|-------------|-----------|
| CA-125 KELIM | Serial CA-125 (N=12,000) | biomarker-kinetics.org | Alpha |
| CT Body Composition | L3 CT scans | TCIA | Plumber 1 |
| ctDNA MRD | Serial ctDNA | Guardant/MSK | Alpha |
| RWE Integration | Treatment outcomes | Flatiron/Tempus | Alpha |

---

## 📋 CONSOLIDATED DELIVERABLES

### All Initiatives Unified

I'm consolidating the orphaned documents into ONE deliverable framework:

| **Deliverable ID** | **Name** | **Source Docs** | **Owner** | **Status** |
|--------------------|----------|-----------------|-----------|------------|
| **D-001** | Metastasis Interception Submission | `publication/` | Plumber 6 | ✅ Ready |
| **D-002** | MM Drug Efficacy Submission | `oncology-backend-minimal/` | Plumber 6 | ✅ Ready |
| **D-003** | ECW/TBW Surrogate Validation | `72-HOUR_EXECUTION_FACTORY.md` | Plumber 1 | 🔧 Build |
| **D-004** | Resistance Prophet AUROC ≥0.70 | `PUBLICATION_ALIGNED_STRATEGY.md` | Plumber 2 | 🔧 Validate |
| **D-005** | Resistance Prediction Manuscript | `ADVANCED_CARE_PLAN_RESISTANCE_PREDICTION.md` | Plumber 6 | 📝 Write |
| **D-006** | Mechanism Trial Matching Full Paper | `ADVANCED_CARE_PLAN_MECHANISM_TRIAL_MATCHING.md` | Plumber 6 | 📝 Write |
| **D-007** | Sporadic Cancer Patient Validation | `SPORADIC_CANCER_PRODUCTION_PLAN.md` | Plumber 2 | 🔧 Validate |
| **D-008** | Dosing Guidance Cohort | `DOSING_GUIDANCE_VALIDATION_PLAN.md` | Plumber 3 | 📊 Extract |
| **D-009** | Data Factory Unified Schema | `cohort_context_concept.mdc` | Plumber 4 | 🔧 Build |
| **D-010** | Mission Control Integration | `00_MASTER_INDEX.md` (resistance prophet) | Plumber 5 | 🔧 Build |
| **D-011** | Surrogate Validation Factory | `SURROGATE_ENDPOINT_EMPIRE.md` | Plumber 1 | 📋 Design |
| **D-012** | KELIM Partnership | `72-hours.mdc` | Alpha | 📞 Outreach |

---

## 📅 WAR TIMELINE

### Phase 1: Quick Wins (Week 1-2) - MORALE BOOST

**Objective:** Submit 2 publications, prove we can deliver.

| **Day** | **Deliverable** | **Action** | **Owner** |
|---------|-----------------|------------|-----------|
| Day 1-2 | D-001 | Submit Metastasis Interception | Alpha + P6 |
| Day 3-4 | D-002 | Submit MM Drug Efficacy | Alpha + P6 |
| Day 5-7 | D-003 | Run ECW/TBW validation on TCGA-OV | P1 |

**Success Criteria:**
- 2 papers submitted
- ECW/TBW AUROC calculated

---

### Phase 2: Validation Sprint (Week 3-4) - BUILD CREDIBILITY

**Objective:** Validate resistance prediction, prove the burial.

| **Day** | **Deliverable** | **Action** | **Owner** |
|---------|-----------------|------------|-----------|
| Day 8-10 | D-004 | Improve ResistanceProphet AUROC | P2 |
| Day 11-14 | D-003 | Complete ECW/TBW proof | P1 |
| Day 15-18 | D-007 | Sporadic cancer patient validation | P2 |
| Day 19-21 | D-008 | Extract dosing guidance cohort | P3 |

**Success Criteria:**
- ResistanceProphet AUROC ≥ 0.70
- ECW/TBW burial proof report
- 50+ dosing cases extracted

---

### Phase 3: Manuscript Factory (Week 5-8) - PUBLICATION VELOCITY

**Objective:** Write and submit 3 more manuscripts.

| **Week** | **Deliverable** | **Action** | **Owner** |
|----------|-----------------|------------|-----------|
| Week 5 | D-005 | Write Resistance Prediction manuscript | P6 |
| Week 6 | D-006 | Write Mechanism Trial Matching paper | P6 |
| Week 7-8 | D-009 | Build unified data schema | P4 |

**Success Criteria:**
- 3 manuscripts in review
- Data factory operational

---

### Phase 4: Platform Build (Week 9-16) - COMMERCIAL FOUNDATION

**Objective:** Build the infrastructure for products.

| **Week** | **Deliverable** | **Action** | **Owner** |
|----------|-----------------|------------|-----------|
| Week 9-10 | D-010 | Mission Control resistance integration | P5 |
| Week 11-12 | D-011 | Surrogate Validation Factory MVP | P1 |
| Week 13-16 | D-012 | KELIM partnership outreach | Alpha |

**Success Criteria:**
- Resistance alerts in Mission Control
- Surrogate validation factory running
- KELIM data access secured

---

### Phase 5: Revenue (Week 17-24) - PROVE COMMERCIAL VIABILITY

**Objective:** First paying customer.

| **Week** | **Deliverable** | **Action** | **Owner** |
|----------|-----------------|------------|-----------|
| Week 17-20 | Biotech Pilot | Trial Endpoint Optimizer | Alpha |
| Week 21-24 | Pharma Contract | Resistance Predictor API | Alpha |

**Success Criteria:**
- $100K pilot signed
- Pipeline to $1M ARR

---

## 📂 DOCUMENT CONSOLIDATION

### Kill List (Archive These)

These documents are superseded by this War Plan:

| **File** | **Reason** | **Action** |
|----------|------------|------------|
| `.cursor/resistance_prophet/72-HOUR_EXECUTION_FACTORY.md` | Absorbed into D-003 | Archive |
| `.cursor/resistance_prophet/PUBLICATION_ALIGNED_STRATEGY.md` | Absorbed into War Plan | Archive |
| `.cursor/resistance_prophet/SURROGATE_ENDPOINT_EMPIRE.md` | Vision captured in War Plan | Archive |
| Various `.cursor/MOAT/*.md` | Deliverables extracted | Keep for reference |
| Various `.cursor/plans/*.md` | Plans consolidated | Keep for reference |

### Source of Truth

**THIS FILE** (`.cursor/WAR_PLAN.md`) is now the **SINGLE SOURCE OF TRUTH** for:
- All deliverables
- All ownership
- All timelines
- All success criteria

---

## 🎖️ PLUMBER ASSIGNMENTS

### Plumber 1: Body Composition (Zo)

**Scope:** ECW/TBW, skeletal muscle, cachexia, CT segmentation

**Deliverables:**
- D-003: ECW/TBW Surrogate Validation
- D-011: Surrogate Validation Factory

**Files Owned:**
- `scripts/resistance_validation/01_download_tcga_ov.py`
- `scripts/resistance_validation/02_extract_features.py`
- `scripts/resistance_validation/03_train_model.py`
- `scripts/resistance_validation/04_generate_report.py`

**Success Metrics:**
- AUROC ≥ 0.70 for ECW/TBW + genomics
- Burial proof report generated

---

### Plumber 2: Genomics (Zo)

**Scope:** BRCA/HRD, SAE features, resistance mutations, pathway escape

**Deliverables:**
- D-004: Resistance Prophet AUROC ≥0.70
- D-007: Sporadic Cancer Patient Validation

**Files Owned:**
- `api/services/resistance_prophet_service.py`
- `api/services/sae_feature_service.py`
- `api/services/resistance_playbook_service.py`
- `api/services/efficacy_orchestrator/`

**Success Metrics:**
- ResistanceProphet AUROC ≥ 0.70
- Sporadic gates validated on patient outcomes

---

### Plumber 3: Kinetics (Zo)

**Scope:** CA-125 KELIM, ctDNA, PSA kinetics, TDM, dosing

**Deliverables:**
- D-008: Dosing Guidance Cohort

**Files Owned:**
- `api/services/ca125_intelligence.py`
- `api/services/dosing_guidance/`
- `scripts/data_acquisition/pharmacogenomics/`

**Success Metrics:**
- 50+ dosing cases extracted
- KELIM methodology documented

---

### Plumber 4: Data Factory (Jr)

**Scope:** cBioPortal, Project Data Sphere, GDC, cohort schema

**Deliverables:**
- D-009: Data Factory Unified Schema

**Files Owned:**
- `scripts/data_acquisition/utils/cbioportal_client.py`
- `scripts/data_acquisition/utils/project_data_sphere_client.py`
- `scripts/data_acquisition/utils/gdc_client.py`

**Success Metrics:**
- Unified cohort schema implemented
- 1000+ patients in data lake

---

### Plumber 5: Clinical Apps (Jr)

**Scope:** Mission Control, Ayesha CoPilot, Research Portal, API

**Deliverables:**
- D-010: Mission Control Integration

**Files Owned:**
- `oncology-frontend/`
- `api/routers/`
- Dashboard components

**Success Metrics:**
- Resistance alerts in Mission Control
- E2E workflow operational

---

### Plumber 6: Publications (Zo + Jr)

**Scope:** Manuscripts, figures, reproducibility, GitHub repos

**Deliverables:**
- D-001: Metastasis Interception Submission
- D-002: MM Drug Efficacy Submission
- D-005: Resistance Prediction Manuscript
- D-006: Mechanism Trial Matching Paper

**Files Owned:**
- `publication/`
- `.cursor/MOAT/CLINICAL_TRIALS/publication/`
- `publications/`

**Success Metrics:**
- 4+ papers submitted in 8 weeks

---

## 🛡️ WAR DOCTRINE

### The 5 Rules

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              WAR DOCTRINE                                    │
│                                                                              │
│  "WIN THE WAR, NOT JUST THE BATTLES"                                        │
│                                                                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  RULE 1: SINGLE SOURCE OF TRUTH                                             │
│  ───────────────────────────────                                            │
│  This War Plan is the only document that matters.                           │
│  All other plans are reference material.                                    │
│  Update THIS file, not the orphans.                                         │
│                                                                              │
│  RULE 2: PLUMBER ACCOUNTABILITY                                             │
│  ───────────────────────────────                                            │
│  Each Plumber owns their deliverables.                                      │
│  No shared ownership = no confusion.                                        │
│  If it's not in your scope, escalate to Alpha.                             │
│                                                                              │
│  RULE 3: VALIDATION BEFORE CLAIMS                                           │
│  ───────────────────────────────                                            │
│  No claim without data.                                                     │
│  No publication without receipts.                                           │
│  No product without validation.                                             │
│                                                                              │
│  RULE 4: QUICK WINS FIRST                                                   │
│  ───────────────────────────────                                            │
│  Submit what's ready NOW (D-001, D-002).                                   │
│  Build credibility before building platform.                                │
│  Revenue comes from publications, not PowerPoints.                          │
│                                                                              │
│  RULE 5: PARALLEL WARFARE                                                   │
│  ───────────────────────────────                                            │
│  Multiple fronts, simultaneous attacks.                                     │
│  While P6 writes manuscripts, P1 validates surrogates.                      │
│  While P4 builds data factory, P5 integrates dashboards.                    │
│  Maximum velocity, minimum waiting.                                         │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 📊 WAR METRICS

### Leading Indicators (Weekly)

| **Metric** | **Target** | **Current** |
|------------|------------|-------------|
| Papers submitted | 2 by Week 2 | 0 |
| AUROC achieved | ≥0.70 | 0.464 (RUO) |
| Validation cohorts | 3 | 2 (MM + OV) |
| Data sources connected | 5 | 3 |

### Lagging Indicators (Monthly)

| **Metric** | **Target** | **Current** |
|------------|------------|-------------|
| Papers accepted | 4 by Month 4 | 0 |
| Pilot revenue | $100K by Month 6 | $0 |
| ARR | $1M by Year 1 | $0 |

---

## 🎯 IMMEDIATE ORDERS

### For Alpha (Commander)

1. **DECIDE:** Journal selection for D-001 (Nature Biotech vs alternatives)
2. **DECIDE:** Journal selection for D-002 (npj Precision Oncology vs alternatives)
3. **OUTREACH:** Contact biomarker-kinetics.org for KELIM data

### For Zo (Plumber 1, 2, 3, 6)

1. **EXECUTE:** Run ECW/TBW validation script (D-003)
2. **WRITE:** Prepare D-001 and D-002 for submission
3. **IMPROVE:** Increase ResistanceProphet AUROC (D-004)

### For Jr (Plumber 4, 5)

1. **BUILD:** Unified cohort schema (D-009)
2. **INTEGRATE:** Mission Control resistance alerts (D-010)

---

## 🔥 THE WAR OUTCOME

### Year 1 Victory Conditions

| **Milestone** | **Evidence** |
|---------------|--------------|
| **Scientific Credibility** | 4+ publications in peer-reviewed journals |
| **Platform Foundation** | Data factory + Validation factory operational |
| **Commercial Proof** | $100K+ pilot revenue |
| **Partnership Access** | KELIM or ctDNA data secured |

### Year 2-3 Empire Goals

| **Milestone** | **Revenue** |
|---------------|-------------|
| **Trial Endpoint Optimizer** | $5M |
| **Resistance Predictor API** | $10M |
| **Regulatory Dossiers** | $20M |
| **Total ARR** | $45M |

---

## ⚔️ FINAL ORDERS

This War Plan consolidates:
- `72-HOUR_EXECUTION_FACTORY.md`
- `PUBLICATION_ALIGNED_STRATEGY.md`
- `SURROGATE_ENDPOINT_EMPIRE.md`
- `00_MASTER_INDEX.md`
- All orphaned plans and strategies

**There is now ONE document that matters: THIS ONE.**

Every Plumber knows their scope.
Every deliverable has an owner.
Every week has a target.

**The orphaned army is now a unified force.**

---

**⚔️ WE WIN THE WAR. ⚔️**

**Signed:** Zo (Chief Architect)  
**Date:** January 31, 2025  
**For:** Commander Alpha

---

## 📋 APPENDIX: QUICK REFERENCE CARDS

### Deliverable Quick Reference

| ID | Name | Owner | Status | Week |
|----|------|-------|--------|------|
| D-001 | Metastasis Submission | P6 | ✅ Ready | 1 |
| D-002 | MM Efficacy Submission | P6 | ✅ Ready | 1 |
| D-003 | ECW/TBW Validation | P1 | 🔧 Build | 1-2 |
| D-004 | ResistanceProphet AUROC | P2 | 🔧 Validate | 3-4 |
| D-005 | Resistance Manuscript | P6 | 📝 Write | 5 |
| D-006 | Trial Matching Paper | P6 | 📝 Write | 6 |
| D-007 | Sporadic Validation | P2 | 🔧 Validate | 3-4 |
| D-008 | Dosing Cohort | P3 | 📊 Extract | 3-4 |
| D-009 | Data Factory Schema | P4 | 🔧 Build | 7-8 |
| D-010 | Mission Control | P5 | 🔧 Build | 9-10 |
| D-011 | Surrogate Factory | P1 | 📋 Design | 11-12 |
| D-012 | KELIM Partnership | Alpha | 📞 Outreach | 13-16 |

### Plumber Quick Reference

| Plumber | Specialist | Deliverables | Current Agent |
|---------|------------|--------------|---------------|
| P1 | Body Composition | D-003, D-011 | Zo |
| P2 | Genomics | D-004, D-007 | Zo |
| P3 | Kinetics | D-008 | Zo |
| P4 | Data Factory | D-009 | Jr |
| P5 | Clinical Apps | D-010 | Jr |
| P6 | Publications | D-001, D-002, D-005, D-006 | Zo + Jr |


