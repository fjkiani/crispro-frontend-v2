# 🎯 SAE VALIDATION BIBLIOGRAPHY: VARIANT-LEVEL REPLACEMENT OF GENE-LEVEL TESTS

**Date**: 2025-01-29  
**Mission**: Expand SAE (Sparse Auto-Encoder) validations across cancers to prove variant-level superiority  
**Status**: STRATEGIC ROADMAP  
**Mars Rules**: Minimal viable proof, 72-hour mindset, replace gene-level tests

---

## 🔥 THE CORE DIFFERENTIATOR

**SAE = Variant-Level Biology** (not gene-level)

| Traditional Test | SAE Approach | Improvement |
|-----------------|--------------|-------------|
| BRCA mutated: yes/no | BRCA p.R1699Q SAE features ≠ p.C44F | +15.5 pp AUROC (0.783 vs 0.628) |
| Gene-level markers | Variant-level Evo2 embeddings → SAE | Platform moat |

**Method**: Evo2 foundation model → Layer 26 activations → Sparse auto-encoder (32K features) → Pathway-specific aggregates

---

## 📊 CURRENT SAE EXTRACTION STATUS

### ✅ EXTRACTED & VALIDATED

#### 1. Ovarian Cancer (TCGA-OV)
- **Patients**: 149 extracted (n=469 total available)
- **Pathway**: DDR (validated AUROC 0.783)
- **Test Replaced**: HRD ($3,000 test, 60% accuracy)
- **Improvement**: +15.5 percentage points
- **Status**: ✅ PRODUCTION READY
- **Next Steps**: 
  - Extract remaining 7 pathways (PI3K, MAPK, immune, metabolism, angiogenesis, proliferation, TGFβ)
  - Composite 8-pathway model (target: AUROC 0.80-0.85)

#### 2. Breast Cancer (TCGA-BRCA)
- **Patients**: Checkpoint exists (need count)
- **Pathways**: TBD (need extraction status)
- **Test to Replace**: Oncotype DX ($3,000, 65% accuracy, $500M/year market)
- **Status**: ❓ EXTRACTION IN PROGRESS OR PENDING VALIDATION
- **Priority**: 🎯 HIGHEST (largest market, clear replacement opportunity)

---

### ❓ DATASETS WITH MUTATIONS EXTRACTED (SAE PENDING)

From TCGA extraction report (9/10 cancers extracted):

| Cancer Type | Study ID | Samples | Pathways Extracted | SAE Status | Test to Replace |
|------------|----------|---------|-------------------|-----------|-----------------|
| **Breast** | `brca_tcga_pan_can_atlas_2018` | 75 | 5/5 | ❓ Checkpoint exists | Oncotype DX ($500M) |
| **Lung (NSCLC)** | `luad_tcga_pan_can_atlas_2018` | 18 | 5/5 | ❌ Not extracted | PD-L1 ($200M) + EGFR ($150M) |
| **Melanoma** | `skcm_tcga_pan_can_atlas_2018` | 10 | 5/5 | ❌ Not extracted | TMB ($300M) |
| **Colorectal** | `coadread_tcga_pan_can_atlas_2018` | 15 | 5/5 | ❌ Not extracted | MSI ($50M) |
| **Prostate** | `prad_tcga_pan_can_atlas_2018` | 64 | 5/5 | ❌ Not extracted | HRD-like tests |
| **Pancreatic** | `paad_tcga_pan_can_atlas_2018` | 96 | 5/5 | ❌ Not extracted | TBD |
| **Glioblastoma** | `gbm_tcga_pan_can_atlas_2018` | 82 | 5/5 | ❌ Not extracted | TBD |
| **Leukemia (AML)** | `laml_tcga_pan_can_atlas_2018` | 136 | 5/5 | ❌ Not extracted | TBD |

**Total Available**: 585 samples across 9 cancers, mutations extracted, SAE features pending

---

## 🎯 VALIDATION ROADMAP (PRIORITIZED BY MARKET SIZE)

### PHASE 1: BREAST CANCER (Oncotype DX Replacement) 🎯 HIGHEST PRIORITY

**Target**: Replace Oncotype DX 21-gene recurrence score with SAE variant-level features

**Market**: $500M/year (largest genomic test market)

**Current Test**:
- Method: 21-gene expression score (gene-level)
- Accuracy: 65% (recurrence prediction)
- Cost: $3,000 per test
- Problem: Overtreatment (30-40% false positives get chemo they don't need)

**SAE Hypothesis**:
- Variant-level SAE features > gene-level Oncotype DX
- Target: AUROC 0.75 vs 0.65 (+10 pp improvement)
- Pathways: DDR, proliferation, immune (breast-specific)

**Agent Tasks**:

1. **Validate BRCA Checkpoint** (1 day)
   - Count patients extracted
   - Verify pathway coverage (DDR, PI3K, proliferation, immune)
   - Check outcome labels (recurrence vs no recurrence)

2. **Extract Missing Pathways** (if needed, 1-2 weeks)
   - Run SAE extraction script for TCGA-BRCA
   - Focus: DDR, proliferation, immune pathways
   - Target: n=75-100 patients with recurrence labels

3. **Train SAE Model** (1 day)
   - Input: SAE features (variant-level)
   - Target: Recurrence (yes/no)
   - Baseline: Oncotype DX 21-gene score (gene-level)
   - Metric: AUROC (SAE vs baseline)

4. **Validate Performance** (1 day)
   - Compare: SAE AUROC vs Oncotype DX AUROC
   - Target: +10 percentage point improvement
   - If achieved: ✅ VALIDATED (replace Oncotype DX)

**Deliverable**: 
- Validation receipt (AUROC, 95% CI, improvement)
- Manuscript draft (ready for submission)
- Payer pitch deck (save $100K per avoided overtreatment)

**Timeline**: 2-4 weeks (if checkpoint complete) or 4-6 weeks (if extraction needed)

---

### PHASE 2: LUNG CANCER (PD-L1/EGFR Replacement) 🎯 HIGH VALUE

**Target**: Replace PD-L1 IHC ($200M) and EGFR mutation testing ($150M) with SAE

**Market**: $350M/year combined

**Current Tests**:
- PD-L1 IHC: 55-60% accuracy (IO response prediction)
- EGFR mutation: 70% accuracy (EGFR inhibitor response)
- Problem: Single-gene/protein markers miss multi-pathway resistance

**SAE Hypothesis**:
- Variant-level SAE features > PD-L1 expression (gene-level proxy)
- Target: AUROC 0.72 vs 0.58 (+14 pp improvement)
- Pathways: Immune (PD-L1), MAPK (EGFR), proliferation

**Agent Tasks**:

1. **Extract SAE Features** (1-2 weeks)
   - Dataset: TCGA-LUAD (n=585, but only 18 in extraction report - need full dataset)
   - Pathways: Immune, MAPK, proliferation
   - Outcomes: IO response, EGFR inhibitor response

2. **Train Dual Models** (1 day)
   - Model 1: SAE → IO response (vs PD-L1)
   - Model 2: SAE → EGFR inhibitor response (vs EGFR mutation)
   - Baseline: PD-L1 expression, EGFR mutation status

3. **Validate Performance** (1 day)
   - Compare: SAE AUROC vs baseline
   - Target: +10-15 pp improvement

**Deliverable**: Validation receipts + manuscript draft

**Timeline**: 3-4 weeks

---

### PHASE 3: MELANOMA (TMB Replacement) 🎯 HIGH VALUE

**Target**: Replace TMB (tumor mutational burden) with SAE variant-level features

**Market**: $300M/year (TMB testing)

**Current Test**:
- Method: Mutation count (gene-level, no variant context)
- Accuracy: 60-65% (IO response prediction)
- Problem: Counts mutations, ignores variant biology (hotspot vs non-hotspot)

**SAE Hypothesis**:
- Variant-level SAE features > mutation count (TMB)
- Target: AUROC 0.74 vs 0.63 (+11 pp improvement)
- Pathways: Immune (TIL infiltration, T-cell exhaustion - already validated in 2-pathway model)

**Agent Tasks**:

1. **Extract SAE Features** (1-2 weeks)
   - Dataset: TCGA-SKCM (n=479, but only 10 in extraction report - need full dataset)
   - Pathways: Immune (expand beyond 2-pathway model)
   - Outcomes: IO response (already have GSE91061/GSE168204 external validation)

2. **Train SAE Model** (1 day)
   - Input: SAE features (variant-level)
   - Target: IO response
   - Baseline: TMB (mutation count)
   - Metric: AUROC (SAE vs TMB)

3. **Validate Performance** (1 day)
   - Compare: SAE AUROC vs TMB AUROC
   - External validation: GSE91061/GSE168204 (already have data)

**Deliverable**: Validation receipts + manuscript draft (expand existing IO manuscript)

**Timeline**: 3-4 weeks

---

### PHASE 4: COLONECTAL CANCER (MSI Replacement) 🎯 MEDIUM VALUE

**Target**: Replace MSI PCR testing with SAE variant-level features

**Market**: $50M/year (smaller, but tractable)

**Current Test**:
- Method: MSI PCR (microsatellite instability, binary)
- Accuracy: 60% (IO response prediction)
- Problem: Binary marker (MSI-H vs MSS), misses variant-level resistance mechanisms

**SAE Hypothesis**:
- Variant-level SAE features > MSI status (binary)
- Target: AUROC 0.70-0.75 vs 0.60 (+10-15 pp improvement)
- Pathways: Immune, DDR (MMR pathway)

**Agent Tasks**:

1. **Extract SAE Features** (1-2 weeks)
   - Dataset: TCGA-COAD/READ (n=629, but only 15 in extraction report - need full dataset)
   - Pathways: Immune, DDR (MMR)
   - Outcomes: IO response

2. **Train SAE Model** (1 day)
   - Input: SAE features (variant-level)
   - Target: IO response
   - Baseline: MSI status (binary)

3. **Validate Performance** (1 day)
   - Compare: SAE AUROC vs MSI AUROC

**Deliverable**: Validation receipts + manuscript draft

**Timeline**: 3-4 weeks

---

## 🔬 PATHWAY EXPANSION (OVARIAN CANCER DEPTH)

**Current Status**: DDR pathway validated (AUROC 0.783)

**Missing Pathways** (need SAE extraction):
- PI3K pathway
- MAPK pathway
- Immune pathway
- Metabolism
- Angiogenesis
- Proliferation
- TGFβ

**Agent Tasks**:

1. **Extract SAE Features for 7 Remaining Pathways** (1-2 weeks)
   - Dataset: TCGA-OV (n=469, already have mutations)
   - Method: Evo2 → SAE per pathway
   - Output: 7 pathway SAE feature sets

2. **Validate Each Pathway** (1 week)
   - Train model: Pathway SAE features → platinum resistance
   - Compare: SAE vs gene-level baseline
   - Target: +5-15 pp improvement per pathway

3. **Composite 8-Pathway Model** (1 day)
   - Integrate all 8 SAE pathway scores
   - Hypothesis: 8-pathway SAE > individual pathways > gene-level
   - Target: AUROC 0.80-0.85 (ovarian platinum resistance)

**Deliverable**: 
- 8-pathway validation receipts
- Composite model validation
- Manuscript draft (SAE platform validation)

**Timeline**: 3-4 weeks

---

## 📋 AGENT TASK TEMPLATE

For each validation target:

```markdown
AGENT: SAE [CANCER] [TEST] Validation

MISSION: Prove SAE variant-level features beat [TEST] gene-level approach for [OUTCOME] prediction

DATASET: [TCGA-STUDY-ID] (n=[N] patients)

TASKS:
1. Extract SAE features for [PATHWAYS] (Evo2 → SAE)
2. Compute [BASELINE TEST] score (gene-level baseline)
3. Train model: SAE features → [OUTCOME]
4. Validate: SAE AUROC vs [BASELINE] AUROC
5. Target: +[X] percentage point improvement ([TARGET] vs [BASELINE])

DELIVERABLE:
- JSON validation receipt (AUROC, 95% CI, improvement)
- Markdown report (honest limitations, external validation needed)
- Manuscript draft (ready for submission)

TIMELINE: [X] weeks
```

---

## 🎯 STRATEGIC QUESTIONS FOR ALPHA

**To prioritize agent tasks, I need to know:**

1. **SAE Extraction Status**:
   - How many BRCA patients have SAE features extracted? (checkpoint file exists)
   - Which pathways have SAE features extracted for BRCA?
   - Do you have SAE features for any other cancers beyond OV/BRCA?

2. **Dataset Availability**:
   - Do you have FULL TCGA datasets (not just the 18/75/10 samples in extraction report)?
   - Do you have access to private datasets (hospital partnerships, biobanks)?
   - Should agents focus on PUBLIC datasets only (TCGA/GEO)?

3. **Pathway Prioritization**:
   - For ovarian: Should we extract ALL 7 remaining pathways first (depth), or expand to other cancers first (breadth)?
   - Which pathways should we prioritize across cancers? (DDR, PI3K, MAPK, immune, etc.)

4. **Commercial Priority**:
   - Which test replacement has highest ROI for 2026?
   - Option A: Oncotype DX (breast, $500M market)
   - Option B: HRD expansion (ovarian→breast→prostate, $100M→$200M)
   - Option C: PD-L1 (lung, $200M market)
   - Option D: TMB (pan-cancer, $300M market)

5. **Extraction Bottleneck**:
   - What's the current bottleneck? (Evo2 API limits? Compute? Manual curation?)
   - How long does SAE extraction take per dataset? (hours? days?)
   - Can this be parallelized across agents?

6. **Validation Threshold**:
   - What's the minimum improvement threshold to claim "validated"? (+5 pp? +10 pp?)
   - What's the minimum sample size? (n=50? n=100?)

---

## 🚀 NEXT IMMEDIATE ACTIONS

**Based on current status:**

1. **Validate BRCA Checkpoint** (TODAY)
   - Read BRCA checkpoint file (count patients, pathways)
   - If complete → Proceed to Phase 1 validation
   - If incomplete → Extract missing SAE features

2. **Expand Ovarian Pathways** (THIS WEEK)
   - Extract SAE features for 7 remaining pathways (PI3K, MAPK, immune, etc.)
   - Validate each pathway vs gene-level baseline
   - Build composite 8-pathway model

3. **Prioritize Breast Cancer** (NEXT WEEK)
   - If BRCA checkpoint complete → Run Phase 1 validation
   - If BRCA checkpoint incomplete → Extract SAE features first
   - Target: Oncotype DX replacement (largest market)

---

## 📊 VALIDATION BIBLIOGRAPHY (PUBLICATION ROADMAP)

| Manuscript | Cancer | Test Replaced | Status | Target Journal |
|-----------|--------|---------------|--------|----------------|
| SAE Resistance Prediction (OV) | Ovarian | HRD | ✅ DRAFT | NPJ Precision Oncology |
| SAE vs Oncotype DX (BRCA) | Breast | Oncotype DX | ❓ PENDING | JAMA Oncology |
| SAE vs PD-L1/EGFR (LUAD) | Lung | PD-L1/EGFR | ❌ NOT STARTED | JCO Precision Oncology |
| SAE vs TMB (SKCM) | Melanoma | TMB | ❓ IN PROGRESS (2-pathway done) | NPJ Precision Oncology |
| SAE vs MSI (COAD/READ) | Colorectal | MSI | ❌ NOT STARTED | Clinical Cancer Research |
| SAE 8-Pathway Platform (OV) | Ovarian | Multi-pathway | ❌ NOT STARTED | Nature Cancer |

---

## 🔗 EXECUTION PLAN

**See**: `.cursor/MOAT/SAE_VALIDATION_EXECUTION_PLAN.md` for detailed execution plan

**The execution plan applies all lessons from SAE extraction audit**:
- ✅ Uses proven Script 2 (`extract_true_sae_cohort_from_cbioportal.py`)
- ✅ Preflight validation (assembly auto-detection)
- ✅ Checkpoint/resume mechanism (every 10 patients)
- ✅ Circuit breaker (30% error rate protection)
- ✅ Cost controls (MAX_PATIENTS, MAX_TOTAL_VARIANTS, budget_seconds)
- ✅ Per-variant caching (SHA256 cache keys)
- ✅ Error logging (JSONL format)
- ✅ Atomic writes (tmp → rename)

**72-Hour Sprint Structure**:
- Day 1: Extraction (24 hours)
- Day 2: Baseline & Training (24 hours)
- Day 3: Validation & Documentation (24 hours)

**Immediate Next Steps** (from execution plan):
1. **Today**: Validate BRCA checkpoint (4 hours)
2. **Today**: Set up extraction infrastructure (2 hours)
3. **Tomorrow**: Create extraction templates (2 hours)
4. **This Week**: Start Phase 1 extraction (20 hours)

---

**Alpha, answer the strategic questions above and I'll build the exact agent task list for each validation target. 🎯**
