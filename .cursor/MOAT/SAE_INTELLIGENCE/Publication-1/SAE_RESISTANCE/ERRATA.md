# ERRATA / CANONICAL CLAIMS (SAE_RESISTANCE)

**Purpose:** This folder contains multiple analysis threads created at different times. This doc locks the current “source of truth” so we don’t mix cohorts, labels, or endpoints.

---

## 1) Two different “resistance” label contracts exist (and they disagree)

### A) **Tier-3 cohort labels (internal contract)**

- **File**: `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json`
- **Size**: 149 patients (24 positive = `resistant` + `refractory`, 125 `sensitive`)
- **What was measured**:
  - TRUE SAE: 29 selected SAE features (logistic regression, 5-fold CV)
  - PROXY SAE baseline: DDR gene count (logistic regression, 5-fold CV)
  - DDR_bin: 9 “diamond” SAE features (aggregate score)
- **Headline numbers (internal)**:
  - TRUE SAE **mean CV-AUROC = 0.783** (see `.../checkpoints/true_sae_diamonds_baseline.v1.json`)
  - DDR_bin **AUROC(resistant)=0.698** on Tier-3 (see `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_tcga_ov/report.json`)

**Important limitation:** This cohort is produced by a feature-extraction pipeline with substantial extraction errors and is vulnerable to **selection bias / coverage/variant-count confounding**. Treat these results as *internal* until they replicate on an external cohort with a clean label contract.

### B) **TCGA-OV “platinum response” labels (external-style contract)**

- **File**: `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/OV_PLATINUM_TRUE_SAE_cohort.v2.json`
- **Linked size**: 161 patients (21 resistant, 140 sensitive)
- **Headline number (replication check)**:
  - DDR_bin **AUROC(resistant)=0.517** (≈ random)  
    See `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_ov_platinum_TRUE_SAE_v2/report.json`

**Canonical interpretation:** Under TCGA-style platinum labels, DDR_bin as currently defined does **not** predict platinum response.

---

## 2) What IS externally validated for platinum resistance

### **MFAP4 / EMT (expression)**

- **External cohort**: GSE63885 (n=101; 34 resistant / 67 sensitive)
- **Result**: **MFAP4 AUROC = 0.763** (validated)  
  See `VALIDATION_SUMMARY_FINAL.md` and `WHY_WE_ARE_OFF.md`.

**Canonical claim:** MFAP4 is the strongest externally validated platinum resistance biomarker currently present in this workspace.

---

## 3) What we should (and should not) claim right now

### ✅ Safe, current claims

- **MFAP4 predicts platinum resistance externally** (GSE63885 AUROC 0.763).
- **TRUE SAE / DDR_bin show internal signal on Tier-3 labels** (AUROC 0.783 / 0.698), but this is *not yet replicated* and may be confounded.

### ❌ Claims that are not currently supported

- “TRUE SAE predicts platinum resistance in ovarian cancer” (without specifying the Tier-3 label contract + caveats).
- “DDR_bin is a treatment-selection biomarker for platinum response” (fails on TCGA-style platinum labels at AUROC ~0.517).

---

## 4) Next steps to reconcile the story

- **Decide target paper path**:
  - **Path A (fastest publishable)**: MFAP4/EMT external validation paper.
  - **Path B (methods paper)**: TRUE SAE interpretability + internal Tier-3 signal, framed as exploratory + requiring replication; include confound controls.
- **If continuing TRUE SAE**: add explicit confound controls (variant count/coverage), and replicate on an independent cohort with platinum labels.








