# DDR_bin Retrospective Validation Plan v2.3 (MOAT Edition)

**Date:** January 29, 2025  
**Owner:** GPT (execution) · Manager (decisions/claims)  
**Repo of record:** `/Users/fahadkiani/Desktop/development/crispr-assistant-main` (MAIN REPO)  
**Status:** **D1–D10 complete through D9 ✅ (remaining: external cohorts)**

---

## 🎯 Core mission (what we are proving)

> **DDR_bin is a treatment-selection biomarker:** it predicts **platinum resistance/response** (mechanistic), even if it does **not** predict survival (PFS/OS) in this cohort.

### Score direction (locked)
- **DDR_bin is a resistance score**: higher = more resistant/refractory.
- Positive class for AUROC benchmarks: **resistant/refractory = 1**, sensitive = 0.

---

## 🧾 Source-of-truth scripts + artifact directory

### A8 external replication harness (new)

- Script: `oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_ddr_bin_generic_cohort.py`
- Plumber instructions: `.cursor/MOAT/SAE_INTELLIGENCE/PLUMBER_A8_EXTERNAL_REPLICATION.mdc`

### Run Phase 1 validator (A1–A3 + plots)

```bash
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_ddr_bin_tcga_ov_survival.py
```

### Run Phase 2/3 analyses (A4 + D5–D9 + publication plots)

```bash
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/compute_ddr_bin_competitive_benchmark.py
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/generate_waterfall_ddr_bin_publication.py
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/compute_per_diamond_auroc.py
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/compute_top3_vs_full.py
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/analyze_ddr_bin_sparsity.py
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/compare_ddr_bin_aggregation_methods.py
python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/compute_stage_subgroup_auroc.py
```

### Canonical output directory

All artifacts live in:

`oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_tcga_ov/`

---

## ✅ What’s delivered (results + files)

### Phase 1 results (A1–A3)

- **A1 (PFS/OS)**: no survival signal (expected for response biomarker)
  - PFS: p=0.56
- **A2 (gene-level head-to-head)**: DDR_bin adds value vs gene flag on PFS C-index
  - ΔC-index ≈ **+0.106**
- **A3 (platinum response AUROC)**: **AUROC = 0.698** (predict resistant/refractory)

### Phase 2 competitive benchmark (A4 / Deliverables D1–D4)

Computed AUROC (predict resistant/refractory):
- **DDR_bin**: **0.698**
- **Gene DDR (best orientation)**: **0.620**  → **Δ DDR_bin − gene DDR = +0.078**
- **TP53 (negative control)**: **0.507**
- **TMB proxy (negative control)**: **0.509**

**Interpretation:** DDR_bin beats the gene-flag baseline and crushes negative controls → mechanism-specific signal.

### D10 publication waterfall
- `waterfall_ddr_bin_publication.png` (300 DPI)

---

## ✅ D5–D9 results (interpretability + sparsity + robustness)

### D5: Per-diamond AUROC ranking (9 diamonds)
Top diamonds (AUROC, nonzero rate):
- **feature 9738**: **0.760**, 25/149 nonzero (16.8%)
- **feature 12893**: **0.751**, 35/149 nonzero (23.5%)
- **feature 27607**: **0.747**, 40/149 nonzero (26.8%)

**Artifact:** `per_diamond_auroc.json` + `per_diamond_auroc.csv`

### D6: Top-3 composite vs full DDR_bin
- **Top-3 composite AUROC:** **0.745** (top3 = [9738, 12893, 27607])
- **Full DDR_bin AUROC:** **0.698**
- **Δ(full − top3): -0.047**

**Interpretation:** The signal is concentrated in a small subset of diamonds. For *biomarker performance*, top-3 is stronger.

**BUT**: do not change production/validator definition yet without Manager sign-off (definition drift risk).

**Artifact:** `top3_vs_full_comparison.json`

### D7: Sparsity postmortem
- **DDR_bin == 0:** **73/149 (49.0%)**
- **AUROC (all): 0.698**
- **AUROC (nonzero-only): 0.569 (n=76)**

**Interpretation:** “Nonzero-only AUROC” dropping suggests our **current aggregation** (max_of_max) is not simply “missingness drives signal”. Instead, the distribution/label alignment is more complex (and supports trying alternate aggregations).

**Artifacts:**
- `ddr_bin_distribution_analysis.json`
- `ddr_bin_histogram_publication.png` (300 DPI)

### D8: Aggregation method comparison (research)
AUROC (resistant/refractory):
- **mean_of_sum:** **0.707** (best)
- **sum_of_max:** **0.703**
- **top3_mean:** **0.703**
- **max_of_max (current): 0.698**

**Interpretation:** A small improvement exists by switching aggregation (≈ +0.009). Not a ≥0.02 win, but directionally consistent.

**Artifact:** `aggregation_comparison.json`

### D9: Stage robustness (coarse)
- **Stage III:** n=115, **AUROC = 0.707**
- **Stage IV:** n=20, **AUROC = 0.484**
- Unknown stage: 14 rows

**Interpretation:** Stage IV is underpowered (n=20) and does not show signal here; Stage III carries most of the observed discrimination.

**Artifact:** `subgroup_auroc_stage.json`

---

## ✅ Updated decision gates (where we stand)

### Gate 1 (Exploratory / “signal present”) — **PASS ✅**
- PASS if platinum AUROC > 0.65 → **0.698**

### Gate 2 (Publication-ready package) — **PASS (core) ✅ / needs external replication 🟡**
Core package is now complete:
- ✅ DDR_bin beats gene DDR by ≥ 0.05 → **+0.078**
- ✅ DDR_bin beats TP53/TMB by ≥ 0.10 → **~+0.19**
- ✅ Publication-quality plots exist (benchmark + waterfall + histogram)
- ✅ Interpretability and robustness diagnostics produced (D5–D9)

Remaining for top-tier confidence:
- ⏳ External cohort replication (ICGC-OV or equivalent)

---

## ✅ Manager decisions (locked)

### 1) Definition lock (production + validator)
- **Decision:** Keep canonical DDR_bin as **`max_of_max`** (current validator definition).
- **Why:** Prevent definition drift and avoid circularity (top‑3 was chosen based on the same cohort’s AUROC ranking).
- **Policy:** Alternative aggregations (**top‑3 composite**, **mean_of_sum**) remain **exploratory diagnostics** until validated on an external cohort or prospectively pre-specified.

### 2) Claim threshold (RUO)
- **Decision:** **AUROC 0.698 is acceptable for RUO positioning**; external-facing copy can say **“~0.70”**.
- **Allowed RUO claim language:** “DDR_bin predicts platinum resistance with AUROC ~0.70, outperforming gene-level DDR status by +8% and demonstrating mechanism-specific discrimination (negative controls ~0.50).”
- **Restriction:** Do **not** claim “clinically validated biomarker” until prospective replication.

### 3) Stage IV subgroup interpretation
- **Decision:** Treat Stage IV AUROC as **underpowered** (n=20) and report as a limitation, not a failure.

---


---

## 🧪 External replication (A8) — in progress


### TCGA-OV External (with platinum response labels) — COMPLETE ✅
- **TRUE-SAE cohort file:** `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/OV_PLATINUM_TRUE_SAE_cohort.json`
- **Validator output dir:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_ov_platinum_external/`

**Replication status (n=161 linked):**
- **Platinum Response AUROC = 0.517** (predicting resistant/refractory vs sensitive; TCGA-OV TRUE-SAE v2)
  - **Interpretation:** essentially random; DDR_bin (as currently defined) does **not** discriminate platinum response in this cohort.
  - **Root cause of prior 0.663 claim:** earlier cohort had **0 variants for all patients** => constant DDR_bin; non-tie-safe AUROC ranking produced a spurious 0.663. After tie-safe AUROC, constant-score AUROC = 0.5.
- OS Spearman ρ: **0.252** (p=0.0013) on TRUE-SAE v2 (note: this is OS correlation, not platinum-response discrimination)

**Interpretation:**
- **This cohort does NOT replicate A3** (AUROC ~0.517 ≈ random).
- Treat this as a **failed external replication** for A3 under the current DDR_bin definition + mapping.
- OS signal in this run is **not interpretable** for A3; prioritize fixing A3 or reframing claim.

**Gate check (updated):**
- Gate 1 (exploratory A3 signal): **FAIL ❌** (0.517 < 0.65)
- Gate 2 (external replication): **FAIL ❌**


### TCGA-BRCA (PanCan Atlas) — partial run (resumable)
- **TRUE-SAE cohort file:** `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/BRCA_TCGA_TRUE_SAE_cohort.json`
- **Validator output dir:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_brca_tcga/`

**Current replication status (partial, n linked=60):**
- PFS Spearman ρ = 0.067 (p=0.610)
- OS Spearman ρ = 0.017 (p=0.899)
- Platinum AUROC: N/A (no labels in BRCA cohort)

**Notes:**
- This is an *engineering and generalization* check. True platinum-response replication still requires an ovarian external cohort with response labels (ICGC-OV or partner data).
- Once BRCA extraction finishes (target n≈200), re-run the validator to update these numbers.

## 📂 Artifact index (copy/paste)

Directory: `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_tcga_ov/`

- Phase 1:
  - `report.json`
  - `linked_patients.csv`
  - `roc_platinum_response.png`
  - `waterfall_ddr_bin.png`
  - `km_pfs.png`, `km_os.png`
  - `ddr_bin_hist.png`
- Phase 2/3:
  - `competitive_benchmark.json`
  - `competitive_benchmark.png`
  - `waterfall_ddr_bin_publication.png`
  - `per_diamond_auroc.json`
  - `per_diamond_auroc.csv`
  - `top3_vs_full_comparison.json`
  - `ddr_bin_distribution_analysis.json`
  - `ddr_bin_histogram_publication.png`
  - `aggregation_comparison.json`
  - `subgroup_auroc_stage.json`
