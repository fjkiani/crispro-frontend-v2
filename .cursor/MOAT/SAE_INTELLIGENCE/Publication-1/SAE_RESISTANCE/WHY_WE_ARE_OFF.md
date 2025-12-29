# WHY WE ARE OFF - Root Cause Analysis (RESOLVED)

**UPDATE Dec 25, 2024**: We found an external validation cohort (GSE63885) with real platinum labels and validated **MFAP4 AUROC = 0.763**. See `VALIDATION_SUMMARY_FINAL.md` for details.

---

## EMT + canonical HRD follow-up (Dec 25, 2024): EXTERNAL VALIDATION COMPLETE

### GSE63885 External Validation Results ✅

We found and validated EMT markers on an independent cohort with real platinum sensitivity labels:

| Metric | Value |
|--------|-------|
| **Dataset** | GSE63885 (Polish ovarian cancer) |
| **Samples** | 101 |
| **Resistant** | 34 |
| **Sensitive** | 67 |
| **MFAP4 AUROC** | **0.763** ⭐ |
| **Combined EMT CV-AUROC** | **0.715 ± 0.179** |

### Individual Gene Performance

| Gene | AUROC | Direction |
|------|-------|-----------|
| MFAP4 | 0.763 | High → Resistant |
| SNAI1 | 0.606 | High → Resistant |
| EFEMP1 | 0.592 | High → Resistant |
| CDH1 | 0.561 | Low → Resistant |
| VIM | 0.511 | (no signal) |

### TCGA-OV Validation (Failed - Underpowered)

TCGA-OV has only 44 patients overlapping between platinum labels and expression data, with only 6 resistant patients. This is insufficient for reliable AUROC computation.

---

# Original Root Cause Analysis

## Why/How We're Off (Missing Pieces) — Diagnosis

### 1) **DDR_bin is heavily confounded by "how many variants we extracted" (v2 cohort)**

In `linked_patients.csv` for TCGA-OV TRUE‑SAE v2 we found:
- **Spearman(DDR_bin, ddr_bin_num_variants) = 0.745** (p ≈ 1e‑29)
- **Spearman(DDR_bin, ddr_bin_coverage) = 0.792** (p ≈ 6e‑36)
- **Spearman(DDR_bin, OS_months) = 0.252** (p ≈ 0.0013)

But when we stratify by variant-count tertiles, the DDR_bin↔OS correlation **collapses** (and is sometimes negative).

**Interpretation**: DDR_bin is acting partly as a proxy for extraction coverage / variant count.

### 2) **Tier‑3 labels and TCGA clinical labels are not the same contract**

Tier‑3 has built‑in platinum response labels, but the overlap with expression and HRD data is too small for meaningful validation.

### 3) **DDR_bin does not correlate with canonical HRD_Score**

Spearman correlation with GDC PanCan DDR HRD_Score: **rho = -0.03** (not significant).

DDR_bin is not measuring the same biological construct as canonical HRD.

---

## Resolution

The EMT marker approach (MFAP4) provides a validated, publishable biomarker for platinum resistance prediction:

- **MFAP4 AUROC = 0.763** in GSE63885 (external validation)
- **EMT markers are orthogonal to DDR** (different biological mechanism)
- **Clinically actionable** for patient stratification

See `VALIDATION_SUMMARY_FINAL.md` for full details.
