# Multi-Biomarker v3: Expression Signatures — Honest Results

**Date**: December 25, 2024  
**Status**: Completed — No breakthrough

---

## Executive Summary

**Expression signatures did NOT improve platinum response prediction.**

| Model | AUROC | Status |
|-------|-------|--------|
| DDR_bin alone | 0.64 | Baseline |
| Immune score alone | 0.65 | Similar |
| DDR_bin + Immune | **0.67** | **Modest improvement** |
| All features | 0.58 | Worse (overfitting) |

---

## What We Tried

### Expression Data Acquired
- **300 samples** from TCGA-OV PanCancer Atlas
- **61 genes** across 6 signatures (HRD, Immune, Proliferation, MAPK, PI3K, Apoptosis)
- **76 patients** matched to our Tier-3 cohort (DDR_bin + expression)

### Individual Signature Performance (Platinum Response)

| Signature | AUROC | Interpretation |
|-----------|-------|----------------|
| **Immune score** | **0.65** | Best individual |
| **DDR_bin (SAE)** | **0.64** | Second best |
| TP53 status | 0.57 | Weak |
| HRD expression | 0.51 | No signal |
| Proliferation | 0.43 | No signal |
| MAPK expression | N/A | No signal |
| PI3K expression | N/A | No signal |

### Combined Models

| Model | CV-AUROC | Notes |
|-------|----------|-------|
| DDR_bin + Immune | **0.67 ± 0.19** | Best combination |
| DDR_bin + TP53 + all | 0.58 | Overfitting (too many features) |

---

## Unexpected Finding: Immune Paradox

**High immune infiltration correlates with RESISTANCE, not sensitivity.**

| Group | Mean Immune Score |
|-------|-------------------|
| Resistant | +0.23 |
| Sensitive | -0.08 |

**Possible explanations:**
1. **Immune evasion**: High immune = tumor already under immune pressure → developed escape
2. **Inflammation**: Chronic inflammation promotes resistance
3. **Confounding**: Advanced stage tumors have both higher immune infiltration AND worse outcomes

---

## Why Expression Didn't Work

1. **Z-score normalization**: Expression data is already normalized, washing out absolute differences
2. **Pathway-level vs gene-level**: Simple gene mean doesn't capture pathway dysfunction
3. **Sample mismatch**: Only 76/149 Tier-3 patients had expression data
4. **Small resistant group**: Only 13 resistant patients → unstable estimates

---

## What We Learned

### Expression signatures for platinum response are weak
- HRD expression: AUROC 0.51 (random)
- Proliferation: AUROC 0.43 (worse than random)
- These don't add predictive value

### Immune score shows a signal, but counterintuitive
- High immune → Resistant
- Needs biological validation
- May be confounded by stage/grade

### DDR_bin remains the best single feature
- AUROC 0.64 for platinum response (modest)
- AUROC 0.72 for extreme survival (strong)
- Adding expression doesn't help meaningfully

---

## Comparison to Manager's Predictions

| Manager's Prediction | Reality |
|---------------------|---------|
| MAPK_bin: AUROC +0.05-0.10 | ❌ No signal (low mutation rate) |
| PI3K_bin: AUROC +0.03 | ❌ No signal (low mutation rate) |
| HRD expression: AUROC +0.08-0.12 | ❌ No improvement (AUROC 0.51) |
| Immune: AUROC +0.05-0.08 | ✅ Modest (AUROC 0.65, but wrong direction) |
| Combined: AUROC 0.72-0.76 | ❌ Not achieved (AUROC 0.67) |

**Manager's strategy assumed mutations/expression would capture bypass pathways. In ovarian cancer, this doesn't work — the biology is different.**

---

## Final Recommendation

### For Platinum Response Prediction
- **Don't pursue expression signatures** — they don't add value
- **Accept that baseline genomics can't predict acquired resistance**
- **Focus on prospective serial monitoring** (DDR-MONITOR trial)

### For Publication
- **Lead with survival prediction** (DDR_bin AUROC 0.73 for extreme survival)
- **Acknowledge platinum prediction failure** honestly
- **Position as prognostic, not predictive**

### For Future Work
1. **Serial ctDNA** (ΔDDR_bin over time)
2. **Germline integration** (BRCA status)
3. **Functional assays** (HRD testing, RAD51 foci)

---

## Bottom Line

**Expression data was a "fire in the hole" that didn't ignite.**

The multi-pathway strategy sounds logical but doesn't work in ovarian cancer because:
- MAPK/PI3K mutations are rare (<5%)
- Expression signatures are weak for platinum prediction
- Acquired resistance can't be predicted from baseline

**The path forward is prospective monitoring, not more retrospective features.**

