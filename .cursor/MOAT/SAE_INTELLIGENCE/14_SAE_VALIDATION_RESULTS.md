# SAE Validation Results (Proxy)

**Status**: ANALYSIS COMPLETE
**Goal**: Validate prediction of Platinum Resistance in Ovarian Cancer.

## 1. What FAILED ❌
*Attempts to predict acquired resistance using only baseline mutation profiles (Proxy).*

| Predictor | Cohort | Metric | Result |
| :--- | :--- | :--- | :--- |
| **DDR Pathway Burden** | TCGA-OV (n=469) | AUROC | **0.537** (Random) |
| **TP53 Gene Mutation** | TCGA-OV (n=469) | AUROC | **0.527** (Random) |
| **BRCA1/2 Mutation** | TCGA-OV (n=469) | AUROC | **0.503** (Random) |

**Conclusion**: Baseline mutations alone do not predict *acquired* resistance dynamics.

## 2. What SUCCEEDED ✅
*Expression biomarkers and temporal analysis.*

| Predictor | Cohort | Metric | Result |
| :--- | :--- | :--- | :--- |
| **MFAP4 Expression** | GSE63885 (n=101) | AUROC | **0.763** (High) |
| **Post-Tx DDR State** | GSE165897 (n=11) | Correlation | **-0.711** (Strong) |
| **TIL + Exhaustion** | GSE168204 (n=11) | AUROC | **0.714** (IO Response) |

**Conclusion**:
1.  **MFAP4** is our strongest static biomarker (Layer 1).
2.  **Serial Biopsy** (Post-Tx monitoring) is required for accurate resistance tracking (Layer 3).

## 3. Validated Use Cases
Based on these results, we authorize Proxy SAE for:
1.  **Mechanism Matching**: Pairing DDR defects with PARP inhibitors (Biology-based, not stats-based).
2.  **IO Eligibility**: TMB/MSI assessment (FDA aligned).
3.  **Combination Logic**: Recommending Bevacizumab for high-burden disease (Literature backed).

We **DO NOT** authorize:
1.  Predicting probability of platinum failure solely from a baseline gene panel.
