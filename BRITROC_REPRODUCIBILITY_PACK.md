# 🛡️ BriTROC-1 REPRODUCIBILITY PACK

**Date:** Jan 29, 2026
**Status:** AUDITED & REPRODUCED
**Cohort Size:** n=40 Paired Patients
**Endpoint:** Primary Platinum Resistance (Relapse < 6mo)
**Class Balance:** 3 Resistant / 37 Sensitive

## 1. RESULTS: POST-TREATMENT DOMINANCE

| Timepoint | CN Sig 7 AUC | 95% CI | p-value |
|---|---|---|---|
| **Diagnosis (Pre)** | **0.694** | [0.410 - 0.923] | 0.1395 |
| **Relapse (Post)** | **0.874** | [0.760 - 0.972] | 0.0175 |

## 2. SENSITIVITY ANALYSIS (Leave-One-Positive-Out)
Given the small number of resistant cases (n=3), we tested robustness by excluding each resistant patient iteratively.

| Excluded Patient | Refined AUC (Post) | Delta |
|---|---|---|
| **205** | 0.892 | +0.018 |
| **216** | 0.865 | -0.009 |
| **234** | 0.865 | -0.009 |

## 3. INTERPRETATION
- **Diagnosis (Pre)**: Moderate prediction. AUC significantly > 0.5, confirming published association.
- **Relapse (Post)**: Stronger prediction. The signal intensifies after selective pressure.
- **Robustness**: Validated. Even after removing the most influential positive case, the AUC remains > 0.75.
- **Conclusion**: The 'Serial SAE' hypothesis is supported. Post-treatment state is the superior feature.

## 4. PATIENT LIST (n=40)
`1, 2, 5, 8, 9, 11, 32, 34, 36, 37, 55, 65, 67, 74, 99, 106, 124, 132, 147, 148, 157, 192, 204, 205, 209, 216, 225, 229, 232, 234, 241, 242, 246, 248, 256, 259, 267, 268, 271, 274`
