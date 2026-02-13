# SAE Failure Postmortem: True Feature Validation

**Date**: January 28, 2025
**Status**: Critical Failure Analysis (Tier 3)
**Objective**: Analyze why TRUE SAE (32k features) failed to unblock Feature->Pathway mapping despite successful infrastructure.

## 1. Executive Summary
While the infrastructure for extracting 32k SAE features using Evo2 on Modal was successfully built and deployed, the downstream biomarker discovery failed to yield FDR-significant features.

*   **Infrastructure**: ✅ Successful (H100/Modal deployment, 149 patients processed).
*   **Biomarker Signal**: ❌ Failed (0 FDR-significant features at p<0.05).
*   **Result**: Cannot currently rely on TRUE SAE for primary clinical decision making. **Proxy SAE** (gene-level) remains production standard.

## 2. The Artifact Reality
We have three tiers of data, often confused in documentation:
*   **Tier 1 (N=10)**: `sae_features_tcga_ov_platinum.json`. Useless. Too small.
*   **Tier 2 (N=49)**: `Tier2_statistical_minimum.json`. Intermediate.
*   **Tier 3 (N=149)**: `Tier3_validation_cohort.json`. Best available. 125 Sensitive / 24 Resistant.

**Critical Finding**: Even at Tier 3 (N=149), no features survived Benjamini-Hochberg correction regarding platinum resistance.

## 3. Root Causes
1.  **Sample Size/Imbalance**: 24 resistant patients is insufficient for 32k-dimensional discovery. Power analysis suggests ~64/group needed for medium effects.
2.  **Statistical Harshness**: Multiple testing correction (32k tests) requires extremely low p-values ($< 10^{-5}$), which noisy biological data rarely achieves without massive N.
3.  **Engineering Diversions**: Cycles spent on assembly mismatches (GRCh37/38) and Modal payload limits reduced time for cohort expansion.

## 4. The "Diamond" Excavation Plan
Despite 0 FDR-significant features, Tier 3 revealed **11 Large-Effect Features** (|d| > 0.5) that did not pass p-value thresholds but show strong directional signal.

**High-Priority Tasks ("Plumber" Protocol):**
1.  **Map Large-Effect Features**: Manually investigate the 11 "Diamond" candidates (e.g., Feature 27607, d=+0.635). Do they map to known drivers (DIS3, NF1)?
2.  **Head-to-Head**: Run predictive probe (Logistic Regression): TRUE SAE vectors vs PROXY SAE vectors. If TRUE does not add AUC, deprecate until N > 500.
3.  **Canonicalize**: Stop using Tier 1 files. Merge Tier 2/3 into a single canonical dataset.

## 5. Strategic Pivots
*   **Production**: Stick to **Proxy SAE** (Gene-level: DIS3, TP53, etc.) for clinical resistance predictions. It is validated and interpretable.
*   **Research (RUO)**: Keep TRUE SAE as "Steerability V1" experiment. Use "Diamond" features for hybrid bins (e.g., "True SAE DDR Bin") rather than full feature vectors.

## 6. Validated Proxy Wins (The Baseline)
While TRUE SAE struggled, Proxy SAE validated successfully:
*   **DIS3**: RR=2.08, p=0.0145 (MMRF)
*   **NF1**: RR=2.10 (TCGA-OV)
*   **MAPK**: RR=1.97
This confirms the *biological* hypothesis (pathways drive resistance), even if the *feature* hypothesis (SAE captures it better) is technically stalled.

## 7. Validation Incident (Jan 29, 2026)

**The Event**: A claim of "Diamond SAE AUROC 0.783" was widely propagated. Use of a nested cross-validation audit revealed the true Mean AUROC to be **0.555 (Random)**.

**Root Cause Hypothesis (Leakage)**:
The "Diamond" features were likely selected using the *entire* cohort (Feature Selection Bias) before cross-validation. This "peeking" creates a subset of features that perform well on the training data but fail on unseen folds. When a specific fold partition happens to align with the selection, the score spikes (Fold 4 = 0.784), creating a "hallucination" of success.

**Correction Protocol**:
1.  **Strict Nesting**: Feature selection MUST happen *inside* the CV loop.
2.  **Patient-Level Splits**: Ensure variants from the same patient never span train/test.
3.  **Head-to-Head**: Future tests must run TRUE vs PROXY on locked folds.

