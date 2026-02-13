# 00 SAE Master Source of Truth

**Status**: CURRENT (Audited Jan 29, 2026)
**Authority**: Zo (Nyx)

## 1. Executive Summary: The Reality
After a comprehensive audit, the status of SAE Intelligence is:

1.  **Diamond True SAE (FAILED)**: The previous claim of AUROC 0.783 was found to be a **Single-Fold Hallucination** (Best Fold=0.784, Mean=0.555). It is random noise.
2.  **MFAP4 Biomarker (SUCCESS)**: Independently validated (AUROC 0.763) for intrinsic platinum resistance.
3.  **Proxy SAE (ACTIVE)**: Useful for gene-level pathway aggregation and trial matching, but limited predictive power (AUROC 0.628).
4.  **Broad Spectrum True SAE (FAILED)**: Using all 32k features blindly yields random noise (AUROC 0.555).

## 2. The Validated Hierarchy

| Component | Status | AUROC | Role |
| :--- | :--- | :--- | :--- |
| **1. MFAP4 Biomarker** | 🏆 **VALIDATED** | **0.763** | Intrinsic Resistance (Layer 1) |
| **2. Proxy SAE** | ✅ **Active** | **0.628** | Baseline / Trial Matching |
| **3. Diamond True SAE** | ❌ **FAILED** | 0.555 | *Research Only (Deprecated)* |
| **4. Broad True SAE** | ❌ **FAILED** | 0.555 | Deprecated |

## 3. The 9 Diamond Features (Debunked)
*Source: Audit of Tier 3 Cohort (N=149)*

**CRITICAL WARNING**: These features failed nested cross-validation (Mean AUROC 0.46-0.55). They appeared significant only due to overfitting on a single fold.

| Feature Index | Effect Size (d) | Fold 4 AUROC | Mean AUROC | Status |
| :--- | :--- | :--- | :--- | :--- |
| **27607** | +0.635 | 0.78 | 0.55 | ❌ Spurious |
| **26220** | +0.609 | 0.78 | 0.55 | ❌ Spurious |
| **16337** | +0.634 | 0.78 | 0.55 | ❌ Spurious |
| **12893** | +0.597 | 0.78 | 0.55 | ❌ Spurious |
| **6020** | +0.573 | 0.78 | 0.55 | ❌ Spurious |
| **22868** | +0.544 | 0.78 | 0.55 | ❌ Spurious |
| **1407** | +0.537 | 0.78 | 0.46 | ❌ Spurious |
| **31362** | +0.517 | 0.78 | 0.55 | ❌ Spurious |
| **9738** | +0.530 | 0.78 | 0.55 | ❌ Spurious |

**Action**: Do NOT use these for clinical decision making.

## 4. The 4-Layer Resistance Framework

This framework remains valid, but Layer 2 is now powered primarily by **Proxy SAE**, not Diamond SAE.

1.  **Layer 1 (Intrinsic)**: MFAP4 / EMT State.
2.  **Layer 2 (Genetic)**: **Proxy SAE** (Gene Mutations -> Pathway Burden).
3.  **Layer 3 (Adaptive)**: CA-125 Kinetics (Serial Delta).
4.  **Layer 4 (Clearance)**: ABCB1 / Drug Efflux.

## 5. Strategic Directives
1.  **HALT** all "True SAE" integration into clinical UI. It is a "Show, Don't Tell" feature only.
2.  **PIVOT** all validation efforts to MFAP4 (Layer 1) and Serial Kinetics (Layer 3).
3.  **USE** Proxy SAE for trial matching (Mechanism-based), as it is logically sound even if predictive power is moderate.
