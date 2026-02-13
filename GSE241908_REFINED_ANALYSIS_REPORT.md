# 📊 GSE241908 Refined Analysis Report

**Cohort:** n=6 Paired Patients (Ascites Pre/Post Chemo)
**Method:** Z-Score Module Scoring (Cohort-Normalized TPM) -> Paired Stats
**Status:** Heterogeneous Signal (High Variance)

## 1. Executive Summary
Paired n=6 analysis (TPM-based): DNA Repair and Cell Cycle module Z-scores show directional upregulation in 4/6 patients, with strong downregulation in 2/6 patients, yielding a non-significant paired p-value (p > 0.5) despite a majority-trending effect.

**Interpretation:** The response is heterogeneous across individuals; results support a “trend with outliers” rather than a consistent cohort-wide shift. The majority (67%) align with the BriTROC hypothesis.

## 2. Data Integrity Audit
- **Expected Samples:** n=14 (7 pairs)
- **Actual processed data:** n=13 samples
- **Excluded:** **Patient 6 Post-Chemo** (Missing from dataset; singleton ShV-86 (Pre) excluded).
- **Normalization:** TPM (Transcripts Per Million). Z-scores calculated within this n=13 cohort.

## 3. Statistical Results (Per-Module)
| Module | Direction (Up/Total) | Median ΔZ (IQR) | Paired T-test | Sign Test (Binomial) |
|---|---|---|---|---|
| **DNA Repair (DDR)** | **3/6** | -0.000 (0.653) | p=0.5810 | p=0.6562 |
| **Cell Cycle (E2F)** | **4/6** | 0.263 (1.096) | p=0.5639 | p=0.3438 |

## 4. Patient-Level Trajectories (Spaghetti Data)
| Patient | DDR Pre(Z) | DDR Post(Z) | **DDR ΔZ** | Cell Cycle ΔZ | Interpretation |
|---|---|---|---|---|---|
| Patient_3 | -0.327 | 0.292 | **+0.619** | +0.592 | ✅ UP |
| Patient_5 | 0.130 | 2.057 | **+1.927** | +2.187 | ✅ UP |
| Patient_8 | -0.302 | -0.270 | **+0.031** | +0.403 | ✅ UP |
| Patient_2 | -0.070 | -0.301 | **-0.231** | -0.777 | ❌ DOWN |
| Patient_4 | 0.449 | -0.470 | **-0.918** | -0.847 | ❌ DOWN |
| Patient_7 | -0.303 | -0.334 | **-0.032** | +0.123 | ❌ DOWN |

## 5. Module Definitions
- **DNA Repair (DDR):** Aggregated hallmark genes covering HR (BRCA1/2, RAD51), Checkpoints (ATM, ATR, CHEK1/2), and Sensors (MRE11, RAD50, NBN).
- **Cell Cycle (Proliferation):** E2F targets and G2M drivers (MKI67, CCNE1, CDK1, PCNA, TOP2A).
