# SAE Biomarker Discovery Protocol

## 1. System Overview
This protocol defines the methodology for identifying SAE (Sparse Autoencoder) features from the Evo2 DNA language model that correlate with clinical outcomes (specifically platinum chemotherapy response in Ovarian Cancer).

**Objective**: Find "monosemantic" biological signals in the 32,768-dimensional SAE feature space that predict resistance better than gene-level markers.

## 2. Methodology

### Phase 1: Feature Extraction
*   **Model**: Evo2 (7B parameter base model).
*   **Input**: ±4096bp context window around patient variants.
*   **Layer**: Layer 26 activations (optimal for biological features).
*   **SAE Transform**: 1920-dim Evo2 activations -> 32,768-dim sparse features.
*   **Sparsity**: Top-k (k=64) features kept per position.

### Phase 2: Statistical Correlation
For each of the 32,768 features, we compute:
1.  **Pearson Correlation**: Linear relationship with outcome (Sensitive=1.0, Resistant=0.5, Refractory=0.0).
2.  **Effect Size (Cohen's d)**: Magnitude of difference between resistant/sensitive populations.
3.  **Cross-Validation**: Stability of signal across data splits.
4.  **FDR Correction**: Benjamini-Hochberg correction for multiple testing.

## 3. Execution Pipeline

### Scripts
*   **Extraction**: `scripts/sae/extract_sae_features_cohort.py`
*   **Analysis**: `scripts/sae/analyze_biomarkers.py`
*   **Service**: `api/services/biomarker_correlation_service.py`

### Run Command
```bash
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots
```

## 4. Interpretation Standards

### Thresholds for "Diamond" Features
*   **P-value**: < 0.01
*   **Cohen's d**: > 0.3
*   **Stability**: Selected in > 60% of CV folds

### Clinical Mapping (Future State)
Biomarkers are mapped to drug scoring (WIWFM) via:
*   **Positive Correlation**: Boosts drug efficacy score.
*   **Negative Correlation**: Penalizes score (indicates resistance).
*   **Mechanism**: Mapped to biological interpretation (e.g., DNA repair capacity) where possible.

## 5. Current Status
*   **Extraction**: Completed for sub-cohort (69 patients, ~3k variants).
*   **Analysis**: Statistical engine ready.
*   **Next Step**: Full cohort execution and "Diamond" mining (see Postmortem).
