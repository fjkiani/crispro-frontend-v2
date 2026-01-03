# Supplementary Materials (Draft)

## Supplementary Table S1 — Cohort characteristics (GSE63885)
- Total samples: 101
- Labeled platinum sensitivity: 101
  - Resistant: 34
  - Sensitive: 67
- Platform: GPL570

## Supplementary Table S2 — Probe mapping used (GPL570)
- MFAP4: 212713_at
- EFEMP1: 201842_s_at; 201843_s_at; 228421_s_at
- VIM: 1555938_x_at; 201426_s_at
- CDH1: 201130_s_at; 201131_s_at
- SNAI1: 219480_at

## Supplementary Methods — Data processing details

### Data download
- Series matrix: `GSE63885_series_matrix.txt.gz` downloaded from GEO and decompressed.
- Platform annotation: `GPL570.annot.gz` downloaded from GEO FTP and decompressed.

### Label extraction
Labels were parsed from GEO `!Sample_characteristics_ch1` fields:
- `platinium sensitivity (resistant: DFS<180 days; moderately sensitive: 180>DFS>732; highly sensitive: DFS>732 days)`

### Feature construction
- Probe → gene: if multiple probes per gene, compute mean across probes.
- Standardize each gene across samples (z-score).

### Models
- Univariate: AUROC for MFAP4.
- Multivariate: logistic regression on standardized EMT features.

## Supplementary Code & Reproducibility

### Data artifacts
`oncology-coPilot/oncology-backend-minimal/data/external/GSE63885/`
- `GSE63885_series_matrix.txt`
- `sample_annotations.csv`
- `GPL570.annot.txt`
- `emt_probe_mapping.json`
- `emt_platinum_auroc_results.json`

### Figure generation
- Script: `oncology-coPilot/oncology-backend-minimal/scripts/figures/generate_gse63885_figures.py`
- Output directory: `.cursor/MOAT/SAE_INTELLIGENCE/Publication-1/SAE_RESISTANCE/figures/`

### Key figure outputs (this submission)
- `fig_gse63885_roc_mfap4.(png|pdf)`
- `fig_gse63885_roc_emt_score.(png|pdf)`
- `fig_gse63885_box_mfap4_by_platinum.(png|pdf)`
- `fig_gse63885_box_emt_by_platinum.(png|pdf)`
- `fig_gse63885_cohort_flow.(png|pdf)`
