# GSE165897 Processing Status

**Date**: 2026-01-13  
**Dataset**: GSE165897 (GEO)  
**Status**: ✅ **Data Downloaded and Processed**

---

## Dataset Overview

- **Access**: GEO (public, no application required)
- **Cohort**: 11 HGSOC patients, paired pre/post-NACT scRNA-seq
- **Sample Size**: 22 scRNA-seq samples (11 paired)
- **Treatment**: Platinum-based neoadjuvant chemotherapy (NACT)
- **SAE Suitability Score**: 8.7/10

---

## Processing Status

### ✅ Completed

1. **Data Download**
   - Location: `scripts/data_acquisition/sae/`
   - Files:
     - `GSE165897_cellInfo_HGSOC.tsv.gz` (cell metadata)
     - `GSE165897_UMIcounts_HGSOC.tsv.gz` (expression matrix)
     - Series matrix files (GPL16791, GPL20301, GPL24676)

2. **Pathway Kinetics Computation**
   - Script: `scripts/serial_sae/pathway_kinetics_gse165897.py`
   - Results: `data/serial_sae/gse165897/results/`
   - Pathways analyzed: DDR, MAPK, PI3K, VEGF
   - Method: Pseudo-bulk aggregation → pathway scores → kinetics (Δ = post - pre)

3. **Analysis Outputs**
   - `pathway_scores.csv`: Pre/post pathway scores per patient
   - `pathway_kinetics.csv`: Pathway kinetics (Δ values)
   - `kinetic_patterns.csv`: Pattern classifications
   - `pathway_correlations.csv`: Correlation matrix
   - `pathway_kinetics_heatmap.png`: Visualization
   - `pathway_kinetics_report.txt`: Summary report

---

## Key Findings

### Pathway Kinetics Summary (n=11 patients)

| Pathway | Mean Δ | Median Δ | Increases (Δ > 0.05) | Decreases (Δ < -0.05) |
|---------|--------|----------|---------------------|----------------------|
| **DDR** | -0.0140 | -0.0129 | 0/11 | 0/11 |
| **MAPK** | 0.0002 | -0.0007 | 0/11 | 1/11 |
| **PI3K** | -0.0097 | -0.0139 | 1/11 | 0/11 |
| **VEGF** | 0.0268 | 0.0125 | **4/11** | 1/11 |

### Key Observations

1. **VEGF Activation**: 4/11 patients (36%) show VEGF pathway increase (Δ > 0.05) post-NACT
   - Suggests angiogenesis activation as a potential resistance mechanism
   - Largest increase: EOC349 (Δ = 0.122)

2. **DDR Suppression**: Mean DDR decrease (-0.014) suggests DNA repair downregulation post-treatment
   - Largest decrease: EOC1005 (Δ = -0.048)

3. **Pathway Correlations**: 
   - DDR ↔ MAPK: r=0.492
   - DDR ↔ VEGF: r=0.479
   - PI3K ↔ VEGF: r=0.512

---

## Missing Data

### ⚠️ Treatment Response Labels

**Status**: Not available in GEO metadata

**What's Missing**:
- RECIST criteria (complete/partial/stable/progressive disease)
- Platinum sensitivity classification
- Progression-free survival (PFS)
- Overall survival (OS)

**Action Required**:
1. Extract treatment response from publication (Zhang et al., 2022)
2. Create `resistance_labels.json` or `resistance_labels.csv` mapping patient_id → response
3. Re-run analysis with resistance stratification

**Expected Format**:
```json
{
  "EOC1005": "resistant",
  "EOC136": "sensitive",
  ...
}
```

or CSV:
```csv
patient_id,resistance_label
EOC1005,resistant
EOC136,sensitive
...
```

---

## Next Steps

### Immediate

1. **Extract Treatment Response Labels**
   - Search publication (Zhang et al., 2022) for RECIST criteria
   - Map patient IDs to treatment response
   - Create resistance labels file

2. **Resistance-Stratified Analysis**
   - Correlate pathway kinetics with treatment response
   - Identify predictive patterns in pre-treatment scores
   - Compare responders vs non-responders

### Short-term

3. **Integration with TCGA-OV**
   - Compare pathway kinetics patterns
   - Validate findings across datasets
   - Identify consistent resistance mechanisms

4. **Publication Analysis**
   - Review original publication for additional insights
   - Extract clinical outcomes if available
   - Cross-reference with pathway findings

---

## Files and Locations

### Data Files
- `scripts/data_acquisition/sae/GSE165897_*.tsv.gz`
- `scripts/data_acquisition/sae/GSE165897-*_series_matrix.txt.gz`

### Results
- `data/serial_sae/gse165897/results/pathway_kinetics.csv`
- `data/serial_sae/gse165897/results/pathway_scores.csv`
- `data/serial_sae/gse165897/results/pathway_kinetics_report.txt`
- `data/serial_sae/gse165897/results/gse165897_analysis_summary.txt`

### Scripts
- `scripts/serial_sae/pathway_kinetics_gse165897.py` (main processing)
- `scripts/serial_sae/gse165897_summary.py` (summary generation)

---

## Integration with Mission

**Phase 2 Status**: ✅ Complete
- Systematic search: ✅ Complete
- Top datasets identified: MSK_SPECTRUM (9.5/10), GSE165897 (8.7/10)
- GSE165897: ✅ Downloaded and processed
- TCGA-OV: ✅ Correlated with clinical outcomes (426 patients)

**Mission Status**: Phase 2 complete. GSE165897 processed (11 patients, pathway kinetics computed). TCGA-OV correlation complete (426 patients).

---

## References

- **GEO Accession**: GSE165897
- **Publication**: Zhang et al., 2022 (need to extract treatment response)
- **Platform**: 10x Genomics scRNA-seq
- **Analysis Pipeline**: Pseudo-bulk aggregation → pathway scoring → kinetics computation

'll create a comprehensive markdown documentation file for your dataset acquisition process. Here's the complete documentation:

text
# Dataset Acquisition Documentation
**CrisPRO.ai Precision Oncology Project**  
**Date Created:** January 13, 2026  
**Last Updated:** January 13, 2026  
**Author:** Agent D (Dataset Acquisition)

---

## Project Overview

**Mission:** Acquire multi-modal public datasets for ovarian cancer resistance prediction modeling, combining:
- TIL scores
- TCR sequencing (diversity metrics)
- RNA-seq (IO pathway analysis using Sparse Autoencoders)
- PD-L1 status
- IO treatment outcomes (RESPONSE DATA)

**Target Hypothesis:** Combined AUROC >0.75 vs. PD-L1 alone (~0.60-0.65)

**Cancer Types:** Melanoma or NSCLC (with checkpoint inhibitors: PD-1/PD-L1/CTLA-4)

**Scaffolding Data:** TCGA ovarian + MSK SPECTRUM (already validated, NO IO treatment)

---

## Acquired Datasets

### 1. GSE91061 - Melanoma Anti-PD-1 Treatment (✅ BREAKTHROUGH DATA)

**Study:** Riaz et al. (2017) - "Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab"  
**Journal:** Cell 171(4): 934-949.e16  
**DOI:** https://doi.org/10.1016/j.cell.2017.09.028

**Data Type:** RNA-seq (bulk) - Pre/On-treatment paired samples  
**Platform:** Illumina HiSeq 2500  
**Samples:** 110 total
- 51 pre-treatment baseline
- 59 on-treatment biopsies
- **Paired samples:** Available for subset of patients

**Clinical Data:**
- ✅ **Response labels:** Complete/Partial Response, Stable/Progressive Disease
- ✅ **PFS/OS data:** Available
- ✅ **Treatment:** Nivolumab (anti-PD-1)
- ✅ **Cancer type:** Melanoma

**Access:** Public (GEO)

**Download Instructions:**
```bash
# Method 1: GEO FTP (Raw FASTQ files)
wget -r -np -nH --cut-dirs=3 \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/

# Method 2: Processed count matrices
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.scaled.csv.gz

# Method 3: SRA Toolkit
prefetch --option-file SRR_Acc_List.txt
fasterq-dump SRR*
Priority Files:

GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz - Raw count matrix

GSE91061_BMS038109Sample.hg19KnownGene.scaled.csv.gz - Normalized counts

GSE91061_BMS038109Sample.txt.gz - Sample metadata/clinical annotations

Individual FASTQ files (SRA accessions available in Series Matrix)

Key Features:

✅ IO treatment (nivolumab)

✅ RNA-seq for SAE analysis

✅ Response/survival data

⚠️ No TCR-seq (bulk RNA can estimate repertoire)

⚠️ No explicit TIL scores (can estimate from deconvolution)

⚠️ PD-L1 data needs extraction from publication supplements

GEO Accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061

2. GSE78220 - Melanoma Anti-PD-1 Treatment (✅ BREAKTHROUGH DATA)
Study: Hugo et al. (2016) - "Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma"
Journal: Cell 165(1): 35-44
DOI: https://doi.org/10.1016/j.cell.2016.02.065

Data Type: RNA-seq (bulk) - Pre-treatment baseline samples
Platform: Illumina HiSeq 2000
Samples: 38 pre-treatment biopsies

Clinical Data:

✅ Response labels: Responder vs. Non-responder

✅ Treatment: Pembrolizumab (anti-PD-1)

✅ Cancer type: Metastatic melanoma

✅ PD-L1 IHC data: Available in publication

Access: Public (GEO)

Download Instructions:

bash
# Processed expression data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientID_vs_GeneSymbol_HUGO_VoomCounts.txt.gz

# Raw FASTQ via SRA
prefetch --option-file SRR_Acc_List.txt
Key Features:

✅ IO treatment (pembrolizumab)

✅ RNA-seq for SAE analysis

✅ Response data

✅ PD-L1 status available

✅ Gene expression signatures for TME/TIL

⚠️ No TCR-seq

⚠️ Pre-treatment only (no on-treatment pairs)

GEO Accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220

3. GSE165897 - Ovarian Cancer Scaffolding Data (✅ SCAFFOLDING)
Study: Zhang et al. (2022) - "Longitudinal single-cell RNA-seq analysis reveals stress-promoted chemoresistance in metastatic ovarian cancer"
Journal: Science Advances 8(8): eabm1831
DOI: https://doi.org/10.1126/sciadv.abm1831

Data Type: scRNA-seq (10x Genomics Chromium) + Bulk RNA-seq
Platform: Illumina HiSeq 4000, HiSeq 2500, NovaSeq 6000
Samples: 22 paired treatment-naïve and post-NACT samples (11 patients)

Clinical Data - ALL 11 scRNA-seq Patients:

Patient ID	PFI (days)	Resistance Label	CRS	Stage	Age	CA125 TN	CA125 PN
EOC1005	65	Resistant	2	IVA	73	3776	343
EOC136	520	Sensitive	2	IVA	64	2647	212
EOC153	393	Sensitive	2	IVA	78	1063	93
EOC227	230	Resistant	2	IVA	74	445	33
EOC3	14	Resistant	2	IVA	67	821	221
EOC349	36	Resistant	2	IVB	67	2155	67
EOC372	460	Sensitive	1	IIIC	68	3180	334
EOC443	177	Resistant	3	IVA	54	2295	82
EOC540	126	Resistant	2	IIIC	62	155	7
EOC733	83	Resistant	1	IVA	72	22079	3579
EOC87	30	Resistant	1	IIIC	62	998	346
Resistance Classification: PFI <180 days (6 months) = Resistant; ≥180 days = Sensitive
Result: 8 resistant, 3 sensitive

Treatment: NACT (neoadjuvant chemotherapy) - platinum-taxane, NO immunotherapy
Median PFI: 4.2 months (126 days)

Access:

Public: GEO GSE165897 (count matrices, metadata)

Controlled: EGA EGAS00001005010 (raw FASTQ files)

Download Instructions:

bash
# Public metadata and processed counts
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165897/suppl/GSE165897_cellinfo_HGSOC.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165897/matrix/

# GenomeSpy visualization tool (for visual exploration)
# https://csbi.ltdk.helsinki.fi/pub/projects/lahtinen_et_al_2023/

# For raw FASTQ files - requires EGA account and DAC approval
# Submit request at: https://ega-archive.org/datasets/EGAD00001005010
Key Features:

✅ scRNA-seq (93,650 cells QC'd → 51,786 final)

✅ Longitudinal pre/post-chemo pairs

✅ PFI/resistance labels

✅ Stress-associated cell states identified

✅ TME characterization (CAFs, immune cells)

❌ NO IO treatment (platinum-based chemo only)

❌ Ovarian cancer (not melanoma/NSCLC)

Use Case: Scaffolding data for methodological development, NOT for Agent D breakthrough hypothesis

GEO Accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897

4. IMvigor210 - Urothelial Cancer Anti-PD-L1 (✅ BREAKTHROUGH DATA - Alternative)
Study: Mariathasan et al. (2018) - "TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells"
Journal: Nature 554(7693): 544-548
DOI: https://doi.org/10.1038/nature25501

Data Type: RNA-seq (bulk) - Pre-treatment samples
Platform: Illumina HiSeq
Samples: 348 patients

Clinical Data:

✅ Response labels: CR/PR/SD/PD (RECIST)

✅ PFS/OS data: Available

✅ Treatment: Atezolizumab (anti-PD-L1)

✅ PD-L1 IHC: IC0/IC1/IC2+ scoring

✅ TME signatures: TGFβ, T-effector, angiogenesis

⚠️ Cancer type: Urothelial (not melanoma/NSCLC, but IO-treated)

Access: Requires registration at http://research-pub.gene.com/IMvigor210CoreBiologies/

Download Instructions:

text
# Install package
install.packages("http://research-pub.gene.com/IMvigor210CoreBiologies/
  packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz", 
  repos = NULL, type = "source")

# Load data
library(IMvigor210CoreBiologies)
data(cds) # Expression data
Key Features:

✅ IO treatment (atezolizumab)

✅ RNA-seq for SAE analysis

✅ Response/survival data

✅ PD-L1 status

✅ Large cohort (n=348)

⚠️ Urothelial cancer (different biology than melanoma/NSCLC)

⚠️ Pre-treatment only

Controlled-Access Datasets (Pending Approval)
5. BriTROC Shallow WGS - Ovarian Cancer Genomics (PENDING)
Dataset ID: EGAD00001011049
Study: UK Translational Research in Ovarian Cancer (BriTROC) Consortium
Publication: Piskorz, Ennis, Macintyre et al. (Ann. Oncol. 2015)

Data Type: Shallow whole genome sequencing (~0.1-0.5X coverage)
Platform: Illumina HiSeq 2500, HiSeq 4000
Samples: 679 total (diagnosis + post-relapse pairs)

Clinical Data:

✅ Longitudinal design: Diagnosis → Relapse pairs

✅ CNA profiling: Genome-wide copy number alterations

✅ Treatment: Platinum-based chemotherapy

✅ Cancer type: High-grade serous ovarian cancer (HGSOC)

❌ NO IO treatment

Access: Requires EGA DAC approval (EGAC00001000388)

Application Status: Template prepared (see below)

Request Template:

text
Request for Research Access - BriTROC Shallow WGS Data (EGAD00001011049)

Principal Investigator: [Your Name]
Institution: [Your Institution]

Research Project: Multi-modal biomarker development for platinum resistance 
prediction in high-grade serous ovarian cancer (HGSOC)

Purpose:
1. Copy number alteration (CNA) profiling from shallow WGS
2. Clonal evolution analysis (diagnosis vs. post-relapse)
3. Integration with transcriptomic data (TCGA, HERCULES, GSE165897)

Scientific Justification:
Longitudinal design critical for acquired resistance mechanisms. 
Shallow WGS provides cost-effective genome-wide CNA detection.

Data Usage Commitment:
- Solely for stated research purpose
- No re-identification attempts
- No data sharing outside User Institution
- Cite Piskorz, Ennis, Macintyre et al. (Ann. Oncol. 2015)
- Acknowledge BriTROC/Ovarian Cancer Action/Cancer Research UK
- Destroy raw data upon project completion

Project Duration: [Timeline]
EGA Portal: https://ega-archive.org/datasets/EGAD00001011049

Public Reference Datasets (Already Validated)
TCGA Ovarian Cancer (TCGA-OV)
Samples: 271 patients (Grade G2-G4, Stage IIIA-IV)
Data Types:

RNA-seq (RSEM normalized)

Clinical annotations (PFS, stage, grade)

Reverse phase protein array (RPPA)

Somatic mutations

Access: Broad Firehose (https://gdac.broadinstitute.org/)

Download:

