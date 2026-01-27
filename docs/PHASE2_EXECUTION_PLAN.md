# Phase 2 Execution Plan: Serial SAE Dataset Analysis

**Date:** January 13, 2026  
**Status:** ✅ **DATASETS FOUND** - Ready for execution  
**Total Paired Samples:** 68 immediate (57 + 11) + 276 pending (BriTROC-1)

---

## 🎯 IMMEDIATE EXECUTION (This Week)

### Dataset 0: TCGA-OV Direct Download (OPTIONAL - Baseline Validation)
**Access:** ✅ Public (GDC Data Portal)

**Use Case:** Baseline SAE computation for ~300 primary tumors (not serial monitoring)

**Steps:**
1. **Download Data** (1-2 hours)
   ```bash
   python3 scripts/serial_sae/download_tcga_ov_gdc.py
   gdc-client download -m data/serial_sae/tcga_ov/tcga_ov_rnaseq_manifest.csv
   ```

2. **Process RNA-seq** (2-4 hours)
   - Normalize expression data
   - Compute baseline SAE pathway scores
   - Match with clinical outcomes

**Value:** Large cohort for baseline validation, complements serial monitoring

**Timeline:** 3-6 hours total

**Note:** This is for baseline validation, not serial monitoring (TCGA has mostly primary-only samples)

---

### Dataset 1: cBioPortal MSK-SPECTRUM (40 Paired Patients)

**Access:** ✅ Public (cBioPortal API)

**Steps:**
1. **Download Data** (15-30 min manual + 5 min processing)
   ```bash
   # Manual download via cBioPortal web interface
   # Study ID: msk_spectrum_tme_2022
   # Extract: 40 patients with paired primary+recurrent samples
   # See: scripts/serial_sae/MANUAL_DOWNLOAD_INSTRUCTIONS.md
   ```

2. **Process Downloaded Files** (5 minutes)
   - Run: `python3 scripts/serial_sae/process_manual_downloads.py`
   - Extracts mutations/expression for 40 paired patients
   - Saves processed data for SAE computation

3. **Compute SAE Pathway Scores** (2-3 hours)
   - Expression-based: Use mRNA expression data (Z-scores, FPKM, TPM)
   - Mutation-based: Use mutation calls for pathway scoring
   - Compute 7D mechanism vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]

4. **Calculate Pathway Kinetics** (1 hour)
   - ΔDDR_BIN = (DDR_T2 - DDR_T1) / time_interval
   - ΔMAPK_BIN = (MAPK_T2 - MAPK_T1) / time_interval
   - ΔPI3K_BIN = (PI3K_T2 - PI3K_T1) / time_interval

5. **Correlate with Outcomes** (1-2 hours)
   - Time to progression (TTP)
   - Progression-free survival (PFS)
   - Overall survival (OS)
   - Response status

**Deliverable:** 
- `data/serial_sae/cbioportal_tcga_msk_paired_analysis.json`
- Pathway kinetics per patient
- Correlation results (ΔSAE vs outcomes)

**Timeline:** 1-2 hours total (manual download + processing + SAE computation)

---

### Dataset 4: GSE165897 (11 Paired Patients, scRNA-seq)

**Access:** ✅ Public (GEO)

**Steps:**
1. **Download Data** (30 min - 2 hours)
   ```bash
   # Automated: python3 scripts/serial_sae/download_gse165897.py
   # Or manual: Download from GEO (GSE165897)
   # Format: 10x Genomics scRNA-seq
   ```

2. **Process scRNA-seq** (1-2 days)
   - Quality control, normalization
   - Cell type annotation
   - Pathway score computation at single-cell level

3. **Aggregate by Patient/Timepoint** (4-6 hours)
   - Compute pathway scores per sample
   - Compare treatment-naïve vs post-NACT
   - Identify resistant subpopulations

4. **Correlate with Outcomes** (2-3 hours)
   - PFS (stress-high: 14.9 mo vs stress-low: 21.2 mo, P=0.0037)
   - Chemotherapy response
   - Clonal structure changes

**Deliverable:**
- `data/serial_sae/gse165897_scrnaseq_analysis.json`
- Single-cell pathway scores
- Aggregated patient-level pathway kinetics
- Correlation results

**Timeline:** 2-3 days total (scRNA-seq processing is compute-intensive)

---

## 📊 COMBINED ANALYSIS (Week 1)

**Total Samples:** 51 paired patients (40 + 11) for serial monitoring  
**Plus:** ~300 TCGA-OV primary tumors for baseline validation

**Analysis Plan:**
1. **Combine Results** (2 hours)
   - Merge cBioPortal + GSE165897 findings
   - Standardize pathway score computation
   - Create unified dataset

2. **Meta-Analysis** (4-6 hours)
   - Pooled correlation: ΔSAE vs TTP/PFS
   - Stratified by treatment type
   - Sensitivity analysis

3. **Resistance Prediction Model** (4-6 hours)
   - Logistic regression: ΔSAE → resistance probability
   - Threshold optimization
   - Cross-validation

**Deliverable:**
- `data/serial_sae/combined_analysis.json`
- `receipts/serial_sae_validation_pilot.json`
- Figures: Pathway kinetics plots, correlation plots, ROC curves

**Timeline:** 1-2 days

---

## 🎯 MEDIUM-TERM (2-4 Weeks)

### Dataset 2: BriTROC-1 (276 Paired Patients)

**Access:** ⚠️ Controlled (EGA - requires approval)

**Steps:**
1. **Submit EGA Access Request** (Day 1)
   - Institutional approval required
   - Data access committee review
   - Timeline: 2-4 weeks

2. **Download Data** (1-2 days, after approval)
   - Accession: EGAS00001007292
   - Download: sWGS + targeted panel data

3. **Compute SAE** (8-12 hours)
   - Mutation-based pathway scores (targeted panel)
   - Copy number-based pathway scores (sWGS)
   - Combine for comprehensive pathway assessment

4. **Validate on Large Cohort** (2-3 days)
   - 276 paired patients (largest paired HGSOC cohort)
   - Platinum sensitivity correlation
   - HRD status validation

**Deliverable:**
- `data/serial_sae/britroc1_analysis.json`
- Large-cohort validation results
- Publication-ready findings

**Timeline:** 2-4 weeks (approval) + 1 week (analysis)

---

## 📧 LONG-TERM (3-6 Months)

### Dataset 3: Williams et al. Nature 2025

**Access:** ❓ Unknown (data embargoed)

**Action Items:**
1. **Contact Authors** (This Week)
   - Email corresponding author
   - Inquire about data availability timeline
   - Request pre-publication access (if possible)

2. **Monitor Repositories** (Ongoing)
   - Check EGA for deposition
   - Check dbGaP for deposition
   - Monitor paper updates

3. **Prepare Analysis Pipeline** (While Waiting)
   - Design single-cell SAE computation
   - Prepare clonal evolution analysis
   - Test on GSE165897 (similar data type)

**Value:** Exceptional granularity (5-20 samples per patient, clonal dynamics)

---

## 📋 PRIORITY RANKING

### Tier 1: Execute Now ✅
1. **cBioPortal TCGA+MSK** - 57 paired, open access, full SAE compatible
2. **GSE165897** - 11 paired, open access, single-cell resolution

**Total:** 68 paired samples, immediate execution

### Tier 2: Submit Access Request ⏳
3. **BriTROC-1** - 276 paired, controlled access, 2-4 week approval

**Value:** Largest paired HGSOC cohort for validation

### Tier 3: Contact Authors 📧
4. **Williams et al.** - 18 patients (dense sampling), data embargoed

**Value:** Exceptional granularity, clonal dynamics

---

## ✅ SUCCESS METRICS

**Minimum Viable (Week 1):**
- ✅ Download 51 paired samples (40 cBioPortal + 11 GSE165897)
- ✅ Compute SAE pathway scores at both timepoints
- ✅ Calculate pathway kinetics (ΔSAE)
- ✅ Correlate with outcomes (TTP, PFS, response)
- ✅ Proof-of-concept: ΔSAE predicts resistance

**Ideal (Month 1):**
- ✅ Validate on 276 additional patients (BriTROC-1)
- ✅ Establish resistance prediction thresholds
- ✅ Publication-ready findings

---

## 🚀 IMMEDIATE NEXT STEPS

1. **Download cBioPortal Data** (Today - Manual)
   - Instructions: `scripts/serial_sae/MANUAL_DOWNLOAD_INSTRUCTIONS.md`
   - Process: `scripts/serial_sae/process_manual_downloads.py`
   - Extract 40 paired patients
   - Save to `data/serial_sae/cbioportal_paired/`

2. **Download GSE165897** (Today)
   - Script: `scripts/serial_sae/download_gse165897.py`
   - Download scRNA-seq data
   - Save to `data/serial_sae/gse165897/`

3. **Compute SAE Scores** (This Week)
   - Script: `scripts/serial_sae/compute_serial_sae.py`
   - Compute pathway scores at both timepoints
   - Calculate pathway kinetics

4. **Correlate with Outcomes** (This Week)
   - Script: `scripts/serial_sae/correlate_delta_sae_outcomes.py`
   - ΔSAE vs TTP/PFS/response
   - Generate figures

5. **Submit BriTROC-1 Request** (This Week)
   - Prepare EGA access application
   - Submit to data access committee

---

**Status:** ✅ **READY FOR EXECUTION**  
**Next:** Download datasets and begin analysis
