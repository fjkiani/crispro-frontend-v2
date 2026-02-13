# Phase 2: Serial SAE Pathway Kinetics Validation — Single Source of Truth

**Last Updated:** January 13, 2026  
**Status:** Phase 2 complete — datasets identified, GSE165897 + TCGA-OV analyzed  
**Consolidated From:** PHASE2_AGENT_INSTRUCTIONS, PHASE2_DATA_HUNT_SPECIFICATION, PHASE2_EXECUTION_PLAN, PHASE2_EXECUTION_STATUS, PHASE2_FINAL_STRATEGY, PHASE2_FINDINGS_SUMMARY, PHASE2_QUICK_REFERENCE, PHASE2_SERIAL_SAE_VALIDATION_MASTER, PHASE2_SYSTEMATIC_SEARCH_RESULTS

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Mission & Objective](#2-mission--objective)
3. [Agent Instructions & Search Strategy](#3-agent-instructions--search-strategy)
4. [Data Hunt Specification](#4-data-hunt-specification)
5. [Systematic Search Results](#5-systematic-search-results)
6. [Dataset Discovery & Evaluation](#6-dataset-discovery--evaluation)
7. [Execution Status & Final Strategy](#7-execution-status--final-strategy)
8. [GSE165897 Analysis](#8-gse165897-analysis)
9. [TCGA-OV Validation](#9-tcga-ov-validation)
10. [Code Audit & Gaps](#10-code-audit--gaps)
11. [Next Steps & Action Items](#11-next-steps--action-items)
12. [Quick Reference](#12-quick-reference)
13. [Files & Locations](#13-files--locations)
14. [References](#14-references)

---

## 1. Executive Summary

### Mission Objective
Validate SAE pathway kinetics (DDR, MAPK, PI3K, VEGF changes from pre → post-treatment) for mechanism-specific chemotherapy resistance prediction using paired longitudinal samples from ovarian cancer patients.

### Current Status
- **Dataset Discovery:** Complete — 4 primary datasets + TCGA-OV baseline; 362 paired samples total identified.
- **GSE165897 Processing:** Complete — 11 paired patients processed, pathway kinetics computed.
- **TCGA-OV Baseline:** Complete — 426 patients, baseline SAE and outcome correlations.
- **Resistance Analysis:** Partial — 3/11 GSE165897 patients with known PFI; VEGF activation in resistant vs sensitive noted.
- **Code Implementation:** Gaps — expression-based scoring and VEGF pathway need to be added where missing.

### Key Finding
**VEGF pathway activation in resistant patients** (+0.0694 Δ) versus suppression in sensitive patients (-0.0323 Δ), suggesting angiogenesis as a chemotherapy resistance mechanism (n=3; full cohort PFI pending).

### Decisions
- **TCGA-OV:** Primary focus for baseline validation (complete).
- **MSK-SPECTRUM / cBioPortal TCGA+MSK:** Abandoned — small cohort, no expression via API.
- **GSE165897:** Secondary serial proof-of-concept (complete).
- **BriTROC-1:** Pending EGA access (276 paired).
- **Williams et al.:** Contact authors (18 patients, dense sampling).

---

## 2. Mission & Objective

**Find datasets with paired samples (baseline + progression) that allow:**
1. Compute SAE pathway scores at multiple timepoints  
2. Correlate pathway changes (ΔSAE) with clinical outcomes  
3. Validate serial monitoring hypothesis  

**Minimum viable:** n ≥ 20 patients with paired samples (ideal: 3–5 datasets, mix ovarian + other cancers, public where possible).

---

## 3. Agent Instructions & Search Strategy

### Exact Search Strategy

**Step 1: PubMed (30 min)**  
- Query 1 — Ovarian serial: `("ovarian cancer" OR "ovarian carcinoma") AND ("serial sequencing" OR "longitudinal sequencing" OR "paired samples" OR "primary recurrent") AND ("whole exome" OR "WES" OR "RNA-seq")`  
- Query 2 — Ovarian ctDNA + tumor: `("ovarian cancer" OR "ovarian carcinoma") AND ("ctDNA" OR "circulating tumor DNA") AND ("tumor sequencing" OR "tissue sequencing")`  
- Query 3 — Ovarian organoid: `("ovarian cancer" OR "ovarian carcinoma") AND ("organoid" OR "patient-derived organoid") AND ("drug testing" OR "resistance" OR "serial")`  

**Step 2: Specific studies**  
- Goranova et al. (ovarian + ctDNA) — paired tumor?  
- Parkinson et al. (ovarian + serial + plasma) — tumor sequencing?  
- Patch et al. (recurrent ovarian + whole genome) — primary + recurrent pairs?  

**Step 3: Repositories**  
- GEO: "ovarian cancer serial" OR "ovarian cancer longitudinal"; multiple samples per patient.  
- cBioPortal: Ovarian studies with "recurrent" or "progression" samples.  
- dbGaP: Ovarian studies; note serial samples (may require application).  

**Step 4: Document**  
- One section per dataset in deliverable; summary table; priority tiers (Tier 1 = best → Tier 4 = future).

### What to Look For (prioritize)
- Paired samples (baseline + progression or mid-treatment)  
- Genomic data: WES mutations or RNA-seq expression  
- Outcomes: TTP, PFS, response, resistance labels  
- Accessible: public or requestable  
- n ≥ 20 patients with pairs  

### Skip
- Only plasma ctDNA (no tumor sequencing)  
- Single timepoint only  
- &lt;10 paired samples  
- No outcome data  
- Restricted access with no contact  
- Only CNV (no mutations/expression)  

---

## 4. Data Hunt Specification

### Required Information Per Dataset

For each dataset, extract (structure below):

- **Study:** authors, journal, year, PMID, DOI, title  
- **Cohort:** cancer_type, sample_size, serial_timepoints, treatment  
- **Data availability:** status (public/restricted/not available), repository, accession_id, data_types (WES, RNA-seq, WGS, targeted_panel)  
- **SAE compatibility:** can_compute_pathway_scores (yes/no), reasoning, required_data, missing_data, feasibility (high/medium/low)  
- **Outcome data:** available (yes/no), outcomes (OS, PFS, TTP, response, resistance), time_to_progression, resistance_labels  
- **Notes**

### SAE Compatibility Checklist
- [ ] Somatic mutations (WES/panel/MAF) or expression (RNA-seq/array)  
- [ ] Paired samples (baseline + progression)  
- [ ] Outcome data (TTP, PFS, response)  
- [ ] Accessible (public or requestable)  
- [ ] n ≥ 20 for validation  

**Feasibility:** High = mutations + outcomes + pairs + accessible; Medium = expression + outcomes + pairs + accessible; Low = missing critical data or not accessible.

### Reporting Template (per dataset)

```markdown
### Dataset X: [Study Name]
**Citation:** Author et al., Journal Year (PMID: xxxxx)
**Cohort:** Cancer type, n=, timepoints, treatment
**Data Availability:** Status, repository, accession, data types
**SAE Compatibility:** YES/NO, reasoning, feasibility
**Outcome Data:** Available, outcomes listed
**Notes:** ...
```

### Success Criteria
- **Minimum:** 1–2 datasets, n ≥ 20 paired, accessible, SAE computable.  
- **Ideal:** 3–5 datasets, n ≥ 20 paired, mix ovarian + other cancers, publicly accessible.  

---

## 5. Systematic Search Results

**Databases searched:** GENIE-BPC, GEO/SRA, dbGaP, PubMed/PMC, JCI.  
**Datasets evaluated:** 5. **Validation studies reviewed:** 2 (TCGA 2011, Verhaak 2013).

### GENIE-BPC
- **Ovarian cohort:** Not available.  
- **Available:** NSCLC, CRC, Prostate, Breast, Pancreas, Bladder, Melanoma.  
- **Action:** Monitor for future ovarian release.

### Longitudinal Datasets Evaluated (summary)

| Dataset            | n (Paired) | Access     | SAE suitability      | Priority | Status        |
|--------------------|------------|------------|----------------------|----------|---------------|
| MSK_SPECTRUM       | 57/40*     | Public/dbGaP| Full (RNA+Mut)       | HIGH     | Abandoned*    |
| GSE165897          | 11         | Public     | Full (scRNA-seq)     | HIGH     | Complete      |
| GSE241908          | 7          | Public     | Bevacizumab-specific | Secondary| Optional      |
| GSE217177          | —          | Public     | Cross-sectional      | —        | Not suitable  |
| GSE184880          | —          | Public     | Treatment-naive only | —        | Not suitable  |
| BriTROC-1 (EGA)    | 276        | Controlled | Partial (Mut+CN)     | MEDIUM   | Pending       |
| Williams Nature 2025| 18        | Unknown    | Full (scRNA+WGS)     | HIGH*    | Contact authors|
| TCGA-OV (baseline) | ~300       | Public     | Full (RNA-seq)       | HIGH     | Complete      |

*MSK-SPECTRUM abandoned: small cohort, no expression via cBioPortal API.

### TCGA Validation vs Published Studies
- TCGA-OV pathway results aligned with TCGA 2011 (Nature) and Verhaak 2013 (JCI): TP53, RB1, PI3K/RAS, FOXM1, NOTCH, HR deficiency.  
- Mesenchymal subtype: poor prognosis, higher platinum resistance (63% vs 23%), shorter OS.

---

## 6. Dataset Discovery & Evaluation

### Total Paired Samples Identified
- **Immediate (open):** 68 (57 cBioPortal + 11 GSE165897) — cBio arm later abandoned.  
- **Pending:** 294 (276 BriTROC-1 + 18 Williams).  
- **Grand total:** 362 paired samples.

### Tier 1 — Execute Now (done)
1. **GSE165897** — 11 paired, public, full SAE (scRNA-seq); processed.  
2. **TCGA-OV** — ~300 baseline, public; baseline validation complete.

### Tier 2 — Submit Access
3. **BriTROC-1** — 276 paired, EGA; submit application (2–4 weeks).

### Tier 3 — Contact Authors
4. **Williams et al.** — 18 patients, 5–20 samples/patient; data embargoed.

### Tier 4 — Optional / Not Suitable
- GSE241908 (n=7, bevacizumab); GSE217177 (cross-sectional); GSE184880 (treatment-naive only).

---

## 7. Execution Status & Final Strategy

### Completed
- Paired sample identification (e.g. MSK-SPECTRUM 40 paired; script `download_cbioportal_paired.py`; output `paired_patients.json`).  
- Molecular download scripts created; cBioPortal API issues led to pivot.  
- SAE computation script created (placeholder until data).  
- **Decision: Focus on TCGA-OV** — manifest created (434 RNA-seq, 482 mutation, 1,204 clinical); download script ready.  
- **GSE165897:** Downloaded, processed, pathway kinetics computed.  
- **TCGA-OV:** Baseline SAE and outcome correlations completed (426 patients).

### MSK-SPECTRUM Abandoned
- Reasons: 40 paired only, no gene expression via API, limited SAE utility.  
- Fallback: TCGA-OV baseline + GSE165897 serial.

### TCGA-OV Execution Plan (done)
1. Download (GDC manifest + gdc-client).  
2. Process RNA-seq (normalize, match clinical).  
3. Compute baseline SAE (7D mechanism vectors).  
4. Correlate with outcomes (OS, PFS, recurrence).

### Immediate Next Steps (if resuming)
- Option A: Manual cBioPortal download for MSK-SPECTRUM if needed (see `MANUAL_DOWNLOAD_INSTRUCTIONS.md`).  
- Option B: GSE165897 (done).  
- Option C: BriTROC-1 EGA application; Williams et al. contact.

---

## 8. GSE165897 Analysis

### Overview
- **Access:** Public (GEO).  
- **Cohort:** 11 HGSOC, paired pre/post-NACT scRNA-seq (22 samples).  
- **Treatment:** Platinum-based NACT.  
- **Publication:** Zhang et al., Science Advances 8, eabm1831 (2022).

### Processing Status — Complete
- Data: `scripts/data_acquisition/sae/GSE165897_*.tsv.gz` and series matrices.  
- Script: `scripts/serial_sae/pathway_kinetics_gse165897.py`.  
- Results: `data/serial_sae/gse165897/results/` (pathway_scores.csv, pathway_kinetics.csv, kinetic_patterns.csv, correlations, heatmap, report).

### Pathway Kinetics (n=11)

| Pathway | Mean Δ | Median Δ | Increases (Δ>0.05) | Decreases (Δ<-0.05) |
|---------|--------|----------|--------------------|---------------------|
| DDR     | -0.0140| -0.0129  | 0/11               | 0/11                |
| MAPK    | 0.0002 | -0.0007  | 0/11               | 1/11                |
| PI3K    | -0.0097| -0.0139  | 1/11               | 0/11                |
| VEGF    | 0.0268 | 0.0125   | 4/11               | 1/11                |

### Resistance-Stratified (Partial — 3/11 with PFI)
- Resistant: PFI &lt; 180 days; Sensitive: PFI ≥ 180 days.  
- **VEGF:** Resistant Δ = +0.0694 (n=1), Sensitive Δ = -0.0323 (n=2); difference 0.1017 (Cohen’s d 6.99).  
- Interpretation: VEGF activation post-chemotherapy in resistant vs suppression in sensitive — angiogenesis as potential resistance mechanism.  
- **Limitation:** n=3; need full PFI for remaining 8 (Table S1, Zhang et al.).

### Missing Clinical Data
- PFI for 8/11; RECIST, platinum sensitivity, PFS, OS.  
- Action: Extract PFI from Table S1 (Science Advances supp) or GenomeSpy (Lahtinen et al.); create `resistance_labels.json`; re-run stratified analysis.

---

## 9. TCGA-OV Validation

- **Samples:** 434 processed; 426 with complete SAE + clinical.  
- **Pathway–OS correlations (significant):**  
  - HER2: r=-0.182, p=0.0002  
  - PI3K: r=-0.143, p=0.0032  
  - Efflux: r=-0.148, p=0.0022  
  - RAS_MAPK: r=-0.111, p=0.0218  
- **Validation:** Aligned with TCGA 2011 and Verhaak 2013 (TP53, RB1, PI3K, FOXM1, NOTCH, HR deficiency, mesenchymal subtype).

---

## 10. Code Audit & Gaps

### Present
- Pathway gene lists (DDR, MAPK, PI3K) in `biomarker_enriched_cohorts/scripts/compute_pathway_burden_features.py`.  
- Mutation-based pathway scoring.  
- CA-125/KELIM in `api/services/ca125_intelligence.py` (KELIM to be verified).

### Gaps
1. **VEGF pathway genes** — Add VEGF_GENES (e.g. VEGFA, VEGFR1, VEGFR2, HIF1A).  
2. **Expression-based pathway scoring** — Implement scoring from expression matrix (e.g. mean log2(expr+1)) in addition to mutation-based.  
3. **scRNA-seq pipeline** — Pseudo-bulk aggregation and pathway scores (scanpy/anndata) where not yet done.  
4. **KELIM** — Confirm implementation in CA125Intelligence service.

---

## 11. Next Steps & Action Items

### Immediate
1. **Extract remaining PFI** (GSE165897) from Table S1 or GenomeSpy; build full `resistance_labels.json`; re-run resistance-stratified analysis.  
2. **Code fixes:** Add VEGF genes; expression-based pathway scoring; verify KELIM.

### Short-term
3. **BriTROC-1:** Submit EGA access request.  
4. **Williams et al.:** Contact corresponding author for data timeline.

### Long-term
5. BriTROC-1 analysis after approval (276 paired).  
6. Combined meta-analysis (GSE165897 + BriTROC-1) when available.

---

## 12. Quick Reference

- **Need:** Paired samples (baseline + progression), genomic data (WES or RNA-seq), outcomes (TTP/PFS/response), n≥20, accessible.  
- **Search:** PubMed (ovarian serial/ctDNA/organoid); GEO (ovarian serial/longitudinal); cBioPortal (ovarian, recurrent); dbGaP (ovarian).  
- **Per dataset report:** Citation, cohort, data availability, SAE compatible?, outcomes, feasibility.  
- **Deliverable:** Single doc with sections per dataset, summary table, priority tiers (this master replaces previous separate PHASE2_*.md files).  
- **Success:** ≥1–2 datasets with n≥20 paired, accessible, SAE computable; ideal 3–5 datasets, public.

---

## 13. Files & Locations

### Data
- **GSE165897:** `scripts/data_acquisition/sae/GSE165897_*.tsv.gz`, series matrices.  
- **TCGA-OV:** `data/serial_sae/tcga_ov/` (manifest, RNA-seq, mutation, clinical).  
- **cBioPortal paired:** `data/serial_sae/cbioportal_paired/paired_patients.json` (MSK-SPECTRUM, abandoned path).

### Results
- **GSE165897:** `data/serial_sae/gse165897/results/` (pathway_kinetics.csv, pathway_scores.csv, report, heatmap).  
- **TCGA-OV:** `data/serial_sae/tcga_ov/results/` (pathway_correlations, baseline SAE).

### Scripts
- `scripts/serial_sae/pathway_kinetics_gse165897.py`  
- `scripts/serial_sae/download_tcga_ov_gdc.py`, `download_tcga_ov_files.sh`  
- `scripts/serial_sae/download_cbioportal_paired.py`, `download_molecular_data.py`, `compute_serial_sae.py`  
- `scripts/serial_sae/README.md`, `MANUAL_DOWNLOAD_INSTRUCTIONS.md` (if present)

---

## 14. References

- **GSE165897:** Zhang et al. (2022). Science Advances 8(8): eabm1831. DECIDER cohort: Lahtinen et al., Cancer Cell (2023).  
- **TCGA-OV:** TCGA Research Network (2011). Nature 474(7353): 609–615. Verhaak et al. (2013). J Clin Invest 123(1): 517–525.  
- **BriTROC-1:** EGA EGAS00001007292.  
- **MSK-SPECTRUM:** dbGaP phs002857.v3.p1; cBioPortal msk_spectrum_tme_2022.

---

**Doctrine status:** Active.  
**Applies to:** All Phase 2 serial SAE validation work.  
**This document is the single source of truth; previous PHASE2_*.md files have been retired.**
