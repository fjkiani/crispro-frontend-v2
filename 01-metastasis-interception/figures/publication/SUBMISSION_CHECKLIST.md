# 📋 PUBLICATION SUBMISSION CHECKLIST

**Title:** Stage-Aware CRISPR Design for Metastatic Cascade Interception: A Multi-Modal Validation Framework  
**Submission Date:** October 2025  
**Target Journal:** *Nature Methods* / *Nature Biotechnology* / *Cell Systems*

---

## ✅ PRE-SUBMISSION CHECKLIST

### Manuscript Components
- [ ] **Title Page**
  - [ ] Full title
  - [ ] Short title (≤50 characters)
  - [ ] Author list with affiliations
  - [ ] Corresponding author contact
  - [ ] Word count (Abstract: 185, Main text: ~4000)
  - [ ] Conflict of interest statement
  - [ ] Author contributions (CRediT taxonomy)

- [ ] **Abstract** ✅
  - [X] Background (1-2 sentences)
  - [X] Methods (2-3 sentences)
  - [X] Results (3-4 sentences with key metrics)
  - [X] Conclusions (1-2 sentences)
  - [X] ≤250 words
  - [X] RUO disclaimers included

- [ ] **Main Text**
  - [ ] Introduction (motivation, prior work, gap)
  - [X] Methods (complete technical details)
  - [ ] Results (8 validation analyses with figures)
  - [ ] Discussion (limitations, future work, impact)
  - [ ] References (numbered, formatted per journal style)

- [X] **Figures** ✅ **COMPLETE**
  - [X] Figure 1: Metastatic cascade workflow
  - [X] Figure 2A-D: Per-step validation (AUROC/AUPRC)
  - [X] Figure 3: Specificity matrix with enrichment
  - [X] Figure 4: Precision@K curves
  - [X] Figure 5: Ablation study (signal importance)
  - [X] Figure 6: Structural validation (pLDDT, iPTM, per-step, best structure) **NEW**
  - [X] Figure 7: Guide RNA validation distributions
  - [X] Figure 8: Target Lock heatmap
  - [X] All legends include n, RUO notes, statistical tests

- [X] **Tables** ✅ **COMPLETE**
  - [X] Table S2: Comprehensive per-step metrics (CSV + LaTeX)
  - [X] Table S4: Structural validation (15 guides with pLDDT, iPTM, verdict) **NEW**
  - [ ] Table 1 (optional): Summary of key performance metrics

- [ ] **Supplementary Materials**
  - [X] Supplementary Figures S1-S3 (calibration, effect sizes, confounders)
  - [X] Supplementary Methods (extended technical details)
  - [X] Table S2 (comprehensive validation metrics)
  - [ ] Supplementary Note 1: Statistical methodology
  - [ ] Supplementary Note 2: Infrastructure deployment guide

---

## 📊 DATA & CODE AVAILABILITY

### Data Repository (Zenodo/Figshare)
- [X] Raw validation datasets (9 CSV files)
- [X] Ground truth gene sets with citations (JSON)
- [X] Processed metrics (CSV + JSON)
- [X] Figure source data
- [ ] DOI assignment
- [ ] Data availability statement in manuscript

### Code Repository (GitHub)
- [X] Complete analysis pipeline (scripts/metastasis/)
- [X] Reproducibility guide (one-command execution)
- [X] Frozen environment (environment.yml)
- [X] Docker/Modal deployment instructions
- [ ] README with quick start
- [ ] LICENSE file (MIT/Apache 2.0)
- [ ] DOI via Zenodo
- [ ] Code availability statement in manuscript

### Model Weights & Services
- [X] Evo2 proxy endpoints documented
- [X] AlphaFold 3 Server structural validation **NEW**
- [X] Boltz-2 service deployment guide (backup infrastructure)
- [ ] Enformer deployment guide (when deployed)
- [ ] Model version tracking (provenance metadata)

---

## 🔬 RESEARCH INTEGRITY

### Reproducibility
- [X] Random seed fixed (seed=42)
- [X] Bootstrap iterations documented (n=1000)
- [X] All dependencies versioned (requirements.txt)
- [X] Container digests pinned (PUBLICATION_SPEC_LOCKED.md)
- [X] Provenance metadata in all outputs
- [X] One-command reproduction script tested

### Statistical Rigor
- [X] Bootstrap confidence intervals (1000 iterations, stratified)
- [X] Multiple testing corrections (where applicable)
- [X] Effect sizes reported (Cohen's d)
- [X] Power analysis (n=38 genes adequate for detection)
- [X] Outlier analysis (confounder check)
- [X] Calibration assessment (reliability diagrams)

### Research Use Only (RUO) Disclaimers
- [X] Chromatin stubs documented in Abstract
- [X] Chromatin stubs documented in Methods
- [X] Structural validation (AF3) documented in Methods **UPDATED**
- [X] Structural validation results documented (RESULTS_STRUCTURAL.md) **NEW**
- [X] All figure legends include RUO notes
- [X] Limitations section addresses both

---

## 📧 SUBMISSION MATERIALS

### Cover Letter
- [ ] Journal editor addressed by name
- [ ] Brief summary of significance (2-3 sentences)
- [ ] Key findings highlighted (AUROC 0.976, n=38 genes)
- [ ] Novelty statement (first stage-aware CRISPR framework)
- [ ] Suggested reviewers (3-5 names with expertise)
- [ ] Competing interests disclosure
- [ ] Data/code availability confirmation

### Author Information
- [ ] ORCID IDs for all authors
- [ ] Current affiliations verified
- [ ] Corresponding author contact (email, phone)
- [ ] Author contributions (CRediT): conceptualization, methodology, software, validation, formal analysis, investigation, data curation, writing, visualization, supervision, funding

### Supplementary Files
- [ ] Supplementary Information PDF (figures, methods, notes)
- [ ] Supplementary Data 1: All validation datasets (ZIP)
- [ ] Supplementary Data 2: Code repository (GitHub link + Zenodo DOI)
- [ ] Supplementary Data 3: Figure source data
- [ ] Reporting summary (Nature journals)

---

## 🎯 JOURNAL-SPECIFIC REQUIREMENTS

### Nature Methods
- [ ] Reporting summary (automated template)
- [ ] Life Sciences Reporting Summary
- [ ] Statistics checklist
- [ ] Software availability statement
- [ ] Ethics statement (if applicable)
- [ ] Research involving Human Participants statement

### Nature Biotechnology
- [ ] Similar to Nature Methods
- [ ] Industrial/commercial application statement
- [ ] Patent disclosure

### Cell Systems
- [ ] STAR Methods section (detailed technical protocols)
- [ ] Key Resources Table
- [ ] Lead Contact information
- [ ] Materials Availability statement
- [ ] Data and Code Availability statement

---

## ✅ FINAL CHECKS (BEFORE SUBMISSION)

### Manuscript Quality
- [ ] Spell check (US English)
- [ ] Grammar check (Grammarly/language tool)
- [ ] Reference formatting (journal style)
- [ ] Figure quality (≥300 DPI for print)
- [ ] Figure file formats (PDF/EPS for vector, TIFF/PNG for raster)
- [ ] Table formatting (editable, not images)
- [ ] Page numbers on all pages
- [ ] Line numbers for review

### Technical Verification
- [ ] All URLs in manuscript are live
- [ ] All DOIs resolve correctly
- [ ] GitHub repository is public
- [ ] Data repository is public
- [ ] Reproduction script runs successfully
- [ ] All figures regenerate from source data
- [ ] No proprietary/confidential data included

### Submission Portal
- [ ] Journal account created
- [ ] Manuscript uploaded (Word/LaTeX)
- [ ] Figures uploaded (separate files)
- [ ] Tables uploaded (editable format)
- [ ] Supplementary materials uploaded
- [ ] Cover letter uploaded
- [ ] Author information complete
- [ ] Suggested reviewers entered
- [ ] Conflicts of interest declared
- [ ] Copyright/license agreement signed

---

## 📌 POST-SUBMISSION TRACKING

### Timeline
- [ ] Submission confirmation received
- [ ] Editor decision (expect 2-4 weeks)
- [ ] Reviews received (expect 4-8 weeks)
- [ ] Revisions completed (if requested)
- [ ] Resubmission (if applicable)
- [ ] Acceptance notification
- [ ] Proofs reviewed
- [ ] Publication date

### Communication
- [ ] Preprint posted (bioRxiv/medRxiv)
- [ ] Lab website updated
- [ ] Social media announcement prepared
- [ ] Press release prepared (if high impact)

---

## 🎓 ESTIMATED COMPLETION TIME

- **Manuscript finalization:** 2-3 days
- **Data/code repository setup:** 1 day
- **Supplementary materials:** 1-2 days
- **Final checks & submission:** 1 day

**Total:** ~5-7 days to submission

---

**Status:** Ready for final manuscript assembly  
**Next Action:** Complete Introduction, Results, Discussion sections  
**Blocker:** None - all validation complete, methods documented

