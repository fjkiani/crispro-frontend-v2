# Phase 2 Final Strategy: TCGA-OV Focus

**Date:** January 13, 2026  
**Decision:** Focus on TCGA-OV, leave MSK-SPECTRUM

---

## 🎯 DECISION RATIONALE

### MSK-SPECTRUM: Abandoned ❌
**Reasons:**
- **Cohort too small:** 40 paired patients (insufficient for robust analysis)
- **No gene expression data:** Cannot compute expression-based SAE pathway scores
- **Limited utility:** Without expression data, cannot validate pathway kinetics hypothesis

### TCGA-OV: Primary Focus ✅
**Advantages:**
- **Large cohort:** ~300-400 samples with full RNA-seq
- **Comprehensive data:** Expression, mutations, clinical outcomes
- **Public access:** No approval needed, direct download
- **Validated:** Well-curated TCGA data
- **Baseline validation:** Perfect for establishing SAE pathway score distributions

---

## 📊 DATA ASSETS

### TCGA-OV (Primary)
- **434 RNA-seq files** - Full transcriptome expression
- **482 mutation files** - Somatic mutations
- **1,204 clinical files** - Survival, recurrence, treatment
- **~300-400 samples** - Large cohort for validation
- **Status:** ✅ Manifest created, ready to download

### GSE165897 (Secondary - Optional)
- **11 paired patients** - Treatment-naïve vs post-NACT
- **Single-cell RNA-seq** - High-resolution pathway analysis
- **Status:** ⏳ Script ready, can proceed in parallel

---

## 🚀 EXECUTION PLAN

### Phase 1: TCGA-OV Download (1-2 hours)
```bash
# Download files
./scripts/serial_sae/download_tcga_ov_files.sh

# Or manually:
gdc-client download -m data/serial_sae/tcga_ov/tcga_ov_rnaseq_manifest.csv
```

**Output:**
- RNA-seq expression files (~300-400 samples)
- Mutation files (MAF format)
- Clinical data files

### Phase 2: Process RNA-seq Data (2-4 hours)
**Tasks:**
1. Extract expression matrices from GDC files
2. Normalize expression data (TPM, FPKM, or log-normalized)
3. Match samples to clinical data
4. Prepare for SAE computation

**Script:** `scripts/serial_sae/process_tcga_ov_rnaseq.py` (to be created)

### Phase 3: Compute Baseline SAE Scores (2-3 hours)
**Tasks:**
1. Compute pathway scores for all TCGA-OV samples
2. Generate 7D mechanism vectors: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
3. Match with clinical outcomes (survival, recurrence)
4. Validate pathway score distributions

**Script:** `scripts/serial_sae/compute_tcga_ov_sae.py` (to be created)

### Phase 4: Baseline Validation Analysis (2-3 hours)
**Tasks:**
1. Correlate SAE pathway scores with clinical outcomes
2. Stratify by survival, recurrence, treatment response
3. Identify prognostic pathway signatures
4. Generate baseline validation report

**Output:**
- Pathway score distributions
- Correlation with outcomes
- Prognostic pathway signatures
- Validation report

---

## 📋 DELIVERABLES

### Immediate (This Week)
1. ✅ TCGA-OV manifest created
2. ⏳ TCGA-OV files downloaded
3. ⏳ RNA-seq data processed
4. ⏳ Baseline SAE scores computed
5. ⏳ Validation analysis completed

### Documentation
- `docs/TCGA_OV_DOWNLOAD_STRATEGY.md` ✅
- `docs/PHASE2_FINAL_STRATEGY.md` ✅ (this file)
- `scripts/serial_sae/download_tcga_ov_gdc.py` ✅
- `scripts/serial_sae/download_tcga_ov_files.sh` ✅

---

## 🎯 SUCCESS METRICS

### Minimum Viable
- ✅ Download TCGA-OV RNA-seq data (~300 samples)
- ✅ Process expression data
- ✅ Compute baseline SAE pathway scores
- ✅ Correlate with clinical outcomes

### Ideal
- ✅ Validate pathway score distributions
- ✅ Identify prognostic pathway signatures
- ✅ Establish baseline for serial monitoring comparison
- ✅ Generate publication-ready figures

---

## ⚠️ LIMITATIONS & MITIGATION

### Limitation: No Paired Samples in TCGA-OV
**Impact:** Cannot compute serial SAE (ΔSAE) for resistance prediction

**Mitigation:**
- Use TCGA-OV for baseline validation
- Establish pathway score distributions
- Use GSE165897 (11 paired) for serial monitoring proof-of-concept
- Future: BriTROC-1 (276 paired) for large-scale serial validation

### Limitation: Primary Tumors Only
**Impact:** Cannot track pathway changes during progression

**Mitigation:**
- Focus on baseline pathway signatures
- Validate against clinical outcomes
- Use as foundation for serial monitoring studies

---

## 📈 NEXT PHASES

### Phase 3: GSE165897 Serial Analysis (Optional)
- 11 paired patients (treatment-naïve → post-NACT)
- Single-cell resolution
- Compute ΔSAE for resistance prediction

### Phase 4: BriTROC-1 Validation (Future)
- 276 paired patients (largest cohort)
- Controlled access via EGA
- Large-scale serial monitoring validation

---

**Status:** ✅ **STRATEGY FINALIZED**  
**Next:** Execute TCGA-OV download and processing
