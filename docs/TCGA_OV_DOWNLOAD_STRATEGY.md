# TCGA-OV Direct Download Strategy

**Date:** January 13, 2026  
**Source:** TCGA GDC Data Portal  
**Cohort:** TCGA-OV (High-Grade Serous Ovarian Cancer)

---

## 📊 Data Available

### RNA-seq
- **~300 primary tumors** with full RNA-seq coverage
- **Format:** Gene Expression Quantification (HTSeq, FPKM, TPM)
- **Coverage:** Comprehensive transcriptome

### Mutations
- **~300 samples** with WGS/WXS
- **Format:** Masked Somatic Mutations (MAF files)
- **Coverage:** Full exome/genome

### Clinical Data
- **Survival:** Overall survival (OS), progression-free survival (PFS)
- **Recurrence:** Time to recurrence, recurrence status
- **Treatment:** Chemotherapy regimens, response

### Paired Samples
- **Limited:** Mostly primary tumors only
- **Note:** TCGA typically has primary-only samples
- **Use Case:** Baseline SAE computation, not serial monitoring

---

## 🎯 Use Cases

### Option A: Baseline SAE Validation
- Use TCGA-OV primary tumors for baseline SAE computation
- Validate pathway scores against clinical outcomes
- Compare with MSK-SPECTRUM paired samples

### Option B: Large Cohort Analysis
- ~300 samples for robust pathway score validation
- Stratify by clinical outcomes (survival, recurrence)
- Identify pathway signatures associated with outcomes

### Option C: Combine with Other Sources
- TCGA primary tumors as baseline
- MSK-SPECTRUM recurrent samples for progression
- Cross-validate pathway kinetics

---

## 📥 Download Methods

### Method 1: GDC Data Transfer Tool (Recommended)
```bash
# Install
pip install gdc-client

# Download using manifest
gdc-client download -m data/serial_sae/tcga_ov/tcga_ov_rnaseq_manifest.csv
```

### Method 2: GDC Data Portal Web Interface
1. Navigate to: https://portal.gdc.cancer.gov/
2. Select project: **TCGA-OV**
3. Add files:
   - RNA-Seq (Gene Expression Quantification)
   - Mutations (Masked Somatic Mutation)
   - Clinical (Clinical Supplement)
4. Download manifest
5. Use `gdc-client` to download files

### Method 3: Direct API Download (RECOMMENDED - No gdc-client needed)
```bash
# Download files directly via GDC API
python3 scripts/serial_sae/download_tcga_ov_api.py

# Or limit to first N files for testing:
python3 scripts/serial_sae/download_tcga_ov_api.py 10
```

**Advantages:**
- No external tools required (just Python + httpx)
- Direct HTTP download
- Progress tracking
- MD5 verification
- Automatic retry on failure

---

## ⚡ Timeline

- **Query & Manifest:** 5-10 minutes
- **Download:** 1-2 hours (depending on connection)
- **Processing:** 2-4 hours (RNA-seq normalization, mutation processing)
- **Total:** 3-6 hours

---

## ✅ Pros

- **Comprehensive:** ~300 samples with full RNA-seq
- **Public Access:** No approval needed
- **Validated:** Well-curated TCGA data
- **Expression Data:** Full transcriptome coverage
- **Clinical Outcomes:** Survival, recurrence, treatment data

---

## ⚠️ Cons

- **Limited Paired Samples:** Mostly primary tumors only
- **Not Ideal for Serial Monitoring:** Need baseline + progression pairs
- **Large Download:** Several GB of data
- **Processing Time:** RNA-seq files need normalization

---

## 🔄 Integration with Serial SAE Mission

### Current Status
- **MSK-SPECTRUM:** 40 paired patients (primary + recurrent) ✅
- **GSE165897:** 11 paired patients (scRNA-seq) ✅
- **TCGA-OV:** ~300 primary tumors (baseline only) ⏳

### Strategy
1. **Use TCGA-OV for baseline validation:**
   - Compute SAE pathway scores for ~300 primary tumors
   - Validate against clinical outcomes (survival, recurrence)
   - Establish baseline pathway signatures

2. **Use MSK-SPECTRUM for serial monitoring:**
   - 40 paired patients (primary → recurrent)
   - Compute ΔSAE (pathway kinetics)
   - Validate resistance prediction

3. **Combine insights:**
   - Baseline signatures from TCGA-OV
   - Serial kinetics from MSK-SPECTRUM
   - Cross-validate findings

---

## 📋 Next Steps

1. **Download TCGA-OV Data** (1-2 hours)
   - Run: `python3 scripts/serial_sae/download_tcga_ov_gdc.py`
   - Download files using GDC client

2. **Process RNA-seq Data** (2-4 hours)
   - Normalize expression data
   - Compute pathway scores
   - Match with clinical outcomes

3. **Baseline SAE Analysis** (2-3 hours)
   - Compute SAE pathway scores for all samples
   - Correlate with survival/recurrence
   - Identify prognostic pathway signatures

4. **Integrate with Serial Data** (1-2 hours)
   - Compare baseline signatures (TCGA) with serial kinetics (MSK)
   - Validate pathway changes during progression

---

**Status:** ⏳ **READY TO EXECUTE**  
**Priority:** Medium (complements serial monitoring mission)
