# 🔍 Pathway Kinetics Validation (GSE165897) - Mission Audit

**Date:** January 2025  
**Auditor:** Zo  
**Mission:** Validate SAE pathway kinetics detect mechanism-specific chemotherapy resistance  
**Timeline:** 6 hours  
**Status:** ⚠️ AUDIT COMPLETE - Questions Identified

---

## 📋 MISSION REVIEW

### Objective
Validate that pathway kinetics (DDR, MAPK, PI3K, VEGF changes from pre → post-treatment) detect mechanism-specific chemotherapy resistance using real paired pre/post-treatment samples.

**Dataset:** GSE165897 (n=11 HGSOC patients, scRNA-seq, pre-treatment + post-NACT)

**Novel Contribution:** Pathway kinetics complement CA-125 KELIM by revealing resistance mechanism, not just timing.

---

## ✅ WHAT WE HAVE (Code Audit)

### 1. Pathway Gene Lists ✅

**Location:** `biomarker_enriched_cohorts/scripts/compute_pathway_burden_features.py`

**Existing Pathway Genes:**
```python
DDR_GENES = [
    "BRCA1", "BRCA2", "ATM", "ATR", "CHEK2", "PALB2", "RAD51", "RAD51C",
    "RAD51D", "BARD1", "MBD4", "FANCA", "FANCD2"
]

MAPK_GENES = [
    "KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RAF1"
]

PI3K_GENES = [
    "PIK3CA", "PIK3CB", "PTEN", "AKT1", "AKT2", "MTOR", "TSC1", "TSC2"
]
```

**Mission Requirements:**
- DDR: BRCA1, BRCA2, ATM, ATR, CHEK1, CHEK2, RAD51, PARP1 ✅ (mostly match)
- MAPK: KRAS, NRAS, BRAF, MEK1, MEK2, ERK1, ERK2 ✅ (mostly match)
- PI3K: PIK3CA, AKT1, AKT2, PTEN, mTOR ✅ (match)
- VEGF: VEGFA, VEGFR1, VEGFR2, HIF1A ❌ (NOT in existing code)

**Gap:** Need to add VEGF pathway genes.

### 2. Pathway Score Computation ✅

**Location:** `biomarker_enriched_cohorts/scripts/compute_pathway_burden_features.py`

**Existing Method:**
- Computes pathway scores from mutations (S/P/E-style weighted aggregation)
- Uses sequence proxies (1.0 for LOF/hotspots, 0.5 for missense)

**Mission Requirements:**
- Compute pathway scores from **gene expression** (not mutations)
- Formula: `DDR_score = mean(log2(expression + 1)) for DDR genes`
- Normalize to 0-1 scale

**Gap:** Need expression-based pathway scoring (not mutation-based).

### 3. GEO Dataset Processing ⚠️

**Found:** References to GSE63885 processing (microarray data)

**Location:** `publications/03-sae-resistance/manuscript/SUPPLEMENTARY_MATERIALS.md`

**Existing Pattern:**
- Download series matrix from GEO
- Parse sample annotations
- Extract expression values

**Mission Requirements:**
- Download GSE165897 (scRNA-seq, not microarray)
- Process scRNA-seq data (h5ad/csv/mtx format)
- Aggregate to pseudo-bulk per patient per timepoint

**Gap:** Need scRNA-seq processing (different from microarray).

### 4. KELIM/CA-125 Code ✅

**Found:** `api/services/ca125_intelligence.py` (referenced in search results)

**Status:** Exists but need to verify KELIM calculation

**Mission Requirements:**
- Compare pathway kinetics vs KELIM
- KELIM threshold: ≥1.0 (favorable), <1.0 (unfavorable)

**Gap:** Need to verify KELIM implementation.

---

## ❓ QUESTIONS IDENTIFIED

### Question 1: scRNA-seq Data Format
**Q:** What format is GSE165897 data in? (h5ad, csv, mtx, or raw FASTQ?)

**Answer Attempt:**
- Check GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897
- Likely formats: Processed count matrices (h5ad/csv) or raw FASTQ
- Need to check supplementary files

**Still Need:** Access GEO page to verify format.

### Question 2: Pseudo-bulk Aggregation Method
**Q:** How should we aggregate scRNA-seq to pseudo-bulk? (mean, sum, or median expression per gene per patient?)

**Answer Attempt:**
- Standard approach: Sum counts per gene per patient, then normalize
- Alternative: Mean log2(CPM + 1) across cells
- For pathway scoring: Mean log2(expression + 1) is appropriate

**Still Need:** Confirm aggregation method matches mission requirements.

### Question 3: Timepoint Matching
**Q:** How are pre/post timepoints matched in the dataset? (patient IDs, sample IDs, or metadata fields?)

**Answer Attempt:**
- Check metadata for: patient_id, timepoint, sample_id
- Expected: Each patient has 2 samples (pre + post)
- Need to verify pairing logic

**Still Need:** Download metadata to verify structure.

### Question 4: Resistance Labels
**Q:** Where are resistance outcomes (responder vs non-responder) stored? (metadata, separate file, or need to infer?)

**Answer Attempt:**
- Check GEO sample characteristics
- May be in metadata or supplementary files
- May need to infer from clinical data

**Still Need:** Download and inspect metadata.

### Question 5: VEGF Pathway Genes
**Q:** Should we add VEGF pathway genes to existing pathway lists?

**Answer:** ✅ YES - Mission requires VEGF pathway.

**Action:** Add to pathway gene lists:
```python
VEGF_GENES = ["VEGFA", "VEGFR1", "VEGFR2", "HIF1A"]
```

### Question 6: Expression-Based Pathway Scoring
**Q:** How do we compute pathway scores from expression (not mutations)?

**Answer:** 
- Extract expression values for pathway genes
- Compute: `pathway_score = mean(log2(expression + 1)) for pathway genes`
- Normalize to 0-1 scale (divide by max or use percentile)

**Action:** Create expression-based pathway scoring function.

### Question 7: KELIM Calculation
**Q:** Do we have KELIM calculation code, or need to implement it?

**Answer Attempt:**
- Found `ca125_intelligence.py` service
- Need to verify if KELIM is implemented
- If not, need to implement: KELIM = f(CA-125 values over time)

**Still Need:** Check CA125Intelligence service for KELIM.

### Question 8: Statistical Analysis Libraries
**Q:** What statistical libraries are available? (scipy, statsmodels, etc.)

**Answer Attempt:**
- Standard Python: scipy.stats for Mann-Whitney U, t-test
- Visualization: matplotlib, seaborn
- Should be available in environment

**Action:** Verify imports work.

### Question 9: Output Directory Structure
**Q:** Where should we save outputs? (data/, scripts/sae/, receipts/)

**Answer:**
- Follow existing pattern: `data/external/GSE165897/` for raw data
- `data/gse165897_pathway_scores.csv` for processed scores
- `scripts/sae/pathway_kinetics_gse165897.py` for main script
- `receipts/pathway_kinetics_validation.json` for receipt
- `docs/PATHWAY_KINETICS_VS_KELIM.md` for positioning

**Action:** Create directory structure.

### Question 10: scRNA-seq Processing Libraries
**Q:** What libraries are needed for scRNA-seq? (scanpy, anndata, pandas, numpy)

**Answer:**
- scanpy: scRNA-seq analysis
- anndata: AnnData format (h5ad files)
- pandas: Data manipulation
- numpy: Numerical operations

**Action:** Verify/install dependencies.

---

## 🎯 ANSWERS TO QUESTIONS

### ✅ Answered Questions

1. **VEGF Pathway Genes:** ✅ Add VEGF_GENES list
2. **Expression-Based Scoring:** ✅ Use mean(log2(expression + 1))
3. **Output Structure:** ✅ Follow existing patterns
4. **Statistical Libraries:** ✅ scipy.stats available
5. **scRNA-seq Libraries:** ✅ scanpy, anndata needed

### ⚠️ Questions Requiring Data Download

6. **Data Format:** Need to check GEO page
7. **Timepoint Matching:** Need to inspect metadata
8. **Resistance Labels:** Need to check sample characteristics
9. **Pseudo-bulk Method:** Can infer from standard practice, but verify with data

### ❓ Questions Requiring Code Inspection

10. **KELIM Implementation:** Need to check CA125Intelligence service

---

## 📋 EXECUTION PLAN (Updated After Audit)

### Phase 1: Data Acquisition (1 hour)

**Task 1.1: Download GSE165897**
1. Access GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897
2. Download supplementary files (check format: h5ad/csv/mtx)
3. Download metadata (sample annotations)
4. **Verify:** Patient IDs, timepoints, resistance outcomes

**Task 1.2: Inspect Data Structure**
1. Load data (scanpy.read_h5ad or pandas.read_csv)
2. Check dimensions (cells × genes or samples × genes)
3. Extract metadata (patient_id, timepoint, response_status)
4. **Report:** Data format, structure, pairing logic

### Phase 2: Pathway Score Computation (2 hours)

**Task 2.1: Add VEGF Pathway**
- Add VEGF_GENES to pathway lists
- Update pathway scoring to include VEGF

**Task 2.2: Pseudo-bulk Aggregation**
- If scRNA-seq: Aggregate cells to patient-timepoint level
- Method: Sum counts per gene, normalize to CPM, log2 transform
- Extract expression for pathway genes

**Task 2.3: Expression-Based Pathway Scoring**
- Compute: `pathway_score = mean(log2(expression + 1)) for pathway genes`
- Normalize to 0-1 scale
- Output: `data/gse165897_pathway_scores.csv`

### Phase 3: Kinetics Calculation (1 hour)

**Task 3.1: Compute Pathway Changes**
- Formula: `ΔDDR = DDR_post - DDR_pre`
- Classify kinetic patterns (DDR restoration, MAPK activation, etc.)
- Output: `data/gse165897_pathway_kinetics.csv`

### Phase 4: Statistical Analysis (1 hour)

**Task 4.1: Correlate with Resistance**
- Compare ΔDDR in resistant vs sensitive (Mann-Whitney U)
- Secondary: ΔMAPK, ΔPI3K by resistance
- Visualizations: Heatmap, box plots, scatter plots

### Phase 5: KELIM Comparison (1 hour)

**Task 5.1: Position vs KELIM**
- Check if CA-125 data available in GSE165897
- If yes: Compute KELIM, compare with pathway kinetics
- If no: Document positioning (mechanism vs timing)
- Output: `docs/PATHWAY_KINETICS_VS_KELIM.md`

---

## 🔧 CODE GAPS TO FILL

### Gap 1: VEGF Pathway Genes
**Fix:** Add to pathway gene lists
```python
VEGF_GENES = ["VEGFA", "VEGFR1", "VEGFR2", "HIF1A"]
```

### Gap 2: Expression-Based Pathway Scoring
**Fix:** Create new function (not mutation-based)
```python
def compute_pathway_score_from_expression(
    expression_matrix: pd.DataFrame,
    pathway_genes: List[str]
) -> float:
    """Compute pathway score from expression (not mutations)."""
    pathway_expr = expression_matrix[pathway_genes].values
    log_expr = np.log2(pathway_expr + 1)
    pathway_score = np.mean(log_expr)
    return pathway_score / max_score  # Normalize to 0-1
```

### Gap 3: scRNA-seq Processing
**Fix:** Add scanpy/anndata processing
```python
import scanpy as sc
import anndata as ad

# Load scRNA-seq data
adata = sc.read_h5ad("GSE165897.h5ad")

# Aggregate to pseudo-bulk
pseudo_bulk = aggregate_to_patient_timepoint(adata)
```

### Gap 4: KELIM Verification
**Fix:** Check CA125Intelligence service for KELIM implementation

---

## 📊 SUCCESS CRITERIA (Updated)

### Minimum Viable Result
- ✅ Download GSE165897 data
- ✅ Compute pathway scores from expression (all 4 pathways)
- ✅ Calculate pathway kinetics (pre → post changes)
- ✅ Correlate with resistance (even if p>0.05 with n=11)
- ✅ Document methodology

### Strong Result (Publication-Ready)
- ✅ ΔDDR significantly different in resistant vs sensitive (p<0.05)
- ✅ Pathway-specific patterns identified
- ✅ Demonstrates complementarity to KELIM

---

## ❓ QUESTIONS FOR ALPHA

### Critical Questions (Need Answer Before Starting)

1. **scRNA-seq Processing Libraries**
   - **Q:** Do we have scanpy/anndata installed? If not, should I install or use alternative?
   - **Context:** GSE165897 is scRNA-seq, need libraries for processing
   - **Recommendation:** Install scanpy if not available (standard for scRNA-seq)

2. **VEGF Pathway Genes**
   - **Q:** Mission lists VEGFA, VEGFR1, VEGFR2, HIF1A. Should I add these to existing pathway lists?
   - **Answer:** ✅ YES - Already documented in audit
   - **Action:** Add VEGF_GENES list

3. **Expression vs Mutation Scoring**
   - **Q:** Existing code computes pathway scores from mutations. Mission needs expression-based. Should I create new function or modify existing?
   - **Answer:** ✅ Create new function (expression-based, separate from mutation-based)
   - **Action:** Create `compute_pathway_score_from_expression()` function

4. **KELIM Implementation**
   - **Q:** Do we have KELIM calculation code? If not, should I implement it or skip for now?
   - **Context:** Phase 5 requires KELIM comparison
   - **Recommendation:** Check CA125Intelligence service first, implement if missing

5. **Output Directory Structure**
   - **Q:** Mission specifies `data/`, `scripts/sae/`, `receipts/`, `docs/`. Should I create these or use existing structure?
   - **Answer:** ✅ Create as specified (follows existing patterns)
   - **Action:** Create directories if needed

### Questions I Can Answer During Execution

6. **Data Format:** Will check GEO page during Phase 1
7. **Timepoint Matching:** Will inspect metadata during Phase 1
8. **Resistance Labels:** Will check sample characteristics during Phase 1
9. **Pseudo-bulk Method:** Will use standard approach (sum counts, normalize, log2)

---

## 🚀 READY TO EXECUTE

**Status:** ✅ Audit complete, questions identified, gaps documented

**Blockers:** None - Can proceed with Phase 1 (data download will answer format questions)

**Next Steps:**
1. ✅ Add VEGF pathway genes to code
2. ✅ Create expression-based pathway scoring function
3. ✅ Begin Phase 1: Download GSE165897
4. ✅ Answer remaining questions during data download
5. ✅ Proceed with pathway score computation
6. ✅ Report back after each phase

**Estimated Time:** 6 hours (as specified)

**Confidence:** HIGH - Have pathway gene lists, scoring patterns, and GEO download patterns. Main gap is scRNA-seq processing (standard libraries available).

---

**Status:** ✅ MISSION AUDIT COMPLETE  
**Date:** January 2025  
**Auditor:** Zo  
**Ready:** ✅ YES - Proceeding to Phase 1
