# EXTRACTION PIECE 5.7: Additional Critical Discoveries

**Source**: Lines 6000-6500, 15000-15500, 18000-18500, 20000-20500, 25000-25500, 28000-28500 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 🚨 **CRITICAL DISCOVERIES**

### **1. Biomarker Correlation Service Implementation** (Lines 6000-6500)

**Complete Implementation**: Full biomarker correlation service built with:
- Statistical methods: Pearson, Spearman, Chi-square, Cohen's d
- Cross-validation stability testing
- Bootstrap confidence intervals
- Multiple testing correction (Bonferroni/FDR)
- Feature matrix building from SAE cohort data

**Key Code Structure**:
- `BiomarkerCorrelationService` class (688 lines)
- Methods: `load_sae_cohort_data()`, `build_feature_matrix()`, `compute_pearson_correlation()`, `compute_spearman_correlation()`, `compute_chi_square()`, `compute_cohen_d()`, `compute_cv_stability()`, `compute_bootstrap_ci()`, `apply_multiple_testing_correction()`
- Output: `sae_tcga_ov_platinum_biomarkers.json` with top features ranked by statistical significance

**Success Criteria**:
- Process all patients in SAE cohort
- Compute correlation for all 32K SAE features
- Rank by statistical significance + effect size
- Identify top 100 features with p < 0.01
- Bootstrap confidence intervals for top features

---

### **2. API Key Authentication Issues** (Lines 15000-15500)

**Problem**: Hundreds of `403 Forbidden` errors from SAE endpoint

**Root Cause**:
- Backend sending `X-API-Key` header
- Modal SAE service expecting different key or no key configured
- Mismatch between backend env var and Modal container env

**Fix Required**:
- Set `SAE_API_KEY` in both backend and Modal service
- Ensure keys match exactly
- Redeploy Modal service with correct key
- Restart backend with matching key

**Key Insight**: API key authentication must be configured consistently across services.

---

### **3. GRCh37 vs GRCh38 Assembly Mismatch** (Lines 18000-18500)

**Problem**: Ensembl 400 errors when requesting sequences:
```
{"error":"Cannot request a slice whose start (242497721) is greater than 242193529 for 2."}
```

**Root Cause**:
- TCGA-OV mutations from cBioPortal are on **GRCh37/hg19**
- Cohort script was hardcoding `"assembly": "GRCh38"`
- Positions like `2:242501817` are invalid for GRCh38 (beyond chromosome length)
- Ensembl returns 400 → SAE service wraps as 500

**Fix Applied**:
```python
# Changed in scripts/sae/extract_sae_features_cohort.py
payload = {
    "chrom": str(chrom).replace("chr", ""),
    "pos": int(pos),
    "ref": str(ref).upper(),
    "alt": str(alt).upper(),
    "model_id": model_id,
    "assembly": "GRCh37",  # <-- switched from GRCh38
    "window": 8192,
}
```

**Key Insight**: TCGA data uses GRCh37, not GRCh38. Always check data source assembly version.

---

### **4. Modal Architecture Cost Concerns** (Lines 18000-18500)

**User Concern**: "That's a horrible architecture - what if it lingers and keeps running and burns all of our credits?"

**Problem**:
- Modal web endpoint is **public, always-on** with no built-in:
  - Rate limiting
  - Per-job caps
  - Kill switch
  - Separation between RUO "cohort runs" and normal product traffic
- Cohort script blindly hammers endpoint until success/exhaustion
- No gating by "allow long-running RUO job" flag

**Proposed Hardening**:

**At Modal Layer**:
- Disable or gate public web endpoint for SAE
- Use Modal function calls from gated job runner instead of exposed HTTP endpoint
- OR keep web endpoint but add:
  - Simple auth (API key/header)
  - `MAX_REQUESTS_PER_MIN` + 429 responses

**At Backend & Script Layer**:
- Keep SAE cohort extraction **only callable from CLI scripts**, not user-facing APIs
- In `extract_sae_features_cohort.py`:
  - Enforce hard `MAX_PATIENTS` (e.g. 50-100 per run)
  - Require explicit `--allow-full-run` to go beyond
  - Add global cap on total variants per run (e.g. `max_variants = 10_000`)
  - Bail out early if error rate exceeds threshold (>30% 4xx/5xx)
- Add `ENABLE_SAE_COHORT_RUN=0` default, require explicit `=1` to enable

**Key Insight**: Modal architecture requires explicit cost controls and gating for RUO jobs.

---

### **5. WIWFM Integration Architecture** (Lines 25000-25500)

**Complete Integration Plan**: Detailed code examples for integrating SAE biomarkers into WIWFM:

**Step 1: Extract Patient SAE Features**
```python
def extract_patient_sae_features(mutations: List[Dict]) -> np.ndarray:
    # Aggregate SAE features across all mutations
    # Returns: Aggregated SAE feature vector (32K-dim)
```

**Step 2: Drug-Specific Biomarker Mapping**
```python
class SAEBiomarkerDrugMapper:
    def define_drug_mappings(self) -> Dict[str, Dict]:
        # Map drugs to relevant SAE biomarkers
        # Platinum agents: Direct correlation (trained on platinum response)
        # PARP inhibitors: Indirect (DNA repair proxy, same pathway)
```

**Step 3: Compute Drug-Specific SAE Score**
```python
def compute_patient_sae_score(
    drug_name: str,
    patient_sae_features: np.ndarray,
    biomarker_mapper: SAEBiomarkerDrugMapper
) -> float:
    # Returns: SAE score in [-1.0, +1.0] range (0 = neutral)
```

**Step 4: Apply SAE Boost to WIWFM Confidence**
```python
def apply_sae_biomarker_boost(
    drug_name: str,
    base_confidence: float,
    patient_mutations: List[Dict],
    biomarker_mapper: SAEBiomarkerDrugMapper
) -> Tuple[float, Dict]:
    # Returns: (boosted_confidence, provenance)
    # sae_score ∈ [-1, +1] → boost ∈ [-0.15, +0.15]
    # boosted_confidence = min(base_confidence + sae_boost, 0.95)
```

**Integration into Drug Scorer**:
- Extended `score_drug()` method with SAE parameters
- Conditional SAE boost application (gated by `enable_sae` flag)
- Provenance tracking for SAE attribution

**Key Insight**: SAE integration requires drug-specific biomarker mapping and careful confidence modulation.

---

### **6. Modal Payload Size Optimization** (Lines 28000-28500)

**Problem**: Modal SAE service crashing due to massive JSON payload

**Root Cause**:
- Line 346: `"features": features.cpu().numpy().tolist()` was trying to serialize a **32,768-dimensional vector** (potentially across 8,193 sequence positions = **268 million floats**) as JSON
- This is:
  - **Massive payload** (~1-2GB of JSON)
  - **Crashing Modal** (memory/serialization timeout)
  - **Useless for downstream** (we only need the top-k features for biomarker analysis)

**Fix Applied**:
```python
# Removed full features array from response
# Only return top_features (the 64 most active features)
return {
    "top_features": [
        {"index": int(idx), "value": float(val)}
        for idx, val in zip(top_indices, top_values)
    ],
    "stats": {...},
    "provenance": {...}
}
```

**Result**:
- Payload reduced from ~1-2GB to ~1KB
- No crashes, no timeouts
- Downstream analysis only needs top-k anyway

**Key Insight**: Always return only what's needed, not full feature vectors.

---

## 📊 **KEY INSIGHTS**

### **Assembly Version Handling**
1. **TCGA Data**: Always GRCh37/hg19
2. **Clinical NGS**: Usually GRCh38/hg38
3. **Always Specify**: Assembly parameter in API calls
4. **Liftover**: May be needed for cross-assembly comparison

### **Modal Service Architecture**
1. **Warm Containers**: Stay alive with traffic
2. **Cost Controls**: Required for RUO jobs
3. **Gating**: Use feature flags and caps
4. **Manual Stop**: Required to force code updates

### **Biomarker Integration**
1. **Drug-Specific**: Map biomarkers to drugs by mechanism
2. **Confidence Modulation**: SAE boost capped at ±15%
3. **Provenance**: Track all SAE attribution
4. **Gating**: Require explicit flags for production use

---

## 🔗 **CONTEXT & CONNECTIONS**

- **Builds on**: Bug fixes (Piece 3.3, Piece 5.6), Modal deployment (Piece 2.5)
- **Enables**: Biomarker correlation analysis, WIWFM integration
- **Blocks**: Real SAE feature extraction until assembly fix
- **Resolved**: All issues identified and fixes documented

---

## 📝 **NOTES**

- Biomarker correlation service is fully implemented (688 lines)
- Assembly mismatch was critical bug preventing extraction
- Modal architecture requires explicit cost controls
- Payload optimization was essential for service stability
- WIWFM integration architecture is fully designed but not yet implemented

---

## 🎯 **QUESTIONS RESOLVED**

- ✅ How to implement biomarker correlation? → Full service with statistical methods
- ✅ Why Ensembl 400 errors? → GRCh37 vs GRCh38 assembly mismatch
- ✅ How to prevent Modal cost overruns? → Add gating, caps, and kill switches
- ✅ How to integrate SAE into WIWFM? → Drug-specific mapping + confidence boost
- ✅ Why Modal service crashing? → Massive payload from full feature array








