# Enhanced Biomarker Extraction Specification

## Mission

**Goal**: Expand biomarker extraction from mutation data to improve coverage from 4.7% → 15-20% (HRD) and 1.3% → 3-5% (MSI), enabling meaningful sporadic gate application in benchmarks.

**Current State** (After Audit):
- HRD: Only BRCA1/2 → 4.7% coverage (56/1185 patients)
- MSI: Only core MMR → 1.3% coverage (15/1185 patients)
- TMB: Mutation count only (low quality)

**Target State**:
- HRD: Full HRR pathway → 15-20% coverage (~180-240 patients)
- MSI: Expanded MMR panel → 3-5% coverage (~35-60 patients)
- TMB: Improved calculation with better filtering

---

## Task 1: Enhanced HRD Score Estimation

### Current Implementation

**File**: `scripts/benchmark/benchmark_common/utils/biomarker_extractor.py`

**Current Logic** (lines 88-97):
```python
# Only checks BRCA1/2
brca_mutations = [
    mut for mut in mutations
    if mut.get("gene", "").upper() in ["BRCA1", "BRCA2"] 
    and mut.get("variant_type") not in ["Silent", "Intron"]
]
if brca_mutations:
    return 50.0, "estimated_brca"
return None, None
```

**Coverage**: 56/1185 (4.7%)

### Enhanced Implementation

**Expand to Full HRR Pathway**:

```python
# Core HRR pathway genes (High confidence for HRD estimation)
HRR_GENES_CORE = {
    "BRCA1", "BRCA2",           # Current: 4.7%
    "PALB2",                     # BRCA2 binding partner
    "RAD51C", "RAD51D",          # RAD51 paralogs (HR restoration)
    "BRIP1",                     # BRCA1 binding partner
    "BARD1",                     # BRCA1 binding partner
}

# Extended HRR pathway (Medium confidence)
HRR_GENES_EXTENDED = {
    "ATM",                       # DNA damage sensor
    "CHEK2",                     # DNA damage checkpoint
    "FANCA", "FANCC", "FANCD2",  # Fanconi anemia pathway
    "RAD50", "MRE11", "NBN",     # MRN complex
}

# Estimation Logic:
# - Core HRR mutation → HRD = 50-60 (high confidence)
# - Extended HRR mutation → HRD = 40-50 (medium confidence)
# - Multiple HRR mutations → HRD = 60-70 (biallelic loss)
```

**Expected Coverage**: 
- Core HRR: ~10-12% (120-140 patients)
- Extended HRR: ~5-8% (60-95 patients)
- **Total: 15-20%** (~180-240 patients)

### Implementation Details

**Function Signature**:
```python
def extract_hrd_from_patient(
    patient: dict, 
    estimate_hrd: bool = True,  # Enable estimation
    confidence_threshold: str = "medium"  # "high" (core only) or "medium" (core + extended)
) -> Tuple[Optional[float], Optional[str]]:
```

**Logic**:
1. Try direct HRD field first (unchanged)
2. If `estimate_hrd=True`:
   - Check for core HRR mutations → return 55.0, "estimated_hrr_core"
   - Check for extended HRR mutations → return 45.0, "estimated_hrr_extended"
   - Check for multiple HRR mutations → return 65.0, "estimated_hrr_biallelic"
3. Return `None, None` if no HRR mutations found

**Variant Type Filtering**:
- Include: Missense, Nonsense, Frameshift, Splice, InDel
- Exclude: Silent, Intron, Synonymous, UTR

**Confidence Levels**:
- `"high"`: Core HRR only (BRCA1/2/PALB2/RAD51C/RAD51D/BRIP1/BARD1)
- `"medium"`: Core + Extended (adds ATM/CHEK2/FANCA/etc.)

---

## Task 2: Enhanced MSI Status Estimation

### Current Implementation

**Current Logic** (lines 120-140):
```python
# Only checks core MMR genes
mmr_genes = ["MLH1", "MSH2", "MSH6", "PMS2"]
```

**Coverage**: 15/1185 (1.3%)

### Enhanced Implementation

**Expand to Full MMR Panel**:

```python
# Core MMR genes (High confidence for MSI-H)
MMR_GENES_CORE = {
    "MLH1", "MSH2", "MSH6", "PMS2"  # Current: 1.3%
}

# Extended MMR genes (Medium confidence)
MMR_GENES_EXTENDED = {
    "PMS1",                     # MMR pathway
    "MLH3",                     # MLH1 paralog
    "MSH3",                     # MSH2 binding partner
    "EXO1",                     # MMR exonuclease
    "POLD1", "POLE",            # DNA polymerase (proofreading)
}

# Estimation Logic:
# - Core MMR mutation → MSI-H (high confidence)
# - Extended MMR mutation → MSI-H (medium confidence)
# - Multiple MMR mutations → MSI-H (stronger signal)
```

**Expected Coverage**:
- Core MMR: ~1.3% (15 patients) - current
- Extended MMR: ~2-4% (24-48 patients)
- **Total: 3-5%** (~35-60 patients)

### Implementation Details

**Function Signature**:
```python
def extract_msi_from_patient(
    patient: dict,
    estimate_msi: bool = True,  # Enable estimation
    confidence_threshold: str = "medium"  # "high" (core only) or "medium" (core + extended)
) -> Tuple[Optional[Literal["MSI-H", "MSS"]], Optional[str]]:
```

**Logic**:
1. Try direct MSI field first (unchanged)
2. If `estimate_msi=True`:
   - Check for core MMR mutations → return "MSI-H", "estimated_mmr_core"
   - Check for extended MMR mutations → return "MSI-H", "estimated_mmr_extended"
   - If no MMR mutations → return "MSS", "estimated_no_mmr" (optional - may return None)
3. Return `None, None` if no direct field and no MMR mutations

**Variant Type Filtering**:
- Include: Missense (pathogenic), Nonsense, Frameshift, Splice, InDel
- Exclude: Silent, Intron, Synonymous, UTR

**Confidence Levels**:
- `"high"`: Core MMR only (MLH1/MSH2/MSH6/PMS2)
- `"medium"`: Core + Extended (adds PMS1/MLH3/MSH3/EXO1/POLD1/POLE)

---

## Task 3: Improved TMB Calculation

### Current Implementation

**Current Logic** (lines 40-70):
```python
# Simple mutation count
tmb = len(mutations) / 30.0  # Divide by 30 Mb (approximate exome size)
```

**Issues**:
- Includes all mutations (synonymous, intronic, UTR)
- No filtering for pathogenic variants
- Low quality estimates

### Enhanced Implementation

**Better Mutation Filtering**:

```python
def extract_tmb_from_patient(
    patient: dict,
    use_pathogenic_only: bool = True,  # Filter to pathogenic variants
    min_maf: float = 0.01  # Exclude common variants (MAF > 1%)
) -> Tuple[Optional[float], Optional[str]]:
```

**Filtering Logic**:
1. **Variant Type Filtering**:
   - Include: Missense, Nonsense, Frameshift, Splice, InDel
   - Exclude: Silent, Intron, Synonymous, UTR, 3'UTR, 5'UTR

2. **Pathogenicity Filtering** (if `use_pathogenic_only=True`):
   - Include: Variants with ClinVar "Pathogenic" or "Likely Pathogenic"
   - Include: Variants with VEP "HIGH" or "MODERATE" impact
   - Exclude: Variants with ClinVar "Benign" or "Likely Benign"

3. **MAF Filtering** (if `min_maf` provided):
   - Exclude: Variants with gnomAD MAF > `min_maf` (common variants)

4. **Calculation**:
   - TMB = (filtered_mutation_count) / 30.0  # mutations per Mb
   - Source: "estimated_from_mutations_filtered"

**Expected Improvement**:
- Current: All mutations → TMB = 2-5 mutations/Mb (low quality)
- Enhanced: Pathogenic only → TMB = 0.5-2 mutations/Mb (higher quality, more realistic)

---

## Task 4: Update `build_tumor_context()` Function

### Current Implementation

**Current Logic** (lines 150-180):
```python
tmb, tmb_source = extract_tmb_from_patient(patient)
hrd, hrd_source = extract_hrd_from_patient(patient)
msi, msi_source = extract_msi_from_patient(patient)
```

### Enhanced Implementation

**Add Configuration Parameters**:

```python
def build_tumor_context(
    patient: dict,
    estimate_hrd: bool = True,
    estimate_msi: bool = True,
    hrd_confidence: str = "medium",  # "high" or "medium"
    msi_confidence: str = "medium",  # "high" or "medium"
    tmb_pathogenic_only: bool = True,
    tmb_min_maf: float = 0.01
) -> Dict[str, Any]:
```

**Logic**:
1. Call enhanced extraction functions with configuration parameters
2. Build `tumor_context` dict (unchanged structure)
3. Include `biomarker_sources` dict (unchanged)
4. Add `biomarker_confidence` dict (NEW):
   ```python
   biomarker_confidence = {
       "hrd": "high" if hrd_source == "direct" else "estimated",
       "msi": "high" if msi_source == "direct" else "estimated",
       "tmb": "high" if tmb_source == "direct" else "estimated"
   }
   ```

---

## Task 5: Validation Script

### Create Validation Script

**File**: `scripts/benchmark/validate_enhanced_biomarkers.py`

**Purpose**: Validate enhanced extraction improves coverage without introducing errors

**Tests**:
1. **Coverage Test**: Verify HRD coverage increased from 4.7% → 15-20%
2. **MSI Coverage Test**: Verify MSI coverage increased from 1.3% → 3-5%
3. **TMB Quality Test**: Verify TMB values are more realistic (lower, but higher quality)
4. **Source Tracking Test**: Verify all biomarkers have source tracking
5. **Confidence Tracking Test**: Verify confidence levels are tracked correctly

**Expected Output**:
```
Enhanced Biomarker Extraction Validation
========================================
HRD Coverage:
  - Before: 56/1185 (4.7%)
  - After: 180/1185 (15.2%) ✅

MSI Coverage:
  - Before: 15/1185 (1.3%)
  - After: 42/1185 (3.5%) ✅

TMB Quality:
  - Before: Mean 3.2 mutations/Mb (all mutations)
  - After: Mean 1.1 mutations/Mb (pathogenic only) ✅

Source Tracking: 100% ✅
Confidence Tracking: 100% ✅
```

---

## Acceptance Criteria

### Coverage Targets

| Biomarker | Current | Target | Status |
|-----------|---------|--------|--------|
| **HRD** | 4.7% (56/1185) | 15-20% (180-240) | ⏸️ Pending |
| **MSI** | 1.3% (15/1185) | 3-5% (35-60) | ⏸️ Pending |
| **TMB** | 100% (low quality) | 100% (high quality) | ⏸️ Pending |

### Quality Requirements

1. **Source Tracking**: All biomarkers must have `biomarker_sources` dict
2. **Confidence Tracking**: All biomarkers must have `biomarker_confidence` dict
3. **Validation**: Coverage improvements verified by validation script
4. **No Regression**: Direct biomarker fields still work (unchanged)
5. **Backward Compatible**: Default behavior unchanged (estimation disabled by default)

### Code Quality

1. **Type Hints**: All functions must have proper type hints
2. **Docstrings**: All functions must have docstrings explaining logic
3. **Error Handling**: Graceful handling of missing data
4. **Testing**: Validation script passes all tests

---

## Implementation Order

### Step 1: Enhanced HRD Extraction (2 hours)

1. Update `extract_hrd_from_patient()` with expanded gene lists
2. Add confidence threshold parameter
3. Implement multi-mutation detection (biallelic loss)
4. Test with validation script

### Step 2: Enhanced MSI Extraction (1 hour)

1. Update `extract_msi_from_patient()` with expanded gene lists
2. Add confidence threshold parameter
3. Test with validation script

### Step 3: Improved TMB Calculation (1 hour)

1. Update `extract_tmb_from_patient()` with filtering logic
2. Add pathogenic-only and MAF filtering
3. Test with validation script

### Step 4: Update `build_tumor_context()` (30 min)

1. Add configuration parameters
2. Add confidence tracking
3. Test end-to-end

### Step 5: Create Validation Script (1 hour)

1. Create `validate_enhanced_biomarkers.py`
2. Implement all 5 tests
3. Run on full dataset

**Total Time**: ~5-6 hours

---

## Gene Lists Reference

### HRR Pathway Genes (HRD Estimation)

**Core (High Confidence)**:
- BRCA1, BRCA2, PALB2, RAD51C, RAD51D, BRIP1, BARD1

**Extended (Medium Confidence)**:
- ATM, CHEK2, FANCA, FANCC, FANCD2, RAD50, MRE11, NBN

**Sources**:
- COSMIC: HRD gene panel
- TCGA: HRD score associations
- Literature: HRR pathway reviews

### MMR Pathway Genes (MSI Estimation)

**Core (High Confidence)**:
- MLH1, MSH2, MSH6, PMS2

**Extended (Medium Confidence)**:
- PMS1, MLH3, MSH3, EXO1, POLD1, POLE

**Sources**:
- COSMIC: MSI gene panel
- TCGA: MSI status associations
- Literature: MMR pathway reviews

---

## Questions for Manager

1. **Confidence Threshold Default**: Should default be "high" (conservative) or "medium" (more coverage)?
   - **Recommendation**: "medium" (maximize coverage for benchmarks)

2. **TMB Pathogenic Filtering**: Should we filter to pathogenic variants only?
   - **Recommendation**: Yes (improves quality, more realistic TMB values)

3. **MSI MSS Assignment**: Should we assign "MSS" when no MMR mutations found?
   - **Recommendation**: No (return None - we don't know MSI status without direct test)

4. **Biallelic Loss Detection**: Should we detect multiple HRR mutations as biallelic loss?
   - **Recommendation**: Yes (higher HRD score = 65.0 for biallelic loss)

---

## Success Metrics

### Coverage Improvement

- **HRD**: 4.7% → 15-20% (3-4x improvement)
- **MSI**: 1.3% → 3-5% (2-3x improvement)
- **TMB**: Quality improvement (pathogenic filtering)

### Benchmark Impact

- **Sporadic Gates**: Will apply to 15-20% of patients (vs 4.7% before)
- **IO Boost**: TMB-high patients will get 1.35x boost
- **PARP Rescue**: HRD-high patients will get full PARP effect
- **Expected Correlation**: r=0.037 → r=0.10-0.15 (modest improvement, but meaningful)

---

**Status**: ⏸️ **AWAITING AGENT IMPLEMENTATION**

**Next Step**: Agent implements enhanced extraction → Zo validates → Re-run audit → Proceed with benchmarks

