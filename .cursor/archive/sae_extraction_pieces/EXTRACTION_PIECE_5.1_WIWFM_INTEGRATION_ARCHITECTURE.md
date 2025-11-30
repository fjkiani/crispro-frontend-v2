# EXTRACTION PIECE 5.1: WIWFM Integration Architecture

**Source**: Lines 25040-25290 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-XX  
**Status**: ✅ Complete

---

## Overview

This piece documents the complete architecture for integrating SAE biomarkers into the WIWFM (Will It Work For Me) drug efficacy prediction pipeline. This is Phase 4 of the SAE implementation, currently **PENDING** awaiting Manager validation approval.

---

## Current State vs Planned Architecture

### Current State
- SAE features computed but **NOT used** to modulate drug efficacy scores
- SAE is "display only" (shows in UI but doesn't affect recommendations)
- **Manager's Policy**: Wait for validation + written SAE policy before integrating

### Planned Architecture

```
Drug Efficacy Pipeline (WIWFM):
  
  Ayesha's Mutations → Evo2 → Insights (functionality, essentiality, chromatin)
                           ↓
                      Pathway Scores (DDR, MAPK, PI3K, VEGF, HER2)
                           ↓
                      Evidence (literature, trials)
                           ↓
                      S/P/E Base Score (0-1)
                           ↓
                   ╔═══════╩═══════╗
                   ║ SAE MODULE   ║  ← NEW INTEGRATION
                   ║ (Phase 3)    ║
                   ╚═══════╤═══════╝
                           ↓
              Extract Patient SAE Features (32K-dim)
                           ↓
              Map to Drug-Specific Biomarkers
                           ↓
              Compute SAE Boost/Penalty (-0.15 to +0.15)
                           ↓
              Final Confidence = Base + SAE_Boost (capped at 0.95)
```

---

## Step 1: Extract Patient SAE Features

**Function**: `extract_patient_sae_features(mutations: List[Dict]) -> np.ndarray`

**Purpose**: Extract SAE features for all patient mutations and aggregate into a single 32K-dimensional vector.

**Implementation**:
- For each mutation, call Modal SAE service (`/extract_features`)
- Pass: `chrom`, `pos`, `ref`, `alt`, `assembly` (GRCh38 for Ayesha)
- Aggregate using mean pooling across all mutations
- Returns: `(32768,)` numpy array

**File**: `api/services/sae_biomarker_drug_mapper.py` (to be created)

---

## Step 2: Drug-Specific Biomarker Mapping

**Class**: `SAEBiomarkerDrugMapper`

**Purpose**: Map drugs to relevant SAE biomarkers based on mechanism of action.

**Drug Mappings**:
- **Platinum agents** (Carboplatin, Cisplatin): `biomarker_weight=1.0` (direct platinum response correlation)
- **PARP inhibitors** (Olaparib, Niraparib): `biomarker_weight=0.6` (indirect DNA repair proxy)
- **MEK/RAF inhibitors** (Trametinib): `biomarker_weight=0.0` (no platinum correlation)

**Feature Selection**:
- Platinum: Top 20 features weighted by correlation strength
- PARP: Top 15 features with 0.6 weight
- MEK: No features (use hotspot detection instead)

**Method**: `get_top_features(n: int, weight: float) -> Dict[int, float]`
- Extracts top N features from biomarker correlation analysis
- Weights by `pearson_r * weight`

---

## Step 3: Compute Drug-Specific SAE Score

**Function**: `compute_patient_sae_score(drug_name, patient_sae_features, biomarker_mapper) -> float`

**Purpose**: Compute drug-specific SAE score for patient.

**Algorithm**:
1. Get drug mapping (feature weights + biomarker weight)
2. Weighted sum: `score = Σ(patient_features[feat_idx] * feat_weight)`
3. Normalize using `tanh(score / len(feature_weights))` → `[-1, +1]`
4. Apply drug-specific biomarker weight
5. Returns: SAE score in `[-1.0, +1.0]` range (0 = neutral)

---

## Step 4: Apply SAE Boost to WIWFM Confidence

**Function**: `apply_sae_biomarker_boost(drug_name, base_confidence, patient_mutations, biomarker_mapper) -> Tuple[float, Dict]`

**Purpose**: Apply SAE biomarker boost to WIWFM confidence.

**Algorithm**:
1. Extract patient SAE features (Step 1)
2. Compute drug-specific SAE score (Step 3)
3. Convert to confidence boost: `sae_boost = sae_score * 0.15` → `[-0.15, +0.15]`
4. Apply boost: `boosted_confidence = min(base_confidence + sae_boost, 0.95)`
5. Return: `(boosted_confidence, provenance_dict)`

**Provenance Includes**:
- `sae_score`: Raw SAE score
- `sae_boost`: Applied boost amount
- `base_confidence`: Original S/P/E confidence
- `boosted_confidence`: Final confidence
- `biomarker_features_used`: Number of features used

---

## Integration into Drug Scorer

**File**: `api/services/efficacy_orchestrator/drug_scorer.py`

**New Parameters**:
- `patient_mutations: Optional[List[Dict]] = None`
- `sae_biomarker_mapper: Optional[SAEBiomarkerDrugMapper] = None`
- `enable_sae: bool = False`

**Integration Point**:
- After computing base S/P/E confidence
- Before returning `DrugScoreResult`
- Only if `enable_sae=True` and all parameters provided

**Code Pattern**:
```python
# Compute base confidence
confidence = compute_confidence(tier, seq_pct, path_pct, insights_dict, confidence_config)

# NEW: Apply SAE biomarker boost if enabled
if enable_sae and patient_mutations and sae_biomarker_mapper:
    confidence, sae_provenance = apply_sae_biomarker_boost(
        drug["name"],
        confidence,
        patient_mutations,
        sae_biomarker_mapper
    )
    
    # Add SAE provenance to result
    drug_result.provenance["sae_boost"] = sae_provenance

return drug_result
```

---

## Status

**Current Status**: ⏸️ **PENDING** (Awaiting Validation)

**Blocking Factors**:
1. Manager policy: Wait for validation + written SAE policy
2. Validation required: AUROC/AUPRC computed on ≥200 patients
3. Feature flag: `ENABLE_SAE_BIOMARKERS=true` (default: false)

**Design Status**: ✅ **COMPLETE**
- File: `.cursor/rules/AYESHA_SAE_WIWFM_INTEGRATION_PLAN.mdc` (517 lines)
- All 4 steps designed with code examples
- Integration points identified

**Implementation Files** (To be created after validation approval):
- `api/services/sae_biomarker_drug_mapper.py` (~300 lines)
- `api/services/efficacy_orchestrator/sae_booster.py` (~150 lines)
- `src/components/ayesha/SAEBiomarkerCard.jsx` (~100 lines)

---

## Key Design Decisions

1. **Confidence Boost Range**: ±15% maximum (never >95% total)
2. **Drug-Specific Weighting**: Platinum (1.0), PARP (0.6), MEK (0.0)
3. **Feature Selection**: Top-N features weighted by correlation strength
4. **Normalization**: `tanh()` to bound scores in `[-1, +1]`
5. **Provenance**: Full transparency on boost/penalty reasoning

---

## Related Documents

- `.cursor/ayesha/ZO_SAE_WIWFM_INTEGRATION_REVIEW.md` - Current state analysis
- `.cursor/rules/AYESHA_SAE_WIWFM_INTEGRATION_PLAN.mdc` - Complete integration plan
- `.cursor/ayesha/SAE_LIFT_GATE_POLICY_V1.md` - Manager-approved policy

