# ⚔️ PHASE 2: INTEGRATION POINTS ANALYSIS ⚔️

**Date**: January 15, 2025  
**Status**: ✅ **COMPLETE**  
**Objective**: Map all integration points between components with code evidence

---

## Task 2.1: SAE→WIWFM Integration Status

### Current State (Verified in Code)

**File**: `api/services/efficacy_orchestrator/orchestrator.py:342-389`

**SAE Extraction**:
```python
# Line 342: SAE extraction is OPTIONAL (requires include_sae_features flag)
if (request.options or {}).get("include_sae_features"):
    # Line 347: Extract SAE features from real data
    sae_bundle = extract_sae_features_from_real_data(...)
    # Line 378: Add to response (DISPLAY ONLY)
    response.sae_features = sae_features_to_dict(sae_bundle)
    # Line 380-386: Add attribution to confidence_breakdown (DISPLAY ONLY)
    response.provenance["confidence_breakdown"]["sae_attribution"] = {
        "boosting_features": sae_bundle.boosting_features,
        "limiting_features": sae_bundle.limiting_features,
        "overall_impact": sae_bundle.overall_impact
    }
```

**SAE NOT Modulating Confidence**:
- **File**: `api/services/efficacy_orchestrator/drug_scorer.py` (searched entire file)
- **Finding**: NO SAE references found
- **Evidence**: Confidence computed from S/P/E only (line 171: `raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd`)

### Gap Analysis

**Manager's Vision** (from `.cursor/rules/AYESHA_SAE_WIWFM_INTEGRATION_PLAN.mdc`):
- "SAE must live inside S/P/E and modulate confidence"
- "S/P/E base ± SAE lifts/penalties"

**Current Reality**:
- SAE extracted and displayed
- SAE attribution shown in confidence_breakdown
- **BUT**: SAE does NOT modulate confidence scores

**Blocking Factors**:
1. Manager policy: "DO NOT INTEGRATE SAE INTO EFFICACY YET"
2. Validation required: ≥200 TCGA patients, AUROC/AUPRC computed
3. Written SAE policy approval needed

---

## Task 2.2: S/P/E→SAE Data Flow

### Data Flow Mapping (Verified in Code)

**Step 1: Evo2 Sequence Scoring → SAE**
- **File**: `orchestrator.py:95` → `sequence_processor.score_sequences()`
- **Output**: `seq_scores[0].sequence_disruption`, `seq_scores[0].calibrated_seq_percentile`
- **SAE Input**: Lines 350-352 - Uses `sequence_disruption` and `calibrated_seq_percentile`

**Step 2: Pathway Aggregation → SAE**
- **File**: `orchestrator.py:105` → `aggregate_pathways()`
- **Output**: `pathway_scores` (Dict[str, float])
- **SAE Input**: Line 360 - Uses `pathway_disruption=pathway_scores`

**Step 3: Insights Bundle → SAE**
- **File**: `orchestrator.py:161` → `bundle_insights()`
- **Output**: `InsightsBundle` (functionality, chromatin, essentiality, regulatory)
- **SAE Input**: Lines 354-359 - Uses all 4 insights chips

**Step 4: Evidence → SAE**
- **File**: `orchestrator.py:112-150` → `literature()`, `clinvar_prior()`
- **Output**: `evidence_results`, `clinvar_result`
- **SAE Input**: Lines 362-370 - Uses ClinVar data (classification)

**Step 5: SAE Feature Extraction**
- **File**: `api/services/sae_service.py` → `extract_sae_features_from_real_data()`
- **Input**: All S/P/E outputs (Evo2 scores, pathway scores, insights, ClinVar)
- **Output**: `SAEFeatures` bundle (6 core features)

**Step 6: SAE → Confidence (CURRENTLY DISABLED)**
- **File**: `drug_scorer.py` - NO SAE integration
- **Gap**: SAE features not passed to `score_drug()` method
- **Gap**: No SAE confidence modulation logic

### Complete Data Flow Diagram

```
User Input (mutations)
  ↓
Evo2 Scoring (orchestrator.py:95)
  ├─→ sequence_disruption
  ├─→ calibrated_seq_percentile
  └─→ SAE Input (line 350-352)
  ↓
Pathway Aggregation (orchestrator.py:105)
  ├─→ pathway_scores (DDR, MAPK, PI3K, VEGF)
  └─→ SAE Input (line 360)
  ↓
Insights Bundle (orchestrator.py:161)
  ├─→ functionality, chromatin, essentiality, regulatory
  └─→ SAE Input (line 354-359)
  ↓
Evidence Gathering (orchestrator.py:112-150)
  ├─→ literature strength
  ├─→ ClinVar classification
  └─→ SAE Input (line 362-370)
  ↓
SAE Feature Extraction (orchestrator.py:347)
  ├─→ dna_repair_capacity
  ├─→ pathway_burden
  ├─→ hotspot_mutation
  ├─→ essentiality_hrr_genes
  ├─→ exon_disruption_score
  └─→ mechanism_vector
  ↓
SAE Display (orchestrator.py:378)
  ├─→ response.sae_features (DISPLAY ONLY)
  └─→ confidence_breakdown.sae_attribution (DISPLAY ONLY)
  ↓
Drug Scoring (drug_scorer.py:171)
  ├─→ NO SAE INPUT ❌
  ├─→ confidence = compute_confidence(S, P, E, insights)
  └─→ SAE NOT MODULATING ❌
```

---

## Task 2.3: Ayesha Orchestrator Integration

### Complete End-to-End Flow (Verified in Code)

**File**: `api/routers/ayesha_orchestrator_v2.py:324-671`

**Flow 1: Pre-NGS (No Tumor Context)**
```
POST /api/ayesha/complete_care_v2
  ↓
1. Clinical Trials (line 350-360)
   └─→ POST /api/ayesha/trials/search
   └─→ Returns: top 10 trials, SOC recommendation
  ↓
2. CA-125 Intelligence (line 362-370)
   └─→ ca125_intelligence.analyze_ca125()
   └─→ Returns: burden, forecast, resistance signals
  ↓
3. Drug Efficacy (line 372-380)
   └─→ POST /api/efficacy/predict
   └─→ Returns: "awaiting_ngs" message
  ↓
4. NGS Fast-Track (line 382-390)
   └─→ ngs_fast_track.generate_checklist()
   └─→ Returns: ctDNA, HRD, IHC orders
  ↓
Response: Complete care plan (pre-NGS)
```

**Flow 2: Post-NGS (Tumor Context Available)**
```
POST /api/ayesha/complete_care_v2
  ↓
1. Clinical Trials (same as Flow 1)
  ↓
2. CA-125 Intelligence (same as Flow 1)
  ↓
3. Drug Efficacy (line 372-380)
   └─→ POST /api/efficacy/predict
   └─→ WITH tumor_context
   └─→ Returns: Ranked drugs with S/P/E scores
   └─→ SAE features extracted (if include_sae_features=true)
  ↓
4. SAE Features (line 392-400)
   └─→ compute_sae_features() (Phase 2)
   └─→ Returns: 6 core SAE features
  ↓
5. Resistance Playbook (line 402-410)
   └─→ POST /api/care/resistance_playbook
   └─→ WITH tumor_context + SAE features
   └─→ Returns: 5 detection rules, 7 combos, 6 switches
  ↓
6. Resistance Prophet (line 412-420, if enabled)
   └─→ get_resistance_prophet_service().predict_resistance()
   └─→ WITH current_sae_features + baseline_sae_features
   └─→ Returns: Risk level, signals, actions
  ↓
7. Food Validator (line 422-430, if food_query provided)
   └─→ POST /api/hypothesis/validate_food_dynamic
   └─→ Returns: SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED
  ↓
Response: Complete care plan (post-NGS)
```

### Integration Points

**1. WIWFM → SAE**
- **Location**: `orchestrator.py:342-389`
- **Status**: ✅ SAE extracted from WIWFM outputs
- **Gap**: SAE not modulating WIWFM confidence

**2. Resistance Playbook → SAE**
- **Location**: `ayesha_orchestrator_v2.py:402-410`
- **Status**: ✅ SAE features passed to resistance playbook
- **Usage**: DNA repair capacity used for HR restoration detection

**3. CA-125 Intelligence → Drug Recommendations**
- **Location**: `ayesha_orchestrator_v2.py:362-370`
- **Status**: ✅ CA-125 analyzed independently
- **Gap**: CA-125 not directly influencing drug ranking (separate service)

**4. Resistance Prophet → SAE**
- **Location**: `ayesha_orchestrator_v2.py:412-420`
- **Status**: ✅ Resistance Prophet uses SAE features
- **Usage**: DNA repair restoration + pathway escape detection

---

## Summary

### Integration Status Matrix

| Integration Point | Status | Code Evidence | Gap |
|------------------|--------|--------------|-----|
| S/P/E → SAE | ✅ Complete | orchestrator.py:342-389 | None |
| SAE → Confidence | ❌ Not Done | drug_scorer.py (no SAE) | CRITICAL |
| WIWFM → Resistance Playbook | ✅ Complete | ayesha_orchestrator_v2.py:402-410 | None |
| SAE → Resistance Prophet | ✅ Complete | ayesha_orchestrator_v2.py:412-420 | None |
| CA-125 → Drug Ranking | ⚠️ Indirect | Separate service | Minor |
| Food Validator → SAE | ✅ Complete | food_spe_integration.py:458-463 | None |

### Critical Finding

**SAE→WIWFM Integration Gap**: SAE features are extracted from S/P/E outputs but do NOT modulate drug efficacy confidence. This is the primary gap between current state and Manager's vision.




