# Final Comprehensive Plan: Mechanism-Based Trial Matching & Resistance Prediction

**Date**: January 28, 2025  
**Status**: ✅ **APPROVED** - Strategic Review Complete  
**Location**: Main Repository (`/Users/fahadkiani/Desktop/development/crispr-assistant-main/`)  
**Reviewed By**: Manager Agent (STRATEGIC_REVIEW_MECHANISM_TRIAL_MATCHING.md)

---

## Executive Summary

**Mission**: Deliver mechanism-based trial matching AND mechanism-based resistance/sensitivity prediction

**Current State**: 80% complete - All services exist, need wiring and enhancement (VERIFIED)

**Timeline**: 3-4 days (12-16 hours total, includes 30% buffer)

**Independence**: This work is **INDEPENDENT** of Zo2's outcome calibration. Proceed in parallel.

---

## Agent Coordination (Clarified)

| Agent | Task | Status | Dependency |
|-------|------|--------|------------|
| **Zo** | MBD4+TP53 Analysis | ✅ Complete | None |
| **Zo2** | S/P/E Calibration (r=0.037 → r>0.15) | 🔄 Phase 1 | None |
| **This Agent** | Trial Matching + Resistance | 📋 Approved | None (independent) |

**Key Insight**: Mechanism fit doesn't predict outcomes; it ranks trials by pathway alignment. This is different from Zo2's outcome prediction work.

---

## What EXISTS ✅

### Core Services (All Implemented - VERIFIED)

| Service | File | Lines | Status | Notes |
|---------|------|-------|--------|-------|
| **Mechanism Fit Ranker** | `api/services/mechanism_fit_ranker.py` | 275 | ✅ Complete | Formula: 0.7×eligibility + 0.3×mechanism_fit |
| **Resistance Prophet** | `api/services/resistance_prophet_service.py` | 689 | ✅ Complete | 2-of-3 signals, needs mechanism breakdown |
| **Pathway→Vector** | `api/services/pathway_to_mechanism_vector.py` | - | ✅ Complete | 7D vector conversion |
| **Autonomous Trial Agent** | `api/services/autonomous_trial_agent.py` | - | ✅ Complete | 10 query templates, DDR detection |
| **SAE Feature Service** | `api/services/sae_feature_service.py` | - | ✅ Complete | TRUE SAE flag ready |

### Data Resources

| Resource | File | Status | Notes |
|----------|------|--------|-------|
| **MoA Vectors** | `api/resources/trial_moa_vectors.json` | ✅ 47 trials tagged | See MoA Coverage section |
| **SAE Mapping** | `api/resources/sae_feature_mapping.json` | ✅ 288KB, 88 features | |

### Integration Points

| Integration | File | Status |
|-------------|------|--------|
| **Resistance Prophet** | `api/routers/ayesha_orchestrator_v2.py` | ✅ Integrated (opt-in flag) |
| **Trial Search** | `api/routers/trials_agent.py` | ⚠️ MechanismFitRanker not imported |

---

## What NEEDS WORK ⚠️

### Gap 1: Mechanism Fit NOT in Trial Response

**Current**: Trial search returns trials without mechanism fit scores  
**Needed**: Each trial has `mechanism_fit_score`, `combined_score`, `mechanism_alignment`

**File**: `api/routers/trials_agent.py`  
**Work**: 3-4 hours

### Gap 2: Resistance Prophet Needs Mechanism Breakdown

**Current**: Signals track pathway shifts in provenance, but no structured breakdown  
**Needed**: 
- `mechanism_breakdown` with DDR/HRR/exon changes
- `pathway_contributions` showing formula weights
- `escaped_pathways` list
- `mechanism_alignment` dict
- `baseline_source` documentation (when using population average)

**File**: `api/services/resistance_prophet_service.py`  
**Work**: 3-4 hours (adjusted from 2-3h)

### Gap 3: Validation Scripts Missing

**Needed**:
- `scripts/validate_mechanism_trial_matching.py`
- `scripts/validate_mechanism_resistance_prediction.py`
- `scripts/validate_mbd4_tp53_mechanism_capabilities.py`

**Work**: 6-8 hours (adjusted from 4-6h - validation always takes longer)

---

## MoA Coverage Report (NEW - Per Strategic Review)

### Current Coverage: 47 Trials Tagged

**Validation Required**:
```python
# MoA Coverage Check
{
    "total_trials_tagged": 47,
    "trials_with_complete_vectors": "TBD - verify 7D vectors",
    "trials_missing_moa": "TBD - list untagged",
    "pathway_coverage": {
        "DDR": "count",
        "MAPK": "count",
        "PI3K": "count",
        "VEGF": "count",
        "HER2": "count",
        "IO": "count",
        "Efflux": "count"
    }
}
```

**Fallback for Untagged Trials**: Eligibility-only ranking (no mechanism fit score)

**Action**: Include MoA coverage report in validation output

---

## Baseline Handling Documentation (NEW - Per Strategic Review)

### When Baseline SAE is Missing

**Default**: Use population average (0.50) with confidence penalty

**Required Output**:
```json
{
  "baseline_source": "population_average",
  "baseline_value": 0.50,
  "baseline_penalty_applied": true,
  "confidence_cap": "MEDIUM",
  "reason": "No baseline SAE available - using population average with confidence cap"
}
```

**Validation**: Ensure this is tested in `validate_mechanism_resistance_prediction.py`

---

## AYESHA Integration Flow (NEW - Per Strategic Review)

### How This Integrates with MBD4+TP53 Analysis

```
Patient Mutations → S/P/E Efficacy → Pathway Scores
                                          │
                                          ▼
                              Mechanism Vector (7D)
                                          │
                          ┌───────────────┴───────────────┐
                          │                               │
                          ▼                               ▼
                  Trial Matching                  Resistance Prediction
                  (This Agent)                    (This Agent)
                          │                               │
                          ▼                               ▼
                  Ranked Trials                   Risk Signals
                          │                               │
                          └───────────────┬───────────────┘
                                          │
                                          ▼
                              Clinical Decision Support
```

**Integration Test Required**: Validate full flow from mutations → trials + resistance

---

## Implementation Plan (3-4 Days, 12-16h Total)

### Day 1: Wire Mechanism Fit to Trial Search (3-4 hours)

#### Task 1.1: Add Mechanism Vector to Request Schema

**File**: `api/routers/trials_agent.py`

```python
class PatientDataRequest(BaseModel):
    # ... existing fields ...
    mechanism_vector: Optional[List[float]] = Field(None, description="7D mechanism vector [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]")
    efficacy_predictions: Optional[Dict[str, Any]] = Field(None, description="For auto-extraction of mechanism vector")
```

#### Task 1.2: Extract Mechanism Vector (Auto or Manual)

**File**: `api/routers/trials_agent.py`

```python
async def autonomous_trial_search(request: PatientDataRequest):
    # ... existing code ...
    
    # Extract mechanism vector (priority: explicit > efficacy > compute)
    mechanism_vector = request.mechanism_vector
    
    if not mechanism_vector and request.efficacy_predictions:
        # Auto-extract from efficacy predictions
        from api.services.pathway_to_mechanism_vector import convert_pathway_scores_to_mechanism_vector
        pathway_scores = request.efficacy_predictions.get('provenance', {}).get('confidence_breakdown', {}).get('pathway_disruption', {})
        if pathway_scores:
            mechanism_vector = convert_pathway_scores_to_mechanism_vector(
                pathway_scores,
                tmb=request.tumor_context.get('tmb') if request.tumor_context else None,
                msi_status=request.tumor_context.get('msi_status') if request.tumor_context else None
            )
```

#### Task 1.3: Apply Mechanism Fit Ranking

**File**: `api/routers/trials_agent.py`

```python
    # After agent.search_for_patient() call
    trials = results.get("matched_trials", [])
    
    if mechanism_vector and len(trials) > 0:
        from api.services.mechanism_fit_ranker import MechanismFitRanker
        ranker = MechanismFitRanker(alpha=0.7, beta=0.3)
        
        ranked_trials = ranker.rank_trials(
            trials=trials,
            sae_mechanism_vector=mechanism_vector,
            min_eligibility=0.60,
            min_mechanism_fit=0.50
        )
        
        # Mechanism fit scores are already in ranked_trials from ranker
        results["matched_trials"] = ranked_trials
        results["mechanism_fit_applied"] = True
        results["moa_coverage"] = {
            "trials_with_moa": sum(1 for t in ranked_trials if t.get('mechanism_fit_score', 0) > 0),
            "trials_without_moa": sum(1 for t in ranked_trials if t.get('mechanism_fit_score', 0) == 0)
        }
    else:
        results["mechanism_fit_applied"] = False
        if not mechanism_vector:
            results["mechanism_fit_warning"] = "No mechanism vector provided - ranking by eligibility only"
```

**Files to Modify**:
- `api/routers/trials_agent.py`

---

### Day 2: Enhance Resistance Prophet with Mechanism Breakdown (3-4 hours)

#### Task 2.1: Enhance DNA Repair Restoration Signal

**File**: `api/services/resistance_prophet_service.py`

**Current**: Returns `ResistanceSignalData` with provenance  
**Enhancement**: Add `mechanism_breakdown` to provenance

```python
async def _detect_dna_repair_restoration(
    self,
    current_sae: Dict,
    baseline_sae: Dict
) -> ResistanceSignalData:
    # ... existing code ...
    
    # Determine baseline source
    baseline_source = "patient_baseline"
    if not baseline_sae or baseline_sae.get("is_population_average", False):
        baseline_source = "population_average"
        baseline_penalty = True
    else:
        baseline_penalty = False
    
    # ADD: Mechanism breakdown
    ddr_change = current_sae.get("pathway_burden_ddr", 0.0) - baseline_sae.get("pathway_burden_ddr", 0.0)
    hrr_change = current_sae.get("essentiality_hrr", 0.0) - baseline_sae.get("essentiality_hrr", 0.0)
    exon_change = current_sae.get("exon_disruption_score", 0.0) - baseline_sae.get("exon_disruption_score", 0.0)
    
    provenance = {
        # ... existing provenance ...
        "mechanism_breakdown": {
            "ddr_pathway_change": float(ddr_change),
            "hrr_essentiality_change": float(hrr_change),
            "exon_disruption_change": float(exon_change)
        },
        "pathway_contributions": {
            "ddr": 0.60,  # From Manager's C1 formula
            "hrr": 0.20,
            "exon": 0.20
        },
        "baseline_source": baseline_source,
        "baseline_penalty_applied": baseline_penalty,
        "confidence_cap": "MEDIUM" if baseline_penalty else None
    }
    
    return ResistanceSignalData(...)
```

#### Task 2.2: Enhance Pathway Escape Signal

**File**: `api/services/resistance_prophet_service.py`

```python
async def _detect_pathway_escape(
    self,
    current_sae: Dict,
    baseline_sae: Dict,
    drug_class: Optional[str]
) -> ResistanceSignalData:
    # ... existing code ...
    
    # ADD: Identify escaped pathways
    pathway_names = ["DDR", "MAPK", "PI3K", "VEGF", "HER2", "IO", "Efflux"]
    escaped_pathways = []
    mechanism_alignment = {}
    
    for i, pathway in enumerate(pathway_names):
        change = current_vector[i] - baseline_vector[i]
        mechanism_alignment[pathway] = float(change)
        
        # If drug targets this pathway but pathway burden dropped, it's escape
        if self._drug_targets_pathway(drug_class, pathway) and change < -0.15:
            escaped_pathways.append(pathway)
    
    provenance = {
        # ... existing provenance ...
        "escaped_pathways": escaped_pathways,
        "mechanism_alignment": mechanism_alignment
    }
    
    return ResistanceSignalData(...)
```

**Files to Modify**:
- `api/services/resistance_prophet_service.py`

---

### Day 3-4: Create Validation Scripts (6-8 hours)

#### Task 3.1: Mechanism-Based Trial Matching Validation

**File**: `scripts/validate_mechanism_trial_matching.py` (NEW)

**8-Task Verification**:
1. Trial data quality (47 MoA-tagged trials exist)
2. Mechanism vector structure (7D format correct)
3. Mechanism fit computation (cosine similarity correct)
4. Combined score formula (0.7×eligibility + 0.3×mechanism_fit)
5. Ranking accuracy (Top-3 ≥0.70, MRR ≥0.65) - **MVP thresholds**
6. Pathway alignment (DDR trials rank higher for DDR-high patients)
7. Edge cases (all-zero vector, missing MoA vectors)
8. Consistency (deterministic results)
9. **MoA Coverage Report** (NEW - per strategic review)

**Metrics to Report**:
```python
{
    "top3_accuracy": 0.75,  # MVP Target: ≥0.70 (stretch: ≥0.80)
    "mrr": 0.68,  # MVP Target: ≥0.65 (stretch: ≥0.75)
    "mechanism_fit_scores": {
        "mean": 0.72,
        "std": 0.15,
        "min": 0.50,
        "max": 0.95
    },
    "pathway_alignment": {
        "DDR": 0.88,
        "MAPK": 0.75,
        "PI3K": 0.68
    },
    "moa_coverage": {
        "total_trials_tagged": 47,
        "complete_vectors": 45,
        "missing_moa": 2,
        "pathway_breakdown": {...}
    }
}
```

#### Task 3.2: Mechanism-Based Resistance Prediction Validation

**File**: `scripts/validate_mechanism_resistance_prediction.py` (NEW)

**8-Task Verification**:
1. Signal detection logic (DNA repair, pathway escape)
2. Mechanism breakdown accuracy (DDR 60%, HRR 20%, exon 20%)
3. Risk stratification (HIGH/MEDIUM/LOW thresholds)
4. Signal fusion (2-of-3 logic)
5. Pathway escape detection (escaped pathways identified)
6. **Baseline handling** (population average when missing, penalty applied)
7. Confidence modulation (MEDIUM cap when CA-125 missing or baseline missing)
8. Consistency (deterministic)

**Metrics to Report**:
```python
{
    "signal_detection_accuracy": {
        "dna_repair_restoration": 0.80,
        "pathway_escape": 0.72,
        "ca125_kinetics": 0.65
    },
    "risk_stratification": {
        "high_risk_auroc": 0.68,  # MVP Target: ≥0.65 (stretch: ≥0.70)
        "sensitivity": 0.75,
        "specificity": 0.68
    },
    "mechanism_breakdown_accuracy": {
        "ddr_contribution": 0.60,
        "hrr_contribution": 0.20,
        "exon_contribution": 0.20
    },
    "baseline_handling": {
        "population_average_used": 5,
        "patient_baseline_used": 15,
        "confidence_cap_applied": 5
    }
}
```

#### Task 3.3: MBD4+TP53 End-to-End Integration Validation (NEW - Per Strategic Review)

**File**: `scripts/validate_mbd4_tp53_mechanism_capabilities.py` (NEW)

**Full Integration Test**:
1. Start with MBD4+TP53 mutations
2. Run S/P/E efficacy → get pathway scores
3. Convert pathway scores → mechanism vector
4. Run trial matching with mechanism vector
5. Run resistance prediction with SAE features
6. Verify DDR trials rank higher (mechanism fit >0.80)
7. Verify resistance signals detect DDR pathway changes
8. Generate unified validation report

**Output**:
```python
{
    "integration_test": "MBD4+TP53 HGSOC",
    "mutations": ["MBD4 R361*", "TP53 R175H"],
    "mechanism_vector": {"DDR": 0.88, "MAPK": 0.12, ...},
    "trial_matching": {
        "top_5_trials": [...],
        "ddr_trials_in_top_3": True,
        "average_mechanism_fit": 0.85
    },
    "resistance_prediction": {
        "risk_level": "HIGH",
        "ddr_pathway_change": -0.23,
        "signals_detected": 2
    },
    "integration_success": True
}
```

---

## Success Criteria (Adjusted MVP Thresholds - Per Strategic Review)

### Trial Matching:
| Metric | MVP Target | Stretch Target | Measurement |
|--------|------------|----------------|-------------|
| Mechanism fit in response | ✅ | ✅ | Each trial has `mechanism_fit_score` |
| DDR trials rank higher | ✅ | ✅ | Top trials have DDR alignment >0.70 |
| Top-3 accuracy | **≥0.70** | ≥0.80 | Mechanism-matched trials in top-3 |
| MRR | **≥0.65** | ≥0.75 | Mean reciprocal rank |
| Response time | <10s | <10s | For <100 trials |
| MoA coverage report | ✅ | ✅ | Included in validation output |

### Resistance Prediction:
| Metric | MVP Target | Stretch Target | Measurement |
|--------|------------|----------------|-------------|
| Mechanism breakdown present | ✅ | ✅ | DNA repair + pathway escape signals |
| Signal detection accuracy | ≥0.70 | ≥0.75 | True positive rate |
| Risk stratification AUROC | **≥0.65** | ≥0.70 | High risk prediction |
| Mechanism alignment shown | ✅ | ✅ | Pathway-level breakdown |
| Baseline handling documented | ✅ | ✅ | Source + penalty shown |

**Rationale**: This is mechanism alignment, not outcome prediction. Lower thresholds are acceptable for MVP.

---

## What We're NOT Doing (Scope Exclusion)

❌ **NOT outcome prediction** - No PFS/OS correlation (deferred to Zo2)  
❌ **NOT efficacy validation** - Our benchmark shows r=0.037 with PFS  
❌ **NOT TRUE SAE requirement** - Proxy SAE works, TRUE SAE is enhancement  
❌ **NOT drug efficacy modulation** - SAE doesn't change drug confidence  

**We focus ONLY on mechanism-based trial matching and resistance prediction.**

---

## Key Files Reference

| File | Purpose | Status |
|------|---------|--------|
| `mechanism_fit_ranker.py` | Core ranking logic | ✅ Ready (275 lines) |
| `resistance_prophet_service.py` | Resistance prediction | 🔄 Needs mechanism enhancement (689 lines) |
| `autonomous_trial_agent.py` | Query generation | ✅ Already enhanced |
| `trials_agent.py` | Trial search endpoint | 🔄 Needs mechanism fit wiring |
| `pathway_to_mechanism_vector.py` | Pathway→vector conversion | ✅ Ready |
| `sae_feature_service.py` | SAE features | ✅ Ready |

---

## The Deliverables

### For MBD4+TP53 Patient - Trial Matching:

```json
{
  "query": "MBD4+TP53 ovarian cancer trial matching",
  "patient_mechanism_vector": {
    "DDR": 0.88,
    "MAPK": 0.12,
    "PI3K": 0.15,
    "VEGF": 0.10,
    "HER2": 0.05,
    "IO": 0.0,
    "Efflux": 0.0
  },
  "mechanism_fit_applied": true,
  "moa_coverage": {
    "trials_with_moa": 45,
    "trials_without_moa": 2
  },
  "top_trials": [
    {
      "nct_id": "NCT05678901",
      "title": "PARP + ATR Inhibitor in DDR-Deficient Ovarian Cancer",
      "mechanism_fit_score": 0.92,
      "combined_score": 0.87,
      "mechanism_alignment": {
        "DDR": 0.95,
        "MAPK": 0.10,
        "PI3K": 0.05
      },
      "why_matched": "High DDR pathway alignment (0.95) with PARP+ATR mechanism"
    }
  ]
}
```

### For MBD4+TP53 Patient - Resistance Prediction:

```json
{
  "resistance_prediction": {
    "risk_level": "HIGH",
    "probability": 0.72,
    "signal_count": 2,
    "signals": [
      {
        "type": "DNA_REPAIR_RESTORATION",
        "detected": true,
        "mechanism_breakdown": {
          "ddr_pathway_change": -0.23,
          "hrr_essentiality_change": -0.15,
          "exon_disruption_change": 0.05
        },
        "pathway_contributions": {
          "ddr": 0.60,
          "hrr": 0.20,
          "exon": 0.20
        },
        "baseline_source": "patient_baseline",
        "baseline_penalty_applied": false
      },
      {
        "type": "PATHWAY_ESCAPE",
        "detected": true,
        "escaped_pathways": ["DDR"],
        "mechanism_alignment": {
          "DDR": -0.15,
          "MAPK": 0.05
        }
      }
    ],
    "interpretation": "DNA repair capacity dropped (DDR pathway -0.23, HRR -0.15) and DDR pathway escape detected. This suggests treatment failure likely within 3-6 months. This is mechanism alignment, not outcome prediction."
  }
}
```

---

## Bottom Line

**Mission**: Mechanism-based trial matching + mechanism-based resistance prediction  
**Current State**: 80% complete - All services exist (VERIFIED)  
**Remaining Work**: Wire mechanism fit (3-4h) + Enhance resistance (3-4h) + Validation (6-8h)  
**Timeline**: 3-4 days (12-16 hours total, includes 30% buffer)  
**Success**: DDR-targeting trials rank higher, resistance signals show mechanism alignment  
**Independence**: This work is INDEPENDENT of Zo2's outcome calibration

**Status**: ✅ **APPROVED - PROCEED WITH DAY 1 IMPLEMENTATION**

---

## Strategic Review Acknowledgment

This plan was reviewed by Manager Agent on January 28, 2025. The following clarifications were incorporated:

1. ✅ MoA coverage report requirement added
2. ✅ Baseline handling documentation added
3. ✅ Integration test with AYESHA analysis added
4. ✅ Success thresholds adjusted to MVP levels (0.70/0.65/0.65)
5. ✅ Timeline adjusted with 30% buffer (12-16h)
6. ✅ Independence from Zo2 clarified

**Ready to ship!**
