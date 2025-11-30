# Final Comprehensive Plan: Mechanism-Based Trial Matching & Resistance Prediction

**Date**: January 28, 2025  
**Status**: ✅ **REALITY CHECKED** - Ready for Implementation  
**Location**: Main Repository (`/Users/fahadkiani/Desktop/development/crispr-assistant-main/`)

---

## Executive Summary

**Mission**: Deliver mechanism-based trial matching AND mechanism-based resistance/sensitivity prediction

**Current State**: 80% complete - All services exist, need wiring and enhancement

**Timeline**: 3 days (9-13 hours total)

---

## What EXISTS ✅

### Core Services (All Implemented)

| Service | File | Status | Notes |
|---------|------|--------|-------|
| **Mechanism Fit Ranker** | `api/services/mechanism_fit_ranker.py` | ✅ Complete | Formula: 0.7×eligibility + 0.3×mechanism_fit |
| **Resistance Prophet** | `api/services/resistance_prophet_service.py` | ✅ Complete | 2-of-3 signals, needs mechanism breakdown |
| **Pathway→Vector** | `api/services/pathway_to_mechanism_vector.py` | ✅ Complete | 7D vector conversion |
| **Autonomous Trial Agent** | `api/services/autonomous_trial_agent.py` | ✅ Complete | 10 query templates, DDR detection |
| **SAE Feature Service** | `api/services/sae_feature_service.py` | ✅ Complete | TRUE SAE flag ready |

### Data Resources

| Resource | File | Status |
|----------|------|--------|
| **MoA Vectors** | `api/resources/trial_moa_vectors.json` | ✅ 47 trials tagged |
| **SAE Mapping** | `api/resources/sae_feature_mapping.json` | ✅ 288KB, 88 features |

### Integration Points

| Integration | File | Status |
|-------------|------|--------|
| **Resistance Prophet** | `api/routers/ayesha_orchestrator_v2.py` | ✅ Integrated (opt-in flag) |
| **Trial Search** | `api/routers/trials_agent.py` | ⚠️ Needs mechanism fit wiring |

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

**File**: `api/services/resistance_prophet_service.py`  
**Work**: 2-3 hours

### Gap 3: Validation Scripts Missing

**Needed**:
- `scripts/validate_mechanism_trial_matching.py`
- `scripts/validate_mechanism_resistance_prediction.py`
- `scripts/validate_mbd4_tp53_mechanism_capabilities.py`

**Work**: 4-6 hours

---

## Implementation Plan (3 Days)

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
    
    # If still no mechanism vector, could compute from mutations via SAE service
    # (Defer to Day 2 if needed)
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
    else:
        results["mechanism_fit_applied"] = False
        if not mechanism_vector:
            results["mechanism_fit_warning"] = "No mechanism vector provided - ranking by eligibility only"
```

**Files to Modify**:
- `api/routers/trials_agent.py`

---

### Day 2: Enhance Resistance Prophet with Mechanism Breakdown (2-3 hours)

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
        }
    }
    
    return ResistanceSignalData(...)
```

#### Task 2.2: Enhance Pathway Escape Signal

**File**: `api/services/resistance_prophet_service.py`

**Current**: Tracks pathway shifts in provenance  
**Enhancement**: Add structured `escaped_pathways` and `mechanism_alignment`

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

### Day 3: Create Validation Scripts (4-6 hours)

#### Task 3.1: Mechanism-Based Trial Matching Validation

**File**: `scripts/validate_mechanism_trial_matching.py` (NEW)

**8-Task Verification**:
1. Trial data quality (47 MoA-tagged trials exist)
2. Mechanism vector structure (7D format correct)
3. Mechanism fit computation (cosine similarity correct)
4. Combined score formula (0.7×eligibility + 0.3×mechanism_fit)
5. Ranking accuracy (Top-3 ≥80%, MRR ≥0.75)
6. Pathway alignment (DDR trials rank higher for DDR-high patients)
7. Edge cases (all-zero vector, missing MoA vectors)
8. Consistency (deterministic results)

**Metrics to Report**:
```python
{
    "top3_accuracy": 0.85,  # Target: ≥0.80
    "mrr": 0.78,  # Target: ≥0.75
    "mechanism_fit_scores": {
        "mean": 0.72,
        "std": 0.15,
        "min": 0.50,
        "max": 0.95
    },
    "pathway_alignment": {
        "DDR": 0.88,  # DDR trials rank higher for DDR-high patients
        "MAPK": 0.75,
        "PI3K": 0.68
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
6. Baseline handling (population average when missing)
7. Confidence modulation (MEDIUM cap when CA-125 missing)
8. Consistency (deterministic)

**Metrics to Report**:
```python
{
    "signal_detection_accuracy": {
        "dna_repair_restoration": 0.82,
        "pathway_escape": 0.75,
        "ca125_kinetics": 0.68
    },
    "risk_stratification": {
        "high_risk_auroc": 0.72,  # Target: ≥0.70
        "sensitivity": 0.78,  # Target: ≥0.75
        "specificity": 0.71  # Target: ≥0.70
    },
    "mechanism_breakdown_accuracy": {
        "ddr_contribution": 0.60,  # Verified: matches formula
        "hrr_contribution": 0.20,
        "exon_contribution": 0.20
    }
}
```

#### Task 3.3: MBD4+TP53 End-to-End Validation

**File**: `scripts/validate_mbd4_tp53_mechanism_capabilities.py` (NEW)

**Test Both Capabilities Together**:
1. Run trial matching with MBD4+TP53 mechanism vector
2. Run resistance prediction with MBD4+TP53 SAE features
3. Verify DDR trials rank higher (mechanism fit >0.80)
4. Verify resistance signals detect DDR pathway changes
5. Generate unified validation report

---

## Success Criteria

### Trial Matching:
| Metric | Target | Measurement |
|--------|--------|-------------|
| Mechanism fit in response | ✅ | Each trial has `mechanism_fit_score` |
| DDR trials rank higher | ✅ | Top trials have DDR alignment >0.70 |
| Top-3 accuracy | ≥0.80 | Mechanism-matched trials in top-3 |
| MRR | ≥0.75 | Mean reciprocal rank |
| Response time | <10s | For <100 trials |

### Resistance Prediction:
| Metric | Target | Measurement |
|--------|--------|-------------|
| Mechanism breakdown present | ✅ | DNA repair + pathway escape signals |
| Signal detection accuracy | ≥0.75 | True positive rate |
| Risk stratification AUROC | ≥0.70 | High risk prediction |
| Mechanism alignment shown | ✅ | Pathway-level breakdown |

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
| `mechanism_fit_ranker.py` | Core ranking logic | ✅ Ready |
| `resistance_prophet_service.py` | Resistance prediction | 🔄 Needs mechanism enhancement |
| `autonomous_trial_agent.py` | Query generation | ✅ Already enhanced |
| `trials_agent.py` | Trial search endpoint | 🔄 Needs mechanism fit |
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
        }
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
**Current State**: 80% complete - All services exist  
**Remaining Work**: Wire mechanism fit (3-4h) + Enhance resistance (2-3h) + Validation (4-6h)  
**Timeline**: 3 days (9-13 hours total)  
**Success**: DDR-targeting trials rank higher, resistance signals show mechanism alignment  

**Ready to proceed with implementation!**





