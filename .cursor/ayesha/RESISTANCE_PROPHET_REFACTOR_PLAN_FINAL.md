---
title: Resistance Prophet Refactor Plan (Mars Protocol)
author: Zo (Antigravity Agent)
date: 2026-01-28
status: FINAL_PLAN
iterations: 20
---

# 🧬 Resistance Prophet Refactor: The Mars Protocol

## 🔬 The Diagnosis (20-Point Pathology)
The `ResistanceProphetService.py` (1782 lines) is a "Fragile Monolith" containing critical operational risks. 
Our Deep Audit (iterations 1-20) identified:

### 🚨 Critical Blockers (Immediate Failures)
1.  **Signal/Gene Mismatch**: Prophet outputs 'SIGNALS' (e.g., `DNA_REPAIR_RESTORATION`), but Playbook expects 'GENES' (`MBD4`). **Result: Zero actionable intelligence.**
2.  **Inverted Logic**: `_detect_pathway_escape` only checks for *target loss*, ignoring *bypass activation* (the primary resistance mechanism in OV).
3.  **Math Flaws**: Uses naive multiplication for probability updates (`prob * 1.4`), necessitating artificial `.95` caps.
4.  **Baseline Artifacts**: Defaults to `0.50` for all priors, creating inverted deltas for naturally high/low markers.

### ⚠️ Architectural Debt (Future Failures)
5.  **Disease Leakage**: MM-specific code (cytogenetics) is hardcoded into the generic engine.
6.  **Singleton State**: Global instance holds references to services that may stale.
7.  **Manual Serialization**: Brittle dict construction drops fields silently.
8.  **Error Swallowing**: Playbook failures are logged but hidden from the API response.
9.  **Open Loop**: No mechanism to log predictions for ground-truthing (AUROC improvement impossible).

---

## 🛠️ The Solution: "Organ Extraction" Architecture

We will decompose the service into 5 disconnected, testable "Organs".

### 1. The Core Infrastructure
*   `api/services/resistance_prophet/constants.py`: Thresholds (`0.15`), Risk Levels, and static dicts.
*   `api/services/resistance_prophet/schemas.py`: `ResistanceRiskLevel`, `MechanismBreakdown`, `ResistanceSignalData`, `ResistancePrediction`.

### 2. The Logic Organs
*   `api/services/resistance_prophet/baseline.py`: `get_population_baseline` and future baseline logic.
*   `api/services/resistance_prophet/aggregation.py`: `compute_probability` (with Bayesian Odds + Reliability Weights) and `stratify_risk`.
*   `api/services/resistance_prophet/actions.py`: `determine_actions` and rationale builders.

### 3. The Signal Detectors (Disease-Specific)
*   **Key Design**: Playbook lookups keyed by `ResistanceSignal` Enum (e.g., `DNA_REPAIR_RESTORATION`), not fake genes.
*   `api/services/resistance_prophet/signals/ovarian.py`:
    *   `detect_restoration`: Corrected logic (Positive Delta).
    *   `detect_escape`: Focused on **Target Loss** (Drop). *Bypass Logic differed to Phase 2*.
*   `api/services/resistance_prophet/signals/mm.py`: High-risk genes, Cytogenetics, Drug-class mutations.

### 4. Missing Data & Gate Policy
*Location*: `api/services/resistance_prophet/aggregation.py`
*   **Reliability Weights**: `sig.weight = sig.confidence * baseline_reliability` (0.2 for population baseline).
*   **Treatment Naive Gate**: `has_prior_therapy(history)` check.
    *   If `False` (Truly Naive): Return `NOT_APPLICABLE`.
    *   Ayesha (Post-Chemo Maintenance) is **NOT** Naive -> Prediction Runs.

---

## 📅 Execution Plan (Safe Extraction Sequence)

### Phase 1: Safe Structure (Day 1)
1.  **Extract Schemas**: Create `schemas.py` (pure move).
2.  **Extract Constants**: Create `constants.py` (add `ResistanceSignal` keys to Playbook dicts).
3.  **Contract Test**: Create `tests/test_resistance_contract.py`.
4.  **Extract MM**: Move isolated MM block to `signals/mm.py`.

### Phase 2: Logic Sanitation (Day 1)
5.  **Extract Logic**: Move Baseline, Aggregation, and Actions.
    *   *Inject Fix*: Implement Reliability Weighting & Odds Math.
    *   *Inject Gate*: Implement `has_prior_therapy` check.
6.  **Extract Ovarian**: Move `detect_restoration` to `signals/ovarian.py`.
    *   *Inject Fix*: Correct `positive delta` logic.
    *   *Scope Check*: Implement Target Loss only, no Bypass yet.

### Phase 3: Integration (Day 2)
7.  **Service Wrapper**: Rebuild `ResistanceProphetService`.
8.  **Verify**: Run Ayesha's case.

## 🏁 Success Criteria
1.  **Audit Trail**: A doctor can read `knowledge.json` to verify rules.
2.  **Actionability**: Ayesha's dashboard shows *specific* drugs (e.g., "Switch to IO") instead of generic text.
3.  **Math Integrity**: No probabilities > 1.0, no arbitrary caps.
