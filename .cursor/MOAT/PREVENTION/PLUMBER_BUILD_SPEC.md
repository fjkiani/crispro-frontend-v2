# 🔧 Plumber Build Spec (Stage-3 Prevention MVP)

**Date:** December 24, 2025  
**Audience:** Plumber (engineering), Zo (spec owner)  
**Goal:** Ship a demo-credible “Prevention MVP” that turns DDR_bin trend → alert → tiered actions.

---

## 0) Definition of Done (MVP)

### **MVP Output (one patient)**
Given a patient with:
- baseline DDR_bin
- ≥1 follow-up DDR_bin value

System must produce:
1. **Alert level**: NONE / WARNING / ALERT / CRITICAL
2. **Trend**: STABLE / DECLINING / IMPROVING
3. **Tiered actions** (Tier 1/2/3) with evidence levels
4. **Next tests** (targeted panel suggestion)
5. **Provenance** (model/version, inputs used)

### **MVP Acceptance tests**
- Deterministic (same input → same output)
- Unit tests cover threshold logic + edge cases
- API returns valid JSON contract
- UI can render output (even minimal)

---

## 1) Data Contracts

### **1.1 Patient Prevention State (add to patient state)**

```json
{
  "prevention": {
    "baseline": {
      "date": "2025-01-15",
      "ddr_bin": 0.88,
      "source": "tumor|ctdna",
      "variant_count": 52,
      "model_version": "true_sae_diamonds.v1"
    },
    "history": [
      {
        "date": "2025-04-15",
        "ddr_bin": 0.87,
        "source": "ctdna",
        "variant_count": 10,
        "model_version": "true_sae_diamonds.v1"
      }
    ]
  }
}
```

### **1.2 Prevention Output Contract**

```json
{
  "prevention": {
    "alert_level": "NONE|WARNING|ALERT|CRITICAL",
    "trend": "STABLE|DECLINING|IMPROVING|INSUFFICIENT_DATA",
    "delta_from_baseline": -0.06,
    "confidence": "HIGH|MEDIUM|LOW",
    "interpretation": "...",
    "recommended_actions": [
      {
        "tier": 1,
        "action": "...",
        "evidence_level": "A|B|C",
        "safety_notes": ["..."],
        "when": "..."
      }
    ],
    "next_tests": [
      {
        "test": "ctDNA targeted panel (RAD51C/BRCA1/PALB2)",
        "rationale": "Confirm HR restoration mechanism",
        "turnaround_days": 5,
        "cost_usd": 500
      }
    ],
    "provenance": {
      "engine": "prevention_v1_rule_based",
      "threshold_policy": "manager_decisions_v1",
      "inputs_used": {
        "baseline_date": "...",
        "latest_date": "...",
        "variant_count_latest": 10
      }
    }
  }
}
```

---

## 2) Services to Implement

### **2.1 `prevention_alert_service.py` (pure logic)**
Responsibilities:
- compute delta from baseline
- compute trend (latest vs previous)
- map to alert level

Edge cases:
- missing baseline → INSF DATA
- low variant count → LOW confidence / INSUFFICIENT_DATA
- only one datapoint → trend = INSUFFICIENT_DATA

### **2.2 `prevention_recommendation_service.py` (tiered actions)**
Inputs:
- alert_level, bins (DDR now; later MAPK/PI3K), current therapy

Outputs:
- Tier 1/2/3 actions
- safety notes
- next tests

**Important:** No medical “orders”; output phrased as “discuss with oncologist” unless Manager authorizes stronger language.

---

## 3) Endpoints

Pick ONE (Manager decision):

### **Option A: Standalone**
- `POST /api/prevention/evaluate`

### **Option B: Inside Complete Care**
- `POST /api/complete_care/v2?include_prevention=true`

---

## 4) Fixtures

Create 3 deterministic fixtures:

1. **Stable**: baseline 0.88, latest 0.87 → alert NONE
2. **Warning**: baseline 0.88, latest 0.82 → WARNING
3. **Critical**: baseline 0.88, latest 0.70 → ALERT/CRITICAL (per Manager thresholds)

Store in:
- `oncology-coPilot/oncology-backend-minimal/data/patient_states/` (or equivalent)

---

## 5) Tests

### **Unit tests (must-have)**
- threshold mapping
- variant_count guardrails
- trend computation

### **Contract tests (must-have)**
- output JSON schema validates

### **Validation runner integration**
Add a new script:
- `scripts/validation/validate_prevention_mvp.py`

And wire it into:
- `scripts/validation/run_validations.py`

---

## 6) “Hard Experiments” (turn unknowns into numbers)

These are the experiments that convert Manager uncertainty into locked thresholds:

1. **Stability vs variant count**
- subsample Tier-3 variants to n=5/10/20
- measure DDR_bin variance

2. **Separation check (DDR-high should rank DDR trials)**
- reuse mechanism-fit validator logic for trial ranking sanity

3. **Reversion sanity**
- if any RAD51C/BRCA reversion examples exist, check DDR_bin directionality

---

## 7) Known Limitations (must be explicit)

- No longitudinal real-world validation yet (unless we secure cohort)
- No true causal guarantee of interventions
- Must display RUO + provenance

