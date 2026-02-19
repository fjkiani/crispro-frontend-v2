# Holistic Score – Master Reference

**Last updated:** 2026-02-12  
**Status:** Single source of truth (consolidated from 7 docs)

---

## 1. Overview

There are **two distinct systems**:

| System | Scope | Status |
|--------|--------|--------|
| **Trial Holistic Score** | Patient–trial pairs | Implemented |
| **Holistic Clinical Benefit (D-P-M-T-S)** | Patient–regimen pairs | Implemented |

---

## 2. Trial Holistic Score (Patient–Trial)

### Formula
```
Holistic Score = (0.5 × Mechanism Fit) + (0.3 × Eligibility) + (0.2 × PGx Safety)
```

Weights in code: `api/services/holistic_score/models.py`  
- MECHANISM_FIT_WEIGHT = 0.5  
- ELIGIBILITY_WEIGHT = 0.3  
- PGX_SAFETY_WEIGHT = 0.2  

### Interpretation thresholds
`api/services/holistic_score/interpreter.py`  
- HIGH: ≥ 0.8  
- MEDIUM: ≥ 0.6  
- LOW: ≥ 0.4  
- VERY_LOW: &lt; 0.4  

### Code anchors
| Role | Path |
|------|------|
| Service | `oncology-coPilot/oncology-backend-minimal/api/services/holistic_score/service.py` |
| Models | `api/services/holistic_score/models.py` |
| Mechanism fit | `api/services/holistic_score/mechanism_fit.py` |
| Eligibility | `api/services/holistic_score/eligibility_scorer.py` |
| PGx safety | `api/services/holistic_score/pgx_safety.py` |
| Interpreter | `api/services/holistic_score/interpreter.py` |
| Trial integration | `api/services/ayesha_care_plan/trial_service.py` → `_add_holistic_scores()` |
| Frontend component | `oncology-coPilot/oncology-frontend/src/components/trials/HolisticScoreCard.jsx` |
| Trial card | `TrialMatchCard.jsx`, `TrialMatchesCard.jsx` |

### Integration
- `trial_service.py` uses `get_holistic_score_service()` and `_add_holistic_scores()`
- Patient profile built via `_build_patient_profile_for_holistic()`
- Mechanism vector via `_compute_mechanism_vector_from_tumor_context()` or defaults

### Tests
- `oncology-coPilot/oncology-backend-minimal/tests/test_holistic_score_integration.py`
- 5/5 tests passing (service import, single trial, batch, trial integration, mechanism vector)

### Response shape
```json
{
  "holistic_score": 0.94,
  "mechanism_fit_score": 0.881,
  "eligibility_score": 1.0,
  "pgx_safety_score": 1.0,
  "holistic_interpretation": "HIGH",
  "holistic_recommendation": "...",
  "holistic_caveats": []
}
```

---

## 3. Holistic Clinical Benefit Score (Patient–Regimen, D-P-M-T-S)

### Purpose
Score patient–regimen pairs for trial enrollment, next-line, or monitoring.

### Components
- **D:** Diagnostic fit (DDR_bin, disease context)
- **P:** Prognostic risk (PFI, PFS, line of therapy)
- **M:** Predictive mechanism fit (7D cosine similarity; reuses `holistic_score.mechanism_fit`)
- **T:** Therapeutic dynamics (KELIM, CA-125 early response)
- **S:** Safety/tolerability (PGx; reuses `holistic_score.pgx_safety`)

### Formula
```
Overall = w_D × D + w_P × P + w_M × M + w_T × T + w_S × S
```

### Use-case weights
`api/services/resistance/config/holistic_clinical_benefit_config.py`

| Use case | D | P | M | T | S |
|----------|---|---|---|---|---|
| trial_enrollment | 0.20 | 0.10 | 0.45 | 0.00 | 0.25 |
| next_line | 0.10 | 0.25 | 0.35 | 0.20 | 0.10 |
| monitoring | 0.05 | 0.20 | 0.15 | 0.45 | 0.15 |

For trial_enrollment, T is set to 0 and weights are renormalized (D, P, M, S only).

### Code anchors
| Role | Path |
|------|------|
| Service | `api/services/holistic_clinical_benefit/service.py` |
| D | `api/services/holistic_clinical_benefit/diagnostic_fit.py` |
| P | `api/services/holistic_clinical_benefit/prognostic_risk.py` |
| T | `api/services/holistic_clinical_benefit/therapeutic_dynamics.py` |
| Regimen MoA | `api/services/holistic_clinical_benefit/regimen_moa.py` |
| Config | `api/services/resistance/config/holistic_clinical_benefit_config.py` |
| API | `api/routers/resistance.py` → POST `/api/resistance/holistic-clinical-benefit` |
| Frontend | `HolisticClinicalBenefitCard.jsx`, `TrialHolisticScoreCard.jsx` |
| Hook | `useHolisticClinicalBenefit.js` |

### Tests
- `api/services/holistic_clinical_benefit/test_service.py`
- `api/services/holistic_clinical_benefit/test_integration.py`
- `api/services/csi/test_integration.py` (uses holistic_clinical_benefit)

---

## 4. Dependencies

### Trial holistic score
- Mechanism fit: 7D vectors (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
- Eligibility: recruiting, disease, age, location, biomarkers
- PGx: DPYD, TPMT, UGT1A1, CYP2D6, etc.

### Holistic clinical benefit
- Reuses `holistic_score.mechanism_fit` (M)
- Reuses `holistic_score.pgx_safety` (S)
- D: DDR_bin
- P: Timing engine (PFI, PFS, line of therapy)
- T: Kinetic engine (KELIM, CA-125)

---

## 5. Validation notes

- **TOPACIO AUROC 0.714:** Mentioned in prior design docs as trial matching validation. No validation script found in repo; treat as design claim.
- **Trial holistic score:** Formula and thresholds validated by `test_holistic_score_integration.py` (5/5 passing).

---

## 6. Canonical flows

### Trial matching with holistic score
1. Call trial matching (e.g. `/api/ayesha/complete_care_v2`)
2. `trial_service._add_holistic_scores()` computes scores per trial
3. Response includes `holistic_score`, `mechanism_fit_score`, `eligibility_score`, `pgx_safety_score`, `holistic_interpretation`
4. `HolisticScoreCard.jsx` renders breakdown and interpretation

### Holistic clinical benefit
1. POST `/api/resistance/holistic-clinical-benefit` with `patient_id`, `regimen_id`, `use_case`
2. `compute_holistic_clinical_benefit_score()` orchestrates D, P, M, T, S
3. Returns `D`, `P`, `M`, `T`, `S`, `overall`, `weights`, `component_available`, `breakdown`, `interpretation`, `recommendation`, `provenance`

---

## 7. Run tests

```bash
# Trial holistic score
cd oncology-coPilot/oncology-backend-minimal
python3 -m pytest tests/test_holistic_score_integration.py -v -s

# Holistic clinical benefit
python3 -m pytest api/services/holistic_clinical_benefit/test_service.py api/services/holistic_clinical_benefit/test_integration.py -v -s
```
