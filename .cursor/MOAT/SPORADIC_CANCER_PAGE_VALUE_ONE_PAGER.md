# Sporadic Cancer Page — Value (One Source of Truth)

## What the page provides (today)

- **Sporadic workflow UI** for the majority of patients (**germline-negative or unknown**):
  - **Quick Intake (Level 0/1)** and **NGS Upload (Level 2)** via `SporadicWorkflow`.
  - Produces a **TumorContext** (TMB, MSI status, HRD proxy score, completeness score) stored in `SporadicContext`.

- **WIWFM therapy-fit handoff** (demo-ready):
  - CTA routes to `/validate` which renders **`HypothesisValidator`**.
  - `HypothesisValidator` calls `POST /api/efficacy/predict` and injects:
    - `germline_status`
    - `tumor_context`
  - Results include **per-drug sporadic scoring provenance** (`sporadic_gates_provenance`) via `SporadicProvenanceCard`.

- **Explainable scoring adjustments (RUO)**:
  - **PARP gate**: penal vs rescue based on germline status + HRD threshold.
  - **IO gate**: boost based on TMB/MSI.
  - **Confidence caps**: confidence clamping based on completeness tier (L0/L1/L2).

- **Transparency UX on the sporadic page**:
  - **Gates preview**: shows which gates will apply given current TumorContext.
  - **Copy clinician summary**: copies a markdown summary (context + biomarkers + gates) to clipboard.

---

## What the platform provides beyond the page (already built)

- **Clinical trial matching infrastructure** (backend + FE components exist):
  - Germline-required trial exclusion + biomarker boosts (TMB/MSI/HRD) via hybrid search.
  - Operational note: trial indexing may be **partial** until full AstraDB seeding is run; treat as RUO demo mode until verified.

---

## Who this is for (core value + differentiator)

### Oncologist
- **Value**: fast triage from “germline-negative” → “tumor-context-first” → explainable therapy-fit results.
- **Differentiator**: deterministic gates + per-drug prova black-box ranking).

### Patient
- **Value**: stepwise experience to translate “germline negative” into “what tumor data matters next”.
- **Differentiator**: explicit confidence tiering (L0/L1/L2) + RUO guardrails.

### Pharma
- **Value**: biomarker-stratified logic + receipts for cohort-dependence (what validates where).
- **Differentiator**: reproducible validation suite + provenance artifacts (reports/figures/sha256).

---

## What is validated vs what is not (claims discipline)

### ✅ Implementation validated (deterministic correctness)
- Gate math + precedence are unit/integration tested.
- End-to-end UI flow exists: tumor context → efficacy predict → provenance display.

### ✅ Outcome benchmarking validated (where appropriate cohorts exist)
- IO biomarkers validated on biomarker-enriched cohorts (UCEC/COADREAD) in `oncology-coPilot/oncology-backend-minimal/biomarker_enriched_cohorts/`.

### ⚠️ Blocked / must be framed as future work
- **HRD proxy / PARP outcome validation**: bloc assay ground truth.

---

## Key files (source of truth)

### Frontend
- `oncology-coPilot/oncology-frontend/src/pages/SporadicCancerPage.jsx`
- `oncology-coPilot/oncology-frontend/src/components/sporadic/SporadicWorkflow.jsx`
- `oncology-coPilot/oncology-frontend/src/context/SporadicContext.jsx`
- `oncology-coPilot/oncology-frontend/src/pages/HypothesisValidator.jsx`
- `oncology-coPilot/oncology-frontend/src/components/sporadic/SporadicProvenanceCard.jsx`

### Backend (gates)
- `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sporadic_gates.py`

### Outcome-labeled validation package
- `oncology-coPilot/oncology-backend-minimal/biomarker_enriched_cohorts/`

---

## Gaps / next upgrades (high value)

- **Persona toggle** (Clinician / Patient / Pharma) to change copy + panels (same data).
- **Trials CTA**: replace “Coming Soon” with a real minimal trial shortlist (or explicit demo-mode stub).
- **Provider report**: generate a downloadable PDF/Markdown with provenance and RUOermark.

