# Sporadic Cancer Page — Demo Readiness Plan (3 Personas)

## Key discovery (current gap)
- The **"Run Efficacy Prediction (WIWFM)"** CTA on `oncology-coPilot/oncology-frontend/src/pages/SporadicCancerPage.jsx` routes to **`/validate`**, which is currently mapped to **`HypothesisValidator`** (not WIWFM/efficacy) in `oncology-coPilot/oncology-frontend/src/App.jsx`.
- Result: the sporadic workflow is **not demo-correct** for the intended “tumor-context → efficacy” story.

## Product frame: what this page should do
- Turn **germline-negative or unknown** patients into **tumor-context-first** analysis.
- Make the next action obvious and safe: **"Run therapy fit (WIWFM)"**, **"Trial match"**, **"Next tests"**, with **transparent provenance** and **RUO posture**.

---

## Persona 1: Oncologist (Clinical decision support; RUO posture)

### Core value
- **Fast, structured triage*poradic case into actionable next steps:
  - *What is missing?* (L0/L1/L2 completeness)
  - *What should be ordered next?* (ctDNA/HRD/MSI/TMB where applicable)
  - *What therapy classes are plausible right now?* (IO gate visibility; PARP penalty/rescue logic)
  - *What should be done this week?* (trial shortlist, monitoring plan)

### Differentiator
- **Deterministic gates + provenance** (no black box):
  - PARP penalty / rescue is explainable.
  - IO boosts are explainable.
  - Confidence caps are explainable via completeness tiers.
  - Validation receipts exist (UCEC/COADREAD + missingness simulation + scripts).

### Demo-ready build (what we ship)
- Replace current WIWFM CTA with a real efficacy flow (see Build Tasks).
- Add a **Clinician Summary card**:
  - “Data Level: L0/L1/L2”
  - “Key biomarkers: TMB, MSI, HRD proxy (if present)”
  - “Gates that will apply: PARP, IO, confidence cap (with reasons)”
  - “Next tests recommended (RUO)”
- Add **Copy-to-clipboard provider note** (Markdown)+ top 5 drug classes + top 5 trials.

---

## Persona 2: Patient (Actionable clarity; plain language)

### Core value
- Make the workflow understandable and reduce anxiety:
  - “You’re germline-negative — this is common. Next we look at tumor signals.”
  - “Here’s what we can say now vs what requires NGS.”
  - “Here are the next concrete steps you can take this week.”

### Differentiator
- **Transparent guardrails**:
  - Confidence is not a vibe; it’s tied to data completeness.
  - “Recommendation” vs “Prediction” language enforced.

### Demo-ready build (what we ship)
- Add a **Patient mode** toggle:
  - Hide technical jargon (HRD proxy, etc.) behind tooltips.
  - Show a “What to ask your oncologist” checklist.
  - Show “Why NGS matters” + timeline panel.

---

## Persona 3: Pharma (Trial enrichment & biomarker strategy)

### Core value
- Show how gates support **trial enrichment / segmentation** decisions:
  - Sensitivity to TMB cutoffs (UCEC sweep).
  - Baseline stratevs OR gate).
  - Cohort dependence receipts (COADREAD transparency).
  - Missingness stress-test for operational realism.

### Differentiator
- **Receipts + reproducibility**:
  - One-command validation suite.
  - SHA256 provenance on cohort JSONs.
  - Figures inventory + claims map.

### Demo-ready build (what we ship)
- Add a **Pharma mode** toggle:
  - Surface links to: `oncology-coPilot/oncology-backend-minimal/biomarker_enriched_cohorts/docs/*`
  - Show “Validation artifacts” panel: figures + reports + cohort paths.

---

## Build tasks (to make page demo-ready)

### P0 — Fix the WIWFM action (must-do)
- Update the CTA on `src/pages/SporadicCancerPage.jsx` so it routes to a real efficacy experience.
  - Option A (fastest): route to `MutationExplorer` and ensure it consumes `useSporadic().getEfficacyPayload()` when calling `/api/efficacy/predict`.
  - Option B: create a dedicated route `/wiwfm` (new page) that only runs efficacy with the current tumor context.

### P0 — Ensure tumor context is aused in efficacy calls
- Confirm the frontend code that calls `/api/efficacy/predict` merges:
  - `germline_status: useSporadic().germlineStatus`
  - `tumor_context: useSporadic().tumorContext`

### P0 — Demo-level Trial CTA
- Replace “Coming Soon” with a live call:
  - Use existing trials endpoint (or minimal stub) to display a list of 5 trials.
  - If endpoint isn’t available, show “Demo mode” with clear provenance.

### P1 — Persona mode toggle
- Add a segmented control (Clinician / Patient / Pharma) that toggles:
  - copy, visible panels, and action CTAs (same data)

### P1 — Provider export
- “Copy summary” button that produces clinician-facing markdown.

