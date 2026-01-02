# Errata — SPORADIC_CANCER_PRODUCTION_PLAN.md (Jan 2, 2026)

## Critical gaps vs “production ready” claim

### 1) WIWFM CTA is miswired
- `oncology-coPilot/oncology-frontend/src/pages/SporadicCancerPage.jsx` routes the main CTA to **`/validate`**.
- `oncology-coPilot/oncology-frontend/src/App.jsx` maps **`/validate` → `HypothesisValidator`** (not efficacy/WIWFM).
- Impact: the core demo story “sporadic tumor context → efficacy” is **not executable** from the sporadic page today.

**Fix plan**: see `.cursor/MOAT/SPORADIC_CANCER_PAGE_DEMO_PLAN.md`.

### 2) The plan references non-existent or incorrect frontend APIs
- The plan references `useSporadicContext`, but the actual hook is `useSporadic()` from `src/context/SporadicContext.jsx`.

### 3) Embedded gate-math examples are inconsistent
- In the plan’s example `validate_sporadic_gates.py`, TMB bool reference `0.60 → 0.78`, but current gate logic is `TMB ≥20` ⇒ **1.35x**, so `0.60 → 0.81`.
- The “TMB + MSI double boost” framing is misleading: implementation uses an if/elif chain, so **TMB ≥20 takes precedence** over MSI-H.

### 4) E2E is “script exists” but not “run + receipted”
- A script can exist and still be broken by routing / payload wiring.
- For production-readiness, we need:
  - E2E run against a live backend
  - captured receipts (request/response JSON + logs)

## What is missing for demo readiness (minimum)
- Fix WIWFM CTA routing and ensure the efficacy caller includes `germline_status` + `tumor_context`.
- Render `sporadic_gates_provenance` in the efficacy results.
- Replace “Clinical Trials Coming Soon” with either a real trial query or explicit demo-mode stub.

