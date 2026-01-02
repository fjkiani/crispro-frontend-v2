# Sporadic Cancer Page — Audit Summary (Main Branch)

**Date**: Jan26  
**Branch**: `main`  

## What was true vs what was wrong

- **P0 wiring** (route + payload + provenance): already present in the codebase via `HypothesisValidator` + `useSporadic().getEfficacyPayload()`.
- **The page itself** (`SporadicCancerPage.jsx`) was temporarily broken by bad edits (unclosed JSX / malformed template literal). This has now been corrected in `main`.

## Current state (should be true in `main` now)

- **WIWFM CTA**: `/sporadic-cancer` CTA routes to `/validate`.
- **WIWFM implementation**: `/validate` renders `HypothesisValidator` which calls `POST /api/efficacy/predict`.
- **Tumor context injection**: `HypothesisValidator` uses `getEfficacyPayload()` to include `germline_status` + `tumor_context`.
- **Provenance display**: `SporadicProvenanceCard` renders per-drug when `sporadic_gates_provenance` is present.
- **UX additions on Sporadic page**:
  - Gates preview chips (PARP / IO / confidence cap)
  - “Copy Clinician Summary” (clipboard markdown)

## Key files

- `oncology-coPilocology-frontend/src/pages/SporadicCancerPage.jsx`
- `oncology-coPilot/oncology-frontend/src/pages/HypothesisValidator.jsx`
- `oncology-coPilot/oncology-frontend/src/context/SporadicContext.jsx`
- `oncology-coPilot/oncology-frontend/src/components/sporadic/SporadicProvenanceCard.jsx`

