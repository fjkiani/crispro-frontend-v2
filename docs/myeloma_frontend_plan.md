# Myeloma Digital Twin – Frontend Implementation Plan

## Goals
- Make live Evo2 scoring obvious (mode, upstream, timing) and resilient to long runs
- Provide clean inputs with validation and presets
- Display interpretable outputs (zeta, impact, pathways) with clear next actions

## Phase 1 – UX foundation (this sprint)
- Components to add
  - `LiveJobBanner`: shows mode (live), upstream service, elapsed timer, tips about cold start
  - `VariantInputList`: dynamic mutation rows (gene, hgvs_p, variant_info, build) with inline validation and presets
- Integrations
  - Update `MyelomaDigitalTwin.jsx` to render `LiveJobBanner` above current ToolRunner
  - Keep API calls via ToolRunner; no backend changes required
- Result view
  - Enhance `MyelomaResponseDisplay` to surface `mode`, `upstream_service`, `zeta_score` and errors

## Phase 2 – Interpretability
- Components
  - `DeltaSummaryCard`: final verdict, min delta, window used (once backend returns it)
  - `PathwayScoresCard`: compact gauges for RAS/MAPK & TP53
  - `VariantDetailCard`: per-variant zeta, impact, ClinVar link, CTA to CRISPR designer
- Hooks
  - `useApiClient`: centralized fetch with 10‑min timeout and abort support
  - `useResultCache`: memoize results keyed by variant signature

## Phase 3 – Advanced controls (when backend supports it)
- `ScoringSettings`: window sizes (1k/2k/4k/8k), strict centering, transcript/exon mode
- `DeltaProfileChart`: ±100 bp local delta profile sparkline

## UX rules
- Inputs must validate `variant_info` (chr:pos REF>ALT) and trim whitespace
- Show a “Warming model (cold start)” banner after 30s
- If |zeta| < 0.5: display “Near‑neutral in current genomic context” helper text
- On error from backend: show `UpstreamErrorAlert` with exact message; never fabricate outputs

## Files to create/update
- Add
  - `oncology-coPilot/oncology-frontend/src/components/myeloma/LiveJobBanner.jsx`
  - `oncology-coPilot/oncology-frontend/src/components/myeloma/VariantInputList.jsx`
- Update
  - `oncology-coPilot/oncology-frontend/src/pages/MyelomaDigitalTwin.jsx`
  - `oncology-coPilot/oncology-frontend/src/components/myeloma/MyelomaResponseDisplay.jsx`

## Definition of done (Phase 1)
- Page shows Live banner with API base and elapsed time
- Default inputs are valid and analyzable in one click
- Results show mode, upstream service, zeta per variant, and any Evo2 errors 

## Progress & Breakthroughs (current status)

- What we shipped
  - Live-only Evo2 scoring end-to-end (no mocks), 10‑minute timeout, redirect follow
  - Strict reference validation (hard 400 on REF mismatch)
  - Multi-window scoring (1k/2k/4k/8k) → returns min_delta and window_used
  - Exon-tight scoring (±600 bp) → returns exon_delta
  - Frontend UX foundation: LiveJobBanner + VariantInputList; results now surface mode, upstream, zeta, min_delta, window_used, exon_delta

- What we uncovered
  - Context dilution is real: 8kb windows can hide signal (KRAS G12D ~ −0.0084)
  - Multi-window and exon view increase magnitude (e.g., min_delta ~ −0.0521 @ 1kb; exon_delta ~ −0.0873) but some sites remain near‑neutral in Evo2
  - Coordinate hygiene is critical: BRAF V600E differs between hg19 vs hg38 (140453136 vs 140753336); strict REF checks prevented silent errors
  - Cold starts: first call ~1–2 minutes; warm calls ~5–10s

- Breakthroughs
  - We now return only verifiable, model‑grounded results with transparent error messages (no fabricated fallbacks)
  - We can quantify “context sensitivity” via min_delta/window_used and exon_delta, making interpretation more trustworthy

- Where we’re going next
  - Add model selector (1B/7B/40B) with backend routing by model_id to tune latency/cost
  - Add local ±100 bp delta profile (sparkline) and 3‑alt sensitivity probe to expose sharp local peaks
  - Add interpretability badges and confidence tags (window consistency, exon corroboration)
  - Expand UI: DeltaSummaryCard, PathwayScoresCard, VariantDetailCard (Phase 2)

- Operating guidance (for users)
  - Always enter exact hg build and coordinates (REF>ALT); the app now fails fast on mismatches
  - If |zeta|, |min_delta| and |exon_delta| are all < 0.5 → likely neutral in Evo2’s view; consider orthogonal models (protein/splice) for confidence
  - For known positives (e.g., BRAF V600E), use the correct build mapping to avoid REF errors

- KPIs we’re tracking
  - Time to first live result (cold) and subsequent warm latency
  - Fraction of inputs hitting strict REF errors (data quality)
  - Distribution of |min_delta| across tested variants; shift after exon and profile features 