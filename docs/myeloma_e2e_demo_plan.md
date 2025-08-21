## Myeloma Digital Twin – End-to-End Demo Plan

### Purpose
- Establish a seamless, investor-ready demo that proves end-to-end capability: live Evo2 scoring, interpretable evidence, pathway-level decisioning, and trustworthy provenance.
- Unify conversational guidance, agents, backend orchestration, and a clear UX into one cohesive experience.

### Definition of “Solved” (Demo Acceptance)
- Live-only results (no mocks) with model switching (1B/7B/40B), batch inputs, warmup, and robust timeouts.
- For each variant: zeta_score, min_delta + window_used, exon_delta, local ±100bp delta profile, 3-alt sensitivity probe, and a confidence score with reasons.
- Pathway-level decision: RAS/MAPK and TP53 aggregation → final “Likely Resistant/Sensitive”, with thresholds documented and versioned.
- Evidence and provenance panel: exact inputs, Ensembl URLs, fetched REF base, model/version/endpoint, request/response JSON download, and run signature.
- Reproducibility: one-click “Reproduce run” re-executes with pinned model and identical inputs, producing matching outputs.

## End-to-End Workflow (User Journey)
1) Configure and Warm Up
   - Select Evo2 model (1B/7B/40B), click “Confirm & Warm Up.”
   - Live banner shows upstream URL, elapsed time, and readiness.
2) Input Variants (Batch)
   - Add KRAS/NRAS/BRAF/TP53 variants via presets or manual entry.
   - Inline validation for chr:pos REF>ALT and genome build.
3) Run Analysis
   - Backend orchestrates: /score_variant → /score_variant_multi → /score_variant_exon → /score_variant_profile → /score_variant_probe.
4) Interpret Results
   - For each variant: zeta/min_delta/exon_delta, profile sparkline, probe table, confidence with reasons.
   - Pathway aggregation surfaces final verdict and shows contribution per gene.
5) Trust & Export
   - Evidence panel shows provenance, Ensembl window and REF check, upstream endpoint, and a signed JSON for audit.
   - “Reproduce” re-runs with pinned versions.

## Architecture Overview

### Backend (Orchestrator)
- File: `oncology-coPilot/oncology-backend-minimal/api/index.py`
- Responsibilities:
  - Single public API facade (Vercel or local uvicorn), long timeouts, follow redirects.
  - Model routing: 1B/7B/40B via env `EVO_URL_1B|7B|40B`.
  - Orchestrated calls: `/score_variant`, `/score_variant_multi`, `/score_variant_exon`, `/score_variant_profile`, `/score_variant_probe`.
  - Confidence derivation: effect size (min_delta), window consistency, exon corroboration.
  - Pathway decision: RAS/MAPK + TP53 aggregation with documented thresholds.
  - Provenance: return mode, upstream_service, selected_model, and errors transparently.

### Evo2 Services (Modal)
- File(s): `src/services/evo_service/main.py` (40B + 7B-in-main endpoints)
- Endpoints:
  - `/score_delta` – basic ref vs alt scoring
  - `/score_variant` – Ensembl window fetch + delta
  - `/score_variant_multi` – windows [1k,2k,4k,8k] → min_delta, window_used
  - `/score_variant_exon` – tight window (±600 bp) → exon_delta
  - `/score_variant_profile` – local ±100bp profile
  - `/score_variant_probe` – 3-alt sensitivity probe
- Constraints:
  - Strict REF validation (no silent mismatches)
  - Clear error messages; no fallbacks

### Frontend (Next.js/React or current React app)
- Files (current foundation):
  - Live banner: `oncology-coPilot/oncology-frontend/src/components/myeloma/LiveJobBanner.jsx`
  - Variant input list: `oncology-coPilot/oncology-frontend/src/components/myeloma/VariantInputList.jsx`
  - Results display: `oncology-coPilot/oncology-frontend/src/components/myeloma/MyelomaResponseDisplay.jsx`
  - Page: `oncology-coPilot/oncology-frontend/src/pages/MyelomaDigitalTwin.jsx`
- Components to add:
  - ModelSelector: chooses `evo2_1b|evo2_7b|evo2_40b`
  - DeltaProfileChart: ±100bp sparkline with peak marker
  - SensitivityProbeTable: alt→delta table with most disruptive ALT
  - EvidencePanel: upstream URL, model version, Ensembl links, REF check, JSON download
  - UpstreamErrorAlert: exact upstream exception messages

## Agents (Specialized AI Assistant Model as Orchestrated Agents)
- Conversation Orchestrator Agent
  - Guides user, translates goals to backend requests, maintains session state.
- Variant QA Agent
  - Validates inputs (build, chr, pos, REF>ALT), proposes corrections (e.g., TP53 coordinate), links Ensembl.
- Evidence Synthesizer Agent
  - Summarizes per-variant signals (zeta/min_delta/exon_delta/profile/probe) and external annotations (ClinVar, COSMIC, OncoKB) into a human-readable rationale.
- Threshold Calibrator Agent
  - Monitors cohort/batch distributions and recommends threshold updates; writes config diffs for review.
- Provenance Notary Agent
  - Signs runs, stores artifacts (request/response JSON, hashes, timestamps), and prepares audit packages.
- Benchmark & Regression Agent
  - Runs nightly panel, checks AUROC/calibration, blocks deploy on regression, posts summary to dashboard.

## Backend API Contract (Key Routes)
- POST `/api/predict/myeloma_drug_response`
  - Request: `{ model_id, gene|hgvs_p|variant_info|build OR mutations[] }`
  - Response: `{ prediction, pathway_scores, detailed_analysis[], mode, upstream_service }`
  - detailed_analysis[i].evo2_result includes: `zeta_score, min_delta, window_used, exon_delta, confidence_score, confidence_reason`
- POST `/api/evo/score_variant_profile`
  - Request: `{ assembly, chrom, pos, ref, alt, flank?, radius? }`
  - Response: `{ profile: [{offset, delta}], peak_delta, peak_offset }`
- POST `/api/evo/score_variant_probe`
  - Request: `{ assembly, chrom, pos, ref }`
  - Response: `{ probes: [{alt, delta}], top_alt, top_delta }`
- POST `/api/evo/warmup`
  - Request: `{ model_id }`
  - Response: `{ status, selected_model, upstream_service, elapsed_sec }`

## Evidence & Trust (What We Display)
- Inputs: assembly, chr, pos, REF>ALT (trimmed & validated)
- Ensembl URLs: region fetch link, window sizes, exon flank
- REF check: fetched base at index; mismatch → fail with guidance
- Model & endpoint: model_id, Modal URL, version/hash if available
- Metrics: zeta_score, min_delta/window_used, exon_delta, profile peak, top_alt from probe
- Confidence: numeric [0-1] + reason text (effect size, window consistency, exon corroboration)
- Provenance: downloadable JSON of full inputs/outputs, timestamp, optional signature

## Demo Script (Condensed)
1) Warm up model (show banner → ready)
2) Add default 4 variants (KRAS/NRAS/BRAF/TP53); fix TP53 coordinate if needed via Variant QA agent
3) Run analysis; watch live statuses
4) Show per-variant cards with evidence; open delta profile and probe
5) Show pathway aggregation → Likely Resistant; rationale side panel
6) Download provenance JSON; click Reproduce → matching results

## Validation & Benchmarks
- Initial Panel: ~12–20 myeloma-relevant variants (KRAS/NRAS codon 12/13/61, BRAF V600E, TP53 hotspots with correct hg38)
- Metrics: AUROC, AUPRC, ECE/calibration, cross-model agreement rates, REF mismatch rate
- Process: nightly automated run, dashboard, deploy block on regression

## Operations & Reliability
- Timeouts: 10 min total, 30 s connect, retry on transient 5xx with backoff
- Cold-start: pre-demo warmup; proactive warming before batch
- Logging: info logs per endpoint (start, params summary, outcome), errors include upstream detail
- Caching: memoize results by variant signature; surface cache hits in UI

## Milestones & DoD
- M1: 40B fully wired (done), 7B routed (done), warmup (done), batch (done)
- M2: Profile + Probe wired and visible; EvidencePanel added; ModelSelector in UI
- M3: Benchmark panel live; confidence calibration set; provenance signing enabled
- DoD: Reproducible demo run with all evidence artifacts, passing nightly benchmark, and no schema mismatches

## Dependencies & Config
- Backend: FastAPI, httpx; Vercel deploy (`api/index.py`), ENV: `EVO_URL_1B|7B|40B`
- Modal services: `src/services/evo_service/main.py` (40B + 7B app), GPU H100/A10G
- Frontend: React/MUI components; long-timeout fetch; error surfacing and abort controls

## Future Extensions
- Add 1B service for ultra-low latency
- Integrate AlphaMissense/ESM scores for orthogonal effect validation
- Expand pathway model (PI3K/AKT, NF-κB) and therapy-specific classifiers
- Add cohort-based calibration and individualized treatment recommendations 