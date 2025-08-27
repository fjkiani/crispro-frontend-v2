# Myeloma Digital Twin — Agentic Pipeline Plan and Documentation

## 1) Executive Summary
We deployed a live, interpretable Evo2-powered variant-effect pipeline and a Myeloma Digital Twin UI that converts patient mutations into pathway-level decisions with full provenance. This document consolidates the current state, architecture, evidence protocol, benchmarking approach, analytics plan, and a concrete roadmap to agentic autonomy.

## 2) Achievements to Date (Live)
- Live Evo2 endpoints (7B) proxied by a thin backend with strict REF-allele validation and explicit upstream error surfacing
- Full evidence suite per variant: multi-window zeta, exon-tight corroboration, local delta profile, 3‑alt sensitivity probe
- Confidence heuristic (effect size + window consistency + exon corroboration) and per‑variant call chips (Disruptive/Neutral/Unknown)
- Myeloma Digital Twin UI: model selector, warmup, batch scoring, evidence actions (Run All Profiles/Probes), provenance downloads
- Dual-model comparison (optional): 7B↔40B per-variant agreement rate with per-variant comparisons
- ClinVar-driven benchmark generator (high-confidence reviews) + chunked benchmarking
- Live benchmark (7B, n=160): AUROC 0.9709, AUPRC 0.9614; calibration table
- Supabase logging (optional): per-run and per-variant inserts; ready for dashboarding

## 3) System Architecture
- Vercel thin proxy (`oncology-backend-minimal`):
  - `/api/predict/myeloma_drug_response` (batch scoring with non‑fatal per‑variant errors, optional `dual_compare`)
  - `/api/evo/score_variant*` passthrough (delta, multi, exon, profile, probe, warmup)
  - Provenance in every response (mode, selected_model, upstream_service, run_signature)
  - Async Supabase logging: runs + variants (+ alt_model agree_rate)
- Modal Evo2 services (`src/services/evo_service/main.py`):
  - 7B and 40B endpoints with parity; strict REF base checks and identical-sequence guard
- Frontend (`oncology-frontend`):
  - `MyelomaDigitalTwin.jsx`, `MyelomaResponseDisplay.jsx`, `ModelSelector.jsx`, `VariantInputList.jsx`, `LiveJobBanner.jsx`
  - `ProgressFlow.jsx` revamped: dynamic pill/timeline with icons and custom steps

## 4) Evidence‑First Protocol (VUS Doctrine — Evidence‑or‑Abstain)
- Signals per variant:
  - minΔ from multi-window; window consistency; exon corroboration; local profile peak; 3‑alt probe; (optional) 7B vs 40B corroboration
- Confidence (baseline): `0.5*s1 + 0.3*s2 + 0.2*s3`
- Decision policy:
  - Likely Disruptive if: |minΔ| ≥ P90_pathogenic AND s3 ≥ 0.7 AND exon corroborates AND |peakΔ| large
  - Likely Neutral if: |minΔ| ≤ P50_benign AND |peakΔ| small AND no corroboration
  - Unknown if: evidence weak/contradictory or 7B vs 40B disagree or confidence < 0.4
- Strict hygiene: REF base validation; transparent upstream errors; never fabricate

References: `/.cursor/rules/vus_testing_doctrine.mdc`

## 5) API Contracts (Key)
- POST `/api/predict/myeloma_drug_response`
  - Input: `{ model_id, mutations: [{ gene, hgvs_p, variant_info: "chr:pos REF>ALT", build }], dual_compare? }`
  - Output: `{ prediction, pathway_scores, detailed_analysis: [ per-variant evidence ], provenance..., (dual_compare)? }`
- POST `/api/evo/score_variant`, `/score_variant_multi`, `/score_variant_exon`, `/score_variant_profile`, `/score_variant_probe`, `/evo/warmup`

## 6) Frontend UX (Myeloma Digital Twin)
- Pre‑flight: add REF check badges (via `/api/evo/refcheck`), model warmup ETA; scenario presets (RAS‑heavy, TP53‑heavy, VUS)
- Orchestration: single CTA runs staged pipeline (warm → score chunked → evidence for top‑K → compare → aggregate → persist)
- Evidence: per‑variant chips + mini sparkline profile, probe table, ClinVar/Ensembl links
- Dual-model: toggle + agreement rate; highlight disagreements
- Exports: request.json, response.json, results.csv, “reproduce run” curl
- Job log ribbon: time-stamped stage events; chunk progress

## 7) Benchmarking & Data
- `scripts/generate_clinvar_panel.py`: high-confidence GRCh38 SNVs; VCF alleles; optional sample Ensembl validation; balanced per gene; provenance CSV
- `scripts/run_myeloma_benchmark.py`: chunked 20/req; computes AUROC/AUPRC and calibration; writes raw/metrics/summary artifacts
- Nightly workflow ready to upload artifacts; extend to 40B and bootstrapped CIs

## 8) Analytics (Supabase)
- Tables:
  - `mdt_runs(run_signature, model_id, prediction, ras_sum, tp53_sum, num_variants, upstream, alt_model, agree_rate, created_at)`
  - `mdt_run_variants(id, run_signature, gene, chrom, pos, zeta, min_delta, exon_delta, confidence, call, raw, created_at)`
  - (plan) `mdt_events(run_signature, stage, message, t)` for job logs
- Dashboard (plan): time‑series AUROC/AUPRC; latency histograms; error rates; agreement rate; run drill‑down

## 9) Patent & Compliance Notes (high level)
- Claim the evidence‑first pipeline (methods, systems), confidence calculus, provenance protocol, and serverless orchestration; not the OSS model
- Include formulas, parameter ranges, and artifacts for enablement

References: `/.cursor/rules/patent_claim_doctrine.mdc`

## 10) Partial / Missing (Agentic Gaps)
- Autonomous planner across tools (hypothesis → scoring → design → report)
- Active learning loop (auto-select uncertain/VUS, schedule re-scoring, track drift)
- Cross-model auto-corroboration (7B vs 40B) with disagreement handling in UI/logs (basic version exists; needs automation and UI surfacing)
- Continuous evaluation dashboard (time-series metrics, latency, error rates) from Supabase

## 11) Next 1–2 Sprints to Become “Agentic”
- Background worker (planner)
  - Poll Supabase for new cases → warm models → run chunked scoring → compute calls → export reports
  - Emit job events to `mdt_events`
- Disagreement policy
  - Auto-run 40B on high-risk or low-confidence 7B calls; flag discrepancies in UI/logs
- Auto-benchmark nightly
  - Run fixed 7B & 40B panels; compute AUROC/AUPRC + bootstrapped 95% CIs; persist to Supabase; add simple dashboard view
- Triage loop
  - Generate `vus_panel.json`; score; rank by composite risk; queue for wet-lab or deeper profiling; store outcomes and close loop

## 12) Detailed Task Breakdown
- Backend
  - [x] `/api/evo/refcheck` endpoint calling Ensembl once per variant
  - [x] `/api/twin/run` orchestration endpoint with staged execution and polling handle
  - [x] Supabase runs/variants inserts; [ ] Supabase `mdt_events` dashboard wiring
  - [ ] Evidence archive export (S3/GCS) packing profiles/probes artifacts
- Frontend
  - [ ] Pre‑flight REF badges in `VariantInputList`
  - [ ] Job log ribbon with stage timestamps; chunk progress indicators
  - [ ] Dual‑model comparison UI (agreement bar; disagreement list anchor to variant rows)
  - [ ] Export buttons (CSV/JSON) and “Reproduce run” snippet in `EvidencePanel`
  - [ ] Presets panel (RAS, TP53, VUS)
- Analytics
  - [ ] Nightly 7B & 40B; compute CIs; push to Supabase
  - [ ] Minimal dashboard page: show last N runs metrics + agree_rate trend

## 13) Risks & Mitigations
- External API limits → caching, retries, fast mode; limited pre-flight sample validation
- Serverless timeouts → chunked requests, long read/write timeouts, streaming logs
- Overfitting to ClinVar → hold-out sets; gene-level cross-validation; add orthogonal datasets
- Calibration → add ECE and bootstrapped CIs; avoid overclaiming confidence as probability

## 14) Environment & Config
- Backend env: EVO URLs; SUPABASE_URL/KEY; tables; timeouts; chunk sizes
- Frontend env: VITE_API_ROOT

## 15) Ownership & Timeline (suggested)
- Sprint 1: Backend orchestrator + refcheck ✅; Frontend job log + exports; Nightly 40B; Supabase events
- Sprint 2: Dashboard; dual-model automation; VUS triage worker; presets; calibration CIs

---
This plan operationalizes an agentic, evidence-first Digital Twin capable of end-to-end, audit-ready R&D analysis with continuous evaluation. It defines concrete endpoints, UI/UX, analytics, and sprint-ready tasks to reach autonomous orchestration in 2 sprints. 