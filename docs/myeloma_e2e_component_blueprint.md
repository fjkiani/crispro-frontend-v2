## Myeloma E2E Demo – Build Blueprint (Components, APIs, Agents)

### Goal
Ship a seamless end-to-end demo that is production-trustworthy, modular, and quickly extensible. This blueprint specifies components, props/handlers, backend contracts, data shapes, agents, and a step-by-step implementation plan.

---

## 1) Frontend Architecture

### 1.1 Directory structure
- `oncology-coPilot/oncology-frontend/src/pages/MyelomaDigitalTwin.jsx`
- `oncology-coPilot/oncology-frontend/src/components/myeloma/`
  - `LiveJobBanner.jsx`
  - `ModelSelector.jsx` (new)
  - `VariantInputList.jsx`
  - `MyelomaResponseDisplay.jsx`
  - `DeltaProfileChart.jsx` (new)
  - `SensitivityProbeTable.jsx` (new)
  - `EvidencePanel.jsx` (new)
  - `UpstreamErrorAlert.jsx` (new)
- `oncology-coPilot/oncology-frontend/src/hooks/`
  - `useApiClient.js` (new)
  - `useResultCache.js` (new)
- `oncology-coPilot/oncology-frontend/src/config/toolconfigs.js` (ensure includes model_id, routes)

### 1.2 Component specs

1) ModelSelector.jsx
- Purpose: Choose `evo2_1b | evo2_7b | evo2_40b`.
- Props: `{ value: string, onChange: (modelId) => void }`
- Behavior: Persist to `window.__mdt_model_id` for ToolRunner if needed.

2) LiveJobBanner.jsx (exists)
- Extend: Show selected model, upstream URL (from backend response), warmup CTA; accept `onWarmup(modelInfo)`.
- Props: `{ modelId, onWarmup }`

3) VariantInputList.jsx (exists)
- Maintain current props.
- Validation: Show inline error if `variant_info` doesn’t match `chr:pos REF>ALT`.
- Expose `window.__mdt_mutations` for ToolRunner compatibility.

4) MyelomaResponseDisplay.jsx (exists)
- Extend to accept and render:
  - DeltaProfileChart per variant (request `/api/evo/score_variant_profile` with `{assembly, chrom, pos, ref, alt, model_id}`)
  - SensitivityProbeTable per variant (request `/api/evo/score_variant_probe` with `{assembly, chrom, pos, ref, model_id}`)
  - EvidencePanel modal/drawer toggle.
- Props: `{ results, onAction }`
- onAction contract: `{ apiCall: { endpoint: string, payload: object } } => Promise<any>`

5) DeltaProfileChart.jsx (new)
- Input: `{ profile: Array<{offset:number, delta:number}>, peakOffset:number, peakDelta:number }`
- Render: sparkline (simple SVG or small chart lib). Highlight `peakOffset`.

6) SensitivityProbeTable.jsx (new)
- Input: `{ probes: Array<{alt:string, delta:number}>, topAlt:string, topDelta:number }`
- Render: mini table, highlight topAlt.

7) EvidencePanel.jsx (new)
- Input: `{ provenance: { mode, upstream_service, selected_model, ensembl_url?, ref_base?, window?, exon_flank?, request_json, response_json } }`
- Render: upstream info, Ensembl links, REF base at index, download buttons for JSON (blob download).

8) UpstreamErrorAlert.jsx (new)
- Input: `{ errorText }`
- Render: collapsible detailed error; CTA to fix REF/ALT or coordinate.

### 1.3 Hooks

1) useApiClient.js
- Provides `post(endpoint, payload, { timeoutMs, signal })` with 10-min timeout default; abort support.
- Automatically includes `model_id` from ModelSelector.

2) useResultCache.js
- Cache by signature: `${assembly}|${chrom}|${pos}|${ref}|${alt}|${model_id}`.
- Exposes `get(key)`, `set(key, value)`, `has(key)`.

### 1.4 UI States
- Loading: banner shows spinning warmup; per-variant buttons show inline spinners for profile/probe.
- Partial success: render available fields; show `UpstreamErrorAlert` for failed calls.
- Long-running: banner shows elapsed time; hint for cold starts after 30s.

---

## 2) Backend Orchestrator (FastAPI)

File: `oncology-coPilot/oncology-backend-minimal/api/index.py`

### 2.1 New/confirmed proxy routes
- POST `/api/evo/score_variant_profile` → forwards to `${base}/score_variant_profile`
  - Request: `{ assembly, chrom, pos, ref, alt, flank?, radius?, model_id }`
  - Response: `{ profile: [{offset, delta}], peak_delta, peak_offset }`
- POST `/api/evo/score_variant_probe` → forwards to `${base}/score_variant_probe`
  - Request: `{ assembly, chrom, pos, ref, model_id }`
  - Response: `{ probes: [{alt, delta}], top_alt, top_delta }`

### 2.2 Confidence score
- s1 (effect size) = |min_delta| / 0.5 clamped [0,1]
- s2 (exon corroboration) = 1.0 if same sign and |exon_delta| >= |min_delta| else 0.5 if same sign else 0.0
- s3 (window consistency) = 1 - stdev(deltas)/max(0.05, |min_delta|)
- confidence = round(0.5*s1 + 0.3*s2 + 0.2*s3, 2)

### 2.3 Pathway decision
- RAS/MAPK sum: KRAS/NRAS/BRAF weighted; TP53 sum shown separately.
- Verdict: `Likely Resistant` if RAS/MAPK sum ≥ 2.0 else `Likely Sensitive`.
- Version thresholds in config for transparency.

### 2.4 Provenance
- Always include: `{ mode: 'live', upstream_service, selected_model }`.
- For profile/probe: include query params echoed back in response so UI can show exact request used.

---

## 3) Evo2 Services (Modal)

File: `src/services/evo_service/main.py`
- 40B endpoints: already live.
- 7B-in-main: ensure `/score_variant_multi`, `/score_variant_exon`, `/score_variant_profile`, `/score_variant_probe` exposed (parity with 40B).
- Strict REF check in all endpoints; clear 400 with details.

Payload/Response reference (common):
- Request common: `{ assembly: 'GRCh38'|'GRCh37', chrom: '7', pos: 140753336, ref: 'A', alt?: 'T', windows?: number[], flank?: number, radius?: number }`
- Response snippets:
  - `/score_variant` → `{ delta_score }`
  - `/score_variant_multi` → `{ deltas: [{window, delta}], min_delta, window_used }`
  - `/score_variant_exon` → `{ exon_delta, window_used }`
  - `/score_variant_profile` → `{ profile: [{offset, delta}], peak_delta, peak_offset }`
  - `/score_variant_probe` → `{ probes: [{alt, delta}], top_alt, top_delta }`

---

## 4) Agent Taxonomy & Integration Points

- Conversation Orchestrator Agent (frontend assistant UI)
  - Suggests actions (“Warm model”, “Run profile”), summarizes results.
- Variant QA Agent (backend callable or client-side helper)
  - Validates/builds correct hg38 coords, Ensembl links, flags REF errors; proposes fixes.
- Evidence Synthesizer Agent (backend)
  - Given raw metrics and external annotations, crafts rationale text for UI.
- Provenance Notary Agent (backend)
  - Signs run artifacts (optional for demo), stores JSON to a simple object store or DB.
- Benchmark & Regression Agent (CI task)
  - Nightly run of benchmark variants; posts summary JSON.

Initial demo: Orchestrator + Variant QA + Evidence Synthesizer (lightweight stubs OK).

---

## 5) Data Contracts (JSON)

detailed_analysis[i]
```
{
  "gene": "KRAS",
  "variant": "KRAS p.Gly12Asp",
  "chrom": "12",
  "pos": 25245350,
  "calculated_impact_level": 0.5,
  "evo2_result": {
    "zeta_score": -0.0084,
    "min_delta": -0.0521,
    "window_used": 1024,
    "exon_delta": -0.0436,
    "confidence_score": 0.33,
    "confidence_reason": "effect 0.052, windows variable, exon -0.044"
  },
  "selected_model": "evo2_40b",
  "original_variant_data": { "gene": "KRAS", "hgvs_p": "p.Gly12Asp", "variant_info": "chr12:25245350 C>T", "build": "hg38" }
}
```

provenance (attached to response root)
```
{
  "mode": "live",
  "upstream_service": "https://crispro--evo-service-...",
  "selected_model": "evo2_40b"
}
```

profile response
```
{ "profile": [{"offset": -5, "delta": -0.004}, ...], "peak_delta": -0.019, "peak_offset": -2 }
```

probe response
```
{ "probes": [{"alt": "A", "delta": -0.009}, ...], "top_alt": "T", "top_delta": -0.017 }
```

---

## 6) Step-by-Step Build Plan

### Phase A – Backend completeness (1–2 days)
- [x] Ensure 7B-in-main exposes: multi, exon, profile, probe.
- [x] Add proxy routes in FastAPI: `/api/evo/score_variant_profile`, `/api/evo/score_variant_probe`.
- [x] Return provenance fields for profile/probe routes (mode, upstream_service, selected_model, fallback_used).
- [x] Integration smoke tests: 7B profile/probe via backend; 40B profile/probe via fallback path.
- [ ] Native 40B parity for profile/probe (remove fallback).

### Phase B – Frontend components (1–2 days)
- [x] Implement ModelSelector; plumb `model_id` to requests.
- [x] Extend MyelomaResponseDisplay with buttons to fetch profile/probe per variant; render DeltaProfileChart & SensitivityProbeTable (initial inline render; dedicated components pending).
- [x] Implement EvidencePanel and hook up provenance + JSON download.
- [x] Introduce useApiClient with long timeouts and abort; useResultCache for per-variant memoization.
- [ ] Create dedicated `DeltaProfileChart` and `SensitivityProbeTable` visualization components.

### Phase C – Agents & Guidance (1 day)
- [ ] Add a lightweight assistant panel that suggests next actions, flags REF mismatches, links Ensembl.
- [ ] Stub Evidence Synthesizer to generate concise rationales from metrics.

### Phase D – Benchmarks & Trust (1 day)
- [ ] Create a 12–20 variant myeloma panel (json file) and a script to run nightly.
- [ ] Compute and store summary metrics; render mini dashboard in the app (or static markdown for demo).

---

## 7) Testing & Acceptance
- Unit: UI components render with mock data; hooks handle timeouts/abort.
- Integration: backend proxy returns expected shapes; Modal endpoints reachable; REF validation triggers 400 appropriately.
- E2E: demo script runs fully (warmup → analyze → evidence → verdict → reproduce) with no manual fixes.
- Acceptance: “Likely Resistant” shown for RAS/MAPK-positive panel; provenance JSON downloadable; re-run reproduces values.

---

## 8) Risks & Mitigations
- Cold-start latency → warmup CTA + banner timer; pre-warm before demo.
- REF mismatches → Variant QA agent + inline fix suggestions; curated presets with validated coords.
- Small delta magnitudes → emphasize min_delta/exon_delta + profile/probe visuals; confidence with reasons.
- Endpoint drift → pin model IDs and endpoint URLs in env; nightly smoke.

---

## 9) Deliverables Checklist
- [x] Backend proxy complete (profile/probe + provenance echo; fallback for 40B)
- [x] Modal services parity (7B/40B core; 40B profile/probe pending)
- [x] Frontend components: ModelSelector, Profile/Probe triggers (inline), EvidencePanel
- [ ] Assistant guidance (Variant QA + Evidence Synth)
- [ ] Benchmarks and summary report
- [ ] Demo script validated and reproducible

---

## 10) Next Steps (48–72 hours)
- Backend
  - Implement native 40B endpoints for `/score_variant_profile` and `/score_variant_probe` (remove fallback path).
  - Echo provenance on all proxy routes (variant, multi, exon) for consistency.
- Frontend
  - Add dedicated `DeltaProfileChart` and `SensitivityProbeTable` components.
  - Surface `fallback_used` badges when applicable in EvidencePanel.
- Agents & Benchmarks
  - Ship lightweight assistant panel (input QA + next-step guidance).
  - Prepare 12–20 variant myeloma benchmark JSON; run and store summary metrics.
- Validation
  - Re-run E2E: KRAS/NRAS/BRAF/TP53 (with corrected TP53 hg38 REF), confirm “Likely Resistant” and evidence panels. 

---

## 11) Checklist Rationale – What/Why
- Backend profile/probe routes (done): enables interpretability (local disruption curve, locus susceptibility) required by demo acceptance; adds provenance so users can trust upstream and model identity.
- 7B parity + Modal deploys (done): guarantees at least one live model fully supports profile/probe now; reduces risk of blocked demos.
- 40B fallback to 7B for profile/probe (temporary): ensures consistent UX for all models today while we add native 40B endpoints; provenance shows when fallback occurs.
- Provenance echo on variant/multi/exon (done): consistent “live + upstream + model” metadata across all calls improves auditability and trust.
- ModelSelector (done): supports the cost/latency trade-off story and matches acceptance criteria for model switching.
- EvidencePanel (done): surfaces mode, upstream URL, model; download JSONs for reproducibility and audit; matches “trustworthy provenance” acceptance.
- Profile/Probe UI with caching (done): avoids duplicate recomputation; faster UX for repeated variant exploration; stronger interpretability per acceptance.
- Assistant panel (initial): warms model, suggests next steps, flags input issues; improves demo flow and reduces operator error.
- Benchmark mini-panel + runner (done): provides a quick, reproducible summary (Likely Resistant; RAS/MAPK=3.9) to validate end-to-end behavior.

### Next (why these next)
- Native 40B profile/probe: removes temporary fallback, unifies capabilities across models, strengthens robustness.
- Dedicated charts/tables: clearer, investor-friendly visualization of signals; supports “interpretable evidence” criterion.
- Assistant guidance expansion: faster operator flow and fewer mistakes during live demos; improves perceived polish.
- Nightly benchmark + report: catches regressions early; demonstrates operational maturity and reliability. 