# Agent Handoff — **Backend Core (A–Y)** for Tests (Axes Coverage) + CoPilot (RUO)

## What exists today (receipts in code)

### Core backend endpoints in scope
- **Therapy-Fit scenarios**: `GET /api/ayesha/therapy-fit/scenarios`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py`
  - Preview computation: `oncology-coPilot/oncology-backend-minimal/api/services/ayesha_fit/precompute.py`
- **Therapy-Fit analyze**: `POST /api/ayesha/therapy-fit/analyze?level=all`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py`
  - Pipeline execution: `oncology-coPilot/oncology-backend-minimal/api/services/ayesha_fit/service.py` → `EfficacyOrchestrator`
- **Ayesha trials**: `POST /api/ayesha/trials/search`
  - Router exists (legacy): `oncology-coPilot/oncology-backend-minimal/api/routers/_ayesha_trials_legacy.py`
- **Ayesha complete care (compat)**: `POST /api/ayesha/complete_care_plan`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha.py`
- **Trials agent (generic)**: `POST /api/trials/agent/search`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/trials_agent.py`

## Backend workplan (A–Y). Leave **Z** for frontend agent.

### ✅ Receipt Pack (required evidence)
**Rule:** Every “✅ implemented” claim must include at least one of:
- (a) `openapi.json` proof the path exists, or
- (b) a grep/rg receipt in code, or
- (c) a captured JSON payload snippet.

#### Receipt 1 — OpenAPI path existence
Proof (via `/openapi.json`) that these paths exist in the running app:
- `/api/ayesha/therapy-fit/analyze` ✅
- `/api/ayesha/therapy-fit/scenarios` ✅
- `/api/ayesha/therapy-fit/health` ✅
- `/api/ayesha/trials/search` ✅
- `/api/trials/agent/search` ✅

#### Receipt 2 — Real `analyze?level=all` JSON (canonical keys + unknown encoding)
Captured snippet from `POST /api/ayesha/therapy-fit/analyze?level=all` (clipped to axis keys + meta):

```json
{
  "L1": {
    "pathway_scores": { "ddr": null, "mapk": null, "pi3k": null, "vegf": null, "her2": null, "io": null, "efflux": null },
    "pathway_scores_meta": { "normalized": true, "defaulted": true, "defaulted_keys": ["ddr","mapk","pi3k","vegf","her2","io","efflux"], "normalization_version": "v1-null-defaults" }
  },
  "L2": {
    "pathway_scores": { "ddr": 0.003, "mapk": null, "pi3k": null, "vegf": null, "her2": null, "io": null, "efflux": null },
    "pathway_scores_meta": { "normalized": true, "defaulted": true, "defaulted_keys": ["mapk","pi3k","vegf","her2","io","efflux"], "normalization_version": "v1-null-defaults" }
  }
}
```

#### Receipt 3 — Real `scenarios` JSON (mechanism panel presence)
Captured snippet from `GET /api/ayesha/therapy-fit/scenarios` (clipped):
- Top-level keys include: `l2_scenarios`, `l3_scenarios`, `preview_cache`
- Sample L2 card includes `mechanism_panel.mechanism_axes` (7) and `mechanism_panel.current_mechanism_vector` (7)

#### Receipt 4 — Real `health` JSON (preview cache metadata)
Captured snippet from `GET /api/ayesha/therapy-fit/health` (clipped):

```json
{
  "status": "ok",
  "preview_cache": {
    "status": "stale",
    "generated_at": "2026-02-11T03:00:41Z",
    "ttl_seconds": 3600,
    "scenario_version_hash": "0c02636456d3593a65b18d1d3a3173f0",
    "pipeline_version_hash": "34a6fd72a54619a8dc7de453f7281106",
    "persisted": true
  }
}
```


### A) Canonical endpoint map (declare the source of truth)
**Goal:** Remove ambiguity for CoPilot + FE by explicitly declaring canonical endpoints and keeping aliases stable.

- **Therapy-Fit (Tests page data)**:
  - Canonical: `GET /api/ayesha/therapy-fit/scenarios`, `POST /api/ayesha/therapy-fit/analyze`
- **Ayesha trials**:
  - Canonical: `POST /api/ayesha/trials/search`
  - Keep also: `POST /api/trials/agent/search` (generic)
- **Ayesha complete care**:
  - Keep stable: `POST /api/ayesha/complete_care_plan` (compat)
  - Also exists: `POST /api/ayesha/complete_care_v2`, `POST /api/complete_care/v2`

**Acceptance:** A short backend-facing doc comment in:
- `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha.py`
- `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py`
Naming canonical paths and what aliases exist for compatibility.

**Status:** ⏳ Pending (documentation-only; no runtime change).

### B) Guarantee Therapy-Fit scenario payload contains mechanism panel vectors
**Goal:** The FE must always be able to compute axis status from structured fields when previews exist.

In `GET /api/ayesha/therapy-fit/scenarios`, ensure each scenario card surfaces:
- `mechanism_panel.mechanism_axes` (canonical ordering)
- `mechanism_panel.current_mechanism_vector` (0–1 floats)

**Receipt anchor:** `api/services/ayesha_fit/precompute.py:compute_mechanism_panel()` already returns these keys, and `api/routers/ayesha_therapy_fit.py:/scenarios` already exposes `mechanism_panel`.

**Acceptance:** No scenario card returns `mechanism_panel=null` unless its preview run status is explicitly failed/degraded, in which case:
- keep `mechanism_panel=null`
- include a deterministic `preview_status` and/or `error` string (verbatim).
**Status:** ✅ Mechanism panel vectors are present; ⏳ `preview_status` enum is not standardized yet (optional).

### C) Guarantee L3 preview entries can carry a “tripwire error” verbatim
**Goal:** Enable Alert D (“Expression Coverage Too Sparse”) without FE inference.

In L3 matrix entries (`preview_matrix.by_l2[l2_id]`), keep:
- `error: string` when preview computation triggers expression sparsity tripwire(s)
- OR `mechanism_panel.provenance.tripwire_error: string`

**Acceptance:** FE can display the exact string (no paraphrase), and the string originates from backend only.
**Status:** ✅ `error` slot present in L3 matrix entries (verbatim string contract).

### D) Guarantee analyze returns `inputs_used` exactly (for completeness gates)
**Goal:** Missing-data gates must be computed from structured `inputs_used`.

In `POST /api/ayesha/therapy-fit/analyze`, ensure each returned level includes:
- `inputs_used.mutations[]`
- `inputs_used.tumor_context{ ... }`

**Receipt anchor:** `api/routers/ayesha_therapy_fit.py` already exposes `inputs_used` per level.
**Status:** ✅ Implemented.

### E) Stabilize `pathway_scores` axis key naming across levels
**Goal:** Avoid FE special-casing by making axis keys consistent.

For `pathway_scores` returned under `analyze.*.pathway_scores`:
- prefer keys: `ddr, mapk, pi3k, vegf, her2, io, efflux` (lowercase)
- if upstream provides mixed keys, adapt in backend (non-breaking) by adding normalized duplicates.

**Acceptance (updated; “no hallucination” safe):**
- A single run of `analyze?level=all` always includes these seven keys.
- Canonical key values are **`number|null`**:
  - `number` means known/computed upstream
  - `null` means unknown/not computed/invalid
- A `pathway_scores_meta` object is returned with:
  - `normalized: true`
  - `defaulted: boolean`
  - `defaulted_keys: string[]` (the canonical keys that are `null`)
  - `upstream_keys_present: string[]` (for debugging)
  - `normalization_version: string`

✅ Implemented (non-breaking): `api/routers/ayesha_therapy_fit.py` now normalizes pathway scores with **null defaults + meta**:
- Router helper: `_normalize_pathway_scores_with_meta(...)`
- Version tag: `normalization_version = "v1-null-defaults"`

**Code receipts (grep):**
- `_normalize_pathway_scores_with_meta` is defined in `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py`
- `pathway_scores_meta` is returned from both `/bundle` and `/analyze` payload assembly via `**_normalize_pathway_scores_with_meta(...)`

⚠️ Frontend follow-up (Z): current FE numeric coercion treats `null` as `0` in some paths; FE must treat `null` as unknown and gray out axes listed in `pathway_scores_meta.defaulted_keys`.

### F) Preserve “no hallucination” by never adding derived fields that look like measurements
**Goal:** Avoid backend “helpfulness” that becomes frontend overclaim.

Allowed backend additions:
- `provenance.*`
- `*_status` fields that are **explicit**, bounded enums (e.g., `ok|degraded|failed`)
- `error` fields that are verbatim strings

Not allowed (for this contract):
- inferred axis statuses (“High/Low”) unless they are direct translations of a numeric axis field present in payload.
**Status:** ✅ In scope changes comply (we only added meta + nulls; no derived biological claims).

### G) CoPilot backend compatibility shims (only if needed)
**Goal:** If CoPilot is routing to any endpoint that isn’t present in `api/main.py`, add a thin alias router.

Verify these exist in OpenAPI and `api/main.py` router includes:
- `/api/ayesha/complete_care_plan` (present)
- `/api/ayesha/trials/search` (present in legacy trials router)
- `/api/ayesha/therapy-fit/scenarios` and `/api/ayesha/therapy-fit/analyze` (present)
- `/api/evidence/literature`, `/api/evidence/deep_analysis`, `/api/evidence/rag-query` (present)
- `/api/evo/score_variant_profile` (present)

**Acceptance:** CoPilot can call listed endpoints without 404s in a default dev boot.

✅ Bug fix implemented (blocking): `api/services/ayesha_fit/service.py` was throwing `NameError: _norm_drug_name not defined` in the MBD4 “biohacking” injection path, causing `therapy-fit/analyze` to 500. Added a deterministic local `_norm_drug_name()` helper to unblock contract.
**Status:** ✅ Implemented.

### H) Add health endpoints (optional but recommended)
**Goal:** Make integration debugging easy and deterministic.

Add (or confirm existing):
- `GET /api/ayesha/therapy-fit/health` → includes preview cache status, scenario hash, pipeline hash
- `GET /api/ayesha/complete_care_v2/health` (already exists in `ayesha_orchestrator_v2.py`)
**Status:** ✅ Implemented (`api/routers/ayesha_therapy_fit.py`).

### I) Confirm preview cache persistence defaults are safe
**Goal:** Preview scenarios shouldn’t thrash on every server restart.

Preview cache controls live in:
- `api/services/ayesha_fit/precompute.py` (env vars: `AYESHA_THERAPY_FIT_PREVIEW_CACHE_*`)

**Acceptance:** Defaults keep cache persisted (unless explicitly disabled), and cache metadata returns:
- `scenario_version_hash`
- `pipeline_version_hash`
- `generated_at`, `ttl_seconds`
**Status:** ✅ Verified (code + runtime receipts)
- Code defaults in `api/services/ayesha_fit/precompute.py`:
  - `AYESHA_THERAPY_FIT_PREVIEW_CACHE_PERSIST` defaults to `"1"` (enabled)
  - `AYESHA_THERAPY_FIT_PREVIEW_CACHE_DISABLE` defaults to `"0"` (enabled)
- Runtime receipt: `GET /api/ayesha/therapy-fit/health` returns `persisted: true` and the expected hash+TTL fields.

### J–Y) Reserved backend slots (fill only if required)
These are intentionally reserved to avoid scope creep. Only implement if integration breaks:
- J) Add explicit `contract_version` bump if response keys change
- K) Add server-side validation errors for malformed `scenario_id` / `l3_scenario_id`
- L) Add deterministic error codes (in addition to strings) for UI rendering
- M–Y) (unused)

**Status:** Not needed yet.

---

## Confidence / completeness (honest)

### Do I feel “100%”?
- **Backend contract stability for pathway scores:** **Yes** (keys always present; unknowns are `null`; meta tells FE what was defaulted).
- **End-to-end UI correctness without frontend changes:** **No** until FE updates its numeric parsing to treat `null` as unknown (Z).

### Manager questions (only if we want to extend backend beyond current scope)
These are not blockers for today’s backend contract, but would reduce future drift:
1. **Normalization naming policy**: keep canonical keys lowercase (`ddr`) vs also provide canonical display keys (`DDR`). (Today: lowercase + preserve all upstream keys.)
2. **Contract versioning**: should we bump `contract_version` (e.g., `v2.1`) now that we added `pathway_scores_meta`?
3. **Health endpoint requirement**: do we want `GET /api/ayesha/therapy-fit/health` mandated for demo readiness (yes/no)?
4. **Error coding**: should tripwire errors get a stable `error_code` in addition to verbatim string for UI logic?
5. **Scope of meta propagation**: should the same `*_meta` pattern apply to other “coverage vs value” fields beyond pathway scores (e.g., IO expression signatures), or keep limited to Tests page only?

---

## Z) Frontend-only (handled by frontend agent)
All UI/UX polish, deep-links, scenario selector improvements, and CoPilot chat rendering are explicitly out-of-scope for this backend handoff.

