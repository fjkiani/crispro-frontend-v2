# 🔧 Plumber Backlog: Production Resistance Capability

This is the low-level execution backlog to make the resistance capability **production-grade**, **testable**, and **non-drifting**.

**Default decisions are provided below for Q2/Q3** so the plumber can proceed immediately. Manager can override later if desired (these choices affect every endpoint contract).

**Execution order:** P0 → P1 → P2. Within P0, contract unification (Q1-Q4) must happen before validation runner work, because validation scripts need stable contracts to test against.

---

## 🎯 Next 15 Deliverables (Plumber’s focus; ignore the rest until these are green)

**Rule:** each deliverable must end with a **receipt** (Ring‑0 tests or Ring‑1 `report.json`) + **one copy/paste command** to reproduce.

### Task A — Scaffolding + Validation Hardening (Deliverables 1–8)
1. **`requirements-dev.txt`** (or single documented install path) exists for Ring‑0
2. **Ring‑0 canonical command** documented + not silently skipped
3. **Fixtures folder** created (3–5 JSON fixtures for L0/L1/L2, CA‑125 present/missing, HRD present/missing, expression present/missing)
4. **Fixture-driven Ring‑1 validator** created (consumes fixture pack; emits `out/<contract_id>/report.json`)
5. **Copy‑on‑write receipts** implemented across Ring‑1 validators (`report_<timestamp>.json` + copy to `report.json`)
6. **Ring‑1 suite stays green** after copy‑on‑write change
7. **Ring‑1 suite “strict mode”**: fails if pytest missing *when requested* (no silent success)
8. **One “fixture pack README”** explaining each fixture + expected behaviors

### Task B — Reliability refactor (Deliverables 9–12)
9. **Remove internal localhost HTTP** from `complete_care_universal.py::_extract_insights_bundle()` (direct service calls)
10. **Ring‑1 validator proving “no localhost HTTP required”** (fixture-based)
11. **Provenance versioning fields present** (`code_version`, `contract_version`, `inputs_snapshot_hash`) in Ring‑1 receipt
12. **Centralized config of timeouts/flags** for orchestrator/internal clients (no magic numbers scattered)

### Task C — End‑to‑end resistance detection MVP (Deliverables 13–15)
13. **E2E OV “`Ayesha`-like” fixture** (mutations-only required; optional CA‑125/HRD; expression optional)
14. **Ring‑1 E2E receipt** emitted at:
    - `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/resistance_e2e_ov_ayesha_v1/report.json`
15. **MVP output contract verified**: Risk + Mechanisms + Actions (no trials/SL in MVP), with explicit missing-input flags + caps

---

## ✅ Plumber Execution Contract (A / B / C)

**Purpose:** Stop “code for code’s sake.” Every task must end with **proof artifacts** (Ring‑0/Ring‑1 receipts) and a single command that reproduces them.

### What Zo will accept on every PR (non‑negotiable)
- **PR description includes**: Why / What changed / How to verify
- **Receipt attached** (one of):
  - Ring‑0: tests
  - Ring‑1: `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/<contract_id>/report.json`
- **One copy‑paste command** that reproduces the receipt locally.
- **No new endpoints** unless explicitly required by the task.

---

### Task A (Plumber): Scaffolding + Validation Hardening (must finish before any more features)
**Goal:** Make the current resistance system **provable** and **non‑drifting** with fixtures and stable receipts.

**Scope (must-do):**
- **A1 — Ring‑0 runnable everywhere**
  - Create `oncology-coPilot/oncology-backend-minimal/requirements-dev.txt` (or a documented single install path) so Ring‑0 is not environment-dependent.
  - Provide the canonical Ring‑0 command in this doc.
- **A2 — Fixture pack for “end‑to‑end resistance detection”**
  - Add a minimal fixture set (3–5 JSON fixtures) that cover:
    - L0 vs L1 vs L2 completeness behavior (caps + flags)
    - “2-of-3 trigger” edge cases
    - CA‑125 history present vs missing
    - HRD score present vs missing (optional)
    - expression present vs missing (RUO; **must not be required for MVP**)
  - At least one Ring‑1 validator must consume these fixtures and emit `report.json`.
- **A3 — Report stability (copy-on-write)**
  - Preserve history while keeping a stable canonical pointer:
    - Write `report_<timestamp>.json`, then copy to `report.json`
  - Apply to all Ring‑1 validators.

**Acceptance (Task A is “done” only if all are true):**
- Ring‑1 suite stays green:
  - `python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/run_resistance_validation_suite.py --ring1`
- Ring‑0 runs via a single documented command (not silently skipped unless explicitly requested).
- A new fixture-driven Ring‑1 validator emits a canonical `out/<contract_id>/report.json`.

**Verification commands (plumber must keep these working):**
- Ring‑1:
  - `python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/run_resistance_validation_suite.py --ring1`
- Ring‑0 (canonical target; update if different):
  - `cd oncology-coPilot/oncology-backend-minimal && PYTHONPATH=. python3 -m pytest -q`

---

### Task B (Plumber): Reliability Refactor (reduce brittleness, keep contracts stable)
**Goal:** Remove hidden failure modes so “end‑to‑end resistance detection” doesn’t depend on localhost HTTP or fragile wiring.

**Scope (must-do):**
- **B1 — Remove internal `http://localhost:8000` calls where possible**
  - Start with `complete_care_universal.py::_extract_insights_bundle()` and move to direct service calls.
  - Keep response shape identical (additive provenance only).
- **B2 — Centralize provenance versioning**
  - Ensure `ResistanceContract.provenance` includes stable version fields:
    - `code_version` (commit hash or build-time string)
    - `contract_version` (semantic or timestamp+hash)
    - `inputs_snapshot_hash`

**Acceptance:**
- Ring‑1 suite stays green.
- A Ring‑1 validator proves the tested orchestrator path does not require localhost HTTP.
- Provenance fields are present and stable in a fixture-based receipt.

---

### Task C (Plumber): End‑to‑End Resistance Detection Capability (product path)
**Goal:** A single, reproducible “`Ayesha`-like” run that returns: **risk → mechanisms → actions**, with evidence tiers and provenance.

**Scope (must-do):**
- **C0 — Canonical product surface (non-negotiable)**
  - **Use**: `POST /api/complete_care/v2` (Universal Orchestrator) as the patient-facing surface.
  - `/api/resistance/predict` remains **dev/testing only**.
- **C1 — One canonical end‑to‑end run (fixture-driven, mutations-only MVP)**
  - Create a deterministic “`Ayesha`-like OV” fixture input with:
    - **Required**: tumor mutations (VCF/MAF-parsed list)
    - **Optional**: CA‑125 history, HRD score
    - **Not required for MVP**: expression data (MFAP4/EMT gated off if absent)
  - Run through `POST /api/complete_care/v2` and emit a Ring‑1 `report.json` with:
    - **Risk**: platinum resistance probability (0–1) + category + confidence cap
    - **Mechanisms**: MAPK / PI3K / DDR / CA‑125 (if present) with receipts + limitations
    - **Actions**: monitoring plan + treatment actions (no trials/SL in MVP)
- **C2 — Gold label + validation contract (Ring‑2 target)**
  - **Gold label**: Platinum-Free Interval (PFI) < 6 months
  - Endpoint name (canonical): `platinum_resistance_pfi_lt_6mo`
  - Cohort: TCGA‑OV (n≈469)
  - Metrics: AUROC + sensitivity + specificity
  - Receipt target path: `data/validation/reports/ov_platinum_pfi_validation.json`
- **C3 — MVP output scope (do not expand)**
  - MVP output is **Risk + Mechanisms + Actions** only.
  - **Do not include** trials or synthetic lethality in MVP output (Week‑2+).

**Acceptance:**
- One command produces a deterministic Ring‑1 receipt for an end‑to‑end “`Ayesha`-like” case under:
  - `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/resistance_e2e_ov_ayesha_v1/report.json`
- Receipt explicitly shows missing inputs + caps confidence appropriately (no overclaiming).
- The run uses the **universal orchestrator** (`/api/complete_care/v2`) as the canonical surface.
- “Expression missing” is treated as a **data gap**, not an error.

---

## 🧭 Oversight loop (Zo → Plumber)

**Daily status update (required):**
- **What changed** (1–3 bullets)
- **Receipts produced** (paths to `report.json` and/or tests run)
- **Command to reproduce** (exact)
- **What’s next** (1–3 bullets)

**Stop condition (no more work until resolved):**
- Ring‑1 is red, or
- receipt cannot be reproduced with the documented command, or
- a change adds endpoints/contracts outside Task A/B/C scope.

### Answer Status Summary

| Question | Status | Action Needed |
|----------|--------|---------------|
| **Q2** (Evidence tier mapping) | ✅ **DEFAULT SET** | Manager only needs to override if we want `TREND` treated as Tier 3 (default: Tier 4) |
| **Q3** (Disease normalization) | ✅ **DEFAULT SET** | Manager only needs to override if we want strict 422 everywhere (default: safe degradation, no risky fallbacks) |
| **Q1** (Contract unification) | ✅ **RECOMMENDATION PROVIDED** | Option A (keep all, unify schema) — can proceed |
| **Q4** (Backward compatibility) | ✅ **PARTIALLY ANSWERED** | Option B (version endpoints) recommended — audit callers first |
| **Q5** (Test dependencies) | ✅ **ANSWERED** | pytest already in requirements.txt — low priority cleanup |
| **Q6** (Report stability) | ✅ **PARTIALLY ANSWERED** | Current approach works but loses history — medium priority enhancement |
| **Q7** (Orchestrator refactoring) | ✅ **ANSWERED** | Option B (incremental) — start with `_extract_insights_bundle()` refactor |
| **Q8** (Data plumbing priority) | ✅ **ANSWERED** | OV expression first (confirmed by strategy docs) |
| **Q9** (Fixture scope) | ✅ **RECOMMENDATION PROVIDED** | Option C (iterative) — start minimal, expand |
| **Q10** (Provenance versioning) | ✅ **RECOMMENDATION PROVIDED** | Option C (timestamp + hash) with multi-field approach |

### Quick reference: which questions block which tasks

| Question | Blocks | Priority | Status |
|----------|--------|----------|--------|
| **Q2** (Evidence tier mapping) | All P0 contract work | 🔴 CRITICAL | ✅ Default set |
| **Q3** (Disease normalization) | All P0 contract work | 🔴 CRITICAL | ✅ Default set |
| **Q1** (Contract unification strategy) | Output schema unification | 🟡 HIGH | ✅ Can proceed |
| **Q4** (Backward compatibility) | Output schema unification | 🟡 HIGH | ✅ Can proceed (audit first) |
| **Q5** (Test dependencies) | Ring-0 validation runner | 🟢 MEDIUM | ✅ Answered |
| **Q6** (Report stability) | Report output stabilization | 🟢 MEDIUM | ✅ Can proceed |
| **Q7** (Orchestrator refactoring) | P1 orchestrator hardening | 🟢 MEDIUM | ✅ Can proceed |
| **Q8** (Data plumbing priority) | P1 data plumbing | 🟢 MEDIUM | ✅ Answered |
| **Q9** (Fixture scope) | Fixture creation | 🟢 LOW | ✅ Can proceed |
| **Q10** (Provenance versioning) | Provenance centralization | 🟢 LOW | ✅ Can proceed |

---

## P0 — Contract & duplication cleanup (must-do)

- **Unify “resistance output schema”**
  - Create one canonical schema for:
    - `endpoint`, `risk`, `mechanisms[]`, `actions[]`, `receipts[]`, `provenance`
  - Apply to:
    - `POST /api/care/resistance_playbook`
    - `POST /api/resistance/predict`
    - `POST /api/complete_care/universal` (when `include_resistance*` toggles are enabled)

- **Unify disease normalization**
  - Today there are multiple disease labels (`ovarian`, `ov`, `ovarian_cancer_hgs`, `multiple_myeloma`, `mm`, etc.).
  - Make a single normalizer and remove risky defaults (e.g., “unknown → ovarian” fallback should become explicit error or “unknown”).

- **Unify evidence-tier vocabulary**
  - Map/merge:
    - `EvidenceLevel` in `resistance_playbook_service.py`
    - evidence fields in `resistance_prophet_service.py`
    - manager tiers (Tier 1–5 in `Manager.mdc`)
  - Output must never contradict itself (e.g., “VALIDATED” vs “LITERATURE_ONLY” vs “TREND”).

---

## P0 — Make validations runnable (Ring 0/1)

- **Tooling: make Ring-0 runnable everywhere**
  - Ensure dev/CI has a standard way to install test deps (at minimum: `pytest`).
  - Add a documented command in `oncology-coPilot/oncology-backend-minimal/requirements.txt` or a `requirements-dev.txt` (preferred) so “Ring-0” is not environment-dependent.

- **Add a deterministic validation runner CLI**
  - Runs without server + without network.
  - Executes:
    - unit tests (resistance playbook, biomarker intelligence, safety)
    - deterministic validators that emit `report.json` (e.g., `validate_ov_nf1_playbook.py`)
  - Produces a single summary JSON and a stable location for canonical reports.

- **Stabilize report outputs**
  - Timestamped reports are good for auditing, but we also need a stable “latest” pointer (or a canonical report file per contract).
  - Add: `out/<contract_id>/report.json` as the canonical file.

- **Add fixtures for ResistanceProphet + ResistanceDetectionService**
  - Small JSON fixtures that cover:
    - L0/L1/L2 completeness behavior
    - “2-of-3 trigger” edge cases
    - confidence caps + baseline penalty behavior

---

## P1 — Orchestrator hardening (product reliability)

- **Remove internal `http://localhost:8000` calls where possible**
  - Orchestrators currently call other endpoints via HTTP for “insights extraction”.
  - Refactor toward direct function/service calls for:
    - determinism
    - testability
    - performance

- **Centralize provenance**
  - Ensure every surface returns:
    - `run_id`
    - `versions` (code version / contract version)
    - `flags`
    - `inputs_snapshot_hash`

---

## P1 — Data plumbing (to unlock stronger validation)

- **OV expression ingestion path**
  - Define where expression comes from (file upload / integration) and how it’s stored.
  - Needed for MFAP4/EMT rules to be usable in product, not only in scripts.

- **MM CNV/cytogenetics integration**
  - Define input schema for:
    - CNV calls
    - cytogenetics (del(17p), t(4;14), 1q gain, t(11;14))
  - Add fixture datasets and validation scripts.

---

## P2 — Nice-to-haves (after MVP is stable)

- **CI wiring**
  - Ring-0 on every PR (fast)
  - Ring-1 on promotion PRs
  - Ring-2 scheduled (network/cBioPortal refresh)

- **Docs**
  - “How to add a new endpoint contract”
  - “How to add a new biomarker rule with receipts”
  - “How to interpret report.json”

---

---

## ❓ Manager Questions (blocking execution decisions)

### Q1: Contract unification strategy
**Context:** We have 3 resistance endpoints (`/api/care/resistance_playbook`, `/api/resistance/predict`, orchestrators). They overlap but serve different use cases.

**Question:** Should we:
- **Option A:** Keep all 3, but make them share a single output schema (backward-compatible wrapper)?
- **Option B:** Deprecate one endpoint (which one?) and redirect to a canonical endpoint?
- **Option C:** Create a new unified endpoint and mark the old ones as legacy (with migration path)?

**Recommendation:** Option A (keep all, unify schema) — least breaking, fastest to ship.

---

### Q2: Evidence tier mapping (critical for product honesty) ✅ **DEFAULT SET**
**Context:** We currently have multiple evidence vocabularies across the stack.

**Default decision (so plumber can proceed):**
- **Resistance outputs use Manager evidence tiers (Tier 1–5) as the single product-facing vocabulary.**
- **Drug efficacy uses its own evidence banding (`supported/consider/insufficient`) and must NOT be mapped into Tier 1–5.**

**Why:** Tier 1–5 is the language we use to communicate clinical certainty; efficacy banding is an internal S/P/E confidence gate and answers a different question.

**Schema consequence (prevents drift):**
- Resistance mechanisms/actions must use: `evidence_tier: TIER_1|TIER_2|TIER_3|TIER_4|TIER_5`
- Drug efficacy items must use: `evidence_band: supported|consider|insufficient` (separate field name)

**Default mapping (apply immediately):**
- `STANDARD_OF_CARE` → **Tier 1**
- `CLINICAL_TRIAL` → **Tier 2**
- `VALIDATED` (p<0.05, pinned receipt) → **Tier 3**
- `TREND` (p<0.15 OR underpowered) → **Tier 4** (explicitly *not validated*)
- `LITERATURE_BASED` / `LITERATURE_ONLY` → **Tier 4**
- `PRECLINICAL` / `LOW_EVIDENCE` / `EXPERT_OPINION` → **Tier 5**

**Manager override (only if desired):**
- Do we ever allow `TREND` to be treated as Tier 3? **Default: NO**. If you want exceptions, require a pinned receipt + pre-registered cohort/endpoint contract and promote via the validation gate instead.

---

### Q3: Disease normalization — unknown disease handling ✅ **DEFAULT SET**
**Context:** Multiple disease normalizers exist today and at least one uses a risky fallback (unknown → ovarian).

**Default decision (so plumber can proceed safely):**
- **Remove ALL risky disease fallbacks** (no more `unknown → ovarian` or `unknown → mm` inside normalizers).
- Canonical normalizer returns `(is_valid, normalized_disease)` and for unknown returns: `(False, 'unknown')`.

**Behavior by surface (explicit + safe):**
- **Explicit resistance endpoints** (`/api/resistance/predict`, `/api/care/resistance_playbook`): if `is_valid=False` → return **422** with supported diseases list.
- **Orchestrators** (`/api/complete_care/v2`): if `is_valid=False` → **graceful degradation** (skip disease-specific modules like SOC/resistance; include a warning + provenance flag).

**Why:** This avoids silently producing ovarian outputs for non-ovarian patients, while still keeping the universal orchestrator usable for partial/unknown profiles.

**Implementation consequence:**
- Create one shared normalizer module and delete/stop using the other ad-hoc normalizers (`complete_care_universal/config.py`, `therapy_fit/config.py`, `extraction_agent.py`).
- Add `provenance.flags += ['invalid_disease_input']` and `normalized_disease='unknown'` when invalid.

**Manager override (only if desired):**
- If you want strict mode everywhere (including orchestrator), switch orchestrator behavior to 422 as well. Default keeps orchestrator graceful.

---

### Q4: Output schema backward compatibility ✅ **PARTIALLY ANSWERED**
**Context:** Unifying output schema may change response shapes for existing endpoints.

**Finding:** 
- `complete_care_universal.py` already uses `/v2` endpoint pattern (line 626: `@router.post("/v2")`)
- Ayesha orchestrator also uses `/complete_care_v2` endpoint
- **Versioning pattern exists** — suggests breaking changes are acceptable if versioned

**Question:** Do we need to maintain backward compatibility, or can we version the API?
- **Option A:** Maintain backward compatibility (add new fields, don't remove old ones) — **safest, but adds technical debt**
- **Option B:** Version endpoints (`/api/v2/resistance/predict`) and deprecate v1 with migration guide — **clean, follows existing pattern**
- **Option C:** Breaking change is acceptable if we document migration path (who uses these endpoints today?) — **risky, need to audit callers first**

**Recommendation:** Option B (version endpoints) — **follows existing pattern in codebase**. However, we should **audit who's calling `/api/care/resistance_playbook` and `/api/resistance/predict`** before making breaking changes. If no external callers, Option C is acceptable.

---

### Q5: Test dependencies (tooling) ✅ **ANSWERED**
**Context:** Ring-0 validation requires `pytest`, but it's not in `requirements.txt`.

**Finding:** `pytest>=7.4.0` and `pytest-asyncio>=0.21.0` are **already in `requirements.txt`** (lines 24-25). However, the validation runner script gracefully skips pytest if not installed (for environments without it).

**Question:** Should we move pytest to `requirements-dev.txt` to keep production lean?
- **Option A:** Move pytest to `requirements-dev.txt` (standard practice, keeps production image smaller)
- **Option B:** Keep in `requirements.txt` (simpler, but pytest is lightweight anyway)
- **Option C:** Keep current behavior (pytest in main requirements, runner skips if missing)

**Recommendation:** Option A (move to `requirements-dev.txt`) — cleaner separation, but current setup works. **Low priority** — this is already functional.

---

### Q6: Report stability (canonical vs timestamped) ✅ **PARTIALLY ANSWERED**
**Context:** Validators emit timestamped reports (`report_20251226_063623.json`), but we also need a stable "latest" pointer.

**Finding:** Current validators write to `out/<contract_id>/report.json` (e.g., `out/ddr_bin_ov_platinum_TRUE_SAE_v2/report.json`). **No timestamped files found** — reports are written directly to canonical paths. However, there's no mechanism to preserve historical runs or link to "latest" if we add timestamping later.

**Question:** Should we add timestamping + canonical linking, or keep current direct-write approach?
- **Option A:** Keep current (direct write to canonical path) — simple, but overwrites previous runs
- **Option B:** Add timestamping + copy-on-write: write `report_<timestamp>.json`, then copy to `report.json` (preserves history, works on Windows)
- **Option C:** Add versioned paths: `out/<contract_id>/v<version>/report.json` + symlink `latest.json` (explicit versioning, Unix-friendly)

**Recommendation:** Option B (copy-on-write) — preserves audit trail, works everywhere. **Medium priority** — current approach works but loses history.

---

### Q7: Orchestrator refactoring (internal HTTP calls) ✅ **ANSWERED**
**Context:** Orchestrators call `http://localhost:8000/api/insights/*` endpoints internally. This is brittle (requires server, not testable, slower).

**Finding:** 
- `complete_care_universal.py` has `_extract_insights_bundle()` that makes 4 HTTP calls to insights endpoints (lines 131-274)
- Pattern document (`.cursor/ayesha/Deliverables/Iterations/I8_DATA_FLOW.md`) states: **"Direct Import: When: Same backend process"** — this is the recommended pattern
- There's an `insights/bundle_client.py` package that could be used directly (found in `api/services/insights/README.md`)

**Question:** Refactoring strategy?
- **Option A:** Big-bang refactor: replace all HTTP calls with direct service imports (risky, but cleanest)
- **Option B:** Incremental: add direct service calls alongside HTTP, feature-flag, migrate one endpoint at a time
- **Option C:** Keep HTTP but add integration test harness that mocks the server (less ideal, but safer)

**Recommendation:** Option B (incremental) — safer, can ship value while refactoring. **Start by refactoring `_extract_insights_bundle()` to use `insights.bundle_client` directly instead of HTTP calls.**

---

### Q8: Data plumbing priority (OV expression vs MM CNV) ✅ **ANSWERED**
**Context:** Both are P1, but which unlocks more validation value first?

**Finding:**
- `Manager.mdc` explicitly states **"OV first (Ayesha track)"** as the strategy (line 5 in `RESISTANCE_VALIDATION_PLAN.md`)
- MFAP4 has **external validation** (GSE63885, AUROC 0.763) — stronger receipts than MM CNV
- MM resistance playbook already works at gene-level (DIS3, TP53 mutations validated) — CNV adds completeness but not new validation

**Question:** Which should we prioritize?
- **Option A:** OV expression (MFAP4/EMT rules) — needed for Ayesha track, has external validation (GSE63885)
- **Option B:** MM CNV/cytogenetics — needed for MM resistance playbook completeness, but validation exists at gene-level already

**Recommendation:** Option A (OV expression) — **confirmed by strategy docs**. Aligns with Ayesha-first strategy, has stronger validation receipts, unlocks production-ready OV resistance capability.

---

### Q9: Fixture scope for ResistanceProphet/ResistanceDetectionService
**Context:** Need fixtures to test L0/L1/L2 behavior, 2-of-3 triggers, confidence caps.

**Additional context:**
- `resistance_detection_service.py` implements "2-of-3 trigger" logic (HRD drop, DNA repair capacity drop, CA-125 inadequate response)
- `resistance_prophet_service.py` has L0/L1/L2 progressive enhancement with confidence caps
- Existing unit tests in `tests/test_resistance_playbook.py` cover some edge cases (HR restoration, ABCB1 upregulation, MAPK/PI3K activation)

**Question:** How comprehensive should fixtures be?
- **Option A:** Minimal (3–5 fixtures covering happy path + 1 edge case) — **fast to implement, can expand**
- **Option B:** Comprehensive (10–15 fixtures covering all edge cases + regression cases) — **prevents regressions, but slower**
- **Option C:** Start minimal, expand as bugs are found (iterative) — **pragmatic, but risks missing edge cases**

**Recommendation:** Option C (iterative) — **start with Option A** (happy path + L0/L1/L2 + 2-of-3 trigger edge cases), then expand as we find bugs or add new features. **Balance between coverage and velocity.**

---

### Q10: Provenance versioning strategy
**Context:** Need to track `code_version`, `contract_version`, `flags`, `inputs_snapshot_hash`.

**Additional context:**
- `complete_care_universal.py` already generates `run_id` with timestamp: `f"complete_care_universal_v2_{int(datetime.utcnow().timestamp())}"` (line 1020)
- No explicit contract versioning found in current code
- `Testing.mdc` suggests versioned outputs with metadata (commit hash, python version, data checksums) — **pattern exists**

**Question:** How should we version contracts?
- **Option A:** Semantic versioning (`contract_version: "1.0.0"`) — explicit, but requires discipline
- **Option B:** Git commit hash (`code_version: "abc123"`) — automatic, but not human-readable
- **Option C:** Timestamp + hash (`contract_version: "2025-01-26_abc123"`) — both human-readable and unique

**Recommendation:** Option C (timestamp + hash) — **best of both worlds, matches existing `run_id` pattern**. However, we should also include:
- `code_version`: git commit hash (automatic, for reproducibility)
- `contract_version`: semantic version (manual, for breaking changes)
- `run_id`: timestamp-based (for traceability)
- `inputs_snapshot_hash`: SHA256 of normalized inputs (for cache invalidation)

---

## References (canonical docs)
- Validation strategy: `.cursor/MOAT/RESISTANCE_VALIDATION_PLAN.md`
- Manager blueprint: `.cursor/MOAT/SAE_INTELLIGENCE/Manager.mdc`
- `Ayesha`/MOAT state of build: `.cursor/ayesha/MOAT/COMPREHENSIVE_AYESHA_ARCHIVE_AUDIT.md`


