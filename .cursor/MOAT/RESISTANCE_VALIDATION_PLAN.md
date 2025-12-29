# 🧪 Resistance + Cancer-Weakness Validation Plan (MOAT-wide)

**Canonical manager blueprint (why/what):** `.cursor/MOAT/SAE_INTELLIGENCE/Manager.mdc`  
**Canonical execution plan (how/when):** **this file**  
**Canonical weakness-exploitation pillar:** `.cursor/ayesha/SYNTHETIC_LETHALITY_COMPLETE.md` + `oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/`

**Goal:** validate a *production MOAT capability* for cancer resistance **and** cancer weakness exploitation as a **system** (not just SAE), with pinned contracts, deterministic validators, and CI gates.

**Scope order:** OV first (Ayesha track) → MM → cancer-agnostic onboarding via the same contracts.

---

## What “validated” means (promotion gate)

We only call something **validated for endpoint X** when all are true:

1. **Endpoint contract is pinned** (labels + inclusion/exclusion + timepoint)
2. **Inputs contract is pinned** (required vs optional signals; L0/L1/L2)
3. **Deterministic validator exists** (script path pinned)
4. **Validator emits `report.json`** (metrics + CI + cohort metadata)
5. **Confound checks computed** (at minimum: coverage/variant-count sensitivity when relevant)

Everything else must be shown as **Tier 4–5 evidence** with explicit caveats.

---

## Current code surfaces (what exists today)

We keep multiple entrypoints, but standardize contracts across them:

- **Resistance playbook (action utility, rule-based)**
  - Endpoints:
    - Legacy: `POST /api/care/resistance_playbook` (compat)
    - Strict: `POST /api/care/resistance_playbook_v2` (explicit disease)
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/care.py`
  - Logic: `oncology-coPilot/oncology-backend-minimal/api/services/resistance_playbook_service.py`

- **Resistance prediction (disease-agnostic wrapper)**
  - Endpoint: `POST /api/resistance/predict`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/resistance.py`
  - Calls: `ResistanceProphetService` and/or `ResistancePlaybookService`

- **Orchestrators (product surface)**
  - Universal: `POST /api/complete_care/v2` and alias `/api/complete_care/universal`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/complete_care_universal.py`

- **Cancer weakness exploitation (Synthetic Lethality / Essentiality)**
  - Service: `oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/`
  - Canonical doc: `.cursor/ayesha/SYNTHETIC_LETHALITY_COMPLETE.md`
  - Pilot benchmark: **50% drug match, 100% Evo2 usage** (in doc)

---

## Contracts (pin these first)

### Input completeness contract (progressive enhancement)
We validate (and cap claims) by input completeness:

- **L0**: no tumor NGS (disease/stage + baseline labs/clinical)
- **L1**: partial tumor metrics (e.g., HRD score or TMB only)
- **L2**: full tumor NGS (mutations + CNV + key scores; optional expression; optional serial markers)

### Endpoint contracts (initial pinned set)

- **`ov_platinum_response_tcga_style`**
  - Label: resistant/refractory vs sensitive
  - Timepoint: baseline
  - Receipt: DDR_bin is ~random under this contract
    - `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_ov_platinum_TRUE_SAE_v2/report.json`

- **`ov_early_progression_tcga`**
  - Label: early progression (per validator)
  - Receipt: EMT + canonical HRD validators exist
    - `oncology-coPilot/oncology-backend-minimal/scripts/validation/test_emt_hrd_early_progression_tcga_ov.py`

- **`ov_platinum_resistance_expression_external`**
  - Label: platinum resistant vs sensitive (external)
  - Receipt: MFAP4 story in publication
    - `.cursor/MOAT/SAE_INTELLIGENCE/Publication-1/SAE_RESISTANCE/`

- **`sl_pilot_drug_match_v1`** (synthetic lethality)
  - Label: binary match/no-match vs expected drug class/therapy in curated pilot cases
  - Receipt: pilot benchmark results + corrected benchmark script
    - `.cursor/ayesha/SYNTHETIC_LETHALITY_COMPLETE.md`

- **`mm_drug_class_resistance`** (planned)
  - Labels: per-line PI / IMiD / anti-CD38 response
  - Inputs: mutations + CNV/cytogenetics + treatment line

---

## Validation tracks (MOAT-wide)

### Track A — Predictive resistance (baseline)
- Metrics: AUROC + AUPRC; calibration (ECE)
- Acceptance (V1 “predictive”): AUROC ≥ 0.65 with CI not crossing 0.55 on ≥1 external cohort (or pre-registered internal cohort with strict confounds)

### Track B — Monitoring & early warning (on-therapy)
- Metrics: lead-time gain (weeks), sensitivity @ fixed FPR
- Acceptance (V1 “early warning”): deterministic correctness on fixtures + stable FP behavior; cohort evaluation when time-series available

### Track C — Action utility (recommendations)
- Metrics: deterministic correctness + contraindication rate ~0 + provenance completeness

### Track D — Safety (toxicity/PGx)
- Metrics: fixture correctness + “no missed high-risk” invariant; RUO until outcomes validation

### Track E — Cancer weakness exploitation (Synthetic Lethality)
- Metrics: drug match accuracy (binary + Top-N once API supports ranked list), essentiality correlation targets (DepMap proxy), Evo2 usage rate
- Acceptance (V1): deterministic benchmark runner passes + reports metrics; “rules-only bypass” must be measured and controlled (see SL doc)

---

## Sprint plan (deliverables in sprints)

We run 1-week sprints for plumbing + validation deliverables.

### Sprint 0 (DONE) — Contract foundation + safety guardrails
**Goal:** remove drift sources (disease defaults, evidence tier contradictions, localhost coupling).

**Delivered (code receipts):**
- Shared disease normalization (no risky defaults)
  - `api/services/disease_normalization.py`
- Resistance evidence tier mapping (manager Tier 1–5)
  - `api/services/resistance_evidence_tiers.py`
- Canonical resistance contract models + builder
  - `api/contracts/resistance_contract.py`
  - `api/contracts/resistance_builders.py`
- `/api/resistance/predict` emits additive canonical `contract`
- Orchestrator no longer calls localhost for resistance playbook (direct service call)
- Added strict `POST /api/care/resistance_playbook_v2` and stabilized legacy endpoint

**Acceptance:** Ring-1 suite passes (NF1 validator + receipt parse).

---

### Sprint 1 (NOW) — Unify product surface outputs + provenance scoreboard
**Goal:** UI + orchestrator consume **one canonical resistance object** consistently.

**Deliverables:**
- Orchestrator: emit `results["resistance_contract"]` (single canonical object) and deprecate ad-hoc fields over time
- Add deterministic contract validator (Ring-1): asserts required keys + provenance flags + evidence tier presence
- Add weekly scoreboard stub in `report.json` (per Track A–E)

**Acceptance:**
- Ring-1 includes new contract validator that passes on pinned fixtures
- Orchestrator response includes canonical `resistance_contract` with provenance + evidence tiers

---

### Sprint 2 — OV MVP hardening (Ayesha) ✅ DONE
**Goal:** ship OV resistance MVP that is honest and measurable.

**Deliverables shipped:**
- ✅ **L0/L1/L2 caps + missing-input warnings**
  - Logic: `oncology-coPilot/oncology-backend-minimal/api/services/input_completeness.py`
  - Wired into:
    - `POST /api/resistance/predict` (`api/routers/resistance.py`)
    - `POST /api/complete_care/v2` (and `/api/complete_care/universal`) (`api/routers/complete_care_universal.py`)
- ✅ **OV expression ingestion interface (RUO)**
  - Endpoint: `POST /api/expression/ingest`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/expression.py`
  - Storage path: `oncology-coPilot/oncology-backend-minimal/data/expression/`
- ✅ **Serial marker ingestion contract for CA-125 (time-series)**
  - Plumbed through universal orchestrator → ResistanceProphet via `biomarker_history.ca125_history`
  - Router: `oncology-coPilot/oncology-backend-minimal/api/routers/complete_care_universal.py`

**Acceptance (verified):**
- ✅ Ring‑1 validation scripts pass
  - Re-run command:
    - `python3 oncology-coPilot/oncology-backend-minimal/scripts/validation/run_resistance_validation_suite.py --ring1`
- ✅ Contract shows missing inputs explicitly and enforces caps (L0=0.4, L1=0.6, L2=0.8)

---

### Sprint 3 — Synthetic Lethality: validation harness + contract alignment
**Goal:** treat SL as a first-class “weakness exploitation” pillar with the same rigor.

**Deliverables:**
- Add a deterministic SL benchmark runner that emits `report.json` for `sl_pilot_drug_match_v1`
- Add a contract mapping so SL output can be shown as MOAT actions (treatments/tests) with evidence tiers
- Add a guardrail check to surface GUIDANCE_FAST / rules-only bypass status in reports

**Acceptance:**
- `sl_pilot_drug_match_v1/report.json` exists and is stable
- Evo2 usage rate is reported explicitly

**Sprint 3 implementation pointers (shipped):**
- Ring‑1 validator (deterministic; no network):
  - `oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_synthetic_lethality_pilot_benchmark.py`
  - Output: `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/sl_pilot_drug_match_v1/report.json`
- Pinned pilot receipt source (Evo2-backed benchmark artifact):
  - `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/results/benchmark_efficacy_20251203_210907.json`
- Contract alignment helper (SL → ActionCards with manager-facing evidence tiers):
  - `oncology-coPilot/oncology-backend-minimal/api/contracts/resistance_builders.py::action_cards_from_synthetic_lethality_result`

---

### Sprint 4 — MM resistance (drug-class endpoints)
**Goal:** pin MM contracts and produce first receipts per drug class.

**Deliverables:**
- MM endpoint contract(s): PI/IMiD/anti-CD38
- Minimal dataset pinned (even small) + first `report.json`
- CNV/cytogenetics input schema + fixtures

**Acceptance:**
- At least one MM drug-class endpoint has pinned cohort + deterministic validator

---

## Registry (what we’ve pinned)

- **Ring-1 suite runner:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/run_resistance_validation_suite.py`
- **NF1 action-utility validation script:** `.../scripts/validation/validate_ov_nf1_playbook.py`
- **DDR_bin receipts (anti-overclaim):** `.../scripts/validation/out/ddr_bin_ov_platinum_TRUE_SAE_v2/report.json`
- **Synthetic lethality canonical doc:** `.cursor/ayesha/SYNTHETIC_LETHALITY_COMPLETE.md`
- **Synthetic lethality pilot receipt (Ring‑1):** `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/sl_pilot_drug_match_v1/report.json`

---

## How this avoids "flying blind"

- Every sprint deliverable must ship with:
  - pinned contract,
  - deterministic validation script,
  - stable `report.json`,
  - provenance fields + flags,
  - and explicit evidence tiering.

That’s how we make the MOAT system safer, more reproducible, and incrementally closer to real patient impact.
