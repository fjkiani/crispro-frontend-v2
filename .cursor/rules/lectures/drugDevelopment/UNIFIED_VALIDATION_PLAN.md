# 🎯 Unified Validation Plan: Evidence-First (Claim → Receipt → Next Receipt)

**Purpose:** Establish a single source of truth for what CrisPRO has **actually validated** (with reproducible receipts) vs what is **not yet validated** (and must not be claimed).

**Scope (tight):** This plan focuses only on validated claims across:
- `publications/sporadic_cancer/`
- `publications/synthetic_lethality/`
- `publications/01-metastasis-interception/`

We explicitly avoid noise from other folders unless a claim is backed by a reproducible artifact.

---

## Definitions (non-negotiable)

- **Validated**: backed by a reproducible artifact (saved report/receipt) and a runnable script/command.
- **Benchmarked (demo-safe)**: deterministic sanity checks (synthetic tests, offline validators). Useful, but **not** outcome evidence.
- **Projected**: ROI / business impact numbers. Allowed only if explicitly labeled as projections and assumptions are listed.

---

## ✅ Validated Claims Inventory (with receipts)

### 1) Sporadic Cancer — deterministic tumor-context gates + provenance (AACR package)

**Validated (deterministic behavior, not outcomes):**
- Sporadic gate logic passes deterministic validation (6/6).
- E2E workflow produces provenance-bearing outputs.

**Primary receipts:**
- `publications/sporadic_cancer/submission_aacr/receipts/validate_sporadic_gates.txt` (6/6 PASS)
- `publications/sporadic_cancer/submission_aacr/receipts/validate_sporadic_gates_report.json` (status PASS)
- `publications/sporadic_cancer/submission_aacr/SUPPLEMENT.md` (claim → code → receipt table)
- `publications/sporadic_cancer/submission_aacr/MANUSCRIPT_DRAFT.md` (results section enumerates receipts)

**Reproduce:**
- `oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_sporadic_gates.py` (writes `.../out/sporadic_gates/report.json`)

**Explicit non-claims (forbidden here):
- clinical outcome benefit (OS/PFS/response) from sporadic gates.

---

### 2) Synthetic Lethality — publication suite + benchmark composition

**Validated (in this publication bundle):**
- Benchmark composition is explicitly tracked.
- A publication suite runner produces saved results packs.

**Primary receipts:**
- `publications/synthetic_lethality/docs/results_pack.md`
- `publications/synthetic_lethality/docs/baseline_comparison.md` (e.g., SL-positive n=70, SL-negative n=30)
- `publications/synthetic_lethality/results/publication_suite_20251230_*.md` (full run receipts)

**Reproduce:**
- `publications/synthetic_lethality/code/run_publication_suite.py` (writes JSON + MD receipts)
- `publications/synthetic_lethality/code/validate_test_cases.py` (writes `results/validation_report.json`)

**Explicit non-claims:**
- Any survival/response improvement claims without outcome-labeled cohorts.

---

### 3) Metastasis Interception — Target-Lock validation package

**Validated (in this publication bundle):**
- Per-step AUROC with bootstrap CIs on a defined dataset/seed.
- Precision@K, enrichment tests, effect sizes, and calibration curves are reported.

**Primary receipts:**
- `publications/01-metastasis-interception/figures/LEGENDS.md`
  - Dataset: n=38 primary genes, 8 steps (304 points), 50 positive labels
  - Bootstrap: B=1000, seed=42
  - Mean AUROC = 0.976 ± 0.035
  - Macro P@3 = 1.000 ± 0.000
  - RUO: chromatin is a deterministic stub (Enformer-ready)
- `publications/01-metastasis-interception/REPRODUCIBILITY.md`

**Explicit non-claims:**
- Claims that the chromatin component is validated (it is explicitly stubbed RUO).

---

### 4) PGx Dosing Guidance — safety foundation (used in platform story)

**Validated (strongest clinical-adjacent receipts we have):**
- Cohort exists and is curated.
- Sensitivity on known toxicity cases is explicitly reported.
- CPIC mapping/logic is documented.

**Primary receipts:**
- `oncology-coPilot/oncology-backend-minimal/dosing_guidance_validation/docs/SME_REVIEW_PACKAGE.md (N=59; sensitivity 6/6 flagged)
- `oncology-coPilot/oncology-backend-minimal/dosing_guidance_validation/docs/CPIC_ALIGNMENT_SUMMARY.md`

**Explicit non-claims:**
- “MedWatch reduction” or real-world adverse event reduction without outcome receipts.

---

## 🚫 Forbidden Claims (until outcome receipts exist)

Do **not** claim any of the following without an outcome-labeled cohort + statistical analysis receipts:
- “Escape from the 71% Phase 2 failure zone”
- “Improves Phase 2 success rate to X%”
- “95% MedWatch reduction”
- “Responder identification = 96.6%” (unless backed by a published, reproducible evaluation artifact)
- “Top-5 accuracy 17/17 patients” (unless backed by a published, reproducible evaluation artifact)

---

## Validation Roadmap (receipt-driven, build-last)

### Tier 0 (Immediate): Claim hygiene + receipts ledger

**Deliverable (must exist before anything else):**
- `VALIDATED_CLAIMS_LEDGER.md` (one table: claim → dataset → script → output paths → reproduce command)

Acceptance criteria:
- Every “validated” line has an on-disk receipt.
- Every “projected” line is explicitly labeled projected.

---

### Tier 1 (Next): Cross-capability *measurement* on existing cases (no new product code)

Goal: demonstrate we can compute multiple signals **on the same cases** and store receipts.

**Inputs:** Start with the strongest curated cohorts we already have:
- dosing guidance cohort (PGx)
- sporadic gate scenario suite receipts
- synthetic lethality benchmark suite

**Outputs (must be saved):**
- `reports/cross_capability_metrics.json`
- `reports/cross_capability_metrics.md`

Must include:
- coverage/missingness table per feature
- per-case provenance snapshots

Acceptance criteria:
- deterministic rerun: same inputs → same outputs
- coverage quantified (not implied)

---

### Tier 2 (Then): Define unified score *evaluation spec* before implementing a unified score

Use the critique in `ec2.mdc`:
- define what “risk–benefit” means (utility function, veto rules)
- define calibration expectations (e.g., ECE target)) for any probability-like outputs

Deliverable:
- `UNIFIED_SCORE_EVALUATION_SPEC.md` (acceptance tests + datasets required)

---

### Tier 3 (Outcome validation): only after Tier 1–2 receipts exist

- Pick one outcome-labeled dataset per claim type (efficacy vs safety vs enrollment).
- Pre-register the analysis plan inside the repo (thresholds, endpoints, stats).

Deliverable:
- `OUTCOME_VALIDATION_PROTOCOL.md`

---

## Bottom line

We already have three publication bundles with real receipts. This plan forces the platform narrative to be earned by:
1) freezing validated claims,
2) forbidding outcome/business claims until outcome receipts exist, and
3) requiring each next step to generate new receipts before any “platform” code is written.
