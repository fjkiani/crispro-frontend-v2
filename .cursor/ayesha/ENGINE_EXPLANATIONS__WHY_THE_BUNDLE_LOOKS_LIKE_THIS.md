---
title: "ENGINE EXPLANATIONS — Why the Bundle Looks Like This (Engineering/RUO)"
ruo: true
contract:
  - "This document is engineering: code paths, flags, equations, and re-run receipts."
  - "Nothing here should be quoted as payload-derived unless it cites a bundle JSON path."
  - "All runtime claims must cite an artifact path (receipt) and/or a code file path + line numbers."
related_quote_safe_packet:
  - ".cursor/ayesha/TUMOR_BOARD_PACKET__BUNDLE_ONLY.md"
mission_queue:
  - ".cursor/ayesha/OPEN_QUESTIONS__MISSION_QUEUE.md"
---

## 0) Purpose

The tumor board packet must remain quote-safe and payload-bounded. This file exists to answer the separate engineering question:

> “Why does the JSON look like this, and which backend levers produced it?”

This file may include:
- re-run receipts
- backend code anchors
- controlled ablations
- known drift between payload snapshots

It must **not** be merged into the bundle-only packet.

## 1) Drift exists → we version artifacts

The bundle-only packet is bounded to the captured payload (see `artifacts/payload.extracted.json` referenced in `.cursor/ayesha/TUMOR_BOARD_PACKET__BUNDLE_ONLY.md`).

Separately, we have observed “payload drift” when re-running endpoints later (e.g., variant exclusion reasons and sequence scoring counts can change as normalization logic changes). That is expected during active development and is why this file exists as an engineering companion rather than mixing with tumor-board quoting material.

## 2) Dominant backend levers (code anchors)

### 2.1 Therapy-Fit L1 forces fast mode and disables evidence in the default path

In the Therapy-Fit service, L1 request options are set explicitly:

```125:141:oncology-coPilot/oncology-backend-minimal/api/services/ayesha_fit/service.py
options={
    "include_low_confidence": True,
    "fast": True,
    "allow_heuristic_sequence": True,
    "include_evidence": False,
    "include_insights": True,
}
```

Engineering implication:
- A bundle can show `evidence_enabled=true` in flags yet still produce citation-empty drug rows if `include_evidence` is forced off in the caller layer, or if evidence-gathering is gated by fast-mode.

### 2.2 Evidence gathering is gated off in fast mode

Evidence gathering eligibility is computed as:

```37:40:oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/helpers/evidence_gatherer.py
gather_evidence = (not fast_mode) and evidence_enabled_flag and primary_gene and primary_variant.get("hgvs_p")
```

Engineering implication:
- If `fast=true`, evidence collection is **not executed**, even if the run flags say “evidence enabled.”

### 2.3 Efficacy score formula (why “bands” can exist)

Raw efficacy (likelihood-of-benefit) uses:

```286:297:oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py
raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd + clinvar_prior
if tier != "insufficient":
    lob = raw_lob
else:
    # Default behavior is to suppress efficacy to 0 for insufficient evidence (clinical-safe),
    # but research/validation runs may request a continuous scalar even for low-confidence drugs.
    if getattr(confidence_config, "allow_insufficient_efficacy", False):
        lob = raw_lob * 0.5
    else:
        lob = 0.0
```

Engineering implication:
- “Score bands” (e.g., 0.32 vs 0.12) are often the formula behaving consistently given a small number of non-zero terms (e.g., some drugs get `path_pct=1.0`, others `path_pct=0.0`, while evidence stays at 0.0).

## 3) Mission receipts (completed) — the “why” proofs

These are the manager’s highest-yield clarifications. Each mission must have its own proof artifact (receipt) or code+data snippet.

### Mission 4 — DepMap “terrain bonus” verification (COMPLETED)

**Proof target:** Show the Ovary/Fallopian lineage inputs and the numeric rule that yields the “+0.15 confidence boost” language in SL `essential_pathways`.

- **Dataset**: `publications/synthetic_lethality/data/depmap_essentiality_by_context.json`
- **Rule (code)**: DepMap boost injected when `depmap_mean_effect < -0.5`, granting +0.15.
  - Code anchor: `oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/dependency_identifier.py` (depmap boost logic; see scribe mission section in legacy file for exact lines if needed)

**Receipt-backed facts (from dataset)**:
- ATR, Ovary/Fallopian: `depmap_mean_effect = -1.05`
- WEE1, Ovary/Fallopian: `depmap_mean_effect = -2.77`
- PARP1, Ovary/Fallopian: `depmap_mean_effect = -0.21`

Interpretation (engineering-only):
- ATR and WEE1 qualify for the +0.15 boost; PARP1 does not.

### Mission 3 — SL confidence 0.85 decomposition (COMPLETED)

**Proof target:** Show exactly how SL `recommended_drugs[].confidence = 0.85` is computed (base + boosts).

- Code anchor: `oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/drug_recommender.py` (confidence composition)
- Drug catalog anchor: `oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/constants.py` (drug indications, pathways)

Engineering decomposition for Ceralasertib (ATR) in ovarian:
- base = 0.40
- ovarian indication boost = +0.20
- pathway alignment (ATR present) = +0.10
- DepMap terrain boost (ATR essential) = +0.15
- Total = 0.85

### Mission 1 — Pathway alignment trace for efficacy (COMPLETED)

**Proof target:** Explain (with trace) why paclitaxel/adavosertib get pathway percentile 1.0 while olaparib gets 0.0 in this payload family.

Code anchors:
- `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py` (`path_pct` computation)
- `oncology-coPilot/oncology-backend-minimal/api/services/pathway/drug_mapping.py` (`get_pathway_weights_for_drug`)
- `oncology-coPilot/oncology-backend-minimal/api/services/pathway/panel_config.py` (ovarian drug panel weights)
- `oncology-coPilot/oncology-backend-minimal/api/services/pathway/aggregation.py` (pathway aggregation + hotspot lifts)

Engineering explanation (high level):
- If the run’s `pathway_scores` contains a pathway a drug heavily weights, `path_pct` rises; if not, `path_pct` can be 0.

### Mission 5 — MBD4 asymmetry fix/decision (COMPLETED → GAP IDENTIFIED)

**Proof target:** Prove why SL can normalize the MBD4 indel but efficacy excludes it (intentional vs gap).

Code anchors:
- WIWFM exclusion path: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sequence_processor.py` (excludes non-SNVs with `non_snv_requires_normalization`)
- SL normalization path: `oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/sl_agent.py` (HGVS→GRCh38 resolution + left-pad single-base deletions)

Engineering conclusion:
- This is a **GAP_TO_FIX** if the product goal is “same normalization policy across SL and WIWFM.”

### Mission 2 — Evidence Engine autopsy (COMPLETED → BROKEN)

**Proof target:** Explain why `citations_count=0` and why PubMed queries yield zero.

- Receipt: `artifacts/mission2_evidence_autopsy_2026-02-13.json`
- Code anchor: `oncology-coPilot/oncology-backend-minimal/api/routers/evidence/literature.py`

Engineering verdict:
- The literature path returns `agent_fallback=true` due to an import error: `cannot import name 'genai' from 'google' (unknown location)` (see receipt), resulting in `top_results_count=0` and no citations.

### Mission 6 — MFAP4 & “Diamond SAE feature 27607” wiring (COMPLETED)

**Proof target:** Show where these are defined and whether they are wired into Therapy-Fit by default.

Code anchors:
- Definitions: `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet/constants.py`
- Detectors: `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet/signals/ovarian.py`
- Orchestration: `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet_service.py`
- Therapy-Fit default path: `oncology-coPilot/oncology-backend-minimal/api/services/ayesha_fit/service.py` and `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py`

Engineering conclusion:
- MFAP4 + DIAMOND features are wired into Resistance Prophet ovarian signals, but are **not** included in the Therapy-Fit bundle by default (only present behind debug-only paths).

