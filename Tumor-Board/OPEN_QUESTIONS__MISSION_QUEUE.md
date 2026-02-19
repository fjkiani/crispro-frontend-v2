---
title: "OPEN QUESTIONS — Mission Queue (Receipts Required)"
ruo: true
contract:
  - "This file is the ONLY place to add new mission orders."
  - "One mission must produce one proof artifact (receipt) or a code+data snippet."
  - "Every acceptance criterion must be satisfiable by: (a) JSON field path, (b) file path + line numbers, or (c) dataset snippet."
  - "Do not write re-run/code-derived facts into the bundle-only packet."
related_docs:
  - "TUMOR_BOARD_PACKET__BUNDLE_ONLY.md"
  - "ENGINE_EXPLANATIONS__WHY_THE_BUNDLE_LOOKS_LIKE_THIS.md"
---

## 0) How to use this queue

Each mission below is a “tumor-board justification clarifier” or an “engine hardening truth test.”

Template:
- **Mission ID**
- **Question**
- **Owner agent**
- **Inputs**
- **Exact artifact to produce**
- **Acceptance criteria**
- **Status**
- **Notes / next action**

## 1) Completed missions (manager’s high-yield set)

### M4 — DepMap “terrain bonus” verification
- **Question**: What exact DepMap lineage numbers and rule yield the “+0.15 confidence boost” language for ATR/WEE1 essentiality?
- **Owner agent**: Completed (Memento)
- **Inputs**: `publications/synthetic_lethality/data/depmap_essentiality_by_context.json`
- **Artifact produced**: Included in `ENGINE_EXPLANATIONS__WHY_THE_BUNDLE_LOOKS_LIKE_THIS.md` (Mission 4 section)
- **Acceptance**: Must show ATR/WEE1/PARP1 mean effects for Ovary/Fallopian + cite the boost rule location in code.
- **Status**: ✅ Completed

### M3 — SL confidence 0.85 decomposition
- **Question**: How is `recommended_drugs[].confidence = 0.85` computed?
- **Owner agent**: Completed (Memento)
- **Inputs**: SL drug catalog + recommender logic
- **Artifact produced**: Included in `ENGINE_EXPLANATIONS__WHY_THE_BUNDLE_LOOKS_LIKE_THIS.md` (Mission 3 section)
- **Acceptance**: Must decompose into base + boosts and cite the line(s) producing the field.
- **Status**: ✅ Completed

### M1 — Pathway alignment trace for efficacy
- **Question**: Why do paclitaxel/adavosertib get pathway percentile 1.0 while olaparib gets 0.0 (in this payload family)?
- **Owner agent**: Completed (Memento)
- **Inputs**: Pathway scoring + ovarian panel drug weights
- **Artifact produced**: Included in `ENGINE_EXPLANATIONS__WHY_THE_BUNDLE_LOOKS_LIKE_THIS.md` (Mission 1 section)
- **Acceptance**: Must cite code path for `path_pct`, the drug weights source, and the pathway aggregation source.
- **Status**: ✅ Completed

### M5 — MBD4 asymmetry: SL normalizes indel but efficacy excludes it
- **Question**: Is the asymmetry intentional, or a gap?
- **Owner agent**: Completed (Memento) — gap identified
- **Inputs**: WIWFM sequence processor + SL HGVS resolution path
- **Artifact produced**: Included in `ENGINE_EXPLANATIONS__WHY_THE_BUNDLE_LOOKS_LIKE_THIS.md` (Mission 5 section)
- **Acceptance**: Must cite both code paths and clearly state the gap/decision.
- **Status**: ✅ Completed (Gap logged)

### M2 — Evidence engine autopsy (why zero citations)
- **Question**: Are citations zero because evidence was skipped, or because providers returned zero?
- **Owner agent**: Completed (Memento) — broken dependency
- **Inputs**: `/api/evidence/literature` runtime + evidence router code
- **Artifact produced**: `artifacts/mission2_evidence_autopsy_2026-02-13.json`
- **Acceptance**: Must show runtime response flags + error and cite code path.
- **Status**: ✅ Completed (BROKEN: `genai` import error)

### M6 — MFAP4 + “Diamond SAE feature 27607” wiring audit
- **Question**: Are MFAP4/27607 wired into the Ayesha Therapy-Fit bundle by default?
- **Owner agent**: Completed (Memento)
- **Inputs**: Resistance Prophet constants + ovarian detectors + therapy-fit service
- **Artifact produced**: Included in `ENGINE_EXPLANATIONS__WHY_THE_BUNDLE_LOOKS_LIKE_THIS.md` (Mission 6 section)
- **Acceptance**: Must cite definition + detector + orchestration and answer wiring status.
- **Status**: ✅ Completed

## 2) Open missions (next)

### O1 — Fix Evidence Engine dependency (restore citations deterministically)
- **Question**: What is the minimal code change to make `/api/evidence/literature` return deterministic PMIDs again in this environment?
- **Owner agent**: Unassigned
- **Inputs**: `oncology-coPilot/oncology-backend-minimal/api/routers/evidence/literature.py`, requirements/env, `artifacts/mission2_evidence_autopsy_2026-02-13.json`
- **Artifact to produce**:
  - `artifacts/receipt__evidence_literature_fixed_<date>.json` showing non-empty `pmids` and `top_results_count>0` for a stable query
  - Patch with code+line anchors
- **Acceptance criteria**:
  - Must remove/replace the failing `genai` import path or guard it cleanly
  - Must show runtime receipt with `agent_fallback=false` OR successful fallback to E-utilities returning results
- **Status**: ⏳ Pending

### O2 — Unify indel normalization across SL and WIWFM (MBD4 parity)
- **Question**: Should WIWFM accept the same normalized allele-complete variant objects that SL uses (HGVS→GRCh38 + left-padding), and if yes, what exact adapter layer change is required?
- **Owner agent**: Unassigned
- **Inputs**: `sequence_processor.py` exclusion logic, SL normalization helpers
- **Artifact to produce**:
  - A before/after bundle (or minimal receipt) showing MBD4 included in WIWFM `variants_sent_to_engine` rather than excluded
  - Patch with code+line anchors
- **Acceptance**:
  - Must demonstrate that a non-SNV indel no longer becomes “unscored by design” for WIWFM
  - Must include safeguards to avoid invalid VCF representations
- **Status**: ⏳ Pending

### O3 — Decide where Resistance Prophet signals belong in Therapy-Fit bundle (non-debug)
- **Question**: Should MFAP4/DIAMOND features appear as a separate section of the Therapy-Fit bundle by default (RUO), and what minimal safe contract should be returned?
- **Owner agent**: Unassigned
- **Inputs**: `resistance_prophet_service.py`, `ayesha_fit/service.py`, therapy-fit router
- **Artifact to produce**:
  - Proposed schema snippet + one runtime receipt showing the new field is present (with empty-safe defaults when missing expression data)
- **Acceptance**:
  - Must not change existing top-level keys in a breaking way
  - Must be explicit about “Unknown” when expression/SAE features absent
- **Status**: ⏳ Pending

