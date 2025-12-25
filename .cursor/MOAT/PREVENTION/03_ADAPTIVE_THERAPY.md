# 03 — Adaptive Therapy (DDR_bin-Guided)

**Purpose:** Operationalize “dose modulation” as a deterministic playbook.  
**Status:** Spec-only (requires Manager thresholds).

---

## What Adaptive Therapy Is

Adaptive therapy is **intentional dose reduction / holiday** to reduce selection pressure so sensitive cells suppress resistant clones (fitness cost).

---

## Inputs

- Baseline DDR_bin
- Latest DDR_bin
- Trend slope (2+ datapoints)
- Current regimen class (PARP, platinum, ATR, etc.)
- Tolerance (toxicity grade, labs) — hard gate

---

## Rule (Template — Manager must lock)

### Trigger
- If ΔDDR_bin ≤ **T_WARN** (from baseline): enter adaptive mode

### Action
- Reduce dose by **X%** OR hold for **Y days**

### Monitoring
- Re-check DDR_bin in **Z weeks**

### Success criterion
- DDR_bin rebounds by ≥ **R%** OR slope turns positive

### Failure criterion
- DDR_bin continues falling OR CA-125 rises sharply OR clinical progression

---

## Safety Guardrails (Non-negotiable)

- Do not propose adaptive therapy if:
  - Grade ≥3 toxicity
  - significant cytopenias / renal impairment
  - provider override

---

## Outputs

- Action: “Consider reducing PARP dose by X%”
- Rationale: “selection pressure reduction”
- Evidence tier: **B/C** (depends on indication)
- Provenance + RUO

