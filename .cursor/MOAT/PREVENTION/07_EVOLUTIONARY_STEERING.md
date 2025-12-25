# 07 — Evolutionary Steering (Therapy Cycling)

**Purpose:** Encode the “PARP → ATR → PARP” cycling concept as a guarded, RUO-only advanced option.

---

## Core Idea

Switch therapies to kill what the previous therapy selected for.

Example hypothesis:
- PARP selects for HR restoration
- HR-restored cells may become ATR-dependent
- ATR inhibitor can reduce HR-restored population
- return to PARP when DDR_bin rebounds

---

## Inputs

- Longitudinal DDR_bin (≥2–3 points)
- Therapy history
- Tolerance constraints

---

## Guardrails

- RUO-only language
- Trial-first if ATR inhibitor not standard
- Clear “unknowns” surfaced (failure modes)

---

## Outputs

- Cycle plan template:
  - Switch condition (DDR_bin < T1)
  - Return condition (DDR_bin > T2)
  - Monitoring cadence

---

## Manager Decisions Needed

- Are we allowed to mention ATR cycling in demo?
- If yes, what wording (options-only vs recommendation)?

