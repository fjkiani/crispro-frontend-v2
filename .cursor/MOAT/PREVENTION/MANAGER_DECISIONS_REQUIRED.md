# �� Manager Decisions Required (Prevention Build)

**Date:** December 24, 2025  
**Audience:** Manager (decision-maker), Plumber (builder)  
**Purpose:** Convert PREVENTION vision into a buildable spec by locking assumptions.

---

## How to Answer
For each item, please provide:
- **Decision** (pick one / provide numbers)
- **Confidence**: HIGH / MEDIUM / LOW
- **Rationale + Source**: (guideline, paper, clinical precedent, expert judgment)
- **If unknown**: approve the **default** + tell us what experiment would de-risk it

---

## A) Claims We Are Allowed To Make (Demo / Product)

### **A1. What can we say publicly about DDR_bin today?**
Pick ONE:
- **(1) Observational only**: “DDR_bin is a research biomarker that correlates with resistance risk.”
- **(2) Early-warning claim (conservative)**: “May flag resistance earlier than imaging in some patients (RUO).”
- **(3) Timing claim (aggressive)**: “Detects resistance 3–6 months early.”

**Decision:** ____

### **A2. Are we allowed to recommend therapy changes, or only present ‘options’?**
Pick ONE:
- **(1) Options-only (RUO)**: “Consider discussing X with your oncologist.”
- **(2) Recommendation (clinical decision support)**: “We recommend X as Tier 1.”

**Decision:** ____

---

## B) Intervention Triggers (Numbers We Must Lock)

### **B1. DDR_bin trigger type**
Pick ONE:
- **Absolute threshold** (e.g., DDR_bin < 0.80)
- **Relative drop from baseline** (e.g., ΔDDR_bin ≤ -0.06)
- **Hybrid** (either absolute or relative)

**Decision:** ____

### **B2. Trigger thresholds (numbers)**
Fill in the exact values:
- **WARNING** when: ____
- **ALERT** when: ____
- **CRITICAL** when: ____

### **B3. Minimum confirmation**
Pick ONE:
- **Single sample triggers**
- **Requires 2 consecutive samples**
- **2-of-3 rule** (DDR_bin + CA-125 + ctDNA VAF change)

**Decision:** ____

### **B4. Monitoring cadence**
Pick ONE:
- q **4 weeks**
- q **8 weeks**
- q **12 weeks**
- **adaptive** (depends on trend)

**Decision:** ____

---

## C) Which Intervention Strategy Is Tier-1 (Default)

### **C1. At first WARNING (clone likely 3–10%), what is Tier-1?**
Pick ONE:
- **Adaptive therapy** (dose reduce)
- **Combination therapy** (add agent)
- **Trial-first** (enroll)

**Decision:** ____

### **C2. If combination therapy is Tier-1: which combination is our canonical story?**
Pick ONE:
- PARP + **low-dose carboplatin** (AUC 2–3)
- PARP + **bevacizumab**
- PARP + **ATR inhibitor**
- Other: ____

**Decision:** ____

### **C3. If adaptive therapy is Tier-1: what is the dosing rule?**
Provide:
- **Dose reduction %**: ____
- **Stop criteria** (if any): ____
- **Restart criteria**: ____

---

## D) Safety Guardrails (Hard Rules)

### **D1. “Do not escalate” conditions**
List the red-line conditions where our system must refuse to propose escalation (examples):
- Grade ≥3 toxicity
- Cytopenias below thresholds
- Renal impairment
- Pregnancy
- Provider override

**Guardrails:** ____

### **D2. Evidence tiers**
Confirm the tiers:
- **Tier 1**: RCT / guideline
- **Tier 2**: prospective / strong observational
- **Tier 3**: mechanistic / preclinical / hypothesis

**Decision:** ____

---

## E) Data Requirements (Reality Check)

### **E1. Do we require ctDNA to compute DDR_bin?**
Pick ONE:
- **ctDNA required**
- **tissue acceptable**
- **either** (prefer ctDNA)

**Decision:** ____

### **E2. Minimum variant count for DDR_bin stability**
Provide:
- **Minimum variants**: ____
- **If below minimum**: fallback behavior (show “insufficient data” / proxy signal)

---

## F) Validation Strategy (What “Done” Means)

### **F1. What is our MVP validation bar for Stage-3 demo?**
Pick ONE:
- **(1) Offline simulation + deterministic validator only**
- **(2) Retrospective cohort (single timepoint) + plausibility**
- **(3) Longitudinal retrospective evidence**
- **(4) Prospective pilot (clinical partner)**

**Decision:** ____

### **F2. Longitudinal data source**
Pick ONE:
- Partner clinic (prospective)
- Public cohort (name: ____)
- Pharma dataset (partner: ____)

**Decision:** ____

---

## G) Product Shape (What We Actually Build First)

### **G1. MVP surface**
Pick ONE:
- **(1) Dashboard only** (monitoring + alerts)
- **(2) Dashboard + recommendations** (tiered interventions)
- **(3) Full orchestration** (complete-care integrates prevention)

**Decision:** ____

### **G2. Where should prevention live in the product?**
Pick ONE:
- Inside **Complete Care** response
- Separate **/prevention** endpoint + UI

**Decision:** ____

---

## Output Requested From Manager
Please return answers in the same order as sections A–G.

