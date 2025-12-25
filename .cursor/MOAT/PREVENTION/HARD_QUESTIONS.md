# 🔥 HARD QUESTIONS: What We Don't Know Yet

**Date:** December 24, 2025  
**Purpose:** Questions for Manager to answer + Plumber to investigate  
**Status:** ⚠️ **BLOCKERS & UNKNOWNS**

---

## 🎯 CATEGORY 1: DDR_bin THRESHOLDS & TIMING

### **Q1.1: What DDR_bin threshold triggers intervention?**

We proposed:
- DDR_bin ≥ 0.85 → Continue monotherapy
- DDR_bin 0.80-0.85 → Consider intervention
- DDR_bin 0.70-0.80 → Escalate intervention
- DDR_bin < 0.70 → Switch therapy

**HARD QUESTION:** These thresholds are GUESSES. Where does the data come from?

- **For Manager:** Do we have longitudinal DDR_bin data from any cohort? Can we correlate DDR_bin values with time-to-progression?
- **For Plumber:** Can we retrospectively compute DDR_bin trajectories from existing TCGA-OV data (multiple timepoints)?


**MANAGER ANSWER:**

**Answer:** Use **relative drop from baseline (ΔDDR_bin)**, not absolute thresholds.
- **WARNING**: ΔDDR_bin ≤ **-0.05**
- **ALERT**: ΔDDR_bin ≤ **-0.08**
- **CRITICAL**: ΔDDR_bin ≤ **-0.12**
Always show **RUO** + "discuss with oncologist".

**Confidence:** MEDIUM

**Source:** General biomarker monitoring practice (delta-from-baseline) + need for calibration given variant count noise; not a published DDR_bin standard yet.

**Plumber Action:** Implement delta-based alerting + run stability experiments (subsample variants n=5/10/20) to see if -0.05 is above noise floor.

---

### **Q1.2: How fast does DDR_bin drop before progression?**

We said "3-6 months before clinical progression" but:

- Is the drop LINEAR or EXPONENTIAL?
- Does 6% drop (0.88 → 0.82) always mean 5-8% resistant clone?
- How do we convert DDR_bin delta → resistant clone percentage?

**For Manager:** What's the mathematical relationship between DDR_bin and resistant clone fraction?

**For Plumber:** Can we model this from first principles (if resistant clone has feature X, and sensitive has feature Y, then DDR_bin = weighted average)?


**MANAGER ANSWER:**

**Answer:** We **cannot claim "3-6 months"** as a universal timing law yet. For demo, phrase as: "**Potential early signal** that may precede clinical progression; timing varies."

**Confidence:** LOW

**Source:** Lack of longitudinal DDR_bin ground truth.

**Plumber Action:** Build a simple "lead time estimator" that outputs **UNKNOWN** unless ≥3 timepoints exist; do not hardcode months.

---

### **Q1.3: How often should we measure DDR_bin?**

We proposed "every 3 months" but:

- Is that too infrequent? Could resistance explode between measurements?
- Is that too frequent? ctDNA costs money.
- Does measurement frequency depend on baseline DDR_bin?

**For Manager:** What's the optimal monitoring interval based on tumor kinetics?


**MANAGER ANSWER:**

**Answer:** **q8-12 weeks** (aligns with many oncology follow-up cadences), with "adaptive tightening" if trending down:
- Stable: q12 weeks
- WARNING: q6-8 weeks
- ALERT: q4-6 weeks (if feasible)

**Confidence:** MEDIUM

**Source:** Practical clinical cadence constraints + ctDNA cost/turnaround realities.

**Plumber Action:** Build cadence recommendation logic (rule-based) + surface cost/turnaround in UI.

---

## 🎯 CATEGORY 2: INTERVENTION EFFICACY (THE BIGGEST UNKNOWN)

### **Q2.1: Does combination therapy actually suppress resistant clone?**

We proposed "add carboplatin AUC 2-3 when DDR_bin drops."

**HARD QUESTION:** This is based on PAOLA-1 (bevacizumab, not carboplatin). Do we have evidence that:

1. Low-dose carboplatin kills HR-proficient (resistant) cells preferentially?
2. Adding carboplatin at 5-8% resistant clone prevents dominance?
3. DDR_bin stabilizes after adding carboplatin?

**For Manager:** Is there ANY clinical data on DDR_bin (or similar marker) response to combination therapy?

**For Plumber:** Can we simulate this in silico? (Model resistant clone dynamics under different intervention scenarios)


**MANAGER ANSWER:**

**Answer:** We must **not claim suppression is proven**. For MVP we present **tiered options**:
- **Tier 1 (highest real-world precedent in OV)**: **PARP + bevacizumab** (PAOLA-1 is real RCT precedent; not DDR_bin triggered but combo exists).
- **Tier 2**: PARP + ATR (trial-first framing).
- **Tier 3**: PARP + low-dose carboplatin (hypothesis; toxicity-flagged).

**Confidence:** HIGH on "don't overclaim"; MEDIUM on tier ordering.

**Source:** Combination regimens exist; DDR_bin-triggered combo does not yet have direct evidence.

**Plumber Action:** Recommendation engine must label tiers + evidence level and avoid causal language ("will suppress").

---

### **Q2.2: Does adaptive therapy work in ovarian cancer?**

Gatenby's work is in PROSTATE cancer. We're assuming it translates to ovarian.

**HARD QUESTIONS:**
- Do ovarian cancer resistant clones have the same fitness cost as prostate?
- Is the sensitive cell regrowth rate fast enough to outcompete resistant cells?
- What's the optimal dose reduction (50%? 75%? Complete holiday?)?

**For Manager:** Has adaptive therapy been tested in any gynecologic malignancy?


**MANAGER ANSWER:**

**Answer:** Treat as **experimental / hypothesis** in OV. Include it as **Tier 3** unless Manager later upgrades based on literature review.

**Confidence:** LOW-MEDIUM

**Source:** Adaptive therapy evidence strongest outside OV; translation uncertain.

**Plumber Action:** Keep adaptive therapy template, but default to "experimental / trial context".

---

### **Q2.3: When does evolutionary steering work vs. fail?**

We proposed: PARP → ATR inhibitor → PARP cycling.

**HARD QUESTIONS:**
- What if the resistant clone develops ATR independence?
- What if the sensitive clone doesn't regrow during ATR phase?
- How many cycles can you do before both populations are exhausted?

**For Manager:** What's the failure mode of evolutionary steering? When does it stop working?


**MANAGER ANSWER:**

**Answer:** For MVP: **RUO-only, trial-first**. Do not imply it is standard or reliably effective.
- Failure modes to explicitly surface: pathway redundancy, ATR-independence, toxicity limits, no sensitive regrowth.

**Confidence:** MEDIUM (on failure modes), LOW (on efficacy frequency).

**Source:** General evolutionary therapy concepts; limited OV-specific proof.

**Plumber Action:** UI must show a "Known unknowns / failure modes" box whenever steering is displayed.

---

## 🎯 CATEGORY 3: BIOLOGICAL ASSUMPTIONS

### **Q3.1: Is DDR_bin actually measuring HR restoration?**

Our 9 diamond features all map to DDR pathway (TP53 dominant).

**HARD QUESTION:** But TP53 is not the same as RAD51C/BRCA1 reversion.

- Does DDR_bin drop when RAD51C reverts? (We assume yes, but haven't tested)
- Does DDR_bin drop when BRCA1 reverts?
- Does DDR_bin drop for OTHER resistance mechanisms (e.g., drug efflux)?

**For Manager:** What exactly is DDR_bin measuring at the molecular level?

**For Plumber:** Can we test this by:
1. Finding patients with known RAD51C reversions → check if DDR_bin is low
2. Finding patients with known drug efflux → check if DDR_bin is unchanged


**MANAGER ANSWER:**

**Answer:** Right now, we should describe DDR_bin as: "**DDR-pathway signal correlated with resistance**" **not** "HR restoration confirmed."

**Confidence:** MEDIUM

**Source:** Your own mapping shows **TP53 dominance**; HR restoration requires specific reversion evidence.

**Plumber Action:** Add a "mechanism confirmation" step: if RAD51C/BRCA reversion detected → upgrade mechanism label to "HR restoration supported."

---

### **Q3.2: Are there multiple resistance mechanisms at once?**

We modeled "single resistant clone emerges." But in reality:

- Could a patient have BOTH RAD51C reversion AND ABCB1 upregulation?
- If DDR_bin drops AND Efflux_bin rises, which intervention wins?
- Do we need a multi-bin decision tree?

**For Manager:** How often do patients have multiple concurrent resistance mechanisms?


**MANAGER ANSWER:**

**Answer:** Assume **yes** (real tumors are messy). MVP should support multi-signal output:
- DDR_bin drop + MAPK rise + Efflux rise can coexist; actions must prioritize safety + evidence.

**Confidence:** HIGH

**Source:** Broad oncology reality (poly-clonal evolution).

**Plumber Action:** Multi-bin decision engine with tie-break rules (e.g., safety/contraindications > evidence tier > mechanism strength).

---

### **Q3.3: What about non-DDR resistance (MAPK, PI3K)?**

We focused on DDR_bin, but the Manager's `idea.mdc` mentioned MAPK_bin and PI3K_bin.

**HARD QUESTIONS:**
- Do we have validated MAPK_bin features? (We have MAPK pathway genes, but no TRUE SAE validation)
- What interventions work for MAPK escape? (MEK inhibitors? Which ones?)
- How do we combine DDR_bin + MAPK_bin + PI3K_bin in decision logic?

**For Plumber:** We need to:
1. Extract MAPK_bin features from Tier-3 cohort (if signal exists)
2. Validate PI3K_bin features
3. Build multi-bin decision engine


**MANAGER ANSWER:**

**Answer:** For MVP, **don't block on TRUE-SAE MAPK/PI3K bins**. Use **proxy pathway burdens** (gene/pathway aggregation) if TRUE bins aren't validated.

**Confidence:** HIGH (pragmatic).

**Source:** Avoid waiting on new TRUE-SAE validations.

**Plumber Action:** Implement MAPK/PI3K as proxy bins first; keep TRUE-bin hooks as optional later.

---

## 🎯 CATEGORY 4: CLINICAL FEASIBILITY

### **Q4.1: Can oncologists actually use this?**

We built a complex decision system:
- DDR_bin monitoring every 3 months
- Multiple intervention tiers
- Evolutionary steering (cycling therapies)

**HARD QUESTION:** Is this realistic in clinical practice?

- Do oncologists have time to interpret DDR_bin dashboards?
- Will patients tolerate frequent ctDNA draws?
- Are ATR inhibitors actually accessible? (Most are Phase 1)

**For Manager:** What's the minimum viable product (MVP) for clinical adoption?


**MANAGER ANSWER:**

**Answer:** MVP = **Dashboard + alerting + options-only tiered actions + next-test suggestions**. No complex cycling plans as "recommendations."

**Confidence:** HIGH

**Source:** Adoption reality; reduce cognitive load.

**Plumber Action:** Build the simplest UI that shows: baseline, trend, alert, top 3 actions, next test.

---

### **Q4.2: What's the regulatory pathway?**

We mentioned "FDA companion diagnostic" but:

- Is DDR_bin a diagnostic, prognostic, or predictive biomarker? (Different pathways)
- Do we need a prospective clinical trial to validate?
- What's the timeline and cost for FDA clearance?

**For Manager:** What's the realistic regulatory strategy?


**MANAGER ANSWER:**

**Answer:** Start as **RUO / decision support**. Treat "companion diagnostic" as **future** (requires prospective evidence).

**Confidence:** HIGH

**Source:** Typical SaMD/CDx requirements.

**Plumber Action:** Hardcode RUO banners + provenance in every prevention output.

---

### **Q4.3: Who pays for this?**

- ctDNA every 3 months = $500-1000/test × 4/year = $2000-4000/patient/year
- DDR_bin computation = software cost
- Clinical decision support = integration cost

**For Manager:** What's the reimbursement model? Will payers cover prevention-based monitoring?


**MANAGER ANSWER:**

**Answer:** For now: piggyback on existing **ctDNA reimbursement**; we are the software layer. Position as **avoiding late progression costs** (hospitalizations, failed lines).

**Confidence:** MEDIUM

**Source:** Payer logic (avoid high downstream costs), but needs real-world economics later.

**Plumber Action:** Surface "estimated cost / cadence" in outputs, but don't promise reimbursement.

---

## 🎯 CATEGORY 5: DATA & VALIDATION GAPS

### **Q5.1: Do we have longitudinal data to validate?**

Our Tier-3 cohort (149 patients) is SINGLE TIMEPOINT.

**HARD QUESTION:** We can't validate DDR_bin trajectory without longitudinal data.

- Does TCGA-OV have multiple timepoints per patient? (Probably not)
- Is there a public dataset with serial ctDNA samples + outcomes?
- Do we need to partner with a clinical site for prospective collection?

**For Manager:** Where can we get longitudinal DDR_bin data?

**For Plumber:** Research available datasets:
- MMRF CoMMpass (has longitudinal for MM)
- GRAIL CCGA (ctDNA cohort)
- POG570 (BC Cancer longitudinal)


**MANAGER ANSWER:**

**Answer:** TCGA-OV is mostly single timepoint. We likely need: **clinical partner prospective pilot** or a known longitudinal research cohort.

**Confidence:** HIGH

**Source:** TCGA structure.

**Plumber Action:** Implement pipeline so it works with longitudinal when available; don't wait for it to build MVP.

---

### **Q5.2: Sample size for intervention validation?**

If we want to run a clinical trial:
- How many patients needed to show PFS improvement?
- What's the effect size we expect? (50% PFS improvement = how many patients?)
- Can we do a biomarker-stratified trial (enroll only DDR_bin-drop patients)?

**For Manager:** What's the power calculation for a DDR_bin-guided intervention trial?


**MANAGER ANSWER:**

**Answer:** **Pilot**: 30-60 triggered patients (feasibility, biomarker trajectory, safety)
**Efficacy RCT**: likely **200-400+** depending on effect size/trigger rate.

**Confidence:** LOW-MEDIUM (needs statistician + event rate assumptions).

**Source:** General survival trial sizing intuition.

**Plumber Action:** None required for MVP; keep trial design doc as placeholder until stats consult.

---

### **Q5.3: How do we handle ctDNA-negative patients?**

Not all patients shed ctDNA.

- If ctDNA is undetectable, can we still compute DDR_bin?
- Do we fall back to tissue-based SAE?
- What's the false negative rate for ctDNA-based monitoring?

**For Manager:** What's the backup plan for ctDNA-negative patients?


**MANAGER ANSWER:**

**Answer:** Must support "**insufficient data**" mode + fallback:
- Tissue-based re-biopsy if clinically reasonable, or proxy signals.

**Confidence:** HIGH

**Source:** ctDNA shedding variability.

**Plumber Action:** If variant_count < minimum → set confidence LOW + recommend alternative test.

---

## 🎯 CATEGORY 6: TECHNICAL IMPLEMENTATION

### **Q6.1: How do we compute DDR_bin from ctDNA?**

Current TRUE SAE pipeline:
1. Get variant list from ctDNA
2. Extract reference + alt sequences (101bp)
3. Run through Evo2 → activations
4. Run through SAE → 32K sparse features
5. Aggregate features into DDR_bin

**HARD QUESTIONS:**
- ctDNA has LOW variant counts (maybe 5-10 variants vs. 50-100 in tumor)
- Is DDR_bin stable with low variant counts?
- What's the minimum variant count for reliable DDR_bin?

**For Plumber:** Test DDR_bin stability:
1. Subsample Tier-3 variants (simulate ctDNA-like sparsity)
2. Recompute DDR_bin
3. Measure variance (is it still discriminative?)


**MANAGER ANSWER:**

**Answer:** **Minimum variants for "usable"**: **≥10** (below that = LOW confidence)
**Strong confidence**: **≥20**

**Confidence:** MEDIUM (must be empirically checked).

**Source:** Variance/noise considerations.

**Plumber Action:** Run subsampling stability study and set thresholds from measured CV.

---

### **Q6.2: How fast can we compute DDR_bin?**

Current pipeline runs in Modal (cloud GPU).

- What's the turnaround time from ctDNA result → DDR_bin score?
- Can we batch process? (Multiple patients at once)
- What's the cost per patient?

**For Plumber:** Benchmark the pipeline:
1. Time from variant input → DDR_bin output
2. Cost per patient (Modal GPU cost)
3. Optimize if needed


**MANAGER ANSWER:**

**Answer:** For MVP, target **< 1 hour** compute after variants available; hard upper bound **<24h**. Cost per patient must be tracked but not promised.

**Confidence:** MEDIUM

**Source:** Operational feasibility.

**Plumber Action:** Benchmark + log runtime and estimated cost in provenance.

---

### **Q6.3: How do we version and track DDR_bin over time?**

Patient gets DDR_bin at Month 0, 3, 6, 9...

- How do we store this in patient state?
- How do we visualize trends?
- How do we handle recalibration (if model improves)?

**For Plumber:** Design the data schema:
```json
{
  "patient_id": "AK_001",
  "ddr_bin_history": [
    {"date": "2025-01-15", "value": 0.88, "variant_count": 52, "model_version": "v1.0"},
    {"date": "2025-04-15", "value": 0.87, "variant_count": 48, "model_version": "v1.0"},
    {"date": "2025-07-15", "value": 0.82, "variant_count": 55, "model_version": "v1.0"}
  ]
}
```


**MANAGER ANSWER:**

**Answer:** Always store: model_version, variant_count, source (ctDNA/tissue), timestamp, baseline used.

**Confidence:** HIGH

**Source:** Auditability requirement.

**Plumber Action:** Implement state schema exactly like in `PLUMBER_BUILD_SPEC.md`.

---

## 🎯 CATEGORY 7: STRATEGIC QUESTIONS

### **Q7.1: Should we publish DDR_bin before clinical validation?**

We have TRUE SAE AUROC 0.783. We could publish now.

**TRADE-OFFS:**
- Publish now → Establish priority, attract collaborators
- Wait for clinical validation → Stronger paper, less risk of failure

**For Manager:** What's the publication strategy?


**MANAGER ANSWER:**

**Answer:** Yes, publish as **retrospective / method + association**, with explicit limitations (no longitudinal validation yet).

**Confidence:** MEDIUM

**Source:** Standard path: publish signal → attract collaborators → run pilot.

**Plumber Action:** Keep outputs reproducible + versioned to support paper artifacts.

---

### **Q7.2: Who are our competitors?**

- Foundation Medicine (ctDNA-based monitoring)
- Guardant Health (LUNAR, ctDNA for MRD)
- Grail (CCGA, multi-cancer early detection)

**HARD QUESTION:** Do any of them have pathway-bin-level monitoring?

**For Manager:** What's our competitive moat vs. Foundation/Guardant?


**MANAGER ANSWER:**

**Answer:** Competitors do ctDNA/MRD; most **don't** provide a **mechanism-bin + action framework** tied to interpretable features + trial matching.

**Confidence:** MEDIUM (needs market scan).

**Source:** General landscape knowledge.

**Plumber Action:** None.

---

### **Q7.3: Should we partner with pharma?**

ATR inhibitors (ceralasertib) are owned by AstraZeneca.
RAD51 inhibitors (CYT-0851) are owned by Cyteir.

**HARD QUESTION:** Do we need pharma partnerships to enable evolutionary steering?

**For Manager:** What's the partnership strategy?

## 📊 PRIORITY RANKING

### **Must Answer Before Building:**
1. **Q1.1** - DDR_bin thresholds (what triggers intervention?)
2. **Q2.1** - Does combination therapy actually work?
3. **Q3.1** - Is DDR_bin measuring HR restoration?
4. **Q5.1** - Do we have longitudinal data?

### **Must Answer Before Clinical Trial:**
5. **Q4.2** - Regulatory pathway
6. **Q5.2** - Sample size calculation
7. **Q2.2** - Does adaptive therapy work in OV?

### **Must Answer Before Scale:**
8. **Q4.3** - Reimbursement model
9. **Q6.1** - ctDNA-based DDR_bin stability
10. **Q7.2** - Competitive positioning

---

## 🛠️ PLUMBER INVESTIGATION TASKS

### **Week 1: Data Investigation**
1. Check TCGA-OV for longitudinal samples (probably none)
2. Research POG570, GRAIL CCGA for longitudinal ctDNA
3. Test DDR_bin stability with subsampled variants

### **Week 2: Modeling**
1. Model DDR_bin → resistant clone relationship
2. Simulate intervention scenarios (combo vs. adaptive)
3. Build multi-bin decision tree (DDR + MAPK + PI3K)

### **Week 3: Technical**
1. Benchmark DDR_bin computation time/cost
2. Design patient state schema for longitudinal DDR_bin
3. Build trend analysis algorithm

---

## �� MANAGER RESPONSE TEMPLATE

For each question, we need:
1. **Answer** (if known)
2. **Confidence** (HIGH/MEDIUM/LOW/UNKNOWN)
3. **Source** (literature, expert opinion, first principles)
4. **Action** (what Plumber should do next)


**MANAGER ANSWER:**

**Answer:** Use **relative drop from baseline (ΔDDR_bin)**, not absolute thresholds.
- **WARNING**: ΔDDR_bin ≤ **-0.05**
- **ALERT**: ΔDDR_bin ≤ **-0.08**
- **CRITICAL**: ΔDDR_bin ≤ **-0.12**
Always show **RUO** + "discuss with oncologist".

**Confidence:** MEDIUM

**Source:** General biomarker monitoring practice (delta-from-baseline) + need for calibration given variant count noise; not a published DDR_bin standard yet.

**Plumber Action:** Implement delta-based alerting + run stability experiments (subsample variants n=5/10/20) to see if -0.05 is above noise floor.

---


**MANAGER ANSWER:**

**Answer:** We **cannot claim "3-6 months"** as a universal timing law yet. For demo, phrase as: "**Potential early signal** that may precede clinical progression; timing varies."

**Confidence:** LOW

**Source:** Lack of longitudinal DDR_bin ground truth.

**Plumber Action:** Build a simple "lead time estimator" that outputs **UNKNOWN** unless ≥3 timepoints exist; do not hardcode months.

---


**MANAGER ANSWER:**

**Answer:** **q8-12 weeks** (aligns with many oncology follow-up cadences), with "adaptive tightening" if trending down:
- Stable: q12 weeks
- WARNING: q6-8 weeks
- ALERT: q4-6 weeks (if feasible)

**Confidence:** MEDIUM

**Source:** Practical clinical cadence constraints + ctDNA cost/turnaround realities.

**Plumber Action:** Build cadence recommendation logic (rule-based) + surface cost/turnaround in UI.

---


**MANAGER ANSWER:**

**Answer:** We must **not claim suppression is proven**. For MVP we present **tiered options**:
- **Tier 1 (highest real-world precedent in OV)**: **PARP + bevacizumab** (PAOLA-1 is real RCT precedent; not DDR_bin triggered but combo exists).
- **Tier 2**: PARP + ATR (trial-first framing).
- **Tier 3**: PARP + low-dose carboplatin (hypothesis; toxicity-flagged).

**Confidence:** HIGH on "don't overclaim"; MEDIUM on tier ordering.

**Source:** Combination regimens exist; DDR_bin-triggered combo does not yet have direct evidence.

**Plumber Action:** Recommendation engine must label tiers + evidence level and avoid causal language ("will suppress").

---


**MANAGER ANSWER:**

**Answer:** Treat as **experimental / hypothesis** in OV. Include it as **Tier 3** unless Manager later upgrades based on literature review.

**Confidence:** LOW-MEDIUM

**Source:** Adaptive therapy evidence strongest outside OV; translation uncertain.

**Plumber Action:** Keep adaptive therapy template, but default to "experimental / trial context".

---


**MANAGER ANSWER:**

**Answer:** For MVP: **RUO-only, trial-first**. Do not imply it is standard or reliably effective.
- Failure modes to explicitly surface: pathway redundancy, ATR-independence, toxicity limits, no sensitive regrowth.

**Confidence:** MEDIUM (on failure modes), LOW (on efficacy frequency).

**Source:** General evolutionary therapy concepts; limited OV-specific proof.

**Plumber Action:** UI must show a "Known unknowns / failure modes" box whenever steering is displayed.

---


**MANAGER ANSWER:**

**Answer:** Right now, we should describe DDR_bin as: "**DDR-pathway signal correlated with resistance**" **not** "HR restoration confirmed."

**Confidence:** MEDIUM

**Source:** Your own mapping shows **TP53 dominance**; HR restoration requires specific reversion evidence.

**Plumber Action:** Add a "mechanism confirmation" step: if RAD51C/BRCA reversion detected → upgrade mechanism label to "HR restoration supported."

---


**MANAGER ANSWER:**

**Answer:** Assume **yes** (real tumors are messy). MVP should support multi-signal output:
- DDR_bin drop + MAPK rise + Efflux rise can coexist; actions must prioritize safety + evidence.

**Confidence:** HIGH

**Source:** Broad oncology reality (poly-clonal evolution).

**Plumber Action:** Multi-bin decision engine with tie-break rules (e.g., safety/contraindications > evidence tier > mechanism strength).

---


**MANAGER ANSWER:**

**Answer:** For MVP, **don't block on TRUE-SAE MAPK/PI3K bins**. Use **proxy pathway burdens** (gene/pathway aggregation) if TRUE bins aren't validated.

**Confidence:** HIGH (pragmatic).

**Source:** Avoid waiting on new TRUE-SAE validations.

**Plumber Action:** Implement MAPK/PI3K as proxy bins first; keep TRUE-bin hooks as optional later.

---


**MANAGER ANSWER:**

**Answer:** MVP = **Dashboard + alerting + options-only tiered actions + next-test suggestions**. No complex cycling plans as "recommendations."

**Confidence:** HIGH

**Source:** Adoption reality; reduce cognitive load.

**Plumber Action:** Build the simplest UI that shows: baseline, trend, alert, top 3 actions, next test.

---


**MANAGER ANSWER:**

**Answer:** Start as **RUO / decision support**. Treat "companion diagnostic" as **future** (requires prospective evidence).

**Confidence:** HIGH

**Source:** Typical SaMD/CDx requirements.

**Plumber Action:** Hardcode RUO banners + provenance in every prevention output.

---


**MANAGER ANSWER:**

**Answer:** For now: piggyback on existing **ctDNA reimbursement**; we are the software layer. Position as **avoiding late progression costs** (hospitalizations, failed lines).

**Confidence:** MEDIUM

**Source:** Payer logic (avoid high downstream costs), but needs real-world economics later.

**Plumber Action:** Surface "estimated cost / cadence" in outputs, but don't promise reimbursement.

---


**MANAGER ANSWER:**

**Answer:** TCGA-OV is mostly single timepoint. We likely need: **clinical partner prospective pilot** or a known longitudinal research cohort.

**Confidence:** HIGH

**Source:** TCGA structure.

**Plumber Action:** Implement pipeline so it works with longitudinal when available; don't wait for it to build MVP.

---


**MANAGER ANSWER:**

**Answer:** **Pilot**: 30-60 triggered patients (feasibility, biomarker trajectory, safety)
**Efficacy RCT**: likely **200-400+** depending on effect size/trigger rate.

**Confidence:** LOW-MEDIUM (needs statistician + event rate assumptions).

**Source:** General survival trial sizing intuition.

**Plumber Action:** None required for MVP; keep trial design doc as placeholder until stats consult.

---


**MANAGER ANSWER:**

**Answer:** Must support "**insufficient data**" mode + fallback:
- Tissue-based re-biopsy if clinically reasonable, or proxy signals.

**Confidence:** HIGH

**Source:** ctDNA shedding variability.

**Plumber Action:** If variant_count < minimum → set confidence LOW + recommend alternative test.

---


**MANAGER ANSWER:**

**Answer:** **Minimum variants for "usable"**: **≥10** (below that = LOW confidence)
**Strong confidence**: **≥20**

**Confidence:** MEDIUM (must be empirically checked).

**Source:** Variance/noise considerations.

**Plumber Action:** Run subsampling stability study and set thresholds from measured CV.

---


**MANAGER ANSWER:**

**Answer:** For MVP, target **< 1 hour** compute after variants available; hard upper bound **<24h**. Cost per patient must be tracked but not promised.

**Confidence:** MEDIUM

**Source:** Operational feasibility.

**Plumber Action:** Benchmark + log runtime and estimated cost in provenance.

---


**MANAGER ANSWER:**

**Answer:** Always store: model_version, variant_count, source (ctDNA/tissue), timestamp, baseline used.

**Confidence:** HIGH

**Source:** Auditability requirement.

**Plumber Action:** Implement state schema exactly like in `PLUMBER_BUILD_SPEC.md`.

---


**MANAGER ANSWER:**

**Answer:** Yes, publish as **retrospective / method + association**, with explicit limitations (no longitudinal validation yet).

**Confidence:** MEDIUM

**Source:** Standard path: publish signal → attract collaborators → run pilot.

**Plumber Action:** Keep outputs reproducible + versioned to support paper artifacts.

---


**MANAGER ANSWER:**

**Answer:** Competitors do ctDNA/MRD; most **don't** provide a **mechanism-bin + action framework** tied to interpretable features + trial matching.

**Confidence:** MEDIUM (needs market scan).

**Source:** General landscape knowledge.

**Plumber Action:** None.

---


**MANAGER ANSWER:**

**Answer:** Needed for **Stage-4 steering** (ATR/RAD51 access), not needed for **Stage-3 MVP demo**.

**Confidence:** HIGH

**Source:** Drug access realities.

**Plumber Action:** Keep trial-first hooks; don't hardcode specific pharma dependencies.

---

---

*Document Owner: Zo*  
*Last Updated: December 24, 2025*  
*Status: ✅ MANAGER ANSWERS INTEGRATED (December 24, 2025)*


**MANAGER ANSWER:**

**Answer:** Needed for **Stage-4 steering** (ATR/RAD51 access), not needed for **Stage-3 MVP demo**.

**Confidence:** HIGH

**Source:** Drug access realities.

**Plumber Action:** Keep trial-first hooks; don't hardcode specific pharma dependencies.

---
