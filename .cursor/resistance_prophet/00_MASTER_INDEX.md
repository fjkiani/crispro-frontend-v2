# ⚔️ RESISTANCE PROPHET: MASTER INDEX ⚔️

**Date:** January 14, 2025  
**Status:** APPROVED - Ready for Implementation  
**Timeline:** 3 weeks to production  
**Owner:** Zo (prediction engine) + Jr (dashboard integration)

---

## 📋 **NAVIGATION**

### **Core Documentation**
- **[Architecture Overview](architecture/00_ARCHITECTURE_OVERVIEW.md)** - System design & philosophy
- **[Implementation Guide](implementation/00_IMPLEMENTATION_GUIDE.md)** - Step-by-step build plan
- **[Integration Guide](integration/00_INTEGRATION_GUIDE.md)** - Mission Control integration
- **[Validation Guide](validation/00_VALIDATION_GUIDE.md)** - Testing & metrics

---

## 🎯 **THE VISION: ONE INTEGRATED SYSTEM**

### **What We're Building:**

**Resistance Prophet** = The prediction engine that forecasts treatment failure 3-6 months early

**Mission Control** = The dashboard that surfaces predictions and enables immediate action

**Together:** A complete system that predicts resistance → displays insights → enables intervention

---

## 🧠 **THE BREAKTHROUGH INSIGHT**

### **What the Manager Revealed:**

> "Instant Oracle overlaps with existing Co-Pilot endpoint - enhance it instead of creating new one"

**Translation:**
- ❌ **Don't build:** Separate "Instant Oracle" agent
- ✅ **Do enhance:** `/api/ayesha/complete_care_v2` with Resistance Prophet predictions
- ✅ **Mission Control:** Calls enhanced endpoint, displays predictions

### **What Jr's Data Hunt Revealed:**

**161 patients with:**
- ✅ Genomic data (mutations, SAE features)
- ✅ Platinum response outcomes (sensitive/resistant/refractory)
- ✅ Clinical data (stage, OS, treatment history)

**This IS the resistance prediction training set.**

---

## 📋 **QUICK START**

### **For Zo (Backend Implementation):**
1. Read [Architecture Overview](architecture/00_ARCHITECTURE_OVERVIEW.md)
2. Follow [CA-125 Service Implementation](implementation/01_CA125_SERVICE.md)
3. Follow [DNA Repair Detection](implementation/02_DNA_REPAIR_DETECTION.md)
4. Follow [Pathway Escape Detection](implementation/03_PATHWAY_ESCAPE_DETECTION.md)
5. Follow [Backend Integration](implementation/04_BACKEND_INTEGRATION.md)
6. Execute [Validation Plan](validation/01_RETROSPECTIVE_VALIDATION.md)

### **For Jr (Dashboard Integration):**
1. Read [Integration Overview](integration/00_INTEGRATION_GUIDE.md)
2. Implement [Critical Alerts](integration/01_CRITICAL_ALERTS.md)
3. Implement [Resistance Playbook](integration/02_RESISTANCE_PLAYBOOK.md)
4. Implement [Next Best Actions](integration/03_NEXT_BEST_ACTIONS.md)
5. Execute [E2E Testing](validation/02_E2E_TESTING.md)

---

## 🎯 **SUCCESS METRICS**

### **Resistance Prophet Validation (Week 1):**
- ✅ Test on Jr's 161 patients (retrospective validation)
- ✅ Target: >70% accuracy for resistance prediction
- ✅ Target: 3-6 weeks earlier detection vs imaging alone

### **Mission Control Integration (Week 3):**
- ✅ Resistance predictions display in real-time
- ✅ Critical alerts fire when resistance detected (>70% probability)
- ✅ Next best actions include resistance recommendations
- ✅ Complete care plan includes resistance playbook

---

## 📋 **COMPLETE TIMELINE**

### **Week 1: Resistance Prophet Core**
- [ ] Day 1-2: Build `ResistanceProphetService` (CA-125 kinetics)
- [ ] Day 3: Add DNA repair restoration detection
- [ ] Day 4: Add pathway escape detection
- [ ] Day 5: Test on Jr's 161 patients (validation)

### **Week 2: Backend Integration**
- [ ] Day 6-7: Integrate into `/api/ayesha/complete_care_v2`
- [ ] Day 8: Add resistance prediction endpoint
- [ ] Day 9: Test end-to-end (API → predictions)
- [ ] Day 10: Document + provenance

### **Week 3: Mission Control Integration (Jr)**
- [ ] Day 11-12: Critical Alerts integration
- [ ] Day 13-14: Resistance Playbook card
- [ ] Day 15: Next Best Actions integration
- [ ] Day 16-17: Testing + polish

---

## ⚔️ **WHAT MAKES THIS REVOLUTIONARY**

### **Before:**
- ❌ Resistance detected AFTER imaging shows progression
- ❌ Oncologists react to failure (too late to intervene)
- ❌ Patients lose 3-6 months of treatment window

### **After (With Resistance Prophet + Mission Control):**
- ✅ Resistance predicted 3-6 WEEKS before imaging
- ✅ Oncologists switch therapy PROACTIVELY (before failure)
- ✅ Patients preserve treatment window (better outcomes)

### **The Competitive Advantage:**
**NOBODY ELSE PREDICTS RESISTANCE BEFORE PROGRESSION.**

- Tempus: Reports mutations, doesn't predict resistance
- Foundation Medicine: Reports biomarkers, doesn't predict resistance
- Guardant: Reports ctDNA changes, doesn't predict resistance

**We integrate CA-125 + SAE + Treatment Line Intelligence → Predict failure BEFORE it happens.**

---

## 🚨 **MANAGER QUESTIONS: PREVENT HALLUCINATION** 🚨

**Date:** January 14, 2025  
**Purpose:** Get explicit manager approval before implementation to prevent deviation

---

### **SECTION 1: DATA & VALIDATION SCOPE**

**Q1 (Validation Dataset):**
- Jr found 161 patients with platinum response data (sensitive/resistant/refractory).
- Should I use ALL 161 patients for validation, or stratify by stage/treatment line?
- Minimum N required: 40 (per manager's previous guidance) - we have 161 (4x target). Proceed?

**Q2 (Platinum Response Mapping):**
- Jr's data has platinum response categories (sensitive/resistant/refractory).
- Should I treat this as:
  - **Binary**: sensitive vs. resistant/refractory (simpler, higher N)
  - **3-class**: sensitive vs. resistant vs. refractory (more granular, lower N per class)
- Manager's preference?

**Q3 (CA-125 Validation Data):**
- Current validation dataset: Jr's 161 patients (genomics + platinum response).
- **CRITICAL GAP**: Do we have CA-125 history (time-series values) for these 161 patients?
- If NOT available from GDC/cBio, should I:
  - **Option A**: Validate DNA repair + pathway escape ONLY (skip CA-125 for now)
  - **Option B**: Use CA-125 kinetics as "prospective validation" (test on Ayesha live)
  - **Option C**: Find another data source for CA-125 histories
- Manager's decision?

**Q4 (Success Metrics - Retrospective):**
- For retrospective validation (Jr's 161 patients):
  - **Primary Metric**: AUROC for predicting platinum response (target ≥0.70)?
  - **Secondary Metrics**: Sensitivity, specificity, PPV?
  - **Subgroup Analysis**: Stage IV only vs. IIIC+IV?
- Confirm metrics before I run validation?

---

### **SECTION 2: ARCHITECTURE & INTEGRATION**

**Q5 (CA-125 Service Integration):**
- Current: `ca125_intelligence.py` exists and is operational.
- Plan: `ResistanceProphetService` calls `CA125Intelligence.analyze()` for kinetics.
- **Question**: Should I extend `ca125_intelligence.py` with resistance detection methods, or keep all resistance logic in `ResistanceProphetService`?
- Manager's preference for separation of concerns?

**Q6 (SAE Feature Source):**
- Current: `sae_feature_service.py` exists.
- Plan: `ResistanceProphetService` calls `SAEFeatureService.extract_features()` for DNA repair capacity and mechanism vector.
- **Question**: Should I add new methods to `sae_feature_service.py` (e.g., `detect_dna_repair_restoration()`), or keep detection logic in `ResistanceProphetService`?
- Manager's preference?

**Q7 (Instant Oracle Endpoint):**
- Manager said: "Instant Oracle overlaps with existing Co-Pilot endpoint - enhance it instead."
- Current endpoint: `/api/ayesha/complete_care_v2`
- **Plan**: Add `resistance_prediction` key to existing response schema.
- **Question**: Should resistance prediction be:
  - **Always computed** (default behavior for all `/complete_care_v2` calls)
  - **Opt-in** (new query parameter like `include_resistance_prediction=true`)
  - **Conditional** (only if patient has CA-125 history OR baseline SAE features)
- Manager's decision?

**Q8 (Mission Control Integration - Jr's Work):**
- Jr will integrate Resistance Prophet predictions into Mission Control dashboard.
- **Question**: Should resistance predictions be:
  - **Displayed immediately** (in Critical Alerts banner on page load)
  - **On-demand** (user clicks "Check Resistance Risk" button)
  - **Proactive** (auto-fire alert if resistance probability >70%)
- Manager's UX preference?

---

### **SECTION 3: CLINICAL DECISION SUPPORT**

**Q9 (Resistance Probability Thresholds):**
- Plan: Use 2-of-3 signal detection (CA-125 + DNA repair + pathway escape).
- **Question**: What probability thresholds should trigger different urgency levels?
  - **HIGH urgency** (immediate action): ≥70% probability?
  - **MEDIUM urgency** (discuss at next visit): 50-70% probability?
  - **LOW urgency** (monitor): <50% probability?
- Manager's clinical guidance?

**Q10 (Recommendation Actions):**
- Plan: Generate actionable recommendations based on resistance signals.
- **Question**: What actions should be recommended for HIGH urgency (≥70% probability)?
  - **ESCALATE_IMAGING**: Order CT/PET within 1 week?
  - **CONSIDER_SWITCH**: Discuss therapy switch within 2 weeks?
  - **REVIEW_RESISTANCE_PLAYBOOK**: Prepare next-line options within 1 week?
- Manager's clinical workflow alignment?

**Q11 (Next-Line Options Source):**
- Plan: When resistance detected, suggest next-line therapies (e.g., Niraparib + Bevacizumab).
- **Question**: Should next-line options come from:
  - **Static rules** (hardcoded in `ResistanceProphetService`)
  - **Treatment Line Intelligence** (call `TreatmentLineService.get_next_line_options()`)
  - **Resistance Playbook** (call `ResistancePlaybookService.get_combos()`)
- Manager's preference for source of truth?

---

### **SECTION 4: TIMELINE & SCOPE GATES**

**Q12 (3-Week Timeline Feasibility):**
- Proposed timeline:
  - **Week 1**: Build `ResistanceProphetService` core (CA-125 + DNA repair + pathway escape)
  - **Week 2**: Backend integration into `/api/ayesha/complete_care_v2` + validation on Jr's 161 patients
  - **Week 3**: Jr integrates into Mission Control dashboard (Critical Alerts, Resistance Playbook, Next Best Actions)
- **Question**: Is this realistic, or should I adjust scope/timeline?
- Manager's timeline approval?

**Q13 (Validation Requirements Before Production):**
- Plan: Validate on Jr's 161 patients (retrospective validation).
- **Question**: What validation metrics are required BEFORE production deployment?
  - **Minimum AUROC**: ≥0.70 for platinum response prediction?
  - **Minimum Sensitivity**: ≥0.75 for detecting resistance?
  - **Minimum Specificity**: ≥0.70 to avoid false alarms?
- Manager's go/no-go criteria?

**Q14 (Phase 2 Scope - ctDNA Integration):**
- Current plan focuses on CA-125 + SAE features (Phase 1 MVP).
- Phase 2 would add ctDNA monitoring (VAF trends, reversion mutations).
- **Question**: Should I design `ResistanceProphetService` architecture to accommodate ctDNA later (e.g., abstract signal detection interface)?
- Or focus purely on MVP and refactor later?
- Manager's architectural guidance?

---

### **SECTION 5: EDGE CASES & ERROR HANDLING**

**Q15 (Missing CA-125 Data):**
- **Scenario**: Patient has genomic data BUT no CA-125 history.
- **Question**: Should resistance prediction:
  - **Skip CA-125 signal** (use only DNA repair + pathway escape, 2-of-2 detection)
  - **Return "INSUFFICIENT_DATA"** (require all 3 signals)
  - **Lower confidence** (e.g., cap at MEDIUM if only 2 signals available)
- Manager's decision?

**Q16 (Baseline SAE Features Missing):**
- **Scenario**: Patient has current mutations BUT no baseline SAE features (pre-treatment).
- **Question**: Should resistance prediction:
  - **Use population average** as baseline (e.g., DNA repair capacity = 0.5)
  - **Return "INSUFFICIENT_DATA"** (require baseline)
  - **Estimate baseline** from current mutations (reverse-engineer)
- Manager's preference?

**Q17 (False Positive Risk):**
- **Scenario**: Resistance Prophet predicts HIGH risk (≥70%), but patient actually responds well.
- **Question**: How should we handle false positives?
  - **Adjust thresholds** (e.g., require 3-of-3 signals for HIGH urgency)?
  - **Add disclaimer** ("Research use only - discuss with oncologist")?
  - **Track feedback** (oncologist can report false alarms)?
- Manager's risk mitigation strategy?

---

### **SECTION 6: INTEGRATION WITH EXISTING SYSTEMS**

**Q18 (Resistance Playbook Service):**
- Current: `resistance_playbook_service.py` exists (combo strategies + next-line switches).
- Plan: `ResistanceProphetService` calls `ResistancePlaybookService.get_recommendations()`.
- **Question**: Should resistance predictions automatically trigger Resistance Playbook display in Mission Control?
- Or should it be a separate user action ("View Resistance Playbook")?
- Manager's UX guidance?

**Q19 (Treatment Line Intelligence):**
- Current: `treatment_line_service.py` exists (line appropriateness, cross-resistance).
- Plan: `ResistanceProphetService` uses treatment line context for pathway escape detection.
- **Question**: Should Treatment Line Intelligence be consulted:
  - **Before resistance prediction** (filter signals based on line appropriateness)
  - **After resistance prediction** (use predictions to update line recommendations)
  - **Both** (bidirectional integration)
- Manager's integration strategy?

**Q20 (Co-Pilot Conversational Access):**
- Manager said: "Enhance existing Co-Pilot endpoint."
- **Question**: Should resistance predictions be accessible via Co-Pilot natural language queries?
  - Example: "What's my risk of resistance?"
  - Should this call `/api/ayesha/complete_care_v2` with `include_resistance_prediction=true`?
- Manager's conversational UX preference?

---

## 🎯 **MANAGER: PLEASE ANSWER ALL 20 QUESTIONS ABOVE**

**Purpose:** Prevent hallucination, scope creep, and architectural mismatch.

**Format:** Please answer inline (e.g., "Q1: Use ALL 161 patients, stratify by stage for subgroup analysis").

**Timeline:** I will wait for your answers before proceeding with implementation.

**Accountability:** Your answers will be copied into each implementation file as header comments to prevent deviation.

---

## ✅ **MANAGER ANSWERS (APPROVED – Jan 14, 2025)**\n
Q1: Use ALL 161 patients for primary analysis; report subgroup analyses for Stage IIIC+IV. Ensure N per subgroup ≥40; if not, collapse to “advanced (III–IV)” and note power caveat.\n
Q2: Primary = Binary (sensitive vs resistant+refractory). Secondary (exploratory) = 3-class only if per-class N≥30; otherwise skip.\n
Q3: CA‑125 histories are not guaranteed in this cohort. Phase 1 (retrospective): proceed without CA‑125 (signals = DNA repair + pathway escape). Phase 1b (prospective): use CA‑125 kinetics for Ayesha live. If a reliable external CA‑125 source appears, consider Phase 1c.\n
Q4: Metrics for retrospective validation: Primary AUROC (target ≥0.70). Secondary: sensitivity, specificity, PPV/NPV, and calibration curve. Include subgroup metrics for Stage IIIC+IV.\n
\n
Q5: Keep CA‑125 computation in `ca125_intelligence.py`. Keep resistance logic (signal fusion, thresholds, probability) in `ResistanceProphetService` for separation of concerns.\n
Q6: Keep detection logic in `ResistanceProphetService`. If needed, add minimal getters to `sae_feature_service.py` (no complex logic there).\n
Q7: Make resistance prediction opt‑in via `include_resistance_prediction=true`, and compute conditionally (only when minimal inputs present). Default = off to control latency/cost.\n
Q8: Mission Control behavior: show Critical Alerts automatically only when HIGH risk (≥70% and ≥2 signals). Otherwise provide an on‑demand “Check Resistance Risk” action. Also enable proactive toast if probability crosses 70% during session.\n
\n
Q9: Thresholds: HIGH if probability ≥0.70 and ≥2 signals; MEDIUM if 0.50–0.69 or exactly 1 signal; LOW if <0.50. Display signal count with the probability.\n
Q10: HIGH urgency actions: (1) ESCALATE_IMAGING within 1 week, (2) CONSIDER_SWITCH within 2 weeks, (3) REVIEW_RESISTANCE_PLAYBOOK within 1 week. MEDIUM: MONITOR_WEEKLY x4 weeks and re‑assess; LOW: routine monitoring.\n
Q11: Use `ResistancePlaybookService` as the source of truth for next‑line options; consult `TreatmentLineService` for line appropriateness and cross‑resistance. Avoid hardcoded static rules in Prophet.\n
\n
Q12: Timeline is feasible with scope guard: Week 1 (Prophet core without CA‑125), Week 2 (integration + validation on 161), Week 3 (Mission Control wiring). Defer retrospective CA‑125 to Phase 2.\n
Q13: Go/No‑Go: AUROC ≥0.70 overall AND sensitivity ≥0.75 AND specificity ≥0.70 on the 161‑patient set. If unmet, ship as RUO with thresholds raised to minimize false positives and add clear disclaimers.\n
Q14: Design Prophet with a pluggable “signal provider” interface so ctDNA (VAF trends, reversions) can be added in Phase 2 without refactor.\n
\n
Q15: When CA‑125 missing, skip that signal; cap confidence at MEDIUM unless ≥2 non‑CA‑125 signals present. Flag “INSUFFICIENT_CA125_DATA” in provenance.\n
Q16: If baseline SAE missing, use population average baseline (0.50) with a confidence penalty and provenance flag; do not attempt to back‑infer baseline from current mutations.\n
Q17: Mitigate false positives by: RUO disclaimer, clinician feedback loop, optional “strict mode” requiring 3‑of‑3 signals for HIGH, and periodic threshold tuning based on outcomes.\n
\n
Q18: Auto‑open Resistance Playbook card when HIGH risk; otherwise present a “View Resistance Playbook” button. Always log the trigger in provenance.\n
Q19: Consult Treatment Line Intelligence both before (filter/appropriateness) and after (update recommendations) Prophet computation.\n
Q20: Yes. Support Co‑Pilot queries like “What’s my resistance risk?” by calling `/api/ayesha/complete_care_v2?include_resistance_prediction=true`. Provide concise, clinician‑friendly responses with links to details.\n
\n
---

## 🔥 **FOR AYESHA'S LIFE** ⚔️

**Week 1-2: Zo builds Resistance Prophet prediction engine**
- CA-125 kinetics analysis
- DNA repair restoration detection
- Pathway escape detection
- Validated on Jr's 161 patients

**Week 3: Jr integrates into Mission Control dashboard**
- Critical Alerts banner
- Resistance Playbook card
- Next Best Actions
- Complete end-to-end workflow

**Result:** A complete system that predicts resistance → displays insights → enables action

