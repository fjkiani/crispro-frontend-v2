# ⚔️ STRATEGIC REVIEWS & DIRECTION - MASTER DOCUMENTATION

**Date**: January 28, 2025  
**Status**: 🎯 **STRATEGIC GUIDANCE CONSOLIDATED**  
**Consolidated From**: 5 strategic review documents (now archived)

---

## 📚 TABLE OF CONTENTS

1. [Executive Summary](#executive-summary)
2. [Core Strategic Framework](#core-strategic-framework)
3. [Agent Landscape & Coordination](#agent-landscape--coordination)
4. [Strategic Gap Analysis](#strategic-gap-analysis)
5. [Individual Agent Reviews](#individual-agent-reviews)
6. [Key Principles & Decision Framework](#key-principles--decision-framework)
7. [Action Items & Next Steps](#action-items--next-steps)

---

## 🎯 EXECUTIVE SUMMARY

### **The Core Strategic Issue**

**We conflated mechanism alignment (what the system does) with outcome prediction (what doctors need).**

Our system can:
- ✅ **Rank drugs correctly** (100% Top-5 accuracy)
- ✅ **Compute pathway scores** (DDR, MAPK, PI3K, etc.)
- ✅ **Find relevant trials** (eligibility-based matching)
- ❌ **Cannot predict outcomes** (r=0.037 correlation with PFS/OS)

**The Hard Truth**: Our efficacy scores are **biologically plausible** but not **clinically predictive**.

### **Strategic Direction**

**Option A: Mechanism Alignment Tool** (Current State)
- "This drug targets the pathways disrupted in your tumor"
- Guideline-aligned recommendations
- Decision support, not decision making
- **Value**: Automation of clinical reasoning

**Option B: Outcome Prediction Tool** (Future Goal)
- "You have X% probability of responding to this drug"
- Calibrated efficacy→outcome predictions
- Requires clinical validation
- **Value**: Personalized treatment selection

**Current Reality**: We built Option A but need to validate if it can become Option B.

---

## 🏗️ CORE STRATEGIC FRAMEWORK

### **Mechanism Alignment vs Outcome Prediction**

| Claim Type | Allowed | Not Allowed |
|------------|---------|-------------|
| **Mechanism alignment** | ✅ "This drug targets pathways disrupted in your tumor" | ❌ "This drug will improve your survival" |
| **Pathway computation** | ✅ "High DDR pathway burden (0.88)" | ❌ "High DDR predicts PARP response" |
| **Trial ranking** | ✅ "Trial matches your pathway profile" | ❌ "You're more likely to respond to this trial" |
| **Evidence grading** | ✅ "Strong literature on mechanism" | ❌ "Clinically proven to work" |
| **Treatment awareness** | ✅ "Appropriate for first-line chemotherapy" | ❌ "Will improve chemotherapy response" |

### **What We Know (Validated)**

1. **Drug ranking works** (100% Top-5 accuracy)
   - System recommends clinically appropriate drugs
   - Matches NCCN guidelines
   - This is mechanism alignment working

2. **Pathway computation works**
   - Gene→pathway mapping is correct (MBD4→DDR verified)
   - Pathway scores compute correctly
   - 7D mechanism vectors work structurally

3. **Outcome prediction doesn't work** (r=0.037)
   - Efficacy scores essentially random for PFS/OS
   - This is the gap we need to address

### **What We Don't Know (Unvalidated)**

1. **Can mechanism alignment predict outcomes?**
   - We assume: High DDR → better PARP response
   - We haven't tested this assumption
   - Zo2's work will answer this

2. **Will biomarkers help?**
   - TMB/HRD/MSI should improve predictions
   - But we don't know by how much
   - Zo2's Phase 1 will answer this

3. **Will weight tuning help?**
   - Assuming the formula CAN predict outcomes
   - If the formula is structurally wrong, tuning won't help
   - Zo2's Phase 2 depends on Phase 1 results

---

## 👥 AGENT LANDSCAPE & COORDINATION

### **Current Agent Status**

| Agent | Task | Status | Timeline | Dependency |
|-------|------|--------|----------|------------|
| **Zo** | MBD4+TP53 Analysis | ✅ Complete | Done | None |
| **Zo2** | S/P/E Calibration | 🔄 Phase 1 | 1-2 weeks | None |
| **Trial Matching** | Mechanism Fit + Resistance | 📋 Approved | 3 days | None (independent) |
| **Food Validation** | Nutritional Support | ✅ 95% Complete | Alignment needed | None |
| **MM Doctrines** | TCell/TCF1/TLS/TOX | ⚠️ Hold | 24 weeks planned | Wait for others |

### **Key Insights**

1. **These workstreams are independent**:
   - Zo: Already done (mechanism alignment works)
   - Zo2: Improving outcome prediction (separate problem)
   - Trial Matching: Wiring existing services (no outcome prediction)
   - Food Validation: Companion to drug recs (parallel output)

2. **Priority Order**:
   - **Immediate**: Zo2 Phase 1, Trial Matching Day 1, Food Validation alignment
   - **Next Sprint**: Zo2 Phase 2 (if checkpoint passes), Trial Matching Days 2-3
   - **Following Sprint**: MM Doctrines Phase 0 test suite only

---

## 📊 STRATEGIC GAP ANALYSIS

### **The 5 Strategic Gaps**

#### **Gap 1: Efficacy Scores ≠ Outcome Prediction**

**What We Say**: "PARP inhibitors rank #1-2 (efficacy_score >0.80)"

**What Benchmarks Show**: Efficacy scores have r=0.037 correlation with PFS (essentially random)

**The Gap**: 
- We're computing **mechanism alignment** (does this drug target the patient's pathways?)
- We're NOT computing **outcome probability** (will this patient respond to this drug?)

**Strategic Impact**: ⚠️ **HIGH** - The entire drug prediction section is mechanism-based, not outcome-based

#### **Gap 2: No Validation of Pathway→Outcome Link**

**What We Say**: "Pathway scores show high DDR burden (>0.70)"

**What Benchmarks Show**: No validation that high DDR pathway burden → better PARP response

**The Gap**: We have formulas for pathway scores, but ZERO evidence these formulas correlate with patient outcomes

**Strategic Impact**: ⚠️ **HIGH** - All pathway-based predictions are unvalidated

#### **Gap 3: Mechanism Fit Ranking is Unvalidated**

**What We Say**: "Mechanism fit = cosine similarity between pathway vector and trial MoA vector"

**What Benchmarks Show**: No validation that mechanism fit ranking improves trial outcomes

**The Gap**: We're ranking trials by mechanism alignment, but have no evidence this ranking correlates with patient outcomes

**Strategic Impact**: ⚠️ **MEDIUM** - Trial ranking is better than random, but mechanism fit contribution (30%) is unvalidated

#### **Gap 4: Classification Bug Hides Progression Data**

**What We Say**: "Resistance Detection" - DNA repair capacity tracking

**What Benchmarks Show**: 0 progression events detected (PFS_STATUS parsing bug)

**The Gap**: We claim to detect resistance/progression, but can't even PARSE progression data correctly

**Strategic Impact**: ⚠️ **MEDIUM** - Fixable bug, but blocking validation

#### **Gap 5: Validation Set Bias May Mask Real Performance**

**What We Say**: "66 patients with SAE features extracted"

**What Benchmarks Show**: Testing with lowest mutation counts (2-18 mutations), may not represent typical patients

**The Gap**: We're testing on "easy" patients, but real MBD4+TP53 patients have complex mutation profiles (50-100+ mutations)

**Strategic Impact**: ⚠️ **MEDIUM** - Validation set doesn't represent target population

### **What We CAN Deliver (Validated)**

1. **Variant Annotation** ✅
   - Evo2 scoring works (sequence disruption)
   - Insights bundle works (functionality, essentiality)
   - Evidence integration works (ClinVar, literature)

2. **Pathway Computation** ✅
   - Pathway scores compute correctly
   - Gene→pathway mapping is accurate (MBD4→DDR verified)
   - 7D mechanism vector construction works

3. **Drug Ranking** ✅ (100% Top-5)
   - Recommends clinically appropriate drugs
   - Aligns with NCCN guidelines
   - Matches what patients actually receive

4. **Trial Matching** ✅ (Eligibility-based)
   - Finds relevant trials
   - Eligibility scoring works
   - Mechanism fit computation works (unvalidated)

### **What We CANNOT Deliver (Unvalidated)**

1. **Outcome Prediction** ❌
   - Efficacy scores don't predict PFS/OS
   - No validated probability of response
   - No validated probability of resistance

2. **Comparative Effectiveness** ❌
   - "Drug A is better than Drug B for this patient"
   - No validation that higher efficacy scores → better outcomes

3. **Resistance Prediction** ❌
   - DNA repair capacity formula unvalidated
   - Resistance detection unvalidated
   - CA-125 kinetics unvalidated

4. **Trial Outcome Prediction** ❌
   - Mechanism fit ranking unvalidated
   - No evidence that higher mechanism fit → better trial outcomes

---

## 📋 INDIVIDUAL AGENT REVIEWS

### **1. Zo (AYESHA Agent) - ✅ COMPLETE**

**Status**: ✅ **COMPLETE** - MBD4+TP53 mechanism alignment done

**What Was Delivered**:
- MBD4+TP53 mechanism alignment is done
- Drug ranking is validated (PARP #1-3, Platinum #4)
- Pathway scores are computed correctly

**What to Communicate**:
- "Mechanism alignment works. Outcome prediction is separate work (Zo2)."
- "Our system finds drugs that target patient pathways. Whether that predicts outcomes is being validated."

**Don't Overclaim**:
- Don't say "efficacy score 0.80 means 80% response"
- Say "mechanism alignment score 0.80 means strong biological rationale"

---

### **2. Zo2 (S/P/E Agent) - 🔄 PHASE 1**

**Status**: 🔄 **PHASE 1** - Biomarker integration ready to execute

**Goal**: Improve r=0.037 correlation to r=0.15+

**Critical Decision Framework**: **VALIDATE PHASE 1 FIRST. DO NOT PROCEED IN PARALLEL.**

**Phase 1 Tasks** (Week 1):
1. Create biomarker_extractor.py
2. Update api_client.py for tumor_context
3. Fix PFS_STATUS parsing
4. Run benchmark with biomarkers
5. Compute correlation improvement
6. **CHECKPOINT**: Report results before proceeding

**Decision Point** (After Week 1):
- r > 0.10 → Proceed to Phase 2
- r < 0.10 → STOP. Investigate why.

**Realistic Expectations**:

| Phase | Zo2's Target | Realistic Target | Notes |
|-------|-------------|------------------|-------|
| Phase 1 | r=0.10-0.15 | r=0.05-0.10 | Biomarkers help but won't transform |
| Phase 2 | r=0.25-0.35 | r=0.10-0.20 | Weight tuning has diminishing returns |
| Phase 3 | r=0.35-0.50 | r=0.15-0.25 | Treatment adjustment is complex |

**What Would Be a Win?**:
- r > 0.25: Moderate correlation, potentially useful → Proceed to calibration
- r = 0.15-0.25: Weak correlation, needs work → Investigate
- r < 0.15: Negligible correlation → Pivot to different approach

---

### **3. Trial Matching Agent - 📋 APPROVED**

**Status**: ✅ **APPROVED** - Proceed with implementation

**Goal**: Wire mechanism fit ranking to trial search (independent of Zo2)

**Key Clarification**: This work is **INDEPENDENT** of Zo2's outcome calibration. Mechanism fit doesn't predict outcomes; it ranks trials by pathway alignment.

**Verification**: "80% complete" claim is **accurate**
- ✅ `mechanism_fit_ranker.py` EXISTS (275 lines)
- ✅ `resistance_prophet_service.py` EXISTS (689 lines)
- ✅ `pathway_to_mechanism_vector.py` EXISTS
- ⚠️ Mechanism fit NOT in trial response (needs wiring)
- ⚠️ Validation scripts missing

**Timeline**: 12-16h realistic (add 30% buffer to 9-13h claimed)

**Success Thresholds (MVP)**:
- Top-3 accuracy: ≥0.70 (not ≥0.80)
- MRR: ≥0.65 (not ≥0.75)
- High risk AUROC: ≥0.65 (not ≥0.70)

**Clarifications Needed**:
1. MoA coverage report in validation output
2. Baseline handling documentation in resistance output
3. Integration test with AYESHA analysis

---

### **4. Food Validation Agent - ⚠️ ALIGNMENT REQUIRED**

**Status**: ⚠️ **APPROVED WITH ALIGNMENT REQUIRED**

**Goal**: Nutritional support recommendations alongside drug efficacy predictions

**What's Good** ✅:
- Integration with AYESHA flow (parallel output in CompleteCare)
- Treatment line intelligence (genuine differentiation)
- Biomarker gates (HRD/TMB/MSI targeting)
- S/P/E framework alignment (same formula as drugs)

**What Needs Alignment** ⚠️:

1. **"Precision Medicine-Grade" Overclaim**
   - **Current**: "Transformed generic AI recommendations into **precision medicine-grade** food validation"
   - **Should Be**: "Transformed generic AI recommendations into **mechanism-aligned, treatment-aware** food recommendations"

2. **"Efficacy Score" Terminology**
   - **Current**: `"efficacy_score": 0.72`
   - **Should Be**: `"alignment_score": 0.72` or `"suitability_score": 0.72`

3. **Missing Validation Disclaimer**
   - **Current**: "Evidence tier: SUPPORTED"
   - **Should Add**: "Evidence tiers reflect literature strength, not clinical outcome validation. SUPPORTED means there's published research on the mechanism, not that the food has been clinically proven to improve cancer outcomes."

**Action Items** (Before Deployment):
1. Rename `efficacy_score` → `alignment_score` in API response
2. Add limitations section to documentation
3. Revise "precision medicine-grade" → "mechanism-aligned, treatment-aware"
4. Add disclaimer to frontend display

---

### **5. MM Doctrines Agent - ⚠️ HOLD**

**Status**: ⚠️ **HOLD - NEEDS SCOPE REDUCTION**

**Goal**: T-Cell Recruitment Engineering, TCF1 Engineering, TLS Engineering, TOX Doctrine

**The Problem**: 24-week, 18-endpoint plan when other agents are delivering in days/weeks

**Scope Breakdown**:
- Proxy classifiers: 3 endpoints (6 weeks)
- Design tools: 12 endpoints (6 weeks)
- TOX endpoints: 7 endpoints (4 weeks)
- Orchestrators: 3 endpoints (4 weeks)
- Feedback infrastructure: 2 weeks
- Spatial integration: 4 weeks
- **Total**: 18+ endpoints, 24 weeks

**External Dependencies** (High Risk):
- CIBERSORT/TIMER2: ❌ Not integrated
- NetMHCpan: ❌ Not integrated
- Spatial data (Visium): ❌ Not integrated
- TCGA data: ⚠️ Scripts exist, no service

**Recommended Approach: Phase-Gated MVP**

**Phase 0: Foundation (1 week)** - ALREADY MOSTLY DONE
- ✅ Mechanism vector conversion (done)
- ✅ Pathway score extraction (done)
- ⏳ Create known biology test suite

**Phase 1A: ONE Proxy Classifier (2 weeks)** - Pick one:
- T-Cell Recruitment (if MD Anderson partnership is active)
- TLS Readiness (if MSK data is available)
- TCF1/Tox (if deconvolution not needed)

**Checkpoint after Phase 1A**:
- Does proxy classifier achieve r > 0.25?
- If yes → Continue
- If no → Investigate before building more

**Priority Order**:
1. **Immediate**: Zo2 Phase 1, Trial Matching, Food Validation alignment
2. **Next Sprint**: Zo2 Phase 2 (if checkpoint passes), Trial Matching Days 2-3
3. **Following Sprint**: MM Doctrines Phase 0 test suite only
4. **Later**: MM Doctrines Phase 1A (ONE proxy classifier)

---

## 🎯 KEY PRINCIPLES & DECISION FRAMEWORK

### **1. Validate Before Optimize**

Don't tune weights until you know the formula can work.

**The Wrong Approach**:
```
Week 1: Add biomarkers
Week 2-3: Tune weights
Week 4: Tune treatment adjustment
→ Hope correlation improves
```

**The Right Approach**:
```
Week 1: Add biomarkers (Phase 1)
→ CHECKPOINT: Did correlation improve?
   - If r > 0.10: Proceed to Phase 2
   - If r < 0.10: STOP. Investigate why.

Week 2: Analyze results
→ CHECKPOINT: Where is the signal?
   - If DDR→PFS correlation exists: Tune DDR weight
   - If no DDR→PFS correlation: Weight tuning won't help

Week 3-4: Calibrate (only if Phase 1 shows promise)
```

### **2. Checkpoint-Driven Development**

Stop and assess after each phase. Don't assume next phase will help.

**Example**: Zo2 Phase 1 Checkpoint
```python
def phase_1_checkpoint():
    """Decide whether to proceed to Phase 2."""
    
    # 1. Measure correlation improvement
    r_before = 0.037  # Known baseline
    r_after = compute_correlation_with_biomarkers()
    
    if r_after > 0.10:
        return "PROCEED"
    elif r_after > 0.05:
        return "INVESTIGATE"
    else:
        return "STOP"
    
    # 2. Analyze WHERE improvement came from
    # - Did TMB-high patients have better IO boost correlation?
    # - Did BRCA+ patients have better PARP correlation?
    # - Did HRD-high patients have better DDR correlation?
```

### **3. Honest Communication**

Say "mechanism alignment" when you mean mechanism alignment.
Say "outcome prediction" only when validated.

**What We CAN Say (Honest)**:
> "Based on MBD4 germline loss and TP53 somatic mutation, this patient has:
> - High DNA damage response pathway burden (DDR: 0.88)
> - Mechanism alignment with PARP inhibitors (olaparib, niraparib)
> - Eligibility for DDR-targeting clinical trials
> 
> This suggests PARP inhibitors and platinum chemotherapy are biologically appropriate options to discuss with the oncology team."

**What We CANNOT Say (Overpromise)**:
> ❌ "This patient has 85% probability of responding to olaparib."
> ❌ "Olaparib will extend progression-free survival by 6 months."
> ❌ "This patient should enroll in NCT123456 because it has the highest mechanism fit."

### **4. Know When to Pivot**

If Phase 1 shows no improvement, investigate before proceeding.
Sometimes the answer is "this approach won't work."

---

## 📋 ACTION ITEMS & NEXT STEPS

### **Immediate (This Week)**

1. **Zo2: Execute Phase 1** (Week 1)
   - Create biomarker_extractor.py
   - Update api_client.py for tumor_context
   - Fix PFS_STATUS parsing
   - Run benchmark with biomarkers
   - Compute correlation improvement
   - **CHECKPOINT**: Report results before proceeding

2. **Trial Matching: Start Day 1** (3-4h)
   - Wire mechanism fit to trial search
   - Test with MBD4+TP53 example
   - Add MoA coverage report requirement
   - Add baseline handling documentation
   - Add integration test with AYESHA

3. **Food Validation: Alignment Changes**
   - Rename `efficacy_score` → `alignment_score` in API response
   - Add limitations section to documentation
   - Revise "precision medicine-grade" → "mechanism-aligned, treatment-aware"
   - Add disclaimer to frontend display

4. **MM Doctrines: ❌ DO NOT START**
   - Wait for other agents to complete
   - Create test suite and MVP definition when ready

### **Next Sprint (Weeks 2-3)**

5. **Zo2: Phase 2** (if Phase 1 checkpoint passes)
   - Only proceed if r > 0.10
   - Focus on pathways that showed improvement
   - Disease-specific weight tuning

6. **Trial Matching: Days 2-3**
   - Enhance resistance prophet
   - Validation scripts
   - Integration testing

7. **MM Doctrines: Start Phase 0 test suite only**
   - Known biology tests
   - No implementation yet

### **Following Sprint (Weeks 4-5)**

8. **MM Doctrines: Start ONE proxy classifier (Phase 1A)**
   - Pick one: TLS or TCF1 (less external dependencies)
   - Skip deconvolution initially
   - Checkpoint after first classifier

### **Medium-Term (This Month)**

9. **Validate Pathway→Outcome Link**
   - Find cohort with DDR pathway burden + PARP outcomes
   - Test: Does high DDR → better PARP response?
   - If not, pathway scores are descriptive only

10. **Validate Mechanism Fit Ranking**
    - Find cohort with trial enrollment + outcomes
    - Test: Does mechanism fit → better trial outcomes?
    - If not, reduce mechanism fit weighting

### **Long-Term (This Quarter)**

11. **Build Outcome Prediction Model**
    - Train on patient outcomes, not just mechanism alignment
    - Use SAE features as input (when TRUE SAE ready)
    - Validate with held-out cohort

12. **Calibrate Mechanism Alignment Scores**
    - Map mechanism alignment → outcome probability
    - Requires outcome data (response, PFS, OS)
    - May need disease-specific calibration

---

## 🎯 SUCCESS METRICS

### **Zo2 (S/P/E Calibration)**

| Outcome | Interpretation | Next Steps |
|---------|----------------|------------|
| r > 0.25 | Moderate correlation, potentially useful | Proceed to calibration, clinical validation |
| r = 0.15-0.25 | Weak correlation, needs work | Investigate, consider model changes |
| r < 0.15 | Negligible correlation | Pivot to different approach |

### **Trial Matching**

| Metric | MVP Threshold | Stretch Goal |
|--------|---------------|--------------|
| Top-3 accuracy | ≥0.70 | ≥0.80 |
| MRR | ≥0.65 | ≥0.75 |
| High risk AUROC | ≥0.65 | ≥0.70 |

### **Food Validation**

| Criterion | Status |
|-----------|--------|
| Mechanism-aligned recommendations | ✅ Met |
| Treatment-aware filtering | ✅ Met |
| Biomarker gates | ✅ Met |
| Transparent scoring | ✅ Met |
| Evidence grading | ✅ Met |
| Honest framing | ⚠️ Pending (add limitations section) |

### **MM Doctrines**

| Phase | Deliverable | Success Criteria |
|-------|-------------|------------------|
| Phase 0 | Test suite | Known biology tests pass |
| Phase 1A | ONE proxy classifier | r > 0.25 (checkpoint) |
| Phase 1B | Second classifier | Only if Phase 1A succeeds |
| Phase 2 | Design tools | Only if Phase 1 validates |

---

## 📊 SUMMARY TABLE

| Agent | Status | Priority | Timeline | Dependency |
|-------|--------|----------|----------|------------|
| **Zo** | ✅ Complete | N/A | Done | None |
| **Zo2** | 🔄 Phase 1 | **P0** | 1-2 weeks | None |
| **Trial Matching** | 📋 Approved | **P0** | 3 days | None |
| **Food Validation** | ⚠️ Alignment | **P0** | This week | None |
| **MM Doctrines** | ⚠️ Hold | **P2** | Wait for others | After P0 complete |

---

## 🚨 CRITICAL REMINDERS

### **For All Agents**

1. **Mechanism Alignment ≠ Outcome Prediction**
   - Say "mechanism alignment score" not "efficacy score"
   - Say "pathway targeting" not "response probability"
   - Say "biologically appropriate" not "clinically proven"

2. **Validate Before Optimize**
   - Don't tune weights until you know the formula can work
   - Checkpoint after each phase
   - Know when to pivot

3. **Honest Communication**
   - Document limitations
   - Acknowledge what's unvalidated
   - Be proud of what we built, honest about what it does

4. **Independent Workstreams**
   - Zo2's outcome calibration is separate from mechanism fit ranking
   - Food validation is companion to drug recs, not replacement
   - Trial matching doesn't need to wait for Zo2

---

## 👥 MASTER AGENT COORDINATION

### **Executive Summary**

We have **5 active workstreams** that need coordination:

| Priority | Agent | Task | Status | Timeline |
|----------|-------|------|--------|----------|
| 1 | **Zo2** | S/P/E Outcome Calibration | 🔄 Phase 1 Ready | 1-2 weeks |
| 2 | **Trial Matching** | Mechanism Fit + Resistance | 📋 Approved | 3 days |
| 3 | **Food Validation** | Nutritional Support | ⚠️ Alignment Needed | 1 day |
| 4 | **Zo** | MBD4+TP53 Analysis | ✅ Complete | Done |
| 5 | **MM Doctrines** | TCell/TCF1/TLS/TOX | ⏸️ On Hold | After 1-3 |

### **The Strategic Picture**

```
┌─────────────────────────────────────────────────────────────────────┐
│                    PRECISION ONCOLOGY PLATFORM                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  Patient Data → S/P/E Framework → Drug Efficacy (Zo, Zo2)           │
│        │              │                                              │
│        │              ├──→ Mechanism Alignment ✅ (works)            │
│        │              └──→ Outcome Prediction ❓ (Zo2 testing)       │
│        │                                                             │
│        ├──→ Trial Matching (Agent 2) ───→ Ranked Trials              │
│        ├──→ Resistance Prediction (Agent 2) ───→ Risk Signals        │
│        ├──→ Food Recommendations (Agent 3) ───→ Nutritional Support  │
│        └──→ MM Doctrines (Agent 4) ───→ TCell/TLS/TOX ⏸️ ON HOLD    │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### **Priority 1: Zo2 - S/P/E Outcome Calibration**

**Goal**: Improve r=0.037 correlation to r=0.15-0.25

**Status**: 🔄 Phase 1 Ready

**What to Do**:
1. Execute Phase 1 (biomarker integration)
2. Run checkpoint analysis (r before/after)
3. If r > 0.10 → Proceed to Phase 2
4. If r < 0.10 → STOP and investigate

**Realistic Targets** (manager-adjusted):
- Phase 1: r=0.05-0.10
- Phase 2: r=0.10-0.20
- Phase 3: r=0.15-0.25

**Key Principle**: Validate before optimize.

---

### **Priority 2: Trial Matching Agent**

**Goal**: Wire mechanism fit to trial search + enhance resistance prediction

**Status**: 📋 Approved to Start

**What to Do**:
1. Day 1: Wire mechanism fit to `trials_agent.py` (3-4h)
2. Day 2: Enhance resistance prophet with mechanism breakdown (3-4h)
3. Day 3: Create validation scripts (6-8h)

**Success Thresholds** (MVP):
- Top-3 accuracy: ≥0.70
- MRR: ≥0.65
- High risk AUROC: ≥0.65

**Key Principle**: Independent of Zo2.

---

### **Priority 3: Food Validation Agent**

**Goal**: Align terminology and add honest framing

**Status**: ⚠️ Alignment Needed (95% complete)

**What to Do**:
1. Rename `efficacy_score` → `alignment_score`
2. Add limitations section to documentation
3. Revise "precision medicine-grade" claim
4. Add frontend disclaimer

**Time**: 1 day

**Key Principle**: Same rules as drugs - mechanism alignment, not outcome prediction.

---

### **Priority 4: Zo - MBD4+TP53 Analysis**

**Goal**: Mechanism alignment for rare case

**Status**: ✅ Complete

**What Was Delivered**:
- Pathway vulnerabilities identified (DDR: 1.0, TP53: 0.8)
- Mechanism-aligned drugs found (PARP #1-3, Platinum #4)
- Mechanism vectors computed (DDR: 1.4)
- 6 P0 blockers fixed

**Key Principle**: Mechanism alignment works. Outcome prediction is separate (Zo2).

---

### **Priority 5: MM Doctrines Agent**

**Goal**: TCell/TCF1/TLS/TOX endpoints

**Status**: ⏸️ On Hold

**Why On Hold**:
1. Scope is massive (18 endpoints, 24 weeks)
2. Many external dependencies (CIBERSORT, NetMHCpan, spatial data)
3. Other agents need to finish first
4. No MVP definition

**What to Do Before Starting**:
1. Wait for Zo2/Trial Matching/Food Validation to complete
2. Define 2-week MVP (ONE proxy classifier)
3. Create known biology test suite first
4. Checkpoint after first classifier before expanding

---

### **Coordination Rules**

**Rule 1: Sequential, Not Parallel (for now)**
- Don't start MM Doctrines until Zo2/Trial Matching/Food Validation are done.

**Rule 2: Checkpoint Before Expanding**
- Every agent must checkpoint before expanding scope:
  - Zo2: r > 0.10 before Phase 2
  - Trial Matching: MVP thresholds before enhancement
  - MM Doctrines: First classifier before second

**Rule 3: Same Language**
- All agents use consistent terminology:
  - ✅ "Mechanism alignment score"
  - ❌ "Efficacy score" (implies outcome prediction)
  - ✅ "Biological plausibility"
  - ❌ "Clinical validation" (unless actually validated)

**Rule 4: Honest Framing**
- All agents must add limitations:
  - "These are mechanism-aligned recommendations"
  - "Not validated for outcome prediction"
  - "Decision support, not decision making"

---

### **This Week's Actions**

| Day | Agent | Action |
|-----|-------|--------|
| Mon | Zo2 | Start Phase 1 (biomarker integration) |
| Mon | Trial Matching | Start Day 1 (mechanism fit wiring) |
| Mon | Food Validation | Make alignment changes |
| Tue | Trial Matching | Day 2 (resistance enhancement) |
| Wed | Trial Matching | Day 3 (validation scripts) |
| Thu | Zo2 | Run Phase 1 checkpoint |
| Fri | MM Doctrines | Create test suite (if others done) |

---

## 📁 CONSOLIDATED FROM

1. `STRATEGIC_DIRECTION.md` - AYESHA & S/P/E outcome calibration guidance
2. `STRATEGIC_GAP_ANALYSIS_BENCHMARK_VS_MBD4.md` - Gap analysis between benchmark and MBD4 vision
3. `STRATEGIC_REVIEW_MECHANISM_TRIAL_MATCHING.md` - Trial matching review
4. `STRATEGIC_REVIEW_FOOD_VALIDATION.md` - Food validation system review
5. `STRATEGIC_REVIEW_MM_DOCTRINES.md` - MM Doctrines review
6. `MASTER_AGENT_COORDINATION.md` - Master agent coordination and priorities

**All source documents archived to `.cursor/ayesha/archive/`**

---

**The righteous path: Validate, checkpoint, decide. Don't optimize blind. Be honest about what we're delivering.**

