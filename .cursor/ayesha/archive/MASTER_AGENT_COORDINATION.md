# Master Agent Coordination - All Workstreams

**Date**: January 28, 2025  
**Status**: 🎯 **STRATEGIC ALIGNMENT COMPLETE**  
**Purpose**: Single source of truth for all agent priorities and coordination

---

## Executive Summary

We have **5 active workstreams** that need coordination:

| Priority | Agent | Task | Status | Timeline |
|----------|-------|------|--------|----------|
| 1 | **Zo2** | S/P/E Outcome Calibration | 🔄 Phase 1 Ready | 1-2 weeks |
| 2 | **Trial Matching** | Mechanism Fit + Resistance | 📋 Approved | 3 days |
| 3 | **Food Validation** | Nutritional Support | ⚠️ Alignment Needed | 1 day |
| 4 | **Zo** | MBD4+TP53 Analysis | ✅ Complete | Done |
| 5 | **MM Doctrines** | TCell/TCF1/TLS/TOX | ⏸️ On Hold | After 1-3 |

---

## The Strategic Picture

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
│        ├──────────────────────────────────────────────────────────── │
│        │                                                             │
│        ├──→ Trial Matching (Agent 2) ───→ Ranked Trials              │
│        │                                                             │
│        ├──→ Resistance Prediction (Agent 2) ───→ Risk Signals        │
│        │                                                             │
│        ├──→ Food Recommendations (Agent 3) ───→ Nutritional Support  │
│        │                                                             │
│        └──→ MM Doctrines (Agent 4) ───→ TCell/TLS/TOX ⏸️ ON HOLD    │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Priority 1: Zo2 - S/P/E Outcome Calibration

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

**Document**: `zo2.mdc` (916 lines)

---

## Priority 2: Trial Matching Agent

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

**Document**: `FINAL_COMPREHENSIVE_PLAN_ALIGNED.md` (git worktree)

---

## Priority 3: Food Validation Agent

**Goal**: Align terminology and add honest framing

**Status**: ⚠️ Alignment Needed (95% complete)

**What to Do**:
1. Rename `efficacy_score` → `alignment_score`
2. Add limitations section to documentation
3. Revise "precision medicine-grade" claim
4. Add frontend disclaimer

**Time**: 1 day

**Key Principle**: Same rules as drugs - mechanism alignment, not outcome prediction.

**Document**: `MASTER_SUMMARY_FOR_MANAGER.md`

---

## Priority 4: Zo - MBD4+TP53 Analysis

**Goal**: Mechanism alignment for rare case

**Status**: ✅ Complete

**What Was Delivered**:
- Pathway vulnerabilities identified (DDR: 1.0, TP53: 0.8)
- Mechanism-aligned drugs found (PARP #1-3, Platinum #4)
- Mechanism vectors computed (DDR: 1.4)
- 6 P0 blockers fixed

**Key Principle**: Mechanism alignment works. Outcome prediction is separate (Zo2).

**Document**: `mbd4-tp53-hgsoc-analysis-2e21acc4.plan.md`

---

## Priority 5: MM Doctrines Agent

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

**Corrections Applied**:
- Mechanism vector service EXISTS (not empty)
- Pathway score extraction EXISTS
- Phase 0 is mostly done

**Document**: `MM_DOCTRINES_EXECUTION_DIRECTIVE.mdc`

---

## Coordination Rules

### Rule 1: Sequential, Not Parallel (for now)

Don't start MM Doctrines until Zo2/Trial Matching/Food Validation are done.

### Rule 2: Checkpoint Before Expanding

Every agent must checkpoint before expanding scope:
- Zo2: r > 0.10 before Phase 2
- Trial Matching: MVP thresholds before enhancement
- MM Doctrines: First classifier before second

### Rule 3: Same Language

All agents use consistent terminology:
- ✅ "Mechanism alignment score"
- ❌ "Efficacy score" (implies outcome prediction)
- ✅ "Biological plausibility"
- ❌ "Clinical validation" (unless actually validated)

### Rule 4: Honest Framing

All agents must add limitations:
- "These are mechanism-aligned recommendations"
- "Not validated for outcome prediction"
- "Decision support, not decision making"

---

## This Week's Actions

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

## Strategic Reviews Created

| Document | Purpose |
|----------|---------|
| `STRATEGIC_DIRECTION.md` | Master strategy for all agents |
| `STRATEGIC_REVIEW_MECHANISM_TRIAL_MATCHING.md` | Trial matching review |
| `STRATEGIC_REVIEW_FOOD_VALIDATION.md` | Food validation review |
| `STRATEGIC_REVIEW_MM_DOCTRINES.md` | MM doctrines review |
| `MASTER_AGENT_COORDINATION.md` | This document |

---

## Key Decisions Made

### For Zo2

- ✅ Validate Phase 1 first, checkpoint before Phase 2
- ✅ Realistic targets: r=0.15-0.25 (not r=0.50)
- ✅ Independent of other agents

### For Trial Matching

- ✅ Approved to start independently
- ✅ MVP thresholds: 0.70/0.65/0.65 (not 0.80/0.75/0.70)
- ✅ Add MoA coverage report

### For Food Validation

- ✅ Rename efficacy_score → alignment_score
- ✅ Add limitations section
- ✅ Same rules as drugs

### For MM Doctrines

- ⏸️ On hold until other agents done
- ⏸️ Need MVP definition (2-week deliverable)
- ⏸️ Start with ONE proxy classifier
- ✅ Mechanism vector service already exists

---

## The Righteous Path

1. **Finish what's started** before starting new things
2. **Validate before optimizing** - checkpoint at each phase
3. **Be honest about limitations** - mechanism alignment ≠ outcome prediction
4. **Sequence properly** - don't start 6-month projects when 1-week projects are unfinished

---

**Status**: 🎯 **ALL AGENTS HAVE CLEAR DIRECTION**

- Zo: ✅ Done
- Zo2: 🔄 Execute Phase 1
- Trial Matching: 📋 Start Day 1
- Food Validation: ⚠️ Make alignment changes
- MM Doctrines: ⏸️ Wait for others

