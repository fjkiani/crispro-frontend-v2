# Strategic Review: MM Doctrines Execution Directive

**Date**: January 28, 2025  
**Reviewer**: Manager  
**Document Reviewed**: `MM_DOCTRINES_EXECUTION_DIRECTIVE.mdc` (1089 lines)  
**Status**: ⚠️ **NEEDS SCOPE REDUCTION & PRIORITIZATION**

---

## Executive Summary

This is an **ambitious 24-week, 18-endpoint plan** covering four MM doctrines:
- T-Cell Recruitment Engineering (TCell)
- TCF1 Engineering
- TLS Engineering
- TOX Doctrine (Autoimmune Blueprint)

**Verdict**: ⚠️ **GOOD RESEARCH, BUT SCOPE IS TOO LARGE FOR IMMEDIATE EXECUTION**

The plan is well-researched and thorough, but it's a 6-month roadmap when we have 4 other agents already working on nearer-term deliverables. This needs prioritization.

---

## Agent Landscape (Updated with MM Agent)

| Agent | Task | Status | Timeline |
|-------|------|--------|----------|
| **Zo** | MBD4+TP53 Analysis | ✅ Complete | Done |
| **Zo2** | S/P/E Calibration | 🔄 Phase 1 | 1-2 weeks |
| **Trial Matching** | Mechanism Fit + Resistance | 📋 Approved | 3 days |
| **Food Validation** | Nutritional Support | ✅ 95% Complete | Alignment needed |
| **MM Doctrines** | TCell/TCF1/TLS/TOX | 📋 **NEEDS PRIORITIZATION** | 24 weeks planned |

**Key Issue**: MM Doctrines is a 24-week plan while other agents are delivering in days/weeks.

---

## What's Good ✅

### 1. Thorough Research

The plan integrates learnings from four proven analysis plans:
- MBD4+TP53 Analysis Plan
- SOTA Benchmarks Plan
- Ayesha Universalization Plan
- Proxy SAE Validation Plan

### 2. Code-Verified References

The agent verified actual code state:
- ✅ `pathway_to_mechanism_vector.py` EXISTS and is fully implemented (317 lines)
- ✅ `mechanism_fit_ranker.py` EXISTS
- ✅ `sae_feature_service.py` EXISTS with 7D vector structure
- ❌ No MM doctrine endpoints exist yet

### 3. WHY/WHAT/HOW/WHERE Structure

Each doctrine has clear:
- WHY: The clinical problem
- WHAT: What we're building
- HOW: Technical implementation
- WHERE: Code references

### 4. Manager Questions Identified

7 critical questions (Q1-Q7) identified before implementation.

---

## What's Problematic ⚠️

### Problem 1: Scope is Massive

| Component | Count | Est. Time |
|-----------|-------|-----------|
| Proxy classifiers | 3 endpoints | 6 weeks |
| Design tools | 12 endpoints | 6 weeks |
| TOX endpoints | 7 endpoints | 4 weeks |
| Orchestrators | 3 endpoints | 4 weeks |
| Feedback infrastructure | N/A | 2 weeks |
| Spatial integration | N/A | 4 weeks |
| **Total** | **18+ endpoints** | **24 weeks** |

**This is a 6-month project**, not aligned with our current 1-3 week sprints.

### Problem 2: External Dependencies

| Dependency | Status | Risk |
|------------|--------|------|
| CIBERSORT/TIMER2 | ❌ Not integrated | HIGH - Needs external tool |
| NetMHCpan | ❌ Not integrated | HIGH - Needs local install or API |
| Spatial data (Visium) | ❌ Not integrated | HIGH - Needs new parsers |
| TCGA data | ⚠️ Scripts exist, no service | MEDIUM - Needs unification |

### Problem 3: Validation Targets May Be Optimistic

| Target | Claimed | Realistic |
|--------|---------|-----------|
| T-Cell recruitment r | > 0.40 | Unknown (no baseline) |
| TCF1/Tox signature r | > 0.50 | Unknown (no baseline) |
| TLS readiness r | > 0.35 | Unknown (no baseline) |

**Same issue as other agents**: We don't know if mechanism alignment → outcome prediction.

### Problem 4: No MVP Definition

The plan is all-or-nothing. No clear MVP that delivers value in 2-4 weeks.

---

## Corrections to the Plan

### Correction 1: Mechanism Vector Service Status

**Plan Says**: "File `api/services/pathway_to_mechanism_vector.py` exists but is EMPTY - needs implementation"

**Actual State**: ✅ **FULLY IMPLEMENTED** (317 lines)
- `convert_pathway_scores_to_mechanism_vector()` exists (line 185)
- `normalize_pathway_name()` exists (line 61)
- `validate_mechanism_vector()` exists (line 91)
- TP53→DDR 50% contribution implemented (line 232)

**Impact**: Phase 0 Task 0.1 is ALREADY DONE.

### Correction 2: Pathway Score Extraction

**Plan Says**: Need to create `pathway_score_extractor.py`

**Actual State**: `get_mechanism_vector_from_response()` already exists in `pathway_to_mechanism_vector.py` (line 283)

**Impact**: Phase 0 Task 0.2 is PARTIALLY DONE.

---

## Manager Questions Answered

### Q1: TP53 Contribution to DDR

**Question**: Should TP53 be added to DDR (50% weight)?

**Answer**: ✅ **ALREADY IMPLEMENTED** in `pathway_to_mechanism_vector.py:232`

```python
mechanism_vector[ddr_idx] = ddr_score + (tp53_score * 0.5)
```

### Q2: TCGA Service Architecture

**Question**: Create reusable service or keep scripts standalone?

**Answer**: **DEFER** - Use existing scripts for now. Create service when we have 3+ use cases.

### Q3: Deconvolution Service Choice

**Question**: Which tool: CIBERSORTx, TIMER2, quanTIseq?

**Answer**: **DEFER to Phase 1** - Start with bulk RNA proxy only. Add deconvolution when validated.

### Q4: NetMHCpan Integration

**Question**: Local installation or API wrapper?

**Answer**: **DEFER to Phase 3** - Not needed for proxy classifiers. Add when TOX endpoints built.

### Q5: Spatial Data Processing Priority

**Question**: Start with bulk RNA proxy only?

**Answer**: **YES** - Start with bulk RNA (LOW confidence). Add spatial in Phase 6.

### Q6: Pathway Score Extractor

**Question**: Handle all 3 sources or prioritize SAE only?

**Answer**: **ALREADY IMPLEMENTED** - `get_mechanism_vector_from_response()` handles this.

### Q7: Evo2 Integration for MM Doctrines

**Question**: Should MM endpoints use Evo2?

**Answer**: **YES, for design endpoints** - Follow existing `/api/design/*` pattern.

---

## Recommended Approach: Phase-Gated MVP

### Option A: Full 24-Week Plan ❌ NOT RECOMMENDED

Too much scope, too many dependencies, no near-term value.

### Option B: Phase-Gated MVP ✅ RECOMMENDED

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

**Phase 1B: Second Proxy Classifier (2 weeks)** - If Phase 1A succeeds

**Phase 2: Design Tools (4 weeks)** - Only if Phase 1 validates

---

## Strategic Alignment Required

### Same Rule as Other Agents

**Mechanism Alignment vs Outcome Prediction**:

| Claim Type | Allowed | Not Allowed |
|------------|---------|-------------|
| Mechanism alignment | ✅ "This predicts chemokine signature" | ❌ "This predicts checkpoint response" |
| Literature validation | ✅ "Published studies show correlation" | ❌ "We validated r > 0.40" |
| Proxy confidence | ✅ "MED confidence (r = 0.3-0.5)" | ❌ "HIGH confidence (r > 0.50)" |

### Add Honest Framing

The plan should state:

> "These are proxy classifiers based on bulk RNA signatures. They predict mechanism alignment, not clinical outcomes. Validation against published checkpoint trial cohorts is planned but not yet complete."

---

## Recommended Priority Order

Given the current agent landscape:

### Immediate (This Week)
1. **Zo2**: Execute Phase 1 (biomarker integration)
2. **Trial Matching**: Start Day 1 (mechanism fit wiring)
3. **Food Validation**: Make alignment changes
4. **MM Doctrines**: ❌ DO NOT START - Wait for above to complete

### Next Sprint (Weeks 2-3)
5. **Zo2**: Phase 2 (if Phase 1 checkpoint passes)
6. **Trial Matching**: Days 2-3 (resistance, validation)
7. **MM Doctrines**: Start Phase 0 test suite only

### Following Sprint (Weeks 4-5)
8. **MM Doctrines**: Start ONE proxy classifier (Phase 1A)
9. Checkpoint after Phase 1A before continuing

---

## Action Items for MM Agent

### Before Starting Implementation

1. **Update Plan**: Remove Phase 0 Task 0.1 (already done)
2. **Update Plan**: Acknowledge pathway extraction exists
3. **Add MVP Definition**: Define 2-week deliverable
4. **Add Checkpoints**: Add Phase 1A checkpoint before Phase 1B
5. **Add Honest Framing**: Clarify mechanism alignment vs outcome prediction

### When Starting Implementation (After Other Agents)

6. **Pick ONE Proxy Classifier**: Start with TLS or TCF1 (less external dependencies)
7. **Skip Deconvolution Initially**: Use bulk RNA proxy only
8. **Create Test Suite First**: Known biology tests before implementation
9. **Validate Before Expanding**: Checkpoint after first classifier

---

## Summary

| Aspect | Assessment |
|--------|------------|
| Research quality | ✅ **Excellent** (thorough, code-verified) |
| Scope | ⚠️ **Too large** (24 weeks, 18 endpoints) |
| Dependencies | ⚠️ **Many external tools needed** |
| Priority | ⚠️ **Should wait for other agents** |
| MVP definition | ❌ **Missing** (needs 2-week deliverable) |
| Corrections | ⚠️ **Mechanism vector exists, not empty** |

---

## Final Verdict

**Status**: ⚠️ **HOLD - NEEDS SCOPE REDUCTION**

**Recommended Action**:
1. Other agents complete their current work first
2. MM agent creates test suite and MVP definition
3. Start with ONE proxy classifier after other agents done
4. Checkpoint before expanding scope

---

## Documents Created

- `.cursor/ayesha/STRATEGIC_REVIEW_MM_DOCTRINES.md` (this document)

---

**The righteous path: Great research deserves proper sequencing. Don't start a 6-month project when 1-week projects are unfinished.**

