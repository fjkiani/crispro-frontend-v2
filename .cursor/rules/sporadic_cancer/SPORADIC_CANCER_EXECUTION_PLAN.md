# ⚔️ SPORADIC CANCER STRATEGY - EXECUTION PLAN

**Date**: January 5, 2025  
**Commander**: Alpha  
**Executor**: Zo (Agent)  
**Mission**: Build tumor-centric analysis for sporadic (germline-negative) cancer patients

**Related Doctrine**: `.cursor/rules/specialized_systems/archive/sporadic_cancer_strategy_doctrine.mdc`

---

## ✍️ CHANGELOG (Agent Quick Read)

- Model/schema
  - Switched `TumorContext` to Pydantic BaseModel (not dataclass) with enums/validation.
  - Added validation tasks: `msi_status` ∈ {"MSI-H","MSS",null}; clamp numeric fields (`tmb >= 0`, `0 ≤ hrd_score ≤ 100`); compute `completeness_score`.
- Level 0 behavior
  - Do not infer MSI at Level 0: set `"msi_status": null` (no MSI-derived boosts).
  - Added `priors_refresh_date` to provenance when priors are used.
- API contracts
  - Efficacy response: added `confidence_version` and `priors_refresh_date` under `provenance.flags` and top-level `provenance.confidence_version`.
  - Ingest endpoint now accepts optional `"report_json": {}` (preferred when available) to bypass PDF parsing.
  - Added contract-stability note: preserve existing shapes; add only under `provenance`/`provenance.flags`.
- Clinical rules
  - Explicitly stated: “Unknown MSI → do not infer; treat as null (no MSI boosts).”
  - Reaffirmed clamps/caps for boosts/penalties.
- Day 1 tasks & acceptance
  - Added tasks to define Pydantic enums/validation and completeness scoring.
  - Acceptance now requires `confidence_version` and `priors_refresh_date` in outputs.
- Notes and guardrails
  - M1 now includes explicit contract stability warning (no breaking changes).
  - Added preference for JSON ingestion and no PDF storage reiteration.

---

## 🛡️ AGENT EXECUTION SAFETY CHECKLIST (READ FIRST!)

**⚠️ BEFORE YOU START CODING, AGENT:**

### **✅ DO:**
1. **READ EXISTING CODE FIRST** - Always use `codebase_search` or `read_file` before creating new files
2. **EXTEND, DON'T REPLACE** - Add new methods/classes, don't rewrite existing logic
3. **USE EXISTING PATTERNS** - Reference similar files (e.g., `BiomarkerContext` for schemas)
4. **CREATE NEW FILES ONLY** - Don't modify core orchestrators without explicit approval
5. **SIMPLE IMPLEMENTATIONS** - Start with basic logic, iterate based on feedback
6. **TEST AS YOU GO** - Write tests for each module before moving to next
7. **ASK WHEN UNCERTAIN** - If a requirement is ambiguous, ask for clarification

### **❌ DON'T:**
1. **DON'T HALLUCINATE DATA** - Use published TCGA stats or conservative estimates only
2. **DON'T REWRITE CORE LOGIC** - `ayesha_orchestrator.py`, `efficacy_orchestrator` are off-limits
3. **DON'T ADD COMPLEX DEPENDENCIES** - No PDF parsing libraries, no ML models
4. **DON'T CREATE NEW DATABASES** - Use existing logging/provenance infrastructure
5. **DON'T OVERCOMPLICATE** - Simple multipliers (0.5x, 1.3x) not ML models
6. **DON'T MODIFY EXISTING COMPONENTS** - Create new ones in `sporadic/` subdirectory
7. **DON'T SKIP TESTING** - Every module needs at least 1 smoke test

### **🎯 SUCCESS CRITERIA:**
- ✅ All new code in NEW files (not modifications to existing files)
- ✅ All integrations are EXTENSIONS (new methods, not replacements)
- ✅ All data is SOURCED (published stats or conservative estimates)
- ✅ All tests PASS (pytest with golden snapshots)
- ✅ All endpoints DOCUMENTED (docstrings + provenance)

**IF YOU VIOLATE THESE RULES, EXECUTION WILL BE HALTED AND ROLLED BACK!**

---

## ⚠️ CRITICAL QUESTIONS FOR MANAGER (ANSWER BEFORE EXECUTION)

### **Q1: Tumor NGS Parser Strategy** 🤔
**Question**: Foundation Medicine and Tempus reports have different formats. Should I:
- **Option A**: Build parsers for both (2 days work) ✅ **COMPREHENSIVE**
- **Option B**: Start with Foundation Medicine only (1 day) ⚠️ **FASTER BUT INCOMPLETE**
- **Option C**: Create generic JSON schema and require manual conversion ⚠️ **PUNTS PROBLEM**

**My Recommendation**: Option B (Foundation first), then Option A when proven valuable

**Risk**: If Ayesha has Tempus report, we won't support it yet

---

### **Q2: PARP Penalty Logic** 🤔
**Question**: When germline negative, how much should PARP inhibitor score drop?
- **Current**: PARP gets 0.8-0.9 efficacy for BRCA+ patients
- **Proposed**: Drop to 0.3-0.4 for germline-negative?
- **Evidence**: HRD can exist without germline BRCA, but PARP less effective

**My Recommendation**: Drop PARP to 0.4 for germline-negative UNLESS tumor shows HRD signature

**Risk**: Too aggressive penalty might miss sporadic HRD cases

---

### **Q3: Clinical Trials Filtering** 🤔
**Question**: Should germline-negative patients see BRCA-required trials at all?
- **Option A**: Hard exclude (cleaner UX) ✅ **MY PREFERENCE**
- **Option B**: Show but mark "Not Eligible" (transparent but cluttered)
- **Option C**: Show if tumor has HRD signature (conditional)

**My Recommendation**: Option A with explanation text

**Risk**: Might hide trials that accept sporadic HRD

---

### **Q4: Testing Strategy** 🤔
**Question**: How do I test without real tumor NGS reports?
- **Option A**: Mock Foundation Medicine JSON (faster but less realistic)
- **Option B**: Use Ayesha's anonymized report (realistic but privacy concern)
- **Option C**: Download Foundation demo report from their website

**My Recommendation**: Option C if available, else Option A

**Risk**: Mocked tests might miss real-world edge cases

---

### **Q5: Frontend Priority** 🤔
**Question**: Day 5 frontend work - which is P0?
- **Component A**: GermlineStatusBanner (shows negative status) - **AWARENESS**
- **Component B**: TumorNGSUpload (upload .pdf/.json) - **FUNCTIONALITY**
- **Component C**: Updated trial cards (hide BRCA-required) - **UX**

**My Recommendation**: Do all 3 (they're small), but prioritize B → A → C

**Risk**: Upload without banner might confuse users

---

### **Q6: Backward Compatibility** 🤔
**Question**: Should old endpoint behavior (germline-focused) still work?
- **Option A**: Add new fields but keep defaults (backward compatible) ✅ **SAFE**
- **Option B**: Make germline_status required (breaking change) ⚠️ **RISKY**

**My Recommendation**: Option A - make all new fields optional with smart defaults

**Risk**: None if we maintain old behavior as default

---

## 🎯 MANAGER'S DECISIONS (FILL THIS OUT):

```markdown
Q1 (Tumor NGS Parsers): Option B - Foundation first (Tempus in next iteration)
Q2 (PARP Penalty): Apply 0.6x multiplier when germline_negative AND HRD < 42 AND no somatic BRCA biallelic loss.
                      Remove penalty (1.0x) if HRD ≥ 42 OR biallelic BRCA loss OR genomic scar signature present.
Q3 (Trials Filtering): Default Option A (hard exclude germline-required). Auto-include if tumor HRD+;
                      provide “Show excluded trials” toggle in UI for transparency.
Q4 (Testing Strategy): Option C (Foundation demo report) primary; fallback Option A (faithful mock JSON).
                       Option B only with explicit consent and redaction.
Q5 (Frontend Priority): All 3, order: B → A → C (Upload → Banner → Trial cards)
Q6 (Backward Compatibility): Option A - Optional fields with smart defaults; no breaking changes.

Additional Notes:
- Timeline OK to extend to 10 days if Tempus parser is included.
- Prioritize correctness and auditability over breadth; ship Foundation path first.
```

---

## 🤔 ZO'S CRITICAL IMPLEMENTATION QUESTIONS (ANSWER BEFORE DAY 1)

### **Q1: Disease Priors Source & Maintenance** 🔬
**Question**: Where exactly are we getting the disease priors (HRD prevalence, TMB/MSI distributions)?
- **Option A**: Manual curation from literature (static, needs periodic updates)
- **Option B**: Automated extraction from Cohort Lab summaries (dynamic, version-controlled)
- **Option C**: Hybrid (seed from literature, augment from Cohort Lab) ✅ **ZO'S RECOMMENDATION**

**Why it matters**: 
- If priors are wrong, Level 0 outputs are misleading
- Provenance must cite source (`disease_priors_version`)
- How often do we update? Who validates?

**Manager's Decision**: Option C — Hybrid (seed from literature, augment from Cohort Lab). Update monthly; each release tagged via `disease_priors_version` and `priors_refresh_date`.

---

### **Q2: Platinum Response as HRD Proxy - How Reliable?** 🧪
**Concern**: Platinum sensitivity correlates with HRD, but it's not 1:1
- Platinum-sensitive could be: HRD-high OR DNA repair proficient with other sensitivities
- Platinum-resistant could be: HRD-low OR high with reversion mutations

**Question**: What's our false positive/negative rate tolerance for Level 0 HRD estimates?
- **Current plan**: Platinum-sensitive → assume HRD-likely (remove PARP penalty)
- **Risk**: Overestimating HRD → inappropriate PARP recommendations

**Zo's Recommendation**: 
- Keep conservative: Platinum-sensitive → "possible HRD" (confidence 0.3)
- Add rationale: "Platinum response suggests potential HRD; tumor NGS recommended for confirmation"
- **DO NOT** fully remove PARP penalty unless HRD score ≥42

**Manager's Decision**: Adopt conservative policy. Label as “possible HRD” (confidence ~0.3). Do not remove PARP penalty unless HRD ≥ 42 or BRCA biallelic loss or genomic scar present.

---

### **Q3: Co-Pilot Interactive Elevation - Implementation Strategy** 🤖
**Current plan**: Co-Pilot asks targeted questions to elevate Level 0 → 1

**Question**: How does this work in the backend?
- **Option A**: Co-Pilot triggers `/api/tumor/quick_intake` with incremental fields → re-computes each time
- **Option B**: Co-Pilot updates `SessionContext.tumorContext` directly → one final submit ✅ **ZO'S RECOMMENDATION**
- **Option C**: Co-Pilot builds a "pending context" → user reviews → submit

**Manager's Decision**: Option B — Update `SessionContext.tumorContext` and submit once; show diffs in provenance.

---

### **Q4: Confidence Scaling Formula** 📊
**Current plan**: "Confidence scaled based on field completeness"

**Question**: What's the exact formula?

**Zo's Proposed Formula** (Additive with caps):
```python
base_confidence = 0.3  # Level 0 baseline

# Additive boosts
if tmb_present: base_confidence += 0.1
if msi_present: base_confidence += 0.1
if hrd_present: base_confidence += 0.1
if somatic_mutations_present: base_confidence += 0.1

# Caps
# Level 0: 0.3-0.4 (priors only)
# Level 1: 0.4-0.6 (partial data)
# Level 2: 0.6-0.9 (full report)
```

**Manager's Decision**: Approve proposed formula and caps. Expose `confidence_cap` and `completeness_score` in `provenance.flags`.

---

### **Q5: M6 Priors Resource - File Structure** 📁
**Question**: What does `api/resources/disease_priors.json` look like?

**Zo's Proposed Structure**:
```json
{
  "version": "v1.0",
  "last_updated": "2025-01-05",
  "diseases": {
    "ovarian_hgs": {
      "prevalence": {
        "tp53_mutation": 0.95,
        "hrd_high": 0.48,
        "msi_high": 0.02
      },
      "distributions": {
        "tmb": {"median": 5.2, "q1": 3.1, "q3": 8.7, "high_cutoff": 10},
        "hrd": {"median": 32, "q1": 18, "q3": 48, "high_cutoff": 42}
      },
      "platinum_response": {
        "sensitive_hrd_correlation": 0.68
      },
      "sources": ["TCGA-OV", "cBioPortal ovarian_tcga_pan_can_atlas_2018"]
    }
  }
}
```

**Manager's Decision**: Approve structure. Require `version`, `last_updated`, per‑disease `sources`. Add quarterly archive snapshots.

---

### **Q6: PARP Penalty Math - Edge Cases** 🧮
**Current rule**: Apply 0.6x multiplier when germline_negative AND HRD < 42

**Question**: What if we're in Level 0 (no HRD score at all)?
- **Option A**: Apply penalty (conservative)
- **Option B**: No penalty (optimistic)
- **Option C**: Apply smaller penalty (0.8x) with rationale "HRD unknown" ✅ **ZO'S RECOMMENDATION**

**Zo's Proposed Logic**:
```python
if germline_negative:
    if hrd_score is None:  # Level 0
        parp_score *= 0.80
        rationale += "⚠️ PARP efficacy reduced (germline negative, HRD unknown). Tumor NGS recommended."
    elif hrd_score < 42:  # Level 1/2
        parp_score *= 0.60
        rationale += "⚠️ PARP efficacy reduced (germline negative, HRD < 42)"
    else:  # HRD ≥42
        # No penalty
        rationale += "✅ PARP eligible (somatic HRD-high)"
```

**Manager's Decision**: Option C — Apply 0.80x at Level 0 with explicit rationale; adopt full logic block above.

---

### **Q7: TMB/MSI Boosts - Interaction Logic** 💉
**Current rule**: Apply boost if tmb >= 10 OR msi_status == "MSI-high"

**Question**: What if BOTH are true (TMB-high AND MSI-high)?
- **Option A**: Apply boost twice (1.25 * 1.25 = 1.56x) → too aggressive?
- **Option B**: Cap at single boost (1.25x max)
- **Option C**: Stronger boost for both (1.35x) with different rationale ✅ **ZO'S RECOMMENDATION**

**Zo's Proposed Logic**:
```python
if msi_high:
    io_score *= 1.30
    rationale += "✅ Immunotherapy boost (MSI-high)"
elif tmb >= 20:
    io_score *= 1.35
    rationale += "✅ Immunotherapy boost (TMB very high ≥20)"
elif tmb >= 10:
    io_score *= 1.25
    rationale += "✅ Immunotherapy boost (TMB-high ≥10)"

# Cap at 1.0
io_score = min(io_score, 1.0)
```

**Manager's Decision**: Option C — Use stronger boost for both (1.35x for ≥20; 1.30x for MSI-H) with cap at 1.0; adopt logic as proposed.

---

### **Q8: Quick Intake Frontend UX** 🎨
**Question**: Where does `TumorQuickIntake.jsx` appear in the UI?
- **Option A**: Modal on Co-Pilot page
- **Option B**: Dedicated "/quick-intake" route
- **Option C**: Inline in ResearchPortal above WIWFM
- **Option D**: Smart banner "No tumor NGS? Quick intake here" + inline ✅ **ZO'S RECOMMENDATION**

**Manager's Decision**: Option D — Smart banner + inline in ResearchPortal; consider dedicated route in v2.

---

### **Q9: Provenance Flags - Complete List** 📝
**Question**: What's the COMPLETE list of provenance flags for sporadic logic?

**Zo's Proposed Complete List**:
```python
provenance.flags = {
    # Level/Mode
    "level": "L0" | "L1" | "L2",
    "no_report_mode": bool,
    "tumor_context_source": "Foundation" | "Tempus" | "Manual" | "Quick Intake" | None,
    
    # Priors
    "disease_priors_used": bool,
    "disease_priors_version": "v1.0",
    "platinum_proxy_used": bool,
    
    # Gating Logic
    "germline_negative": bool,
    "parp_penalty_applied": bool,
    "parp_penalty_factor": 0.60 | 0.80 | 1.0,
    "immunotherapy_boost_applied": bool,
    "immunotherapy_boost_factor": 1.25 | 1.30 | 1.35,
    
    # Biomarkers
    "hrd_score_available": bool,
    "tmb_available": bool,
    "msi_available": bool,
    
    # Confidence
    "confidence_cap": 0.4 | 0.6 | None,
    "completeness_score": 0.0-1.0
}
```

**Manager's Decision**: Approve list. Add `confidence_version` and `priors_refresh_date` fields; enforce presence in responses.

---

## 🎯 MISSION CONTEXT

### **THE PROBLEM:**
- **85-90% of cancers are sporadic** (not hereditary)
- Ayesha's germline testing: **NEGATIVE** (38 genes, CustomNext-Cancer®)
- Traditional platforms focus on hereditary cancers (only 10-15% of patients)
- **Our platform needs to shift from germline → tumor-centric analysis**

### **THE OPPORTUNITY:**
- Address the **majority** of cancer patients (sporadic cases)
- Integrate tumor NGS + treatment history + clinical trials
- Competitive advantage: Most platforms ignore sporadic cases

### **AYESHA'S CASE (REFERENCE):**
- **Germline Status**: Negative (not hereditary)
- **Cancer Type**: High-grade serous ovarian carcinoma
- **Stage**: IIIC-IV (peritoneal carcinomatosis, ascites)
- **Treatment Line**: 3 (post-platinum progression)
- **Need**: Tumor NGS to find somatic drivers (TP53, possible HRD)

---

## 📋 WHAT WE WILL BUILD

## 🧱 MODULAR ARCHITECTURE & HIERARCHY (BUILD-BY-MODULE)

We will ship a composable system with clear module boundaries, contracts, and flags. Each module can be developed, tested, and deployed independently.

### Hierarchy of Operation Modes
- Level 0: No-report mode (Quick Intake → disease priors + proxies) → low-confidence outputs
- Level 1: Partial data mode (hand-entered TMB/MSI/HRD/mutations) → medium confidence
- Level 2: Full report mode (Foundation/Tempus parsers) → highest confidence

### Module Map

M1. Intake & Orchestration
- Responsibility: Collect inputs (Quick Intake, Upload, Manual), drive Level 0/1/2 flow
- Endpoints:
  - POST `/api/tumor/quick_intake` (Level 0)
  - POST `/api/tumor/ingest_ngs` (Level 2)
  - POST `/api/efficacy/predict` (dispatch with `germline_status`, `tumor_context`)
- Inputs: patient context, germline status, intake fields or file
- Outputs: `TumorContext`, normalized payload to scoring
- Flags: `ENABLE_QUICK_INTAKE`, `ENABLE_REPORT_UPLOAD`
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **NO EXISTING TUMOR CODE** - confirmed via codebase search, safe to build fresh
  - ✅ **PatientContext exists** (`api/schemas/ayesha.py`) - extend it, don't replace
  - ✅ **Create NEW file:** `api/schemas/tumor_context.py` (not in existing files)
  - ⚠️ **DO NOT MODIFY** `ayesha_orchestrator.py` scoring logic directly - extension only
  - ⚠️ **CONTRACT**: Preserve existing API response shapes; add only under `provenance`/`provenance.flags`

M2. Normalization Layer (TumorContext)
- Responsibility: Normalize parsed/entered data into a single schema with QC + provenance
- Interfaces: `TumorContext` Pydantic BaseModel + validators
- Inputs: parsed metrics, manual fields
- Outputs: structured `TumorContext` with completeness score
- Flags: `STRICT_SCHEMA_VALIDATION`, `MISSING_FIELD_WARN`
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **Use Pydantic BaseModel**, not dataclass (matches existing schema patterns)
  - ✅ **Reference existing pattern:** See `BiomarkerContext` in `api/schemas/ayesha.py`
  - ⚠️ **Schema fields:** Use Optional[] for all fields with Field(default=None)
  - ⚠️ **Completeness score:** Simple count(non-None) / total_fields, don't overcomplicate

M3. Scoring Engine (Sporadic Logic)
- Responsibility: Apply germline gating, PARP penalty, IO boosts, confidence scaling
- Inputs: `germline_status`, `TumorContext` (TMB/MSI/HRD/BRCA)
- Outputs: adjusted `efficacy_score`, `confidence`, `rationale`, `provenance.flags`
- Flags: `ENABLE_PARPTUNE`, `ENABLE_IO_TMB`, `CONFIDENCE_CAP_NO_REPORT`
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **Integration point:** Extend `EfficacyOrchestrator.predict_drug_efficacy()`
  - ✅ **Existing file:** `api/services/efficacy_orchestrator/orchestrator.py`
  - ⚠️ **DO NOT REWRITE** scoring logic - add new methods for sporadic adjustments only
  - ⚠️ **PARP penalty logic:** `if germline_status == "negative" and drug_class == "PARP": score *= 0.5`
  - ⚠️ **TMB boost logic:** `if tumor_context.tmb_high and drug_class in ["checkpoint_inhibitor"]: score *= 1.3`
  - ⚠️ **Confidence cap:** `if tumor_context is None: confidence *= 0.7`

M4. Trials Module (Sporadic-Aware Search)
- Responsibility: Filter/rank trials using sporadic logic and tumor biomarkers
- Endpoints: POST `/api/clinical_trials/search`
- Inputs: `germline_status`, `TumorContext`, line/stage
- Outputs: ranked trials with biomarker badges and exclusion counts
- Flags: `HARD_EXCLUDE_GERMLINE_REQUIRED`, `SHOW_EXCLUDED_TOGGLE`
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **Existing service:** `AutonomousTrialAgent` in `api/services/autonomous_trial_agent.py`
  - ✅ **Integration:** Add filters to existing search, don't rewrite from scratch
  - ⚠️ **Germline exclusion:** Check trial description/criteria for "BRCA", "germline", "hereditary"
  - ⚠️ **TMB/MSI boosting:** Simple keyword match on "biomarker", "TMB", "MSI" in trial criteria

M5. Frontend UX Module
- Responsibility: UI for Banner, Quick Intake, Upload, Trials, Confidence/Provenance
- Components:
  - `GermlineStatusBanner.jsx`
  - `TumorQuickIntake.jsx`
  - `TumorNGSUpload.jsx`
  - Trials result badges + sporadic-aware labels
- State: `SessionContext` holds `germline_status`, `tumorContext`, `run_id`
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **Location:** `oncology-frontend/src/components/` (create new subdirectory `sporadic/`)
  - ✅ **Existing pattern:** Reference `CoPilot/` directory structure for modular components
  - ⚠️ **DO NOT MODIFY** existing components - create new ones only
  - ⚠️ **Banner placement:** Add to `ResearchPortal` page top, not global App.jsx
  - ⚠️ **File upload:** Use `<input type="file" accept=".pdf,.json">` + FormData, NO PDF parsing library

M6. Priors & Knowledge (Resources)
- Responsibility: Provide disease priors (HRD prevalence, MSI/TMB distributions, driver frequencies)
- Files: `api/resources/disease_priors.json` (versioned), loader with TTL cache
- Inputs: `cancer_type/subtype`, platinum response proxy
- Outputs: prior-informed estimates for Level 0/1
- Flags: `PRIORS_VERSION`, `ENABLE_PRIORS_WARNINGS`
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **Source:** Use TCGA published stats from NCI (https://portal.gdc.cancer.gov/)
  - ✅ **Format:** Simple JSON with disease_type → { hrd_prevalence, tmb_median, msi_prevalence }
  - ⚠️ **DO NOT hallucinate numbers** - use conservative estimates if exact data unavailable
  - ⚠️ **Ovarian HRD:** ~50% (published), Breast HRD: ~20%, Colorectal MSI: ~15%
  - ⚠️ **TMB medians:** Lung: 8-10 mut/Mb, Melanoma: 15-20 mut/Mb, Ovarian: 3-5 mut/Mb

M7. Privacy & Security (Policy)
- Responsibility: PII/retention controls; hashed artifacts; audit fields
- Policies: no raw PDF retention by default; encrypt if retained; 7-day TTL; redaction
- Logs: `run_id`, `report_hash`, `parser_version`, timestamps only
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **NO PDF STORAGE** - read, parse, discard immediately
  - ✅ **Hash only:** Store SHA256 hash of report for audit, not raw file
  - ⚠️ **DO NOT implement encryption** - complexity not needed for Phase 1
  - ⚠️ **Logs:** Use existing logging infrastructure, don't create new DB tables

M8. QA, Testing & Observability
- Responsibility: Unit tests per module; golden snapshots; smoke tests; metrics/logging
- Tests:
  - M1 Quick Intake → L0 TumorContext
  - M2 schema validation + completeness scoring
  - M3 gates math (PARP, IO boost) under different contexts
  - M4 filters/ranking with biomarker badges
- Telemetry: latency, error rates, confidence distribution, priors usage counters
- **⚠️ AGENT EXECUTION NOTES:**
  - ✅ **Test location:** `oncology-backend-minimal/tests/test_sporadic_*.py`
  - ✅ **Golden snapshots:** Use `pytest` with `--snapshot-update` for baseline
  - ⚠️ **DO NOT test PDF parsing** - mock the parser output, test normalization only
  - ⚠️ **Metrics:** Use simple counters in provenance, don't add Prometheus/Grafana yet

---

## 🎯 VALUE & OUTCOMES (WHY THIS MATTERS)

- **Serve the 85–90% majority**: Deliver high‑quality, tumor‑centric analysis for germline‑negative patients, not just hereditary cases.
- **Operate without friction**: Work even when no NGS report is available (Level 0), then progressively improve as data arrives.
- **Faster decisions**: Immediate, transparent, research‑grade recommendations with confidence and provenance so providers can act.
- **One experience, not silos**: Connect sporadic analysis to our WIWFM, Evidence, Cohort Lab, and Dossier flows—no dead ends.
- **Audit‑ready**: Every recommendation comes with run_id, inputs, policy flags, and thresholds so it’s defensible.

Outcomes for Ayesha (example):
- Today: quick intake → conservative guidance and trial options with low confidence (explicitly labeled).
- Tomorrow: upload Foundation report → re‑compute with full HRD/TMB/MSI → higher confidence, refined trials and therapy fit.
- Always: traceable provenance, printable provider report, and next‑step prompts.

---

## 🔗 INTEGRATION WITH EXISTING PLATFORM CAPABILITIES

This plan reuses and elevates what we already built:

- **Efficacy Orchestrator (S/P/E)**
  - Reused by M3 Scoring Engine. We add sporadic gates (PARP penalty, IO boosts), but keep the same contracts.
  - Insights bundle continues to flow into WIWFM; rationale/provenance extended with sporadic flags.

- **Insights Endpoints (Functionality/Chromatin/Essentiality/Regulatory)**
  - Optional enrichment for variants present in `TumorContext`; chips appear in Dossier and VUS flows.
  - Confidence modulation leverages insights lifts as before.

- **Fusion Engine (AlphaMissense)**
  - When missense GRCh38 coverage exists, maintain coverage checks and fused‑S contribution.
  - Display Fusion eligibility in Dossier/VUS; non‑blocking in sporadic pipeline.

- **Evidence Services (ClinVar/Literature)**
  - ClinVar priors surfaced when variants overlap; literature remains best‑effort with clear fallbacks.
  - Trial ranking reasons include biomarker evidence strings (e.g., “TMB ≥10”).

- **Cohort Lab (Datasets/Benchmarks)**
  - Level 0 priors seeded from Cohort Lab summaries (disease‑level TMB/MSI/HRD prevalence).
  - Future: Show cohort overlays (e.g., HRD outcomes) in the Dossier for sporadic context.

- **Sessions & Caching**
  - Quick Intake/Upload outputs persisted in `SessionContext` and Redis; single‑flight prevents duplicate parsing.
  - Provenance shows cache hits for reproducibility and speed.

- **Co‑Pilot**
  - Drives elevation from Level 0 → 1 by asking targeted questions (TMB, MSI, HRD, BRCA).
  - Can trigger WIWFM, Trials, and Dossier generation with the current `TumorContext`.

- **IND Package / Dossier**
  - Provider report enhanced with sporadic fields (HRD/TMB/MSI, penalties/boosts applied, trial badges).
  - Same document generation path; new sections added for tumor context and confidence rationale.

- **Design & Safety (Demo)**
  - CRISPR Readiness panel continues to consume chromatin/essentiality where present; remains clearly RUO/demo.

Integration Map (Module → Existing Capability):
- M1/M2 → Efficacy Orchestrator (input contracts extended), Sessions
- M3 → WIWFM (sporadic gates), Insights lifts, Fusion eligibility
- M4 → Trials router/service (filters + biomarker badges), Evidence
- M5 → Frontend components (Banner/Upload/Trials/Dossier/VUS), SessionContext
- M6 → Datasets (Cohort Lab priors), versioned resources
- M7 → Existing privacy/policy framework (hashing, TTLs, provenance)
- M8 → Test harness (golden snapshots), client logger, health metrics

---

### **Core Capabilities:**

1. **Germline-Status Gating** ✅ Backend
   - Down-weight hereditary-targeted drugs (e.g., PARP) when germline negative
   - Unless somatic HRD is high (tumor-based)
   - Boost tumor-agnostic therapies (immunotherapy if TMB-high)

2. **TumorContext Support** ✅ Backend
   - Accept somatic mutations (TP53, KRAS, PIK3CA, etc.)
   - TMB (tumor mutational burden)
   - MSI (microsatellite instability)
   - HRD score (somatic homologous recombination deficiency)
   - Copy number alterations (amplifications/deletions)

3. **Tumor NGS Parsers** ✅ Backend
   - Foundation Medicine report parser
   - Tempus report parser
   - Output: Structured `TumorContext` JSON

4. **Sporadic-Aware Clinical Trials** ✅ Backend
   - Exclude "BRCA-required" trials
   - Prioritize tumor-agnostic trials (TMB-high, MSI-high)
   - Prioritize somatic biomarker trials

5. **Frontend UX** ✅ Frontend
   - Germline Status Banner
   - Tumor NGS upload/parse
   - Sporadic-aware trial cards

6. **Provenance & Caching** ✅ Backend
   - Persist `germline_status`, `tumor_context_source`
   - Cache parsed NGS reports
   - Track trial queries with TTL

---

## 🧭 REPORT-AGNOSTIC & PARTIAL-DATA CAPABILITY (NO NGS REPORT MODE)

We will support three intake levels with progressive enhancement and transparent confidence:

- Level 0 (No Report – Minimal Intake):
  - Inputs: cancer type/subtype, stage, treatment line, prior platinum response (sens/res), basic labs if available.
  - Heuristics:
    - Disease priors (e.g., HGSOC: ~95% TP53 mutation; HRD prevalence ~40-50%).
    - Platinum sensitivity as HRD proxy (sens → HRD-likely).
    - MSI/TMB disease prevalence distributions (from cBioPortal summaries).
  - Output: Conservative recommendations with low confidence; explicit RUO and “data needed” prompts.
  - Confidence: capped at 0.3-0.4; explicit provenance “no_report_mode”.

- Level 1 (Partial Structured Data):
  - Inputs: manually entered mutations (gene + hgvs_p), TMB number (if known), MSI status (IHC/PCR), HRD score (if provided by lab).
  - Output: Full sporadic logic (penalties/boosts) applied on available fields; unknowns set to “insufficient data”.
  - Confidence: scaled based on field completeness (e.g., +0.1 for each of TMB/MSI/HRD present).

- Level 2 (Full Report Ingestion):
  - As already specified (parsers for Foundation/Tempus).
  - Confidence: standard computation with complete provenance.

Interactive Completion (Co‑Pilot):
- The Co‑Pilot will ask targeted questions to elevate from Level 0 → Level 1:
  - “Do you have a TMB value on the pathology report?”
  - “Is MSI-H reported (IHC/PCR)?”
  - “Any known tumor variants (e.g., TP53, BRCA1/2)?”
  - “Response to last platinum therapy?”

Confidence & RUO:
- Always display a confidence bar linked to completeness.
- Level 0 outputs labeled “preliminary – report recommended”.
- Provenance flags: `no_report_mode`, `disease_priors_used`, `platinum_proxy_used`.

Data Sources:
- Disease-level distributions and priors from cBioPortal/GDC snapshots (bundled JSON in `api/resources/disease_priors.json`).
- Maintain versioned provenance for priors.

---

## 🏗️ STEP-BY-STEP IMPLEMENTATION (7 DAYS)

### **✅ DAY 1-2: Backend Foundation**

**Tasks:**
- [ ] Create `api/schemas/tumor_context.py`
- [ ] Extend `predict_drug_efficacy` with `germline_status` and `tumor_context`
- [ ] Implement PARP penalty logic:
  - If `germline_status == "negative"` AND (`hrd_score` missing OR `hrd_score < 42`)
  - Then: Penalize PARP class efficacy_score by 0.3-0.5
  - Rationale: "PARP inhibitors less effective without germline BRCA or high somatic HRD"
- [ ] Implement tumor-agnostic boosts:
  - If `tmb >= 10` OR `msi_status == "MSI-high"`
  - Then: Boost checkpoint inhibitor class by 0.2-0.3
  - Rationale: "Immunotherapy effective in TMB-high/MSI-high tumors"
- [ ] Surface decisions in `rationale[]` and `provenance.flags`
- [ ] Add reportless heuristics in orchestrator:
  - If Level 0: apply disease priors and platinum proxy with conservative weights; cap confidence.
  - Record `no_report_mode` and `disease_priors_version` in provenance.
- [ ] Define `TumorContext` as Pydantic BaseModel with enums/validation:
  - `msi_status`: `"MSI-H" | "MSS" | null` (unknown stays null; do not infer)
  - Clamp numeric fields: `tmb >= 0`, `0 <= hrd_score <= 100` (reject outliers)
  - Provide `completeness_score` = count(non-null tracked fields)/total_tracked_fields

**Files to Create/Modify:**
- `oncology-coPilot/oncology-backend-minimal/api/schemas/tumor_context.py` (NEW)
- `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py` (MODIFY)
- `oncology-coPilot/oncology-backend-minimal/api/routers/efficacy.py` (MODIFY - add params)

**Acceptance:**
- PARP penalty works: `germline_status="negative"` + no HRD → PARP score drops
- Immunotherapy boost works: `tmb=18` → checkpoint inhibitor score increases
- Confidence/provenance include `confidence_version` and `priors_refresh_date`

---

### **✅ DAY 3: Tumor NGS Parsers**

**Tasks:**
- [ ] Create router: `POST /api/tumor/ingest_ngs` (accepts PDF/JSON, returns `TumorContext`)
- [ ] Implement Foundation Medicine parser:
  - Extract TMB, MSI status, HRD score
  - Extract variant table (gene, hgvs_p, VAF)
  - Extract copy number alterations (amplifications/deletions)
- [ ] Implement Tempus parser:
  - Parse comparable fields
  - Normalize to shared `TumorContext` schema
- [ ] Store minimal cache keyed by report hash
- [ ] Include `tumor_context_source` in provenance

**Files to Create:**
- `oncology-coPilot/oncology-backend-minimal/api/routers/tumor.py` (NEW)
- `oncology-coPilot/oncology-backend-minimal/api/services/tumor_ngs_parser.py` (NEW)
- `oncology-coPilot/oncology-backend-minimal/api/resources/disease_priors.json` (NEW)

**Acceptance:**
- Foundation Medicine report → valid `TumorContext` JSON
- Tempus report → valid `TumorContext` JSON
- Cached by file hash (subsequent uploads instant)

---

### **✅ DAY 4: Clinical Trials Filtering**

**Tasks:**
- [ ] Extend `POST /api/clinical_trials/search` to accept:
  - `germline_status` (negative/positive/unknown)
  - `tumor_context` hints (TMB, MSI, HRD, mutations)
- [ ] Implement filters:
  - Exclude trials requiring germline BRCA/Lynch
  - Include tumor-agnostic trials (TMB/MSI biomarker-based)
  - Include somatic HRD trials
- [ ] Rank by match strength:
  - Biomarker match (TMB/MSI/HRD) > Mechanism match > Line eligibility

**Files to Modify:**
- `oncology-coPilot/oncology-backend-minimal/api/routers/clinical_trials.py` (MODIFY)
- `oncology-coPilot/oncology-backend-minimal/api/services/hybrid_trial_search_service.py` (MODIFY)

**Acceptance:**
- BRCA-required trials excluded when `germline_status="negative"`
- TMB-high trials prioritized when `tmb >= 10`
- Trial cards show biomarker match badges

---

### **✅ DAY 5: Frontend Wiring**

**Tasks:**
- [ ] Create `GermlineStatusBanner.jsx` component:
  - Shows when `germline_status === 'negative'`
  - Message: "Sporadic Cancer: Germline testing negative. Analysis focused on tumor genomics and treatment history."
- [ ] Create NGS upload component:
  - File upload (PDF/JSON)
  - Calls `POST /api/tumor/ingest_ngs`
  - Stores `tumorContext` in `SessionContext`
  - Displays parsed TMB, MSI, HRD, mutations
- [ ] Update clinical trial results:
  - Show "Sporadic-aware" label
  - Hide hereditary-only trials
  - Highlight biomarker-matched trials (TMB/MSI/HRD badges)
- [ ] Add “Quick Intake” form (no report):
  - Fields: cancer type/subtype, stage, line, prior platinum response, optional TMB/MSI/HRD numbers.
  - Calls `/api/tumor/quick_intake` to generate Level 0/1 `TumorContext`.
  - Render confidence band and “what to add next” hints.

**Files to Create/Modify:**
- `oncology-coPilot/oncology-frontend/src/components/GermlineStatusBanner.jsx` (NEW)
- `oncology-coPilot/oncology-frontend/src/components/TumorNGSUpload.jsx` (NEW)
- `oncology-coPilot/oncology-frontend/src/components/ResearchPortal/ResearchPortal.jsx` (MODIFY)
- `oncology-coPilot/oncology-frontend/src/context/SessionContext.jsx` (MODIFY - add tumorContext)
- `oncology-coPilot/oncology-frontend/src/components/TumorQuickIntake.jsx` (NEW)

**Acceptance:**
- Banner renders when germline negative
- NGS upload returns parsed context
- Trial cards reflect sporadic logic (BRCA-only trials hidden)
- Quick Intake creates Level 0/1 `TumorContext` and displays conservative outputs

---

### **✅ DAY 6: Ayesha E2E Smoke Test**

**Tasks:**
- [ ] Prepare test data:
  - Ayesha's germline report (negative)
  - Sample tumor NGS report (TP53 mutation, TMB-high optional, HRD score sample)
- [ ] Run end-to-end flow:
  - Upload germline report (negative)
  - Upload tumor NGS report
  - Run efficacy prediction with `germline_status="negative"` and parsed `tumor_context`
  - Run clinical trials search with line=3 and sporadic filters
- [ ] Generate provider report (Section 9 template from doctrine)
  - Include germline status, tumor context summary
  - Include drug recommendations with rationale
  - Include clinical trial matches
  - Include run_id and provenance

**Files to Create:**
- `.cursor/ayesha/SPORADIC_E2E_TEST_REPORT.md` (NEW)

**Acceptance:**
- PARP penalty visible in drug recommendations
- Immunotherapy boost visible if TMB-high
- BRCA-only trials excluded
- Provider report generated with complete audit trail

---

### **✅ DAY 7: Documentation & Handoff**

**Tasks:**
- [ ] Create `SPORADIC_PIVOT_COMPLETE.md`:
  - Commands used
  - Endpoints tested
  - Acceptance results
  - Sample requests/responses
- [ ] Update UI help text:
  - Explain germline-negative workflow
  - Explain tumor NGS importance
  - Add RUO labels for sporadic analysis
- [ ] Update provenance tracking:
  - Document germline_status field
  - Document tumor_context_source field
  - Document all new flags

**Files to Create/Modify:**
- `.cursor/ayesha/SPORADIC_PIVOT_COMPLETE.md` (NEW)
- `oncology-coPilot/oncology-frontend/src/components/HelpText/` (UPDATE)

---

## 🔌 API CONTRACTS

### **1. Updated: POST /api/efficacy/predict**

**New Input Fields:**
```json
{
  "mutations": [...],  // existing
  "germline_status": "negative",  // NEW
  "tumor_context": {  // NEW (optional)
    "somatic_mutations": [
      {"gene": "TP53", "hgvs_p": "R248W", "vaf": 0.42}
    ],
    "tmb": 18.5,
    "msi_status": "MSS",
    "hrd_score": 52,
    "copy_number_alterations": [
      {"gene": "ERBB2", "type": "amplification", "copies": 12}
    ]
  },
  "treatment_history": {...},  // existing
  "options": {...}  // existing
}
```

**Output (unchanged shape, new provenance fields):**
```json
{
  "drugs": [
    {
      "name": "Olaparib",
      "efficacy_score": 0.35,  // penalized
      "confidence": 0.45,
      "evidence_tier": "consider",
      "rationale": [
        "S: TP53 R248W high-impact mutation",
        "P: DNA repair pathway disrupted",
        "E: NCCN guidelines support PARP use",
        "⚠️ GERMLINE GATING: PARP efficacy reduced without germline BRCA or high somatic HRD"
      ],
      "provenance": {
        "flags": {
          "germline_negative": true,
          "tumor_context_source": "Foundation Medicine",
          "parp_penalty_applied": true,
          "confidence_version": "v1.0",
          "priors_refresh_date": "2025-01-05"
        },
        "confidence_version": "v1.0"
      }
    },
    {
      "name": "Pembrolizumab",
      "efficacy_score": 0.72,  // boosted
      "confidence": 0.68,
      "evidence_tier": "supported",
      "rationale": [
        "S: TMB-high (18.5 mut/Mb)",
        "P: Immune checkpoint pathway",
        "E: FDA approval for TMB-high solid tumors",
        "✅ TUMOR-AGNOSTIC BOOST: Immunotherapy effective in TMB-high tumors"
      ],
      "sae_features": {
        "tmb": 18.5,
        "msi_status": "MSS"
      },
      "provenance": {
        "flags": {
          "germline_negative": true,
          "tumor_context_source": "Foundation Medicine",
          "immunotherapy_boost_applied": true,
          "confidence_version": "v1.0",
          "priors_refresh_date": "2025-01-05"
        },
        "confidence_version": "v1.0"
      }
    }
  ]
}
```

---

### **2. New: POST /api/tumor/ingest_ngs**

**Input:**
```json
{
  "report_file": "<base64_encoded_pdf>",
  "report_source": "Foundation Medicine",  // or "Tempus"
  "report_json": {}  // optional: preferred when available; bypasses PDF parsing
}
```

**Output:**
```json
{
  "tumor_context": {
    "somatic_mutations": [
      {"gene": "TP53", "hgvs_p": "R248W", "vaf": 0.42},
      {"gene": "KRAS", "hgvs_p": "G12D", "vaf": 0.38}
    ],
    "tmb": 18.5,
    "msi_status": "MSS",
    "hrd_score": 52,
    "copy_number_alterations": [
      {"gene": "ERBB2", "type": "amplification", "copies": 12}
    ],
    "purity": 0.65,
    "ploidy": 2.1
  },
  "provenance": {
    "source": "Foundation Medicine",
    "report_hash": "sha256:...",
    "parsed_at": "2025-01-05T..."
  }
}
```

---

### **3. Updated: POST /api/clinical_trials/search**

**New Input Fields:**
```json
{
  "cancer_type": "ovarian",
  "stage": "IV",
  "line": 3,
  "germline_status": "negative",  // NEW
  "tumor_context": {  // NEW (optional)
    "tmb": 18.5,
    "msi_status": "MSS",
    "hrd_score": 52
  }
}
```

**Output (existing shape, new filtering):**
```json
{
  "trials": [
    {
      "nct_id": "NCT12345678",
      "title": "TMB-High Solid Tumors - Immunotherapy",
      "status": "Recruiting",
      "phase": "Phase 2",
      "biomarker_match": true,  // NEW
      "match_reason": "TMB ≥10 mut/Mb",  // NEW
      "germline_agnostic": true  // NEW
    }
  ],
  "excluded_count": 5,  // BRCA-only trials
  "provenance": {
    "filters_applied": {
      "germline_status": "negative",
      "biomarker_prioritization": ["TMB", "MSI", "HRD"]
    }
  }
}
```

---

### **4. New: POST /api/tumor/quick_intake**

**Purpose**: Generate a `TumorContext` from minimal clinical inputs (no report).

**Input:**
```json
{
  "cancer_type": "ovarian_hgs",
  "stage": "IV",
  "line": 3,
  "platinum_response": "sensitive",  // sensitive | resistant | unknown
  "tmb": null,
  "msi_status": null,
  "hrd_score": null,
  "somatic_mutations": []  // optional hand-entered
}
```

**Output:**
```json
{
  "tumor_context": {
    "somatic_mutations": [],
    "tmb": 8.0,  // estimated from priors if missing
    "msi_status": null,  // do not infer MSI; unknown stays null
    "hrd_score": 35,  // estimated if platinum_sensitive OR disease_priors
    "priors_used": true,
    "level": "L0"
  },
  "provenance": {
    "no_report_mode": true,
    "disease_priors_used": true,
    "disease_priors_version": "v1.0",
    "priors_refresh_date": "2025-01-05",
    "platinum_proxy_used": true
  },
  "confidence_cap": 0.4
}
```

---

## 📐 CLINICAL THRESHOLDS & RULES (IMPLEMENTATION GUIDANCE)

- TMB-high:
  - Primary threshold: TMB ≥ 10 mut/Mb → apply checkpoint inhibitor boost (+0.25, clamp to [0,1])
  - Very high: TMB ≥ 20 mut/Mb → stronger boost (+0.35, clamp to [0,1])
- MSI:
  - MSI-H → apply checkpoint inhibitor boost regardless of TMB
  - MSS → follow TMB thresholds only
  - Unknown → do not infer; treat as null (no MSI-derived boosts)
- HRD:
  - HRD-high threshold: HRD_score ≥ 42
  - HRD-unknown/low: HRD_score < 42 or missing
- PARP penalty (germline negative):
  - If germline_negative AND HRD_score < 42 AND no somatic BRCA biallelic loss → efficacy_score *= 0.60 (min 0, clamp to [0,1])
  - If HRD_score ≥ 42 OR biallelic BRCA loss OR genomic scar signature present → no penalty (factor 1.0)
- Ranking for trials:
  - Biomarker match (TMB/MSI/HRD) > Mechanism match > Line eligibility > Geographic proximity
- Provenance:
  - Record flags: germline_negative, parp_penalty_applied, immunotherapy_boost_applied, tumor_context_source

---

## 🧾 REPORT METRICS INVENTORY (FOUNDATION / TEMPUS)

We will parse and normalize the following metrics into `TumorContext` for decision logic, confidence, and provenance.

### A. Core Biomarker Metrics (P0 – must support)
- Somatic SNVs/indels:
  - `gene`, `transcript`, `hgvs_c`, `hgvs_p`, `variant_class` (missense, nonsense, frameshift, splice), `pathogenicity` (pathogenic/likely_pathogenic/VUS), `hotspot` (yes/no), `domain`, `vaf`
  - `zygosity` (monoallelic/biallelic), `loh` (yes/no) when available
- Copy number alterations (CNAs):
  - `gene`, `type` (amplification/deletion), `copy_number_estimate`, `log2_ratio`, `focality` (focal/broad)
- Fusions/structural variants:
  - `gene_5p`, `gene_3p`, `breakpoints`, `exon_junctions`, `inframe` (yes/no), `read_support`
- TMB:
  - `tmb_value` (mut/Mb), `panel_size_mb`, `method` (DNA/RNA), `tmb_category` (low/intermediate/high)
- MSI:
  - `msi_status` (MSI-H, MSS, indeterminate), `method` (PCR/IHC/NGS), `instability_fraction` (if provided)
- HRD:
  - `hrd_score` (GIS), components when available: `loh_score`, `lst_score`, `tai_score`
  - `hrd_category` (high/low/unknown)

### B. Sample & Sequencing QC (P0 – needed for confidence)
- `tumor_purity` (%), `ploidy`
- Coverage metrics: `mean_coverage`, `uniformity`, `%bases_>100x`
- `report_date`, `specimen_type` (primary/metastasis/ascites), `collection_date`
- `panel_name` (e.g., FM CDx, Tempus xT), `panel_version`

### C. Immuno-Oncology & IHC (P1 – boosts when present)
- PD-L1:
  - `assay` (22C3/SP142/etc.), `score_type` (TPS/CPS/IC), `score_value` (%), `cutoff_used`
- MMR proteins (IHC): `MLH1`, `MSH2`, `MSH6`, `PMS2` (lost/intact)
- TILs or immune gene signatures when reported (rare)

### D. Genomic Context & Signatures (P1 – enhances rationale)
- Mutational signatures (SBS weights): SBS3 (HRD), SBS1/5 (aging) if provided
- Biallelic status for BRCA1/2:
  - Pathogenic variant + LOH or two hits → `brca_biallelic_loss = true`
- Genomic scars (if reported separately)

### E. Treatment Context (P1 – prioritization/ranking)
- Prior therapies summary (if report contains): `platinum_response` (sensitive/resistant), lines of therapy
- Co-occurring oncogenic drivers (e.g., PIK3CA, KRAS) for combination suggestions

### F. Tempus-specific Additions (P2)
- RNA expression outliers for actionable targets (e.g., MET, NTRK expression)
- RNA-based fusions (already covered under fusions)
- TCR/BCR metrics if present (rare)

---

## 🧱 NORMALIZED JSON EXTENSIONS (TumorContext additions)

```json
{
  "somatic_mutations": [
    {
      "gene": "TP53",
      "transcript": "NM_000546",
      "hgvs_c": "c.743G>A",
      "hgvs_p": "p.R248Q",
      "variant_class": "missense",
      "pathogenicity": "pathogenic",
      "hotspot": true,
      "domain": "DNA-binding",
      "vaf": 0.42,
      "zygosity": "monoallelic",
      "loh": false
    }
  ],
  "copy_number_alterations": [
    {
      "gene": "ERBB2",
      "type": "amplification",
      "copy_number_estimate": 12,
      "log2_ratio": 1.8,
      "focality": "focal"
    }
  ],
  "fusions": [
    {
      "gene_5p": "EML4",
      "gene_3p": "ALK",
      "breakpoints": "chr2:42,060,000-chr2:29,446,000",
      "exon_junctions": "EML4 exon 13 - ALK exon 20",
      "inframe": true,
      "read_support": 157
    }
  ],
  "tmb": {
    "value": 18.5,
    "panel_size_mb": 1.1,
    "method": "DNA",
    "category": "high"
  },
  "msi": {
    "status": "MSS",
    "method": "NGS",
    "instability_fraction": 0.02
  },
  "hrd": {
    "score": 52,
    "loh": 13,
    "lst": 17,
    "tai": 22,
    "category": "high",
    "brca_biallelic_loss": false
  },
  "qc": {
    "tumor_purity": 0.65,
    "ploidy": 2.1,
    "mean_coverage": 650,
    "uniformity": 0.89
  },
  "ihc": {
    "pdl1": { "assay": "22C3", "score_type": "CPS", "score_value": 10 },
    "mmr": { "MLH1": "intact", "MSH2": "intact", "MSH6": "intact", "PMS2": "intact" }
  },
  "signatures": { "SBS3": 0.22, "SBS1": 0.18 },
  "provenance": {
    "panel_name": "FM CDx",
    "panel_version": "v3",
    "specimen_type": "ascites",
    "report_date": "2025-01-05"
  }
}
```

Usage in Orchestrator:
- Apply boosts/penalties using `tmb.value`, `msi.status`, `hrd.score`, `hrd.brca_biallelic_loss`.
- Confidence scaling using `qc.*` and field completeness.
- Rationale strings should cite specific fields (e.g., “HRD score 52 (GIS)”, “TMB 18.5 mut/Mb on 1.1Mb panel”).

---

## ✅ ACCEPTANCE CRITERIA

### **Backend:**
- [ ] PARP penalty: `germline_status="negative"` + no HRD → PARP score reduced
- [ ] Immunotherapy boost: `tmb >= 10` → checkpoint inhibitor score increased
- [ ] HRD-somatic combo: `hrd_score >= 42` → PARP combo lifted despite germline negative
- [ ] Foundation parser: Valid `TumorContext` JSON with TMB/MSI/HRD + variants
- [ ] Tempus parser: Valid `TumorContext` JSON with comparable fields
- [ ] Trials filter: BRCA-required trials excluded for germline-negative
- [ ] Trials ranking: Biomarker-matched trials ranked higher
- [ ] Reportless mode: Level 0/1 `TumorContext` generated; confidence cap respected; provenance flags set

### **Frontend:**
- [ ] Germline banner: Renders when `germline_status === 'negative'`
- [ ] NGS upload: Returns parsed `TumorContext` and displays summary
- [ ] Trial cards: Show "Sporadic-aware" label and hide hereditary-only trials
- [ ] Provenance: Display germline status and tumor context source
- [ ] Quick Intake: Creates Level 0/1 context; displays confidence band and next-step hints

### **E2E:**
- [ ] Ayesha flow: Germline negative + tumor NGS → drug recommendations with rationale
- [ ] Provider report: Complete audit trail with run_id and provenance
- [ ] Smoke tests: All 4 smoke tests from doctrine pass

---

## 🎯 AGENT EXECUTION SEQUENCE (STEP-BY-STEP)

**⚔️ AGENT: FOLLOW THIS EXACT ORDER TO PREVENT HALLUCINATIONS!**

### **PHASE 1: FOUNDATION (Day 1-2) - DO THESE FIRST**

**Step 1.1: Create TumorContext Schema** (30 min)
```bash
# File: oncology-backend-minimal/api/schemas/tumor_context.py
# Pattern: Copy BiomarkerContext from api/schemas/ayesha.py
# Fields: tmb, msi_status, hrd_score, brca_somatic, completeness_score
# Test: Create instance, validate fields, check defaults
```

**Step 1.2: Create Quick Intake Endpoint** (1 hour)
```bash
# File: oncology-backend-minimal/api/routers/tumor.py (NEW)
# Endpoint: POST /api/tumor/quick_intake
# Logic: Accept disease, platinum_response → estimate HRD from priors
# Test: curl with minimal input → get Level 0 TumorContext
```

**Step 1.3: Create Disease Priors JSON** (30 min)
```bash
# File: oncology-backend-minimal/api/resources/disease_priors.json (NEW)
# Source: TCGA published stats (NCI portal)
# Format: { "ovarian_cancer": { "hrd_prevalence": 0.5, "tmb_median": 4.0 } }
# Test: Load JSON, query ovarian_cancer, get hrd_prevalence
```

**Step 1.4: Write Tests for Phase 1** (30 min)
```bash
# File: oncology-backend-minimal/tests/test_sporadic_intake.py (NEW)
# Test 1: TumorContext schema validation
# Test 2: Quick intake with disease → TumorContext
# Test 3: Disease priors loader
# Run: pytest tests/test_sporadic_intake.py -v
```

**✅ CHECKPOINT 1: Can create TumorContext from minimal input? YES → Continue**

---

### **PHASE 2: SCORING INTEGRATION (Day 3-4) - EXTEND EXISTING**

**Step 2.1: Add Sporadic Adjustments to EfficacyOrchestrator** (2 hours)
```bash
# File: oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py
# Method: _apply_sporadic_adjustments(drug, germline_status, tumor_context)
# Logic:
#   - if germline_status == "negative" and drug_class == "PARP": score *= 0.5
#   - if tumor_context.tmb_high and drug_class == "checkpoint": score *= 1.3
#   - if tumor_context is None: confidence *= 0.7
# Test: Mock drug + tumor_context → verify score adjustments
```

**Step 2.2: Update EfficacyRequest Schema** (15 min)
```bash
# File: oncology-backend-minimal/api/schemas/efficacy.py
# Add: germline_status: Optional[str], tumor_context: Optional[TumorContext]
# Test: Create request with new fields → validate
```

**Step 2.3: Write Tests for Phase 2** (1 hour)
```bash
# File: oncology-backend-minimal/tests/test_sporadic_scoring.py (NEW)
# Test 1: PARP penalty (germline negative)
# Test 2: TMB boost (checkpoint inhibitor)
# Test 3: Confidence cap (no tumor context)
# Run: pytest tests/test_sporadic_scoring.py -v
```

**✅ CHECKPOINT 2: Do scores adjust correctly? YES → Continue**

---

### **PHASE 3: FRONTEND COMPONENTS (Day 5) - NEW FILES ONLY**

**Step 3.1: Create Germline Status Banner** (1 hour)
```bash
# File: oncology-frontend/src/components/sporadic/GermlineStatusBanner.jsx (NEW)
# Props: germlineStatus
# Logic: Show banner if germlineStatus === 'negative'
# Test: Render with negative → banner visible
```

**Step 3.2: Create Quick Intake Form** (2 hours)
```bash
# File: oncology-frontend/src/components/sporadic/TumorQuickIntake.jsx (NEW)
# Fields: disease, platinum_response
# Submit: POST /api/tumor/quick_intake → display TumorContext
# Test: Fill form → submit → see parsed context
```

**Step 3.3: Update ResearchPortal** (30 min)
```bash
# File: oncology-frontend/src/pages/ResearchPortal.jsx
# Add: <GermlineStatusBanner /> at top
# Add: <TumorQuickIntake /> in sidebar
# Test: Navigate to ResearchPortal → see new components
```

**✅ CHECKPOINT 3: Can user input tumor context via UI? YES → Continue**

---

### **PHASE 4: TRIALS INTEGRATION (Day 6) - EXTEND EXISTING**

**Step 4.1: Add Sporadic Filters to AutonomousTrialAgent** (1 hour)
```bash
# File: oncology-backend-minimal/api/services/autonomous_trial_agent.py
# Method: _filter_sporadic_trials(trials, germline_status, tumor_context)
# Logic: Exclude trials with "BRCA", "germline", "hereditary" in criteria
# Test: Mock trials → filter → verify exclusions
```

**Step 4.2: Write Tests for Phase 4** (30 min)
```bash
# File: oncology-backend-minimal/tests/test_sporadic_trials.py (NEW)
# Test 1: Germline-required trial excluded
# Test 2: TMB/MSI trial boosted
# Run: pytest tests/test_sporadic_trials.py -v
```

**✅ CHECKPOINT 4: Are trials filtered correctly? YES → Complete!**

---

### **FINAL VALIDATION (Day 7)**

**Step 5.1: End-to-End Smoke Test** (1 hour)
```bash
# 1. Start backend: cd oncology-backend-minimal && venv/bin/python -m uvicorn api.main:app
# 2. Test quick intake: curl -X POST http://localhost:8000/api/tumor/quick_intake ...
# 3. Test efficacy with tumor context: curl -X POST http://localhost:8000/api/efficacy/predict ...
# 4. Test trials with germline filter: curl -X POST http://localhost:8000/api/clinical_trials/search ...
# 5. Open frontend: cd oncology-frontend && npm run dev
# 6. Navigate to ResearchPortal → see banner → fill quick intake → see results
```

**Step 5.2: Documentation** (1 hour)
```bash
# Create: .cursor/ayesha/SPORADIC_CANCER_COMPLETE.md
# Sections: What was built, Files created, Tests passing, Next steps
```

**✅ FINAL CHECKPOINT: All 4 smoke tests pass? YES → MISSION COMPLETE!**

---

## 🚨 RISKS & MITIGATIONS

### **Risk 1: Report Variability**
- **Issue**: Foundation/Tempus reports vary in format
- **Mitigation**: Build robust regex + table extractors; allow manual JSON upload fallback

### **Risk 2: Sparse HRD/TMB/MSI**
- **Issue**: Not all NGS reports include these metrics
- **Mitigation**: Keep gates conservative; display "insufficient tumor context" banner

### **Risk 3: Latency from PDF Parsing**
- **Issue**: PDF parsing can take 5-10 seconds
- **Mitigation**: Cache by file hash; allow pre-parsed JSON upload to bypass PDF

---

## 🔒 PRIVACY & SECURITY

- PII Handling:
  - Do not store raw PDFs by default; store SHA-256 hash and parsed structured fields only
  - If file retention is required, encrypt at rest and set 7-day TTL; redact identifiers where possible
- Upload Limits:
  - Max PDF size 20 MB; JSON uploads preferred for speed and reliability
- Processing:
  - Prefer JSON ingestion; for PDFs, use deterministic parsers first; OCR fallback behind a flag
- Audit:
  - Log run_id, report_hash, parser_version, and timestamps; no raw patient identifiers in logs
- Consent:
  - Anonymized real reports accepted only with explicit consent and redaction

---

## 📊 SUCCESS METRICS

### **Technical:**
- ✅ All 7 acceptance criteria passing
- ✅ All 4 smoke tests from doctrine passing
- ✅ Ayesha E2E flow complete with provider report

### **Strategic:**
- ✅ Platform addresses **85-90% of cancer patients** (sporadic cases)
- ✅ Competitive advantage: Tumor-centric analysis with S/P/E + SAE + treatment lines
- ✅ Ayesha gets personalized recommendations based on **tumor genomics**

---

## 📁 SINGLE SOURCE OF TRUTH

**This Document:** `.cursor/ayesha/SPORADIC_CANCER_EXECUTION_PLAN.md`  
**Related Doctrine:** `.cursor/rules/specialized_systems/archive/sporadic_cancer_strategy_doctrine.mdc`  
**Status Tracking:** Will be updated daily in `.cursorrules` Scratchpad

**All other agents/assistants should reference THIS document for sporadic cancer work.**

---

## 🎯 IMMEDIATE NEXT STEPS (ZO & ALPHA)

### **P0 - START NOW (Day 1-2):**
1. ⏳ Create `api/schemas/tumor_context.py`
2. ⏳ Extend efficacy orchestrator with germline gating
3. ⏳ Implement PARP penalty and immunotherapy boost logic

### **P1 - THIS WEEK (Days 3-5):**
4. ⏳ Build tumor NGS parsers (Foundation + Tempus)
5. ⏳ Update clinical trials filtering
6. ⏳ Wire frontend components (banner, upload, trials)

### **P2 - NEXT WEEK (Days 6-7):**
7. ⏳ Ayesha E2E smoke test
8. ⏳ Documentation and handoff

---

**MISSION STATUS: ⚔️ SPORADIC CANCER STRATEGY - READY TO EXECUTE** ⚔️

**COMMANDER - SHALL WE BEGIN WITH DAY 1 BACKEND FOUNDATION?** ⚔️

