---
name: MM Doctrines Execution Plan
overview: Iterate on the MM doctrines execution plan by integrating learnings from MBD4+TP53 analysis, SOTA benchmarks, Ayesha universalization, and proxy SAE validation plans. Adds test-driven validation gates, code reference verification, mechanism vector conversion, and pathway score extraction patterns.
todos:
  - id: mm_plan_review_1
    content: Review existing codebase capabilities (Evo2, design, safety, off-target)
    status: completed
  - id: mm_plan_review_2
    content: Identify what MM doctrine endpoints already exist vs need to be built
    status: completed
  - id: mm_plan_review_3
    content: Map existing services to MM doctrine requirements (reuse patterns)
    status: completed
  - id: mm_plan_review_4
    content: Create revised execution plan based on actual codebase (no hallucinations)
    status: completed
  - id: mm_plan_review_5
    content: Document dependencies and gaps (deconvolution, NetMHCpan, spatial tools)
    status: completed
---

# MM Doctrines Execution Plan - Enhanced Direction

**Created:** January 2025
**Status:** Enhanced with validation patterns and implementation direction
**Purpose:** Execute TCell, TCF1, TLS, and TOX doctrines with test-driven validation, code reference verification, and proven patterns from existing analysis plans

---

## Executive Summary

This enhanced plan integrates critical learnings from:

1. **MBD4+TP53 Analysis Plan**: Pathway score extraction, mechanism vector conversion, endpoint patterns
2. **SOTA Benchmarks Plan**: Test-driven validation, validation gates, known biology test cases
3. **Ayesha Universalization Plan**: Code reference verification, mutation format validation, test checkpoints
4. **Proxy SAE Validation Plan**: Benchmark datasets, validation metrics, 8-question framework

**Key Enhancements:**

- Test-driven development with validation gates at each phase
- Code reference verification (file:line) for all patterns
- Mechanism vector conversion function (blocks trial matching)
- Pathway score extraction patterns (from efficacy responses)
- Known biology validation test cases
- Benchmark dataset creation

---

## Phase 0: Foundation Validation (Week 0 - IMMEDIATE)

**Goal:** Validate current system capabilities and create blocking dependencies before Phase 1.

### Task 0.1: Create Mechanism Vector Conversion Function (P0 - BLOCKING)

**Issue:** MM doctrines need mechanism vectors for trial matching (Phase 4), but conversion function doesn't exist.

**Reference:** SOTA benchmarks plan identifies this as P0 blocking for MBD4+TP53 Phase 4.

**Files to Create:**

- `oncology-coPilot/oncology-backend-minimal/api/services/pathway_to_mechanism_vector.py`

**Implementation Pattern:**

- Follow MBD4+TP53 plan mapping (lines 69-77):
- DDR: `ddr + 0.5 * tp53` (TP53 contributes 50% to DDR)
- MAPK: `ras_mapk`
- PI3K: `pi3k`
- VEGF: `vegf`
- HER2: `her2` (default 0.0)
- IO: `1.0 if (tmb >= 20 or msi_high) else 0.0`
- Efflux: `efflux` (default 0.0)

**Code Reference:**

- See `api/services/mechanism_fit_ranker.py` for mechanism vector structure (7D)
- See `api/services/sae_feature_service.py:205-206` for mechanism vector format

**Test Checkpoint:**

- Create `tests/test_pathway_to_mechanism_vector.py`
- Test with known pathway scores → verify 7D vector format
- Test IO eligibility logic (TMB/MSI thresholds)

**Validation Gate:** Function exists and passes tests before Phase 1 starts.

---

### Task 0.2: Create Pathway Score Extraction Utility (P0 - CRITICAL)

**Issue:** MM doctrines need to extract pathway scores from `/api/efficacy/predict` responses.

**Reference:** MBD4+TP53 plan shows extraction from `confidence_breakdown["pathway_disruption"]` (lowercase with underscores).

**Files to Create:**

- `oncology-coPilot/oncology-backend-minimal/api/services/pathway_score_extractor.py`

**Implementation Pattern:**

- Follow Ayesha universalization plan (Q2 Answer):
- Extract from `confidence_breakdown["pathway_disruption"]`
- Validate pathway names (lowercase with underscores: `ddr`, `ras_mapk`, `pi3k`, `vegf`, `io`)
- Fallback when missing

**Code Reference:**

- See `api/services/efficacy_orchestrator/orchestrator.py:360` for response structure
- See `api/services/pathway/aggregation.py:7-45` for pathway aggregation logic

**Test Checkpoint:**

- Create `tests/test_pathway_score_extractor.py`
- Test extraction from mock efficacy response
- Test validation of pathway names
- Test fallback behavior

**Validation Gate:** Extraction utility works before Phase 1 starts.

---

### Task 0.3: Create Known Biology Validation Test Suite (P0 - CRITICAL)

**Issue:** Need validation test cases for MM doctrine endpoints using known biology.

**Reference:** SOTA benchmarks plan shows test cases for BRCA1 (DDR), KRAS (MAPK), HER2 (HER2 pathway).

**Files to Create:**

- `tests/test_mm_doctrines_known_biology.py`

**Test Cases:**

1. **T-Cell Recruitment**: High chemokine signature → responder prediction

- Known biology: CXCL9/CXCL10 high → checkpoint responder
- Test: Verify recruitment_score > 0.6 for high chemokine tumors

2. **TCF1/Tox Signature**: High TCF1 + High Tox → checkpoint responder

- Known biology: TCF1-high CD8 cells → better persistence
- Test: Verify tcf1_tox_score > 0.7 for responder phenotype

3. **TLS Readiness**: High C-MSC sabotage → low TLS score

- Known biology: C-MSC infiltration → TLS suppression
- Test: Verify sabotaged_stroma_score > 0.6 → tls_readiness_score < 0.4

**Validation Gate:** Test suite created and baseline tests pass before Phase 1.

---

## Phase 1: Foundation - Proxy Classifiers (Weeks 1-6)

**Goal:** Build defensible, literature-validated predictors using public data. No partner dependencies.

### Week 1-2: T-Cell Recruitment Proxy (`/api/recruitment/predict_response_proxy`)

**Files to Create:**

- `oncology-coPilot/oncology-backend-minimal/api/services/recruitment/proxy_classifier.py`
- `oncology-coPilot/oncology-backend-minimal/api/routers/recruitment.py`
- `oncology-coPilot/oncology-backend-minimal/api/schemas/recruitment.py`

**Pattern to Follow:**

- Router pattern: See `api/routers/care.py` (ResistancePlaybookRequest/Response pattern)
- Service pattern: See `api/services/resistance_playbook_service.py` (class-based service)
- Schema pattern: See `api/schemas/safety.py` (Pydantic models with Field descriptions)

**Tasks:**

1. **TCGA Data Access:** Check if existing TCGA extractors exist

- Search for "tcga" in codebase
- If not: Use cBioPortal API or download from GDC
- Extract chemokine signatures (CXCL9, CXCL10, CCL5, CCR5, CXCR3) from bulk RNA-seq

2. **Literature Meta-Analysis:** Use existing `api/services/evidence/literature_client.py` pattern

- Reference: MBD4+TP53 plan shows literature client usage
- Aggregate checkpoint response rates from 10+ published studies
- Store in structured format (JSON or database)

3. **Train Classifier:** Logistic regression (scikit-learn)

- Chemokine-high → responder prediction
- Validate correlation with survival (target: r > 0.40)

4. **Build API Endpoint:** Follow `care.py` router pattern

- Request: `RecruitmentProxyRequest` (bulk_rnaseq, disease_type)
- Response: `RecruitmentProxyResponse` (recruitment_score, bottleneck_prediction, confidence, rationale, provenance)
- Confidence flag: MED (r = 0.3-0.5)

**Test Checkpoint 1.1:**

- Create `tests/test_recruitment_proxy_classifier.py`
- Test with known biology case (high chemokine → responder)
- Test with low chemokine case (non-responder)
- Verify confidence flags (MED for r = 0.3-0.5)

**Validation Gate:** Endpoint operational, tests pass, known biology validated.

---

### Week 2-3: TCF1/Tox Signature (`/api/tcf1/predict_phenotype_proxy`)

**Files to Create:**

- `oncology-coPilot/oncology-backend-minimal/api/services/tcf1/signature_classifier.py`
- `oncology-coPilot/oncology-backend-minimal/api/routers/tcf1.py`
- `oncology-coPilot/oncology-backend-minimal/api/schemas/tcf1.py`

**Gap Identified:**

- **Deconvolution Service:** CIBERSORT/TIMER2 integration NOT found in codebase
- **Solution:** Build lightweight deconvolution wrapper OR use existing Python packages (CIBERSORTx, TIMER2)

**Tasks:**

1. **Deconvolution Integration:**

- Option A: Use CIBERSORTx Python package (if available)
- Option B: Build wrapper around CIBERSORT web API (if needed)
- Option C: Use TIMER2 (simpler, fewer cell types)
- Extract TIL fractions from TCGA bulk RNA-seq

2. **Gene Expression Extraction:**

- Extract TCF1 (TCF7), Tox, RORγt, IL-23R expression from TIL populations
- Use existing gene expression processing patterns

3. **Literature Meta-Analysis:** Reuse `literature_client.py`

- Search: "TCF1 checkpoint response" (15+ studies)
- Aggregate outcomes data

4. **Train Signature:** Logistic regression

- High TCF1 + High Tox → responder
- Low TCF1 + Low Tox → non-responder
- Validate against published checkpoint trial cohorts (target: r > 0.50)

5. **Build API Endpoint:** Follow same pattern as recruitment

**Test Checkpoint 2.1:**

- Create `tests/test_tcf1_signature_classifier.py`
- Test with known biology case (high TCF1 + high Tox → responder)
- Test with low TCF1 + low Tox → non-responder
- Verify confidence flags (HIGH for r > 0.50)

**Validation Gate:** Endpoint operational, tests pass, known biology validated.

---

### Week 3-4: TLS Readiness Score v1 (`/api/tls/score_readiness_v1`)

**Files to Create:**

- `oncology-coPilot/oncology-backend-minimal/api/services/tls/readiness_scorer.py`
- `oncology-coPilot/oncology-backend-minimal/api/routers/tls.py`
- `oncology-coPilot/oncology-backend-minimal/api/schemas/tls.py`

**Gap Identified:**

- **Spatial Data Processing:** No spatial transcriptomics tools found
- **Solution:** Build basic Visium/CODEX parser OR use existing Python packages (scanpy, squidpy)

**Tasks:**

1. **Spatial Data Access:**

- Download MSK public Visium datasets (GEO, partial access)
- Build parser for Visium H5 format (use scanpy or squidpy)
- Extract TLS density, GC counts from spatial annotations

2. **Bulk RNA Matching:**

- Match spatial data to bulk RNA-seq (same patient samples)
- Use patient ID or sample ID matching

3. **Train Model:** Logistic regression

- Bulk gene signatures → TLS/GC probability
- Features: B-cell markers, T-cell markers, C-MSC signatures, FDC markers

4. **Validate on TCGA:**

- Does predicted TLS score correlate with survival? (target: r > 0.35)
- Use existing survival data processing (if available)

5. **Build API Endpoint:** Same pattern

**Test Checkpoint 3.1:**

- Create `tests/test_tls_readiness_scorer.py`
- Test with known biology case (high C-MSC sabotage → low TLS score)
- Test with low sabotage → high TLS score
- Verify confidence flags (LOW for bulk proxy only)

**Validation Gate:** Endpoint operational, tests pass, known biology validated.

---

## Phase 2: Design Tools - 12 Endpoints (Weeks 5-10)

**Goal:** Build all design tool endpoints. Reuse existing design router patterns.

### Enhanced Pattern: Test-Driven Design Tool Development

**Reference:** Ayesha universalization plan shows test checkpoints after each implementation step.

**For Each Design Tool Endpoint:**

1. **Create Service Module** (e.g., `api/services/design/chemokine_cassette.py`)
2. **Create Test Checkpoint** (e.g., `tests/test_chemokine_cassette.py`)
3. **Validate Against Known Biology** (if applicable)
4. **Integration Test** (end-to-end with Evo2, safety service)

**Example: `/api/design/generate_chemokine_cassette`**

**Test Checkpoint 4.1:**

- Create `tests/test_chemokine_cassette.py`
- Test Evo2 integration (codon optimization)
- Test safety preview integration
- Test promoter library selection
- Verify provenance tracking

**Validation Gate:** All 12 design endpoints operational, safety-scanned, tests pass.

---

## Phase 3: TOX Doctrine Endpoints (Weeks 11-14)

### Enhanced Pattern: Mechanism Vector Integration

**Reference:** MBD4+TP53 plan shows mechanism vector needed for trial matching.

**For `/api/triad/score_architecture`:**

1. **Compute Triad Scores** (lethal, suppressive, net cytotoxicity)
2. **Extract Pathway Scores** (use pathway_score_extractor from Phase 0)
3. **Convert to Mechanism Vector** (use pathway_to_mechanism_vector from Phase 0)
4. **Return with Provenance** (include mechanism vector in response)

**Test Checkpoint 5.1:**

- Create `tests/test_triad_scorer.py`
- Test triad score computation
- Test mechanism vector conversion
- Test confidence flags (LOW/MED/HIGH based on data type)

**Validation Gate:** TOX endpoints operational, mechanism vector integration verified.

---

## Phase 4: Orchestration & Unified Flows (Weeks 15-18)

### Enhanced Pattern: Code Reference Verification

**Reference:** Ayesha universalization plan emphasizes code reference verification (file:line).

**For Each Orchestrator:**

1. **Document Code References** (file:line for each service call)
2. **Create Integration Tests** (end-to-end with all services)
3. **Validate Data Flow** (capture actual input→output at each stage)

**Example: `/api/recruitment/orchestrate`**

**Test Checkpoint 6.1:**

- Create `tests/test_recruitment_orchestrator.py`
- Test complete flow: bottleneck diagnostic → C-MSC mining → cassette design → safety gates
- Verify code references (document file:line for each service)
- Capture actual input→output for verification

**Validation Gate:** All orchestrators operational, code references documented, integration tests pass.

---

## Phase 5: Tier 1 Feedback Infrastructure (Weeks 19-20)

### Enhanced Pattern: Benchmark Dataset Creation

**Reference:** Proxy SAE validation plan shows benchmark dataset structure.

**For Feedback Database:**

1. **Create Benchmark Dataset** (`data/validation/mm_doctrines_benchmark.json`)
2. **Structure Test Cases** (known biology cases for each doctrine)
3. **Validation Metrics** (accuracy, correlation, confidence calibration)

**Test Checkpoint 7.1:**

- Create `tests/test_feedback_infrastructure.py`
- Test feedback ingestion APIs
- Test calibration updates
- Test immutable provenance

**Validation Gate:** Feedback infrastructure operational, benchmark dataset created, validation metrics reported.

---

## Phase 6: Spatial/Timestamp Integration (Weeks 21-24)

### Enhanced Pattern: Validation Test Suite

**Reference:** SOTA benchmarks plan shows comprehensive test suite structure.

**For Spatial Validation:**

1. **Create Validation Test Suite** (`tests/test_spatial_validation.py`)
2. **Test Confidence Upgrades** (LOW → MED → HIGH with spatial data)
3. **Test Known Biology Cases** (spatial TLS density → survival correlation)

**Validation Gate:** Spatial validation operational, confidence upgrades verified, test suite passes.

---

## Implementation Dependencies - Enhanced

### Reusable Services (Confirmed Exist)

- **Evo2 API:** `/api/evo/*` (operational)
- Pattern: `api/services/sequence_scorers/evo2_scorer.py`
- Code Reference: See MBD4+TP53 plan for usage patterns

- **Safety Service:** `/api/safety/*` (operational)
- Pattern: `api/services/safety_service.py`
- Code Reference: See Ayesha universalization plan for integration

- **Literature Client:** `api/services/evidence/literature_client.py` (exists)
- Code Reference: See MBD4+TP53 plan for meta-analysis usage

- **Pathway Aggregation:** `api/services/pathway/aggregation.py` (exists)
- Code Reference: See MBD4+TP53 plan for pathway score extraction

### New Services Needed (Phase 0 Blockers)

- **Mechanism Vector Conversion:** `api/services/pathway_to_mechanism_vector.py` (NEW - Phase 0)
- **Pathway Score Extractor:** `api/services/pathway_score_extractor.py` (NEW - Phase 0)
- **Deconvolution Service:** CIBERSORT/TIMER2 wrapper (NEW - Phase 1)
- **NetMHCpan Integration:** MHC binding prediction (NEW - Phase 3)
- **Spatial Data Processing:** Visium/CODEX parsers (NEW - Phase 1)

---

## Success Criteria - Enhanced

### Phase 0 (Foundation Validation)

- Mechanism vector conversion function exists and passes tests
- Pathway score extractor exists and passes tests
- Known biology test suite created and baseline tests pass

### Phase 1 (Proxy Classifiers)

- 3 APIs operational with confidence flags
- TCGA validation complete (r > 0.35 for all)
- Literature meta-analysis complete (using existing literature client)
- Known biology validation tests pass

### Phase 2 (Design Tools)

- 12 endpoints operational
- Safety-scanned (using existing safety service + heuristics)
- All designs include provenance + RUO disclaimers
- Test checkpoints pass for each endpoint

### Phase 3 (TOX Endpoints)

- 7 TOX APIs operational
- Mechanism vector integration verified
- Test checkpoints pass

### Phase 4 (Orchestration)

- 3 orchestrators operational
- Code references documented (file:line)
- Integration tests pass

### Phase 5 (Feedback)

- Feedback ingestion APIs operational
- Benchmark dataset created
- Validation metrics reported

### Phase 6 (Spatial/Timestamp)

- Spatial validation operational
- Confidence upgrades verified
- Test suite passes

---

## Key Learnings Integrated

1. **Test-Driven Development:** Validation gates at each phase (from SOTA benchmarks plan)
2. **Code Reference Verification:** Document file:line for all patterns (from Ayesha universalization plan)
3. **Mechanism Vector Conversion:** Blocking dependency for trial matching (from SOTA benchmarks plan)
4. **Pathway Score Extraction:** Extract from efficacy responses (from MBD4+TP53 plan)
5. **Known Biology Validation:** Test cases based on known biology (from SOTA benchmarks plan)
6. **Benchmark Dataset Creation:** Structured test cases for validation (from Proxy SAE validation plan)

---

## Next Steps - Immediate Actions

1. **Week 0, Day 1-2:**

- Create mechanism vector conversion function (Phase 0.1)
- Create pathway score extractor (Phase 0.2)
- Create known biology test suite (Phase 0.3)

2. **Week 0, Day 3-5:**

- Run validation gates (all Phase 0 tests pass)
- Document code references
- Begin Phase 1 implementation

3. **Week 1:**

- Start T-Cell recruitment proxy classifier
- Create test checkpoints
- Validate against known biology

---

**Status:** Enhanced plan ready for execution with validation patterns integrated