<!-- 89fd491b-919d-4434-ae5c-80fc25d31024 05b705cf-9d10-46a0-9a4b-33a3a27e83a9 -->
# Advanced Trial Query System for Complex Clinical Questions

## Executive Summary

Enhance the existing clinical trials system to answer complex multi-criteria queries (e.g., "MBD4 + DNA repair + basket trials + PARP inhibitors") by combining semantic search improvements with direct ClinicalTrials.gov API query building. Leverage existing scripts and services without reinventing.

**Manager Policy Compliance**: This plan implements Manager's P3 (Gemini trial tagging) and P4 (Mechanism fit ranking) from `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md`.

**Critical Manager Requirements**:

- **P3**: Gemini tagging OFFLINE ONLY (never runtime), batch tag 200 trials, human review 30, ≥90% accuracy required
- **P4**: Mechanism fit formula α=0.7 (eligibility) + β=0.3 (mechanism_fit), thresholds: eligibility ≥0.60, mechanism_fit ≥0.50
- **C7**: Patient mechanism vector = [DDR, MAPK, PI3K, VEGF, IO, Efflux] (6D, not 7D - note: HER2 not in manager spec)
- **C7**: Fallback: if patient vector all zeros, disable mechanism_fit (β=0), show "awaiting NGS" message
- **P4**: If mechanism_fit <0.50, show trial but without mechanism boost + add "low mechanism fit" warning
- **P4**: Never suppress SOC; SOC card remains first-class
- **P4**: Provide "Show all trials" toggle for clinician control

## Current Capabilities Analysis

### Existing Components

1. **Autonomous Trial Agent** (`api/services/autonomous_trial_agent.py`)

   - Generates 1-3 simple queries (disease + biomarker)
   - Uses semantic search via HybridTrialSearchService
   - Limited to basic patterns: "disease biomarker trial", "disease category clinical trial"

2. **ClinicalTrials.gov API Integration** (`api/routers/clinical_trials.py`)

   - Basic search: `search_clinical_trials(gene, cancer_type, max_results)`
   - Simple query building: `query.term = gene + cancer_type`
   - Status filter: `filter.overallStatus = RECRUITING`

3. **Extraction Scripts**

   - `extract_fresh_recruiting_trials.py`: Hardcoded ovarian cancer, RECRUITING only
   - `seed_trials_table.py`: Uses `fetch_ovarian_trials()` with fixed filters
   - `reconnaissance_ovarian_trials.py`: Count queries for ovarian trials

4. **Trial Intelligence Pipeline** (`api/services/trial_intelligence_universal/`)

   - 6-stage filtering (hard filters, trial type, location, eligibility, LLM analysis, dossier)
   - Works on pre-fetched candidates
   - Not designed for initial query generation

## Requirements

Answer queries like:

- "MBD4-associated neoplasia + DNA repair deficiency + basket trials"
- "TP53-mutant ovarian cancer + HRD-positive + BRCA-wildtype"
- "Platinum-sensitive/resistant + rare germline DNA repair mutations"
- "PARP inhibitors + ATR/ATM/DNA-PK inhibitors + checkpoint inhibitors"
- "Basket trials + rare mutation registries + precision medicine protocols"

Output: Trial name, phase, therapy type, NCT ID, location, enrollment criteria, PI contact.

## Implementation Plan

### Phase 1: Enhance Autonomous Agent Query Generation (2-3 hours)

**File**: `api/services/autonomous_trial_agent.py`

**Changes**:

1. **Expand `generate_search_queries()` method**:

   - Add support for multiple biomarkers (not just first one)
   - Generate queries for: basket trials, rare disease registries, precision medicine
   - Add intervention-specific queries (PARP, ATR/ATM, checkpoint inhibitors)
   - Add DNA repair pathway queries
   - Generate 5-10 queries instead of 3

2. **Add query templates**:
   ```python
   QUERY_TEMPLATES = {
       'basket_trial': "{condition} basket trial tumor agnostic",
       'rare_mutation': "{gene} mutation rare disease registry",
       'dna_repair': "{condition} DNA repair deficiency syndrome",
       'parp_inhibitor': "{condition} PARP inhibitor",
       'checkpoint_inhibitor': "{condition} PD-1 PD-L1 checkpoint inhibitor",
       'precision_medicine': "{condition} precision medicine protocol",
       'synthetic_lethal': "{gene} synthetic lethal targeted agent"
   }
   ```

3. **Enhance `extract_patient_context()`**:

   - Detect DNA repair pathway mutations (MBD4, BRCA1/2, TP53, HRD)
   - Extract intervention preferences (PARP, ATR/ATM, immunotherapy)
   - Identify rare mutation status
   - Extract platinum sensitivity status
   - **NEW**: Efficacy prediction integration
     - Extract pathway scores from efficacy_predictions (if available)
     - Convert pathway scores to mechanism vector (6D or 7D) for mechanism fit ranking
     - Extract top-ranked drugs from efficacy_predictions
     - Auto-infer interventions from predicted-effective drugs (use DRUG_MECHANISM_DB)
   - **NEW**: Sporadic cancer support (already exists in autonomous_trial_agent.py)
     - Extract germline_status and tumor_context (if provided)
     - Use for filtering and biomarker boosting

**Key Functions to Modify**:

- `generate_search_queries()`: Lines 64-93
- `extract_patient_context()`: Lines 26-62
- Add new method: `_generate_advanced_queries(patient_context) -> List[str]`

### Phase 2: Create Direct API Query Builder (3-4 hours)

**New File**: `api/services/ctgov_query_builder.py`

**Purpose**: Build precise ClinicalTrials.gov API v2 queries with multiple filters

**Key Features**:

1. **Query Builder Class**:
   ```python
   class CTGovQueryBuilder:
       def __init__(self):
           self.params = {}
       
       def add_condition(self, condition: str, operator: str = "AND")
       def add_intervention(self, intervention: str, operator: str = "OR")
       def add_status(self, status: List[str])  # RECRUITING, NOT_YET_RECRUITING
       def add_phase(self, phases: List[str])  # PHASE1, PHASE2, PHASE3
       def add_study_type(self, study_type: str)  # INTERVENTIONAL
       def add_geo(self, country: str = "United States")
       def add_keyword(self, keyword: str)  # "basket", "rare disease", "precision medicine"
       def build(self) -> Dict[str, Any]
   ```

2. **Specialized Query Methods**:

   - `build_dna_repair_query(conditions, mutations, interventions)`
   - `build_basket_trial_query(conditions, mutations)`
   - `build_rare_mutation_query(gene, condition)`
   - `build_immunotherapy_query(condition, dna_repair_allowed=True)`

3. **Query Execution**:

   - `async def execute_query(builder: CTGovQueryBuilder, max_results: int = 1000) -> List[Dict]`
   - Handle pagination (page tokens)
   - Rate limiting (2 req/sec)
   - Deduplication by NCT ID

**Integration Points**:

- Use existing `CLINICAL_TRIALS_BASE_URL` from `api/routers/clinical_trials.py`
- Reuse pagination logic from `extract_fresh_recruiting_trials.py`

### Phase 3: Parameterize Extraction Scripts (2-3 hours)

**Files to Modify**:

1. `scripts/extract_fresh_recruiting_trials.py`
2. `scripts/seed_trials_table.py`

**Changes**:

1. **Add CLI arguments**:
   ```python
   parser.add_argument('--condition', type=str, help='Condition query (e.g., "ovarian cancer")')
   parser.add_argument('--intervention', type=str, help='Intervention query (e.g., "PARP inhibitor")')
   parser.add_argument('--status', nargs='+', default=['RECRUITING'], help='Status filters')
   parser.add_argument('--phase', nargs='+', help='Phase filters (PHASE1, PHASE2, PHASE3)')
   parser.add_argument('--study-type', type=str, default='INTERVENTIONAL', help='Study type')
   parser.add_argument('--keyword', type=str, help='Additional keyword (e.g., "basket", "rare disease")')
   parser.add_argument('--limit', type=int, default=1000, help='Max trials to fetch')
   ```

2. **Refactor `fetch_recruiting_ovarian_trials()`**:

   - Rename to `fetch_trials_by_criteria()`
   - Accept query builder parameters
   - Remove hardcoded "ovarian cancer" filter

3. **Update `extract_and_seed()`**:

   - Use parameterized query builder
   - Support multiple conditions (comma-separated or list)
   - Support intervention filters

**Usage Examples**:

```bash
# DNA repair trials
python scripts/extract_fresh_recruiting_trials.py \
  --condition "ovarian cancer" \
  --intervention "PARP inhibitor" \
  --keyword "DNA repair" \
  --status RECRUITING NOT_YET_RECRUITING \
  --phase PHASE1 PHASE2 PHASE3 \
  --limit 500

# Basket trials
python scripts/extract_fresh_recruiting_trials.py \
  --condition "ovarian cancer" \
  --keyword "basket trial tumor agnostic" \
  --status RECRUITING \
  --limit 200
```

### Phase 4: Create Advanced Query Endpoint (2-3 hours)

**New File**: `api/routers/advanced_trial_queries.py`

**Purpose**: REST endpoint for complex multi-criteria queries

**Endpoint**: `POST /api/trials/advanced-query`

**Request Schema**:

```python
class AdvancedTrialQueryRequest(BaseModel):
    conditions: List[str]  # ["ovarian cancer", "high-grade serous"]
    mutations: List[str] = []  # ["MBD4", "TP53"]
    interventions: List[str] = []  # ["PARP inhibitor", "checkpoint inhibitor"]
    status: List[str] = ["RECRUITING", "NOT_YET_RECRUITING"]
    phases: List[str] = ["PHASE1", "PHASE2", "PHASE3"]
    study_type: str = "INTERVENTIONAL"
    keywords: List[str] = []  # ["basket trial", "rare disease", "DNA repair"]
    geo: Optional[str] = "United States"
    max_results: int = 100
    use_semantic_search: bool = True  # Also use AstraDB semantic search
    enable_mechanism_fit: bool = True  # Rank by mechanism alignment (pathway-based)
    mechanism_vector: Optional[List[float]] = None  # 6D or 7D pathway vector [DDR, MAPK, PI3K, VEGF, IO, Efflux] or [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux] ⚠️ NEED CLARIFICATION: Manager C7 says 6D, plan uses 7D
    # NEW: Efficacy prediction integration
    efficacy_predictions: Optional[Dict[str, Any]] = None  # S/P/E framework output from /api/efficacy/predict
    pathway_scores: Optional[Dict[str, float]] = None  # Pathway disruption scores (e.g., {"ddr": 0.85, "ras_mapk": 0.20})
    # NEW: Sporadic cancer support (already in autonomous_trial_agent.py)
    germline_status: Optional[str] = None  # "positive", "negative", "unknown" (for sporadic filtering)
    tumor_context: Optional[Dict[str, Any]] = None  # {tmb, hrd_score, msi_status} (for biomarker boost)
    auto_infer_interventions: bool = True  # Auto-add interventions from top-ranked drugs in efficacy_predictions
```

**Response Schema**:

```python
class AdvancedTrialQueryResponse(BaseModel):
    success: bool
    total_found: int
    trials: List[TrialDetail]
    query_method: str  # "api_direct", "semantic", "hybrid"
    provenance: Dict
```

**Implementation**:

1. Use `CTGovQueryBuilder` for direct API queries
2. Use enhanced `AutonomousTrialAgent` for semantic search (with germline_status and tumor_context)
3. Merge and deduplicate results
4. Apply trial intelligence pipeline filtering (optional)
5. **NEW**: Efficacy prediction integration

   - If efficacy_predictions provided:
     - Extract top-ranked drugs (e.g., olaparib, niraparib)
     - Look up MoA from DRUG_MECHANISM_DB (api/services/client_dossier/dossier_generator.py)
     - Map MoA to intervention keywords (e.g., "PARP inhibitor", "checkpoint inhibitor")
     - Auto-add to interventions list if auto_infer_interventions=True
   - If pathway_scores provided:
     - Convert to 6D or 7D mechanism vector (see Phase 4.5 for conversion logic)
     - Use for mechanism fit ranking

6. **NEW**: Apply mechanism fit ranking if mechanism_vector provided (see Phase 4.5)
7. Return structured results with PI contact info

**Trial Detail Schema**:

```python
class TrialDetail(BaseModel):
    nct_id: str
    title: str
    phase: str
    status: str
    therapy_types: List[str]  # ["PARP inhibitor", "checkpoint inhibitor"]
    locations: List[LocationInfo]
    enrollment_criteria: str
    genetic_requirements: List[str]
    principal_investigator: Optional[PIInfo]
    site_contact: Optional[ContactInfo]
    source_url: str
    mechanism_fit_score: Optional[float] = None  # Mechanism alignment score (0-1)
    combined_score: Optional[float] = None  # 0.7 × eligibility + 0.3 × mechanism_fit
    mechanism_alignment: Optional[Dict[str, float]] = None  # Per-pathway alignment breakdown
    low_mechanism_fit_warning: Optional[bool] = None  # True if mechanism_fit <0.50 (Manager P4)
    mechanism_boost_applied: Optional[bool] = None  # True if mechanism_fit ≥0.50
```

### Phase 4.5: Integrate Mechanism Fit Ranking (1-2 hours)

**File**: `api/routers/advanced_trial_queries.py` (add to Phase 4)

**Purpose**: Rank trials by mechanism alignment using pathway-based mechanism vectors

**Manager P4 Compliance**:

- Formula: combined_score = (0.7 × eligibility) + (0.3 × mechanism_fit)
- Minimum eligibility to enter top-10: ≥0.60
- Minimum mechanism_fit for mechanism-gated display: ≥0.50
- If mechanism_fit <0.50: show trial but without mechanism boost + add "low mechanism fit" warning
- Never suppress SOC; SOC card remains first-class
- Provide "Show all trials" toggle for clinician control

**Manager C7 Compliance**:

- Fallback: if patient vector all zeros/unknown, disable mechanism_fit (β=0)
- Show explanation: "awaiting NGS; eligibility-only ranking shown"
- Show breakdown: "DDR 0.82 × PARP+ATR → 0.95 fit"

**Integration Points**:

1. **Extract Mechanism Vector from Patient Data**:

   - Option 1: From /api/efficacy/predict response
     - Extract pathway_scores from `provenance["confidence_breakdown"]["pathway_disruption"]`
     - Map to 6D or 7D mechanism vector: [DDR, MAPK, PI3K, VEGF, IO, Efflux] or [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
     - IO = 1.0 if TMB ≥20 OR MSI-High, else 0.0
   - Option 2: Compute from mutations directly
     - Use pathway/aggregation.py to compute pathway burden
     - Build mechanism vector from gene→pathway mapping
   - Option 3: Accept mechanism_vector in request (if already computed)

2. **Apply Mechanism Fit Ranking**:
   ```python
   from api.services.mechanism_fit_ranker import MechanismFitRanker
   
   # Check if mechanism vector is all zeros (Manager C7 fallback)
   if not mechanism_vector or all(v == 0.0 for v in mechanism_vector):
       # Disable mechanism fit (β=0), use eligibility-only ranking
       ranked_trials = rank_trials_by_eligibility_only(trials_from_search)
       mechanism_fit_disabled = True
       mechanism_fit_message = "awaiting NGS; eligibility-only ranking shown"
   else:
       ranker = MechanismFitRanker(alpha=0.7, beta=0.3)  # Manager's P4 formula
       ranked_trials = ranker.rank_trials(
           trials=trials_from_search,
           sae_mechanism_vector=mechanism_vector,  # 6D or 7D pathway vector
           min_eligibility=0.60,  # Manager's P4 threshold
           min_mechanism_fit=0.50  # Manager's P4 threshold
       )
       mechanism_fit_disabled = False
   
   # Add "low mechanism fit" warning for trials with mechanism_fit <0.50 (Manager P4)
   for trial in ranked_trials:
       if trial.get('mechanism_fit_score', 0) < 0.50:
           trial['low_mechanism_fit_warning'] = True
           trial['mechanism_boost_applied'] = False
   ```

3. **Update Trial Results**:

   - Add mechanism_fit_score to each trial
   - Add combined_score (0.7 × eligibility + 0.3 × mechanism_fit)
   - Add mechanism_alignment (per-pathway breakdown)
   - Re-sort trials by combined_score (descending)

**Code References**:

- Mechanism Fit Ranker: `api/services/mechanism_fit_ranker.py`
- Existing Integration: `api/routers/ayesha_trials.py` (lines 597-631)
- Manager Policy: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md` (P4, C7)
- Pathway to Mechanism Vector: `api/services/pathway_to_mechanism_vector.py` (NEW - to be created)

**Benefits**:

- Better trial ranking for MBD4+TP53 cases (high DDR pathway → better PARP trial matching)
- Pathway-based mechanism vectors work now (no SAE required)
- Future: Automatically uses TRUE SAE vectors when Feature→Pathway Mapping complete

### Phase 5: Enhance Trial Data Extraction (2-3 hours)

**File**: `api/services/client_dossier/trial_scraper.py` or new `api/services/trial_data_enricher.py`

**Purpose**: Extract PI contact info, enrollment criteria, genetic requirements from trial data

**Key Functions**:

1. **Extract PI Information**:

   - Parse `contactsLocationsModule` from ClinicalTrials.gov API response
   - Extract: PI name, email, institution, phone
   - Priority: Overall PI > Study Director > Study Chair

2. **Extract Enrollment Criteria**:

   - Parse `eligibilityModule`
   - Extract inclusion/exclusion criteria
   - Identify genetic requirements (BRCA, HRD, MBD4, etc.)
   - Extract biomarker requirements

3. **Extract Therapy Types**:

   - Parse `interventionsModule`
   - Classify: PARP inhibitor, ATR/ATM inhibitor, checkpoint inhibitor, etc.
   - Use existing `DRUG_MECHANISM_DB` from `dossier_generator.py`
   - **NEW**: Mechanism Vector Tagging (required for MechanismFitRanker)
     - **⚠️ MANAGER P3 COMPLIANCE**: Gemini tagging OFFLINE ONLY (never runtime)
       - Batch tag 200 ovarian trials → human spot-review 30 diverse trials
       - Accept batch if ≥90% tag accuracy; otherwise adjust prompt and re-tag
       - Persist metadata: model, version, parsed_at, reviewed_by, source_checksum
       - Update cadence: weekly diff for new/changed trials
       - Uncertain tags default to neutral vector; never force a mechanism label
     - Extract MoA from interventionsModule (runtime fallback only if Gemini tag missing)
     - Map to mechanism vector: [DDR, MAPK, PI3K, VEGF, IO, Efflux] (6D per Manager C7) or [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux] (7D - ⚠️ NEED CLARIFICATION)
     - Use DRUG_MECHANISM_DB for drug → MoA mapping (runtime fallback)
     - Store moa_vector in trial data with versioning metadata (required for MechanismFitRanker.rank_trials())
     - Priority: 1) Pre-tagged Gemini vectors (offline), 2) Runtime keyword matching (fallback only)

4. **Extract Location Details**:

   - Parse `contactsLocationsModule.locations`
   - Extract: facility, city, state, country
   - Format for display

**Integration**:

- Use existing `get_scraped_data_from_sqlite()` pattern
- Cache enriched data in SQLite `scraped_data_json` field
- Reuse parsing logic from `agent_1_seeding/parsers/study_parser.py`

### Phase 6: Testing & Validation (2-3 hours)

**Test Cases**:

1. **MBD4 Query**:
   ```python
   query = {
       "conditions": ["ovarian cancer", "high-grade serous"],
       "mutations": ["MBD4"],
       "keywords": ["DNA repair", "basket trial"],
       "interventions": ["PARP inhibitor"],
       "status": ["RECRUITING", "NOT_YET_RECRUITING"]
   }
   ```

2. **TP53 + HRD Query**:
   ```python
   query = {
       "conditions": ["ovarian cancer"],
       "mutations": ["TP53"],
       "keywords": ["HRD-positive", "BRCA-wildtype"],
       "interventions": ["checkpoint inhibitor"],
       "phases": ["PHASE2", "PHASE3"]
   }
   ```

3. **Platinum Sensitivity Query**:
   ```python
   query = {
       "conditions": ["ovarian cancer"],
       "keywords": ["platinum-sensitive", "platinum-resistant"],
       "interventions": ["PARP inhibitor", "ATR inhibitor"],
       "status": ["RECRUITING"]
   }
   ```

4. **Efficacy Prediction Integration Test**:
   ```python
   query = {
       "conditions": ["ovarian cancer"],
       "mutations": ["BRCA1"],
       "efficacy_predictions": {
           "drugs": [
               {"name": "olaparib", "efficacy": 0.85, "confidence": 0.75},
               {"name": "niraparib", "efficacy": 0.80, "confidence": 0.70}
           ],
           "provenance": {
               "confidence_breakdown": {
                   "pathway_disruption": {
                       "ddr": 0.85,
                       "ras_mapk": 0.20,
                       "pi3k": 0.10
                   }
               }
           }
       },
       "auto_infer_interventions": True,
       "enable_mechanism_fit": True
   }
   ```

5. **Sporadic Cancer Test**:
   ```python
   query = {
       "conditions": ["ovarian cancer"],
       "mutations": ["TP53"],
       "germline_status": "negative",  # Sporadic cancer
       "tumor_context": {
           "tmb": 10.5,
           "hrd_score": 0.65,
           "msi_status": "MSS"
       }
   }
   ```


**Validation**:

- Compare results from direct API vs semantic search
- Verify PI contact info extraction
- Validate genetic requirements parsing
- Test pagination for large result sets
- Measure performance (queries per second, total time)
- **NEW**: Mechanism fit ranking validation
  - Verify top trials have high mechanism_fit_score ≥0.50
  - Verify combined_score = (0.7 × eligibility) + (0.3 × mechanism_fit)
  - Verify mechanism_alignment breakdown per pathway
  - Verify "low mechanism fit" warning for trials with mechanism_fit <0.50
  - Verify fallback logic when mechanism vector is all zeros
- **NEW**: Efficacy prediction integration validation
  - Validate auto-inference of interventions from top-ranked drugs
  - Verify interventions auto-added from efficacy_predictions (e.g., "PARP inhibitor" from olaparib)
  - Verify pathway score conversion to mechanism vector
- **NEW**: Sporadic cancer validation
  - Validate sporadic cancer filtering (verify germline_status and tumor_context used)
  - Verify biomarker boosting for high TMB/MSI-H cases

## File Changes Summary

### New Files

1. `api/services/ctgov_query_builder.py` - Query builder for ClinicalTrials.gov API
2. `api/routers/advanced_trial_queries.py` - Advanced query endpoint
3. `api/services/trial_data_enricher.py` - PI contact and criteria extraction
4. `api/services/pathway_to_mechanism_vector.py` - Convert pathway_scores to 6D or 7D mechanism vector (NEW)

### Modified Files

1. `api/services/autonomous_trial_agent.py` - Enhanced query generation + efficacy integration + sporadic cancer support
2. `scripts/extract_fresh_recruiting_trials.py` - Parameterized extraction
3. `scripts/seed_trials_table.py` - Parameterized seeding
4. `api/main.py` - Register new router
5. `api/services/mechanism_fit_ranker.py` - Verify integration points (existing, may need minor updates)

### Dependencies

- Existing: `httpx`, `requests`, `astrapy`, `google-generativeai`
- No new dependencies required

## Success Criteria

1. Can answer complex queries like "MBD4 + DNA repair + basket trials"
2. Returns structured data with PI contact info
3. Supports both direct API queries and semantic search
4. Extraction scripts support custom parameters
5. Mechanism fit ranking improves trial matching (PARP trials rank higher for MBD4+TP53 cases)
6. Performance: < 10 seconds for queries returning < 100 trials
7. Handles pagination for large result sets (> 1000 trials)
8. **NEW**: Efficacy prediction integration

   - Auto-inference of interventions from top-ranked drugs works
   - Pathway score conversion to mechanism vector works correctly

9. **NEW**: Sporadic cancer support

   - germline_status and tumor_context filtering works
   - Biomarker boosting for high TMB/MSI-H cases works

10. **NEW**: Manager P3/P4 compliance

    - Gemini tagging OFFLINE ONLY (never runtime)
    - Mechanism fit formula and thresholds match Manager P4
    - "Low mechanism fit" warning for trials with mechanism_fit <0.50
    - Fallback logic when mechanism vector is all zeros

## Implementation Order

1. **Phase 1** (2-3h): Enhance autonomous agent - enables better semantic search
2. **Phase 2** (3-4h): Create query builder - enables direct API queries
3. **Phase 3** (2-3h): Parameterize scripts - enables flexible extraction
4. **Phase 4** (2-3h): Create endpoint - exposes functionality via API
5. **Phase 4.5** (1-2h): Integrate mechanism fit ranking - improves trial matching (Manager P4 compliance)
6. **Phase 5** (2-3h): Enhance data extraction - enriches results (Manager P3 compliance)
7. **Phase 6** (2-3h): Testing - validates functionality

**Total Estimated Time**: 14-21 hours (2-3 days)

## Notes

- Leverage existing `agent_1_seeding` parsers and API client patterns
- Reuse pagination logic from `extract_fresh_recruiting_trials.py`
- Use existing trial intelligence pipeline for post-filtering (optional)
- Maintain backward compatibility with existing endpoints
- Cache API responses to reduce rate limiting issues

## Manager Policy Compliance

### Manager P3: Gemini Trial Tagging

- **OFFLINE ONLY**: Never use Gemini in runtime paths
- **Validation Protocol**: Batch tag 200 ovarian trials → human spot-review 30 diverse trials
- **Acceptance Criteria**: ≥90% tag accuracy; otherwise adjust prompt taxonomy and re-tag
- **Metadata Persistence**: Store `model`, `version`, `parsed_at`, `reviewed_by`, `source_checksum` with each record
- **Update Cadence**: Weekly diff for new/changed trials
- **Uncertain Tags**: Default to neutral vector; never force a mechanism label

### Manager P4: Mechanism Fit Ranking

- **Formula**: combined_score = (0.7 × eligibility) + (0.3 × mechanism_fit)
- **Thresholds**: 
  - Minimum eligibility to enter top-10: ≥0.60
  - Minimum mechanism_fit for mechanism-gated display: ≥0.50
- **Low Mechanism Fit**: If mechanism_fit <0.50, show trial but without mechanism boost + add "low mechanism fit" warning
- **SOC Card**: Never suppress SOC; SOC card remains first-class
- **UI Control**: Provide "Show all trials" toggle for clinician control

### Manager C7: SAE-Aligned Trial Ranking

- **Mechanism Vector**: Patient `sae_mechanism_vector` = [DDR, MAPK, PI3K, VEGF, IO, Efflux] (6D) ⚠️ **NOTE**: Plan uses 7D [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux] - NEED CLARIFICATION
- **Trial MoA Vector**: From offline MoA tagging; store versioned; neutral if unknown
- **Fallback**: If patient vector all zeros/unknown, disable mechanism_fit (β=0) and explain "awaiting NGS; eligibility-only ranking shown"
- **Explanation**: Show breakdown ("DDR 0.82 × PARP+ATR → 0.95 fit")
- **Wrong MoA Handling**: Human review gate; uncertain trials remain neutral

## Critical Manager Policy Gaps Requiring Clarification

### Gap M1: Mechanism Vector Dimension Mismatch

**Issue**: Manager C7 specifies 6D vector `[DDR, MAPK, PI3K, VEGF, IO, Efflux]` but plan uses 7D `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`.

**Questions**:

1. Should we use 6D (per Manager C7) or 7D (per current plan)?
2. Is HER2 pathway excluded from mechanism fit ranking?
3. Should `MechanismFitRanker` accept 6D or 7D vectors?
4. What about trials with HER2-targeted therapies - how do we rank them?

**Impact**: High - affects all mechanism vector conversion and ranking logic

### Gap M2: Gemini Tagging Runtime vs Offline

**Issue**: Phase 5 says "If Gemini tagging exists, use it; otherwise compute from interventions" - but Manager P3 says Gemini is OFFLINE ONLY.

**Questions**:

1. Is runtime keyword matching acceptable as fallback for trials without Gemini tags?
2. Or should we ONLY use pre-tagged Gemini vectors (no runtime computation)?
3. What's the fallback strategy for new trials not yet tagged?

**Impact**: Medium - affects Phase 5 implementation approach

### Gap M3: Missing Manager P3 Validation Protocol

**Issue**: Plan doesn't specify the validation protocol required by Manager P3.

**Questions**:

1. Should validation protocol be part of this plan or separate?
2. Who performs the human spot-review (30 trials)?
3. Where do we store the metadata (model, version, parsed_at, reviewed_by, source_checksum)?

**Impact**: Medium - affects Phase 5 implementation completeness

### Gap M4: Missing Manager P4 UI Requirements

**Issue**: Plan doesn't mention Manager P4 UI requirements.

**Questions**:

1. Are UI requirements out of scope for this backend plan?
2. Or should we add API fields to support these UI requirements (e.g., `low_mechanism_fit_warning`, `soc_card_priority`)?
3. Should we add a `show_all_trials` parameter to the API?

**Impact**: Low-Medium - affects API response schema

### Gap M5: Missing Manager C7 Fallback Logic

**Issue**: Plan mentions mechanism fit ranking but doesn't fully implement Manager C7 fallback.

**Questions**:

1. Should fallback logic be in Phase 4.5 or separate?
2. What's the exact format for breakdown explanation?
3. Should breakdown be per-trial or per-pathway?

**Impact**: Medium - affects Phase 4.5 implementation

## Integration Notes

### Efficacy Prediction Integration

- Accept efficacy_predictions (S/P/E framework output) in request
- Extract top-ranked drugs (e.g., olaparib, niraparib)
- Look up MoA from DRUG_MECHANISM_DB (api/services/client_dossier/dossier_generator.py)
- Map MoA to intervention keywords (e.g., "PARP inhibitor", "checkpoint inhibitor")
- Auto-add to interventions list if auto_infer_interventions=True
- Pathway scores location: `provenance["confidence_breakdown"]["pathway_disruption"]`

### Pathway Score to Mechanism Vector Conversion

- Accept pathway_scores (Dict[str, float]) in request
- Map pathway names to mechanism vector indices:
  - 6D (Manager C7): DDR = 0, MAPK = 1, PI3K = 2, VEGF = 3, IO = 4, Efflux = 5
  - 7D (Current Plan): DDR = 0, MAPK = 1, PI3K = 2, VEGF = 3, HER2 = 4, IO = 5, Efflux = 6
  - ⚠️ NEED CLARIFICATION: Which dimension should we use?
- Normalize pathway scores: **DO NOT normalize** - use raw pathway scores directly (cosine similarity is scale-invariant)
- Create mechanism vector: [ddr_score, mapk_score, pi3k_score, vegf_score, io_score, efflux_score] (6D) or [ddr_score, mapk_score, pi3k_score, vegf_score, her2_score, io_score, efflux_score] (7D)
- IO = 1.0 if TMB ≥20 OR MSI-High, else 0.0
- Use for MechanismFitRanker.rank_trials()

### Sporadic Cancer Support

- Already implemented in autonomous_trial_agent.py (germline_status, tumor_context)
- Use for filtering and biomarker boosting
- Pass through to HybridTrialSearchService.search_optimized()

### Existing Integration Points

- MechanismFitRanker: `api/services/mechanism_fit_ranker.py` (existing, use as-is)
- Pathway Aggregation: `api/services/pathway/aggregation.py` (produces pathway_scores)
- Efficacy Orchestrator: `api/services/efficacy_orchestrator/orchestrator.py` (produces pathway_scores)
- DRUG_MECHANISM_DB: `api/services/client_dossier/dossier_generator.py` (drug → MoA mapping)
- Autonomous Trial Agent: `api/services/autonomous_trial_agent.py` (sporadic cancer support)

## Priority Actions Before Implementation

### ⚠️ **BLOCKER**: Manager Policy Clarifications Required

**Must resolve before starting implementation**:

1. **Gap M1**: Mechanism vector dimension (6D vs 7D)

   - Manager C7 says 6D `[DDR, MAPK, PI3K, VEGF, IO, Efflux]`
   - Plan uses 7D `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
   - **Decision needed**: Which dimension to use?

2. **Gap M2**: Gemini tagging runtime fallback

   - Manager P3: Gemini OFFLINE ONLY
   - Plan Phase 5: Runtime keyword matching as fallback
   - **Decision needed**: Is runtime fallback acceptable?

3. **Gap M3**: Manager P3 validation protocol

   - Plan doesn't specify validation workflow
   - **Decision needed**: Include in this plan or separate?

4. **Gap M4**: Manager P4 UI requirements

   - Plan doesn't mention SOC card, "Show all trials" toggle
   - **Decision needed**: Backend API support or frontend-only?

5. **Gap M5**: Manager C7 fallback logic

   - Plan mentions but doesn't fully implement
   - **Decision needed**: Include in Phase 4.5 or separate?

### ✅ **READY TO PROCEED** (After Manager Clarification):

- Phase 1: Enhance Autonomous Agent Query Generation
- Phase 2: Create Direct API Query Builder
- Phase 3: Parameterize Extraction Scripts
- Phase 4: Create Advanced Query Endpoint (partial - mechanism fit needs clarification)
- Phase 6: Testing & Validation

### ⏸️ **PENDING CLARIFICATION**:

- Phase 4.5: Mechanism Fit Ranking (needs dimension clarification)
- Phase 5: Trial Data Extraction (needs Gemini tagging strategy)

---

## Implementation Approach & Strategy

### My Approach to Tackle This

**1. Start with Non-Blocking Phases**

- Begin with Phase 1, 2, 3 (query generation, API builder, parameterized scripts)
- These don't depend on manager clarifications
- Can proceed in parallel while waiting for manager answers

**2. Build Incrementally**

- Phase 1 → Phase 2 → Phase 3 → Phase 4 (core functionality)
- Phase 4.5 and Phase 5 can be added after manager clarifications
- Each phase is independently testable

**3. Use Existing Patterns**

- Follow `ayesha_trials.py` pattern for mechanism fit integration
- Reuse `HybridTrialSearchService` for sporadic cancer support (already works!)
- Leverage `extract_fresh_recruiting_trials.py` pagination logic

**4. Error Handling Strategy**

- Graceful degradation: If mechanism fit fails, fall back to eligibility-only ranking
- If Gemini tags missing, use runtime keyword matching (if manager approves)
- If Neo4j unavailable, use AstraDB-only results
- If efficacy predictions unavailable, skip auto-inference

**5. Testing Strategy**

- Unit tests for each phase independently
- Integration tests for end-to-end flow
- Validation tests for manager policy compliance
- Performance tests for large result sets

### Key Implementation Decisions

**Decision 1: Mechanism Vector Dimension**

- **Temporary Solution**: Support both 6D and 7D, detect dimension automatically
- **Conversion Logic**: If 7D provided, extract first 6 dimensions for Manager C7 compliance
- **Future**: Once manager clarifies, standardize on one dimension

**Decision 2: Gemini Tagging Fallback**

- **Temporary Solution**: Runtime keyword matching as fallback (document as "pending manager approval")
- **Storage**: Store both Gemini tags and runtime fallback with metadata flag
- **Future**: Remove runtime fallback if manager requires OFFLINE ONLY

**Decision 3: Eligibility Score Source**

- **Primary**: Use `_composite_score` from trial intelligence pipeline if available
- **Fallback**: Use `match_score` from soft boost ranking
- **Default**: Conservative 0.7 if neither available

**Decision 4: Mechanism Vector Conversion**

- **Create**: `api/services/pathway_to_mechanism_vector.py` as standalone service
- **Support**: Both 6D and 7D output (configurable)
- **Handle**: Missing pathways default to 0.0
- **IO Calculation**: Separate function for TMB/MSI-based IO score

### Critical Implementation Details

**1. Pathway Name Normalization**

- System uses lowercase with underscores: `"ddr"`, `"ras_mapk"`, `"tp53"`
- Need mapping function: `normalize_pathway_name(pathway: str) -> str`
- Handle variations: `"DNA Repair"` → `"ddr"`, `"RAS/MAPK"` → `"ras_mapk"`

**2. MoA Vector Format Conversion**

- Storage: Dict format `{"ddr": 0.9, "mapk": 0.0, ...}` (human-readable)
- Ranker Input: List format `[0.9, 0.0, ...]` (efficient)
- Conversion: `convert_moa_dict_to_vector()` function needed
- Handle both `"ras_mapk"` and `"mapk"` keys (backward compatibility)

**3. Efficacy Predictions Structure**

- Extract from: `response.provenance["confidence_breakdown"]["pathway_disruption"]`
- Format: `Dict[str, float]` with lowercase pathway names
- Top drugs: Sort by `efficacy_score` descending, take top 5-10
- MoA lookup: Use `DRUG_MECHANISM_DB.get_drug_mechanism(drug_name)`

**4. Sporadic Cancer Filtering**

- Already implemented in `HybridTrialSearchService.search_optimized()`
- Pass `germline_status` and `tumor_context` directly
- Filtering happens in `_requires_germline()` method
- Biomarker boost in `_apply_biomarker_boost()` method

**5. Mechanism Fit Ranking Integration**

- Check for all-zero vector BEFORE calling ranker (Manager C7)
- Convert Dict → List BEFORE calling ranker (fix existing bug)
- Add "low mechanism fit" warning AFTER ranking (Manager P4)
- Store breakdown explanation in trial metadata

### Performance Considerations

**1. Query Execution**

- ClinicalTrials.gov API: Rate limit 2 req/sec (enforce with `asyncio.sleep(0.5)`)
- Pagination: Use `pageToken` from response, handle up to 1000 trials per query
- Caching: Cache API responses for 24 hours (reduce rate limiting)

**2. Semantic Search**

- AstraDB: Batch queries when possible (5-10 queries at once)
- Neo4j: Use connection pooling, timeout after 5 seconds
- Fallback: If Neo4j unavailable, continue with AstraDB-only

**3. Mechanism Fit Ranking**

- Pre-filter trials: Only rank trials with `eligibility_score ≥ 0.60`
- Batch processing: Process 100 trials at a time
- Vector normalization: Cache normalized vectors (L2 norm computation)

**4. Data Extraction**

- Parallel processing: Extract PI info, enrollment criteria, therapy types in parallel
- Caching: Store enriched data in SQLite `scraped_data_json` field
- Incremental updates: Only re-extract if trial data changed

### Error Handling & Edge Cases

**1. Missing Data**

- No pathway scores: Skip mechanism fit ranking, use eligibility-only
- No Gemini tags: Use runtime keyword matching (if approved) or neutral vector
- No eligibility score: Default to 0.7 (conservative)
- No MoA vector: Use neutral vector `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`

**2. API Failures**

- ClinicalTrials.gov API down: Return cached results with warning
- AstraDB unavailable: Fall back to direct API queries only
- Neo4j unavailable: Use AstraDB-only results (already handled)

**3. Invalid Input**

- Empty conditions list: Return error with helpful message
- Invalid mechanism vector dimension: Auto-detect and convert or return error
- Invalid pathway scores format: Validate and return error with expected format

**4. Performance Degradation**

- Large result sets (>1000 trials): Apply pagination, return top 100
- Slow mechanism fit ranking: Add timeout (5 seconds), fall back to eligibility-only
- Slow semantic search: Add timeout (10 seconds), fall back to direct API

### Additional Questions & Considerations

**Q1: Mechanism Vector Dimension (6D vs 7D)**

- **Current State**: Manager C7 says 6D, but existing code uses 7D
- **Question**: Should we support both dimensions temporarily?
- **Recommendation**: Yes - auto-detect dimension, support both until manager clarifies

**Q2: Gemini Tagging Runtime Fallback**

- **Current State**: Manager P3 says OFFLINE ONLY, but plan suggests runtime fallback
- **Question**: Can we implement runtime fallback with a feature flag?
- **Recommendation**: Yes - add `allow_runtime_moa_tagging` flag, default to False

**Q3: Eligibility Score Calculation**

- **Current State**: Multiple sources (pipeline composite, soft boost, default)
- **Question**: Should we always run trial intelligence pipeline for eligibility scoring?
- **Recommendation**: Optional - if pipeline not run, use soft boost or default

**Q4: Mechanism Fit Ranking Performance**

- **Current State**: No performance benchmarks for large trial sets
- **Question**: What's the acceptable latency for mechanism fit ranking?
- **Recommendation**: Target <2 seconds for 100 trials, timeout at 5 seconds

**Q5: Trial MoA Vector Storage**

- **Current State**: Stored in JSON file `trial_moa_vectors.json`
- **Question**: Should we migrate to database for better querying?
- **Recommendation**: Keep JSON for now, add database metadata table later

**Q6: Pathway Score Normalization**

- **Current State**: Plan says "DO NOT normalize" but pathway scores may be in different ranges
- **Question**: Should we validate pathway score ranges before conversion?
- **Recommendation**: Yes - log warnings if scores outside expected range (0.0 to ~0.5)

**Q7: Integration with Trial Intelligence Pipeline**

- **Current State**: Pipeline is optional in Phase 4
- **Question**: Should pipeline be required for mechanism fit ranking?
- **Recommendation**: Optional but recommended - mechanism fit uses `_composite_score` from pipeline

**Q8: SOC Card Handling**

- **Current State**: Manager P4 says "never suppress SOC"
- **Question**: How do we identify SOC trials in the results?
- **Recommendation**: Add `is_soc_trial` flag to trial data, always include in top results

**Q9: "Show All Trials" Toggle**

- **Current State**: Manager P4 requires this UI control
- **Question**: Should backend support this via API parameter?
- **Recommendation**: Yes - add `show_all_trials: bool = False` parameter, bypass mechanism fit filtering when True

**Q10: Breakdown Explanation Format**

- **Current State**: Manager C7 says "DDR 0.82 × PARP+ATR → 0.95 fit"
- **Question**: What's the exact format for per-pathway breakdown?
- **Recommendation**: `{"DDR": 0.82, "MAPK": 0.12, ...}` dict + formatted string for display

### Last Minute Enhancements

**Enhancement 1: Add Pathway Name Normalization Function**

- Create `normalize_pathway_name()` to handle variations
- Support: "DNA Repair" → "ddr", "RAS/MAPK" → "ras_mapk", etc.
- Location: `api/services/pathway_to_mechanism_vector.py`

**Enhancement 2: Add MoA Vector Conversion Utilities**

- Create `convert_moa_dict_to_vector()` for Dict → List conversion
- Create `convert_vector_to_moa_dict()` for List → Dict conversion (for storage)
- Handle both 6D and 7D formats
- Location: `api/services/pathway_to_mechanism_vector.py`

**Enhancement 3: Add Mechanism Vector Validation**

- Validate dimension (6 or 7)
- Validate range (0.0 to 1.0)
- Validate sum (should be reasonable, not all zeros)
- Location: `api/services/pathway_to_mechanism_vector.py`

**Enhancement 4: Add Performance Monitoring**

- Log query execution time
- Log mechanism fit ranking time
- Log API call counts and rate limiting
- Location: Add to `advanced_trial_queries.py` endpoint

**Enhancement 5: Add Caching Strategy**

- Cache ClinicalTrials.gov API responses (24 hours)
- Cache mechanism fit rankings (1 hour, keyed by mechanism vector hash)
- Cache pathway score conversions (indefinite, keyed by pathway scores hash)
- Location: Add caching layer to `advanced_trial_queries.py`

**Enhancement 6: Add Comprehensive Error Messages**

- User-friendly error messages for API failures
- Detailed error messages for invalid input
- Helpful suggestions for missing data
- Location: Error handling in `advanced_trial_queries.py`

**Enhancement 7: Add Provenance Tracking**

- Track which queries were used
- Track which services were called (AstraDB, Neo4j, ClinicalTrials.gov API)
- Track mechanism fit ranking provenance (vector source, dimension, thresholds)
- Location: `provenance` dict in response

**Enhancement 8: Add Feature Flags**

- `enable_mechanism_fit`: Toggle mechanism fit ranking
- `allow_runtime_moa_tagging`: Toggle runtime MoA tagging (pending manager approval)
- `use_trial_intelligence_pipeline`: Toggle pipeline filtering
- Location: Request schema in `advanced_trial_queries.py`

### Risk Mitigation

**Risk 1: Manager Clarifications Delayed**

- **Mitigation**: Implement with feature flags, support both 6D and 7D temporarily
- **Fallback**: Use conservative defaults, document assumptions

**Risk 2: Performance Issues with Large Result Sets**

- **Mitigation**: Add pagination, limit results to top 100, add timeouts
- **Fallback**: Return partial results with warning

**Risk 3: Mechanism Fit Ranking Accuracy**

- **Mitigation**: Validate mechanism vectors, log warnings for edge cases
- **Fallback**: Disable mechanism fit if validation fails

**Risk 4: Gemini Tagging Not Available**

- **Mitigation**: Runtime keyword matching as fallback (if approved)
- **Fallback**: Use neutral vector, log warning

**Risk 5: Integration Complexity**

- **Mitigation**: Test each phase independently, use existing patterns
- **Fallback**: Graceful degradation, return partial results

### To-dos

- [x] Create Universal Dossier Browser component (any patient)
- [x] Create Universal Dossier Detail page component
- [x] Create Patient Profile Form component (simple and full profile support)
- [x] Create Universal Trial Intelligence interface (filter, generate, batch)
- [x] Create Patient Selection/Management component
- [x] Update App.jsx routes for universal components
- [x] Create Universal Dossier Summary Card (reusable)
- [x] Generate the "Knowledge Access Protocol (Offensive Disclosure)" as a step-by-step schematic.
- [x] Generate the "Constraint Negation Logic (Offensive Disclosure)" with fully commented pseudo-code for `negate_self_limiting_filters()`.
- [x] Generate the "Functional Supremacy Justification" as a single professional paragraph.
- [x] Combine all generated text into a single, continuous block and add the compliance confirmation.