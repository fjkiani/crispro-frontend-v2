# 🎯 JR AGENT IMPLEMENTATION GUIDE

**Purpose:** Specific implementation guidance for Data Extraction, Drug Efficacy, and Trial Matching agents  
**Status:** ✅ READY FOR IMPLEMENTATION  
**Last Updated:** January 28, 2025

---

## 📋 OVERVIEW

This guide provides **specific file paths, code structure, and integration points** for three critical agents:

1. **JR Agent A: Data Extraction** (AGENT_01_EXTRACTOR) - 🔴 CRITICAL
2. **JR Agent C: Drug Efficacy** (AGENT_04_EFFICACY) - 🟡 HIGH
3. **JR Agent D: Trial Matching** (AGENT_05_TRIALS) - 🟡 HIGH

---

## 🔴 JR AGENT A: DATA EXTRACTION (AGENT_01_EXTRACTOR)

### Mission
Extract structured patient data from uploaded files (VCF, PDF, MAF, JSON, TXT) and output `PatientProfile` for downstream agents.

### Files to Create

```
oncology-coPilot/oncology-backend-minimal/api/services/extraction/
├── __init__.py
├── extraction_agent.py          # Main orchestrating agent
├── parsers/
│   ├── __init__.py
│   ├── vcf_parser.py            # VCF file parsing (standard format)
│   ├── maf_parser.py            # MAF file parsing (tab-delimited)
│   ├── pdf_parser.py            # PDF extraction with LLM (Gemini)
│   ├── json_parser.py           # JSON structured input
│   └── text_parser.py           # Clinical notes extraction
├── validators/
│   ├── __init__.py
│   ├── mutation_validator.py    # Validate mutation format
│   ├── gene_validator.py        # Validate gene names (HGNC)
│   └── quality_checker.py       # Data quality assessment
├── normalizers/
│   ├── __init__.py
│   ├── hgvs_normalizer.py       # Normalize HGVS notation
│   ├── gene_normalizer.py       # Normalize gene symbols
│   └── disease_normalizer.py    # Normalize disease terms
└── tests/
    ├── test_vcf_parser.py
    ├── test_pdf_parser.py
    ├── test_mutation_validator.py
    └── test_integration.py
```

### Key Implementation Details

#### 1. `extraction_agent.py` - Main Agent

**Location:** `api/services/extraction/extraction_agent.py`

**Key Responsibilities:**
- Orchestrate parsing → validation → normalization → quality checking
- Build `PatientProfile` from extracted data
- Track provenance (source, method, confidence)
- Handle errors gracefully (log warnings, return partial data)

**Integration Points:**
- **Input:** File upload via `/api/agents/extract` endpoint
- **Output:** `PatientProfile` dataclass (defined in `api/services/orchestrator/state.py`)
- **Dependencies:** None (foundation module)

**Reference:** `.cursor/MOAT/orchestration/01_DATA_EXTRACTION_AGENT.mdc` (lines 158-314)

#### 2. `parsers/vcf_parser.py` - VCF Parser

**Location:** `api/services/extraction/parsers/vcf_parser.py`

**Key Features:**
- Parse VCF 4.1, 4.2, 4.3 formats
- Extract mutations from INFO field (GENE, HGVSc, HGVSp)
- Extract VAF from FORMAT fields (AF, AD)
- Handle gzipped VCF files
- Support multi-sample VCFs (use first sample)

**Reference:** `.cursor/MOAT/orchestration/01_DATA_EXTRACTION_AGENT.mdc` (lines 447-556)

#### 3. `parsers/pdf_parser.py` - LLM-Based PDF Extraction

**Location:** `api/services/extraction/parsers/pdf_parser.py`

**Key Features:**
- Use PyMuPDF (`fitz`) for text extraction
- Use Gemini API (`google-genai`) for structured extraction
- Support Foundation Medicine, Tempus, Guardant Health reports
- Extract mutations, clinical data, biomarkers, demographics
- Return JSON with structured data

**Dependencies:**
```python
PyMuPDF>=1.23.0          # PDF text extraction
google-genai>=0.3.0      # LLM for PDF extraction
```

**Reference:** `.cursor/MOAT/orchestration/01_DATA_EXTRACTION_AGENT.mdc` (lines 316-445)

#### 4. `validators/gene_validator.py` - Gene Name Validation

**Location:** `api/services/extraction/validators/gene_validator.py`

**Key Features:**
- Validate gene names against HGNC (Human Gene Nomenclature Committee)
- Normalize aliases to official symbols (e.g., "BRCA-1" → "BRCA1")
- Flag invalid gene names for manual review
- Use existing HGNC API or local database

**Integration:** Called by `extraction_agent.py` after parsing

#### 5. `normalizers/hgvs_normalizer.py` - HGVS Normalization

**Location:** `api/services/extraction/normalizers/hgvs_normalizer.py`

**Key Features:**
- Normalize HGVS notation (coding and protein)
- Validate HGVS format
- Use `hgvs` Python library for parsing
- Handle edge cases (ambiguous notation, legacy formats)

**Dependencies:**
```python
hgvs>=1.5.4              # HGVS notation parsing
```

### API Endpoint

**Endpoint:** `POST /api/agents/extract`

**Location:** `api/routers/agents.py` (add new route)

**Request:**
```python
{
    "file": "binary_file_content",  # multipart/form-data
    "file_type": "vcf",             # vcf, maf, pdf, json, txt
    "metadata": {
        "patient_id": "optional",
        "lab": "Foundation Medicine",
        "report_date": "2025-01-15"
    }
}
```

**Response:**
```python
{
    "patient_profile": PatientProfile,
    "extraction_time_ms": 1234,
    "warnings": ["ca125_missing", "hrd_pending"]
}
```

### Acceptance Criteria

- [ ] Can parse VCF files with >95% accuracy
- [ ] Can parse MAF files with >95% accuracy
- [ ] Can extract mutations from PDF reports with >85% accuracy
- [ ] All gene names normalized to HGNC
- [ ] All HGVS notation validated
- [ ] Data quality flags generated for missing data
- [ ] Provenance tracked for all extractions
- [ ] Processing time <10 seconds for typical files
- [ ] Unit test coverage >80%

### Estimated Time: 2-3 days

---

## 🟡 JR AGENT C: DRUG EFFICACY (AGENT_04_EFFICACY)

### Mission
Rank drugs by efficacy using S/P/E framework (Sequence/Pathway/Evidence) with validated formula: `efficacy = 0.3*S + 0.4*P + 0.3*E + clinvar_prior`

### Files to Create

```
oncology-coPilot/oncology-backend-minimal/api/services/efficacy/
├── __init__.py
├── drug_efficacy_agent.py       # Main orchestrating agent
├── sequence_scorer.py            # S: Evo2 integration
├── pathway_scorer.py             # P: Pathway alignment
├── evidence_scorer.py            # E: Literature + ClinVar
├── drug_ranker.py                # Combine S/P/E, rank drugs
├── tier_classifier.py            # Assign tiers (I/II/III/Research)
├── drug_catalog.py              # Drug metadata (name, class, MoA)
├── constants.py                  # Weights, thresholds, mappings
└── tests/
    ├── test_sequence_scorer.py
    ├── test_pathway_scorer.py
    ├── test_evidence_scorer.py
    ├── test_drug_ranker.py
    └── test_integration.py
```

### Key Implementation Details

#### 1. `drug_efficacy_agent.py` - Main Agent

**Location:** `api/services/efficacy/drug_efficacy_agent.py`

**Key Responsibilities:**
- Orchestrate S/P/E scoring pipeline
- Rank drugs by efficacy score
- Build 7D mechanism vector for trial matching
- Generate rationale for each drug
- Assign evidence tiers

**Integration Points:**
- **Input:** `PatientProfile` (from Module 01), `BiomarkerProfile` (from Module 02), `ResistancePrediction` (from Module 03)
- **Output:** `DrugEfficacyResult` with ranked drugs
- **Dependencies:** Modules 01, 02, 03

**Reference:** `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` (lines 149-354)

#### 2. `sequence_scorer.py` - Evo2 Integration

**Location:** `api/services/efficacy/sequence_scorer.py`

**Key Features:**
- Call existing `/api/evo/score_variant_multi` endpoint
- Filter mutations to drug-relevant genes (target_genes)
- Normalize Evo2 delta scores to 0-1 range
- Handle Evo2 API errors gracefully (fallback to 0.5)
- Track provenance (which mutations scored, delta values)

**Existing Endpoints to Use:**
- `POST /api/evo/score_variant_multi` - Multi-window Evo2 scoring
- `POST /api/evo/score_variant_exon` - Exon-context scoring

**Reference:** `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` (lines 356-449)

#### 3. `pathway_scorer.py` - Pathway Alignment

**Location:** `api/services/efficacy/pathway_scorer.py`

**Key Features:**
- Map drug MoA to pathways (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
- Calculate pathway burden from mutations or resistance predictions
- Score alignment: drug pathways vs patient disrupted pathways
- Use existing pathway mappings from `api/services/pathway/`

**Integration:**
- Use `ResistancePrediction.pathways` if available (from Module 03)
- Otherwise infer from mutations using gene-pathway mappings

**Reference:** `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` (lines 451-512)

#### 4. `evidence_scorer.py` - Literature + ClinVar

**Location:** `api/services/efficacy/evidence_scorer.py`

**Key Features:**
- Call existing `/api/evidence/literature` endpoint
- Search PubMed, OpenAlex, S2 for drug-variant-disease evidence
- Lookup ClinVar for variant pathogenicity
- Score evidence quality (Phase III > Phase II > case reports)
- Extract citations for rationale

**Existing Services to Use:**
- `api/services/evidence/literature_client.py` - Literature search
- `api/services/evidence/clinvar_client.py` - ClinVar lookup

**Reference:** `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` (lines 514-597)

#### 5. `drug_catalog.py` - Drug Metadata

**Location:** `api/services/efficacy/drug_catalog.py`

**Key Features:**
- Define drug catalog with name, class, target_genes, pathways, mechanism
- Map drugs to diseases (indications)
- Track biomarker requirements (e.g., PARP requires HRD+)
- Include FDA approval status

**Example Drugs:**
- PARP inhibitors: Olaparib, Niraparib, Rucaparib
- Platinum: Carboplatin, Cisplatin
- Checkpoint inhibitors: Pembrolizumab, Nivolumab
- Targeted: Vemurafenib (BRAF), Trastuzumab (HER2)
- Anti-VEGF: Bevacizumab

**Reference:** `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` (lines 514-597)

#### 6. `tier_classifier.py` - Evidence Tier Assignment

**Location:** `api/services/efficacy/tier_classifier.py`

**Key Features:**
- Assign tiers based on efficacy score, evidence quality, FDA approval
- **Tier I:** On-label, guideline-backed, Phase III evidence
- **Tier II:** Strong evidence, off-label, Phase II/III
- **Tier III:** Mechanistic support, limited evidence
- **Research:** Preclinical or early-phase only

**Reference:** `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` (lines 267-272)

### Existing Code to Reuse

**✅ Already Built:**
- `api/services/efficacy_orchestrator/` - Existing S/P/E implementation
- `api/routers/efficacy.py` - Existing efficacy endpoint
- `api/services/sequence_scorers/evo2_scorer.py` - Evo2 integration
- `api/services/pathway/` - Pathway aggregation
- `api/services/evidence/` - Literature and ClinVar clients

**⚠️ Refactor Needed:**
- Extract reusable components from `efficacy_orchestrator/`
- Create clean agent interface matching MOAT structure
- Separate concerns (scoring vs ranking vs tier classification)

### API Endpoint

**Endpoint:** `POST /api/agents/drugs`

**Location:** `api/routers/agents.py` (add new route)

**Request:**
```python
{
    "patient_id": "string",
    "mutations": Mutation[],
    "disease": "ovarian_cancer",
    "biomarker_profile": BiomarkerProfile,  # optional
    "resistance_prediction": ResistancePrediction,  # optional
    "drug_classes": ["PARP_inhibitor", "platinum"]  # optional filter
}
```

**Response:**
```python
{
    "drug_efficacy_result": DrugEfficacyResult,
    "calculation_time_ms": 2345
}
```

### Acceptance Criteria

- [ ] S/P/E formula correctly implemented (0.3*S + 0.4*P + 0.3*E + clinvar_prior)
- [ ] Evo2 integration for sequence scoring (S component)
- [ ] Pathway alignment calculated (P component)
- [ ] Evidence from literature integrated (E component)
- [ ] Drugs ranked by efficacy score
- [ ] Tiers assigned (I/II/III/Research)
- [ ] 7D mechanism vector built for trial matching
- [ ] Rationale generated for each drug
- [ ] Unit test coverage >80%

### Estimated Time: 3-4 days

---

## 🟡 JR AGENT D: TRIAL MATCHING (AGENT_05_TRIALS)

### Mission
Match patients to clinical trials using mechanism-based ranking with 7D vector cosine similarity.

### Files to Enhance/Create

```
oncology-coPilot/oncology-backend-minimal/api/services/trials/
├── __init__.py
├── trial_matching_agent.py      # Main agent (NEW - wire existing services)
├── query_generator.py            # Generate search queries (NEW)
├── clinicaltrials_client.py     # ClinicalTrials.gov API (NEW)
├── mechanism_ranker.py          # 7D vector ranking (ENHANCE existing)
├── eligibility_checker.py       # Check eligibility (ENHANCE existing)
├── trial_parser.py              # Parse trial data (NEW)
├── moa_extractor.py             # Extract trial MoA vector (NEW)
├── constants.py                 # Query templates, weights (NEW)
└── tests/
    ├── test_query_gen.py
    ├── test_matching.py
    └── test_eligibility.py
```

### Existing Code to Wire

**✅ Already Built:**
- `api/services/autonomous_trial_agent.py` - Trial search service
- `api/services/mechanism_fit_ranker.py` - 7D vector ranking
- `api/services/ayesha_trial_matching/` - Eligibility filters, scoring engine
- `api/services/trial_intelligence/` - Trial intelligence services

**⚠️ Task: Wire to Orchestrator**
- Create `trial_matching_agent.py` that wraps existing services
- Match MOAT interface (input: `PatientProfile`, output: `TrialMatchingResult`)
- Integrate with orchestrator message bus

### Key Implementation Details

#### 1. `trial_matching_agent.py` - Main Agent (NEW)

**Location:** `api/services/trials/trial_matching_agent.py`

**Key Responsibilities:**
- Generate disease/mutation-specific queries
- Search ClinicalTrials.gov via existing services
- Extract MoA vectors for trials
- Rank by mechanism fit (cosine similarity)
- Check eligibility criteria
- Return ranked matches with rationale

**Integration Points:**
- **Input:** `PatientProfile`, `BiomarkerProfile`, `mechanism_vector` (7D from Module 04)
- **Output:** `TrialMatchingResult` with ranked `TrialMatch[]`
- **Dependencies:** Modules 01, 02, 04

**Wire Existing Services:**
```python
from api.services.autonomous_trial_agent import AutonomousTrialAgent
from api.services.mechanism_fit_ranker import MechanismRanker
from api.services.ayesha_trial_matching.eligibility_filters import EligibilityFilters
```

**Reference:** `.cursor/MOAT/orchestration/05_TRIAL_MATCHING_AGENT.mdc` (lines 177-419)

#### 2. `query_generator.py` - Query Generation (NEW)

**Location:** `api/services/trials/query_generator.py`

**Key Features:**
- Generate 5-10 disease/mutation-specific queries
- Query types: disease, gene, basket, DDR, PARP, IO, TMB-H, MSI-H, HRD
- Context-aware query selection (detect DDR mutations, IO eligibility)
- Priority-based ordering (disease first, then gene-specific, then biomarkers)

**Query Templates:**
```python
QUERY_TEMPLATES = {
    'disease': "{disease} cancer clinical trial",
    'gene': "{gene} mutation {disease} trial",
    'basket': "{disease} basket trial tumor agnostic",
    'ddr': "{disease} DNA repair deficiency trial",
    'parp': "{disease} PARP inhibitor trial",
    'io': "{disease} immunotherapy checkpoint inhibitor trial",
    'tmb_high': "TMB-high tumor agnostic trial",
    'msi_high': "MSI-high dMMR clinical trial",
    'hrd': "{disease} HRD homologous recombination deficiency trial"
}
```

**Reference:** `.cursor/MOAT/orchestration/05_TRIAL_MATCHING_AGENT.mdc` (lines 421-549)

#### 3. `clinicaltrials_client.py` - ClinicalTrials.gov API (NEW)

**Location:** `api/services/trials/clinicaltrials_client.py`

**Key Features:**
- Implement ClinicalTrials.gov v2 API client
- Handle pagination (fetch all results)
- Caching (24h TTL for trial data)
- Rate limiting (respect API limits)
- Filter by status (RECRUITING, NOT_YET_RECRUITING)

**API Endpoint:**
```
https://clinicaltrials.gov/api/v2/studies
```

**Or Reuse Existing:**
- Check `api/services/clinical_trial_search_service.py`
- Check `api/services/hybrid_trial_search.py`
- May already have ClinicalTrials.gov integration

#### 4. `mechanism_ranker.py` - 7D Vector Ranking (ENHANCE)

**Location:** `api/services/trials/mechanism_ranker.py`

**Existing:** `api/services/mechanism_fit_ranker.py`

**Enhancement:**
- Wrap existing `MechanismRanker` class
- Add combined scoring: `0.7*eligibility + 0.3*mechanism_fit`
- Add threshold filtering (min_eligibility=0.5, min_mechanism_fit=0.3)
- Return sorted trials with scores

**Reference:** `.cursor/MOAT/orchestration/05_TRIAL_MATCHING_AGENT.mdc` (lines 551-630)

#### 5. `moa_extractor.py` - Extract Trial MoA Vector (NEW)

**Location:** `api/services/trials/moa_extractor.py`

**Key Features:**
- Extract MoA from trial descriptions/interventions
- Map interventions to 7D vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
- Handle unknown drugs (default to 0.5 for all pathways)
- Use drug catalog from Module 04 to map drugs to pathways

**Example:**
```python
# Trial intervention: "Olaparib"
# Extract: PARP inhibitor → DDR pathway → vector[0] = 1.0

# Trial intervention: "Pembrolizumab"
# Extract: Checkpoint inhibitor → IO pathway → vector[5] = 1.0
```

#### 6. `eligibility_checker.py` - Eligibility Checking (ENHANCE)

**Location:** `api/services/trials/eligibility_checker.py`

**Existing:** `api/services/ayesha_trial_matching/eligibility_filters.py`

**Enhancement:**
- Parse eligibility criteria from trial data
- Match against patient data (mutations, biomarkers, demographics)
- Return `EligibilityCriteria` with matched/unmatched/uncertain lists
- Handle uncertain criteria (data not available)

**Reference:** `.cursor/MOAT/orchestration/05_TRIAL_MATCHING_AGENT.mdc` (lines 661-672)

### API Endpoint

**Endpoint:** `POST /api/agents/trials`

**Location:** `api/routers/agents.py` (add new route)

**Request:**
```python
{
    "patient_id": "string",
    "mutations": Mutation[],
    "disease": "ovarian_cancer",
    "mechanism_vector": [0.8, 0.2, 0.1, 0.0, 0.0, 0.5, 0.0],  # 7D from Module 04
    "biomarker_profile": BiomarkerProfile,  # optional
    "max_results": 10
}
```

**Response:**
```python
{
    "trial_matching_result": TrialMatchingResult,
    "search_time_ms": 3456
}
```

### Acceptance Criteria

- [ ] Generate 5-10 disease/mutation-specific queries
- [ ] Search ClinicalTrials.gov successfully
- [ ] Extract MoA vectors for trials
- [ ] Rank by mechanism fit (cosine similarity)
- [ ] Check eligibility criteria
- [ ] Return top 10 matches with rationale
- [ ] Include per-pathway alignment
- [ ] Search time <10 seconds
- [ ] Unit test coverage >80%

### Estimated Time: 3 days (mostly wiring existing services)

---

## 🔗 INTEGRATION WITH ORCHESTRATOR

### Message Bus Integration

All agents should publish results to the orchestrator message bus:

```python
from api.services.orchestrator.message_bus import MessageBus

# After processing
message_bus.publish(
    event_type="patient_profile_extracted",
    payload=patient_profile,
    source="data_extraction_agent"
)
```

### State Store Integration

All agents should read/write to state store:

```python
from api.services.orchestrator.state_store import StateStore

# Read patient state
state = state_store.get_patient_state(patient_id)

# Write results
state_store.update_patient_state(
    patient_id=patient_id,
    updates={"drug_efficacy_result": result}
)
```

### API Router Registration

Register all agent endpoints in `api/routers/agents.py`:

```python
from fastapi import APIRouter
from api.services.extraction.extraction_agent import DataExtractionAgent
from api.services.efficacy.drug_efficacy_agent import DrugEfficacyAgent
from api.services.trials.trial_matching_agent import TrialMatchingAgent

router = APIRouter(prefix="/api/agents", tags=["agents"])

@router.post("/extract")
async def extract_data(...):
    agent = DataExtractionAgent()
    return await agent.extract(...)

@router.post("/drugs")
async def rank_drugs(...):
    agent = DrugEfficacyAgent()
    return await agent.rank_drugs(...)

@router.post("/trials")
async def match_trials(...):
    agent = TrialMatchingAgent()
    return await agent.match(...)
```

---

## 📚 REFERENCE DOCUMENTS

| Document | Purpose |
|----------|---------|
| `.cursor/MOAT/orchestration/01_DATA_EXTRACTION_AGENT.mdc` | Complete Data Extraction spec |
| `.cursor/MOAT/orchestration/04_DRUG_EFFICACY_AGENT.mdc` | Complete Drug Efficacy spec |
| `.cursor/MOAT/orchestration/05_TRIAL_MATCHING_AGENT.mdc` | Complete Trial Matching spec |
| `.cursor/MOAT/orchestration/00_MASTER_INDEX.mdc` | Navigation hub, dependencies |
| `.cursor/MOAT/orchestration/10_STATE_MANAGEMENT.mdc` | Orchestrator integration |
| `.cursor/MOAT/orchestration/11_API_CONTRACTS.mdc` | API schema definitions |

---

## ✅ IMPLEMENTATION CHECKLIST

### JR Agent A: Data Extraction
- [ ] Create `api/services/extraction/` directory structure
- [ ] Implement `extraction_agent.py` with orchestration logic
- [ ] Implement `vcf_parser.py` with VCF 4.x support
- [ ] Implement `pdf_parser.py` with LLM extraction
- [ ] Implement `gene_validator.py` with HGNC validation
- [ ] Implement `hgvs_normalizer.py` with validation
- [ ] Create API endpoint `/api/agents/extract`
- [ ] Write unit tests (>80% coverage)
- [ ] Integration test with sample files

### JR Agent C: Drug Efficacy
- [ ] Create `api/services/efficacy/` directory structure
- [ ] Implement `drug_efficacy_agent.py` with S/P/E orchestration
- [ ] Implement `sequence_scorer.py` with Evo2 integration
- [ ] Implement `pathway_scorer.py` with pathway alignment
- [ ] Implement `evidence_scorer.py` with literature/ClinVar
- [ ] Implement `drug_catalog.py` with drug metadata
- [ ] Implement `tier_classifier.py` with tier assignment
- [ ] Create API endpoint `/api/agents/drugs`
- [ ] Write unit tests (>80% coverage)
- [ ] Integration test with sample patient data

### JR Agent D: Trial Matching
- [ ] Create `api/services/trials/` directory structure
- [ ] Implement `trial_matching_agent.py` (wire existing services)
- [ ] Implement `query_generator.py` with query templates
- [ ] Implement `clinicaltrials_client.py` (or reuse existing)
- [ ] Enhance `mechanism_ranker.py` (wrap existing)
- [ ] Implement `moa_extractor.py` for trial MoA extraction
- [ ] Enhance `eligibility_checker.py` (wrap existing)
- [ ] Create API endpoint `/api/agents/trials`
- [ ] Write unit tests (>80% coverage)
- [ ] Integration test with sample patient data

---

## 🎯 SUCCESS METRICS

| Agent | Metric | Target |
|-------|--------|--------|
| **Data Extraction** | VCF parsing accuracy | >95% |
| **Data Extraction** | PDF extraction accuracy | >85% |
| **Data Extraction** | Processing time | <10 seconds |
| **Drug Efficacy** | S/P/E formula correctness | 100% |
| **Drug Efficacy** | Evo2 integration | 100% usage |
| **Drug Efficacy** | Tier classification accuracy | >90% |
| **Trial Matching** | Query generation | 5-10 queries |
| **Trial Matching** | Search time | <10 seconds |
| **Trial Matching** | Mechanism fit accuracy | >80% |

---

**Document Status:** ✅ COMPLETE  
**Next Steps:** Agents can start implementation using this guide  
**Questions?** Review the referenced MDC files for detailed specifications

