# 📊 ITERATION 7: RESEARCH & DESIGN SYSTEMS

**Status**: ✅ **COMPLETE**  
**Duration**: 2-3 hours  
**Created**: January 14, 2025  
**Cycle**: I7 (Cycle 7)

---

## **7.1 EXECUTIVE SUMMARY**

### **7.1.1 Research & Design Systems Overview**

The platform includes comprehensive research and design capabilities beyond clinical decision support:

1. **VUS Explorer**: Turn "unknown" variants into "understood and actionable" insights
2. **Metastasis Interception**: Design CRISPR interception weapons for metastatic cascade
3. **Hypothesis Validator**: Universal hypothesis testing (food/supplements, compounds)
4. **CRISPR Design**: Guide RNA design with structural validation
5. **Evidence Systems**: RAG, extraction, background jobs
6. **Cohort Intelligence**: Extract, label, and benchmark datasets

---

## **7.2 RESEARCH & DESIGN SYSTEMS**

#### **7.2.1 Metastasis Interception** (`routers/metastasis_interception.py`):

**Purpose**: Design CRISPR interception weapons for metastatic cascade steps

**Endpoint**: `POST /api/metastasis/intercept`

**Workflow**:
1. **Target Lock**: Identify vulnerable gene for mission step
2. **Design**: Generate guide RNA candidates
3. **Safety**: Preview off-target risks (heuristic)
4. **Score**: Rank candidates by assassin_score

**Service**: `services/metastasis_interception_service.py`
- **Target Lock Scoring**: Functionality (0.35) + Essentiality (0.35) + Chromatin (0.15) + Regulatory (0.15)
- **Guide Generation**: Evo2-powered guide RNA design
- **Off-Target Safety**: Hierarchical alignment (minimap2/BLAST)
- **Assassin Score**: Composite ranking (efficacy + safety + mission fit)

**Health Check**: `GET /api/metastasis/intercept/health`
- Returns ruleset version, mission steps configured, gene sets

---

#### **7.2.2 CRISPR Design Router** (`routers/design.py`):

**Purpose**: Generative design endpoints and spacer efficacy prediction

**Endpoints**:
- `POST /api/design/predict_crispr_spacer_efficacy`: Predict on-target efficacy using Evo2 delta scoring

**Method** (`predict_crispr_spacer_efficacy`):
1. **Context Determination**:
   - If `target_sequence` provided: use directly (should be guide + ±50bp flanks = 120bp)
   - Else if `chrom/pos/ref/alt` provided: fetch from Ensembl (±window_size bp flanks)
   - Else: fallback to guide-only context (low confidence)
2. **Evo2 Scoring**: Call `/api/evo/score` to get delta log-likelihood
3. **Sigmoid Transformation**: `efficacy = 1 / (1 + exp(delta / scale_factor))`
   - `scale_factor = 10.0` (as specified in doctrine)
   - More negative delta = more disruptive = higher efficacy
4. **Fallback**: If Evo2 unavailable, use GC-based heuristic
   - GC content penalty: `0.75 - abs(gc - 0.5)`
   - Homopolymer penalty: `-0.1` if AAAA/TTTT/CCCC/GGGG present

**Response**:
- `efficacy_score`: [0,1] efficacy score
- `evo2_delta`: Delta log-likelihood (if available)
- `confidence`: 0.75 (with context) or 0.50 (guide-only) or 0.30 (heuristic)
- `rationale`: Explanation of scoring method
- `provenance`: Method, model_id, context_length, scale_factor

**Feature Flag**: Requires `enable_design_api = True` in config

---

#### **7.2.3 Datasets Router** (`routers/datasets.py`):

**Purpose**: Cohort intelligence - Extract, label, and benchmark datasets from cBioPortal and GDC

**Endpoints**:
- `GET /api/datasets/cbio/studies`: List available cBioPortal studies
- `GET /api/datasets/cbio/mutations`: Extract mutations from cBioPortal study
- `GET /api/datasets/gdc/projects`: List GDC projects
- `GET /api/datasets/gdc/mutations`: Extract mutations from GDC project

**cBioPortal Integration**:
- **Base URL**: `https://www.cbioportal.org/api`
- **Authentication**: Optional Bearer token (`CBIO_TOKEN` env var)
- **Fallback**: Uses `pyBioPortal` library if API unavailable
- **Data Extraction**:
  - Molecular profiles (mutations)
  - Sample lists
  - Mutation data (DETAILED projection)

**GDC Integration**:
- **Base URL**: `https://api.gdc.cancer.gov`
- **Data Extraction**:
  - Projects
  - Cases
  - Mutations (MAF files)

**Use Cases**:
- Extract cohort mutations for benchmarking
- Label datasets with disease annotations
- Generate mutation frequency distributions

---

#### **7.2.4 VUS Explorer** (`components/vus/`, `routers/insights.py`):

**Purpose**: Transform "unknown" variants into "understood and actionable" insights

**Mission**: Fastest way to turn VUS into actionable insights with clear signals, simple story, and one-click next steps.

**Workflow**:
1. **Inputs and Normalization**: HGVS, VCF, coordinates (GRCh38)
2. **Priors Check**: ClinVar lookup, AlphaMissense coverage
3. **Triumvirate Protocol Gate**: Truncation/frameshift detection
4. **S/P/E Core Signals**: Sequence, Pathway, Evidence
5. **Insights Bundle**: Functionality, Chromatin, Essentiality, Regulatory
6. **Fusion Gating**: AlphaMissense integration
7. **SAE Features**: Interpretable features from Evo2
8. **VUS Triage Scoring**: Final verdict with confidence

**Endpoints**:
- `/api/insights/predict_functionality`: Functionality prediction
- `/api/insights/predict_chromatin`: Chromatin accessibility
- `/api/insights/predict_essentiality`: Gene dependency
- `/api/insights/predict_regulatory`: Regulatory impact
- `/api/efficacy/predict`: WIWFM integration

**Frontend**: `oncology-coPilot/oncology-frontend/src/components/vus/`
- MutationTable, AnalysisResults, VEP table
- InsightChips, CoverageChips, ProvenanceBar
- WIWFM integration, CRISPR readiness panel

**Status**: ✅ **FUNCTIONAL** - Needs live API integration

---

#### **7.2.5 Hypothesis Validator** (`routers/hypothesis_validator.py`):

**Purpose**: Universal hypothesis testing for food/supplements and compounds

**Endpoints**:
- `POST /api/hypothesis/validate_food_ab`: A→B dependency validator (no tumor NGS required)
- `POST /api/hypothesis/validate_food_dynamic`: Dynamic compound extraction + S/P/E scoring

**Capabilities**:
1. **Dynamic Target Extraction**:
   - ChEMBL API (2.1M compounds with target data)
   - PubChem API (110M compounds)
   - LLM literature extraction (backup)

2. **A→B Dependency Mapping**:
   - Disease biology → likely A alterations (TP53, HRD, inflammation)
   - A → B dependencies from disease profile
   - Food compound → B target matching

3. **S/P/E Scoring**:
   - Sequence (Evo2) - 30% weight
   - Pathway (weighted aggregation) - 40% weight
   - Evidence (literature/ClinVar) - 30% weight

4. **Treatment Line Intelligence**:
   - Line appropriateness scores
   - Cross-resistance analysis
   - Sequencing fitness

5. **Evidence Mining**:
   - PubMed XML parsing
   - Diffbot full-text extraction
   - LLM paper reading (Gemini/Anthropic/OpenAI)
   - Dosage extraction from papers

**Services**:
- `dynamic_food_extraction.py`: Target extraction (ChEMBL/PubChem/LLM)
- `enhanced_evidence_service.py`: Evidence mining (PubMed/Diffbot)
- `food_spe_integration.py`: S/P/E scoring integration
- `food_treatment_line_service.py`: Treatment line features

**Data Files**:
- `.cursor/ayesha/hypothesis_validator/data/disease_ab_dependencies.json`
- `.cursor/ayesha/hypothesis_validator/data/food_targets.json`
- `.cursor/ayesha/hypothesis_validator/data/cancer_pathways.json`
- `api/resources/universal_disease_pathway_database.json` (50+ diseases)

**Status**: ✅ **PRODUCTION-READY** - Phase 1 (P/E/SAE) complete, Phase 2 (Evo2) experimental

---

#### **7.2.6 Evidence RAG** (`routers/evidence/rag.py`):

**Purpose**: Evidence retrieval augmented generation for literature search

**Endpoints**:
- `POST /api/evidence/rag-query`: Natural language query → evidence-backed answer
- `POST /api/evidence/rag-add-variant`: Add variant to knowledge base
- `GET /api/evidence/rag-stats`: Knowledge base statistics

**Architecture**:
- **Knowledge Base**: Vector embeddings for clinical literature
- **RAG Query Processor**: Natural language → structured queries
- **PubMed Client Enhanced**: E-utils API integration
- **Clinical Insights Processor**: LLM-based answer generation

**LLM Integration**: Google Gemini (requires `GEMINI_API_KEY`)

**Status**: ✅ **OPERATIONAL** - Lazy initialization from `Pubmed-LLM-Agent-main`

---

#### **7.2.7 Evidence Extraction/Jobs** (`routers/evidence/jobs.py`):

**Purpose**: Background job orchestration for evidence processing

**Endpoints**:
- `POST /api/evidence/jobs/crawl`: Crawl article (Diffbot)
- `POST /api/evidence/jobs/summarize`: Summarize article (Google GenAI)
- `POST /api/evidence/jobs/align`: Align article to knowledge base (Google GenAI)
- `GET /api/evidence/jobs/{job_id}`: Get job status

**Job Types**:
- **Crawl**: Full-text extraction using Diffbot
- **Summarize**: LLM-based summarization
- **Align**: Knowledge base alignment

**Job Storage**: In-memory (`JOBS` dict) or Supabase

**Services**:
- `job_service.py`: Background job management and execution

**Status**: ✅ **OPERATIONAL** - Diffbot + Google GenAI integration

---

### **7.3 CLINICAL GENOMICS COMMAND CENTER**

#### **7.3.1 Backend Router** (`routers/clinical_genomics.py`):

**Purpose**: Unified endpoint for comprehensive variant analysis

**Endpoint**: `POST /api/clinical_genomics/analyze_variant`

**Request Schema**:
```python
class AnalyzeVariantRequest:
    mutations: List[Dict[str, Any]]
    disease: Optional[str] = None
    profile: str = "baseline"  # baseline | richer | fusion
    include: List[str] = []  # optional: ["acmg", "pgx"]
    germline_variants: Optional[List[Dict]] = None
    guides: Optional[List[str]] = None
```

**Response Structure**:
- `efficacy`: S/P/E analysis (wraps `/api/efficacy/predict`)
- `toxicity`: None (SLICE 3 - future)
- `off_target`: None (SLICE 3 - future)
- `kg_context`: None (SLICE 3 - future)
- `provenance`: Run ID, methods, model IDs

**Profile Modes**:
- **baseline**: SP mode (Sequence + Pathway), fast path (skip evidence/insights/calibration)
- **richer**: SPE mode (full S/P/E), force exon scan (multi-window)
- **fusion**: SPE mode, force exon scan, includes fusion scoring

**Integration**:
- Calls `efficacy_orchestrator` directly (avoids nested HTTP, reduces latency)
- Includes SAE features when available
- Includes cohort signals when available
- Includes calibration snapshot when available

---

### **7.4 AGENTS SYSTEM**

#### **7.4.1 Agent Router** (`routers/agents.py`):

**Purpose**: API endpoints for agent management, execution, and results

**Endpoints**:
- `POST /api/agents`: Create new agent
- `GET /api/agents`: List all agents for current user
- `GET /api/agents/{agent_id}`: Get agent details
- `PUT /api/agents/{agent_id}`: Update agent
- `DELETE /api/agents/{agent_id}`: Delete agent
- `POST /api/agents/{agent_id}/execute`: Execute agent
- `GET /api/agents/{agent_id}/results`: Get agent execution results

**Agent Types**:
- `pubmed_sentinel`: Monitors PubMed for new papers
- `trial_scout`: Searches for new clinical trials
- (Additional types defined in agent manager)

**Services**:
- `agent_manager.py`: CRUD operations for agents
- `agent_executor.py`: Agent execution logic
- `agent_scheduler.py`: Scheduled agent execution

**Authentication**: Requires user authentication (via `get_optional_user` middleware)

---

### **7.5 KNOWLEDGE BASE (KB) SYSTEM**

#### **7.5.1 KB Router** (`routers/kb/router.py`):

**Purpose**: Modular knowledge base system for storing and retrieving clinical knowledge

**Architecture**: Modular design with separated concerns:
- **Items Router** (`/api/kb/items`): List and retrieve KB items
- **Search Router** (`/api/kb/search`): Keyword and vector search
- **Admin Router** (`/api/kb/admin`): Cache reload, statistics
- **Validation Router** (`/api/kb/validate`): Schema validation
- **Client Router** (`/api/kb/client`): Client-specific endpoints

**Services**:
- `kb_client.py`: KB client interface
- `kb_store.py`: KB storage backend
- `kb_validator.py`: Schema validation

**Utilities**:
- `rate_limiter.py`: In-memory rate limiting
- `client_extractor.py`: Client IP extraction

---

### **7.6 KNOWLEDGE GRAPH (KG) ROUTER**

#### **7.6.1 KG Context Router** (`routers/kg.py`):

**Purpose**: Minimal KG context stub - returns ClinVar + AlphaMissense coverage flags and pathway context

**Endpoint**: `POST /api/kg/context`

**Request**:
```python
class KGCtxRequest:
    mutations: List[Dict[str, Any]]  # [{gene, chrom, pos, ref, alt, consequence, build}]
```

**Response**:
- `coverage`: Dict mapping gene → {clinvar: bool, alphamissense: bool}
- `pathways`: Dict mapping gene → List[pathway_names]
- `provenance`: Method, timestamp

**Current Implementation**:
- **Minimal stub**: Returns placeholder data
- **Coverage flags**: Assumes ClinVar lookup possible, AlphaMissense for missense variants
- **Pathway hints**: Hardcoded for BRAF/KRAS/NRAS (RAS/MAPK), TP53 (TP53)

**Status**: ⚠️ **STUB** - Needs full KG integration

---

### **7.7 CLIENT DOSSIER SYSTEM**

#### **7.7.1 Client Dossier Service** (`services/client_dossier/`):

**Purpose**: Generate trial dossiers for clients

**Components**:
- `dossier_generator.py`: Main dossier generation logic
- `dossier_renderer.py`: Dossier rendering/formatting
- `trial_filter.py`: Trial filtering logic
- `trial_querier.py`: Trial querying
- `trial_scraper.py`: Trial data scraping

**Status**: ⚠️ **NEEDS EXPLORATION** - Directory exists but not yet documented

---

### **7.8 TRIAL INTELLIGENCE PIPELINE**

#### **7.8.1 Trial Intelligence Services** (`services/trial_intelligence/`):

**Purpose**: Multi-stage trial processing pipeline

**Stages**:
1. **Stage 1 - Hard Filters** (`stage1_hard_filters/`): Disease, stage, line, status
2. **Stage 2 - Trial Type** (`stage2_trial_type/`): Trial type classification
3. **Stage 3 - Location** (`stage3_location/`): Location filtering
4. **Stage 4 - Eligibility** (`stage4_eligibility/`): Eligibility criteria matching
5. **Stage 5 - LLM Analysis** (`stage5_llm_analysis/`): LLM-based trial fit analysis
6. **Stage 6 - Dossier** (`stage6_dossier/`): Dossier generation

**Pipeline** (`pipeline.py`): Orchestrates all stages

**Config** (`config.py`): Pipeline configuration

**Status**: ⚠️ **NEEDS EXPLORATION** - Directory exists but not yet documented

---

### **7.9 ADDITIONAL SERVICES** (Quick Reference)

#### **7.9.1 Safety Validator** (`services/safety_validator.py`):

**Purpose**: Validate therapeutic sequences for biological safety

**Safety Checks**:
1. **Viral Content Blocking**: Block HIV, SARS, Ebola, Influenza sequences
2. **GC Extreme Filtering**: Filter extreme GC content
3. **Homopolymer Filtering**: Filter long homopolymers
4. **Known Toxic Sequences**: Block known toxic sequences

**Safety Levels**:
- **SAFE**: Passes all checks
- **WARNING**: Minor concerns
- **BLOCKED**: Fails critical checks

**Status**: ✅ **OPERATIONAL**

---

#### **7.9.2 Guidance Router** (`routers/guidance.py`):

**Purpose**: Clinical gating facade over efficacy orchestrator

**Endpoints**:
- `POST /api/guidance/chemo`: Chemotherapy recommendations
- `POST /api/guidance/radonc`: Radiation oncology recommendations
- `POST /api/guidance/synthetic_lethality`: Synthetic lethality analysis

**Logic**:
- Evidence strength classification
- Tier determination (I/II/III/research)
- Resistance/sensitivity marker detection
- Fusion Engine integration for fused S scores

**Status**: ✅ **OPERATIONAL**

---

#### **7.9.3 Myeloma Digital Twin** (`routers/myeloma.py`):

**Purpose**: Myeloma-specific drug response prediction

**Endpoint**: `POST /api/predict/myeloma_drug_response`

**Features**:
- Evo2 live scoring
- Preflight validation (format, REF-check, duplicate collapse)
- Confidence calculation
- Policy-based interpretation
- Supabase logging
- Dual-model comparison

**Status**: ✅ **OPERATIONAL**

---

## **7.10 SUMMARY**

### **Research Systems**:
1. **VUS Explorer**: Variant interpretation with S/P/E + SAE
2. **Hypothesis Validator**: Universal compound testing (50+ diseases, 110M+ compounds)
3. **Evidence RAG**: Conversational literature search
4. **Evidence Jobs**: Background processing pipeline
5. **Cohort Intelligence**: Dataset extraction and benchmarking

### **Design Systems**:
1. **Metastasis Interception**: CRISPR guide design for metastasis prevention
2. **CRISPR Design Router**: Guide RNA efficacy prediction
3. **Safety Validator**: Biological safety validation
4. **Structural Validation**: AlphaFold 3 integration (pLDDT ≥70)

### **Integration Points**:
- **VUS Explorer** → **WIWFM**: One-click therapy fit
- **Hypothesis Validator** → **S/P/E Framework**: Unified scoring
- **Evidence RAG** → **Knowledge Base**: Vector embeddings
- **CRISPR Design** → **Evo2**: Sequence-level scoring
- **Metastasis** → **Evo2 + AlphaFold 3**: Multi-modal validation

---

**Status**: ✅ **CYCLE 7 COMPLETE** - Research & Design Systems  
**Next**: Cycle 8 (SC-I4) - Execution Plans & Case Studies

---