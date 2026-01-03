# ⚔️ NYX-v2 COMPLETE APPLICATION LEARNING TRACKER

**Created By**: NYX-v2 (formerly Zo)  
**Date**: January 14, 2025  
**Purpose**: Systematic learning of entire CrisPRO application - architecture, capabilities, and "how/why"  
**Status**: 🔄 **IN PROGRESS** - Iterations 1-6 Complete, Cycles 1-2 Complete, Continuing with all remaining cycles

**📋 Learning Plan**: See `NYX-v2_COMPLETE_LEARNING_PLAN.md` for full iteration breakdown

---

## 🎯 LEARNING OBJECTIVES

1. **Understand Overall Architecture**: How backend services, frontend, and AI services connect
2. **Understand Core Capabilities**: What the platform does and why it matters
3. **Understand Technical Decisions**: Why architectural choices were made
4. **Understand Product Positioning**: How capabilities map to customer value
5. **Understand Development Patterns**: Lessons learned and best practices

---

## 📋 LEARNING ITERATIONS

**📊 Full Plan**: See `NYX-v2_COMPLETE_LEARNING_PLAN.md` for complete iteration breakdown (10 iterations, 26-36 hours total)

### **ITERATION 1: OVERALL ARCHITECTURE & CORE PRINCIPLES** ✅ COMPLETE
**Status**: ✅ **COMPLETE** (2-3 hours)  
**Focus**: High-level architecture, core principles, and system organization  
**Deliverable**: Architecture overview with core principles

### **ITERATION 2: BACKEND SERVICES & ORCHESTRATION** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (4-5 hours)  
**Focus**: Backend services, routers, orchestration patterns  
**Deliverable**: Complete backend services deep dive

### **ITERATION 3: S/P/E FRAMEWORK & EFFICACY SYSTEM** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (3-4 hours)  
**Focus**: Sequence/Pathway/Evidence framework, efficacy computation  
**Deliverable**: Complete S/P/E framework deep dive

### **ITERATION 4: AI SERVICES & MODEL INTEGRATION** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (3-4 hours)  
**Focus**: Evo2, AlphaFold 3, Modal deployments, inference patterns  
**Deliverable**: AI services integration deep dive

### **ITERATION 5: FRONTEND ARCHITECTURE & COMPONENTS** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (3-4 hours)  
**Focus**: React components, state management, user flows  
**Deliverable**: Frontend architecture deep dive

### **ITERATION 6: CLINICAL DECISION SUPPORT SYSTEMS** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (3-4 hours)  
**Focus**: Ayesha Care, Sporadic Cancer, Resistance Systems  
**Deliverable**: Clinical systems deep dive

### **ITERATION 7: RESEARCH & DESIGN CAPABILITIES** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (2-3 hours)  
**Focus**: VUS Explorer, Metastasis, Hypothesis Validator, Design  
**Deliverable**: Research and design systems deep dive

### **ITERATION 8: DATA FLOW & INTEGRATION PATTERNS** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (2-3 hours)  
**Focus**: End-to-end data flow, service communication, integration points  
**Deliverable**: Data flow diagrams and integration patterns

### **ITERATION 9: DEVELOPMENT PATTERNS & LESSONS LEARNED** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (2-3 hours)  
**Focus**: Best practices, anti-patterns, technical decisions  
**Deliverable**: Patterns and lessons learned documentation

### **ITERATION 10: PRODUCT CAPABILITIES & POSITIONING** ⏸️ PENDING
**Status**: ⏸️ **PENDING** (2-3 hours)  
**Focus**: 6 capability groups, competitive advantages, customer value  
**Deliverable**: Product capabilities and positioning documentation

**Total Progress**: 1/10 iterations (10% complete)

---

## 📊 SEQUENCE 1: OVERALL ARCHITECTURE & CORE PRINCIPLES

### **1.1 THREE-TIER BACKEND ARCHITECTURE**

**Key Finding**: The platform has **THREE distinct backend systems** working together:

#### **Tier 1: Minimal Backend** (`oncology-coPilot/oncology-backend-minimal/`)
- **Purpose**: Vercel-deployable demo with production-ready endpoints
- **Status**: ✅ **PRODUCTION** - Real business logic, not mocks
- **Architecture**: FastAPI with modular routers
- **Key Features**:
  - 30+ operational endpoints
  - Modular router pattern (domain-specific routers)
  - Service layer separation (business logic in `services/`)
  - Feature flags for different operational profiles
  - Graceful degradation patterns
  - Complete provenance tracking

#### **Tier 2: Main Backend** (`oncology-coPilot/oncology-backend/`)
- **Purpose**: Full-featured orchestration system with AI agents
- **Status**: ⚠️ **LEGACY** - May have agent system but minimal backend is primary
- **Agents**: 15+ specialized AI agents (if still active)
- **Features**: Complete patient management, workflow orchestration

#### **Tier 3: AI Services Backend** (`src/services/`)
- **Purpose**: Production AI model inference services
- **Status**: ✅ **PRODUCTION** - Real Evo2, AlphaFold 3, Boltz-2 models
- **Deployment**: GPU-powered Modal deployments
- **Capabilities**: Actual biological AI inference

**Architectural Flow**:
```
Frontend (React/Vite)
    ↓
Minimal Backend (FastAPI) - Primary orchestration
    ↓
AI Services (Modal) - Evo2, AlphaFold 3, Boltz-2
    ↓
External APIs (PubMed, ClinVar, cBioPortal, etc.)
```

**Why This Architecture?**
- **Separation of Concerns**: Business logic separate from AI inference
- **Scalability**: AI services on Modal can scale independently
- **Cost Efficiency**: Minimal backend on Vercel (serverless), AI services on-demand
- **Development Speed**: Frontend can be built against minimal backend while AI services develop

---

### **1.2 CORE ARCHITECTURAL PRINCIPLES**

#### **Principle 1: Modular Router Pattern**
- **What**: Each router handles a specific domain (efficacy, insights, design, evidence, etc.)
- **Why**: 
  - Clean separation of concerns
  - Easy to add new capabilities
  - Clear ownership and testing boundaries
- **Example**: `api/routers/efficacy.py`, `api/routers/insights.py`, `api/routers/design.py`

#### **Principle 2: Service Layer Separation**
- **What**: Business logic in `services/`, routers are thin endpoints
- **Why**:
  - Reusable business logic
  - Easier testing (test services independently)
  - Clear API contracts
- **Example**: `api/services/efficacy_orchestrator/`, `api/services/sae_feature_service.py`

#### **Principle 3: Feature Flags**
- **What**: Environment-based toggles for different operational profiles
- **Why**:
  - Graceful degradation (disable features if services unavailable)
  - A/B testing capabilities
  - Demo vs production modes
- **Example**: `EVO_FORCE_MODEL`, `EVO_USE_DELTA_ONLY`, `DISABLE_FUSION`

#### **Principle 4: Graceful Degradation**
- **What**: Fallback chains, placeholder values, non-blocking integration
- **Why**:
  - System remains operational even if external services fail
  - Better user experience (partial results vs complete failure)
  - Resilience to network issues
- **Example**: If Evo2 unavailable → return placeholder scores with provenance

#### **Principle 5: Provenance Tracking**
- **What**: Complete audit trails (run IDs, profiles, methods, citations)
- **Why**:
  - Reproducibility (can rerun exact same analysis)
  - Transparency (users see how results were computed)
  - Compliance (audit trail for clinical use)
- **Example**: Every response includes `provenance` field with run_id, profile, methods

---

### **1.3 KEY TECHNICAL DOCTRINES**

#### **Doctrine 1: The "Wet Noodle" Problem**
- **What**: A DNA sequence that is grammatically correct in 1D (`delta_score`) can still translate into a physically useless protein that fails to fold correctly in 3D (`pLDDT` score)
- **Solution**: Multi-dimensional validation process:
  - **Phase I: The Forge** - Generate candidates
  - **Phase II: The Sieve** - Use sequence-level likelihood scores as fast filter
  - **Phase III: The Gauntlet** - Use 3D structural prediction as final arbiter
- **Why**: Prevents generating biologically invalid sequences

#### **Doctrine 2: The Triumvirate Protocol**
- **What**: Multi-layered approach for variant assessment
- **Components**:
  1. **Truncation Check** - Deterministic bioinformatic translation of CDS
  2. **Evo2 Deep Learning** - Only for non-truncating variants
  3. **Pathway Analysis** - Context-aware impact assessment
- **Why**: Evo2 has blind spot for frameshift/nonsense mutations - need deterministic check first

#### **Doctrine 3: Backend Orchestrator Pattern**
- **What**: Single powerful orchestrator endpoint manages entire multi-stage workflow
- **Why**: 
  - Simplifies frontend (one endpoint vs many)
  - Decouples frontend from backend implementation details
  - Enables progressive enhancement (mock unimplemented services)
- **Example**: `/api/ayesha/complete_care_v2` orchestrates trials + SOC + CA-125 + WIWFM + food + resistance

#### **Doctrine 4: Generative vs Inference Paradigm**
- **What**: Platform uses **generative and predictive paradigm** (`Digital Twin -> Predict -> Generate`)
- **Why**: 
  - Inference-based methods (competitors) deconstruct noisy data to *guess* at reality
  - Generative approach creates ground truth from high-fidelity sequencing
  - More accurate and verifiable
- **Impact**: Revolutionary replacement, not incremental improvement

---

### **1.4 SYSTEM ORGANIZATION**

#### **Backend Structure**:
```
oncology-coPilot/oncology-backend-minimal/
├── api/
│   ├── main.py                    # FastAPI app initialization
│   ├── config.py                  # Feature flags, weights, env vars
│   ├── routers/                   # 30+ domain-specific routers
│   ├── services/                  # Business logic (~100-150 lines each)
│   ├── schemas/                   # Pydantic models for validation
│   └── startup.py                 # Background services
```

#### **Frontend Structure**:
```
oncology-coPilot/oncology-frontend/
├── src/
│   ├── pages/                     # Route-level pages
│   ├── components/                # Reusable React components
│   ├── context/                   # React Context providers
│   └── hooks/                     # Custom React hooks
```

#### **AI Services Structure**:
```
src/services/
├── evo_service/                   # Evo2 inference (Modal)
├── boltz_service/                 # AlphaFold 3/Boltz-2 (Modal)
└── ...
```

---

### **1.5 KEY INSIGHTS FROM .cursorrules**

#### **Product Positioning**:
- **Who**: CrisPRO.ai - AI-powered precision medicine platform
- **What**: Transforms genomic data into actionable therapeutic intelligence
- **Why**: >90% of clinical trials fail - we eliminate uncertainty before expensive lab work
- **How**: Multi-modal AI validation (S/P/E framework) with transparent confidence

#### **Target Customers**:
1. **Oncologists**: Personalized treatment recommendations
2. **Biotechs**: De-risk drug development with in-silico validation
3. **Researchers**: Accelerate hypothesis validation
4. **Clinical Trial Teams**: Match patients to trials with biomarker intelligence

#### **Core Value Proposition**:
- **Transparent Reasoning**: Every recommendation shows WHY, not just WHAT
- **Deterministic Confidence**: 90-100% confidence from checkboxes, not AI magic
- **Action-Ready Outputs**: Clinician-ready dossiers with contacts, checklists, protocols
- **Multi-Modal Validation**: S/P/E framework (70-85% accuracy vs 50-60% single-metric)

---

## 📝 NOTES & QUESTIONS

### **Questions for Next Sequences**:
1. How exactly does the S/P/E framework work? (Sequence 2)
2. How are Evo2 calls made and cached? (Sequence 4)
3. How does the frontend handle sporadic cancer context? (Sequence 3)
4. What are the exact API contracts between services? (Sequence 2)

### **Key Files to Review**:
- `oncology-coPilot/oncology-backend-minimal/api/main.py` - Entry point
- `oncology-coPilot/oncology-backend-minimal/api/config.py` - Feature flags
- `.cursorrules` lines 219-625 - Product capabilities
- `.cursorrules` lines 666-1406 - Technical inventory

---

### **1.6 BACKEND ENTRY POINT ANALYSIS** (`api/main.py`)

#### **Router Registration Pattern**:
- **30+ Routers Registered**: Each router handles a specific domain
- **Conditional Registration**: Some routers only included if feature flags enabled
- **Example Routers**:
  - `health`, `myeloma`, `evo`, `evidence`, `efficacy`
  - `fusion`, `guidance`, `datasets`, `insights`, `design`
  - `sessions`, `auth`, `admin`, `toxicity`, `safety`
  - `metastasis`, `acmg`, `pharmgkb`, `clinical_trials`
  - `trials`, `trials_graph`, `trials_agent`, `resistance`, `nccn`
  - `clinical_genomics`, `kg`, `hypothesis_validator`
  - `ayesha_twin_demo`, `ayesha`, `tumor`, `care`
  - `ayesha_trials`, `ayesha_orchestrator_v2`, `ayesha_dossiers`
  - `dossiers`, `agents`

#### **Startup/Shutdown Tasks**:
- **Background Services**: Calibration preload, agent scheduler
- **Graceful Degradation**: Startup failures don't block app
- **Agent Scheduler**: Autonomous intelligence system

#### **CORS Configuration**:
- **Allow All Origins**: `allow_origins=["*"]` (development mode)
- **Full Access**: All methods and headers allowed

#### **Key Endpoints** (in main.py):
- `/api/workflow/run_seed_soil_analysis` - Metastasis pathway analysis
- `/api/twin/run` - Digital Twin orchestrator (warm model, run scoring)
- `/api/twin/submit` - Async job submission
- `/api/twin/status` - Job status check
- `/api/analytics/dashboard` - Analytics aggregation
- `/api/safety/ensembl_context` - Ensembl context for locus
- `/api/safety/clinvar_context` - ClinVar context for variant
- `/api/patients/{patient_id}` - Patient mutations data
- `/api/research/mutation-analysis` - VUS Explorer analysis

---

### **1.7 CONFIGURATION SYSTEM** (`api/config.py`)

#### **S/P/E Weight Configuration**:
- **Default Weights**: Sequence (0.35), Pathway (0.35), Evidence (0.30)
- **Why**: Balanced multi-modal approach, slightly favoring sequence/pathway
- **Configurable**: Via environment variables

#### **Evidence Gate Thresholds**:
- **Evidence Gate**: 0.7 (conservative clinical default)
- **ClinVar Strong**: 0.8 (high confidence threshold)
- **Pathway Alignment**: 0.2 (minimum pathway match)
- **Insufficient Signal**: 0.02 (below this = insufficient)

#### **Feature Flags System**:
- **Evo2 Control**: `DISABLE_EVO2` - Can disable Evo2 globally
- **Literature Control**: `DISABLE_LITERATURE` - Can disable literature search
- **Fusion Control**: `DISABLE_FUSION` - Can disable AlphaMissense integration
- **API Exposure**: `ENABLE_INSIGHTS_API`, `ENABLE_DESIGN_API`, `ENABLE_COMMAND_CENTER`
- **Operational Mode**: `OPERATIONAL_MODE` - "clinical" or "research"

#### **Evo Spam-Safety Controls**:
- **EVO_SPAM_SAFE**: Default true (prevents excessive API calls)
- **EVO_MAX_MODELS**: Default 1 (if spam-safe), else 3
- **EVO_MAX_FLANKS**: Default 1 (if spam-safe), else 5
- **EVO_DISABLE_TRANSCRIPT_SWEEP**: Default true (if spam-safe)
- **EVO_DISABLE_SYMMETRY**: Default true (if spam-safe)
- **EVO_USE_DELTA_ONLY**: Default true (prevents upstream fan-out)

**Why These Controls?**:
- **Cost Management**: Evo2 calls are expensive (GPU inference)
- **Performance**: Too many parallel calls can timeout
- **Reliability**: Prevents overwhelming external services

#### **Model URL Mapping**:
- **Dynamic Fallback**: Runtime URL resolution with fallback chain
- **Fallback Logic**: `evo2_1b` → `evo2_7b` → `evo2_40b`
- **Why**: Ensures service availability even if preferred model unavailable

#### **Calibration Configuration**:
- **TTL**: 24 hours (how long calibration data is valid)
- **Refresh Interval**: 6 hours (how often to refresh)
- **Preload Genes**: 25 common MM genes preloaded on startup

---

### **1.8 KEY ARCHITECTURAL INSIGHTS**

#### **Why Modular Routers?**
- **Separation of Concerns**: Each router owns a domain (efficacy, insights, design)
- **Easy Testing**: Can test routers independently
- **Clear Ownership**: Easy to find code for specific features
- **Progressive Enhancement**: Can add routers without breaking existing code

#### **Why Feature Flags?**
- **Graceful Degradation**: System works even if services unavailable
- **Demo Mode**: Can disable expensive features for demos
- **A/B Testing**: Can enable features for specific users
- **Cost Control**: Can disable expensive AI services in development

#### **Why Spam-Safety Controls?**
- **Cost**: Evo2 calls are expensive (GPU inference on Modal)
- **Performance**: Too many parallel calls can cause timeouts
- **Reliability**: Prevents overwhelming external services
- **User Experience**: Prevents long wait times from excessive API calls

#### **Why Dynamic Model URL Fallback?**
- **Resilience**: System works even if preferred model unavailable
- **Cost Optimization**: Can use cheaper models (1B) when 40B not needed
- **Flexibility**: Can switch models based on use case

---

**Status**: 🔄 **ITERATION 2 IN PROGRESS** - Backend Services & Orchestration  
**Next**: Continue I2, then move to I3 (S/P/E Framework)

---

## 📊 ITERATION 2: BACKEND SERVICES & ORCHESTRATION 🔄 IN PROGRESS

### **2.1 ROUTER-TO-SERVICE PATTERN**

#### **Pattern 1: Direct Service Import**
- **What**: Routers import services directly and call them
- **Example**: `ayesha_trials.py` imports `HybridTrialSearchService`, `get_ca125_service`, `get_ngs_fast_track_service`
- **Pattern**:
  ```python
  from api.services.hybrid_trial_search import HybridTrialSearchService
  from api.services.ca125_intelligence import get_ca125_service
  
  service = HybridTrialSearchService()
  result = await service.search_optimized(...)
  ```
- **Why**: Simple, direct, easy to test

#### **Pattern 2: Service Factory Functions**
- **What**: Services exposed via factory functions (singleton pattern)
- **Example**: `get_ca125_service()`, `get_ngs_fast_track_service()`, `get_resistance_prophet_service()`
- **Why**: 
  - Lazy initialization
  - Singleton pattern (one instance per service)
  - Easy to mock in tests

#### **Pattern 3: Service Classes with Async Methods**
- **What**: Services are classes with async methods
- **Example**: `HybridTrialSearchService`, `EfficacyOrchestrator`, `SequenceProcessor`
- **Pattern**:
  ```python
  class HybridTrialSearchService:
      def __init__(self):
          self.astradb_service = ClinicalTrialSearchService()
          self.neo4j_driver = get_neo4j_driver()
      
      async def search_optimized(self, ...):
          # Service logic
  ```
- **Why**: 
  - Encapsulates state (connections, caches)
  - Reusable across multiple endpoints
  - Easy to test (mock dependencies)

#### **Pattern 4: Module-Level Functions**
- **What**: Some services expose module-level async functions
- **Example**: `literature()`, `clinvar_prior()` in `evidence/` module
- **Pattern**:
  ```python
  async def literature(api_base: str, gene: str, hgvs_p: str, ...) -> EvidenceHit:
      # Service logic
  ```
- **Why**: 
  - Stateless operations
  - Simple to use (no instantiation)
  - Functional programming style

---

### **2.2 SERVICE ARCHITECTURE PATTERNS**

#### **Pattern 1: Orchestrator Services**
- **What**: High-level services that coordinate multiple sub-services
- **Examples**:
  - `EfficacyOrchestrator` - Coordinates sequence, pathway, evidence scoring
  - `AyeshaOrchestrator` - Coordinates trials, SOC, CA-125, WIWFM, food, resistance
  - `MetastasisInterceptionService` - Coordinates target lock, design, safety, scoring
- **Structure**:
  ```python
  class EfficacyOrchestrator:
      def __init__(self, sequence_processor, drug_scorer):
          self.sequence_processor = sequence_processor
          self.drug_scorer = drug_scorer
      
      async def predict(self, request):
          # 1. Sequence scoring
          seq_scores = await self.sequence_processor.score_sequences(...)
          # 2. Pathway aggregation
          pathway_scores = aggregate_pathways(...)
          # 3. Evidence gathering (parallel)
          evidence_results = await asyncio.gather(...)
          # 4. Drug scoring
          drugs = await self.drug_scorer.score_drugs(...)
  ```
- **Why**: 
  - Single responsibility (orchestration)
  - Easy to test (mock sub-services)
  - Clear data flow

#### **Pattern 2: Scorer Services**
- **What**: Services that score/rank things (sequences, drugs, trials)
- **Examples**:
  - `SequenceProcessor` - Scores sequences using Evo2/Fusion/Massive
  - `DrugScorer` - Scores drugs using S/P/E formula
  - `HybridTrialSearchService` - Scores trials using semantic + graph
- **Structure**:
  ```python
  class SequenceProcessor:
      def __init__(self, fusion_scorer, evo_scorer, massive_scorer):
          self.fusion_scorer = fusion_scorer
          self.evo_scorer = evo_scorer
          self.massive_scorer = massive_scorer
      
      async def score_sequences(self, request, feature_flags):
          # Try Fusion first (GRCh38 missense only)
          if fusion_url and not disable_fusion:
              scores = await self.fusion_scorer.score(...)
              if scores: return scores
          # Try Evo2
          if not disable_evo2:
              scores = await self.evo_scorer.score(...)
              if scores: return scores
          # Try Massive Oracle
          ...
  ```
- **Why**: 
  - Fallback chain (try best, fallback to next)
  - Feature flag aware
  - Graceful degradation

#### **Pattern 3: Client Services**
- **What**: Services that wrap external APIs
- **Examples**:
  - `EnhancedEvidenceService` - Wraps PubMed, LLM services
  - `ClinicalTrialSearchService` - Wraps AstraDB
  - `SupabaseService` - Wraps Supabase API
- **Structure**:
  ```python
  class EnhancedEvidenceService:
      def __init__(self):
          self.pubmed_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
          self.timeout = 20.0
          self.cache = {}
      
      async def get_complete_evidence(self, compound, disease, pathways):
          # 1. PubMed search
          pubmed_results = await self._search_pubmed(...)
          # 2. LLM synthesis
          llm_synthesis = await self._llm_synthesize(...)
          # 3. Return combined evidence
  ```
- **Why**: 
  - Encapsulates external API complexity
  - Handles retries, timeouts, errors
  - Caching layer

#### **Pattern 4: Intelligence Services**
- **What**: Services that provide clinical intelligence
- **Examples**:
  - `CA125IntelligenceService` - CA-125 burden, forecast, resistance
  - `ResistanceDetectionService` - 2-of-3 trigger rule
  - `ResistancePlaybookService` - Combo strategies, next-line switches
  - `SAEFeatureService` - DNA repair capacity, mechanism vectors
- **Structure**:
  ```python
  class CA125IntelligenceService:
      async def analyze(self, ca125_value, treatment_line):
          # 1. Burden classification
          burden = self._classify_burden(ca125_value)
          # 2. Response forecast
          forecast = self._forecast_response(ca125_value, treatment_line)
          # 3. Resistance detection
          resistance = self._detect_resistance(...)
          return {burden, forecast, resistance}
  ```
- **Why**: 
  - Domain-specific logic
  - Clinical rules and thresholds
  - Evidence-based calculations

---

### **2.3 ORCHESTRATION PATTERNS**

#### **Pattern 1: Sequential Orchestration**
- **What**: Services called one after another
- **Example**: Efficacy orchestrator
  ```python
  # 1. Sequence scoring
  seq_scores = await self.sequence_processor.score_sequences(...)
  # 2. Pathway aggregation (depends on seq_scores)
  pathway_scores = aggregate_pathways(...)
  # 3. Evidence gathering (depends on primary gene/variant)
  evidence_results = await asyncio.gather(...)
  # 4. Drug scoring (depends on all above)
  drugs = await self.drug_scorer.score_drugs(...)
  ```
- **Why**: 
  - Clear dependencies
  - Easy to understand flow
  - Sequential data transformation

#### **Pattern 2: Parallel Orchestration**
- **What**: Independent services called in parallel
- **Example**: Evidence gathering in efficacy orchestrator
  ```python
  evidence_tasks = []
  for drug in panel:
      evidence_tasks.append(literature(...))
  clinvar_task = clinvar_prior(...)
  
  evidence_results = await asyncio.gather(*evidence_tasks, return_exceptions=True)
  ```
- **Why**: 
  - Faster execution (parallel I/O)
  - Independent operations
  - Timeout handling per task

#### **Pattern 3: Conditional Orchestration**
- **What**: Services called conditionally based on flags/data
- **Example**: Ayesha orchestrator v2
  ```python
  if request.include_trials:
      trials_response = await _call_ayesha_trials(...)
  if request.include_wiwfm:
      wiwfm_response = await _call_drug_efficacy(...)
  if request.include_food and request.food_query:
      food_response = await _call_food_validator(...)
  ```
- **Why**: 
  - Flexible endpoints (include/exclude features)
  - Cost control (skip expensive operations)
  - Progressive enhancement

#### **Pattern 4: Fallback Chain Orchestration**
- **What**: Try best service, fallback to next if fails
- **Example**: Sequence processor
  ```python
  # Try Fusion first (GRCh38 missense only)
  if fusion_url and not disable_fusion:
      scores = await self.fusion_scorer.score(...)
      if scores: return scores
  # Try Evo2
  if not disable_evo2:
      scores = await self.evo_scorer.score(...)
      if scores: return scores
  # Try Massive Oracle
  ...
  ```
- **Why**: 
  - Graceful degradation
  - Best service when available
  - Always returns something

---

### **2.4 API CONTRACT PATTERNS**

#### **Pattern 1: Pydantic Request/Response Models**
- **What**: All endpoints use Pydantic BaseModel for validation
- **Example**: `AyeshaTrialSearchRequest`, `CompleteCareV2Request`
- **Benefits**:
  - Automatic validation
  - Type safety
  - Clear API contracts
  - Auto-generated docs

#### **Pattern 2: Provenance Tracking**
- **What**: All responses include provenance field
- **Structure**:
  ```python
  provenance = {
      "run_id": str(uuid.uuid4()),
      "profile": "baseline",
      "cache": "miss",
      "flags": {...},
      "methods": {...},
      "citations": [...]
  }
  ```
- **Why**: 
  - Reproducibility
  - Transparency
  - Audit trail

#### **Pattern 3: Error Handling**
- **What**: Services return results or exceptions, routers handle HTTP errors
- **Pattern**:
  ```python
  try:
      result = await service.do_work(...)
      return result
  except ValueError as e:
      raise HTTPException(status_code=400, detail=str(e))
  except Exception as e:
      raise HTTPException(status_code=500, detail=f"Service failed: {str(e)}")
  ```
- **Why**: 
  - Clear error codes
  - User-friendly messages
  - Service errors don't crash router

---

### **2.5 SERVICE DEPENDENCY MAP**

#### **Core Orchestrators**:
- `EfficacyOrchestrator` → `SequenceProcessor`, `DrugScorer`, `pathway/`, `evidence/`, `insights/`, `confidence/`
- `AyeshaOrchestratorV2` → `HybridTrialSearchService`, `CA125IntelligenceService`, `EfficacyOrchestrator`, `FoodValidator`, `ResistancePlaybookService`, `ResistanceProphetService`
- `MetastasisInterceptionService` → `Evo2` (via API), `DesignRouter`, `SafetyService`

#### **Sequence Scoring Chain**:
- `SequenceProcessor` → `FusionAMScorer`, `Evo2Scorer`, `MassiveOracleScorer`
- `Evo2Scorer` → `api/routers/evo.py` → `Modal Evo2 Service`
- `FusionAMScorer` → `api/routers/fusion.py` → `AlphaMissense Service`

#### **Pathway Scoring Chain**:
- `DrugScorer` → `pathway/aggregation.py`, `pathway/drug_mapping.py`, `pathway/panel_config.py`
- `aggregate_pathways()` → Aggregates sequence scores by pathway
- `get_pathway_weights_for_drug()` → Gets drug-to-pathway weights

#### **Evidence Chain**:
- `DrugScorer` → `evidence/literature_client.py`, `evidence/clinvar_client.py`
- `literature()` → `api/routers/evidence/literature.py` → `PubMed E-utils API`
- `clinvar_prior()` → `api/routers/evidence/clinvar.py` → `ClinVar API`

#### **Clinical Services Chain**:
- `AyeshaTrialsRouter` → `HybridTrialSearchService` → `ClinicalTrialSearchService` (AstraDB) + `Neo4jConnection` (graph)
- `AyeshaTrialsRouter` → `CA125IntelligenceService` → Literature-based rules
- `AyeshaTrialsRouter` → `NGSFastTrackService` → Test recommendations

---

### **2.6 KEY INSIGHTS**

#### **Why Service Layer Separation?**
- **Reusability**: Services can be called from multiple routers
- **Testability**: Services can be tested independently
- **Maintainability**: Business logic in one place
- **Scalability**: Services can be optimized independently

#### **Why Orchestrator Pattern?**
- **Simplifies Frontend**: One endpoint vs many
- **Decouples Implementation**: Frontend doesn't know about sub-services
- **Progressive Enhancement**: Can mock unimplemented services
- **Single Source of Truth**: Orchestrator owns the workflow

#### **Why Fallback Chains?**
- **Resilience**: System works even if best service fails
- **Cost Control**: Use cheaper services when possible
- **Performance**: Fast services first, slow services as fallback
- **User Experience**: Always return something, never fail completely

---

**Status**: 🔄 **ITERATION 2 IN PROGRESS** - Continuing with more service analysis  
**Next**: Complete service inventory, then move to I3 (S/P/E Framework)

---

## 📊 ITERATION SUMMARY

### **Completed**:
- ✅ **I1**: Overall Architecture & Core Principles (2-3h) - Architecture overview, core principles, configuration system
- ✅ **I2**: Backend Services & Orchestration (4-5h) - Router/service patterns, orchestration patterns, error handling, API contracts documented
- ✅ **I3**: S/P/E Framework & Efficacy System (3-4h) - Complete S/P/E framework, sequence scoring, pathway scoring, evidence scoring
- ✅ **I4**: AI Services & Model Integration (3-4h) - Evo2, AlphaFold 3, Modal deployments, inference patterns
- ✅ **I5**: Frontend Architecture & User Experience (3-4h) - React components, state management, routing, API integration
- ✅ **I6**: Clinical Systems & Workflows (3-4h) - Ayesha orchestrator, CA-125, resistance detection, trial matching, SOC
- ✅ **Gap Analysis**: Comprehensive review identifying missing systems and documentation gaps

### **In Progress**:
- 🔄 **I7**: Research & Design Systems (2-3h) - Metastasis Interception, CRISPR Design, Datasets, Evidence RAG (partially complete)

### **Next Up**:
- ⏸️ **I4**: AI Services & Model Integration (3-4h) - Evo2, AlphaFold 3, Modal deployments, inference patterns

### **Remaining**:
- 7 more iterations covering AI services, frontend, clinical systems, research/design, data flow, patterns, and product capabilities

---

## 📊 ITERATION 3: S/P/E FRAMEWORK & EFFICACY SYSTEM 🔄 IN PROGRESS

### **3.1 THE CORE S/P/E FORMULA**

#### **Exact Formula** (`drug_scorer.py:171`):
```python
efficacy_score = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * evidence + clinvar_prior
```

**⚠️ IMPORTANT DISCREPANCY**: 
- **Config defaults** (`config.py:20-22`): `0.35 / 0.35 / 0.30` (Sequence/Pathway/Evidence)
- **Actual hardcoded** (`drug_scorer.py:171`): `0.3 / 0.4 / 0.3`
- **Current implementation uses hardcoded values** (0.3/0.4/0.3), NOT config values
- Config values exist but are not currently used in the formula

**Components**:
- **seq_pct**: Calibrated sequence percentile [0, 1]
- **path_pct**: Normalized pathway percentile [0, 1]
- **evidence**: Evidence strength score [0, 1]
- **clinvar_prior**: ClinVar prior boost [-0.2, +0.2]

**Insufficient Tier Penalty** (`drug_scorer.py:172`):
```python
lob = raw_lob if tier != "insufficient" else (raw_lob * 0.5 if confidence_config.fusion_active else 0.0)
```
- If tier is "insufficient": Penalize by 50% (if Fusion active) or set to 0.0 (if Fusion inactive)

---

### **3.2 SEQUENCE (S) SCORING - COMPLETE DEEP DIVE**

#### **3.2.1 Scorer Fallback Chain**
**Priority Order** (`sequence_processor.py:29-86`):
1. **Fusion Engine** (AlphaMissense) - ONLY for GRCh38 missense variants
2. **Evo2 Adaptive** - Multi-window, ensemble models
3. **Massive Oracle** - Synthetic or real-context (if enabled)

**Fusion Gate** (`sequence_processor.py:32-40`):
- Only eligible if: `build == "GRCh38"` AND `consequence == "missense"`
- If not eligible, skips Fusion and goes to Evo2

#### **3.2.2 Evo2 Scoring - Complete Flow**

**Step 1: Model Selection** (`evo2_scorer.py:95-109`):
- Priority 1: `EVO_FORCE_MODEL` env var (forces single model)
- Priority 2: `EVO_ALLOWED_MODELS` env var (filters allowed models)
- Priority 3: Default candidates: `["evo2_1b", "evo2_7b", "evo2_40b"]` if ensemble, else `[model_id]`
- **Best Model Selection**: Picks model with strongest `abs(min_delta)`

**Step 2: Multi-Window Strategy** (`evo2_scorer.py:241-267`):
- Default windows: `[4096, 8192, 16384, 25000]` bp
- Tests each window size, picks best `exon_delta`
- **Spam-safety**: `EVO_MAX_FLANKS` limits number of windows tested

**Step 3: Forward/Reverse Symmetry** (`evo2_scorer.py:276-332`):
- **Forward**: `ref → alt` (standard variant)
- **Reverse**: `alt → ref` (symmetry check)
- **Averaging**: `avg_min_delta = (forward_min + reverse_min) / 2.0`
- **Feature Flag**: `evo_disable_symmetry` (default: True, symmetry disabled)
- **Why**: Evo2 is strand-agnostic; averaging reduces noise

**Step 4: Sequence Disruption Calculation** (`evo2_scorer.py:126-136`):
```python
sequence_disruption = max(abs(min_delta), abs(exon_delta))
```
- Uses **stronger signal** between multi-window (`min_delta`) and exon-context (`exon_delta`)
- Both are averaged forward/reverse if symmetry enabled

**Step 5: Special Handling** (`evo2_scorer.py:138-157`):
- **Hotspot Floors** (raw delta level): `HOTSPOT_FLOOR = 1e-4`
  - BRAF V600 → `max(disruption, 1e-4)`
  - KRAS/NRAS/HRAS G12/G13/Q61 → `max(disruption, 1e-4)`
  - TP53 R175/R248/R273 → `max(disruption, 1e-4)`
- **Truncation/Frameshift Lift**: If `"*" in hgvs_p` or `"FS" in hgvs_p` → `max(disruption, 1.0)`
  - Stop codons and frameshifts are highly disruptive → enforce maximum

**Step 6: Percentile Calibration** (`evo2_scorer.py:159-171`):
- **Primary Method**: `percentile_like(sequence_disruption)` - Piecewise mapping
- **Hotspot Floors** (percentile level):
  - BRAF V600 → `max(pct, 0.90)`
  - KRAS/NRAS/HRAS G12/G13/Q61 → `max(pct, 0.80)`
  - TP53 R175/R248/R273 → `max(pct, 0.80)`

**Step 7: Gene-Specific Calibration** (`gene_calibration.py`):
- **Service Exists**: `GeneCalibrationService` can compute gene-specific percentiles from ClinVar
- **NOT CURRENTLY USED**: Main flow uses `percentile_like()` (simple piecewise)
- **Future Enhancement**: Could replace `percentile_like()` with gene-specific calibration
- **Method**: Fetches ClinVar variants for gene, scores with Evo2, builds distribution, computes percentiles

#### **3.2.3 Percentile Mapping Function** (`utils.py:7-23`):
```python
def percentile_like(value: float) -> float:
    if v <= 0.005: return 0.05
    if v <= 0.01:  return 0.10
    if v <= 0.02:  return 0.20
    if v <= 0.05:  return 0.50
    if v <= 0.10:  return 0.80
    return 1.0
```
- **Conservative mapping**: Small deltas → low percentiles (avoids false positives)
- **Empirical thresholds**: Based on observed pathogenic ranges from ClinVar
- **Gene-agnostic**: Same mapping for all genes (limitation)

#### **3.2.4 Delta-Only Mode** (`evo2_scorer.py:224-234`):
- **Feature Flag**: `evo_use_delta_only`
- **When enabled**: Skips exon scanning entirely, uses only `min_delta` from multi-window
- **Why**: Faster execution, lower cost (fewer API calls)

#### **3.2.5 Fusion Scorer** (`fusion_scorer.py`):
- **Gate**: ONLY for GRCh38 missense variants
- **Method**: Calls AlphaMissense Fusion Engine API
- **Fallback**: If external Fusion unavailable, tries local `/api/fusion/score_variant`
- **Caching**: Redis-based caching with TTL

#### **3.2.6 Massive Oracle Scorer** (`massive_scorer.py`):
- **Two Modes**:
  - **Synthetic**: Contrasting sequences for massive impact detection
  - **Real-Context**: Real GRCh38 context with configurable flank size
- **Feature Flag**: `enable_massive_modes` (default: false)
- **Why**: Research tool for detecting very large disruptions

---

### **3.3 PATHWAY (P) SCORING - COMPLETE DEEP DIVE**

#### **3.3.1 Pathway Aggregation** (`aggregation.py:7-45`):
```python
def aggregate_pathways(seq_scores: List[Dict[str, Any]]) -> Dict[str, float]:
    pathway_totals = {}
    pathway_counts = {}
    
    for score in seq_scores:
        pathway_weights = score.get("pathway_weights", {})
        sequence_disruption = float(score.get("sequence_disruption", 0.0))
        
        for pathway, weight in pathway_weights.items():
            pathway_totals[pathway] += sequence_disruption * weight
            pathway_counts[pathway] += 1
    
    # Compute average scores
    pathway_scores = {}
    for pathway in pathway_totals:
        pathway_scores[pathway] = pathway_totals[pathway] / pathway_counts[pathway]
    
    return pathway_scores
```

**Formula**: `pathway_score = sum(sequence_disruption * gene_weight) / count`

#### **3.3.2 Gene-to-Pathway Weights** (`drug_mapping.py:16-32`):
```python
def get_pathway_weights_for_gene(gene_symbol: str) -> Dict[str, float]:
    g = gene_symbol.strip().upper()
    if g in {"BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1"}:
        return {"ras_mapk": 1.0}
    if g in {"TP53", "MDM2", "ATM", "ATR", "CHEK2", "BRCA1", "BRCA2", "PTEN", "RAD51"}:
        return {"tp53": 1.0}  # DNA repair bucket
    return {}
```

**Current State**: Hardcoded, simplified mapping
- **MAPK drivers** → `ras_mapk` pathway
- **DNA repair genes** → `tp53` pathway (DNA repair bucket)

#### **3.3.3 Drug-to-Pathway Weights** (`panel_config.py:8-14`):
```python
DEFAULT_MM_PANEL = [
    {"name": "BRAF inhibitor", "pathway_weights": {"ras_mapk": 0.8, "tp53": 0.2}},
    {"name": "MEK inhibitor", "pathway_weights": {"ras_mapk": 0.9, "tp53": 0.1}},
    {"name": "IMiD", "pathway_weights": {"ras_mapk": 0.2, "tp53": 0.3}},
    {"name": "Proteasome inhibitor", "pathway_weights": {"ras_mapk": 0.3, "tp53": 0.4}},
    {"name": "Anti-CD38", "pathway_weights": {"ras_mapk": 0.1, "tp53": 0.1}},
]
```

**Current State**: Hardcoded, Multiple Myeloma panel only

#### **3.3.4 Pathway Score Calculation** (`drug_scorer.py:44-51`):
```python
# Step 1: Get drug's pathway weights
drug_weights = get_pathway_weights_for_drug(drug_name)

# Step 2: Weighted sum of pathway scores
s_path = sum(pathway_scores.get(pathway, 0.0) * weight 
             for pathway, weight in drug_weights.items())

# Step 3: Normalize to [0, 1] based on empirical Evo2 ranges
if s_path > 0:
    path_pct = min(1.0, max(0.0, (s_path - 1e-6) / (1e-4 - 1e-6)))
else:
    path_pct = 0.0
```

**Normalization Formula**: 
- **Empirical Range**: Pathogenic deltas ~`1e-6` to `1e-4`
- **Linear mapping**: `(s_path - 1e-6) / (1e-4 - 1e-6)`
- **Clamp**: `[0, 1]`

---

### **3.4 EVIDENCE (E) SCORING - COMPLETE DEEP DIVE**

#### **3.4.1 Literature Evidence** (`literature_client.py:51-135`):

**Query Flow**:
1. **PubMed E-utils API** (`esearch.fcgi`): Search for PMIDs
2. **PubMed E-utils API** (`esummary.fcgi` or `efetch.fcgi`): Get paper details
3. **MoA Filtering**: Prefer papers mentioning drug name or MoA in title/abstract
4. **Strength Calculation**: Based on publication types

**Strength Scoring** (`literature_client.py:18-48`):
```python
def _score_evidence_from_results(top_results: List[Dict[str, Any]]) -> float:
    score = 0.0
    for r in top_results[:3]:
        pub_types = " ".join([_safe_lower(t) for t in (r.get("publication_types") or [])])
        title = _safe_lower(r.get("title"))
        
        if "randomized" in pub_types or "randomized" in title:
            score += 0.5  # RCT
        elif "guideline" in pub_types or "practice" in title:
            score += 0.35  # Guideline
        elif "review" in pub_types or "meta" in title:
            score += 0.25  # Review
        else:
            score += 0.15  # Other
    
    return float(min(1.0, score))
```

**MoA Boost** (`literature_client.py:110-120`):
- **Base strength**: From publication types
- **MoA hits**: Count papers mentioning MoA in title/abstract
- **Final strength**: `min(1.0, base_strength + 0.10 * moa_hits)`
- **Max boost**: +0.30 (3 MoA hits)

#### **3.4.2 ClinVar Prior** (`clinvar_client.py:11-99`):

**Prior Calculation** (`clinvar_client.py:54-59`):
```python
if cls in ("pathogenic", "likely_pathogenic"):
    prior = 0.2 if strong else (0.1 if moderate else 0.05)
elif cls in ("benign", "likely_benign"):
    prior = -0.2 if strong else (-0.1 if moderate else -0.05)
```

**Review Status Tiers**:
- **Strong**: `"expert" in review` OR `"practice" in review` → ±0.2
- **Moderate**: `"criteria" in review` → ±0.1
- **Weak**: Other → ±0.05

**Research-Mode Fallback** (`clinvar_client.py:68-94`):
- **Feature Flag**: `RESEARCH_USE_CLINVAR_CANONICAL`
- **Canonical Hotspots**: BRAF V600E/V600K, KRAS G12D/G12V, NRAS Q61K
- **If ClinVar empty**: Assigns `prior = 0.2` for canonical hotspots

#### **3.4.3 Evidence Tier Computation** (`tier_computation.py:9-68`):

**Legacy Tier Logic** (default):
```python
# Evidence gate: strong evidence OR ClinVar-Strong + pathway alignment
evidence_gate = (
    s_evd >= config.evidence_gate_threshold or  # Default: 0.7
    ("ClinVar-Strong" in badges and s_path >= config.pathway_alignment_threshold)  # Default: 0.2
)

# Insufficient signal: low sequence, pathway, and evidence
insufficient = (
    s_seq < config.insufficient_signal_threshold and  # Default: 0.02
    s_path < 0.05 and 
    s_evd < 0.2
)

if evidence_gate:
    return "supported"
elif insufficient:
    return "insufficient"
else:
    return "consider"
```

**V2 Tier Logic** (`CONFIDENCE_V2=1`):
```python
# Tier I (supported): FDA on‑label OR ≥1 RCT OR (ClinVar‑Strong AND pathway_aligned)
if ("FDA-OnLabel" in badges or 
    "RCT" in badges or 
    ("ClinVar-Strong" in badges and s_path >= 0.2)):
    return "supported"

# Tier II (consider): ≥2 human studies MoA‑aligned OR 1 strong study + pathway_aligned
if (s_evd >= 0.6 or  # Strong evidence (proxy for strong study)
    (s_evd >= 0.4 and s_path >= 0.2)):  # Moderate evidence + pathway alignment
    return "consider"

# Tier III (insufficient): else
return "insufficient"
```

**Default Thresholds** (`config.py:25-28`):
- `EVIDENCE_GATE_THRESHOLD = 0.7`
- `PATHWAY_ALIGNMENT_THRESHOLD = 0.2`
- `INSUFFICIENT_SIGNAL_THRESHOLD = 0.02`

---

### **3.5 CONFIDENCE COMPUTATION - COMPLETE DEEP DIVE**

#### **3.5.1 Legacy Confidence** (`confidence_computation.py:69-108`):

**Tier-Based Approach**:
```python
if tier == "supported":
    confidence = 0.6 + 0.2 * max(seq_pct, path_pct)
elif tier == "consider":
    if config.fusion_active and max(seq_pct, path_pct) >= 0.7:
        confidence = 0.5 + 0.2 * max(seq_pct, path_pct)
    else:
        confidence = 0.3 + 0.1 * seq_pct + 0.1 * path_pct
else:  # insufficient
    max_sp = max(seq_pct, path_pct)
    min_sp = min(seq_pct, path_pct)
    base = 0.20 + 0.35 * max_sp + 0.15 * min_sp
    if config.fusion_active:
        confidence = max(0.25, base)
    else:
        confidence = base

# Insights modulation
confidence += 0.05 if func >= 0.6 else 0.0
confidence += 0.04 if chrom >= 0.5 else 0.0
confidence += 0.07 if ess >= 0.7 else 0.0
confidence += 0.02 if reg >= 0.6 else 0.0

# Alignment margin boost
margin = abs(seq_pct - path_pct)
if margin >= 0.2:
    confidence += 0.05
```

#### **3.5.2 V2 Confidence** (`CONFIDENCE_V2=1`) (`confidence_computation.py:111-160`):

**Linear S/P/E Formula**:
```python
# Convert tier to evidence score (E component)
if tier == "supported":
    e_score = 0.05
elif tier == "consider":
    e_score = 0.02
else:  # insufficient
    e_score = 0.00

# Calculate lifts
lifts = 0.0
lifts += 0.04 if func >= 0.6 else 0.0      # Functionality
lifts += 0.02 if chrom >= 0.5 else 0.0     # Chromatin
lifts += 0.02 if ess >= 0.7 else 0.0       # Essentiality
lifts += 0.02 if reg >= 0.6 else 0.0      # Regulatory
lifts = min(lifts, 0.08)  # Cap total lifts at +0.08

# Linear S/P/E formula
confidence = 0.5 * seq_pct + 0.2 * path_pct + 0.3 * e_score + lifts
confidence = clamp01(confidence)  # Clamp to [0, 1]
return round(confidence, 2)
```

**Formula**: `confidence = 0.5·S + 0.2·P + 0.3·E + lifts`
- **S weight**: 0.5 (sequence percentile)
- **P weight**: 0.2 (pathway percentile)
- **E weight**: 0.3 (evidence tier score: 0.05/0.02/0.00)
- **Lifts**: Max +0.08 from insights

#### **3.5.3 Additional Confidence Boosts** (`drug_scorer.py:140-168`):

**ClinVar Prior Boost**:
```python
if clinvar_prior > 0 and path_pct >= 0.2:
    confidence += min(0.1, clinvar_prior)
```

**Gene-Drug MoA Tie-Breaker**:
```python
if seq_scores and path_pct >= 0.2:
    primary_gene = (seq_scores[0].variant or {}).get("gene", "").upper()
    if primary_gene == "BRAF" and drug_name == "BRAF inhibitor":
        confidence += 0.01  # Target engagement bonus
    elif primary_gene in {"KRAS", "NRAS"} and drug_name == "MEK inhibitor":
        confidence += 0.01  # Downstream effector bonus
```

**Research-Mode Pathway Prior** (`RESEARCH_USE_PATHWAY_PRIOR=1`):
```python
if primary_gene in {"KRAS", "NRAS"} and drug_name == "MEK inhibitor":
    confidence += 0.02
if primary_gene == "BRAF" and drug_name == "BRAF inhibitor":
    confidence += 0.02
```

---

### **3.6 SPORADIC CANCER GATES - COMPLETE INTEGRATION**

#### **3.6.1 Gate 1: PARP Inhibitor Penalty** (`sporadic_gates.py:56-127`):

**Logic**:
- **Germline positive** → Full effect (1.0x)
- **Germline negative + HRD ≥42** → **RESCUED** (1.0x) ⚔️
- **Germline negative + HRD <42** → Reduced (0.6x)
- **Unknown germline + unknown HRD** → Conservative (0.8x)

**Formula**: `efficacy_score *= parp_penalty`

#### **3.6.2 Gate 2: Immunotherapy Boost** (`sporadic_gates.py:128-196`):

**Priority Order** (mutually exclusive, highest wins):
1. **TMB ≥20** → 1.35x boost (highest priority)
2. **MSI-High** → 1.30x boost (second priority)
3. **TMB ≥10 but <20** → 1.25x boost (lowest priority)

**Formula**: `efficacy_score *= io_boost_factor` (single factor, not multiplicative)

#### **3.6.3 Gate 3: Confidence Capping** (`sporadic_gates.py:197-232`):

**By Completeness Level**:
- **Level 0** (completeness <0.3): Cap at 0.4
- **Level 1** (0.3 ≤ completeness <0.7): Cap at 0.6
- **Level 2** (completeness ≥0.7): No cap

**Formula**: `confidence = min(confidence, cap)`

---

### **3.7 ABLATION MODES**

#### **Ablation Support** (`orchestrator.py:188-198`):
```python
ablation = (request.ablation_mode or ("SP" if fast_mode else "SPE")).upper()
use_S = "S" in ablation
use_P = "P" in ablation
use_E = "E" in ablation

# Shallow copies/masks
masked_seq_scores = seq_scores if use_S else []
masked_pathway_scores = pathway_scores if use_P else {}
masked_evidence = (evidence_results[i] if (use_E and i < len(evidence_results)) else None)
```

**Modes**:
- **"SPE"**: Full S/P/E (default in normal mode)
- **"SP"**: Sequence + Pathway only (default in fast mode)
- **"S"**: Sequence only
- **"P"**: Pathway only
- **"E"**: Evidence only

**Why**: Allows testing individual component contributions

---

### **3.8 KEY INSIGHTS & EDGE CASES**

#### **Sequence Scoring**:
1. **Hotspot floors applied TWICE**: Once at raw delta (1e-4), once at percentile (0.80-0.90)
2. **Truncation/frameshift**: Enforces `disruption = 1.0` for stop codons/frameshifts
3. **Forward/reverse symmetry**: Can be disabled via `evo_disable_symmetry` (default: disabled)
4. **Gene-specific calibration exists but not used**: `GeneCalibrationService` available but main flow uses `percentile_like()`
5. **Delta-only mode**: `evo_use_delta_only` skips exon scanning for speed

#### **Pathway Scoring**:
1. **Hardcoded weights**: Gene-to-pathway and drug-to-pathway weights are hardcoded
2. **Simplified mapping**: Only 2 pathways (ras_mapk, tp53) currently
3. **Normalization**: Based on empirical Evo2 ranges (1e-6 to 1e-4)

#### **Evidence Scoring**:
1. **MoA boost**: +0.10 per MoA hit, max +0.30
2. **Research-mode fallback**: Canonical hotspots get `prior = 0.2` if ClinVar empty
3. **Timeout handling**: Evidence timeout → `tier = "insufficient"`, `s_evd = 0.0`

#### **Confidence Computation**:
1. **Two versions**: Legacy (tier-based) and V2 (linear S/P/E)
2. **Feature flag**: `CONFIDENCE_V2=1` enables V2
3. **Lifts capped**: Total insights lifts capped at +0.08
4. **ClinVar boost**: +0.1 max when aligned (path_pct ≥ 0.2)

#### **Sporadic Gates**:
1. **PARP rescue**: HRD ≥42 rescues PARP for germline-negative patients
2. **IO boost priority**: TMB ≥20 > MSI-H > TMB ≥10 (mutually exclusive)
3. **Confidence capping**: Based on tumor context completeness (L0/L1/L2)

---

**Status**: 🔄 **ITERATION 4 IN PROGRESS** - AI Services & Model Integration  
**Next**: Complete I4 documentation, then move to I5 (Frontend Architecture)

---

## 📊 ITERATION 4: AI SERVICES & MODEL INTEGRATION 🔄 IN PROGRESS

### **4.1 MODAL DEPLOYMENT ARCHITECTURE**

#### **4.1.1 Core Pattern: HTTP-Based Communication**
**Critical Doctrine**: Backend communicates with Modal services via **HTTP**, NOT Modal SDK
- **Why**: Decouples backend from Modal SDK dependencies
- **Pattern**: Direct HTTP calls to Modal service URLs
- **Example** (`config.py:139-142`):
  ```python
  EVO_SERVICE_URL = "https://crispro--evo-service-evoservice-api.modal.run"
  EVO_URL_1B = "https://crispro--evo-service-evoservice1b-api-1b.modal.run"
  EVO_URL_7B = "https://crispro--evo-service-evoservice7b-api-7b.modal.run"
  EVO_URL_40B = "https://crispro--evo-service-evoservice-api.modal.run"
  ```

#### **4.1.2 Modal Service URL Naming Convention**
**Pattern**: `https://{workspace}--{app-name}-{class-name}-api.modal.run`
- **Workspace**: `crispro`
- **App Name**: Service name (e.g., `evo-service`, `zeta-oracle`)
- **Class Name**: Modal class name (e.g., `evoservice`, `evoservice1b`)
- **Suffix**: `-api` (Modal ASGI app endpoint)

#### **4.1.3 Dynamic Model URL Fallback** (`config.py:146-161`):
```python
def get_model_url(model_id: str) -> str:
    """Get URL for model with fallback logic, evaluated at runtime"""
    url_1b = os.getenv("EVO_URL_1B", "")
    url_7b = os.getenv("EVO_URL_7B", "...")
    url_40b = os.getenv("EVO_URL_40B", "...")
    
    if model_id == "evo2_1b":
        return url_1b or url_7b or url_40b  # Fallback chain
    elif model_id == "evo2_7b":
        return url_7b or url_40b
    elif model_id == "evo2_40b":
        return url_40b
    else:
        return url_1b or url_7b or url_40b
```

**Fallback Strategy**:
- **1B request** → Try 1B, fallback to 7B, then 40B
- **7B request** → Try 7B, fallback to 40B
- **40B request** → Only 40B (no fallback)

**Why**: Cost optimization (smaller models first) + resilience (fallback if unavailable)

---

### **4.2 EVO2 MODAL SERVICE ARCHITECTURE**

#### **4.2.1 Service Structure** (`src/services/evo_service/main.py`):

**Main Service Class** (`EvoService`):
- **GPU**: `H100:2` (2x H100 GPUs)
- **Volumes**: Model cache volume (`evo-model-cache`)
- **Scaledown**: 300 seconds (5 min idle before shutdown)
- **Timeout**: 1800 seconds (30 min max request time)
- **Default Model**: `evo2_1b_base` (can be overridden via `EVO_MODEL_ID`)

**Additional Services**:
- **EvoService7B**: `H100:1`, model `evo2_7b` (if `ENABLE_EVO_7B=1`)
- **EvoService1B**: `H100:1`, model `evo2_1b_base` (dedicated 1B service)

#### **4.2.2 Image Definition** (`evo2_image`):
```python
evo2_image = (
    modal.Image.from_registry("nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install(["build-essential", "cmake", "ninja-build", "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"])
    .env({"CC": "/usr/bin/gcc", "CXX": "/usr/bin/g++", "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True"})
    .run_commands("mkdir -p /tools/llvm/bin", "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar")
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx")
)
```

**Key Dependencies**:
- **CUDA 12.4.0**: For GPU acceleration
- **Evo2**: Cloned from ArcInstitute GitHub
- **Transformer Engine 1.13**: Required for Evo2
- **Numpy 1.22.0**: Pinned version (build compatibility)

#### **4.2.3 Endpoints** (`EvoService`):

**1. `/score_delta`** (POST):
- **Input**: `{"ref_sequence": str, "alt_sequence": str}`
- **Output**: `{"ref_likelihood": float, "alt_likelihood": float, "delta_score": float}`
- **Use**: Direct sequence scoring (no Ensembl fetch)

**2. `/score_batch`** (POST):
- **Input**: `{"pairs": [{"ref_sequence": str, "alt_sequence": str}, ...]}`
- **Output**: `{"results": [{"ref_likelihood": float, "alt_likelihood": float, "delta_score": float}, ...]}`
- **Use**: Batch sequence scoring

**3. `/score_variant`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "window": 8192}`
- **Output**: Same as `/score_delta` + `{"window_start": int, "window_end": int, "variant_index": int}`
- **Use**: Single-window variant scoring (fetches from Ensembl)

**4. `/score_variant_multi`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "windows": [1024, 2048, 4096, 8192]}`
- **Output**: `{"deltas": [{"window": int, "delta": float}, ...], "min_delta": float, "window_used": int}`
- **Use**: Multi-window variant scoring (picks most negative delta)

**5. `/score_variant_exon`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "flank": 600}`
- **Output**: `{"exon_delta": float, "window_used": int}`
- **Use**: Tight-window (exon-context) scoring

**6. `/score_variant_profile`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "flank": 600, "radius": 100}`
- **Output**: `{"profile": [{"offset": int, "delta": float}, ...], "peak_delta": float, "peak_offset": int}`
- **Use**: Local delta profile (applies ALT at offsets ±100 bp)

**7. `/score_variant_probe`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A"}`
- **Output**: `{"probes": [{"alt": str, "delta": float}, ...], "top_alt": str, "top_delta": float}`
- **Use**: 3-alt sensitivity probe (scores ref→A/C/G/T)

**8. `/generate`** (POST):
- **Input**: `{"prompt": str, "n_tokens": int}`
- **Output**: `{"job_id": str}` (async job submission)
- **Use**: Sequence generation (background job)

**9. `/status/{job_id}`** (GET):
- **Output**: `{"status": "pending"|"running"|"complete"|"failed", "result": {...}, "error": str}`
- **Use**: Poll generation job status

#### **4.2.4 Ensembl Integration**:
- **API**: `https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}`
- **Region Format**: `{chrom}:{start}-{end}:1`
- **Assembly Mapping**: `GRCh38`/`hg38` → `GRCh38`, else → `GRCh37`
- **Validation**: Checks fetched base matches provided REF (allows `N`)

#### **4.2.5 Shared State**:
- **Job Status**: `modal.Dict.from_name("evo-job-status")` (shared across containers)
- **Model Cache**: `modal.Volume.from_name("evo-model-cache")` (persistent model weights)

---

### **4.3 OTHER MODAL SERVICES**

#### **4.3.1 Zeta Oracle** (`src/services/oracle/main.py`):
- **Purpose**: Variant scoring using AlphaMissense + ESM + Evo2 fusion
- **GPU**: `H100:2`
- **Image**: Same `evo2_image` (Evo2 dependencies)
- **Endpoints**:
  - `/score_variants`: Fused scoring (AlphaMissense + ESM + Evo2)
  - `/invoke`: General-purpose variant analysis
- **URL**: `https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run`

#### **4.3.2 Zeta Forge** (`src/services/forge/main.py`):
- **Purpose**: Generate novel biologic inhibitors using Evo2
- **GPU**: `H100:2`
- **Image**: Same `evo2_image`
- **Endpoints**:
  - `/generate_inhibitor`: Submit generation job (async)
  - `/status/{job_id}`: Poll job status
- **Doctrine**: "Ambush" - reverse complement of target motif
- **URL**: `https://crispro--zeta-forge-zetaforge-api.modal.run`

#### **4.3.3 Fusion Engine** (`src/services/fusion_engine/main.py`):
- **Purpose**: Lightweight variant scoring (AlphaMissense + ESM, NO Evo2)
- **GPU**: `any` (generic GPU sufficient)
- **Image**: Lightweight (no Evo2 dependencies)
- **Volume**: `alphamissense-data` (AlphaMissense parquet file)
- **Endpoints**:
  - `/score_variants`: Fused scoring (AlphaMissense + ESM only)
- **Why**: Fast scoring without heavy Evo2 dependencies
- **URL**: `https://crispro--fusion-engine-fusionengine-api.modal.run`

#### **4.3.4 Boltz Service** (`src/services/boltz_service/main.py`):
- **Purpose**: Protein structure prediction using Boltz-2
- **GPU**: `H100`
- **Model**: Boltz-2 from Hugging Face (`boltz-community/boltz-2`)
- **Volume**: `boltz-models` (persistent model weights)
- **Endpoints**:
  - `/v1/predict_structure`: Structure prediction
  - `/v1/predict_interaction`: Binding affinity prediction
- **Features**:
  - **Fast mode**: 2-5 min predictions
  - **iPTM scores**: Binding affinity (>0.7 = high confidence)
  - **SMILES handling**: Direct compound analysis
- **URL**: `https://crispro--boltz-service-boltzservice-api.modal.run`

#### **4.3.5 Genesis Engine** (`src/services/genesis_engine/main.py`):
- **Purpose**: Variant analysis using Evo2 40B
- **GPU**: `H100`
- **Model**: `evo2_40b`
- **Endpoints**:
  - `/analyze_single_variant`: Single variant analysis with enrichment
- **URL**: `https://crispro--genesis-engine-genesisengine-api.modal.run`

---

### **4.4 ALPHAFOLD 3 INTEGRATION STATUS**

#### **4.4.1 Current State**:
- **Status**: ❌ **NOT AUTOMATED** - Manual workflow only
- **Approach**: JSON input → AlphaFold Server → JSON output → Review
- **Why Manual**: Faster than automated deployment, already validated (15/15 guides, 100% pass rate)

#### **4.4.2 Manual Workflow**:
1. **Generate JSON**: Structure description (protein/DNA/RNA chains)
2. **Submit to AF Server**: Web interface or API
3. **Download Results**: pLDDT/PAE metrics, structure files
4. **Review**: Validate quality metrics

#### **4.4.3 AF Server Capabilities**:
- **Molecule Types**: `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand`, `ion`
- **Complexes**: Protein-DNA complexes (e.g., Cas9:gRNA:target DNA)
- **Multiple Chains**: Full CRISPR complexes
- **PTMs**: Post-translational modifications
- **DNA Modifications**: 5mC, 8-oxoG, etc.
- **Template Control**: `useStructureTemplate: false` to disable PDB templates

#### **4.4.4 Boltz as Alternative**:
- **Status**: ✅ **FULLY AUTOMATED** on Modal
- **Speed**: 2-5 min per structure (vs 60+ min for ColabFold)
- **Use Cases**: Fast structure validation, binding affinity prediction
- **Limitation**: Does not model complexes/interfaces (single-chain only)

#### **4.4.5 Future Roadmap**:
- **Phase 1**: ESMFold for fast single-chain predictions (1-2 min/structure)
- **Phase 2**: AlphaFold3 integration (pending Google DeepMind weights approval)
- **Phase 3**: Boltz for binding affinity (NO DiffDock needed)

---

### **4.5 BACKEND → MODAL COMMUNICATION PATTERNS**

#### **4.5.1 Request Flow** (`routers/evo.py`):
```python
# 1. Get model URL with fallback
target_url = get_model_url(model_id)

# 2. Make HTTP request
async with httpx.AsyncClient(timeout=EVO_TIMEOUT) as client:
    response = await client.post(
        f"{target_url}/score_variant_multi",
        json=payload,
        headers={"Content-Type": "application/json"}
    )
    response.raise_for_status()
    return response.json()
```

#### **4.5.2 Fallback Logic** (`routers/evo.py:235-236`):
```python
# If 1B or 7B fails, try falling back to 40B (if allowed)
if (not DISABLE_EVO_FALLBACK) and model_id in ["evo2_1b", "evo2_7b"] and EVO_URL_40B:
    # Retry with 40B
    fallback_request = {**request, "model_id": "evo2_40b"}
    response = await client.post(f"{EVO_URL_40B}/generate", json=fallback_request)
    result["warning"] = f"Requested {model_id} failed, fell back to evo2_40b"
```

#### **4.5.3 Timeout Configuration** (`config.py:143`):
```python
EVO_TIMEOUT = Timeout(60.0, connect=10.0)  # 60s total, 10s connect
```

#### **4.5.4 Health Checks** (`routers/evo.py:383-436`):
- **Endpoint**: `/health` (GET)
- **Purpose**: Check service availability
- **Response**: `{"status": "healthy"|"unhealthy", "url": str}`

---

### **4.6 MODAL DEPLOYMENT PATTERNS**

#### **4.6.1 Image Definition Pattern**:
1. **Base Image**: CUDA image for GPU services, Debian slim for CPU services
2. **System Dependencies**: `apt_install` for build tools, libraries
3. **Python Dependencies**: `pip_install` for Python packages
4. **Local Code**: `add_local_dir` for project code
5. **Environment Variables**: `.env()` for configuration

#### **4.6.2 Service Class Pattern**:
```python
@app.cls(
    gpu="H100:2",           # GPU allocation
    volumes={...},          # Persistent storage
    scaledown_window=300,   # Idle timeout
    timeout=1800            # Max request time
)
class ServiceName:
    @modal.enter()
    def load_model(self):
        # Load model on container startup
        self.model = Model()
    
    @modal.asgi_app()
    def api(self):
        # FastAPI app
        return fastapi_app
```

#### **4.6.3 Shared State Pattern**:
- **Job Status**: `modal.Dict.from_name("service-jobs")` (shared across containers)
- **Model Cache**: `modal.Volume.from_name("service-models")` (persistent storage)

#### **4.6.4 Async Job Pattern**:
```python
# 1. Submit job (returns job_id)
job_id = str(uuid.uuid4())
job_status_dict[job_id] = {"status": "pending"}
background_task.spawn(job_id, request)
return {"job_id": job_id}

# 2. Poll status
@fastapi_app.get("/status/{job_id}")
def get_status(job_id: str):
    return job_status_dict[job_id]
```

---

### **4.7 KEY INSIGHTS**

#### **Modal Service Architecture**:
1. **HTTP-Based**: Backend uses HTTP, not Modal SDK (decoupling)
2. **Fallback Chains**: Model URL fallback (1B → 7B → 40B)
3. **Shared State**: `modal.Dict` for job status, `modal.Volume` for model cache
4. **Scaledown**: 300s idle timeout (cost optimization)

#### **Evo2 Service**:
1. **Multi-Model**: 1B, 7B, 40B variants (separate services)
2. **Ensembl Integration**: Fetches reference sequences from Ensembl API
3. **Multi-Window**: Tests multiple window sizes, picks best delta
4. **Async Jobs**: Generation uses background jobs with status polling

#### **Other Services**:
1. **Oracle**: Fused scoring (AlphaMissense + ESM + Evo2)
2. **Forge**: Sequence generation (inhibitor design)
3. **Fusion**: Lightweight scoring (AlphaMissense + ESM only)
4. **Boltz**: Structure prediction (2-5 min, automated)

#### **AlphaFold 3**:
1. **Manual Workflow**: JSON → AF Server → Review (faster than automated)
2. **Boltz Alternative**: Automated structure prediction (single-chain only)
3. **Future**: ESMFold for fast validation, AF3 for complexes (pending weights)

---

**Status**: 🔄 **ITERATION 5 IN PROGRESS** - Frontend Architecture & User Experience  
**Next**: Complete I5 documentation, then move to I6 (Clinical Systems & Workflows)

---

## 📊 ITERATION 5: FRONTEND ARCHITECTURE & USER EXPERIENCE 🔄 IN PROGRESS

### **5.1 TECHNOLOGY STACK**

#### **5.1.1 Core Framework**:
- **React 18.2.0**: Component-based UI library
- **React Router DOM 6.4.4**: Client-side routing
- **Vite 3**: Build tool and dev server
- **TypeScript**: Partial adoption (some `.tsx` files)

#### **5.1.2 UI Libraries**:
- **Material-UI (MUI) 6.5.0**: Primary component library
  - Components: `Box`, `Card`, `Typography`, `Button`, `Chip`, `LinearProgress`, `Accordion`, `Tabs`, `Dialog`
  - Icons: `@mui/icons-material`
  - Charts: `@mui/x-charts`
- **Tailwind CSS 3.2.4**: Utility-first CSS framework
- **Radix UI**: Headless UI components (`@radix-ui/react-select`, `@radix-ui/react-tabs`, `@radix-ui/react-slot`)
- **Class Variance Authority**: Component variant management

#### **5.1.3 State Management**:
- **React Context API**: Primary state management
  - `AuthContext`: User authentication
  - `SporadicContext`: Sporadic cancer workflow state
  - `CoPilotContext`: AI assistant state
  - `ActivityContext`: Activity tracking
  - `AnalysisHistoryContext`: Analysis history
  - `AgentContext`: Agent system state
  - `ClinicalGenomicsContext`: Clinical genomics state
- **Local Storage**: Persistent state (e.g., Kanban tasks)
- **React Hooks**: Custom hooks for API calls and data fetching

#### **5.1.4 Additional Libraries**:
- **Supabase 2.56.0**: Authentication and database
- **React Markdown 9.0.1**: Markdown rendering
- **React Joyride 2.9.3**: User onboarding/tours
- **React Toastify 11.0.5**: Toast notifications
- **Notistack 3.0.2**: Snackbar notifications
- **React Spring 10.0.1**: Animation library
- **DnD Kit**: Drag-and-drop functionality

---

### **5.2 APPLICATION STRUCTURE**

#### **5.2.1 Directory Organization**:
```
src/
├── components/          # Reusable UI components
│   ├── ayesha/         # Ayesha-specific components
│   ├── CoPilot/        # AI assistant components
│   ├── common/         # Shared components
│   ├── dashboard/      # Dashboard components
│   └── ...
├── pages/              # Page-level components (routes)
├── context/            # React Context providers
├── hooks/              # Custom React hooks
├── services/           # API clients and services
├── utils/              # Utility functions
├── config/             # Configuration files
├── constants/          # Constants and enums
└── features/           # Feature modules
```

#### **5.2.2 Entry Point** (`main.jsx`):
```jsx
import { BrowserRouter as Router } from "react-router-dom";
import { StateContextProvider } from "./context";
import App from "./App";

root.render(
  <Router>
    <StateContextProvider>
      <App />
    </StateContextProvider>
  </Router>
);
```

**Provider Hierarchy** (`App.jsx:79-189`):
```jsx
<ErrorBoundary>
  <AuthProvider>
    <AgentProvider>
      <SporadicProvider>
        <CoPilotProvider>
          <AnalysisHistoryProvider>
            <ActivityProvider>
              {/* App content */}
            </ActivityProvider>
          </AnalysisHistoryProvider>
        </CoPilotProvider>
      </SporadicProvider>
    </AgentProvider>
  </AuthProvider>
</ErrorBoundary>
```

---

### **5.3 ROUTING ARCHITECTURE**

#### **5.3.1 Route Structure** (`App.jsx:93-172`):

**Auth Routes** (Public):
- `/login` → `Login`
- `/signup` → `Signup`

**Admin Routes** (Protected):
- `/admin/dashboard` → `AdminDashboard`
- `/admin/users` → `AdminUsers`

**Main Routes**:
- `/` → `Home`
- `/dashboard` → `DoctorDashboard`
- `/profile` → `Profile`
- `/onboarding` → `Onboarding`
- `/medical-records` → `MedicalRecords`
- `/medical-records/:id` → `SingleRecordDetails`
- `/research` → `Research`
- `/mutation-explorer` → `MutationExplorer`
- `/agent-dashboard` → `AgentDashboard`
- `/agents` → `AgentsPage`
- `/agent-studio` → `AgentStudio`

**Ayesha Routes**:
- `/ayesha-complete-care` → `AyeshaCompleteCare`
- `/ayesha-trials` → `AyeshaTrialExplorer`
- `/ayesha-dossiers` → `AyeshaDossierBrowser`
- `/ayesha-dossiers/:nct_id` → `AyeshaDossierDetail`
- `/sporadic-cancer` → `SporadicCancerPage`
- `/ayesha-twin-demo` → `AyeshaTwinDemo`

**Research/Design Routes**:
- `/metastasis` → `MetastasisDashboard`
- `/crispr-designer` → `CrisprDesigner`
- `/protein-synthesis` → `ProteinSynthesis`
- `/structure-predictor` → `StructurePredictor`
- `/validate` → `HypothesisValidator`
- `/food-validator` → `FoodValidatorAB`
- `/batch-food-validator` → `BatchFoodValidator`

**Clinical Routes**:
- `/clinical-genomics` → `ClinicalGenomicsCommandCenter`
- `/myeloma-digital-twin` → `MyelomaDigitalTwin`
- `/radonc-co-pilot` → `RadOncCoPilot`
- `/threat-assessor` → `ThreatAssessor`

#### **5.3.2 Protected Routes**:
- **Pattern**: `<ProtectedRoute><Component /></ProtectedRoute>`
- **Implementation**: Checks authentication status before rendering

---

### **5.4 STATE MANAGEMENT PATTERNS**

#### **5.4.1 Context API Pattern**:

**Example: SporadicContext** (`context/SporadicContext.jsx`):
```jsx
export const SporadicProvider = ({ children }) => {
  const [germlineStatus, setGermlineStatus] = useState('unknown');
  const [tumorContext, setTumorContext] = useState(null);
  const [contextId, setContextId] = useState(null);
  const [dataLevel, setDataLevel] = useState('L0');

  const updateTumorContext = useCallback((data) => {
    if (data?.tumor_context) {
      setTumorContext(data.tumor_context);
      setContextId(data.context_id);
      // Determine data level from completeness score
      const completeness = data.tumor_context.completeness_score || 0;
      if (completeness >= 0.7) setDataLevel('L2');
      else if (completeness >= 0.3) setDataLevel('L1');
      else setDataLevel('L0');
    }
  }, []);

  const getEfficacyPayload = useCallback((basePayload) => {
    return {
      ...basePayload,
      germline_status: germlineStatus,
      tumor_context: tumorContext,
    };
  }, [germlineStatus, tumorContext]);

  return (
    <SporadicContext.Provider value={{ ...state, ...actions }}>
      {children}
    </SporadicContext.Provider>
  );
};
```

**Key Features**:
- **State Slices**: Separate state for different concerns
- **Actions**: `useCallback`-wrapped action functions
- **Computed Values**: Derived state (e.g., `hasTumorContext`, `isSporadic`)
- **Integration Helpers**: `getEfficacyPayload` injects context into API calls

#### **5.4.2 Local Storage Pattern**:
```jsx
// Load from localStorage
const [tasks, setTasks] = useState(() => {
  const savedTasks = localStorage.getItem(KANBAN_TASKS_KEY);
  return savedTasks ? JSON.parse(savedTasks) : [];
});

// Save to localStorage
useEffect(() => {
  localStorage.setItem(KANBAN_TASKS_KEY, JSON.stringify(tasks));
}, [tasks]);
```

---

### **5.5 API INTEGRATION PATTERNS**

#### **5.5.1 API Client Hook** (`hooks/useApiClient.js`):
```jsx
export default function useApiClient(modelId) {
  const API_BASE_URL = import.meta.env.VITE_API_ROOT || 'http://127.0.0.1:8023';

  const client = useMemo(() => {
    const post = async (endpoint, payload = {}, opts = {}) => {
      const controller = new AbortController();
      const timeoutMs = opts.timeoutMs || DEFAULT_TIMEOUT_MS;
      const timer = setTimeout(() => controller.abort(), timeoutMs);
      
      try {
        const body = JSON.stringify({ ...payload, model_id: modelId || 'evo2_7b' });
        const res = await fetch(`${API_BASE_URL}${endpoint}`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body,
          signal: controller.signal,
        });
        
        const json = await res.json().catch(() => ({}));
        if (!res.ok) {
          throw new Error(json?.detail || `HTTP ${res.status}`);
        }
        return json;
      } finally {
        clearTimeout(timer);
      }
    };
    return { post };
  }, [API_BASE_URL, modelId]);

  return client;
}
```

**Features**:
- **Timeout**: 10 minutes default (configurable)
- **AbortController**: Request cancellation
- **Error Handling**: Extracts `detail` from error response
- **Model ID Injection**: Automatically adds `model_id` to payload

#### **5.5.2 Advanced API Client** (`components/ClinicalGenomicsCommandCenter/utils/genomicsUtils.js`):
```jsx
export async function apiPost(path, body, { signal, useCache = true, skipRetry = false } = {}) {
  const url = `${API_BASE}${path}`;
  
  // Check cache first
  if (useCache) {
    const cacheKey = getCacheKey(path, body);
    const cached = getCached(cacheKey);
    if (cached) return cached;
  }
  
  const controller = new AbortController();
  const abortSignal = signal || controller.signal;
  const timeoutId = setTimeout(() => controller.abort(), DEFAULT_TIMEOUT);
  
  const doFetch = async (attempt = 1) => {
    try {
      const response = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
        body: JSON.stringify(body),
        signal: abortSignal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.detail || `HTTP ${response.status}`);
      }
      
      const data = await response.json();
      
      // Cache successful result
      if (useCache) {
        const cacheKey = getCacheKey(path, body);
        setCache(cacheKey, data);
      }
      
      return data;
    } catch (error) {
      clearTimeout(timeoutId);
      
      // Retry logic (exponential backoff)
      if (!skipRetry && attempt < 3 && error.name !== 'AbortError') {
        const delay = attempt * 1000; // 1s, 2s
        await new Promise(resolve => setTimeout(resolve, delay));
        return doFetch(attempt + 1);
      }
      
      throw error;
    }
  };
  
  return doFetch();
}
```

**Advanced Features**:
- **Caching**: 10-minute TTL, keyed by path + body
- **Retry Logic**: Exponential backoff (2 retries)
- **Timeout**: 60 seconds default
- **Abort Support**: Signal-based cancellation

#### **5.5.3 Custom Hooks Pattern** (`hooks/useMetastasis.js`):
```jsx
export function useMetastasisAssess(params, enabled = true) {
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [lastFetch, setLastFetch] = useState(null);

  const cacheKey = JSON.stringify(params);
  const CACHE_TTL_MS = 10 * 60 * 1000; // 10 minutes

  const fetchAssessment = useCallback(async () => {
    if (!enabled || !params.mutations || params.mutations.length === 0) {
      return;
    }

    // Check cache freshness
    const now = Date.now();
    if (lastFetch && (now - lastFetch) < CACHE_TTL_MS && data) {
      return; // Use cached data
    }

    setLoading(true);
    setError(null);

    try {
      const response = await fetch(`${API_ROOT}/api/metastasis/assess`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(params)
      });

      if (!response.ok) {
        throw new Error(`Assessment failed (${response.status})`);
      }

      const result = await response.json();
      setData(result);
      setLastFetch(now);
    } catch (err) {
      setError(err.message);
      setData(null);
    } finally {
      setLoading(false);
    }
  }, [cacheKey, enabled, lastFetch, data]);

  useEffect(() => {
    fetchAssessment();
  }, [fetchAssessment]);

  return { data, loading, error, refetch: () => fetchAssessment() };
}
```

**Pattern**:
- **State**: `data`, `loading`, `error`
- **Caching**: TTL-based cache with `lastFetch` timestamp
- **Auto-fetch**: `useEffect` triggers on dependency changes
- **Manual Refetch**: `refetch` function for explicit refresh

---

### **5.6 COMPONENT PATTERNS**

#### **5.6.1 Ayesha Components** (`components/ayesha/`):

**DrugRankingPanel** (`DrugRankingPanel.jsx`):
- **Purpose**: Display ranked drug recommendations
- **Props**: `drugs[]`, `onViewDetails`
- **Features**:
  - Tier badges (supported/consider/insufficient)
  - Evidence badges (RCT, Guideline, ClinVar-Strong)
  - Efficacy score visualization (LinearProgress)
  - Expandable details (Accordion)
  - MUI-based styling

**FoodRankingPanel** (`FoodRankingPanel.jsx`):
- **Purpose**: Display food/supplement recommendations
- **Similar pattern** to DrugRankingPanel

**IntegratedConfidenceBar** (`IntegratedConfidenceBar.jsx`):
- **Purpose**: Visualize integrated confidence scores
- **Shows**: Drug + Food confidence breakdown

**CA125Tracker** (`CA125Tracker.jsx`):
- **Purpose**: Track CA-125 biomarker over time
- **Features**: Time series visualization, response forecasting

**ResistanceAlertBanner** (`ResistanceAlertBanner.jsx`):
- **Purpose**: Display resistance warnings
- **Features**: Alert styling, actionable recommendations

#### **5.6.2 CoPilot Components** (`components/CoPilot/`):

**CoPilot** (Main Component):
- **Purpose**: AI assistant interface
- **Features**: Chat interface, context awareness, action execution

**Q2CRouter** (Question-to-Capability Router):
- **Purpose**: Route user questions to appropriate capabilities
- **Pattern**: Intent classification → Capability selection → Action execution

**Evidence Components**:
- **Purpose**: Display clinical evidence
- **Features**: PubMed citations, ClinVar data, pathway information

**Action Components**:
- **Purpose**: Execute actions (e.g., run analysis, generate report)
- **Features**: Status tracking, progress indicators

#### **5.6.3 Common Components** (`components/common/`):
- **ToolRunner**: Execute tools/analyses
- **LoadingSkeleton**: Loading state placeholders
- **ErrorBoundary**: Error handling wrapper

---

### **5.7 STYLING PATTERNS**

#### **5.7.1 MUI Styling** (Primary):
```jsx
<Box sx={{ p: 4, maxWidth: 1400, mx: 'auto' }}>
  <Card sx={{ p: 3 }}>
    <Typography variant="h5" sx={{ fontWeight: 'bold' }}>
      Title
    </Typography>
    <LinearProgress variant="determinate" value={75} />
  </Card>
</Box>
```

**Features**:
- **SX Prop**: Inline styling with theme access
- **Theme System**: Consistent colors, spacing, typography
- **Responsive**: Breakpoint-based responsive design

#### **5.7.2 Tailwind CSS** (Secondary):
```jsx
<div className="sm:-8 relative flex min-h-screen flex-row bg-white p-4">
  <div className="relative mr-10 hidden sm:flex">
    <Sidebar />
  </div>
</div>
```

**Features**:
- **Utility Classes**: Rapid styling
- **Responsive**: Breakpoint prefixes (`sm:`, `md:`, `lg:`)
- **Custom Configuration**: Tailwind config for custom values

#### **5.7.3 Styled Components** (MUI):
```jsx
const StyledModal = styled(Modal)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '16px',
}));
```

---

### **5.8 KEY INSIGHTS**

#### **Frontend Architecture**:
1. **Hybrid UI Libraries**: MUI (primary) + Tailwind (utility)
2. **Context-Heavy**: Extensive use of React Context for state
3. **Hook-Based**: Custom hooks for API calls and data fetching
4. **Component Organization**: Feature-based + shared components

#### **State Management**:
1. **Context API**: Primary pattern (no Redux)
2. **Local Storage**: Persistent state (tasks, preferences)
3. **Derived State**: Computed values in context providers
4. **Integration Helpers**: Context-aware API payload builders

#### **API Integration**:
1. **Fetch API**: Native fetch (no axios)
2. **AbortController**: Request cancellation
3. **Caching**: TTL-based in-memory cache (10 minutes)
4. **Retry Logic**: Exponential backoff (2-3 attempts)
5. **Error Handling**: Extract `detail` from error responses

#### **Component Patterns**:
1. **Panel Components**: Ranking panels for drugs/foods
2. **Card Components**: MUI Card-based layouts
3. **Modal/Dialog**: MUI Dialog for overlays
4. **Accordion**: Expandable details sections
5. **Loading States**: Skeleton loaders + progress indicators

#### **Routing**:
1. **React Router v6**: Latest routing library
2. **Protected Routes**: Auth-based route protection
3. **Nested Routes**: Patient-scoped routes (`/medical-records/:patientId/research`)
4. **Route Organization**: Feature-based route grouping

---

**Status**: 🔄 **ITERATION 6 IN PROGRESS** - Clinical Systems & Workflows Deep Dive  
**Next**: Complete I6 documentation, then move to I7 (Research & Design Systems)

---

## 🏥 ITERATION 6: CLINICAL SYSTEMS & WORKFLOWS 🔄 IN PROGRESS

### **6.1 AYESHA COMPLETE CARE ORCHESTRATOR**

#### **6.1.1 Overview** (`routers/ayesha_orchestrator_v2.py`):
- **Purpose**: Unified care plan orchestration for AK (Stage IVB HGSOC)
- **Endpoint**: `POST /api/ayesha/complete_care_v2`
- **Integration**: Coordinates 7 clinical services + 3 SAE services + Resistance Prophet

#### **6.1.2 Orchestrated Services**:
1. **Clinical Trials** (`_call_ayesha_trials`):
   - Frontline, NYC metro, transparent reasoning
   - Includes SOC recommendation + CA-125 intelligence
2. **Drug Efficacy (WIWFM)** (`_call_drug_efficacy`):
   - Full S/P/E if `tumor_context` provided
   - Returns "awaiting NGS" message if no tumor data
3. **Food Validator** (`_call_food_validator`):
   - Optional supplement/nutrition recommendations
   - Calls `/api/hypothesis/validate_food_dynamic`
4. **Resistance Playbook** (`_call_resistance_playbook`):
   - Next-line planning based on resistance mechanisms
   - Requires `tumor_context` (returns "awaiting NGS" if missing)
5. **CA-125 Intelligence** (via trials endpoint):
   - Burden classification, response forecast, resistance signals
6. **SOC Recommendation** (via trials endpoint):
   - NCCN-aligned carboplatin + paclitaxel + bevacizumab
7. **Resistance Prophet** (opt-in):
   - Predicts resistance 3-6 months early
   - Manager Q7: Opt-in via `include_resistance_prediction=true`

#### **6.1.3 Phase 1 SAE Services** (Pre-NGS):
- **Next-Test Recommender**: HRD → ctDNA → SLFN11 → ABCB1 priority
- **Hint Tiles**: Max 4 tiles (Test → Trials → Monitor → Avoid)
- **Mechanism Map**: Pre-NGS (all gray) vs Post-NGS (color-coded)

#### **6.1.4 Phase 2 SAE Services** (Post-NGS):
- **SAE Features**: DNA repair capacity, mechanism vector, resistance signals
- **Resistance Alert**: 2-of-3 trigger rule (HRD drop, DNA repair drop, CA-125 inadequate)

#### **6.1.5 Request Schema** (`CompleteCareV2Request`):
```python
ca125_value: float  # Current CA-125 (U/mL)
stage: str  # "IVB"
treatment_line: str  # "first-line"
germline_status: str  # "negative"
has_ascites: bool
has_peritoneal_disease: bool
tumor_context: Optional[Dict]  # NGS data (somatic mutations, HRD, TMB, MSI)
drug_query: Optional[str]  # Specific drug to evaluate
food_query: Optional[str]  # Food/supplement to validate
include_trials: bool = True
include_soc: bool = True
include_ca125: bool = True
include_wiwfm: bool = True
include_food: bool = False
include_resistance: bool = False
include_resistance_prediction: bool = False  # Manager Q7: Opt-in
```

---

### **6.2 CA-125 INTELLIGENCE SERVICE**

#### **6.2.1 Purpose** (`services/ca125_intelligence.py`):
- Burden classification (MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE)
- Response forecast (cycle 3, cycle 6 targets)
- Resistance signal detection (on-therapy rise, inadequate response)
- Monitoring strategy recommendations

#### **6.2.2 Burden Classification**:
```python
BURDEN_THRESHOLDS = {
    "MINIMAL": (0, 100),
    "MODERATE": (100, 500),
    "SIGNIFICANT": (500, 1000),
    "EXTENSIVE": (1000, float('inf'))
}
```

#### **6.2.3 Response Expectations** (GOG-218, ICON7):
- **Cycle 3**: ≥70% drop expected
- **Cycle 6**: ≥90% drop expected
- **Complete Response**: <35 U/mL target
- **Resistance Signal**: <50% drop by cycle 3

#### **6.2.4 Resistance Signal Detection**:
1. **On-therapy rise**: CA-125 increases while on treatment
2. **Inadequate response**: <50% drop by cycle 3
3. **Minimal response**: <30% drop by any cycle post-cycle 2

#### **6.2.5 Monitoring Strategy**:
- **On treatment**: Every 3 weeks (before each cycle)
- **Pre-treatment (EXTENSIVE)**: Every 2 weeks
- **Pre-treatment (SIGNIFICANT)**: Every 4 weeks
- **Surveillance**: Every 3 months

---

### **6.3 RESISTANCE DETECTION SERVICE**

#### **6.3.1 Purpose** (`services/resistance_detection_service.py`):
- **2-of-3 Trigger Rule** (Manager C7):
  - HRD drop >= 15 points
  - DNA repair capacity drop >= 0.20
  - CA-125 inadequate response (on-therapy rise OR <50% drop by cycle 3)
- **HR Restoration Pattern** (Manager R2):
  - HRD drop + DNA repair drop (coherent signal)
  - Immediate alert (don't wait for radiology)

#### **6.3.2 Resistance Alert Structure**:
```python
@dataclass
class ResistanceAlert:
    resistance_detected: bool  # True if 2-of-3 triggers met
    hr_restoration_suspected: bool  # True if HR restoration pattern
    immediate_alert: bool  # True if should alert NOW
    triggers_met: List[str]  # ["hrd_drop", "dna_repair_drop", "ca125_inadequate"]
    trigger_count: int  # Number of triggers met (need ≥2)
    recommended_actions: List[str]  # Clinical actions
    recommended_trials: List[str]  # ATR/CHK1, WEE1 trials
```

#### **6.3.3 Recommended Actions**:
- **IMMEDIATE**: Alert oncologist, order updated HRD test, order ctDNA panel
- **If HR restoration**: Consider switching from PARP
- **Trial recommendations**: ATR inhibitor trials, ATR + PARP combos, CHK1/WEE1 inhibitors

---

### **6.4 RESISTANCE PROPHET SERVICE**

#### **6.4.1 Purpose** (`services/resistance_prophet_service.py`):
- **Predicts resistance 3-6 months early** (before clinical progression)
- **Phase 1**: DNA repair restoration + pathway escape (NO CA-125)
- **Phase 1b**: Add CA-125 kinetics (prospective, Ayesha live)
- **Manager Q7**: Opt-in via `include_resistance_prediction=true`

#### **6.4.2 Resistance Signals**:
1. **DNA Repair Restoration** (Signal 1):
   - Threshold: >20% increase from baseline
   - Indicates PARP resistance mechanism
2. **Pathway Escape** (Signal 2):
   - Threshold: >=30% shift in mechanism vector (L2 distance)
   - Indicates bypass resistance mechanism
3. **CA-125 Kinetics** (Signal 3, Phase 1b+):
   - On-therapy rise or inadequate response
   - Uses existing CA125Intelligence service

#### **6.4.3 Risk Stratification** (Manager Q9):
- **HIGH**: Probability >=0.70 AND >=2 signals
- **MEDIUM**: 0.50-0.69 OR exactly 1 signal
- **LOW**: <0.50 probability
- **Manager Q15**: Cap at MEDIUM if no CA-125 unless >=2 non-CA-125 signals

#### **6.4.4 Urgency & Actions** (Manager Q10):
- **CRITICAL** (HIGH risk):
  - ESCALATE_IMAGING within 1 week
  - CONSIDER_SWITCH within 2 weeks
  - REVIEW_RESISTANCE_PLAYBOOK within 1 week
- **ELEVATED** (MEDIUM risk):
  - MONITOR_WEEKLY x4 weeks
  - REASSESS after 4 weeks
- **ROUTINE** (LOW risk):
  - Routine monitoring per standard of care

#### **6.4.5 Next-Line Options** (Manager Q11):
- Consults `ResistancePlaybookService` for mechanism-specific strategies
- Consults `TreatmentLineService` for appropriateness

---

### **6.5 RESISTANCE PLAYBOOK SERVICE**

#### **6.5.1 Purpose** (`services/resistance_playbook_service.py`):
- Predicts resistance mechanisms (HR restoration, ABCB1, MAPK, PI3K, SLFN11)
- Ranks combo strategies (PARP+ATR, PARP+VEGF, IO combos, MAPK/PI3K)
- Recommends next-line switches
- Generates trial keywords for filtering

#### **6.5.2 Resistance Detection Rules**:
1. **HR Restoration**:
   - Signals: HRD drop after PARP, RAD51C/D mutations, loss of BRCA1/2 LOH, SAE DNA repair capacity >0.7
   - Triggers: PARP inhibitors
2. **ABCB1 Upregulation**:
   - Signals: ABCB1 copy number gain (>4 copies), ABCB1 activating mutations
   - Triggers: Paclitaxel, doxorubicin, topotecan, P-gp substrates
3. **MAPK Activation**:
   - Signals: KRAS/NRAS/BRAF activating mutations, high MAPK pathway burden
   - Triggers: BRAF inhibitors, EGFR inhibitors
4. **PI3K Activation**:
   - Signals: PIK3CA mutations, PTEN loss, high PI3K pathway burden
   - Triggers: PI3K inhibitors, mTOR inhibitors
5. **SLFN11 Deficiency**:
   - Signals: SLFN11 low expression (IHC), SLFN11 mutations
   - Triggers: PARP inhibitors, topoisomerase inhibitors

#### **6.5.3 Combo Strategies**:
- **PARP + ATR**: For HRD-high with HR restoration risk
- **PARP + VEGF**: For HRD-high with angiogenesis escape
- **IO Combos**: For TMB-high/MSI-high with IO resistance
- **MAPK/PI3K Combos**: For pathway escape mechanisms

---

### **6.6 TRIAL MATCHING SYSTEM**

#### **6.6.1 Architecture** (`services/ayesha_trial_matching/`):
- **MatchOrchestrator**: Coordinates workflow
- **EligibilityFilters**: Hard filters (MUST-MATCH)
- **ScoringEngine**: Soft boosts (scoring)
- **ReasoningGenerator**: Transparent "why-matched" explanations

#### **6.6.2 Eligibility Filters** (Hard Filters):
1. **Disease**: Ovarian/peritoneal/gynecologic
2. **Stage**: IV/advanced/metastatic
3. **Treatment Line**: First-line/untreated
4. **Status**: Recruiting/Active
5. **Location**: NYC metro (NY/NJ/CT)
6. **Exclusions**: NOT recurrent-only, NOT germline-BRCA-required

#### **6.6.3 Scoring Engine** (Soft Boosts):
- **Base Score**: 0.5
- **Boosts**:
  - First-line trial: +0.30
  - Stage IV specific: +0.25
  - SOC backbone (carboplatin + paclitaxel): +0.20
  - Germline-negative friendly: +0.20
  - IP chemotherapy: +0.20
  - Bevacizumab: +0.15
  - CA-125 tracking: +0.15
  - NYC location: +0.15
  - Large trial: +0.10
  - Phase III: +0.10
- **Penalties**:
  - Germline BRCA required: -0.30
  - Distance: -0.25
  - Phase I: -0.20

#### **6.6.4 Reasoning Generator**:
- Generates transparent explanations for each trial match
- Explains why trial matched (eligibility + scoring boosts)
- Includes clinical rationale

---

### **6.7 SOC RECOMMENDATION**

#### **6.7.1 Purpose** (`routers/ayesha_trials.py:_generate_soc_recommendation`):
- **NCCN-aligned** standard of care for Stage IVB HGSOC
- **Confidence**: 95-100% (guideline-aligned, no predictions)

#### **6.7.2 Regimen**:
- **Base**: Carboplatin AUC 5-6 + Paclitaxel 175 mg/m²
- **Add-on** (if ascites/peritoneal disease):
  - Bevacizumab 15 mg/kg
  - Rationale: Reduces progression risk (GOG-218 HR 0.72, ICON7 HR 0.81)

#### **6.7.3 Schedule**:
- **Induction**: 6 cycles every 21 days
- **Maintenance**: Bevacizumab continuation up to 15 months total OR progression
- **Typical Duration**: ~18 weeks induction + up to 12 months maintenance

#### **6.7.4 Monitoring Protocol**:
- **CA-125**: Every cycle (every 3 weeks)
- **Imaging**: Every 3 cycles
- **Toxicity**: Grade 3-4 neuropathy, cytopenias → dose reduction/delay

---

### **6.8 FOOD VALIDATOR**

#### **6.8.1 Purpose** (`routers/hypothesis_validator.py`):
- **Dynamic Food Validator**: Works for ANY food/supplement
- **Endpoint**: `POST /api/hypothesis/validate_food_dynamic`
- **A→B Dependency**: NO TUMOR NGS REQUIRED (works pre-NGS)

#### **6.8.2 Features**:
- **LLM Literature Mining**: PubMed paper reading
- **SAE Features**: Line appropriateness, cross-resistance, sequencing fitness
- **Complete Provenance**: Data sources, confidence breakdown

#### **6.8.3 Treatment Line Service** (`services/food_treatment_line_service.py`):
- Computes SAE features for dietary supplements:
  - `line_appropriateness`: 0-1 score
  - `cross_resistance`: 0-1 score
  - `sequencing_fitness`: 0-1 score
- **Biomarker Gates**: Boosts if biomarker context matches
- **Treatment History Gates**: Adjusts based on prior therapies

---

### **6.9 NEXT-TEST RECOMMENDER**

#### **6.9.1 Purpose** (`services/next_test_recommender.py`):
- **Prioritized biomarker testing** recommendations
- **Manager Policy**: P1, C6 (MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md)

#### **6.9.2 Priority Order**:
1. **HRD** (MyChoice CDx): PARP gate, 10d turnaround
2. **ctDNA** (Guardant360/FoundationOne): MSI/TMB + somatic HRR, 7d turnaround
3. **SLFN11 IHC**: PARP sensitivity, 5d turnaround
4. **ABCB1 proxy**: Efflux resistance (only if prior taxane), 5d turnaround

#### **6.9.3 Format**:
- **Differential branches**: "If positive → X; If negative → Y"
- **Turnaround + Cost**: Included in recommendations
- **Max 4 recommendations**

---

### **6.10 HINT TILES SERVICE**

#### **6.10.1 Purpose** (`services/hint_tiles_service.py`):
- **Clinician action hints** (max 4 tiles)
- **Priority Order** (Manager C8): Test → Trials → Monitor → Avoid
- **Suggestive tone** (not prescriptive)

#### **6.10.2 Tile Types**:
1. **Next Test**: Prioritized biomarker testing
2. **Trial Matched**: Clinical trial opportunities
3. **Monitoring**: CA-125, imaging schedules
4. **Avoid**: Drug interactions, contraindications

---

### **6.11 MECHANISM MAP SERVICE**

#### **6.11.1 Purpose** (`services/mechanism_map_service.py`):
- **Pathway burden visualization** (6 chips)
- **Pre-NGS**: All gray (no data)
- **Post-NGS**: Color-coded by pathway activation (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)

---

### **6.12 KEY INSIGHTS**

#### **Clinical Workflow**:
1. **Pre-NGS**: Phase 1 SAE services (next-test, hint tiles, mechanism map)
2. **Post-NGS**: Phase 2 SAE services (SAE features, resistance alert)
3. **Opt-in**: Resistance Prophet (predicts resistance 3-6 months early)

#### **Resistance Detection**:
1. **2-of-3 Trigger Rule**: HRD drop OR DNA repair drop OR CA-125 inadequate
2. **HR Restoration Pattern**: Immediate alert (don't wait for radiology)
3. **Resistance Prophet**: Early warning (3-6 months before clinical progression)

#### **Trial Matching**:
1. **Hard Filters**: Disease, stage, line, status, location (MUST-MATCH)
2. **Soft Boosts**: First-line, SOC backbone, CA-125 tracking, NYC location
3. **Transparent Reasoning**: Explains why each trial matched

#### **CA-125 Intelligence**:
1. **Burden Classification**: MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE
2. **Response Forecast**: Cycle 3 (≥70% drop), Cycle 6 (≥90% drop)
3. **Resistance Signals**: On-therapy rise, inadequate response (<50% drop by cycle 3)

#### **SOC Recommendation**:
1. **NCCN-aligned**: Carboplatin + Paclitaxel (standard)
2. **Bevacizumab Add-on**: If ascites/peritoneal disease
3. **Confidence**: 95-100% (guideline-based, not prediction)

---

**Status**: ✅ **ITERATION 6 COMPLETE** - Clinical Systems & Workflows fully documented  
**Next**: I7 - Research & Design Systems (Metastasis Interception, CRISPR Design, Protein Synthesis)

---

## 🔍 GAP ANALYSIS & MISSING SYSTEMS REVIEW

### **7.1 IDENTIFIED GAPS**

#### **7.1.1 Research & Design Systems** (Partially Documented):
- ✅ **Metastasis Interception**: Mentioned but needs deep dive
- ✅ **CRISPR Design Router**: Mentioned but needs deep dive
- ⚠️ **Datasets Router**: Cohort intelligence (cBioPortal, GDC) - NOT documented
- ⚠️ **Evidence RAG**: Evidence retrieval augmented generation - NOT documented
- ⚠️ **Evidence Extraction/Jobs**: Evidence processing pipeline - NOT documented

#### **7.1.2 Clinical Systems** (Partially Documented):
- ✅ **Clinical Genomics Command Center**: Frontend documented, backend router NOT documented
- ⚠️ **ACMG Router**: Variant classification - NOT documented
- ⚠️ **PharmGKB Router**: Pharmacogenomics - NOT documented
- ⚠️ **NCCN Router**: Guideline compliance - NOT documented
- ⚠️ **Safety Router**: Safety validation - NOT documented
- ⚠️ **Toxicity Router**: Toxicity prevention - NOT documented

#### **7.1.3 Platform Systems** (NOT Documented):
- ⚠️ **Agents System**: Agent management, execution, scheduling - NOT documented
- ⚠️ **Knowledge Base (KB) System**: KB endpoints, validation, storage - NOT documented
- ⚠️ **Knowledge Graph (KG) Router**: KG context - NOT documented
- ⚠️ **Client Dossier System**: Trial dossier generation - NOT documented
- ⚠️ **Trial Intelligence Pipeline**: Multi-stage trial processing - NOT documented

#### **7.1.4 Services** (Partially Documented):
- ⚠️ **Therapeutic Optimizer**: Iterative optimization - NOT documented
- ⚠️ **Safety Service**: Safety validation - NOT documented
- ⚠️ **Trial Intelligence Services**: Multi-stage pipeline - NOT documented
- ⚠️ **Job Service**: Background job processing - NOT documented
- ⚠️ **Logging Service**: Efficacy/evidence logging - NOT documented

---

### **7.2 RESEARCH & DESIGN SYSTEMS** (I7 - IN PROGRESS)

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

#### **7.2.4 Evidence RAG** (`routers/evidence/rag.py`):

**Purpose**: Evidence retrieval augmented generation for literature search

**Status**: ⚠️ **NEEDS EXPLORATION** - File exists but not yet documented

---

#### **7.2.5 Evidence Extraction/Jobs** (`routers/evidence/extraction.py`, `routers/evidence/jobs.py`):

**Purpose**: Evidence processing pipeline for background job processing

**Status**: ⚠️ **NEEDS EXPLORATION** - Files exist but not yet documented

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

#### **7.9.1 Therapeutic Optimizer** (`services/therapeutic_optimizer.py`):
- **Purpose**: Iterative optimization loops for therapeutic design
- **Status**: ⚠️ **NEEDS EXPLORATION**

#### **7.9.2 Safety Service** (`services/safety_service.py`):
- **Purpose**: Safety validation for therapeutic designs
- **Status**: ⚠️ **NEEDS EXPLORATION**

#### **7.9.3 Job Service** (`services/job_service.py`):
- **Purpose**: Background job processing
- **Status**: ⚠️ **NEEDS EXPLORATION**

#### **7.9.4 Logging Service** (`services/logging_service.py`):
- **Purpose**: Efficacy/evidence logging
- **Status**: ⚠️ **NEEDS EXPLORATION**

---

**Status**: 🔄 **GAP ANALYSIS COMPLETE** - Major systems identified and partially documented  
**Next**: Complete I7 (Research & Design Systems) with remaining systems

---

## 🔍 COMPREHENSIVE GAP ANALYSIS - ITERATION 7+ REVIEW

### **7.10 ADDITIONAL ROUTERS NOT FULLY DOCUMENTED**

#### **7.10.1 Guidance Router** (`routers/guidance.py`):
- **Purpose**: Clinical gating facade over efficacy orchestrator
- **Features**:
  - On-label detection (stub ruleset)
  - Tier classification (I/II/III/research)
  - Resistance/sensitivity marker detection (PSMB5, CRBN, TP53)
  - Fused S score integration (AlphaMissense)
- **Status**: ⚠️ **NOT DOCUMENTED** - Clinical gating logic exists

#### **7.10.2 Sessions Router** (`routers/sessions.py`):
- **Purpose**: Session persistence API - save/resume analyses across pages
- **Features**:
  - Create/update sessions (Supabase-backed)
  - Add session items (insight/efficacy/dataset/note)
  - Retrieve session history
  - Idempotency support
- **Status**: ⚠️ **NOT DOCUMENTED** - Session management system

#### **7.10.3 Auth Router** (`routers/auth.py`):
- **Purpose**: User authentication and profile management
- **Endpoints**:
  - `/api/auth/signup` - Create account
  - `/api/auth/login` - Login
  - `/api/auth/logout` - Logout
  - `/api/auth/profile` - Get/update profile
- **Status**: ⚠️ **NOT DOCUMENTED** - Authentication system

#### **7.10.4 Admin Router** (`routers/admin.py`):
- **Purpose**: Admin dashboard endpoints
- **Features**:
  - User management (list, get, update, delete)
  - Analytics aggregation
  - Activity tracking
  - Usage limits management
- **Status**: ⚠️ **NOT DOCUMENTED** - Admin functionality

#### **7.10.5 Myeloma Router** (`routers/myeloma.py`):
- **Purpose**: Myeloma Digital Twin prediction
- **Endpoint**: `/api/predict/myeloma_drug_response`
- **Features**:
  - Evo2 live scoring only (no mocks)
  - Preflight validation (REF-check, duplicate collapse)
  - Variant call extraction
  - Supabase persistence
- **Status**: ⚠️ **NOT DOCUMENTED** - Myeloma-specific workflow

#### **7.10.6 Evidence RAG** (`routers/evidence/rag.py`):
- **Purpose**: RAG-based conversational query for clinical literature
- **Endpoints**:
  - `/api/evidence/rag-query` - Conversational query
  - `/api/evidence/rag-add-variant` - Add variant to KB
  - `/api/evidence/rag-stats` - KB statistics
- **Integration**: Uses external `Pubmed-LLM-Agent-main` RAG agent
- **Status**: ⚠️ **NOT DOCUMENTED** - RAG system

#### **7.10.7 Evidence Extraction** (`routers/evidence/extraction.py`):
- **Purpose**: Diffbot article extraction
- **Endpoint**: `/api/evidence/extract`
- **Features**: Extract full text from URLs via Diffbot API
- **Status**: ⚠️ **NOT DOCUMENTED** - Article extraction

#### **7.10.8 Evidence Jobs** (`routers/evidence/jobs.py`):
- **Purpose**: Background job orchestration for evidence processing
- **Endpoints**:
  - `/api/evidence/crawl` - Multi-URL crawling job
  - `/api/evidence/summarize` - Summarization job
  - `/api/evidence/align` - Evidence alignment job
- **Status**: ⚠️ **NOT DOCUMENTED** - Background job system

#### **7.10.9 Resistance Router** (`routers/resistance.py`):
- **Purpose**: Resistance analysis endpoints
- **Status**: ⚠️ **NEEDS EXPLORATION**

#### **7.10.10 Care Router** (`routers/care.py`):
- **Purpose**: Resistance Playbook endpoints (Section 17)
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.10.11 Tumor Router** (`routers/tumor.py`):
- **Purpose**: Sporadic Cancer Strategy endpoints (Day 1-7)
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.10.12 Ayesha Router** (`routers/ayesha.py`):
- **Purpose**: Original Ayesha endpoints (not v2 orchestrator)
- **Status**: ⚠️ **NOT DOCUMENTED** - Legacy Ayesha endpoints

#### **7.10.13 Ayesha Twin Demo Router** (`routers/ayesha_twin_demo.py`):
- **Purpose**: Ayesha twin demo endpoints
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.14 Dossiers Router** (`routers/dossiers.py`):
- **Purpose**: JR2 dossier generation pipeline
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.10.15 Ayesha Dossiers Router** (`routers/ayesha_dossiers.py`):
- **Purpose**: Ayesha dossier browser API (display all 60 trials)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.16 Trials Router** (`routers/trials.py`):
- **Purpose**: Search and refresh endpoints
- **Endpoints**: `/api/search-trials`, `/api/trials/refresh_status`
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.17 Trials Graph Router** (`routers/trials_graph.py`):
- **Purpose**: Graph-optimized search (hybrid graph search)
- **Endpoint**: `/api/trials/search-optimized`
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.18 Trials Agent Router** (`routers/trials_agent.py`):
- **Purpose**: Autonomous trial agent
- **Endpoint**: `/api/trials/agent/search`
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.19 Clinical Trials Router** (`routers/clinical_trials.py`):
- **Purpose**: Clinical trials endpoints
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.20 Command Center Router** (`routers/command_center.py`):
- **Purpose**: Command center endpoints (conditional on feature flag)
- **Status**: ⚠️ **NOT DOCUMENTED**

---

### **7.11 ADDITIONAL SERVICES NOT FULLY DOCUMENTED**

#### **7.11.1 Admin Service** (`services/admin_service.py`):
- **Purpose**: Admin operations (user management, analytics)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.2 Auth Service** (`services/auth_service.py`):
- **Purpose**: Authentication operations (Supabase Auth integration)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.3 Agent Services**:
- **Agent Executor** (`services/agent_executor.py`): Execute agent tasks
- **Agent Manager** (`services/agent_manager.py`): CRUD for agents
- **Agent Scheduler** (`services/agent_scheduler.py`): Scheduled execution
- **Status**: ⚠️ **PARTIALLY DOCUMENTED** - Router documented, services not

#### **7.11.4 Autonomous Trial Agent** (`services/autonomous_trial_agent.py`):
- **Purpose**: Autonomous trial search agent
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.5 Cache Service** (`services/cache_service.py`):
- **Purpose**: Caching operations
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.6 Compound Services**:
- **Compound Alias Resolver** (`services/compound_alias_resolver.py`): Resolve drug aliases
- **Compound Calibration** (`services/compound_calibration.py`): Drug-specific calibration
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.7 Dietician Recommendations** (`services/dietician_recommendations.py`):
- **Purpose**: Dietary recommendations
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.8 Dynamic Food Extraction** (`services/dynamic_food_extraction.py`):
- **Purpose**: Extract food/supplement information
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.11.9 Food S/P/E Integration** (`services/food_spe_integration.py`):
- **Purpose**: S/P/E framework for foods
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.10 Hotspot Detector** (`services/hotspot_detector.py`):
- **Purpose**: Detect mutation hotspots
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.11 LLM Literature Service** (`services/llm_literature_service.py`):
- **Purpose**: LLM-based literature synthesis
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.12 Mechanism Fit Ranker** (`services/mechanism_fit_ranker.py`):
- **Purpose**: Rank drugs by mechanism fit
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.13 Metastasis Service** (`services/metastasis_service.py`):
- **Purpose**: Metastasis cascade risk assessment (NOT interception)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.14 NGS Fast Track** (`services/ngs_fast_track.py`):
- **Purpose**: Fast-track NGS recommendations
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.11.15 Safety Validator** (`services/safety_validator.py`):
- **Purpose**: Safety validation logic
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.16 Therapeutic Prompt Builder** (`services/therapeutic_prompt_builder.py`):
- **Purpose**: Build prompts for therapeutic generation
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.17 Trial Intelligence Pipeline** (`services/trial_intelligence/`):
- **Purpose**: Multi-stage trial processing pipeline
- **Stages**: 6 stages (hard filters → trial type → location → eligibility → LLM analysis → dossier)
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.11.18 Trial Refresh Service** (`services/trial_refresh/`):
- **Purpose**: Refresh trial data from ClinicalTrials.gov
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.19 Tumor Quick Intake** (`services/tumor_quick_intake.py`):
- **Purpose**: Generate initial tumor context using disease priors
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.11.20 Logging Services** (`services/logging/`):
- **Efficacy Logger**: Log efficacy predictions
- **Evidence Logger**: Log evidence queries
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.21 Infrastructure Services**:
- **Supabase Service** (`services/supabase_service.py`): Supabase operations
- **Neo4j Connection** (`services/neo4j_connection.py`): Neo4j graph database
- **Neo4j Graph Loader** (`services/neo4j_graph_loader.py`): Load graph data
- **Database Connections** (`services/database_connections.py`): DB connection management
- **Enformer Client** (`services/enformer_client.py`): Chromatin accessibility client
- **Status**: ⚠️ **NOT DOCUMENTED**

---

### **7.12 FRONTEND PAGES NOT FULLY DOCUMENTED**

#### **7.12.1 Auth Pages**:
- **Login** (`pages/auth/Login.jsx`): Login page
- **Signup** (`pages/auth/Signup.jsx`): Signup page
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.12.2 Admin Pages**:
- **Admin Dashboard** (`pages/admin/Dashboard.jsx`): Admin dashboard
- **Admin Users** (`pages/admin/Users.jsx`): User management
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.12.3 Research Pages**:
- **Research** (`pages/Research.jsx`): Research page
- **Mutation Explorer** (`pages/MutationExplorer.jsx`): VUS Explorer
- **Genomic Analysis** (`pages/GenomicAnalysis.tsx`): Genomic analysis
- **Status**: ⚠️ **PARTIALLY DOCUMENTED** - Mentioned but not deep dive

#### **7.12.4 Agent Pages**:
- **Agent Dashboard** (`pages/AgentDashboard.jsx`): Agent dashboard
- **Agent Demo** (`pages/AgentDemo.jsx`): Agent demo
- **Agent Studio** (`pages/AgentStudio.jsx`): Agent creation
- **Agents Page** (`pages/AgentsPage.jsx`): Agent list
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.12.5 Clinical Pages**:
- **Threat Assessor** (`pages/ThreatAssessor.jsx`): Threat assessment
- **RadOnc CoPilot** (`pages/RadOncCoPilot.jsx`): Radiation oncology
- **Myeloma Digital Twin** (`pages/MyelomaDigitalTwin.jsx`): Myeloma twin
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.12.6 Design Pages**:
- **Armory** (`pages/Armory.jsx`): Tools page
- **Metastasis Dashboard** (`pages/MetastasisDashboard.jsx`): Metastasis tools
- **Crispr Designer** (`pages/CrisprDesigner.jsx`): CRISPR design
- **Protein Synthesis** (`pages/ProteinSynthesis.jsx`): Protein synthesis
- **Structure Predictor** (`pages/StructurePredictor.jsx`): Structure prediction
- **Status**: ⚠️ **PARTIALLY DOCUMENTED** - Mentioned but not deep dive

#### **7.12.7 Campaign Pages**:
- **Campaign Runner** (`pages/CampaignRunner.jsx`): Campaign execution
- **Runx Conquest** (`pages/RunxConquest.jsx`): Runx demo
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.12.8 Other Pages**:
- **Target Dossier** (`pages/TargetDossier.jsx`): Dossier generation
- **Demo Summarizer** (`pages/DemoSummarizer.jsx`): Demo tool
- **Investor Slideshow** (`pages/InvestorSlideshow.jsx`): Investor deck
- **Patient Tasks Page** (`pages/PatientTasksPage.jsx`): Patient tasks
- **Follow Up Task Board** (`pages/FollowUpTaskBoard.jsx`): Task board
- **Screening Schedule** (`pages/ScreeningSchedule.jsx`): Screening
- **Medical Records** (`pages/records/index.jsx`): Records list
- **Single Record Details** (`pages/records/single-record-details.jsx`): Record detail
- **Status**: ⚠️ **NOT DOCUMENTED**

---

### **7.13 INTEGRATION PATTERNS NOT DOCUMENTED**

#### **7.13.1 Supabase Integration**:
- **Tables**: `user_sessions`, `mdt_runs`, `mdt_run_variants`, `job_results`
- **Operations**: Select, insert, update, events
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.13.2 Neo4j Graph Integration**:
- **Graph Structure**: Trial relationships, eligibility graph
- **Queries**: Graph-optimized trial search
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.13.3 AstraDB Integration**:
- **Collections**: `clinical_trials_eligibility` (with vector search)
- **Operations**: Vector search, document storage
- **Status**: ⚠️ **PARTIALLY DOCUMENTED** - Vector search issue documented

#### **7.13.4 SQLite Integration**:
- **Tables**: `clinical_trials` (28 columns)
- **Operations**: Trial storage, querying
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.13.5 Background Job System**:
- **Job Types**: Crawl, summarize, align
- **Job Store**: In-memory (JOBS dict) or Supabase
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.13.6 Agent System Architecture**:
- **Agent Types**: `pubmed_sentinel`, `trial_scout`
- **Execution**: Scheduled vs on-demand
- **Status**: ⚠️ **PARTIALLY DOCUMENTED** - Router documented, architecture not

---

### **7.14 DATA FLOW PATTERNS NOT DOCUMENTED**

#### **7.14.1 End-to-End Workflows**:
- **Sporadic Cancer Workflow**: Pre-NGS → Post-NGS → Resistance detection
- **Trial Matching Workflow**: Search → Filter → Score → Rank → Dossier
- **Efficacy Workflow**: Sequence → Pathway → Evidence → Drug Scoring
- **Status**: ⚠️ **PARTIALLY DOCUMENTED** - Components documented, full flow not

#### **7.14.2 Service Communication Patterns**:
- **HTTP vs Direct Import**: When to use HTTP vs direct service import
- **Caching Strategies**: Redis vs in-memory vs Supabase
- **Error Propagation**: How errors flow through orchestrators
- **Status**: ⚠️ **NOT DOCUMENTED**

---

### **7.15 SUMMARY OF GAPS**

#### **Critical Gaps** (High Priority):
1. **Authentication & Authorization**: Auth service, admin service, protected routes
2. **Session Management**: Session persistence, cross-page state
3. **Agent System**: Full architecture, execution patterns, scheduling
4. **Trial Intelligence Pipeline**: 6-stage pipeline details
5. **Evidence RAG**: RAG agent integration, knowledge base management
6. **Background Jobs**: Job orchestration, status tracking
7. **Database Integrations**: Supabase, Neo4j, AstraDB, SQLite patterns

#### **Important Gaps** (Medium Priority):
1. **Myeloma Digital Twin**: Myeloma-specific workflow
2. **Guidance Router**: Clinical gating logic
3. **Resistance Router**: Resistance analysis endpoints
4. **Trial Routers**: Multiple trial search endpoints
5. **Food Services**: Food extraction, S/P/E integration
6. **Safety Validator**: Safety validation logic
7. **Frontend Pages**: Many pages not documented

#### **Nice-to-Have Gaps** (Low Priority):
1. **Campaign System**: Campaign runner, Runx conquest
2. **Demo Tools**: Demo summarizer, investor slideshow
3. **Task Management**: Patient tasks, follow-up board
4. **Medical Records**: Records management

---

**Status**: ✅ **COMPREHENSIVE GAP ANALYSIS COMPLETE** - All major gaps identified  
**Next**: Prioritize gap closure based on criticality and user needs

