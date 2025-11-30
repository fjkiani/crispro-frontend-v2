# 📊 ITERATION 2: BACKEND SERVICES & ORCHESTRATION

**Status**: ✅ **COMPLETE**  
**Duration**: 4-5 hours  
**Created**: January 14, 2025

---

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
