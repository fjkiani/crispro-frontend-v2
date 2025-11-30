# 🔄 ITERATION 8: DATA FLOW & INTEGRATION PATTERNS

**Status**: ✅ **COMPLETE**  
**Duration**: 2-3 hours  
**Created**: January 14, 2025  
**Cycle**: I8 (Cycle 4)

---

## **8.1 COMPLETE USER FLOW: WIWFM (Will It Work For Me)**

### **8.1.1 Frontend → Backend Flow**

#### **Step 1: User Input** (`MechanisticEvidenceTab.jsx`)
```javascript
// User clicks "Run Deep Analysis"
handleDeepAnalysis() {
  await predict(
    patientProfile.mutations,  // [{gene, chrom, pos, ref, alt, ...}]
    patientProfile.cancer_type,
    profile  // "baseline" | "richer" | "fusion"
  );
}
```

#### **Step 2: API Client** (`useApiClient.js`)
```javascript
// Custom hook for API communication
const client = useApiClient(modelId);
await client.post('/api/clinical_genomics/analyze_variant', {
  mutations: [...],
  disease: "ovarian_carcinoma",
  profile: "baseline",
  include: ["sae_features"]
});
```

**Features**:
- **Timeout**: 10 minutes default (`DEFAULT_TIMEOUT_MS`)
- **AbortController**: Cancels request on timeout
- **Error Extraction**: Extracts `detail` from error responses

#### **Step 3: Backend Router** (`clinical_genomics.py`)
```python
@router.post("/analyze_variant")
async def analyze_variant(request: AnalyzeVariantRequest):
    # Generate unified run_id
    run_id = str(uuid.uuid4())
    
    # Call efficacy orchestrator directly (not HTTP)
    orchestrator = create_efficacy_orchestrator(api_base="http://127.0.0.1:8000")
    efficacy_request = EfficacyRequest(
        mutations=request.mutations,
        model_id="evo2_1b",
        options={
            "profile": request.profile,
            "fast": is_baseline,  # Skip evidence for baseline
            "include_sae_features": True
        }
    )
    efficacy_response = await orchestrator.predict(efficacy_request)
```

**Key Decision**: Direct orchestrator call (not HTTP) to reduce latency

---

### **8.1.2 Backend → AI Services Flow**

#### **Step 4: Efficacy Orchestrator** (`orchestrator.py`)
```python
async def predict(self, request: EfficacyRequest):
    # 1. Sequence Scoring (Evo2/Fusion)
    seq_scores = await sequence_processor.process(
        mutations=request.mutations,
        model_id=request.model_id,
        options=request.options
    )
    
    # 2. Pathway Aggregation
    pathway_scores = await pathway_service.aggregate(
        mutations=request.mutations,
        panel=panel
    )
    
    # 3. Evidence Gathering (if not fast mode)
    if not request.options.get("fast"):
        evidence_scores = await evidence_service.gather(
            mutations=request.mutations,
            disease=request.disease
        )
    
    # 4. Drug Scoring (S/P/E formula)
    drugs = await drug_scorer.score(
        seq_scores=seq_scores,
        pathway_scores=pathway_scores,
        evidence_scores=evidence_scores
    )
```

#### **Step 5: Evo2 Service Call** (`evo2_scorer.py`)
```python
async def score(self, mutation, model_id, options):
    # Check cache first
    cache_key = self._get_cache_key(mutation, model_id, window_flanks, ensemble)
    cached_result = await get_cache(cache_key)
    if cached_result:
        return cached_result
    
    # Call Evo2 via HTTP
    base_url = get_model_url(model_id)  # Modal service URL
    async with httpx.AsyncClient(timeout=60.0) as client:
        response = await client.post(
            f"{base_url}/score_variant_multi",
            json={
                "chrom": mutation["chrom"],
                "pos": mutation["pos"],
                "ref": mutation["ref"],
                "alt": mutation["alt"],
                "model_id": model_id
            }
        )
        result = response.json()
    
    # Cache result
    await set_cache(cache_key, result, ttl=3600)
    return result
```

**Modal Service URL Pattern**:
- **Format**: `https://crispro--evo-service-evoservice{model}-api-{model}.modal.run`
- **1B**: `https://crispro--evo-service-evoservice1b-api-1b.modal.run`
- **7B**: `https://crispro--evo-service-evoservice7b-api-7b.modal.run`
- **40B**: `https://crispro--evo-service-evoservice-api.modal.run`

**Fallback Logic** (`config.py:146-161`):
```python
def get_model_url(model_id: str) -> str:
    if model_id == "evo2_1b":
        return url_1b or url_7b or url_40b  # Fallback chain
    elif model_id == "evo2_7b":
        return url_7b or url_40b
    elif model_id == "evo2_40b":
        return url_40b
```

---

### **8.1.3 External API Integration**

#### **Step 6: PubMed Literature** (`literature_client.py`)
```python
async def search_literature(gene, hgvs_p, disease, max_results=10):
    # Check cache first
    cache_key = literature_cache_key(gene, hgvs_p, disease, max_results)
    cached = await get_cache(cache_key)
    if cached:
        return cached
    
    # Call PubMed E-utils API
    async with httpx.AsyncClient(timeout=60.0) as client:
        # E-search
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            "db": "pubmed",
            "term": f"{gene}[Gene] AND {disease}[MeSH Terms]",
            "retmode": "json",
            "retmax": max_results
        }
        search_response = await client.get(search_url, params=search_params)
        pmids = search_response.json()["esearchresult"]["idlist"]
        
        # E-fetch
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml"
        }
        fetch_response = await client.get(fetch_url, params=fetch_params)
        articles = parse_pubmed_xml(fetch_response.text)
    
    # Cache result
    await set_cache(cache_key, articles, ttl=LITERATURE_TTL)
    return articles
```

#### **Step 7: ClinVar Prior** (`clinvar_client.py`)
```python
async def get_clinvar_prior(mutation):
    # Check cache first
    cache_key = f"clinvar:{mutation['gene']}:{mutation.get('hgvs_p', '')}"
    cached = await get_cache(cache_key)
    if cached:
        return cached
    
    # Call ClinVar API (or local database)
    # ... ClinVar lookup logic ...
    
    # Cache result (24 hour TTL)
    await set_cache(cache_key, prior, ttl=86400)
    return prior
```

---

## **8.2 DATA TRANSFORMATION AT EACH LAYER**

### **8.2.1 Frontend → Backend**

**Input** (React):
```javascript
{
  mutations: [
    {
      gene: "BRAF",
      chrom: "7",
      pos: 140753336,
      ref: "T",
      alt: "A",
      build: "GRCh38",
      consequence: "missense_variant"
    }
  ],
  disease: "ovarian_carcinoma",
  profile: "baseline"
}
```

**Transformation**: JSON serialization → HTTP POST → FastAPI Pydantic model

**Output** (FastAPI):
```python
AnalyzeVariantRequest(
    mutations=[...],  # List[Dict[str, Any]]
    disease="ovarian_carcinoma",
    profile="baseline"
)
```

---

### **8.2.2 Backend → Evo2 Service**

**Input** (Backend):
```python
{
    "chrom": "7",
    "pos": 140753336,
    "ref": "T",
    "alt": "A",
    "assembly": "GRCh38",
    "model_id": "evo2_1b"
}
```

**Transformation**: HTTP POST → Modal service → Evo2 model inference

**Output** (Evo2):
```python
{
    "min_delta": -2.34,
    "window_used": 8192,
    "deltas": [-2.34, -2.12, -1.98]
}
```

**Backend Processing**:
```python
# Convert to SeqScore
seq_score = SeqScore(
    mutation=mutation,
    delta_score=-2.34,
    sequence_disruption=0.85,  # Computed from delta
    calibrated_seq_percentile=0.92  # Percentile mapping
)
```

---

### **8.2.3 Sequence → Pathway Aggregation**

**Input** (Sequence Scores):
```python
[
    SeqScore(gene="BRAF", delta_score=-2.34, ...),
    SeqScore(gene="KRAS", delta_score=-1.89, ...)
]
```

**Transformation**: Gene-to-pathway mapping → Pathway aggregation

**Pathway Mapping** (`drug_mapping.py`):
```python
GENE_TO_PATHWAY = {
    "BRAF": {"ras_mapk": 1.0},
    "KRAS": {"ras_mapk": 1.0},
    "NRAS": {"ras_mapk": 1.0}
}
```

**Aggregation** (`aggregation.py`):
```python
pathway_scores = {
    "ras_mapk": 0.85,  # Average of BRAF + KRAS scores
    "pi3k_akt": 0.42,
    "dna_repair": 0.78
}
```

---

### **8.2.4 Pathway → Drug Scoring**

**Input** (Pathway Scores):
```python
{
    "ras_mapk": 0.85,
    "pi3k_akt": 0.42
}
```

**Transformation**: Drug-to-pathway alignment → Drug scoring

**Drug Mapping** (`panel_config.py`):
```python
DRUG_TO_PATHWAY = {
    "trametinib": {"ras_mapk": 0.9},  # MEK inhibitor
    "dabrafenib": {"ras_mapk": 0.95},  # BRAF inhibitor
    "alpelisib": {"pi3k_akt": 0.85}  # PI3K inhibitor
}
```

**Drug Scoring** (`drug_scorer.py`):
```python
# S/P/E formula
efficacy_score = (
    0.3 * seq_pct +      # Sequence percentile
    0.4 * path_pct +    # Pathway alignment
    0.3 * evidence       # Evidence strength
) + clinvar_prior
```

---

## **8.3 INTEGRATION POINTS**

### **8.3.1 AI Services (Modal)**

| Service | Endpoint | Purpose | Fallback |
|---------|----------|---------|----------|
| **Evo2 1B** | `/score_variant_multi` | Fast sequence scoring | → 7B → 40B |
| **Evo2 7B** | `/score_variant_multi` | Higher accuracy | → 40B |
| **Evo2 40B** | `/score_variant_multi` | Highest accuracy | None |
| **Fusion Engine** | `/score` | AlphaMissense + ESM | None |
| **Zeta Oracle** | `/fused_score` | Multi-model fusion | None |
| **Boltz Service** | `/predict_structure` | Structure prediction | None |

**Communication Pattern**: HTTP POST with JSON payload

**Timeout**: 60-180 seconds (service-dependent)

---

### **8.3.2 External APIs**

| API | Endpoint | Purpose | Caching |
|-----|----------|---------|---------|
| **PubMed E-utils** | `esearch.fcgi` | Literature search | 24h TTL |
| **PubMed E-utils** | `efetch.fcgi` | Article retrieval | 24h TTL |
| **ClinVar** | Local DB/API | Variant classification | 24h TTL |
| **Ensembl** | `/sequence/region` | Reference sequences | No cache |
| **Diffbot** | `/article` | Full-text extraction | No cache |
| **Google Gemini** | `/generateContent` | LLM synthesis | No cache |

**Communication Pattern**: HTTP GET/POST with query params or JSON

**Timeout**: 20-60 seconds (API-dependent)

---

### **8.3.3 Database Integrations**

| Database | Purpose | Operations | Caching |
|----------|---------|------------|---------|
| **Supabase** | User sessions, job results | SELECT, INSERT, UPDATE | No cache |
| **AstraDB** | Clinical trials (vector search) | Vector search, document storage | No cache |
| **SQLite** | Clinical trials (local) | SELECT, INSERT | No cache |
| **Neo4j** | Trial relationships (graph) | Graph queries | No cache |
| **Redis** | Caching layer | GET, SET, DELETE | N/A |

**Communication Pattern**: 
- **Supabase**: HTTP REST API
- **AstraDB**: HTTP REST API (vector search)
- **SQLite**: Direct file access
- **Neo4j**: Bolt protocol (graph queries)
- **Redis**: Redis protocol (caching)

---

## **8.4 CACHING STRATEGY**

### **8.4.1 Caching Layers**

#### **Layer 1: Redis (Distributed)**
- **Purpose**: Shared cache across backend instances
- **TTL**: Service-specific (1h-24h)
- **Key Pattern**: `{service}:{identifier}:{profile}`
- **Fallback**: Memory cache if Redis unavailable

#### **Layer 2: In-Memory (Process-Local)**
- **Purpose**: Fast local cache (no network)
- **TTL**: Same as Redis
- **Key Pattern**: Same as Redis
- **Scope**: Single backend process

#### **Layer 3: Frontend (localStorage)**
- **Purpose**: Client-side caching (not implemented yet)
- **TTL**: User-controlled
- **Key Pattern**: `{endpoint}:{payload_hash}`

---

### **8.4.2 Cache Keys by Service**

| Service | Cache Key Pattern | TTL | Example |
|---------|------------------|-----|---------|
| **Evo2** | `evo2:{model_id}:{hgvs_p}:{hash(window_flanks)}` | 1h | `evo2:evo2_1b:V600E:abc123` |
| **Fusion** | `fusion_am:{hash(tuple(mutation_keys))}` | 1h | `fusion_am:def456` |
| **Insights** | `insights:{variant_key}:{profile}` | 24h | `insights:BRAF_V600E:baseline` |
| **Literature** | `literature:{gene}:{hgvs_p}:{disease}:{max_results}` | 24h | `literature:BRAF:V600E:ovarian:10` |
| **ClinVar** | `clinvar:{gene}:{hgvs_p}` | 24h | `clinvar:BRAF:V600E` |
| **Evidence** | `evidence:{gene}:{disease}` | 30m | `evidence:BRAF:ovarian` |

---

### **8.4.3 Single-Flight Pattern**

**Purpose**: Prevent duplicate concurrent requests for the same resource

**Implementation** (`cache_service.py:59-105`):
```python
async def with_singleflight(key: str, ttl_lock: int, fn: Callable) -> Any:
    # 1. Check cache first
    cached = await get_cache(cache_key)
    if cached:
        return cached
    
    # 2. Try to acquire lock
    lock_acquired = await client.set(lock_key, "1", nx=True, ex=ttl_lock)
    
    if lock_acquired:
        # We got the lock, execute function
        result = await fn()
        await set_cache(cache_key, result, ttl=3600)
        return result
    else:
        # Someone else has the lock, wait and retry
        await asyncio.sleep(0.1)
        cached = await get_cache(cache_key)
        if cached:
            return cached
        # Last resort: execute function anyway
        return await fn()
```

**Use Case**: Multiple users request same variant → Only one Evo2 call, all get same result

---

## **8.5 ERROR PROPAGATION**

### **8.5.1 Error Flow**

```
External API Error (PubMed timeout)
  ↓
Backend Service (literature_client.py)
  ↓
Efficacy Orchestrator (orchestrator.py)
  ↓
Backend Router (clinical_genomics.py)
  ↓
Frontend (useApiClient.js)
  ↓
UI Component (MechanisticEvidenceTab.jsx)
```

### **8.5.2 Error Handling Patterns**

#### **1. Timeout Errors**
```python
# Backend
try:
    result = await asyncio.wait_for(
        call_external_api(),
        timeout=60.0
    )
except asyncio.TimeoutError:
    # Return partial result or default
    return {"evidence_strength": 0.0, "tier": "insufficient"}
```

#### **2. HTTP Errors**
```python
# Backend
try:
    response = await client.post(url, json=payload)
    response.raise_for_status()
except httpx.HTTPStatusError as e:
    # Log error, return fallback
    logger.error(f"HTTP error: {e}")
    return default_result
```

#### **3. Graceful Degradation**
```python
# Orchestrator
try:
    evidence_scores = await evidence_service.gather(...)
except Exception as e:
    logger.warning(f"Evidence gathering failed: {e}")
    # Continue with S/P only (no E)
    evidence_scores = None
```

#### **4. Frontend Error Display**
```javascript
// Frontend
try {
    const result = await client.post(endpoint, payload);
    setResult(result);
} catch (error) {
    // Extract error detail
    const message = error.message || error.detail || "Unknown error";
    setError(message);
    // Show user-friendly error
    showErrorToast(message);
}
```

---

## **8.6 DATA FLOW DIAGRAMS**

### **8.6.1 Complete WIWFM Flow**

```
User Input (React)
  ↓
useApiClient.js (HTTP POST)
  ↓
Backend Router (clinical_genomics.py)
  ↓
Efficacy Orchestrator (orchestrator.py)
  ├─→ Sequence Processor (evo2_scorer.py)
  │     ├─→ Cache Check (Redis/Memory)
  │     ├─→ Evo2 Modal Service (HTTP)
  │     └─→ Cache Store
  ├─→ Pathway Service (aggregation.py)
  │     └─→ Gene-to-Pathway Mapping
  ├─→ Evidence Service (literature_client.py)
  │     ├─→ Cache Check
  │     ├─→ PubMed E-utils API (HTTP)
  │     └─→ Cache Store
  └─→ Drug Scorer (drug_scorer.py)
        └─→ S/P/E Formula
  ↓
Backend Response (JSON)
  ↓
Frontend Display (React Components)
```

---

### **8.6.2 Caching Flow**

```
Request → Cache Check (Redis)
  ├─→ Hit: Return cached result
  └─→ Miss: 
        ├─→ Single-Flight Lock
        ├─→ Execute Function
        ├─→ Cache Result
        └─→ Return Result
```

---

### **8.6.3 Error Propagation Flow**

```
External API Error
  ↓
Service Layer (try/except)
  ├─→ Log Error
  ├─→ Return Fallback/Default
  └─→ Continue Processing
  ↓
Orchestrator (graceful degradation)
  ├─→ Partial Results OK
  └─→ Provenance Tracks Errors
  ↓
Frontend (error display)
  ├─→ Show User-Friendly Message
  └─→ Partial Results Displayed
```

---

## **8.7 SERVICE-TO-SERVICE COMMUNICATION**

### **8.7.1 HTTP vs Direct Import**

**HTTP Communication**:
- **When**: Cross-service boundaries (Modal, external APIs)
- **Why**: Decoupling, scalability, independent deployment
- **Examples**: Evo2 Modal, PubMed API, ClinVar API

**Direct Import**:
- **When**: Same backend process
- **Why**: Lower latency, type safety, shared state
- **Examples**: Orchestrator → Sequence Processor, Drug Scorer

---

### **8.7.2 Async vs Sync**

**Async (Preferred)**:
- **When**: I/O operations (API calls, database queries)
- **Why**: Non-blocking, better concurrency
- **Examples**: `async def predict()`, `await client.post()`

**Sync (Rare)**:
- **When**: CPU-bound operations (calculations, transformations)
- **Why**: Simpler, no async overhead
- **Examples**: `def aggregate_pathways()`, `def compute_percentile()`

---

## **8.8 SUMMARY**

### **Key Patterns**:

1. **Frontend → Backend**: HTTP POST with JSON
2. **Backend → AI Services**: HTTP POST to Modal URLs
3. **Backend → External APIs**: HTTP GET/POST with query params
4. **Backend Internal**: Direct function calls (no HTTP)
5. **Caching**: Redis → Memory fallback, single-flight protection
6. **Error Handling**: Graceful degradation, partial results OK
7. **Timeouts**: Service-specific (20s-180s)

### **Integration Points**:

- **Evo2**: Modal service URLs with fallback chain
- **PubMed**: E-utils API (esearch + efetch)
- **ClinVar**: Local DB or API
- **Supabase**: REST API for sessions/jobs
- **AstraDB**: Vector search for trials
- **Redis**: Caching layer

### **Data Transformations**:

- **Frontend → Backend**: JSON serialization
- **Backend → Evo2**: Variant coordinates → Delta scores
- **Sequence → Pathway**: Gene scores → Pathway aggregation
- **Pathway → Drug**: Pathway alignment → Efficacy scores

---

**Status**: ✅ **CYCLE 4 COMPLETE** - Data Flow & Integration Patterns  
**Next**: Cycle 5 (I9) - Development Patterns & Lessons Learned

---



