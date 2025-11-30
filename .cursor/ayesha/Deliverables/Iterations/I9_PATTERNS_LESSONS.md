# 🎯 ITERATION 9: DEVELOPMENT PATTERNS & LESSONS LEARNED

**Status**: ✅ **COMPLETE**  
**Duration**: 2-3 hours  
**Created**: January 14, 2025  
**Cycle**: I9 (Cycle 5)

---

## **9.1 CORE TECHNICAL DOCTRINES**

### **9.1.1 The "Wet Noodle" Doctrine**

**Problem**: A DNA sequence that is grammatically correct in 1D (`delta_score`) can still translate into a physically useless protein that fails to fold correctly in 3D (`pLDDT` score).

**Why This Matters**:
- Evo2 scores sequences in 1D (sequence likelihood)
- AlphaFold 3 scores sequences in 3D (structural viability)
- **These are NOT equivalent** - high 1D score ≠ high 3D viability

**Solution**: Multi-dimensional validation is mandatory:
1. **Phase I: The Forge** - Generate candidates
2. **Phase II: The Sieve** - Use sequence-level likelihood scores (`delta_score`) as fast filter
3. **Phase III: The Gauntlet** - Use 3D structural prediction (`pLDDT` score) as final arbiter

**Rule**: Any candidate with `pLDDT < 70` is a "wet noodle" → **DISCARD**

**Code Example**:
```python
# ❌ WRONG: Only check delta_score
if delta_score > threshold:
    return candidate  # May be "wet noodle"

# ✅ RIGHT: Multi-dimensional validation
if delta_score > threshold:
    plddt = await alphafold3.predict_structure(candidate)
    if plddt >= 70:  # Structurally viable
        return candidate
    else:
        return None  # "Wet noodle" - discard
```

---

### **9.1.2 The Triumvirate Protocol**

**Problem**: Evo2's `delta_likelihood_score` is fundamentally insensitive to frameshift and nonsense mutations. Relying on it as a single-factor assessment for all variant types is a critical vulnerability.

**Why This Matters**:
- Frameshift mutations cause catastrophic truncation
- Evo2 may score them as "low impact" (blind spot)
- This leads to false negatives (missed pathogenic variants)

**Solution**: Multi-layered approach (Triumvirate Protocol):
1. **Truncation Check** (deterministic) - Bioinformatic translation of CDS to screen for catastrophic mutations
2. **Evo2 Deep Learning** - Only non-truncating variants pass to oracle
3. **Pathway Analysis** - Additional context for complex cases

**Code Example**:
```python
# ❌ WRONG: Direct Evo2 scoring
delta_score = await evo2.score_variant(mutation)
if delta_score > threshold:
    return "pathogenic"

# ✅ RIGHT: Triumvirate Protocol
# Step 1: Truncation check
if is_truncating(mutation):  # Frameshift/nonsense
    return "pathogenic"  # Deterministic - no Evo2 needed

# Step 2: Evo2 for non-truncating variants
if not is_truncating(mutation):
    delta_score = await evo2.score_variant(mutation)
    if delta_score > threshold:
        return "pathogenic"

# Step 3: Pathway analysis (if needed)
pathway_burden = await pathway_service.analyze(mutation)
```

---

### **9.1.3 Backend Orchestrator Pattern**

**Problem**: Legacy applications may attempt to make direct, fine-grained calls to specialized backend services. This is brittle and creates tight coupling.

**Why This Matters**:
- Frontend needs to know about multiple services
- Changes to service APIs break frontend
- Difficult to mock services for testing
- Complex error handling across multiple services

**Solution**: Single, powerful orchestrator endpoint in intermediary backend:
- Frontend calls ONE endpoint
- Orchestrator manages entire multi-stage workflow
- Orchestrator calls Modal services, external APIs, etc.
- Frontend becomes "dumb" client

**Code Example**:
```python
# ❌ WRONG: Frontend calls multiple services
# Frontend
await fetch('/api/evo/score_variant', ...)
await fetch('/api/pathway/aggregate', ...)
await fetch('/api/evidence/literature', ...)

# ✅ RIGHT: Single orchestrator endpoint
# Frontend
await fetch('/api/clinical_genomics/analyze_variant', {
    mutations: [...],
    disease: "ovarian_carcinoma"
})

# Backend orchestrator
@router.post("/analyze_variant")
async def analyze_variant(request):
    # Orchestrator manages all services
    seq_scores = await evo2_service.score(...)
    pathway_scores = await pathway_service.aggregate(...)
    evidence = await evidence_service.gather(...)
    return unified_response
```

**Benefits**:
- Frontend simplicity (one endpoint)
- Service changes don't break frontend
- Easy to mock for testing
- Centralized error handling
- Progressive enhancement (mock unimplemented services)

---

### **9.1.4 Generative vs Inference Paradigm**

**Problem**: The dominant academic paradigm for heterogeneity analysis relies on **inference-based modeling** (e.g., "PhenoPop"). These methods attempt to deconstruct noisy, low-resolution `in vitro` data to *guess* at biological reality.

**Why This Matters**:
- Inference-based methods have "ground truth problem" - cannot validate against patient reality
- Our **Digital Twin** built from high-fidelity sequencing *is* the ground truth
- Generative paradigm (`Digital Twin -> Predict -> Generate`) is revolutionary replacement, not incremental improvement

**Our Approach**:
- **Digital Twin**: High-fidelity sequencing data (ground truth)
- **Predict**: Evo2-powered predictions
- **Generate**: Zeta Forge generates novel therapeutics

**Code Example**:
```python
# ❌ WRONG: Inference-based (guessing from noisy data)
def predict_from_in_vitro_data(noisy_data):
    # Deconstruct noisy data to guess reality
    return guessed_reality

# ✅ RIGHT: Generative paradigm (ground truth → predict → generate)
def predict_from_digital_twin(high_fidelity_sequencing):
    # Ground truth from sequencing
    digital_twin = build_digital_twin(high_fidelity_sequencing)
    
    # Predict from ground truth
    predictions = evo2.predict(digital_twin)
    
    # Generate novel therapeutics
    therapeutics = zeta_forge.generate(predictions)
    return therapeutics
```

---

## **9.2 BEST PRACTICES**

### **9.2.1 Graceful Degradation**

**Pattern**: System remains operational even if external services fail.

**Implementation**:
```python
# ✅ RIGHT: Graceful degradation
try:
    evidence = await literature_service.search(...)
except Exception as e:
    logger.warning(f"Evidence gathering failed: {e}")
    # Continue with S/P only (no E)
    evidence = None
    provenance["evidence_fallback"] = True

# Continue processing with partial results
drugs = score_drugs(seq_scores, pathway_scores, evidence)
```

**Benefits**:
- System remains operational
- Better user experience (partial results vs complete failure)
- Resilience to network issues

---

### **9.2.2 Provenance Tracking**

**Pattern**: Complete audit trails for all operations.

**Implementation**:
```python
# ✅ RIGHT: Complete provenance
provenance = {
    "run_id": str(uuid.uuid4()),
    "profile": "baseline",
    "cache": "miss",  # or "hit"
    "flags": {
        "fusion_active": bool(os.getenv("FUSION_AM_URL")),
        "evo_use_delta_only": bool(os.getenv("EVO_USE_DELTA_ONLY")),
        "evidence_enabled": True
    },
    "sequence_scoring": {
        "method": "evo2_1b",
        "model": "evo2_1b",
        "window_used": 8192
    },
    "fallback": {
        "evidence_timeout": True,
        "clinvar_unavailable": False
    }
}
```

**Benefits**:
- Reproducibility (can rerun exact same analysis)
- Transparency (users see how results were computed)
- Compliance (audit trail for clinical use)
- Debugging (understand why results differ)

---

### **9.2.3 Feature Flags**

**Pattern**: Environment-based toggles for different operational profiles.

**Implementation**:
```python
# ✅ RIGHT: Feature flags
def get_feature_flags():
    return {
        "enable_massive_modes": os.getenv("ENABLE_MASSIVE_MODES", "false") == "true",
        "disable_evo2": os.getenv("DISABLE_EVO2", "false") == "true",
        "disable_fusion": os.getenv("DISABLE_FUSION", "false") == "true",
        "evo_use_delta_only": os.getenv("EVO_USE_DELTA_ONLY", "true") == "true",
        "evo_spam_safe": os.getenv("EVO_SPAM_SAFE", "true") == "true",
        "evo_max_models": int(os.getenv("EVO_MAX_MODELS", "1")),
        "evo_max_flanks": int(os.getenv("EVO_MAX_FLANKS", "1"))
    }

# Usage
flags = get_feature_flags()
if not flags["disable_evo2"]:
    scores = await evo2_service.score(...)
```

**Benefits**:
- Graceful degradation (disable features if services unavailable)
- A/B testing capabilities
- Demo vs production modes
- Spam-safety controls

---

### **9.2.4 Modular Architecture**

**Pattern**: Services are small (~100-150 lines), single responsibility.

**Implementation**:
```python
# ✅ RIGHT: Modular services
# api/services/sequence_scorers/evo2_scorer.py (~150 lines)
class Evo2Scorer:
    async def score(self, mutation, model_id, options):
        # Single responsibility: Evo2 scoring
        pass

# api/services/pathway/aggregation.py (~100 lines)
def aggregate_pathways(mutations, panel):
    # Single responsibility: Pathway aggregation
    pass

# api/services/evidence/literature_client.py (~120 lines)
async def search_literature(gene, hgvs_p, disease):
    # Single responsibility: Literature search
    pass
```

**Benefits**:
- Easy to test (test services independently)
- Clear ownership
- Reusable business logic
- Clear API contracts

---

### **9.2.5 Spam-Safety Controls**

**Pattern**: Prevent excessive and costly API calls.

**Implementation**:
```python
# ✅ RIGHT: Spam-safety controls
EVO_SPAM_SAFE = os.getenv("EVO_SPAM_SAFE", "true") == "true"
EVO_MAX_MODELS = int(os.getenv("EVO_MAX_MODELS", "1" if EVO_SPAM_SAFE else "3"))
EVO_MAX_FLANKS = int(os.getenv("EVO_MAX_FLANKS", "1" if EVO_SPAM_SAFE else "5"))
EVO_DISABLE_SYMMETRY = os.getenv("EVO_DISABLE_SYMMETRY", "true" if EVO_SPAM_SAFE else "false") == "true"

# Usage
if EVO_SPAM_SAFE:
    # Single model, single flank, no symmetry
    scores = await evo2.score(mutation, model_id="evo2_1b", flanks=[4096])
else:
    # Multiple models, multiple flanks, symmetry
    scores = await evo2.score_ensemble(mutations, models=["evo2_1b", "evo2_7b"], flanks=[4096, 8192, 16384])
```

**Benefits**:
- Cost control (prevent excessive API calls)
- Performance (faster responses)
- Resource management (prevent service overload)

---

### **9.2.6 Caching Strategy**

**Pattern**: Redis → Memory fallback with single-flight protection.

**Implementation**:
```python
# ✅ RIGHT: Caching with fallback
async def get_cache(key: str) -> Optional[Any]:
    # Try Redis first
    client = get_redis_client()
    if client:
        cached = await client.get(f"cache:{key}")
        if cached:
            return json.loads(cached)
    
    # Fallback to memory
    if key in _memory_cache:
        timestamp = _cache_timestamps.get(key, 0)
        if time.time() - timestamp < 3600:  # 1 hour TTL
            return _memory_cache[key]
    
    return None

# Single-flight protection
async def with_singleflight(key: str, fn: Callable) -> Any:
    # Prevent duplicate concurrent requests
    lock_acquired = await client.set(f"lock:{key}", "1", nx=True, ex=60)
    if lock_acquired:
        result = await fn()
        await set_cache(f"cache:{key}", result)
        return result
    else:
        # Wait for other request to complete
        await asyncio.sleep(0.1)
        return await get_cache(f"cache:{key}")
```

**Benefits**:
- Performance (faster responses for cached data)
- Cost reduction (fewer API calls)
- Resilience (memory fallback if Redis unavailable)
- Deduplication (single-flight prevents duplicate requests)

---

## **9.3 ANTI-PATTERNS (What to Avoid)**

### **9.3.1 Single-Metric Myopia**

**Problem**: Optimizing for a single metric (e.g., `delta_score`) is insufficient.

**Example**:
```python
# ❌ WRONG: Only optimize delta_score
if delta_score > threshold:
    return candidate  # May be "wet noodle" (structurally non-viable)

# ✅ RIGHT: Multi-dimensional optimization
if delta_score > threshold:
    plddt = await alphafold3.predict_structure(candidate)
    if plddt >= 70:  # Structurally viable
        return candidate
```

**Lesson**: True fitness is multi-dimensional (sequence + structure + function).

---

### **9.3.2 Direct Service Calls from Frontend**

**Problem**: Frontend making direct calls to multiple backend services.

**Example**:
```python
# ❌ WRONG: Frontend calls multiple services
# Frontend
await fetch('/api/evo/score_variant', ...)
await fetch('/api/pathway/aggregate', ...)
await fetch('/api/evidence/literature', ...)

# ✅ RIGHT: Single orchestrator endpoint
# Frontend
await fetch('/api/clinical_genomics/analyze_variant', {...})
```

**Lesson**: Use backend orchestrator pattern for complex workflows.

---

### **9.3.3 Assuming Data Schema**

**Problem**: Assuming structure or content of data source without verification.

**Example**:
```python
# ❌ WRONG: Assume schema
data = json.loads(response.text)
variant = data["variant"]  # May not exist

# ✅ RIGHT: Verify schema first
data = json.loads(response.text)
variant = data.get("variant")  # Safe access
if not variant:
    raise ValueError("Missing variant field")
```

**Lesson**: Always verify schema before coding. Use `head`, `grep` to inspect data first.

---

### **9.3.4 Overly Strict Filters**

**Problem**: Starting with overly strict filters that erroneously discard critical data.

**Example**:
```python
# ❌ WRONG: Overly strict filter
if review_status == "reviewed by expert panel":
    include_variant()  # Misses variants with "criteria provided"

# ✅ RIGHT: Filter loosely, then refine
if "criteria provided" in review_status:
    include_variant()  # Includes both "reviewed by expert panel" and "criteria provided"
```

**Lesson**: Filter loosely, then refine. Start with slightly looser criteria.

---

### **9.3.5 Ignoring Resource Requirements**

**Problem**: Not scaling resources when upgrading data sources.

**Example**:
```python
# ❌ WRONG: Same resources for small and large datasets
@app.cls(cpu=2, memory="4Gi")  # Insufficient for GRCh38
class EvoService:
    def load_model(self):
        # Fails with out-of-memory for large datasets
        pass

# ✅ RIGHT: Scale resources proportionally
@app.cls(cpu=8, memory="64Gi", timeout=1800)  # Scaled for GRCh38
class EvoService:
    def load_model(self):
        # Sufficient resources for large datasets
        pass
```

**Lesson**: When upgrading from small to large datasets, scale resources proportionally.

---

### **9.3.6 Silent Failures**

**Problem**: Services failing silently without logging or fallback.

**Example**:
```python
# ❌ WRONG: Silent failure
try:
    result = await external_api.call()
except Exception:
    pass  # Silent failure - no logging, no fallback

# ✅ RIGHT: Log and fallback
try:
    result = await external_api.call()
except Exception as e:
    logger.warning(f"External API failed: {e}")
    result = default_result  # Fallback
    provenance["fallback"] = True
```

**Lesson**: Always log failures and provide fallbacks. Silent failures are hard to debug.

---

## **9.4 LESSONS LEARNED FROM .CURSORRULES**

### **9.4.1 Modal Deployment Lessons**

**Lesson 1: Correct Deployment Syntax**
- ❌ **Wrong**: `app.cls()(MyClass)` - Extra call after decorator
- ✅ **Right**: `@app.cls(...)` decorator is sufficient

**Lesson 2: Service-to-Service Calls**
- ❌ **Wrong**: HTTP calls for internal Modal services
- ✅ **Right**: Use `modal.Cls.lookup("app-name")` for low-latency communication

**Lesson 3: Resource Scaling**
- ❌ **Wrong**: Same resources for small and large datasets
- ✅ **Right**: Scale resources proportionally (CPU, memory, timeout)

**Lesson 4: Build-Essential Is Critical**
- ✅ **Right**: Always include `.apt_install("build-essential")` for ML libraries

**Lesson 5: Memory Requirements Are Massive**
- ✅ **Right**: Loading 25GB AlphaMissense + 650M ESM requires 64GB+ RAM

---

### **9.4.2 Data Source Lessons**

**Lesson 1: ClinVar Data Sources**
- ❌ **Wrong**: E-utils API for bulk data (rate limiting, unreliable)
- ❌ **Wrong**: `clinvar.vcf.gz` (not comprehensive, missing canonical variants)
- ✅ **Right**: `variant_summary.txt.gz` (most comprehensive and reliable)

**Lesson 2: GTEx Portal API**
- ⚠️ **Warning**: GTEx Portal API is unreliable and documentation is inconsistent
- ✅ **Right**: Consider external dependency as critical risk

**Lesson 3: Verify Schema Before Coding**
- ✅ **Right**: Always perform preliminary inspection (`head`, `grep`) before parsing

**Lesson 4: PyVCF3 Is Unstable**
- ❌ **Wrong**: Relying on PyVCF3 for VCF parsing
- ✅ **Right**: Write custom parser with `pandas` or line-by-line processing

---

### **9.4.3 Evo2 Oracle Lessons**

**Lesson 1: Dual Nature (Sniper vs. Poet)**
- **As Scorer**: Precision instrument for zero-shot variant effect prediction
- **As Generator**: Long-form sequence author (requires >=1000bp context, generates >=500bp)

**Lesson 2: Pathological Prompt Doctrine**
- **High-Quality Context** → High-Quality Output (e.g., BRAF gene)
- **Low-Quality Context** → Pathological Output (e.g., poly-A tracts)

**Lesson 3: The "Wet Noodle" Genesis**
- Previous failures were from violating prompt quality doctrine
- Short-form generative tasks with low-quality bait sequences → junk DNA

---

### **9.4.4 Development Workflow Lessons**

**Lesson 1: Venv Purge & Re-Forge Protocol**
When facing persistent `ModuleNotFoundError`s:
1. **Purge**: `rm -rf venv`
2. **Re-Forge**: `python -m venv venv`
3. **Re-Arm**: `venv/bin/pip install -r requirements.txt`
4. **Install Dev Tools**: Manually install `pytest`, `streamlit`, etc.

**Lesson 2: Service Version Supremacy**
- ✅ **Right**: Always verify latest deployment script before writing client code
- ❌ **Wrong**: Relying on outdated URLs or schemas

**Lesson 3: Client-Side Resilience**
- ✅ **Right**: Test clients must handle self-signed SSL certificates (`verify=False`)

**Lesson 4: Defensive Downstream Communication**
- ✅ **Right**: Never implicitly trust downstream service schema
- ✅ **Right**: Use safe defaults (`response.get("key", default_value)`)

---

## **9.5 CODE ORGANIZATION PATTERNS**

### **9.5.1 Router Pattern**

**Pattern**: Thin endpoints delegating to services.

**Example**:
```python
# ✅ RIGHT: Thin router
@router.post("/predict")
async def predict_efficacy(request: Dict[str, Any]):
    # Delegate to orchestrator
    response = await orchestrator.predict(efficacy_request)
    return response

# ❌ WRONG: Business logic in router
@router.post("/predict")
async def predict_efficacy(request: Dict[str, Any]):
    # Business logic in router (hard to test, not reusable)
    scores = await evo2.score(...)
    pathways = aggregate_pathways(...)
    drugs = score_drugs(...)
    return drugs
```

---

### **9.5.2 Service Layer Pattern**

**Pattern**: Business logic in services, routers are thin.

**Example**:
```python
# ✅ RIGHT: Service layer
# api/services/efficacy_orchestrator/orchestrator.py
class EfficacyOrchestrator:
    async def predict(self, request: EfficacyRequest):
        # Business logic here
        seq_scores = await sequence_processor.process(...)
        pathway_scores = await pathway_service.aggregate(...)
        drugs = await drug_scorer.score(...)
        return EfficacyResponse(drugs=drugs, ...)

# Router delegates
@router.post("/predict")
async def predict_efficacy(request: Dict[str, Any]):
    orchestrator = create_efficacy_orchestrator()
    return await orchestrator.predict(efficacy_request)
```

---

### **9.5.3 Schema Validation Pattern**

**Pattern**: Pydantic models for request/response validation.

**Example**:
```python
# ✅ RIGHT: Pydantic models
class AnalyzeVariantRequest(BaseModel):
    mutations: List[Dict[str, Any]]
    disease: Optional[str] = None
    profile: str = "baseline"
    include: List[str] = []

@router.post("/analyze_variant")
async def analyze_variant(request: AnalyzeVariantRequest):
    # Automatic validation
    if not request.mutations:
        raise HTTPException(status_code=400, detail="mutations required")
    # ...
```

---

## **9.6 ERROR HANDLING PATTERNS**

### **9.6.1 Timeout Handling**

**Pattern**: Use `asyncio.wait_for()` with service-specific timeouts.

**Example**:
```python
# ✅ RIGHT: Timeout handling
try:
    result = await asyncio.wait_for(
        call_external_api(),
        timeout=60.0
    )
except asyncio.TimeoutError:
    logger.warning("External API timeout")
    return default_result
```

---

### **9.6.2 HTTP Error Handling**

**Pattern**: Catch `httpx.HTTPStatusError` and provide fallback.

**Example**:
```python
# ✅ RIGHT: HTTP error handling
try:
    response = await client.post(url, json=payload)
    response.raise_for_status()
except httpx.HTTPStatusError as e:
    logger.error(f"HTTP error: {e}")
    return default_result
```

---

### **9.6.3 Graceful Degradation**

**Pattern**: Continue processing with partial results.

**Example**:
```python
# ✅ RIGHT: Graceful degradation
try:
    evidence_scores = await evidence_service.gather(...)
except Exception as e:
    logger.warning(f"Evidence gathering failed: {e}")
    # Continue with S/P only (no E)
    evidence_scores = None
    provenance["evidence_fallback"] = True

# Continue processing
drugs = score_drugs(seq_scores, pathway_scores, evidence_scores)
```

---

## **9.7 SUMMARY**

### **Key Doctrines**:
1. **"Wet Noodle" Doctrine**: Multi-dimensional validation is mandatory
2. **Triumvirate Protocol**: Truncation check before Evo2
3. **Backend Orchestrator Pattern**: Single orchestrator endpoint
4. **Generative vs Inference Paradigm**: Digital Twin → Predict → Generate

### **Best Practices**:
1. **Graceful Degradation**: System remains operational with partial results
2. **Provenance Tracking**: Complete audit trails for all operations
3. **Feature Flags**: Environment-based toggles for different modes
4. **Modular Architecture**: Small services with single responsibility
5. **Spam-Safety Controls**: Prevent excessive API calls
6. **Caching Strategy**: Redis → Memory fallback with single-flight

### **Anti-Patterns to Avoid**:
1. **Single-Metric Myopia**: Optimizing for single metric only
2. **Direct Service Calls**: Frontend calling multiple services
3. **Assuming Data Schema**: Not verifying data structure
4. **Overly Strict Filters**: Discarding critical data
5. **Ignoring Resource Requirements**: Not scaling resources
6. **Silent Failures**: No logging or fallback

### **Lessons Learned**:
1. **Modal Deployment**: Correct syntax, service-to-service calls, resource scaling
2. **Data Sources**: ClinVar sources, GTEx API, schema verification
3. **Evo2 Oracle**: Dual nature, prompt quality, "wet noodle" genesis
4. **Development Workflow**: Venv purge, service version, client resilience

---

**Status**: ✅ **CYCLE 5 COMPLETE** - Development Patterns & Lessons Learned  
**Next**: Cycle 6 (I10) - Product Capabilities & Positioning

---



