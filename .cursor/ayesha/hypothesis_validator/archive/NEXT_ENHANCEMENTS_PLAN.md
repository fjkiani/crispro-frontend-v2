# ⚔️ NEXT ENHANCEMENTS PLAN - FOOD VALIDATOR

**Status:** 📋 **IMPLEMENTATION GUIDE**  
**Priority:** P1 (Post-MVP Enhancements)  
**Estimated Time:** 8-12 hours total

---

## **🎯 OVERVIEW**

This document provides concrete implementation plans for 4 enhancement areas identified in Priority Fixes completion:

1. **Async LLM Wrapper** - Full LLM-powered timing recommendations
2. **Enhanced Mechanism Extraction** - Sophisticated parsing from abstracts
3. **SAE Rules Expansion** - Systematic compound addition
4. **Performance Optimization** - Caching, timeouts, scaling

---

## **1. ASYNC LLM WRAPPER FOR TIMING RECOMMENDATIONS**

### **Current State:**
- Pattern matching from abstracts (sync, limited)
- No real LLM synthesis for timing

### **Target State:**
- Full LLM-powered extraction with structured output
- Async to avoid blocking
- Fallback to pattern matching if LLM unavailable

### **Implementation:**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/dietician_recommendations.py`

**Step 1: Add Async Method**

```python
async def generate_timing_recommendations_async(
    self, 
    compound: str, 
    evidence: Dict[str, Any],
    use_llm: bool = True
) -> Dict[str, Any]:
    """
    Generate timing recommendations with optional LLM synthesis.
    
    Returns:
        {
            "best_time": "morning with breakfast",
            "with_food": True,
            "timing_rationale": "Based on LLM analysis: ...",
            "meal_suggestions": [...],
            "method": "llm" | "pattern_match" | "hardcoded"
        }
    """
    # Try LLM first if available
    if use_llm and LLM_AVAILABLE and evidence.get("papers"):
        try:
            llm_service = get_llm_service()
            
            # Build context from papers
            papers = evidence.get("papers", [])[:5]
            context = self._build_timing_context(compound, papers)
            
            # Call LLM for structured extraction
            prompt = f"""Extract optimal timing recommendations for {compound} from the following research papers.

Context:
{context}

Provide a JSON response with:
- best_time: Specific time recommendation (e.g., "morning with breakfast", "evening", "with meals")
- with_food: true/false
- timing_rationale: Brief explanation based on the papers
- meal_suggestions: List of meal types or specific foods if mentioned

Papers:
{self._format_papers_for_prompt(papers)}

Respond in JSON format only."""

            # Use LLM service (adjust based on your LLM client)
            llm_response = await llm_service.query(prompt)
            
            # Parse JSON response
            try:
                import json
                result = json.loads(llm_response)
                
                return {
                    "best_time": result.get("best_time", "As directed"),
                    "with_food": result.get("with_food", True),
                    "timing_rationale": result.get("timing_rationale", "Based on literature review"),
                    "meal_suggestions": result.get("meal_suggestions", []),
                    "method": "llm"
                }
            except json.JSONDecodeError:
                # LLM didn't return JSON, fall through to pattern matching
                pass
        except Exception as e:
            print(f"⚠️ LLM timing extraction failed: {e}")
    
    # Fallback to sync pattern matching (existing method)
    return self.generate_timing_recommendations(compound, evidence)
```

**Step 2: Update Router to Use Async**

**File:** `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

```python
@router.post("/api/hypothesis/validate_food_dynamic")
async def validate_food_dynamic(request: Dict[str, Any]):
    # ... existing code ...
    
    # Use async timing if available
    if hasattr(dietician_service, 'generate_timing_recommendations_async'):
        timing = await dietician_service.generate_timing_recommendations_async(
            compound=compound,
            evidence=evidence_result,
            use_llm=True
        )
    else:
        # Fallback to sync
        timing = dietician_service.generate_timing_recommendations(compound, evidence_result)
    
    # ... rest of code ...
```

**Step 3: Add Helper Methods**

```python
def _build_timing_context(self, compound: str, papers: List[Dict]) -> str:
    """Build context string from papers for LLM prompt."""
    contexts = []
    for paper in papers[:3]:
        abstract = paper.get("abstract", "")[:300]  # Limit length
        contexts.append(f"Title: {paper.get('title', '')}\nAbstract: {abstract}")
    return "\n\n".join(contexts)

def _format_papers_for_prompt(self, papers: List[Dict]) -> str:
    """Format papers for LLM prompt."""
    formatted = []
    for i, paper in enumerate(papers[:5], 1):
        formatted.append(
            f"Paper {i}:\n"
            f"  Title: {paper.get('title', '')}\n"
            f"  Abstract: {paper.get('abstract', '')[:400]}"
        )
    return "\n\n".join(formatted)
```

**Estimated Time:** 2-3 hours  
**Testing:** Test with known compounds, verify LLM fallback, test async behavior

---

## **2. ENHANCED MECHANISM EXTRACTION**

### **Current State:**
- Simple keyword matching (`_extract_mechanisms_from_text`)
- Limited to 6 mechanism categories

### **Target State:**
- Named Entity Recognition (NER) for targets
- Pathway relationship extraction
- Confidence scoring per mechanism
- Integration with pathway databases

### **Implementation:**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/enhanced_evidence_service.py`

**Step 1: Enhanced Mechanism Extraction**

```python
def _extract_mechanisms_enhanced(self, text: str, pathways: List[str] = None) -> List[Dict[str, Any]]:
    """
    Enhanced mechanism extraction with confidence scoring.
    
    Returns:
        [
            {
                "mechanism": "anti-inflammatory",
                "confidence": 0.85,
                "targets": ["NF-κB", "COX-2"],
                "evidence_snippet": "...",
                "pathway_alignment": ["inflammation"]
            }
        ]
    """
    mechanisms = []
    text_lower = text.lower()
    
    # Expanded mechanism database with targets and keywords
    mechanism_db = {
        "anti-inflammatory": {
            "keywords": ["inflammation", "nf-kb", "nf-κb", "nf-kappa-b", "cox-2", "prostaglandin", "il-6", "tnf-α"],
            "targets": ["NF-κB", "COX-2", "IL-6", "TNF-α"],
            "pathways": ["inflammation", "nf-kb signaling"]
        },
        "antioxidant": {
            "keywords": ["antioxidant", "oxidative stress", "ros", "glutathione", "superoxide", "free radical"],
            "targets": ["GSH", "SOD", "CAT"],
            "pathways": ["oxidative stress", "redox"]
        },
        "angiogenesis_inhibition": {
            "keywords": ["angiogenesis", "vegf", "vascular", "endothelial", "neovascularization"],
            "targets": ["VEGF", "VEGFR", "FGF"],
            "pathways": ["angiogenesis", "vascular development"]
        },
        "dna_repair": {
            "keywords": ["dna repair", "brca", "parp", "homologous recombination", "hrd", "double-strand break"],
            "targets": ["BRCA1", "BRCA2", "PARP1", "RAD51"],
            "pathways": ["dna repair", "homologous recombination"]
        },
        "apoptosis_induction": {
            "keywords": ["apoptosis", "cell death", "caspase", "bcl-2", "p53", "programmed cell death"],
            "targets": ["Bcl-2", "Bax", "Caspase-3", "p53"],
            "pathways": ["apoptosis", "cell death"]
        },
        "cell_cycle": {
            "keywords": ["cell cycle", "cdk", "cyclin", "checkpoint", "g1/s", "g2/m"],
            "targets": ["CDK1", "CDK2", "Cyclin B", "Cyclin D"],
            "pathways": ["cell cycle"]
        },
        "metabolism": {
            "keywords": ["metabolism", "mtor", "glycolysis", "warburg", "atp", "mitochondria"],
            "targets": ["mTOR", "AMPK", "GLUT1"],
            "pathways": ["metabolism", "warburg effect"]
        },
        "epigenetic": {
            "keywords": ["histone", "dna methylation", "epigenetic", "hdac", "dnmt"],
            "targets": ["HDAC", "DNMT", "HAT"],
            "pathways": ["epigenetic regulation"]
        },
        "immunomodulation": {
            "keywords": ["immune", "t-cell", "nk cell", "cytokine", "interferon", "pd-1"],
            "targets": ["T-cells", "NK cells", "PD-1", "CTLA-4"],
            "pathways": ["immune system", "t-cell activation"]
        },
        "autophagy": {
            "keywords": ["autophagy", "autophagosome", "lc3", "beclin", "mtor"],
            "targets": ["LC3", "Beclin-1", "mTOR"],
            "pathways": ["autophagy"]
        }
    }
    
    # Extract mechanisms with confidence
    for mechanism, data in mechanism_db.items():
        keywords = data["keywords"]
        matches = sum(1 for kw in keywords if kw in text_lower)
        
        if matches > 0:
            # Confidence based on keyword match count and specificity
            confidence = min(0.5 + (matches * 0.1), 0.95)
            
            # Find evidence snippet (sentence containing keywords)
            sentences = text.split('.')
            evidence_snippet = ""
            for sent in sentences:
                if any(kw in sent.lower() for kw in keywords[:3]):  # Check top 3 keywords
                    evidence_snippet = sent[:150].strip()
                    break
            
            # Check pathway alignment
            pathway_alignment = []
            if pathways:
                for pathway in pathways:
                    pathway_lower = pathway.lower()
                    if any(pw in pathway_lower for pw in data.get("pathways", [])):
                        pathway_alignment.append(pathway)
            
            mechanisms.append({
                "mechanism": mechanism,
                "confidence": round(confidence, 2),
                "targets": data.get("targets", []),
                "evidence_snippet": evidence_snippet,
                "pathway_alignment": pathway_alignment,
                "keyword_matches": matches
            })
    
    # Sort by confidence
    mechanisms.sort(key=lambda x: x["confidence"], reverse=True)
    return mechanisms[:10]  # Top 10
```

**Step 2: Integrate with LLM for Advanced Extraction**

```python
async def _extract_mechanisms_llm(self, compound: str, papers: List[Dict]) -> List[Dict[str, Any]]:
    """
    Use LLM to extract mechanisms with structured output.
    
    Returns structured mechanism data with targets and pathways.
    """
    if not LLM_AVAILABLE:
        return []
    
    try:
        llm_service = get_llm_service()
        
        # Build prompt
        papers_text = "\n\n".join([
            f"Title: {p['title']}\nAbstract: {p['abstract'][:500]}"
            for p in papers[:5]
        ])
        
        prompt = f"""Extract mechanisms of action for {compound} from these research papers.

Papers:
{papers_text}

For each mechanism, provide:
- mechanism_name: Brief name (e.g., "anti-inflammatory", "dna_repair")
- targets: List of molecular targets (e.g., ["NF-κB", "COX-2"])
- pathways: List of affected pathways (e.g., ["inflammation", "nf-kb signaling"])
- confidence: 0-1 score based on evidence strength
- evidence: Brief quote from paper supporting this mechanism

Respond as JSON array:
[
  {{
    "mechanism_name": "...",
    "targets": [...],
    "pathways": [...],
    "confidence": 0.85,
    "evidence": "..."
  }}
]"""

        response = await llm_service.query(prompt)
        
        # Parse JSON
        import json
        mechanisms = json.loads(response)
        
        return mechanisms
        
    except Exception as e:
        print(f"⚠️ LLM mechanism extraction failed: {e}")
        return []
```

**Step 3: Update synthesize_evidence_llm to Use Enhanced Extraction**

```python
# In synthesize_evidence_llm method:
# Replace:
mechanisms = self._extract_mechanisms_from_text(summary)

# With:
if LLM_AVAILABLE:
    mechanisms_llm = await self._extract_mechanisms_llm(compound, papers)
    if mechanisms_llm:
        mechanisms = [m.get("mechanism_name", "") for m in mechanisms_llm]
    else:
        mechanisms = self._extract_mechanisms_enhanced(summary, pathways).mechanisms
else:
    mechanisms = self._extract_mechanisms_enhanced(summary, pathways)
```

**Estimated Time:** 3-4 hours  
**Testing:** Test with various compounds, verify target extraction, test confidence scores

---

## **3. SAE RULES EXPANSION (SYSTEMATIC APPROACH)**

### **Current State:**
- 22 compounds in `supplement_treatment_rules.json`
- Manual addition process

### **Target State:**
- 50+ compounds
- Systematic categorization
- Template-based addition
- Validation script

### **Implementation:**

**Step 1: Create Compound Template**

**File:** `.cursor/ayesha/hypothesis_validator/scripts/add_compound_template.json`

```json
{
  "compound_name": "NewCompound",
  "mechanism": "primary_mechanism",
  "high_appropriateness_contexts": [
    "context1",
    "context2"
  ],
  "default_scores": {
    "line_appropriateness": 0.7,
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.75
  },
  "biomarker_gates": {
    "biomarker_name": "expected_value"
  }
}
```

**Step 2: Create Addition Script**

**File:** `.cursor/ayesha/hypothesis_validator/scripts/add_compound_to_sae.py`

```python
#!/usr/bin/env python3
"""
Script to add new compounds to supplement_treatment_rules.json

Usage:
    python add_compound_to_sae.py --name "Coenzyme Q10" --mechanism "mitochondrial_support" --line_app 0.75
"""

import json
import argparse
from pathlib import Path

def add_compound(
    compound_name: str,
    mechanism: str,
    line_appropriateness: float,
    cross_resistance: float = 0.0,
    sequencing_fitness: float = None,
    contexts: list = None,
    biomarker_gates: dict = None
):
    """Add compound to supplement_treatment_rules.json"""
    
    rules_file = Path(__file__).parent.parent / "data/supplement_treatment_rules.json"
    
    # Load existing rules
    with open(rules_file) as f:
        rules_data = json.load(f)
    
    # Default sequencing_fitness
    if sequencing_fitness is None:
        sequencing_fitness = line_appropriateness - 0.05
    
    # Build compound rule
    compound_rule = {
        "mechanism": mechanism,
        "high_appropriateness_contexts": contexts or [],
        "default_scores": {
            "line_appropriateness": line_appropriateness,
            "cross_resistance": cross_resistance,
            "sequencing_fitness": sequencing_fitness
        }
    }
    
    if biomarker_gates:
        compound_rule["biomarker_gates"] = biomarker_gates
    
    # Add to rules
    rules_data["supplement_rules"][compound_name] = compound_rule
    
    # Save
    with open(rules_file, 'w') as f:
        json.dump(rules_data, f, indent=2)
    
    print(f"✅ Added {compound_name} to supplement_treatment_rules.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--mechanism", required=True)
    parser.add_argument("--line_app", type=float, required=True)
    parser.add_argument("--cross_res", type=float, default=0.0)
    parser.add_argument("--seq_fit", type=float)
    parser.add_argument("--contexts", nargs="+", default=[])
    args = parser.parse_args()
    
    add_compound(
        compound_name=args.name,
        mechanism=args.mechanism,
        line_appropriateness=args.line_app,
        cross_resistance=args.cross_res,
        sequencing_fitness=args.seq_fit,
        contexts=args.contexts
    )
```

**Step 3: Create Batch Addition from food_targets.json**

**File:** `.cursor/ayesha/hypothesis_validator/scripts/sync_sae_from_food_targets.py`

```python
#!/usr/bin/env python3
"""
Sync SAE rules from food_targets.json

Automatically adds compounds from food_targets.json to supplement_treatment_rules.json
with reasonable defaults based on mechanisms.
"""

import json
from pathlib import Path

# Mechanism to default scores mapping
MECHANISM_SCORES = {
    "dna_repair_support": {"line_app": 0.9, "seq_fit": 0.85},
    "oxidative_stress_recovery": {"line_app": 1.0, "seq_fit": 0.95},
    "anti_inflammatory": {"line_app": 0.85, "seq_fit": 0.80},
    "nfkb_inhibition": {"line_app": 0.7, "seq_fit": 0.75},
    "antioxidant": {"line_app": 0.7, "seq_fit": 0.70},
    "immune_modulation": {"line_app": 0.8, "seq_fit": 0.85},
    "mitochondrial_support": {"line_app": 0.75, "seq_fit": 0.80},
}

def infer_mechanism_from_targets(targets: list) -> str:
    """Infer mechanism from target list"""
    targets_str = " ".join(targets).lower()
    
    if any(x in targets_str for x in ["brca", "parp", "dna repair"]):
        return "dna_repair_support"
    elif any(x in targets_str for x in ["glutathione", "oxidative", "ros"]):
        return "oxidative_stress_recovery"
    elif any(x in targets_str for x in ["nf-kb", "co and ", "inflammation"]):
        return "anti_inflammatory"
    elif any(x in targets_str for x in ["vegf", "angiogenesis"]):
        return "angiogenesis_inhibition"
    else:
        return "general_support"

def sync_compounds():
    """Sync compounds from food_targets.json"""
    base_path = Path(__file__).parent.parent
    food_targets_file = base_path / "data/food_targets.json"
    rules_file = base_path / "data/supplement_treatment_rules.json"
    
    # Load files
    with open(food_targets_file) as f:
        food_targets = json.load(f)
    
    with open(rules_file) as f:
        rules_data = json.load(f)
    
    added_count = 0
    
    # Process each compound
    for food in food_targets.get("foods", []):
        compound_name = food.get("compound", "")
        if not compound_name:
            continue
        
        # Skip if already exists
        if compound_name in rules_data["supplement_rules"]:
            continue
        
        # Infer mechanism
        targets = food.get("B_targets", [])
        mechanism = infer_mechanism_from_targets(targets)
        
        # Get default scores
        defaults = MECHANISM_SCORES.get(mechanism, {"line_app": 0.7, "seq_fit": 0.70})
        
        # Infer contexts
        contexts = []
        if "HRD" in str(targets) or "BRCA" in str(targets):
            contexts.append("hrd_positive")
        if "oxidative" in str(targets).lower():
            contexts.append("oxidative_stress_high")
        if "inflammation" in str(food.get("mechanisms", [])).lower():
            contexts.append("chronic_inflammation")
        
        # Add compound
        rules_data["supplement_rules"][compound_name] = {
            "mechanism": mechanism,
            "high_appropriateness_contexts": contexts,
            "default_scores": {
                "line_appropriateness": defaults["line_app"],
                "cross_resistance": 0.0,
                "sequencing_fitness": defaults["seq_fit"]
            }
        }
        
        added_count += 1
        print(f"✅ Added {compound_name} (mechanism: {mechanism})")
    
    # Save
    with open(rules_file, 'w') as f:
        json.dump(rules_data, f, indent=2)
    
    print(f"\n🎉 Added {added_count} compounds to SAE rules")

if __name__ == "__main__":
    sync_compounds()
```

**Step 4: Priority Compounds to Add**

```
High Priority (Common supplements):
- Berberine (metabolic support)
- Alpha-lipoic acid (already added)
- Probiotics (already added)
- Magnesium (already added)
- Vitamin B complex (dna synthesis)
- Iron (if anemia present)
- Calcium (bone health post-chemo)
- Chromium (metabolism)
- Manganese (antioxidant cofactor)
- Copper (immune support)

Medium Priority (Evidence-based):
- Astragalus (immune modulation)
- Milk thistle (liver support)
- Milk thistle (liver support)
- Ginseng (energy, immune)
- Echinacea (immune - short term)
- Saw palmetto (hormonal)
- St. John's Wort (depression - interactions!)
- Ginkgo biloba (circulation)
- Valerian (sleep)
```

**Estimated Time:** 2-3 hours  
**Testing:** Run sync script, validate JSON, test new compounds

---

## **4. PERFORMANCE OPTIMIZATION (CACHING & TIMEOUTS)**

### **Current State:**
- Simple in-memory cache in `EnhancedEvidenceService`
- No TTL management
- No timeout configuration
- No request deduplication

### **Target State:**
- Redis-based caching with TTL
- Configurable timeouts
- Request deduplication (single-flight pattern)
- Performance monitoring

### **Implementation:**

**Step 1: Create Cache Service**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/cache_service.py`

```python
"""
Cache Service for Food Validator

Uses Redis if available, falls back to in-memory LRU cache.
"""

import json
import hashlib
from typing import Optional, Any
from functools import lru_cache
import time

# Try Redis
try:
    import redis
    REDIS_AVAILABLE = True
except ImportError:
    REDIS_AVAILABLE = False

# In-memory cache fallback
_memory_cache = {}
_cache_timestamps = {}

class CacheService:
    """Unified caching service with Redis + memory fallback"""
    
    def __init__(self, redis_url: Optional[str] = None, default_ttl: int = 600):
        self.default_ttl = default_ttl  # 10 minutes
        self.redis_client = None
        
        if REDIS_AVAILABLE and redis_url:
            try:
                self.redis_client = redis.from_url(redis_url, decode_responses=True)
                # Test connection
                self.redis_client.ping()
                print("✅ Redis cache initialized")
            except Exception as e:
                print(f"⚠️ Redis unavailable: {e}, using memory cache")
                self.redis_client = None
        else:
            print("⚠️ Redis not configured, using memory cache")
    
    def _make_key(self, prefix: str, *args) -> str:
        """Create cache key from prefix and args"""
        key_str = f"{prefix}:{json.dumps(args, sort_keys=True)}"
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def get(self, prefix: str, *args) -> Optional[Any]:
        """Get from cache"""
        key = self._make_key(prefix, *args)
        
        if self.redis_client:
            try:
                value = self.redis_client.get(key)
                if value:
                    return json.loads(value)
            except Exception as e:
                print(f"⚠️ Redis get error: {e}")
        
        # Fallback to memory
        if key in _memory_cache:
            timestamp = _cache_timestamps.get(key, 0)
            if time.time() - timestamp < self.default_ttl:
                return _memory_cache[key]
            else:
                # Expired
                del _memory_cache[key]
                del _cache_timestamps[key]
        
        return None
    
    def set(self, prefix: str, value: Any, ttl: Optional[int] = None, *args):
        """Set in cache"""
        key = self._make_key(prefix, *args)
        ttl = ttl or self.default_ttl
        
        serialized = json.dumps(value)
        
        if self.redis_client:
            try:
                self.redis_client.setex(key, ttl, serialized)
                return
            except Exception as e:
                print(f"⚠️ Redis set error: {e}")
        
        # Fallback to memory
        _memory_cache[key] = value
        _cache_timestamps[key] = time.time()
    
    def invalidate(self, prefix: str, *args):
        """Invalidate cache entry"""
        key = self._make_key(prefix, *args)
        
        if self.redis_client:
            try:
                self.redis_client.delete(key)
            except:
                pass
        
        if key in _memory_cache:
            del _memory_cache[key]
            del _cache_timestamps[key]
```

**Step 2: Add Single-Flight Pattern**

```python
# In enhanced_evidence_service.py

from functools import lru_cache
import asyncio

# Request deduplication
_active_requests = {}

async def get_complete_evidence_with_cache(
    self,
    compound: str,
    disease: str,
    pathways: List[str] = None,
    cache_service: Optional[CacheService] = None
) -> Dict[str, Any]:
    """
    Get evidence with caching and single-flight deduplication.
    """
    # Check cache
    if cache_service:
        cached = cache_service.get("evidence", compound, disease, pathways)
        if cached:
            return cached
    
    # Single-flight: if request already in progress, wait for it
    request_key = (compound, disease, tuple(pathways or []))
    if request_key in _active_requests:
        return await _active_requests[request_key]
    
    # Start new request
    future = asyncio.create_task(
        self._get_evidence_uncached(compound, disease, pathways)
    )
    _active_requests[request_key] = future
    
    try:
        result = await future
        
        # Cache result
        if cache_service:
            cache_service.set("evidence", result, ttl=600, compound, disease, pathways)
        
        return result
    finally:
        # Clean up
        _active_requests.pop(request_key, None)

async def _get_evidence_uncached(
    self,
    compound: str,
    disease: str,
    pathways: List[str] = None
) -> Dict[str, Any]:
    """Internal method without caching"""
    # Existing get_complete_evidence logic here
    pass
```

**Step 3: Add Timeouts**

```python
# In enhanced_evidence_service.py __init__

import httpx
from httpx import Timeout

def __init__(self):
    # ... existing code ...
    
    # Configurable timeouts
    self.timeout = Timeout(
        connect=5.0,  # Connection timeout
        read=15.0,    # Read timeout
        write=10.0,   # Write timeout
        pool=5.0      # Pool timeout
    )
    
    # Use httpx with timeout
    self.http_client = httpx.AsyncClient(timeout=self.timeout)
```

**Step 4: Performance Monitoring**

```python
import time
from functools import wraps

def monitor_performance(func):
    """Decorator to monitor function performance"""
    @wraps(func)
    async def wrapper(*args, **kwargs):
        start = time.time()
        try:
            result = await func(*args, **kwargs)
            duration = time.time() - start
            
            # Log slow operations (>1s)
            if duration > 1.0:
                print(f"⚠️ Slow operation: {func.__name__} took {duration:.2f}s")
            
            return result
        except Exception as e:
            duration = time.time() - start
            print(f"❌ Error in {func.__name__} after {duration:.2f}s: {e}")
            raise
    return wrapper

# Usage:
@monitor_performance
async def get_complete_evidence(...):
    ...
```

**Estimated Time:** 3-4 hours  
**Testing:** Test caching behavior, verify TTL, test single-flight, measure performance improvements

---

## **📊 IMPLEMENTATION PRIORITY**

### **High Priority (Do First):**
1. **Performance Optimization** - Critical for production readiness
   - **Time:** 3-4 hours
   - **Impact:** User experience, API responsiveness

2. **SAE Rules Expansion** - Quick wins, more compounds
   - **Time:** 2-3 hours
   - **Impact:** Coverage increase

### **Medium Priority:**
3. **Enhanced Mechanism Extraction** - Better insights
   - **Time:** 3-4 hours
   - **Impact:** More accurate mechanism identification

4. **Async LLM Wrapper** - Better timing recommendations
   - **Time:** 2-3 hours
   - **Impact:** More intelligent recommendations

---

## **🧪 TESTING STRATEGY**

### **For Each Enhancement:**

1. **Unit Tests:**
   - Test individual functions
   - Mock external dependencies
   - Verify edge cases

2. **Integration Tests:**
   - Test full workflows
   - Verify cache behavior
   - Test error handling

3. **Performance Tests:**
   - Measure latency improvements
   - Test cache hit rates
   - Load testing

4. **Regression Tests:**
   - Ensure existing functionality works
   - No breaking changes

---

## **🚀 QUICK START**

### **To Implement All 4:**

```bash
# 1. Performance (Critical)
cd oncology-coPilot/oncology-backend-minimal
# Create cache_service.py (3-4 hours)
# Integrate into services (1 hour)

# 2. SAE Expansion (Quick Win)
cd .cursor/ayesha/hypothesis_validator
# Run sync_sae_from_food_targets.py (30 min)
# Add priority compounds manually (1-2 hours)

# 3. Mechanism Extraction (Enhancement)
# Update enhanced_evidence_service.py (3-4 hours)
# Test with various compounds (1 hour)

# 4. Async LLM Wrapper (Polish)
# Update dietician_recommendations.py (2-3 hours)
# Update router (30 min)
```

**Total Estimated Time:** 12-18 hours

---

**⚔️ READY TO IMPLEMENT - ALL ENHANCEMENTS PLANNED**

