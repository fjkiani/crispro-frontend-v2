# 🚀 CHROMADB → ASTRADB MIGRATION PLAN

## **⚔️ MISSION: Eliminate Local File Dependencies for Production Deployment**

**Goal:** Migrate from ChromaDB PersistentClient (local files) to AstraDB (cloud) for serverless compatibility and scalability.

---

## **📊 CURRENT STATE AUDIT**

### **✅ ACTIVE & PRODUCTION-READY:**
1. **Agent 1 Seeding** (`scripts/agent_1_seeding/`)
   - **Vector DB**: ChromaDB PersistentClient (LOCAL FILES) ⚠️
   - **Embeddings**: Google API (`models/embedding-001`) - 768 dimensions ✅
   - **Status**: ACTIVE, but deployment-blocking
   - **Location**: `scripts/agent_1_seeding/database/chromadb_client.py`

2. **ClinicalTrialAgent** (`backend/agents/clinical_trial_agent.py`)
   - **Vector DB**: AstraDB (CLOUD) ✅
   - **Embeddings**: SentenceTransformer (`all-MiniLM-L6-v2`) - 384 dimensions ⚠️
   - **Status**: ACTIVE, but dimension mismatch
   - **Collection**: `clinical_trials` (384-dim)

### **⚠️ LEGACY/MONOLITH (Review & Consolidate):**
3. **load_trials_from_api.py** (`backend/scripts/`)
   - **Vector DB**: AstraDB ✅
   - **Embeddings**: SentenceTransformer (384-dim) ⚠️
   - **Status**: OLD pipeline, uses different schema
   - **Collection**: `clinical_trials` (384-dim)

4. **load_trials_to_astra.py** (`scripts/`)
   - **Vector DB**: AstraDB ✅
   - **Embeddings**: SentenceTransformer (384-dim) ⚠️
   - **Status**: LEGACY, uses Firecrawl (not CT.gov API v2)
   - **Collection**: `clinical_trials` (384-dim)

### **🔧 INFRASTRUCTURE (Ready to Use):**
5. **DatabaseConnections** (`backend/utils/database_connections.py`)
   - **Status**: ✅ Already supports AstraDB!
   - **Method**: `get_vector_db_collection(collection_name)`
   - **Dependencies**: `astrapy` (already in requirements)

---

## **🚨 CRITICAL ISSUES IDENTIFIED**

### **Issue 1: Dimension Mismatch**
- **Agent 1 (Seeding)**: 768-dim (Google Embeddings)
- **ClinicalTrialAgent (Query)**: 384-dim (SentenceTransformer)
- **Problem**: Cannot query 768-dim embeddings with 384-dim queries
- **Impact**: Search won't work after migration

### **Issue 2: Deployment Incompatibility**
- **ChromaDB PersistentClient**: Requires persistent volumes (won't work in serverless)
- **Impact**: Agent 1 seeding fails in Vercel/Modal/Lambda deployments

### **Issue 3: Inconsistent Implementations**
- 3 different AstraDB scripts with different:
  - Schemas (NCI IDs vs NCT IDs)
  - Embedding methods (Google vs SentenceTransformer)
  - Collection names (`clinical_trials` vs `clinical_trials_eligibility`)

---

## **🎯 MIGRATION STRATEGY**

### **Phase 1: Standardize on 768-Dim Google Embeddings** ✅ (Recommended)

**Rationale:**
- Already used by Agent 1 (no seeding change needed)
- Zero cold starts (API-based)
- Consistent with production deployment analysis recommendations
- Cost: ~$0.01 per 1000 embeddings (very affordable)

**Changes:**
1. ✅ Keep Google Embeddings in Agent 1 (no change)
2. ⚠️ Update ClinicalTrialAgent to use Google Embeddings (instead of SentenceTransformer)
3. ⚠️ Re-create AstraDB collection with 768 dimensions
4. ⚠️ Re-seed trials with 768-dim embeddings

### **Phase 2: Migrate Agent 1 to AstraDB** ✅ (Required for Production)

**Replace:**
- `chromadb.PersistentClient` → `DatabaseConnections.get_vector_db_collection()`
- Local file storage → Cloud AstraDB collection

**Benefits:**
- ✅ Serverless compatible (no local files)
- ✅ Auto-scaling (cloud infrastructure)
- ✅ Multi-instance support (shared data)
- ✅ No deployment complexity

### **Phase 3: Consolidate Legacy Scripts** ✅ (Cleanup)

**Actions:**
- Mark `load_trials_to_astra.py` as DEPRECATED (uses Firecrawl, not CT.gov API v2)
- Document `load_trials_from_api.py` as legacy (different schema)
- Consolidate to single seeding pipeline: Agent 1

---

## **📋 DETAILED MIGRATION PLAN**

### **STEP 1: Update Agent 1 - ChromaDB → AstraDB** (2-3 hours)

**File:** `scripts/agent_1_seeding/database/chromadb_client.py`

**Current Implementation:**
```python
chroma_client = chromadb.PersistentClient(path=str(abs_chroma_path))
collection = chroma_client.get_or_create_collection(...)
```

**New Implementation:**
```python
from backend.utils.database_connections import DatabaseConnections

db_manager = DatabaseConnections()
astra_db = db_manager.get_vector_db_connection()
collection = astra_db.get_collection("clinical_trials_eligibility")
```

**Key Changes:**
1. Remove `chromadb` dependency
2. Use `DatabaseConnections` helper (already exists!)
3. Keep Google Embeddings (already correct - 768-dim)
4. Update collection name to `clinical_trials_eligibility` (matches Agent 1 config)

**New File Structure:**
```
scripts/agent_1_seeding/database/
├── chromadb_client.py  → RENAME to astradb_client.py
└── __init__.py         → Update import
```

**Upsert API (Similar Interface):**
```python
# ChromaDB API:
collection.upsert(ids=[...], documents=[...], metadatas=[...])

# AstraDB API (astrapy v2):
collection.insert_many([{
    "_id": source_url,
    "$vector": embedding,  # Generated by Google API
    "eligibility_text": eligibility_text,
    "nct_id": nct_id,
    "title": title,
    "status": status,
    "disease_category": disease_category
}])
```

**Embedding Generation (Keep Google API):**
```python
from chromadb.utils import embedding_functions

google_ef = embedding_functions.GoogleGenerativeAiEmbeddingFunction(
    api_key=os.getenv("GOOGLE_API_KEY"),
    model_name="models/embedding-001"  # 768 dimensions
)

# Generate embeddings for batch
eligibility_texts = [t["eligibility_text"] for t in batch]
embeddings = google_ef(eligibility_texts)  # Returns list of 768-dim vectors
```

---

### **STEP 2: Verify AstraDB Collection Setup** (30 min)

**Check Collection Dimension:**
```python
from backend.utils.database_connections import DatabaseConnections

db_manager = DatabaseConnections()
astra_db = db_manager.get_vector_db_connection()
collection = astra_db.get_collection("clinical_trials_eligibility")

# Check if collection exists and dimension
# If 384-dim, need to re-create with 768-dim
```

**Re-create Collection (if needed):**
```python
from astrapy.constants import VectorMetric

# Drop old collection
astra_db.drop_collection("clinical_trials_eligibility")

# Create new with 768 dimensions
astra_db.create_collection(
    "clinical_trials_eligibility",
    definition={
        "vector": {
            "dimension": 768,  # Google embedding dimension
            "metric": VectorMetric.COSINE
        }
    }
)
```

---

### **STEP 3: Update ClinicalTrialAgent - SentenceTransformer → Google Embeddings** (1 hour)

**File:** `backend/agents/clinical_trial_agent.py`

**Current (384-dim SentenceTransformer):**
```python
from sentence_transformers import SentenceTransformer

self.embedding_model = SentenceTransformer('all-MiniLM-L6-v2')
query_embedding = self.embedding_model.encode(rich_query_text).tolist()
```

**New (768-dim Google Embeddings):**
```python
from chromadb.utils import embedding_functions

google_api_key = os.getenv("GOOGLE_API_KEY")
if not google_api_key:
    raise ValueError("GOOGLE_API_KEY required for query embeddings")

self.google_ef = embedding_functions.GoogleGenerativeAiEmbeddingFunction(
    api_key=google_api_key,
    model_name="models/embedding-001"  # 768 dimensions
)

# In run() method:
query_embedding = self.google_ef([rich_query_text])[0]  # Returns 768-dim vector
```

**Remove SentenceTransformer Dependency:**
```python
# Remove from __init__:
# self.embedding_model = SentenceTransformer(...)

# Update requirements.txt:
# Remove: sentence-transformers
```

---

### **STEP 4: Update Collection Name Consistency** (30 min)

**Standardize Collection Names:**
- **Seeding (Agent 1)**: `clinical_trials_eligibility` (new)
- **Query (ClinicalTrialAgent)**: `clinical_trials_eligibility` (update from `clinical_trials`)

**Update ClinicalTrialAgent:**
```python
# Current:
self.trials_collection = self.db_connections.get_vector_db_collection("clinical_trials")

# New:
self.trials_collection = self.db_connections.get_vector_db_collection("clinical_trials_eligibility")
```

---

### **STEP 5: Re-seed Trials with 768-Dim Embeddings** (1-2 hours for 1000 trials)

**Run Agent 1 with New AstraDB Client:**
```bash
cd oncology-coPilot/oncology-backend
python scripts/agent_1_seeding/main.py --limit 1000
```

**Verify Embeddings:**
```python
# Check collection dimension and count
collection = astra_db.get_collection("clinical_trials_eligibility")
# Should show 768-dim vectors and 1000 trials
```

---

### **STEP 6: Update Tests** (1 hour)

**Files to Update:**
- `tests/agent_1_seeding/test_database.py` - Update ChromaDB mocks → AstraDB mocks
- `tests/agent_1_seeding/test_integration.py` - Update end-to-end tests

**Mock Strategy:**
```python
# Use unittest.mock to mock DatabaseConnections
from unittest.mock import Mock, patch
from backend.utils.database_connections import DatabaseConnections

@patch('scripts.agent_1_seeding.database.astradb_client.DatabaseConnections')
def test_embed_trials_batched(mock_db_class):
    mock_collection = Mock()
    mock_db = Mock()
    mock_db.get_vector_db_connection.return_value = mock_db
    mock_db.get_collection.return_value = mock_collection
    mock_db_class.return_value = mock_db
    
    # Test embedding logic...
```

---

### **STEP 7: Cleanup Legacy Scripts** (30 min)

**Mark as Deprecated:**
- `scripts/load_trials_to_astra.py` - Add deprecation notice
- `scripts/load_trials_from_json.py` - Add deprecation notice

**Add Comments:**
```python
# DEPRECATED: This script uses Firecrawl and SentenceTransformer (384-dim)
# Use scripts/agent_1_seeding/main.py instead (CT.gov API v2 + Google Embeddings 768-dim)
# Migration date: [DATE]
```

---

## **📊 MIGRATION CHECKLIST**

### **Pre-Migration:**
- [ ] Audit existing AstraDB collection dimension (384 vs 768)
- [ ] Verify `GOOGLE_API_KEY` is set in environment
- [ ] Verify `ASTRA_DB_API_ENDPOINT` and `ASTRA_DB_APPLICATION_TOKEN` are set
- [ ] Backup existing ChromaDB data (if needed for rollback)

### **Phase 1: Update Agent 1 (Seeding)**
- [ ] Rename `chromadb_client.py` → `astradb_client.py`
- [ ] Replace `chromadb.PersistentClient` with `DatabaseConnections`
- [ ] Update upsert logic to use `astrapy` API
- [ ] Keep Google Embeddings (768-dim) - no change
- [ ] Update imports in `main.py` and `__init__.py`
- [ ] Remove `chromadb` from requirements.txt (or keep if other scripts use it)

### **Phase 2: Update ClinicalTrialAgent (Query)**
- [ ] Replace SentenceTransformer with Google Embedding API
- [ ] Update query embedding generation (768-dim)
- [ ] Update collection name to `clinical_trials_eligibility`
- [ ] Remove `sentence-transformers` from requirements.txt
- [ ] Update tests

### **Phase 3: Collection Setup**
- [ ] Drop old `clinical_trials` collection (384-dim) OR keep for compatibility
- [ ] Create `clinical_trials_eligibility` collection (768-dim)
- [ ] Verify collection dimension matches embeddings

### **Phase 4: Re-seeding**
- [ ] Run Agent 1 seeding script with new AstraDB client
- [ ] Verify 1000 trials embedded successfully
- [ ] Check embedding dimension (should be 768)
- [ ] Verify metadata fields match expected schema

### **Phase 5: Testing**
- [ ] Unit tests pass (AstraDB mocks)
- [ ] Integration test: Seed → Query → Verify results
- [ ] End-to-end: Search query returns relevant trials
- [ ] Performance test: Query latency <500ms

### **Phase 6: Cleanup**
- [ ] Mark legacy scripts as deprecated
- [ ] Update documentation
- [ ] Remove ChromaDB path from config (optional)
- [ ] Update `.cursorrules` scratchpad

---

## **🔧 IMPLEMENTATION DETAILS**

### **New File: `scripts/agent_1_seeding/database/astradb_client.py`**

```python
"""
AstraDB embedding client for Agent 1 - replaces ChromaDB for production deployment.
Uses Google Embeddings (768-dim) via DatabaseConnections helper.
"""
import os
import asyncio
import logging
from typing import List, Dict, Any, Optional
from pathlib import Path

from chromadb.utils import embedding_functions

# Use centralized database connections
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))
from backend.utils.database_connections import DatabaseConnections

from .. import config
from ..utils.logger import setup_logger
from ..utils.error_handler import should_fail_on_chromadb_errors

logger = setup_logger(__name__)


async def embed_trials_batched(trials: List[Dict[str, Any]], batch_size: Optional[int] = None) -> None:
    """
    Embed trial eligibility text into AstraDB with rate limiting and retries.
    
    Uses Google Embeddings (768-dim) for consistency with Agent 1 seeding.
    Replaces ChromaDB PersistentClient for serverless deployment compatibility.
    
    Args:
        trials: List of parsed trial dictionaries
        batch_size: Batch size for embeddings (defaults to config.CHROMA_BATCH_SIZE)
        
    Raises:
        RuntimeError: If critical failures occur (>50% batches fail)
        ValueError: If GOOGLE_API_KEY or AstraDB credentials are missing
    """
    if not trials:
        logger.warning("No trials to embed")
        return
    
    if batch_size is None:
        batch_size = config.CHROMA_BATCH_SIZE
    
    # Check for Google API key
    google_api_key = os.getenv("GOOGLE_API_KEY")
    if not google_api_key:
        raise ValueError("GOOGLE_API_KEY not found in environment variables")
    
    # Initialize Google embedding function (768-dim)
    google_ef = embedding_functions.GoogleGenerativeAiEmbeddingFunction(
        api_key=google_api_key,
        model_name="models/embedding-001"  # 768 dimensions
    )
    
    # Initialize AstraDB connection
    logger.info("Initializing AstraDB connection...")
    db_manager = DatabaseConnections()
    astra_db = db_manager.get_vector_db_connection()
    
    if not astra_db:
        raise ValueError("Failed to get AstraDB connection. Check ASTRA_DB_API_ENDPOINT and ASTRA_DB_APPLICATION_TOKEN")
    
    collection = astra_db.get_collection("clinical_trials_eligibility")
    logger.info(f"Connected to AstraDB collection: clinical_trials_eligibility")
    
    logger.info(f"Embedding {len(trials)} trials into AstraDB (batch size: {batch_size})...")
    
    total_batches = (len(trials) + batch_size - 1) // batch_size
    failed_batches = 0
    
    for i in range(0, len(trials), batch_size):
        batch = trials[i:i + batch_size]
        batch_num = (i // batch_size) + 1
        
        # Prepare batch data
        eligibility_texts = [t["eligibility_text"] for t in batch if t.get("eligibility_text")]
        
        if not eligibility_texts:
            logger.warning(f"Skipping batch {batch_num} - no eligibility text")
            continue
        
        # Generate embeddings using Google API (768-dim)
        try:
            embeddings = google_ef(eligibility_texts)  # Returns list of 768-dim vectors
        except Exception as e:
            logger.error(f"Failed to generate embeddings for batch {batch_num}: {e}")
            failed_batches += 1
            continue
        
        # Prepare documents for AstraDB
        documents_to_insert = []
        for j, trial in enumerate(batch):
            if not trial.get("eligibility_text"):
                continue
            
            doc = {
                "_id": trial["source_url"],  # Use source_url as document ID
                "$vector": embeddings[j],    # 768-dim vector
                "eligibility_text": trial["eligibility_text"],
                "nct_id": trial.get("nct_id", ""),
                "title": trial.get("title", ""),
                "status": trial.get("status", ""),
                "disease_category": trial.get("disease_category", "")
            }
            documents_to_insert.append(doc)
        
        # Retry logic with exponential backoff
        success = False
        for retry in range(config.CHROMA_MAX_RETRIES):
            try:
                collection.insert_many(documents_to_insert)
                success = True
                embedded_count = min(i + batch_size, len(trials))
                logger.info(
                    f"✅ Embedded batch {batch_num}/{total_batches} "
                    f"({embedded_count}/{len(trials)} trials)"
                )
                break  # Success - exit retry loop
                
            except Exception as e:
                if retry == config.CHROMA_MAX_RETRIES - 1:
                    # Final retry failed
                    logger.error(
                        f"❌ AstraDB batch {batch_num}/{total_batches} failed after "
                        f"{config.CHROMA_MAX_RETRIES} retries: {e}"
                    )
                    failed_batches += 1
                else:
                    # Wait before retry (exponential backoff)
                    wait_time = 2 ** retry
                    logger.warning(
                        f"AstraDB batch {batch_num} failed, retrying in {wait_time}s "
                        f"(attempt {retry + 1}/{config.CHROMA_MAX_RETRIES}): {e}"
                    )
                    await asyncio.sleep(wait_time)
        
        # Rate limiting: sleep between batches (except for last batch)
        if i + batch_size < len(trials) and success:
            await asyncio.sleep(config.CHROMA_RATE_LIMIT)
        
        # Check if should fail (critical failure threshold)
        if should_fail_on_chromadb_errors(
            failed_batches=failed_batches,
            total_batches=batch_num,
            max_failure_rate=config.CHROMA_MAX_FAILURE_RATE
        ):
            raise RuntimeError(
                f"Critical AstraDB failure: {failed_batches}/{batch_num} batches failed "
                f"({failed_batches/batch_num:.1%} > {config.CHROMA_MAX_FAILURE_RATE:.1%} threshold)"
            )
    
    logger.info(f"✅ AstraDB embedding complete: {len(trials)} trials embedded")
```

---

## **🧪 TESTING STRATEGY**

### **Unit Tests:**
```python
# tests/agent_1_seeding/test_astradb_client.py
from unittest.mock import Mock, patch
from scripts.agent_1_seeding.database.astradb_client import embed_trials_batched

@patch('scripts.agent_1_seeding.database.astradb_client.DatabaseConnections')
@patch('scripts.agent_1_seeding.database.astradb_client.embedding_functions.GoogleGenerativeAiEmbeddingFunction')
async def test_embed_trials_batched(mock_google_ef, mock_db_class):
    # Mock AstraDB collection
    mock_collection = Mock()
    mock_db = Mock()
    mock_db.get_collection.return_value = mock_collection
    mock_db_class.return_value.get_vector_db_connection.return_value = mock_db
    
    # Mock Google embeddings
    mock_google_ef.return_value.return_value = [[0.1] * 768 for _ in range(5)]
    
    # Test embedding
    trials = [{"source_url": f"http://test/{i}", "eligibility_text": f"text {i}"} for i in range(5)]
    await embed_trials_batched(trials, batch_size=5)
    
    # Verify upsert called
    assert mock_collection.insert_many.called
```

### **Integration Test:**
```python
# Test: Seed → Query → Verify
# 1. Run Agent 1 seeding (with AstraDB)
# 2. Query ClinicalTrialAgent with test query
# 3. Verify results match expected trials
```

---

## **💰 COST ANALYSIS**

### **AstraDB:**
- **Free Tier**: 20M reads/month, 1M writes/month
- **Paid Tier**: ~$25-50/month for production scale (100K+ trials)
- **Verdict**: Very affordable, fully managed

### **Google Embeddings:**
- **Pricing**: $0.0001 per 1K characters
- **1000 trials**: ~200K characters = $0.02 (one-time)
- **Verdict**: Negligible cost (~$0.01-0.10 per 1000 embeddings)

---

## **🎯 ROLLBACK PLAN**

### **If Migration Fails:**
1. **Keep ChromaDB client** as backup (`chromadb_client.py.backup`)
2. **Keep old collection** (`clinical_trials` 384-dim) for query fallback
3. **Feature flag** to toggle between ChromaDB/AstraDB (if needed)

---

## **✅ SUCCESS CRITERIA**

### **Phase 1 Complete:**
- [ ] Agent 1 seeds 1000 trials to AstraDB successfully
- [ ] Embeddings are 768-dim (verified in collection)
- [ ] No ChromaDB dependency in Agent 1

### **Phase 2 Complete:**
- [ ] ClinicalTrialAgent queries AstraDB successfully
- [ ] Query embeddings are 768-dim (matches seeded data)
- [ ] Search returns relevant results (<500ms latency)
- [ ] No SentenceTransformer dependency

### **Production Ready:**
- [ ] All tests pass
- [ ] Deployment works in serverless (Vercel/Modal)
- [ ] No persistent storage required
- [ ] Auto-scaling verified

---

## **📝 NOTES**

### **Why Keep Google Embeddings:**
- ✅ Already used by Agent 1 (no seeding change)
- ✅ Zero cold starts (API-based)
- ✅ Consistent dimensions (768 across system)
- ✅ Production deployment analysis recommends this approach

### **Why AstraDB over ChromaDB:**
- ✅ Serverless compatible (no local files)
- ✅ Auto-scaling (cloud infrastructure)
- ✅ Multi-instance support (shared data)
- ✅ Fully managed (no maintenance)
- ✅ Already integrated (`DatabaseConnections` helper exists)

---

**STATUS**: 📋 **PLAN COMPLETE** - Ready for implementation  
**ESTIMATED TIME**: 6-8 hours total  
**PRIORITY**: P0 (Required for production deployment)

