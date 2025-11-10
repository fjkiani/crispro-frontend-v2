# 🔄 AGENT 1 RELOCATION PLAN - Main Backend → Minimal Backend

## **MISSION: Migrate Agent 1 Seeding Scripts to Minimal Backend**

### **Current State**
- **Location:** `oncology-coPilot/oncology-backend/scripts/agent_1_seeding/`
- **Database Path:** Points to main backend: `backend/data/clinical_trials.db`
- **ChromaDB Path:** Points to main backend: `backend/data/chroma_data`
- **Status:** ⚠️ Tightly coupled to main backend structure

### **Target State**
- **Location:** `oncology-coPilot/oncology-backend-minimal/scripts/agent_1_seeding/`
- **Database Path:** Point to minimal backend: `data/clinical_trials.db`
- **ChromaDB:** ❌ **DEPRECATED** - Migrate to AstraDB
- **Status:** ✅ Self-contained in minimal backend

---

## **MIGRATION STRATEGY - SIMPLIFIED**

### **Option A: Leave Agent 1 in Main Backend (RECOMMENDED) ⭐**

**Rationale:**
- Agent 1 is a **seeding pipeline**, not a runtime service
- Runs as a **one-time or periodic batch job**, not API endpoint
- Main backend already has Agent 1 complete and functional
- No need to duplicate seeding scripts

**Changes Required:**
1. Update Agent 1 `config.py` database paths to point to **shared minimal backend DB**:
   ```python
   SQLITE_DB_PATH = "../oncology-backend-minimal/data/clinical_trials.db"
   ```
2. Remove ChromaDB seeding - Agent 1 seeds SQLite only
3. Minimal backend `ClinicalTrialSearchService` reads from AstraDB (Agent 1 doesn't touch AstraDB)
4. Keep Agent 1 scripts in main backend as "data pipeline"

**Benefits:**
- ✅ No code duplication
- ✅ Agent 1 remains in main backend (where other agents live)
- ✅ Minimal backend stays focused on production API
- ✅ Clear separation: main backend = data pipelines, minimal backend = API services

---

### **Option B: Move Agent 1 to Minimal Backend (NOT RECOMMENDED)**

**Changes Required:**
1. Copy entire `agent_1_seeding/` directory to minimal backend
2. Update all import paths
3. Update database paths in config
4. Duplicate testing infrastructure

**Why Not:**
- ❌ Code duplication
- ❌ Mixes data pipeline with API code
- ❌ More complex maintenance

---

## **RECOMMENDED IMPLEMENTATION: Option A**

### **Phase 1: Update Agent 1 Config (5 min)**

Edit `oncology-coPilot/oncology-backend/scripts/agent_1_seeding/config.py`:

```python
# Database Configuration
SQLITE_BATCH_SIZE = 100
# Point to minimal backend's shared database
SQLITE_DB_PATH = "../oncology-backend-minimal/data/clinical_trials.db"

# ChromaDB Configuration - DEPRECATED (migrating to AstraDB)
# CHROMA_BATCH_SIZE = 50
# CHROMA_RATE_LIMIT = 1.0
# CHROMA_MAX_RETRIES = 3
# CHROMA_PATH = "backend/data/chroma_data"  # DEPRECATED
# CHROMA_COLLECTION = "clinical_trials_eligibility"  # Now in AstraDB
```

### **Phase 2: Remove ChromaDB Seeding from Agent 1 (10 min)**

Edit `oncology-coPilot/oncology-backend/scripts/agent_1_seeding/main.py`:

- Remove `embed_trials_batched` call (lines ~80-90)
- Remove `--skip-embeddings` flag (no longer needed)
- Update summary report to remove ChromaDB stats

**Rationale:** Embeddings now handled by AstraDB + Google Embedding API in minimal backend's search service.

### **Phase 3: Add AstraDB Seeding Script (Separate, Optional) (30 min)**

Create `oncology-coPilot/oncology-backend-minimal/scripts/seed_astradb_from_sqlite.py`:

```python
"""
Seed AstraDB from SQLite clinical trials database.
Run this AFTER Agent 1 has populated SQLite.
"""
import asyncio
from api.services.database_connections import get_db_connections
from api.services.clinical_trial_search_service import ClinicalTrialSearchService

async def seed_astradb():
    """Read trials from SQLite, embed, and upsert to AstraDB."""
    db = get_db_connections()
    service = ClinicalTrialSearchService()
    
    # Read all trials from SQLite
    conn = db.get_sqlite_connection()
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM clinical_trials")
    trials = cursor.fetchall()
    
    # Upsert to AstraDB
    collection = db.get_vector_db_collection(service.collection_name)
    
    for trial in trials:
        trial_dict = dict(trial)
        embedding = service._generate_embedding(trial_dict['eligibility_text'])
        
        collection.upsert({
            "_id": trial_dict['source_url'],
            "vector": embedding,
            "nct_id": trial_dict['nct_id'],
            "title": trial_dict['title'],
            "status": trial_dict['status'],
            "disease_category": trial_dict['disease_category'],
            "biomarker_requirements": trial_dict['biomarker_requirements'],
            "locations_data": trial_dict['locations_data'],
            "eligibility_text": trial_dict['eligibility_text'],
            "description_text": trial_dict['description_text']
        })
    
    print(f"✅ Seeded {len(trials)} trials to AstraDB")

if __name__ == "__main__":
    asyncio.run(seed_astradb())
```

**Run:**
```bash
cd oncology-coPilot/oncology-backend-minimal
venv/bin/python scripts/seed_astradb_from_sqlite.py
```

---

## **EXECUTION PLAN**

### **Step 1:** Update Agent 1 Config
```bash
# Edit config.py
SQLITE_DB_PATH = "../oncology-backend-minimal/data/clinical_trials.db"
```

### **Step 2:** Run Agent 1 (Seeds SQLite Only)
```bash
cd oncology-coPilot/oncology-backend
venv/bin/python -m scripts.agent_1_seeding.main --limit 1000
```

### **Step 3:** Seed AstraDB (Separate Step)
```bash
cd oncology-coPilot/oncology-backend-minimal
venv/bin/python scripts/seed_astradb_from_sqlite.py
```

### **Step 4:** Verify Search Service
```bash
curl -X POST http://localhost:8000/api/search-trials \
  -H 'Content-Type: application/json' \
  -d '{"query": "ovarian cancer BRCA1 mutation trials"}'
```

---

## **ACCEPTANCE CRITERIA**

✅ Agent 1 seeds SQLite in minimal backend's `data/` directory
✅ Agent 1 no longer touches ChromaDB
✅ AstraDB seeding is a separate, optional script
✅ Search service (`/api/search-trials`) reads from AstraDB
✅ Main backend keeps Agent 1 as "data pipeline"
✅ Minimal backend focuses on API services

---

## **BENEFITS OF THIS APPROACH**

1. **Clear Separation of Concerns:**
   - Main backend = Data pipelines + Batch jobs + Agents
   - Minimal backend = Production API + Services

2. **No Code Duplication:**
   - Agent 1 stays in one place
   - No duplicate testing/maintenance

3. **Flexible Data Flow:**
   - SQLite → AstraDB migration happens separately
   - Easy to re-seed AstraDB without re-running Agent 1

4. **Production-Ready:**
   - Minimal backend has zero dependencies on main backend
   - Agent 1 is a "data preparation" step, not runtime dependency

---

**STATUS:** ✅ **RECOMMENDED STRATEGY DEFINED**
**NEXT:** Update Agent 1 config + implement AstraDB seeding script








