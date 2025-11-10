# ⚔️ IMPLEMENTATION PLAN - STEP BY STEP

## **Implementation Order (Follow Module Dependencies)**

**Dependency Chain:** Config → Utils → API → Parsers → Database → Main → SQL

---

## **PHASE 1: Setup Folder Structure (5 minutes)**

```bash
cd oncology-coPilot/oncology-backend
mkdir -p scripts/agent_1_seeding/{api,parsers,database,utils}
mkdir -p tests/agent_1_seeding
touch scripts/agent_1_seeding/{__init__,main,config}.py
touch scripts/agent_1_seeding/api/{__init__,ctgov_client}.py
touch scripts/agent_1_seeding/parsers/{__init__,study_parser,biomarker_extractor,locations_parser}.py
touch scripts/agent_1_seeding/database/{__init__,migration,sqlite_client,chromadb_client}.py
touch scripts/agent_1_seeding/utils/{__init__,error_handler,logger}.py
touch scripts/migrate_schema_v2.sql
touch tests/agent_1_seeding/{__init__,test_api_client,test_parsers,test_database,test_integration,conftest}.py
```

---

## **PHASE 2: Build Modules (Follow Dependency Order)**

### **Step 1: Config Module (15 min)**
**File:** `scripts/agent_1_seeding/config.py`
- Define all constants from `COMPONENTS/01_config.md`
- Test: Import and verify constants accessible

### **Step 2: Utils Module (30 min)**
**Files:** `utils/logger.py`, `utils/error_handler.py`
- Setup logging
- Implement error handling functions
- Dependencies: Config
- Test: Unit tests for error handler logic

### **Step 3: API Client (1 hour)**
**File:** `api/ctgov_client.py`
- Implement `fetch_ovarian_trials()` with deduplication
- Add error handling integration
- Dependencies: Config, Utils
- Test: Mock API responses, verify deduplication

### **Step 4: Parsers (1 hour)**
**Files:** `parsers/*.py`
- Implement study parser (orchestrator)
- Implement biomarker extractor
- Implement locations parser
- Dependencies: Config
- Test: Mock API study objects, verify parsing

### **Step 5: Database Module (1 hour)**
**Files:** `database/*.py`
- Implement idempotent migration (checks columns first)
- Implement SQLite client with batch commits
- Implement ChromaDB client with rate limiting + retries
- Dependencies: Config, Utils
- Test: Mock database operations, verify idempotency

### **Step 6: Main CLI (30 min)**
**File:** `main.py`
- Implement CLI argument parser
- Implement orchestration function
- Implement summary report generator
- Dependencies: All modules
- Test: End-to-end with mock data

### **Step 7: Migration SQL (15 min)**
**File:** `scripts/migrate_schema_v2.sql`
- Write ALTER TABLE statements
- Write CREATE INDEX statements
- Test: Run migration, verify columns

**Total Build Time: ~4 hours**

---

## **PHASE 3: Testing (1 hour)**

### **Unit Tests (30 min)**
- Test each module independently
- Use mocks for external dependencies

### **Integration Tests (30 min)**
- Test full pipeline with test mode
- Verify end-to-end flow

---

## **PHASE 4: Execution (20 minutes)**

See `EXECUTION/checklist.md` and `EXECUTION/commands.md`

---

## **TOTAL TIME: ~5 hours 20 minutes**

