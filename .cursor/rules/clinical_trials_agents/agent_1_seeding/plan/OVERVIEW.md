# 🤖 AGENT 1: DATA SEEDING AGENT - OVERVIEW 📥

## **⚔️ MISSION**
Bulk-load 1000 ovarian cancer trials from ClinicalTrials.gov API v2 into our database with proper disease hierarchy and embeddings.

---

## **🎯 OBJECTIVES**

### **Primary Goal:**
Populate empty SQLite database with 1000 fresh, recruiting ovarian cancer trials, each tagged with disease category, biomarkers, and location data.

### **Success Criteria:**
- ✅ 1000+ trials inserted into SQLite
- ✅ All trials have `disease_category = "gynecologic_oncology"`
- ✅ All trials have `disease_subcategory = "ovarian_cancer"`
- ✅ ChromaDB contains 1000+ embeddings
- ✅ All trials have `locations_data` JSON with NY locations
- ✅ Script completes in <15 minutes
- ✅ 4/4 unit tests pass

---

## **✅ FINAL DECISIONS (COMMANDER APPROVED)**

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **SQLite Batch Size** | 100 | Balanced (10 commits for 1000 trials) |
| **ChromaDB Batch Size** | 50 | Safe (50 embeddings/min, within API limit) |
| **ChromaDB Retries** | 3 with exponential backoff | Resilience for transient failures |
| **Migration Idempotency** | Check columns first | Prevents crashes on re-run |
| **Error Handling** | Hybrid (fail if >50% fail) | Balance between resilience and data integrity |
| **Test Mode** | Conditional skip in tests | Single test suite, flexible execution |
| **File Paths** | Migration: `scripts/`, DB: relative from root | Clear, maintainable structure |

---

## **📁 FOLDER STRUCTURE**

```
agent_1_seeding/
├── plan/
│   ├── OVERVIEW.md                    # This file - mission, objectives, decisions
│   ├── COMPONENTS/                    # Component-specific documentation
│   │   ├── 01_config.md               # Configuration module
│   │   ├── 02_api_client.md           # ClinicalTrials.gov API fetcher
│   │   ├── 03_parsers.md              # Study/biomarker/locations parsers
│   │   ├── 04_database.md             # Migration, SQLite, ChromaDB
│   │   ├── 05_utils.md                # Error handling, logging
│   │   ├── 06_main_cli.md             # CLI orchestration
│   │   └── 07_tests.md                # Test suite
│   ├── IMPLEMENTATION/
│   │   ├── step_by_step.md            # Implementation order and dependencies
│   │   └── time_estimates.md          # Detailed time breakdown
│   └── EXECUTION/
│       ├── checklist.md               # Pre-flight, execution, verification
│       └── commands.md                # All bash commands
└── AGENT_1_DOCTRINE.md                # Master index (points to all above)
```

---

## **🔍 INFRASTRUCTURE AUDIT SUMMARY**

### **✅ EXISTING ASSETS (95% READY!)**
- SQLite database at `oncology-backend/backend/data/clinical_trials.db`
- ChromaDB with Google Embeddings (production-ready)
- ClinicalTrials.gov API v2 utils
- Frontend integration ready

### **⚠️ GAPS TO CLOSE**
- 5 schema columns missing (disease_category, disease_subcategory, biomarker_requirements, locations_data, last_updated)
- Dedicated ovarian cancer bulk fetcher
- Biomarker keyword extraction
- Locations data JSON formatting

---

## **⚔️ MODULAR ARCHITECTURE**

**7 Modules to Build:**
1. **Config** - Centralized constants
2. **API Client** - ClinicalTrials.gov fetcher with deduplication
3. **Parsers** - Study/biomarker/locations extraction
4. **Database** - Migration, SQLite, ChromaDB
5. **Utils** - Error handling, logging
6. **Main CLI** - Orchestration and CLI arguments
7. **Tests** - Unit and integration tests

**Implementation Order:** Config → Utils → API → Parsers → Database → Main → SQL

---

## **📊 TIME ESTIMATE**

**Total: ~4 hours 20 minutes**
- Schema migration: 30 min
- Code development: 3 hours (7 modules)
- Testing: 30 min
- Smoke test: 5 min
- Full execution: 15 min

---

## **🚀 QUICK START**

1. **Read:** `COMPONENTS/` docs for each module specification
2. **Follow:** `IMPLEMENTATION/step_by_step.md` for build order
3. **Execute:** `EXECUTION/checklist.md` for commands
4. **Verify:** `COMPONENTS/07_tests.md` for test suite

---

## **📚 DOCUMENTATION INDEX**

- **Component Specs:** `COMPONENTS/` (7 files)
- **Build Plan:** `IMPLEMENTATION/step_by_step.md`
- **Execution:** `EXECUTION/checklist.md` + `EXECUTION/commands.md`
- **Full Doctrine:** `AGENT_1_DOCTRINE.md` (detailed reference)

---

**STATUS:** ✅ Ready to build - All decisions finalized, modular structure defined

