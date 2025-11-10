# 🔧 MODULE 1: CONFIGURATION

## **Purpose**
Centralized constants and configuration for Agent 1

## **File Location**
`oncology-coPilot/oncology-backend/scripts/agent_1_seeding/config.py`

## **Constants to Define**

```python
# API Configuration
CTGOV_API_URL = "https://clinicaltrials.gov/api/v2/studies"
CTGOV_RATE_LIMIT = 0.5  # seconds between requests (2 req/sec)
CTGOV_TIMEOUT = 30  # seconds

# Database Configuration
SQLITE_BATCH_SIZE = 100  # Commit every 100 trials (Zo's recommendation)
SQLITE_DB_PATH = "backend/data/clinical_trials.db"  # Relative from oncology-backend/

# ChromaDB Configuration
CHROMA_BATCH_SIZE = 50  # Embed 50 trials per batch
CHROMA_RATE_LIMIT = 1.0  # seconds between batches (50 embeddings/min)
CHROMA_MAX_RETRIES = 3  # Retry failed batches up to 3 times
CHROMA_PATH = "backend/data/chroma_data"  # Relative from oncology-backend/
CHROMA_COLLECTION = "clinical_trials_eligibility"

# Error Handling Configuration
MAX_API_CONSECUTIVE_FAILURES = 5  # Fail if 5 consecutive API failures
MIN_TRIALS_FETCHED = 10  # Fail if <10 trials fetched
CHROMA_MAX_FAILURE_RATE = 0.5  # Fail if >50% ChromaDB batches fail

# Biomarker Keywords
BIOMARKER_KEYWORDS = {
    "BRCA1": ["BRCA1", "BRCA 1"],
    "BRCA2": ["BRCA2", "BRCA 2"],
    "HRD": ["HRD", "homologous recombination deficiency", "HR deficiency"],
    "TP53": ["TP53", "p53"],
    "CCNE1": ["CCNE1", "cyclin E1"],
    "MYC": ["MYC", "c-MYC"]
}

# Disease Classification
DISEASE_CATEGORY = "gynecologic_oncology"
DISEASE_SUBCATEGORY = "ovarian_cancer"

# Migration SQL Path
MIGRATION_SQL_PATH = "scripts/migrate_schema_v2.sql"
```

## **Acceptance Criteria**
- [ ] All constants defined
- [ ] Paths are relative from `oncology-backend/` root
- [ ] Batch sizes match decisions (100 for SQLite, 50 for ChromaDB)
- [ ] Can be imported by all other modules

## **Time Estimate:** 15 minutes

