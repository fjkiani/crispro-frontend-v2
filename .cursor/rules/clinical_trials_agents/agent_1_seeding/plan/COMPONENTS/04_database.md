# 💾 MODULE 4: DATABASE

## **Purpose**
Schema migration, SQLite insertion, ChromaDB embedding

## **File Locations**
- `scripts/agent_1_seeding/database/migration.py` - Idempotent schema migration
- `scripts/agent_1_seeding/database/sqlite_client.py` - Batch SQLite insertion
- `scripts/agent_1_seeding/database/chromadb_client.py` - Rate-limited ChromaDB embedding

## **Key Functions**

### **`run_migration_if_needed(db_path: str, migration_sql_path: str) -> None`** (migration.py)
Idempotent schema migration:
- Checks columns before ALTER TABLE (Option B)
- Only runs migration if columns missing
- Adds 5 new columns + 3 indexes

### **`insert_trials_batched(trials: List[Dict], batch_size: int = 100) -> None`** (sqlite_client.py)
Batch SQLite insertion:
- Commits every 100 trials (from config)
- Progress logging every batch
- Uses `INSERT OR REPLACE` for idempotency

### **`embed_trials_batched(trials: List[Dict], batch_size: int = 50) -> None`** (chromadb_client.py)
Rate-limited ChromaDB embedding:
- Batches of 50 (from config)
- 1 second sleep between batches (50 embeddings/min)
- 3 retries with exponential backoff per batch
- Hybrid error handling (fail if >50% batches fail)

## **Features**
- ✅ Idempotent migration (checks columns before ALTER TABLE)
- ✅ Batch commits every 100 trials
- ✅ Rate limiting (50 embeddings/min) with 3 retries
- ✅ Hybrid error handling (fail if >50% ChromaDB batches fail)

## **Dependencies**
- `config.py` - Batch sizes, paths, retry settings
- `utils/error_handler.py` - `should_fail_on_chromadb_errors()`
- `utils/logger.py` - Logging setup

## **Acceptance Criteria**
- [ ] Migration idempotent (can run multiple times)
- [ ] SQLite batch commits work (100 trials per commit)
- [ ] ChromaDB rate limiting respected (50/min)
- [ ] Retry logic works for transient failures
- [ ] Error handling fails if >50% batches fail

## **Time Estimate:** 1 hour


