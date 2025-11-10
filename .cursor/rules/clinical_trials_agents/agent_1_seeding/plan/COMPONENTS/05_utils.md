# 🛠️ MODULE 5: UTILS

## **Purpose**
Error handling logic and centralized logging

## **File Locations**
- `scripts/agent_1_seeding/utils/error_handler.py` - Hybrid error handling logic
- `scripts/agent_1_seeding/utils/logger.py` - Centralized logging setup

## **Key Functions**

### **`should_fail_on_api_errors(consecutive_failures: int, trials_fetched: int) -> bool`** (error_handler.py)
API failure threshold check:
- Returns `True` if should fail (critical failure)
- Checks: `consecutive_failures >= MAX_API_CONSECUTIVE_FAILURES` AND `trials_fetched < MIN_TRIALS_FETCHED`

### **`should_fail_on_chromadb_errors(failed_batches: int, total_batches: int) -> bool`** (error_handler.py)
ChromaDB failure threshold check:
- Returns `True` if should fail (critical failure)
- Checks: `failed_batches / total_batches > CHROMA_MAX_FAILURE_RATE`

### **`setup_logger() -> logging.Logger`** (logger.py)
Centralized logging setup with consistent format

## **Dependencies**
- `config.py` - Error thresholds (MAX_API_CONSECUTIVE_FAILURES, MIN_TRIALS_FETCHED, CHROMA_MAX_FAILURE_RATE)

## **Acceptance Criteria**
- [ ] Error handler functions work correctly
- [ ] Logging setup consistent across all modules
- [ ] Thresholds match config values

## **Time Estimate:** 30 minutes


