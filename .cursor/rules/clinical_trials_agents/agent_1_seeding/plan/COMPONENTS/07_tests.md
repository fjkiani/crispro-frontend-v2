# 🧪 MODULE 7: TESTS

## **Purpose**
Comprehensive test suite for all modules

## **File Locations**
- `tests/agent_1_seeding/test_api_client.py` - Test API fetcher
- `tests/agent_1_seeding/test_parsers.py` - Test parsers
- `tests/agent_1_seeding/test_database.py` - Test migration, SQLite, ChromaDB
- `tests/agent_1_seeding/test_integration.py` - End-to-end tests
- `tests/agent_1_seeding/conftest.py` - Pytest fixtures (test mode flag support)

## **Test Cases**

### **test_api_client.py**
- `test_api_fetch_small()` - Test API connectivity (10 trials)
- `test_duplicate_detection()` - Verify no duplicate NCT IDs
- `test_rate_limiting()` - Verify rate limiting works

### **test_parsers.py**
- `test_parse_study()` - Test full study parsing
- `test_biomarker_extraction()` - Test biomarker keywords
- `test_locations_parsing()` - Test locations extraction
- `test_error_handling()` - Test missing data handling

### **test_database.py**
- `test_schema_migration()` - Verify new columns exist
- `test_migration_idempotency()` - Verify migration can run twice
- `test_sqlite_batch_commits()` - Verify batch commits work
- `test_chromadb_embedding()` - Test ChromaDB embedding (skip in test mode)

### **test_integration.py**
- `test_full_pipeline_test_mode()` - End-to-end with --test-mode
- `test_full_pipeline_full()` - End-to-end with 1000 trials (skipped by default)

## **Test Mode Support**
- `conftest.py` provides `pytest.mark.test_mode` marker
- Tests check `pytestconfig.getoption("--test-mode")` to skip ChromaDB tests
- Conditional skip logic (Option B)

## **Dependencies**
- All modules (for integration tests)

## **Acceptance Criteria**
- [ ] All unit tests pass
- [ ] Integration tests pass (test mode)
- [ ] Test mode properly skips ChromaDB tests
- [ ] Tests can run independently or together

## **Time Estimate:** 1 hour (30 min unit, 30 min integration)









