# 🌐 MODULE 2: API CLIENT

## **Purpose**
Fetch trials from ClinicalTrials.gov API v2 with duplicate detection and error handling

## **File Location**
`oncology-coPilot/oncology-backend/scripts/agent_1_seeding/api/ctgov_client.py`

## **Key Functions**

### **`fetch_ovarian_trials(limit: int) -> List[Dict]`**
Main fetcher function with:
- Pagination handling (page tokens)
- Duplicate detection by NCT ID
- Progress logging inside loop
- Hybrid error handling (fail if >5 consecutive failures)
- Rate limiting (2 req/sec)

### **`_deduplicate_by_nct_id(trials: List[Dict], new_studies: List[Dict]) -> List[Dict]`**
Helper to deduplicate studies by NCT ID

### **`_handle_api_error(error: Exception, consecutive_failures: int) -> bool`**
Error handling with fail thresholds (from utils/error_handler.py)

## **API Parameters**

```python
params = {
    "query.cond": "ovarian cancer OR peritoneal cancer OR primary peritoneal carcinomatosis",
    "filter.overallStatus": "RECRUITING,NOT_YET_RECRUITING",
    "filter.phase": "PHASE2,PHASE3,PHASE4",
    "filter.geo": "distance(United States, 2000mi)",
    "pageSize": 100,
    "format": "json"
}
```

## **Features**
- ✅ Duplicate detection by NCT ID
- ✅ Progress logging: `"Fetched {len(trials)}/{limit} trials (deduped: {len(new_studies)} new)"`
- ✅ Hybrid error handling (fail if >5 consecutive API failures)
- ✅ Rate limiting: `await asyncio.sleep(CTGOV_RATE_LIMIT)`

## **Dependencies**
- `config.py` - API URL, rate limit, timeout, error thresholds
- `utils/error_handler.py` - `should_fail_on_api_errors()`
- `utils/logger.py` - Logging setup

## **Acceptance Criteria**
- [ ] Fetches 1000 trials in <10 minutes
- [ ] Proper error handling for API failures
- [ ] Rate limiting respected (2 req/sec)
- [ ] No duplicate NCT IDs in final results
- [ ] Progress logged every page

## **Time Estimate:** 1 hour









