# Trial Tagging Infrastructure Improvements - For Plumber

**Date:** January 8, 2025  
**Status:** 📋 **TASK OUTLINE FOR PLUMBER**  
**Purpose:** Document infrastructure improvements needed for robust, production-ready trial tagging

---

## 🎯 **CURRENT STATE AUDIT**

### **What We Have:**
- ✅ Basic batch tagging script (`tag_trials_moa_batch.py`)
- ✅ OpenAI GPT-4o integration
- ✅ Basic rate limiting (1s delay between calls)
- ✅ Progress saving after each batch
- ✅ Automatic LLM confidence-based validation

### **What's Missing:**
- ⚠️ **Robust rate limiting** (exponential backoff, retry logic)
- ⚠️ **Error recovery** (resume from failures)
- ⚠️ **Monitoring/alerting** (track success rates, API costs)
- ⚠️ **Concurrency control** (parallel processing with rate limit respect)
- ⚠️ **Quality metrics** (confidence distribution, validation stats)
- ⚠️ **Cost tracking** (API usage, cost per trial)

---

## 🔧 **INFRASTRUCTURE IMPROVEMENTS NEEDED**

### **1. Enhanced Rate Limiting & Retry Logic** ⭐ **PRIORITY 1**

**Current:** Basic 1s delay, no retry logic for rate limits

**Needed:**
- ✅ **Exponential backoff** for 429 errors (2s, 4s, 8s)
- ✅ **Jitter** to prevent thundering herd
- ✅ **Max retries** configurable (default: 3)
- ✅ **Retry-specific errors** (429, 500, 502, 503, 504)
- ✅ **Circuit breaker** pattern (stop after N consecutive failures)

**Implementation:**
```python
# Exponential backoff with jitter
delay = (2 ** attempt) * INITIAL_BACKOFF
jitter = delay * 0.1 * random.random()
total_delay = delay + jitter
await asyncio.sleep(total_delay)
```

**Status:** ✅ **IMPLEMENTED** in current script

---

### **2. Resume from Failures** ⭐ **PRIORITY 2**

**Current:** Script starts from scratch if interrupted

**Needed:**
- ✅ **Checkpoint file** (save state after each trial)
- ✅ **Resume flag** (`--resume` to continue from last checkpoint)
- ✅ **Failed trials list** (retry only failed trials)
- ✅ **Progress tracking** (trials completed, remaining, ETA)

**Implementation:**
```python
# Checkpoint file: .trial_tagging_checkpoint.json
{
    "last_completed_nct_id": "NCT04284969",
    "completed_trials": ["NCT04284969", "NCT04001023", ...],
    "failed_trials": ["NCT12345678", ...],
    "started_at": "2025-01-08T10:00:00Z",
    "last_updated": "2025-01-08T10:15:00Z"
}
```

**Status:** ⚠️ **NOT IMPLEMENTED** - Needs plumber work

---

### **3. Concurrency Control** ⭐ **PRIORITY 3**

**Current:** Sequential processing (one trial at a time)

**Needed:**
- ✅ **Async semaphore** (limit concurrent API calls)
- ✅ **Rate limit token bucket** (respect API tier limits)
- ✅ **Parallel processing** (process N trials concurrently)
- ✅ **Dynamic throttling** (adjust concurrency based on rate limit errors)

**Implementation:**
```python
# Semaphore for concurrency control
semaphore = asyncio.Semaphore(max_concurrent=5)

async def tag_trial_with_semaphore(trial, semaphore):
    async with semaphore:
        return await tag_trial_with_openai(trial, api_key)

# Process in parallel
tasks = [tag_trial_with_semaphore(trial, semaphore) for trial in trials]
results = await asyncio.gather(*tasks, return_exceptions=True)
```

**Status:** ⚠️ **NOT IMPLEMENTED** - Needs plumber work

---

### **4. Monitoring & Metrics** ⭐ **PRIORITY 4**

**Current:** Basic logging only

**Needed:**
- ✅ **Success rate tracking** (trials tagged / total trials)
- ✅ **Confidence distribution** (histogram of confidence scores)
- ✅ **API call metrics** (total calls, retries, rate limit hits)
- ✅ **Time metrics** (avg time per trial, total elapsed)
- ✅ **Cost tracking** (estimated cost per trial, total cost)

**Implementation:**
```python
# Metrics tracking
metrics = {
    "total_trials": 100,
    "successful": 95,
    "failed": 5,
    "success_rate": 0.95,
    "confidence_distribution": {
        "high (>0.7)": 60,
        "medium (0.5-0.7)": 25,
        "low (<0.5)": 10
    },
    "api_calls": {
        "total": 100,
        "retries": 15,
        "rate_limit_hits": 3
    },
    "cost_estimate": {
        "per_trial": 0.02,  # $0.02 per trial (gpt-4o)
        "total": 2.00
    }
}
```

**Status:** ⚠️ **PARTIALLY IMPLEMENTED** - Basic metrics in logs, needs structured output

---

### **5. Quality Validation** ⭐ **PRIORITY 5**

**Current:** LLM confidence scores only

**Needed:**
- ✅ **Cross-validation** (tag same trial twice, compare results)
- ✅ **Consistency checks** (same intervention → same MoA vector)
- ✅ **Outlier detection** (unusual confidence scores flagged)
- ✅ **Validation report** (summary of tagging quality)

**Implementation:**
```python
# Cross-validation: tag 10% of trials twice
validation_trials = random.sample(trials, int(len(trials) * 0.1))
for trial in validation_trials:
    result1 = await tag_trial_with_openai(trial, api_key)
    await asyncio.sleep(5)  # Wait between calls
    result2 = await tag_trial_with_openai(trial, api_key)
    
    # Compare MoA vectors
    similarity = cosine_similarity(result1['moa_vector'], result2['moa_vector'])
    if similarity < 0.9:
        logger.warning(f"⚠️ Inconsistent tagging for {trial['nct_id']}: {similarity:.2f}")
```

**Status:** ⚠️ **NOT IMPLEMENTED** - Needs plumber work

---

### **6. Cost Optimization** ⭐ **PRIORITY 6**

**Current:** Uses gpt-4o (expensive but accurate)

**Needed:**
- ✅ **Model selection** (gpt-4o for high-confidence, gpt-4o-mini for low-confidence)
- ✅ **Caching** (cache similar trials to avoid duplicate API calls)
- ✅ **Batch API calls** (if OpenAI supports batch endpoint)
- ✅ **Cost tracking** (track actual vs estimated costs)

**Implementation:**
```python
# Model selection based on trial complexity
if is_complex_trial(trial):  # Multiple interventions, unclear MoA
    model = "gpt-4o"  # Use expensive model
else:
    model = "gpt-4o-mini"  # Use cheaper model

# Caching
cache_key = hashlib.md5(f"{trial['title']}{trial['interventions']}".encode()).hexdigest()
if cache_key in cache:
    return cache[cache_key]
```

**Status:** ⚠️ **NOT IMPLEMENTED** - Needs plumber work

---

## 📋 **PLUMBER TASK BREAKDOWN**

### **Task 1: Enhanced Rate Limiting** (2-3 hours)
- [ ] Implement exponential backoff with jitter
- [ ] Add retry logic for rate limit errors (429)
- [ ] Add circuit breaker pattern
- [ ] Test with rate limit simulation

**Files to Modify:**
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/tag_trials_moa_batch.py`

**Status:** ✅ **DONE** - Already implemented in current script

---

### **Task 2: Resume from Failures** (3-4 hours)
- [ ] Create checkpoint file format
- [ ] Save checkpoint after each trial
- [ ] Add `--resume` flag to resume from checkpoint
- [ ] Add `--retry-failed` flag to retry only failed trials
- [ ] Test resume functionality

**Files to Create/Modify:**
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/tag_trials_moa_batch.py`
- `.trial_tagging_checkpoint.json` (checkpoint file)

**Status:** ⚠️ **PENDING** - Needs plumber work

---

### **Task 3: Concurrency Control** (4-5 hours)
- [ ] Implement async semaphore for concurrency control
- [ ] Add rate limit token bucket
- [ ] Process trials in parallel (configurable concurrency)
- [ ] Dynamic throttling based on rate limit errors
- [ ] Test with various concurrency levels

**Files to Modify:**
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/tag_trials_moa_batch.py`

**Status:** ⚠️ **PENDING** - Needs plumber work

---

### **Task 4: Monitoring & Metrics** (2-3 hours)
- [ ] Add structured metrics tracking
- [ ] Generate metrics report (JSON + human-readable)
- [ ] Track success rate, confidence distribution, API calls
- [ ] Cost estimation and tracking
- [ ] Export metrics to file

**Files to Create/Modify:**
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/tag_trials_moa_batch.py`
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/metrics_reporter.py` (new)

**Status:** ⚠️ **PENDING** - Needs plumber work

---

### **Task 5: Quality Validation** (3-4 hours)
- [ ] Implement cross-validation (tag 10% of trials twice)
- [ ] Consistency checks (same intervention → same MoA)
- [ ] Outlier detection (unusual confidence scores)
- [ ] Generate validation report
- [ ] Flag inconsistent tags for review

**Files to Create/Modify:**
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/tag_trials_moa_batch.py`
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/validate_tagging_quality.py` (new)

**Status:** ⚠️ **PENDING** - Needs plumber work

---

### **Task 6: Cost Optimization** (2-3 hours)
- [ ] Model selection based on trial complexity
- [ ] Implement caching for similar trials
- [ ] Track actual vs estimated costs
- [ ] Generate cost report

**Files to Create/Modify:**
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/tag_trials_moa_batch.py`
- `oncology-coPilot/oncology-backend-minimal/scripts/trials/cost_tracker.py` (new)

**Status:** ⚠️ **PENDING** - Needs plumber work

---

## 🎯 **PRIORITY ORDER FOR PLUMBER**

1. **Task 1: Enhanced Rate Limiting** ✅ **DONE** (already implemented)
2. **Task 2: Resume from Failures** ⚠️ **HIGH PRIORITY** (enables long-running batches)
3. **Task 3: Concurrency Control** ⚠️ **MEDIUM PRIORITY** (speeds up processing)
4. **Task 4: Monitoring & Metrics** ⚠️ **MEDIUM PRIORITY** (visibility into process)
5. **Task 5: Quality Validation** ⚠️ **LOW PRIORITY** (nice to have)
6. **Task 6: Cost Optimization** ⚠️ **LOW PRIORITY** (nice to have)

---

## 📊 **EXPECTED IMPROVEMENTS**

### **Before (Current):**
- Sequential processing: 1 trial per second = 100 trials in ~100 seconds
- No retry logic: Failures are permanent
- No resume: Must restart from scratch if interrupted
- Basic logging: Limited visibility

### **After (With Improvements):**
- Parallel processing: 5 trials concurrently = 100 trials in ~20 seconds (5x faster)
- Robust retry logic: 95%+ success rate even with rate limits
- Resume capability: Can resume from any point
- Rich metrics: Full visibility into process, costs, quality

---

## 🔧 **ARCHITECTURE DIAGRAM**

```
┌─────────────────────────────────────────────────────────────┐
│                    Trial Tagging Pipeline                    │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 1. Load Untagged Trials (SQLite)                            │
│    - Filter by status (RECRUITING, ACTIVE)                   │
│    - Exclude already tagged                                  │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. Concurrency Control (Semaphore)                           │
│    - Limit concurrent API calls (default: 5)                 │
│    - Rate limit token bucket                                 │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. Tag Trial (OpenAI GPT-4o)                                │
│    - Extract interventions                                   │
│    - Call OpenAI API with retry logic                        │
│    - Parse MoA vector + confidence                           │
│    - Exponential backoff on rate limits                      │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. Save Progress (Checkpoint + JSON)                          │
│    - Save to trial_moa_vectors.json                          │
│    - Update checkpoint file                                  │
│    - Track metrics                                           │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. Generate Report                                           │
│    - Success rate, confidence distribution                  │
│    - API calls, retries, rate limit hits                    │
│    - Cost estimate                                           │
└─────────────────────────────────────────────────────────────┘
```

---

## 📝 **TESTING REQUIREMENTS**

### **Test 1: Rate Limiting**
- Simulate rate limit errors (429)
- Verify exponential backoff works
- Verify retry logic succeeds

### **Test 2: Resume from Failures**
- Interrupt script mid-run
- Resume from checkpoint
- Verify no duplicate tags

### **Test 3: Concurrency Control**
- Test with various concurrency levels (1, 5, 10)
- Verify rate limits are respected
- Verify no API errors from too many concurrent calls

### **Test 4: Quality Validation**
- Run cross-validation on 10% of trials
- Verify consistency (similarity > 0.9)
- Flag inconsistent tags

---

## 🎯 **SUCCESS CRITERIA**

- ✅ **95%+ success rate** (trials successfully tagged)
- ✅ **<5% rate limit errors** (with retry logic)
- ✅ **Resume capability** (can resume from any point)
- ✅ **5x speedup** (with concurrency control)
- ✅ **Full metrics** (success rate, confidence, costs)
- ✅ **Quality validation** (cross-validation consistency >0.9)

---

**Last Updated:** January 8, 2025  
**Status:** 📋 **READY FOR PLUMBER**  
**Estimated Total Time:** 16-22 hours (all tasks)

