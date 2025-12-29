# ✅ Therapy Fit Demo Script - Verification

**Date:** January 2025  
**Status:** ✅ **VERIFIED - Uses Real API Calls**

---

## 🔍 Verification: No Hard-Coded Values

The demo script (`scripts/demo_therapy_fit.py`) **makes real HTTP API calls** and displays actual responses. Here's how to verify:

### 1. **API Call Verification**

The script makes real HTTP POST requests:

```python
# Line 198 in demo_therapy_fit.py
response = await client.post(url, json=payload, timeout=TIMEOUT)
```

- ✅ Uses `httpx.AsyncClient` for real HTTP calls
- ✅ Calls actual endpoint: `http://127.0.0.1:8000/api/efficacy/predict`
- ✅ Sends real mutation data in payload
- ✅ Waits for actual API response
- ✅ Times the response (shows real latency)

### 2. **Response Data Extraction**

All displayed data comes from the API response:

```python
# Line 246: Extract drugs from response
drugs = response.get("drugs", [])

# Line 104-107: Display actual values from response
name = drug.get("name", "N/A")  # From API response
efficacy = drug.get("efficacy_score", 0.0)  # From API response
confidence = drug.get("confidence", 0.0)  # From API response
tier = drug.get("evidence_tier", "N/A")  # From API response
```

### 3. **How to Verify It's Real**

#### Option A: Run with Debug Mode
```bash
python scripts/demo_therapy_fit.py --test-case AYESHA-001 --debug
```

This will show:
- Response structure (keys in response)
- First drug's actual keys
- Sample of first drug data from API

#### Option B: Check API Logs
When you run the script, check your backend server logs. You should see:
```
INFO: POST /api/efficacy/predict
INFO: Processing mutations: [{'gene': 'MBD4', ...}]
```

#### Option C: Network Inspection
Use a network monitor (like `tcpdump` or browser dev tools) to see the actual HTTP request/response.

#### Option D: Modify Test Case
Change a mutation in the test case - the results should change:
```python
# In TEST_CASES, change MBD4 to a different gene
{"gene": "BRCA1", ...}  # Results will be different
```

### 4. **What Happens If API Fails**

If the API is not running or fails, the script will:
- ❌ Show error message: "API call failed"
- ❌ Display HTTP error code
- ❌ Show response text (first 500 chars)
- ❌ **NOT show any fake/hard-coded data**

### 5. **Response Time Verification**

The script measures actual response time:
```python
start_time = datetime.now()
response = await call_efficacy_predict(...)
elapsed = (datetime.now() - start_time).total_seconds()
```

If it were hard-coded, the response would be instant (< 0.01s). Real API calls take 2-5 minutes for Evo2 scoring.

---

## 🧪 Test It Yourself

### Step 1: Start Backend
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload
```

### Step 2: Run Demo with Debug
```bash
python scripts/demo_therapy_fit.py --test-case AYESHA-001 --debug
```

### Step 3: Verify Output
You should see:
1. "📡 Calling Therapy Fit API..." (real HTTP call)
2. "✅ Response received (XX.XXs)" (actual response time)
3. "🔍 Response Structure (debug):" (shows actual response keys)
4. Drug names, scores, confidence from actual API

### Step 4: Stop Backend and Run Again
```bash
# Stop backend (Ctrl+C)
# Run demo again
python scripts/demo_therapy_fit.py --test-case AYESHA-001
```

You should see:
- ❌ "API call failed"
- ❌ Error message
- ❌ **NO fake data displayed**

---

## 📊 Data Flow

```
Test Case (mutations)
    ↓
HTTP POST to /api/efficacy/predict
    ↓
Backend processes (Evo2, insights, etc.)
    ↓
API Response (JSON)
    ↓
Script extracts data from response
    ↓
Formats and displays (NO hard-coding)
```

---

## ✅ Confirmation

**The demo script:**
- ✅ Makes real HTTP API calls
- ✅ Displays actual API responses
- ✅ Shows real response times
- ✅ Handles API failures gracefully
- ✅ **NO hard-coded drug names**
- ✅ **NO hard-coded scores**
- ✅ **NO hard-coded confidence values**
- ✅ **NO mock data**

All displayed values come directly from the `/api/efficacy/predict` endpoint response.

---

*Verified: January 2025*  
*Script: `scripts/demo_therapy_fit.py`*

