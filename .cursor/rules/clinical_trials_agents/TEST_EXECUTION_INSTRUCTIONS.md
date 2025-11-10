# 🚀 TEST EXECUTION INSTRUCTIONS

**Status:** Code fix applied, backend restart required

---

## **🔧 FIX APPLIED**

✅ Updated `clinical_trial_search_service.py` to support `GEMINI_API_KEY`  
✅ Environment verified: `GEMINI_API_KEY` is set in `.env`  
✅ Code change confirmed: Service now checks both `GEMINI_API_KEY` and `GOOGLE_API_KEY`

---

## **⚠️ REQUIRED ACTION**

**Backend server must be restarted** to load the code fix.

The singleton service (`ClinicalTrialSearchService`) was initialized with old code when the server started. It needs to be recreated with the new code.

---

## **📋 RESTART INSTRUCTIONS**

### **Option 1: If running in terminal**
1. Stop server: Press `Ctrl+C` in the terminal running the backend
2. Restart:
   ```bash
   cd oncology-coPilot/oncology-backend-minimal
   python -m uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
   ```

### **Option 2: If running as background process**
```bash
# Find process
ps aux | grep uvicorn

# Kill it
kill <PID>

# Restart
cd oncology-coPilot/oncology-backend-minimal
python -m uvicorn api.main:app --reload --host 0.0.0.0 --port 8000 &
```

### **Option 3: If running with systemd/docker**
- Restart the service/container

---

## **🧪 AFTER RESTART - RUN TESTS**

### **Automated Test Script:**
```bash
cd oncology-coPilot/oncology-backend-minimal
./test_after_restart.sh
```

### **Manual Tests:**

**Test 1 - Hybrid Search:**
```bash
curl -X POST http://localhost:8000/api/trials/search-optimized \
  -H "Content-Type: application/json" \
  -d '{
    "query": "ovarian cancer BRCA1",
    "patient_context": {
      "condition": "ovarian cancer",
      "location_state": "NY"
    },
    "top_k": 10
  }'
```

**Test 2 - Autonomous Agent:**
```bash
curl -X POST http://localhost:8000/api/trials/agent/search \
  -H "Content-Type: application/json" \
  -d '{
    "mutations": [{"gene": "BRCA1", "hgvs_p": "V600E"}],
    "disease": "ovarian cancer",
    "state": "CA"
  }'
```

---

## **✅ EXPECTED RESULTS**

After restart, you should see:
- ✅ HTTP 200 responses
- ✅ `success: true` in JSON
- ✅ Trial results with `optimization_score` fields
- ✅ Graph-optimized results (site location, sponsor info)

---

## **🔍 VERIFICATION**

Check server logs for:
```
✅ ClinicalTrialSearchService initialized (collection: clinical_trials_eligibility)
✅ Neo4j connection established
```

If you see these, the fix is loaded and tests should work!






