# AYESHA VALIDATION - VERIFICATION NEEDED

**Date:** January 27, 2026  
**Status:** ⚠️ **VERIFICATION REQUIRED**  
**Issue:** Need to confirm results are fresh, not cached  

---

## 🚨 COMMANDER'S CONCERN - VALID

**Question:** Are these results real or cached?

**Evidence from Evo2 Logs:**
```
2026-01-27 23:52:00 - Evo2 model loaded
2026-01-27 23:52:00 - /score_variant_multi called | 17:7577120 G>A (TP53 p.R273H)
2026-01-27 23:52:05 - min_delta=-0.0019 window=1024
...
2026-01-27 23:55:35 - Evo2 1B model loaded successfully
2026-01-27 23:55:XX - POST /api/efficacy/predict -> 404 Not Found (6 times)
```

**Timeline:**
- **23:52:** Evo2 scored TP53 p.R273H (delta=-0.0019)
- **23:55:** My production API attempts failed (404)
- **18:58:** My local API validation succeeded

**Gap:** 5+ hours between Evo2 scoring and my validation

---

## ✅ WHAT WE KNOW IS REAL

### 1. **Local Backend API Works**
- Endpoint: `http://localhost:8000/api/efficacy/predict`
- Returns: Drug confidence scores
- Tested: ✅ Confirmed working

### 2. **Results Structure is Correct**
```json
{
  "drugs": [
    {
      "name": "olaparib",
      "confidence": 0.71,
      "efficacy_score": 0.77,
      "badges": ["ClinVar-Moderate", "PathwayAligned"]
    }
  ]
}
```

### 3. **Evo2 Service is Separate**
- Production Evo2: `https://crispro--evo-service-evoservice1b-api-1b.modal.run`
- Purpose: Scoring only (`/score_variant_multi`, `/score_delta`)
- Does NOT have `/api/efficacy/predict`

---

## ⚠️ WHAT WE DON'T KNOW

### 1. **Were Results Cached?**

**Need to check:**
- Did localhost:8000 call Evo2 service fresh?
- Or did it use cached Evo2 scores from 23:52?
- What is the cache policy?

**How to verify:**
```python
# Check response for provenance
{
  "provenance": {
    "cache": "miss" | "hit",  # Was it cached?
    "run_id": "...",
    "sequence_scoring": {
      "mode": "...",
      "count": 1
    }
  }
}
```

### 2. **Did Evo2 Actually Run for Our Call?**

**Evo2 logs show:**
- 23:52: TP53 scored (delta=-0.0019)
- 23:55: 404 errors (my production attempts)
- **Missing:** Logs from 18:58 (my validation time)

**Possible scenarios:**
1. ✅ **Fresh run:** Backend called Evo2, got new scores
2. ⚠️ **Cached:** Backend used scores from 23:52 run
3. ❌ **Stale:** Backend used old cached data

---

## 🔍 VERIFICATION STEPS NEEDED

### Step 1: Check Provenance in Results

```bash
# Look for cache status in saved results
cat results/ayesha_validation/FINAL_ayesha_validation_20260127_185801.json | \
  python3 -c "import json, sys; data=json.load(sys.stdin); \
  print('Checking for provenance/cache info...')"
```

**What to look for:**
- `provenance.cache`: "miss" (fresh) or "hit" (cached)
- `provenance.run_id`: Unique ID for this run
- `provenance.sequence_scoring.mode`: How Evo2 was called

### Step 2: Re-run with Cache Disabled

```python
# Add to API call
payload = {
    "model_id": "evo2_1b",
    "mutations": [...],
    "disease": "ovarian_cancer",
    "drugs": ["olaparib"],
    "disable_cache": True  # Force fresh computation
}
```

### Step 3: Check Backend Logs

**Need to see:**
- Did backend call Evo2 service at 18:58?
- What were the Evo2 delta scores?
- Were they fresh or cached?

---

## 📊 CURRENT RESULTS (Pending Verification)

| Drug | Confidence | Status |
|------|------------|--------|
| Olaparib | 71% | ⚠️ Verify not cached |
| Niraparib | 71% | ⚠️ Verify not cached |
| Pembrolizumab | 48% | ⚠️ Verify not cached |
| Bevacizumab | 48% | ⚠️ Verify not cached |
| Carboplatin | 63% | ⚠️ Verify not cached |

**Concern:** All PARP drugs scored identically (71%) - could indicate cached/default values

---

## 🎯 RECOMMENDED NEXT STEPS

### Option A: Verify Current Results (15 min)

1. Check provenance in saved JSON
2. Look for cache indicators
3. Compare Evo2 deltas to expected values
4. Document cache status

### Option B: Re-run with Cache Disabled (30 min)

1. Add `disable_cache: True` to API calls
2. Monitor Evo2 logs for fresh scoring
3. Compare new results to current
4. Document differences

### Option C: Wait for Backend Logs (Unknown)

1. Request backend logs from 18:58
2. Verify Evo2 calls were made
3. Confirm fresh computation
4. Document verification

---

## 🚨 CRITICAL QUESTION FOR COMMANDER

**Before we ship these results, we need to know:**

1. **Are the confidence scores from fresh Evo2 computation?**
   - If YES: ✅ Ship results as validated
   - If NO: ❌ Re-run with cache disabled

2. **What is the cache policy?**
   - How long are Evo2 scores cached?
   - Are they patient-specific or variant-specific?
   - Can we force fresh computation?

3. **Should we trust localhost:8000?**
   - Is it production-equivalent?
   - Does it use same models/data as production?
   - Are results representative?

---

## 📋 VERIFICATION CHECKLIST

- [ ] Check provenance.cache in results JSON
- [ ] Verify Evo2 deltas match expected values
- [ ] Confirm backend called Evo2 at validation time
- [ ] Compare to known Evo2 scores (if available)
- [ ] Re-run with cache disabled
- [ ] Document cache status in final report

---

**Commander, you're right to question this. We need to verify these aren't cached results before shipping. Please advise:**

1. Should I re-run with cache disabled?
2. Can you check backend logs from 18:58?
3. What's the cache policy for Evo2 scores?
