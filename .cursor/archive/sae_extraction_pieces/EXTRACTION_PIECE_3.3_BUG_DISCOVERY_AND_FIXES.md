# 📖 EXTRACTION PIECE 3.3: Bug Discovery and Fixes
**Source**: Lines 20000-20500 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents three critical bugs discovered during SAE service deployment and execution, their root causes, fixes implemented, and resolution process. These bugs prevented the SAE service from working correctly.

---

## 🔍 KEY FINDINGS

### **BUG #1: Modal SAE Using Old Code (ModuleList Attribute Error)**

**Error Message:**
```
ModuleList has no attribute `26`
ERROR: SAE extraction failed: ModuleList has no attribute `26`
```

**Root Cause:**
- Modal was still running the **old container** with incorrect code
- Old code used: `self.evo_model.model.blocks[26]` (integer indexing)
- Correct code should use: `self.evo_model.model.blocks._modules['26']` (string key access)
- Modal has a **warm container cache** - even after redeployment, old container kept serving requests

**Why It Happened:**
- Modal keeps warm containers alive as long as they're getting traffic
- Backend + cohort script kept hitting the endpoint → old container never went idle → never died → new code never loaded
- Cold-start time showed `0ms` = warm container = OLD CODE

**Fix Applied:**
1. **Code Fix**: Changed layer access from `blocks[26]` to `blocks._modules['26']`
2. **Deployment Fix**: Redeployed SAE service with corrected code
3. **Container Fix**: Required manual stop of old container via Modal dashboard OR wait 5+ minutes for auto-scale-down

**Resolution Steps:**
1. Stop local backend (halt traffic to Modal)
2. Stop cohort script (prevent new requests)
3. Manually stop Modal app via dashboard OR wait 5 minutes
4. Restart backend and re-run → Modal cold-starts with NEW code

**Key Insight**: Modal's warm container cache requires manual intervention or idle wait to force code updates.

---

### **BUG #2: SAE Weights Not Loading (Checkpoint Format Mismatch)**

**Error Message:**
```
WARNING: Could not load SAE weights: Error(s) in loading state_dict for BatchTopKTiedSAE:
Missing key(s) in state_dict: "W", "b_enc", "b_dec". 
Unexpected key(s) in state_dict: "_orig_mod.W", "_orig_mod.b_enc", "_orig_mod.b_dec".
Using random initialization.
```

**Root Cause:**
- Downloaded SAE checkpoint had `_orig_mod.` prefixes on all keys
- Prefixes come from `torch.compile` when model was saved
- Our `BatchTopKTiedSAE` model expects bare keys: `W`, `b_enc`, `b_dec`
- Key mismatch → `load_state_dict()` failed → model initialized with random weights

**Why It Happened:**
- Hugging Face checkpoint was saved with `torch.compile` enabled
- `torch.compile` adds `_orig_mod.` prefix to all model parameters
- Our loading code didn't account for this prefix

**Fix Applied:**

**Code Change** (`src/services/sae_service/main.py`):
```python
checkpoint = torch.load(sae_weights_path, map_location=device, weights_only=False)

# Strip "_orig_mod." prefix if present (from torch.compile)
if any(k.startswith("_orig_mod.") for k in checkpoint.keys()):
    logger.info("Stripping '_orig_mod.' prefix from checkpoint keys...")
    checkpoint = {k.replace("_orig_mod.", ""): v for k, v in checkpoint.items()}

self.sae_model.load_state_dict(checkpoint)
```

**Result:**
- Checkpoint keys stripped of `_orig_mod.` prefix before loading
- Model loads successfully with correct weights
- SAE features now computed correctly

**Key Insight**: Hugging Face checkpoints may have `torch.compile` prefixes that need stripping.

---

### **BUG #3: Modal App Keeps Running (Warm Container Issue)**

**User Concern:**
- Modal app keeps running indefinitely
- Old container persists even after redeployment
- No auto-stop mechanism

**Root Cause:**
- Modal's architecture keeps warm containers alive as long as they receive traffic
- `scaledown_window=300` (5 minutes) only triggers if container is idle
- Continuous requests prevent idle state → container never scales down

**Current Safeguards (Already in Place):**
- ✅ API key required (blocks unauthorized traffic)
- ✅ Circuit breaker (stops accepting requests after 50% error rate)
- ✅ Local cohort script has caps: `MAX_PATIENTS`, `MAX_TOTAL_VARIANTS`, `ENABLE_SAE_COHORT_RUN`

**What's Still Missing:**
- ❌ Modal app doesn't auto-terminate when idle
- ❌ No timeout-based shutdown

**Real Fix:**
- Use Modal's `scaledown_window` (already set to 300s = 5min)
- **Manual stop required**: Stop app via dashboard when done
- **Alternative**: Wait 5+ minutes without requests → auto-scale-down → next request cold-starts

**Key Insight**: Modal requires manual intervention or idle wait to stop containers. Cost controls prevent runaway costs, but don't auto-stop.

---

### **ADDITIONAL ISSUE: API Key Mismatch**

**Error Message:**
```json
{"detail":"SAE service error: {\"detail\":\"Invalid API key\"}"}
```

**Root Cause:**
- Backend sending: `X-API-Key: dev-unsafe-run`
- Modal service expecting different key (default or env var mismatch)

**Fix:**
- Set same `SAE_API_KEY` value in both:
  - Backend environment: `SAE_API_KEY=some-long-temp-key`
  - Modal container env: `@app.cls(env={"SAE_API_KEY": "some-long-temp-key"})`
- Redeploy Modal service
- Restart backend with matching key

**Key Insight**: API key must match between backend and Modal service.

---

## 📊 KEY INSIGHTS

### **Modal Container Lifecycle**

1. **Warm Containers**: Stay alive as long as receiving traffic
2. **Cold Start**: Only happens when container is idle for `scaledown_window` duration
3. **Code Updates**: Require cold start to take effect
4. **Manual Stop**: Required to force immediate code update

### **Checkpoint Loading Patterns**

1. **torch.compile Prefixes**: Common in Hugging Face checkpoints
2. **Key Stripping**: Necessary before `load_state_dict()`
3. **Detection**: Check for `_orig_mod.` prefix in checkpoint keys
4. **Handling**: Strip prefix before loading

### **Deployment Best Practices**

1. **Stop Traffic First**: Halt backend/scripts before redeployment
2. **Force Cold Start**: Stop Modal app manually or wait for idle
3. **Verify Fixes**: Test with single request before full cohort
4. **Match Keys**: Ensure API keys match between services

---

## 🔗 CONTEXT & CONNECTIONS

- **Builds on**: Deployment instructions (Piece 2.5), Real data extraction (Piece 3.2)
- **Blocks**: Real SAE feature extraction until fixed
- **Resolved**: All three bugs fixed, service operational
- **Key Insight**: Modal deployment requires understanding container lifecycle

---

## 📝 NOTES

- All three bugs were discovered during real deployment
- Fixes were implemented and verified
- Modal warm container cache is a common deployment challenge
- Checkpoint format issues are common with `torch.compile`
- Manual intervention required for Modal container management

---

## 🎯 QUESTIONS RESOLVED

- ✅ Why is Modal using old code? → Warm container cache, requires cold start
- ✅ Why aren't SAE weights loading? → `_orig_mod.` prefix from `torch.compile`
- ✅ Why does Modal keep running? → Warm containers stay alive with traffic
- ✅ How to force code update? → Stop Modal app manually or wait 5+ minutes idle
- ✅ How to fix checkpoint loading? → Strip `_orig_mod.` prefix before loading

