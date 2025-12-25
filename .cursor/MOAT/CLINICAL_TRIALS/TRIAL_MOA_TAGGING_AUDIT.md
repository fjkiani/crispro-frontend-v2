# 🔍 Trial MoA Tagging Deliverable - Comprehensive Audit

**Date:** January 28, 2025  
**Auditor:** Auto  
**Status:** ✅ **AUDIT COMPLETE**  
**Deliverable:** `.cursor/MOAT/CLINICAL_TRIALS/TRIAL_MOA_TAGGING_DELIVERABLE.md`

---

## 📊 Executive Summary

**Current State (VERIFIED):**
- ✅ **Database:** 1,397 trials in `trials` table (SQLite)
- ✅ **Tagged Trials:** **47 trials** in `trial_moa_vectors.json`
  - 5 trials: `manual_intelligence_report` (high confidence, 0.90-0.95)
  - 42 trials: `manual_keyword_matching` (lower confidence, 0.75-0.85)
- ✅ **File Location:** `oncology-coPilot/oncology-backend-minimal/api/resources/trial_moa_vectors.json`
- ✅ **Existing Script:** `scripts/tag_top_ovarian_trials_with_moa.py` (keyword-based tagging)
- ⚠️ **Gap:** `extract_moa_vector_for_trial()` doesn't load from JSON file

**Key Findings:**
1. ✅ Database structure verified (1,397 trials)
2. ✅ MoA vector JSON file exists and has proper structure
3. ⚠️ **CRITICAL GAP:** `extract_moa_vector_for_trial()` function doesn't actually load from JSON file
4. ✅ Gemini API usage patterns exist in codebase (can be reused)
5. ✅ Existing keyword-based tagging script can be enhanced

---

## 🔍 Detailed Findings

### 1. Database Verification ✅

**Location:** `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db`

**Tables Found:**
- `trials` - **1,397 trials** ✅
- `clinical_trials` - (smaller table, ~30 trials)
- `trials_fresh` - Freshly extracted recruiting trials

**Schema (Verified):**
```sql
CREATE TABLE trials (
    id TEXT PRIMARY KEY,           -- NCT ID
    title TEXT,
    status TEXT,
    phases TEXT,
    interventions TEXT,
    conditions TEXT,
    ...
)
```

**Status:** ✅ Database exists and has sufficient trials (1,397 vs 1,397 mentioned in deliverable)

---

### 2. MoA Vector Storage ✅

**Location:** `oncology-coPilot/oncology-backend-minimal/api/resources/trial_moa_vectors.json`

**Structure (Verified):**
```json
{
  "NCT04284969": {
    "moa_vector": {
      "ddr": 0.95,
      "mapk": 0.0,
      "pi3k": 0.0,
      "vegf": 0.0,
      "her2": 0.0,
      "io": 0.0,
      "efflux": 0.0
    },
    "confidence": 0.95,
    "source": "manual_intelligence_report",
    "tagged_at": "2025-01-13T12:00:00Z",
    "reviewed_by": "Zo",
    "provenance": {
      "intelligence_report": "...",
      "extraction_method": "manual_review",
      "primary_moa": "DNA Damage Repair (PARP + ATR inhibitors)"
    }
  }
}
```

**Sources Found:**
- `manual_intelligence_report` - High confidence (0.95), manually reviewed
- `manual_keyword_matching` - Lower confidence (0.75-0.85), keyword-based

**Status:** ✅ File exists with proper structure, 7D vector format correct

---

### 3. Existing Code Review ⚠️

#### 3.1 `extract_moa_vector_for_trial()` Function

**Location:** `oncology-coPilot/oncology-backend-minimal/api/services/trial_data_enricher.py:268`

**Current Implementation:**
```python
def extract_moa_vector_for_trial(
    trial_data: Dict[str, Any],
    use_gemini_tag: bool = False
) -> Tuple[Optional[Dict[str, float]], Dict[str, Any]]:
    # Priority 1: Check for pre-tagged Gemini vectors (if available)
    if use_gemini_tag:
        gemini_vector = trial_data.get('gemini_moa_vector')
        if gemini_vector:
            return gemini_vector.get('vector'), metadata
    
    # Priority 2: Runtime keyword matching (fallback only)
    # ... keyword matching logic ...
```

**⚠️ CRITICAL GAP FOUND:**
- Function does **NOT** load from `trial_moa_vectors.json` file
- Only checks for `gemini_moa_vector` in `trial_data` dict (which is never populated)
- Falls back to runtime keyword matching

**Expected Behavior (Per Deliverable):**
- Should load from `trial_moa_vectors.json` when `use_gemini_tag=True`
- Should check file first, then fall back to keyword matching

**Fix Required:**
```python
# Add this to extract_moa_vector_for_trial():
if use_gemini_tag:
    # Load from JSON file first
    nct_id = trial_data.get('nct_id') or trial_data.get('id')
    if nct_id and nct_id in TRIAL_MOA_VECTORS:
        moa_data = TRIAL_MOA_VECTORS[nct_id]
        metadata['source'] = moa_data.get('source', 'gemini_offline_tagging')
        metadata['version'] = moa_data.get('provenance', {}).get('version', '1.0')
        metadata['parsed_at'] = moa_data.get('tagged_at')
        metadata['reviewed_by'] = moa_data.get('reviewed_by')
        return moa_data.get('moa_vector'), metadata
```

**Status:** ⚠️ **GAP IDENTIFIED** - Function needs update to load from JSON file

---

#### 3.2 `TRIAL_MOA_VECTORS` Loading

**Location:** `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py:44-58`

**Current Implementation:**
```python
TRIAL_MOA_VECTORS = {}
try:
    moa_vectors_path = os.path.join(
        os.path.dirname(__file__), 
        "../resources/trial_moa_vectors.json"
    )
    if os.path.exists(moa_vectors_path):
        with open(moa_vectors_path, "r") as f:
            TRIAL_MOA_VECTORS = json.load(f)
        logger.info(f"✅ Loaded {len(TRIAL_MOA_VECTORS)} trial MoA vectors")
except Exception as e:
    logger.error(f"❌ Failed to load trial MoA vectors: {e}")
    TRIAL_MOA_VECTORS = {}
```

**Status:** ✅ **GOOD** - File is loaded at module initialization in `ayesha_trials.py`

**Issue:** This is only in `ayesha_trials.py`, not in `trial_data_enricher.py` where `extract_moa_vector_for_trial()` is defined

**Recommendation:** Move `TRIAL_MOA_VECTORS` loading to `trial_data_enricher.py` or create a shared module

---

#### 3.3 Existing Tagging Script

**Location:** `scripts/tag_top_ovarian_trials_with_moa.py`

**What It Does:**
- ✅ Queries SQLite database for trials
- ✅ Uses keyword matching (not Gemini)
- ✅ Saves to `trial_moa_vectors.json`
- ✅ Merges with existing vectors

**Limitations:**
- ❌ No Gemini API integration
- ❌ Keyword matching only (lower confidence)
- ❌ No batch processing with rate limiting

**Status:** ✅ **USEFUL REFERENCE** - Can be enhanced for Gemini batch tagging

---

### 4. Gemini API Usage Patterns ✅

**Found Examples:**
1. `api/services/enhanced_evidence_service.py:856` - `_call_gemini_llm()`
   - Uses `google.generativeai` directly
   - Model: `gemini-2.0-flash-exp`
   - API key: `GEMINI_API_KEY` env var

2. `api/services/comprehensive_analysis/llm_explanation_enhancer.py:437`
   - Uses `genai.GenerativeModel()`
   - Model: `gemini-2.0-flash-exp`
   - Async support with `asyncio.to_thread()`

3. `api/services/trial_intelligence/pipeline.py:163`
   - Rate limiting: 30s between calls (free tier = 2/min)
   - Error handling for quota exceeded

**Pattern to Reuse:**
```python
import google.generativeai as genai
import os

genai.configure(api_key=os.getenv("GEMINI_API_KEY"))
model = genai.GenerativeModel("gemini-2.0-flash-exp")
response = model.generate_content(prompt)
```

**Status:** ✅ **PATTERNS FOUND** - Can reuse existing Gemini integration code

---

## 🚨 Critical Issues Found

### Issue 1: `extract_moa_vector_for_trial()` Doesn't Load JSON File ⚠️ **CRITICAL**

**Problem:**
- Function signature says it uses Gemini tags, but doesn't actually load from `trial_moa_vectors.json`
- Only checks for `gemini_moa_vector` in trial_data dict (never populated)

**Impact:**
- Pre-tagged trials in JSON file are **NOT being used**
- All trials fall back to runtime keyword matching (lower accuracy)

**Fix Required:**
1. Load `TRIAL_MOA_VECTORS` in `trial_data_enricher.py`
2. Check JSON file first when `use_gemini_tag=True`
3. Fall back to keyword matching only if not found

---

### Issue 2: Vector Format Mismatch ⚠️ **MEDIUM**

**Problem:**
- JSON file uses dict format: `{"ddr": 0.95, "mapk": 0.0, ...}`
- Some code expects list format: `[0.95, 0.0, ...]`

**Found In:**
- `trial_matching_agent.py:217` - Converts dict to list using `convert_moa_dict_to_vector()`

**Status:** ✅ **HANDLED** - Conversion function exists, but needs to be used consistently

---

## ✅ What's Working Well

1. ✅ **Database Structure** - 1,397 trials available, proper schema
2. ✅ **JSON File Structure** - Proper format with metadata, provenance tracking
3. ✅ **Existing Script** - Keyword-based tagging works, can be enhanced
4. ✅ **Gemini Integration** - Patterns exist in codebase, can be reused
5. ✅ **Loading in ayesha_trials.py** - File is loaded at module init

---

## 📋 Recommendations

### Priority 1: Fix `extract_moa_vector_for_trial()` Function ⚠️ **CRITICAL**

**Action:**
1. Load `TRIAL_MOA_VECTORS` in `trial_data_enricher.py` (or create shared module)
2. Update function to check JSON file first
3. Preserve existing keyword matching as fallback

**Code Changes:**
```python
# In trial_data_enricher.py
import json
import os
from pathlib import Path

# Load at module level
TRIAL_MOA_VECTORS = {}
try:
    moa_vectors_path = Path(__file__).parent.parent / "resources" / "trial_moa_vectors.json"
    if moa_vectors_path.exists():
        with open(moa_vectors_path, "r") as f:
            TRIAL_MOA_VECTORS = json.load(f)
        logger.info(f"✅ Loaded {len(TRIAL_MOA_VECTORS)} trial MoA vectors")
except Exception as e:
    logger.error(f"❌ Failed to load trial MoA vectors: {e}")

# Update extract_moa_vector_for_trial():
def extract_moa_vector_for_trial(...):
    # Priority 1: Check JSON file
    if use_gemini_tag:
        nct_id = trial_data.get('nct_id') or trial_data.get('id')
        if nct_id and nct_id in TRIAL_MOA_VECTORS:
            moa_data = TRIAL_MOA_VECTORS[nct_id]
            metadata.update({
                'source': moa_data.get('source', 'gemini_offline_tagging'),
                'version': moa_data.get('provenance', {}).get('version', '1.0'),
                'parsed_at': moa_data.get('tagged_at'),
                'reviewed_by': moa_data.get('reviewed_by')
            })
            return moa_data.get('moa_vector'), metadata
    
    # Priority 2: Runtime keyword matching (fallback)
    # ... existing keyword matching code ...
```

---

### Priority 2: Create Batch Tagging Script ✅ **ALIGNED WITH DELIVERABLE**

**Action:**
1. Create `scripts/tag_trials_moa_batch.py`
2. Use Gemini API for batch tagging (OFFLINE ONLY)
3. Implement rate limiting (30s between calls)
4. Save to JSON file with proper metadata

**Implementation:**
- Reuse Gemini patterns from `enhanced_evidence_service.py`
- Query database for untagged trials
- Batch process with rate limiting
- Store with Manager P3 compliant metadata

---

### Priority 3: Create Validation Script ✅ **ALIGNED WITH DELIVERABLE**

**Action:**
1. Create `scripts/validate_moa_tagging.py`
2. Select 30 diverse trials for human review
3. Calculate accuracy metrics
4. Generate validation report

---

### Priority 4: Create Update Script ✅ **ALIGNED WITH DELIVERABLE**

**Action:**
1. Create `scripts/update_trial_moa_vectors.py`
2. Diff new/changed trials weekly
3. Tag only new/changed trials
4. Merge with existing tags

---

## 🎯 Questions for Clarification

1. ✅ **Current Tag Count:** **47 trials** (5 manual intelligence, 42 keyword matching) - **VERIFIED**
2. **Gemini API Key:** Is `GEMINI_API_KEY` environment variable set and accessible?
3. **Priority Trials:** Should we prioritize recruiting/active trials first?
4. **Validation Process:** Who will perform the 30-trial human spot-review?
5. **Update Cadence:** Is weekly update cadence acceptable, or should it be more frequent?

---

## 📊 Implementation Readiness Assessment

| Component | Status | Notes |
|-----------|--------|-------|
| Database Access | ✅ Ready | 1,397 trials available |
| JSON File Structure | ✅ Ready | Proper format, metadata tracking |
| Gemini API Integration | ✅ Ready | Patterns exist, can reuse |
| Batch Tagging Script | ⏳ Pending | Needs to be created |
| Validation Script | ⏳ Pending | Needs to be created |
| Update Script | ⏳ Pending | Needs to be created |
| `extract_moa_vector_for_trial()` Fix | ⚠️ **CRITICAL** | Must fix before batch tagging |

---

## ✅ Next Steps

1. **IMMEDIATE:** Fix `extract_moa_vector_for_trial()` to load from JSON file
2. **Phase 1:** Create batch tagging script with Gemini API
3. **Phase 2:** Create validation script
4. **Phase 3:** Create update script
5. **Phase 4:** Human spot-review 30 trials
6. **Phase 5:** Merge and deploy

---

**Audit Complete:** January 28, 2025  
**Status:** ✅ Ready for implementation (after fixing critical gap)

