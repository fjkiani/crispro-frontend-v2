# ⚔️ AYESHA TRIALS SEARCH - VECTOR SUPPORT ISSUE RESOLUTION

**Date**: January 13, 2025  
**Status**: ✅ **RESOLVED** - New collection created with vector support, ready for seeding  
**Collection**: `clinical_trials_eligibility2` (768 dimensions, cosine similarity)

---

## 🔍 **HOW WE DISCOVERED THE ISSUE**

### **Initial Problem**
- Backend `/api/ayesha/trials/search` endpoint returning **0 trials** despite 200+ trials seeded
- Vector search in AstraDB returning **0 results** (no candidates found)

### **Investigation Timeline**

#### **Phase 1: Backend Search Logic Bug**
**Finding**: `HybridTrialSearchService` was incorrectly processing `ClinicalTrialSearchService.search_trials()` response
- **Bug**: Treated dictionary response as a list
- **Fix**: Correctly extract `found_trials` from nested `data` key
- **Result**: Still 0 results (vector search not working)

#### **Phase 2: Vector Field Missing**
**Finding**: AstraDB documents missing `$vector` field
- **Root Cause**: `seed_astradb_from_sqlite.py` using `find_one_and_update` with `update={"$set": document}`
- **Issue**: `$vector` field must be at root level, not nested in `$set`
- **Fix**: Changed to `replace_one` with full document (including `$vector` at root)
- **Result**: Still 0 results (field still not persisting)

#### **Phase 3: Collection Configuration**
**Finding**: Collection `clinical_trials_eligibility` was created **without vector dimensions**
- **Root Cause**: Collection created manually or via script that didn't specify vector dimensions
- **Impact**: AstraDB silently ignores `$vector` field when collection doesn't support vectors
- **Evidence**: 
  - Documents inserted/updated via `astrapy` client did not save `$vector` field
  - `find()` queries did not return `$vector` field
  - Vector search returned 0 results

#### **Phase 4: User Confirmation**
**User Action**: Created new collection `clinical_trials_eligibility2` with:
- **Vector dimensions**: 768 (for Google `text-embedding-004`)
- **Similarity metric**: Cosine
- **Embedding method**: "Bring my own" (we generate embeddings)

**User Observation**: 
- Old collection (`clinical_trials_eligibility`): Some documents showed `{$binary}` for `$vector` (likely manually inserted or from different seeding process)
- New collection (`clinical_trials_eligibility2`): Properly configured with 768 dimensions

---

## 🔧 **FIXES APPLIED**

### **1. Seeding Script Fix (`seed_astradb_from_sqlite.py`)**

**Before**:
```python
collection.find_one_and_update(
    filter={"_id": trial_dict['source_url']},
    update={"$set": document},  # ❌ $vector not at root level
    upsert=True
)
```

**After**:
```python
# Set _id explicitly for upsert
document["_id"] = trial_dict['source_url']
collection.replace_one(
    filter={"_id": trial_dict['source_url']},
    replacement=document,  # ✅ Full document with $vector at root
    upsert=True
)
```

**Key Changes**:
- ✅ Use `replace_one` instead of `find_one_and_update`
- ✅ Include `$vector` at root level of document (not nested in `$set`)
- ✅ Verify embedding generation before upsert (debug logging)
- ✅ Verify `$vector` presence in document before upsert

### **2. String Truncation Fix**

**Issue**: `SHRED_DOC_LIMIT_VIOLATION` (document size > 8000 bytes)

**Fix**: Added `truncate_string_to_bytes()` helper function
```python
def truncate_string_to_bytes(text: str, max_bytes: int = 7500) -> str:
    """Truncate string to fit within byte limit (safe margin for 8000-byte limit)"""
    if not text:
        return text
    
    encoded = text.encode('utf-8')
    if len(encoded) <= max_bytes:
        return text
    
    # Truncate and decode safely
    truncated = encoded[:max_bytes]
    return truncated.decode('utf-8', errors='ignore')
```

**Applied to**:
- `eligibility_text` (truncated to 7500 bytes)
- `description_text` (truncated to 7500 bytes)

### **3. Environment Configuration**

**Updated `.env`**:
```bash
ASTRA_COLLECTION_NAME=clinical_trials_eligibility2  # ✅ New collection with vector support
```

**Backend Restart**: Fully restarted backend to load new environment variable

### **4. Hard Filter Bug Fixes**

**Bug 1: Location Field Mismatch**
- **Issue**: `_apply_ayesha_hard_filters` checking `trial.get("locations", [])`
- **Fix**: Changed to `trial.get("locations_data", [])` (matches AstraDB field)
- **File**: `api/routers/ayesha_trials.py`

**Bug 2: Disease Category Mismatch**
- **Issue**: Search using `disease_category: "ovarian_cancer"`
- **Fix**: Changed to `disease_category: "gynecologic_oncology"` (matches AstraDB values)
- **File**: `api/routers/ayesha_trials.py`

---

## ✅ **CURRENT STATUS**

### **Collection Configuration**
- **Name**: `clinical_trials_eligibility2`
- **Vector Dimensions**: 768 (Google `text-embedding-004`)
- **Similarity Metric**: Cosine
- **Current Documents**: 230 (from initial seeding)

### **Backend Configuration**
- **Environment Variable**: `ASTRA_COLLECTION_NAME=clinical_trials_eligibility2`
- **Backend Status**: ✅ Operational, using new collection
- **Search Endpoint**: `/api/ayesha/trials/search` ✅ Working

### **Test Results**
- **Ayesha Search**: ✅ **1 trial found** (NCT06819007)
  - Passes all hard filters (Recruiting, First-line, NYC metro, Stage IV)
  - Match score: 0.85
  - Phase: PHASE3

---

## 📊 **SEEDING STATISTICS**

### **Initial Seeding (200 trials)**
- **Source**: SQLite `clinical_trials.db`
- **Target**: AstraDB `clinical_trials_eligibility2`
- **Status**: ✅ Complete (230 documents - includes duplicates/updates)
- **Coverage**:
  - ✅ Disease category: 100% (all tagged as `gynecologic_oncology`)
  - ✅ Locations data: 100% (all have `locations_data` array)
  - ⚠️ PI names: 0% (needs investigation - ClinicalTrials.gov API v2 structure)
  - ✅ Embeddings: 100% (all documents have `$vector` field)

### **Next Seeding (500 trials)**
- **Status**: 🔄 **READY TO EXECUTE**
- **Script**: `scripts/seed_astradb_from_sqlite.py`
- **Collection**: `clinical_trials_eligibility2`
- **Expected**: 500 new documents with vector embeddings

---

## 🎯 **LESSONS LEARNED**

### **1. Collection Configuration is Critical**
- **Lesson**: AstraDB collections must be created with vector dimensions from the start
- **Impact**: Documents inserted before vector support cannot be updated with vectors
- **Solution**: Always specify vector dimensions when creating collections

### **2. `$vector` Field Must Be at Root Level**
- **Lesson**: `find_one_and_update` with `$set` does not work for `$vector` field
- **Impact**: Vector field silently ignored if nested in `$set`
- **Solution**: Use `replace_one` with full document (including `$vector` at root)

### **3. Document Size Limits**
- **Lesson**: AstraDB has 8000-byte document size limit
- **Impact**: Large text fields (`eligibility_text`, `description_text`) can exceed limit
- **Solution**: Truncate strings to 7500 bytes (safe margin)

### **4. Field Name Consistency**
- **Lesson**: Backend code must match AstraDB document structure
- **Impact**: Hard filters fail if field names don't match
- **Solution**: Always verify field names in actual documents (not just schema)

### **5. Disease Category Tagging**
- **Lesson**: Use consistent disease category values across system
- **Impact**: Search filters fail if categories don't match
- **Solution**: Standardize on `gynecologic_oncology` (not `ovarian_cancer`)

---

## 📋 **NEXT STEPS**

### **Immediate (Today)**
1. ✅ **Documentation Complete** (this file)
2. 🔄 **Seed 500 More Trials** (execute `seed_astradb_from_sqlite.py` with limit=500)
3. ⏸️ **Verify Seeding** (check document count, vector presence, search functionality)

### **Short-Term (This Week)**
1. **Investigate PI Name Coverage** (0% - needs ClinicalTrials.gov API v2 structure analysis)
2. **Expand Disease Categories** (add more cancer types beyond gynecologic)
3. **Test Search Performance** (with 700+ documents)

### **Long-Term (Next Sprint)**
1. **Automated Collection Creation** (script to create collections with vector support)
2. **Seeding Pipeline** (automated re-seeding when SQLite updates)
3. **Vector Search Optimization** (tune similarity thresholds, top_k values)

---

## 📁 **RELATED DOCUMENTATION**

- **Vector Support Issue**: `.cursor/ayesha/AYESHA_TRIALS_SEARCH_VECTOR_SUPPORT_ISSUE.md`
- **Debug Complete**: `.cursor/ayesha/AYESHA_TRIALS_SEARCH_DEBUG_COMPLETE.md`
- **Root Cause Found**: `.cursor/ayesha/AYESHA_TRIALS_SEARCH_ROOT_CAUSE_FOUND.md`
- **Agent Missions**: `.cursor/ayesha/02_AGENT_MISSIONS_CONSOLIDATED.md`

---

**DOCTRINE STATUS: ACTIVE** ⚔️  
**LAST UPDATED**: January 13, 2025  
**NEXT UPDATE**: After 500-trial seeding completes


