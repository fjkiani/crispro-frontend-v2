# Deliverable: Programmatic Trial MoA Vector Tagging

**Date:** January 28, 2025  
**Status:** 📋 **READY FOR ASSIGNMENT**  
**Priority:** 🟡 **MEDIUM** (blocks mechanism fit ranking expansion)  
**Timeline:** 1-2 weeks  
**Assigned To:** **Trial Tagging Agent** (separate from SAE agent)

---

## 🎯 Executive Summary

**Objective:** Expand trial MoA vector coverage from 47 trials (3.4%) to 200+ trials (14%+) to enable mechanism-based trial matching for a broader set of trials.

**Current State:**
- ✅ 47 trials manually tagged (stored in `trial_moa_vectors.json`)
- ✅ MoA vector structure defined (7D: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
- ✅ Backend infrastructure ready (`extract_moa_vector_for_trial()` function exists)
- ⚠️ **Gap:** Need 153+ more trials tagged programmatically

**Target:**
- **Minimum:** 200 trials (14% coverage)
- **Ideal:** 500+ trials (36%+ coverage)
- **Manager Requirement:** 200+ trials with ≥90% accuracy (Manager P3)

---

## 📊 Current System Understanding

### **Trial Storage**

**Location:** `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db`

**Tables:**
- `trials` - Main trials table (1,397+ trials)
- `trials_fresh` - Freshly extracted recruiting trials

**Schema (relevant fields):**
```sql
CREATE TABLE trials (
    nct_id TEXT PRIMARY KEY,
    title TEXT,
    status TEXT,
    phase TEXT,
    description_text TEXT,
    inclusion_criteria_text TEXT,
    exclusion_criteria_text TEXT,
    objectives_text TEXT,
    metadata_json TEXT,  -- Contains interventions, conditions, etc.
    ...
);
```

### **MoA Vector Storage**

**Location:** `oncology-coPilot/oncology-backend-minimal/api/resources/trial_moa_vectors.json`

**Structure:**
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

### **Existing Functions**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/trial_data_enricher.py`

**Function:** `extract_moa_vector_for_trial(trial_data: Dict, use_gemini_tag: bool = False) -> Tuple[Dict, Dict]`

**Current Behavior:**
- If `use_gemini_tag=True`: Loads from `trial_moa_vectors.json` (if exists)
- If `use_gemini_tag=False`: Uses runtime keyword matching fallback
- Returns: `(moa_vector_dict, moa_metadata)`

**Gap:** No programmatic batch tagging function exists

---

## 🔧 Implementation Approach

### **Option 1: Gemini Batch Tagging (Manager P3 Compliant)** ⭐ **RECOMMENDED**

**Why:** Manager P3 requires Gemini tagging OFFLINE ONLY (never runtime)

**Process:**
1. **Extract Trial Data** from SQLite database
2. **Batch Tag with Gemini** (offline, not runtime)
3. **Human Spot-Review** 30 diverse trials (≥90% accuracy required)
4. **Store in `trial_moa_vectors.json`** with metadata
5. **Update Weekly** for new/changed trials

**Implementation Steps:**

#### **Step 1: Create Batch Tagging Script**

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/tag_trials_moa_batch.py`

**Functionality:**
```python
async def tag_trials_moa_batch(
    trial_nct_ids: List[str],
    batch_size: int = 50,
    use_gemini: bool = True
) -> Dict[str, Dict]:
    """
    Batch tag trials with MoA vectors using Gemini API (offline).
    
    Args:
        trial_nct_ids: List of NCT IDs to tag
        batch_size: Number of trials to process per batch
        use_gemini: Use Gemini API (True) or keyword matching (False)
    
    Returns:
        Dict mapping NCT ID to MoA vector data
    """
    # 1. Load trial data from SQLite
    # 2. For each trial:
    #    - Extract interventions from metadata_json
    #    - Call Gemini API (if use_gemini=True)
    #    - Parse MoA vector from Gemini response
    #    - Store with metadata
    # 3. Return tagged trials
```

**Gemini Prompt Template:**
```
Given the following clinical trial intervention information, 
determine the mechanism of action (MoA) vector.

Trial: {title}
Interventions: {interventions}
Conditions: {conditions}

Return a JSON object with MoA vector (7D):
{
  "ddr": 0.0-1.0,      // DNA Damage Repair
  "mapk": 0.0-1.0,     // RAS/MAPK pathway
  "pi3k": 0.0-1.0,     // PI3K/AKT pathway
  "vegf": 0.0-1.0,     // Angiogenesis
  "her2": 0.0-1.0,     // HER2 pathway
  "io": 0.0-1.0,       // Immunotherapy
  "efflux": 0.0-1.0    // Drug efflux
}

Confidence: 0.0-1.0
Primary MoA: "description"
```

#### **Step 2: Create Validation Script**

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validate_moa_tagging.py`

**Functionality:**
- Human spot-review 30 diverse trials
- Calculate accuracy (≥90% required)
- Generate validation report
- Flag uncertain tags for review

#### **Step 3: Create Update Script**

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/update_trial_moa_vectors.py`

**Functionality:**
- Weekly diff for new/changed trials
- Tag only new/changed trials
- Merge with existing `trial_moa_vectors.json`
- Preserve metadata (model, version, parsed_at, reviewed_by, source_checksum)

---

### **Option 2: Runtime Keyword Matching (Fallback Only)**

**Why:** If Gemini tagging fails or is unavailable

**Process:**
1. Extract interventions from trial data
2. Match keywords to MoA pathways
3. Use `DRUG_MECHANISM_DB` for drug→MoA mapping
4. Store with lower confidence (0.7 vs 0.95)

**Note:** Manager P3 says Gemini OFFLINE ONLY, but runtime fallback acceptable if Gemini tag missing

---

## 📋 Detailed Implementation Plan

### **Phase 1: Database Query & Extraction (2-3 hours)**

**Tasks:**
1. Create function to query SQLite for untagged trials
2. Extract trial data (NCT ID, title, interventions, conditions)
3. Filter to recruiting/active trials (priority)
4. Batch trials for processing

**Deliverable:**
- `get_untagged_trials(limit: int = 200) -> List[Dict]`
- Returns trials not in `trial_moa_vectors.json`

---

### **Phase 2: Gemini Batch Tagging (4-6 hours)**

**Tasks:**
1. Create Gemini API client (offline batch processing)
2. Create prompt template for MoA vector extraction
3. Batch process trials (50 at a time)
4. Parse Gemini responses to MoA vectors
5. Handle errors and retries

**Deliverable:**
- `tag_trials_moa_batch()` function
- Processes 200+ trials offline
- Stores results in `trial_moa_vectors.json`

**Manager P3 Compliance:**
- ✅ OFFLINE ONLY (never runtime)
- ✅ Batch tag 200 trials
- ✅ Human spot-review 30 diverse trials
- ✅ ≥90% accuracy required
- ✅ Metadata persistence (model, version, parsed_at, reviewed_by, source_checksum)

---

### **Phase 3: Validation & Review (2-3 hours)**

**Tasks:**
1. Create validation script
2. Select 30 diverse trials for human review
3. Calculate accuracy metrics
4. Flag uncertain tags
5. Generate validation report

**Deliverable:**
- `validate_moa_tagging.py` script
- Validation report showing ≥90% accuracy
- List of uncertain tags for review

---

### **Phase 4: Integration & Storage (1-2 hours)**

**Tasks:**
1. Merge tagged trials into `trial_moa_vectors.json`
2. Preserve existing 47 manually tagged trials
3. Add metadata (model, version, parsed_at, reviewed_by)
4. Update `extract_moa_vector_for_trial()` to use new tags

**Deliverable:**
- Updated `trial_moa_vectors.json` with 200+ trials
- Backward compatible with existing code
- Metadata tracking complete

---

### **Phase 5: Weekly Update Script (1-2 hours)**

**Tasks:**
1. Create script to diff new/changed trials
2. Tag only new/changed trials
3. Merge with existing tags
4. Preserve all metadata

**Deliverable:**
- `update_trial_moa_vectors.py` script
- Weekly update process documented

---

## 🎯 Success Criteria

**Deliverable Complete When:**
1. ✅ 200+ trials tagged with MoA vectors
2. ✅ Human spot-review completed (30 diverse trials)
3. ✅ Validation report shows ≥90% accuracy
4. ✅ All tags stored in `trial_moa_vectors.json` with metadata
5. ✅ Backward compatible with existing code
6. ✅ Weekly update script functional

**Validation Metrics:**
- **Coverage:** 200+ trials (14%+ of 1,397)
- **Accuracy:** ≥90% (human spot-review)
- **Confidence:** Average confidence ≥0.85
- **Metadata:** All tags have model, version, parsed_at, reviewed_by

---

## 📊 Expected Output

**Before:**
- 47 trials tagged (3.4% coverage)
- All manually tagged

**After:**
- 200+ trials tagged (14%+ coverage)
- 47 manually tagged (preserved)
- 153+ programmatically tagged (Gemini batch)
- All with metadata tracking

**File Structure:**
```json
{
  "NCT04284969": {
    "moa_vector": {...},
    "confidence": 0.95,
    "source": "manual_intelligence_report",
    "tagged_at": "2025-01-13T12:00:00Z",
    "reviewed_by": "Zo",
    ...
  },
  "NCT_NEW_001": {
    "moa_vector": {...},
    "confidence": 0.88,
    "source": "gemini_batch_tagging",
    "tagged_at": "2025-01-28T10:00:00Z",
    "reviewed_by": "TrialTaggingAgent",
    "provenance": {
      "model": "gemini-pro",
      "version": "v1",
      "parsed_at": "2025-01-28T10:00:00Z",
      "reviewed_by": "TrialTaggingAgent",
      "source_checksum": "..."
    }
  }
}
```

---

## 🔗 Integration Points

### **Backend Integration**

**File:** `api/services/trial_data_enricher.py`

**Function:** `extract_moa_vector_for_trial()`

**Current Behavior:**
- Loads from `trial_moa_vectors.json` (if exists)
- Falls back to runtime keyword matching

**After Tagging:**
- More trials will have Gemini tags
- Runtime fallback only for truly untagged trials
- Better mechanism fit accuracy

---

### **Mechanism Fit Ranking**

**File:** `api/services/mechanism_fit_ranker.py`

**Impact:**
- More trials will have MoA vectors
- Mechanism fit ranking works for 200+ trials (vs 47)
- Better trial matching coverage

---

## 🚨 Critical Requirements

### **Manager P3 Compliance (MUST HAVE)**

1. **OFFLINE ONLY:** Never use Gemini in runtime paths
2. **Batch Tagging:** Tag 200 trials offline
3. **Human Review:** Spot-review 30 diverse trials
4. **Accuracy:** ≥90% tag accuracy required
5. **Metadata:** Store model, version, parsed_at, reviewed_by, source_checksum
6. **Update Cadence:** Weekly diff for new/changed trials
7. **Uncertain Tags:** Default to neutral vector; never force a mechanism label

---

## 📝 Deliverables Checklist

- [ ] **Phase 1:** Database query & extraction script
- [ ] **Phase 2:** Gemini batch tagging script
- [ ] **Phase 3:** Validation & review script
- [ ] **Phase 4:** Integration & storage (merge into JSON)
- [ ] **Phase 5:** Weekly update script
- [ ] **Validation:** Human spot-review 30 trials (≥90% accuracy)
- [ ] **Documentation:** Tagging process documented
- [ ] **Testing:** Verify backward compatibility

---

## 🔗 Related Documents

- **Manager Policy:** `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md` (P3)
- **Implementation Review:** `MECHANISM_TRIAL_MATCHING_IMPLEMENTATION_REVIEW.md` (Gap 2)
- **Trial Storage:** `scripts/seed_trials_table.py`
- **MoA Extraction:** `api/services/trial_data_enricher.py`
- **MoA Storage:** `api/resources/trial_moa_vectors.json`

---

## ✅ Assignment Summary

**Assigned To:** **Trial Tagging Agent** (separate from SAE agent)

**Objective:** Expand trial MoA vector coverage from 47 → 200+ trials

**Timeline:** 1-2 weeks

**Priority:** 🟡 MEDIUM (blocks mechanism fit ranking expansion)

**Dependencies:**
- ✅ SQLite database with trials
- ✅ Gemini API access (for batch tagging)
- ✅ Existing `trial_moa_vectors.json` structure

**Blocking:**
- ⚠️ Mechanism fit ranking expansion (needs 200+ trials)
- ⚠️ Publication validation (needs broader trial coverage)

---

*Deliverable Created: January 28, 2025*  
*Status: 📋 READY FOR ASSIGNMENT*  
*Assigned To: Trial Tagging Agent*


