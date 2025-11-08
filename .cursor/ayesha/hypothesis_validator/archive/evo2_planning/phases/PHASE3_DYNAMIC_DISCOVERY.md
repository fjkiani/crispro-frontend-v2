# ⚔️ PHASE 3: DYNAMIC COMPOUND DISCOVERY

**Time Estimate:** 1.5 hours  
**Dependencies:** None (can run in parallel with Phase 1)  
**Status:** Ready to start

---

## **🎯 GOAL**

Allow validation of ANY compound, not just hardcoded list. Extract targets from literature or ChEMBL.

**Core Innovation:** Dynamic target discovery means we can validate compounds not in our database.

---

## **📁 FILES TO CREATE**

- `oncology-coPilot/oncology-backend-minimal/api/services/compound_target_extraction.py` (200+ lines)

---

## **🔧 IMPLEMENTATION**

See [`../services/compound_target_extraction_service.md`](../services/compound_target_extraction_service.md) for complete implementation code.

### **Key Components:**

1. **`CompoundTargetExtractor`** class
   - `extract_targets()` - Main method with priority fallback chain
   - `_check_knowledge_base()` - Check food_targets.json
   - `_query_chembl()` - Query ChEMBL API for targets
   - `_extract_from_literature()` - LLM + PubMed extraction

2. **Priority Chain:**
   ```
   [1] Check cache (if previously extracted)
       ↓
   [2] Check knowledge base (food_targets.json) - FASTEST
       ↓
   [3] Query ChEMBL API - MODERATE SPEED
       ↓
   [4] Extract from PubMed literature (LLM) - SLOWEST
   ```

3. **Output Schema:**
   ```python
   {
       "compound": "Resveratrol",
       "targets": ["SIRT1", "mTOR", "AMPK"],
       "pathways": ["Longevity", "Metabolic regulation"],
       "source": "chembl" | "literature" | "knowledge_base",
       "confidence": 0.75
   }
   ```

---

## **🔗 INTEGRATIONS**

### **Dependencies:**
- `api/services/llm_literature_service.py` (existing) - for literature extraction
- `food_targets.json` (existing) - knowledge base fallback

### **External APIs:**
- ChEMBL REST API: `https://www.ebi.ac.uk/chembl/api/data/`
- PubMed (via LLM service)

---

## **🧪 TESTING**

### **Test Case 1: Known Compound (Knowledge Base)**
```bash
curl -X POST http://127.0.0.1:8000/api/test/extract_targets \
  -H "Content-Type: application/json" \
  -d '{"compound": "Vitamin D", "disease": "ovarian cancer"}'
```

**Expected:** Returns targets from `food_targets.json`, source="knowledge_base"

### **Test Case 2: Unknown Compound (ChEMBL)**
```bash
curl -X POST http://127.0.0.1:8000/api/test/extract_targets \
  -H "Content-Type: application/json" \
  -d '{"compound": "Resveratrol", "disease": "ovarian cancer"}'
```

**Expected:** Queries ChEMBL, extracts SIRT1/mTOR/AMPK, source="chembl"

### **Test Case 3: Very Unknown (Literature)**
```bash
curl -X POST http://127.0.0.1:8000/api/test/extract_targets \
  -H "Content-Type: application/json" \
  -d '{"compound": "Quercetin", "disease": "ovarian cancer"}'
```

**Expected:** Uses LLM to extract from PubMed abstracts, source="literature"

---

## **✅ ACCEPTANCE CRITERIA**

- [ ] Knowledge base check works (food_targets.json)
- [ ] ChEMBL API integration functional (with error handling)
- [ ] LLM literature extraction works (when LLM service available)
- [ ] Fallback chain correct (KB → ChEMBL → Literature)
- [ ] Caching implemented (avoid repeated ChEMBL/literature calls)
- [ ] Error handling graceful (returns empty targets if all fail)
- [ ] Confidence scoring reflects source quality

---

## **📝 NOTES**

- **ChEMBL API:** Free, rate-limited, may be slow
- **LLM Extraction:** Requires Pubmed-LLM-Agent to be available
- **Caching:** Cache extracted targets to avoid repeated API calls
- **Confidence:** knowledge_base > chembl > literature (in that order)

---

**Next Phase:** [`PHASE4_ENHANCED_ENDPOINT.md`](./PHASE4_ENHANCED_ENDPOINT.md) (integrates all 3 services)

