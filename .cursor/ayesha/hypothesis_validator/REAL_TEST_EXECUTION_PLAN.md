# ⚔️ REAL TEST EXECUTION PLAN - NO MOCKS, ACTUAL DATA FLOW

**Commander:** Alpha  
**Architect:** Zo  
**Date:** December 2024  
**Objective:** Test with REAL backends, capture REAL input/output, verify data flow

---

## 🎯 TESTING PHILOSOPHY

**NOT:**
- ❌ Pass/fail boolean tests
- ❌ Mock data
- ❌ Hardcoded expected values
- ❌ Unit tests in isolation

**YES:**
- ✅ Real API calls to running backend
- ✅ Capture actual input → output for each stage
- ✅ Verify data structure and content
- ✅ Show Commander what the system actually produces
- ✅ Modular tests (one component at a time)

---

## 📊 TEST SUITE ARCHITECTURE

### **MODULAR APPROACH:**

```
Real Tests/
├── 01_backend_services/       # Test backend APIs directly
│   ├── test_food_validator.py
│   ├── test_unified_care.py
│   ├── test_sae_features.py
│   └── test_evidence_mining.py
│
├── 02_data_extraction/        # Test data sources
│   ├── test_chembl_targets.py
│   ├── test_pubmed_mining.py
│   └── test_clinvar_lookup.py
│
├── 03_orchestration/          # Test multi-service flows
│   ├── test_parallel_calls.py
│   └── test_fallback_logic.py
│
├── 04_frontend_integration/   # Test React components
│   ├── test_patient_context.py
│   └── test_result_rendering.py
│
└── outputs/                   # Captured real data
    ├── food_validator_vitamin_d.json
    ├── unified_care_ayesha.json
    ├── sae_features_curcumin.json
    └── evidence_omega3.json
```

---

## 🔧 PHASE 1: BACKEND SERVICE TESTS (30 mins)

### **Test 1.1: Food Validator API**
**File:** `01_backend_services/test_food_validator.py`

**What We Test:**
- Hit `/api/hypothesis/validate_food_ab_enhanced` with real inputs
- Capture full JSON response
- Verify structure (not content accuracy)
- Save to `outputs/food_validator_*.json`

**Test Cases:**
```python
TEST_CASES = [
    {
        "name": "vitamin_d_line3",
        "input": {
            "compound": "Vitamin D",
            "disease": "ovarian_cancer_hgs",
            "treatment_line": 3,
            "prior_therapies": ["carboplatin", "paclitaxel"],
            "use_llm": True
        },
        "verify": [
            "response has 'status' key",
            "response has 'sae_features' with 3 sub-keys",
            "response has 'provenance' with run_id",
            "response has 'targets' as list",
            "response has 'dosage' as dict"
        ]
    },
    {
        "name": "curcumin_line2",
        "input": {
            "compound": "Curcumin",
            "disease": "ovarian_cancer_hgs",
            "treatment_line": 2,
            "prior_therapies": ["carboplatin"],
            "use_llm": False  # Test without LLM
        },
        "verify": [
            "response returns within 5 seconds",
            "response has confidence value",
            "sae_features computed (not None)"
        ]
    },
    {
        "name": "green_tea_unknown_disease",
        "input": {
            "compound": "Green Tea Extract",
            "disease": "unknown_disease_xyz",  # Edge case
            "treatment_line": 1,
            "prior_therapies": [],
            "use_llm": False
        },
        "verify": [
            "response status is 'ERROR' or 'UNKNOWN'",
            "response has error message",
            "response still has provenance"
        ]
    }
]
```

**Output Captured:**
```json
// outputs/food_validator_vitamin_d.json
{
  "test_name": "vitamin_d_line3",
  "timestamp": "2024-12-03T15:30:00Z",
  "request_duration_ms": 2341,
  "input": { /* actual input sent */ },
  "output": { /* actual response received */ },
  "structure_verification": {
    "has_status": true,
    "has_sae_features": true,
    "has_provenance": true,
    "sae_keys": ["line_fitness", "cross_resistance", "sequencing_fitness"]
  }
}
```

---

### **Test 1.2: Unified Care API**
**File:** `01_backend_services/test_unified_care.py`

**Test Cases:**
```python
TEST_CASES = [
    {
        "name": "ayesha_hrd_positive",
        "input": {
            "patient_context": {
                "disease": "ovarian_cancer_hgs",
                "treatment_history": [
                    {"line": 1, "drugs": ["Carboplatin", "Paclitaxel"], "outcome": "partial_response"},
                    {"line": 2, "drugs": ["Olaparib"], "outcome": "progression"}
                ],
                "biomarkers": {
                    "brca1_mutant": False,
                    "hrd_positive": True,
                    "tp53_mutant": True
                }
            }
        },
        "verify": [
            "response has 'drug_recommendations' array",
            "response has 'food_recommendations' array",
            "response has 'integrated_confidence' number",
            "each drug has 'sae_features'",
            "provenance tracks both drug and food calls"
        ]
    },
    {
        "name": "brca_germline_positive",
        "input": {
            "patient_context": {
                "disease": "ovarian_cancer_hgs",
                "treatment_history": [],
                "biomarkers": {
                    "brca1_mutant": True,
                    "brca2_mutant": False
                }
            }
        },
        "verify": [
            "drug recommendations include PARP inhibitors",
            "confidence reflects germline status",
            "cross_resistance is LOW (no prior therapies)"
        ]
    }
]
```

**Output Captured:**
```json
// outputs/unified_care_ayesha.json
{
  "test_name": "ayesha_hrd_positive",
  "timestamp": "2024-12-03T15:35:00Z",
  "request_duration_ms": 3456,
  "parallel_execution": true,  // Verify asyncio.gather worked
  "input": { /* actual input */ },
  "output": {
    "drug_recommendations": [ /* actual drugs */ ],
    "food_recommendations": [ /* actual foods */ ],
    "integrated_confidence": 0.XX  // Real computed value
  },
  "structure_verification": {
    "drug_count": 2,
    "food_count": 3,
    "all_drugs_have_sae": true,
    "provenance_has_both_sources": true
  }
}
```

---

### **Test 1.3: SAE Features Computation**
**File:** `01_backend_services/test_sae_features.py`

**What We Test:**
- Call SAE service directly (not through validator)
- Verify null safety
- Capture raw SAE scores

**Test Cases:**
```python
TEST_CASES = [
    {
        "name": "omega3_line3_platinum_resistant",
        "input": {
            "compound": "Omega-3 Fatty Acids",
            "disease_context": "ovarian_cancer_hgs",
            "treatment_history": {
                "current_line": 3,
                "prior_therapies": ["carboplatin", "paclitaxel", "bevacizumab"]
            }
        },
        "verify": [
            "returns line_fitness score (0-1)",
            "returns cross_resistance score (0-1)",
            "returns sequencing_fitness score (0-1)",
            "scores are reasonable (not all 0 or 1)",
            "reasoning strings are present"
        ]
    },
    {
        "name": "vitamin_d_line1_treatment_naive",
        "input": {
            "compound": "Vitamin D",
            "disease_context": "ovarian_cancer_hgs",
            "treatment_history": {
                "current_line": 1,
                "prior_therapies": []
            }
        },
        "verify": [
            "cross_resistance is LOW (no priors)",
            "line_fitness >= 0.7 (appropriate for Line 1)",
            "sequencing_fitness is high"
        ]
    }
]
```

---

## 🔍 PHASE 2: DATA EXTRACTION TESTS (20 mins)

### **Test 2.1: ChEMBL Target Extraction**
**File:** `02_data_extraction/test_chembl_targets.py`

**Test Cases:**
```python
TEST_CASES = [
    {
        "name": "vitamin_d_targets",
        "input": "Vitamin D",
        "verify": [
            "returns VDR as target",
            "returns at least 2 targets",
            "targets have gene symbols",
            "response time < 3 seconds"
        ]
    },
    {
        "name": "curcumin_targets",
        "input": "Curcumin",
        "verify": [
            "returns NF-κB pathway targets",
            "returns COX-2",
            "handles alternative names (turmeric)"
        ]
    },
    {
        "name": "unknown_compound",
        "input": "FakeCompoundXYZ123",
        "verify": [
            "returns empty list or None",
            "does NOT crash",
            "logs warning"
        ]
    }
]
```

**Output:**
```json
// outputs/chembl_vitamin_d.json
{
  "compound": "Vitamin D",
  "chembl_id": "CHEMBL123456",
  "targets_found": ["VDR", "RXR", "CYP24A1"],
  "response_time_ms": 1234,
  "api_status": "success"
}
```

---

### **Test 2.2: PubMed Mining**
**File:** `02_data_extraction/test_pubmed_mining.py`

**Test Cases:**
```python
TEST_CASES = [
    {
        "name": "vitamin_d_ovarian_cancer",
        "input": {
            "compound": "Vitamin D",
            "disease": "ovarian cancer"
        },
        "verify": [
            "returns at least 5 papers",
            "papers have PMIDs",
            "papers have abstracts",
            "papers are relevant (keyword match)",
            "response includes synthesis if LLM enabled"
        ]
    },
    {
        "name": "obscure_compound",
        "input": {
            "compound": "Rare Compound XYZ",
            "disease": "ovarian cancer"
        },
        "verify": [
            "returns 0 papers or very few",
            "does NOT crash",
            "returns 'insufficient evidence' message"
        ]
    }
]
```

**Output:**
```json
// outputs/pubmed_vitamin_d_ovarian.json
{
  "query": "Vitamin D AND ovarian cancer",
  "papers_found": 8,
  "papers": [
    {
      "pmid": "26543123",
      "title": "Vitamin D and ovarian cancer survival",
      "abstract": "...",
      "relevance_score": 0.89
    }
  ],
  "llm_synthesis": "Vitamin D shows moderate evidence..."
}
```

---

## 🔄 PHASE 3: ORCHESTRATION TESTS (15 mins)

### **Test 3.1: Parallel API Calls**
**File:** `03_orchestration/test_parallel_calls.py`

**What We Test:**
- Verify `asyncio.gather()` is actually used
- Measure time savings vs sequential
- Verify both responses captured correctly

**Test:**
```python
async def test_parallel_execution():
    """
    Verify Unified Care calls drug + food APIs in parallel
    """
    start = time.time()
    
    # Make unified care request
    response = await call_unified_care(patient_context)
    
    elapsed = time.time() - start
    
    # Extract timestamps from provenance
    drug_timestamp = response['provenance']['drug_analysis']['timestamp']
    food_timestamp = response['provenance']['food_analysis']['timestamp']
    
    # Calculate time difference (should be < 500ms if parallel)
    time_diff = abs(
        datetime.fromisoformat(drug_timestamp) - 
        datetime.fromisoformat(food_timestamp)
    ).total_seconds()
    
    return {
        "total_time_seconds": elapsed,
        "timestamp_difference_seconds": time_diff,
        "is_parallel": time_diff < 0.5,  # ✅ Confirms parallel execution
        "speedup_vs_sequential": "~2x" if time_diff < 0.5 else "sequential"
    }
```

**Output:**
```json
{
  "total_request_time": 3.2,
  "drug_api_duration": 2.8,
  "food_api_duration": 2.5,
  "timestamp_difference": 0.15,  // ✅ < 0.5s = parallel confirmed
  "parallel_execution_verified": true,
  "speedup": "1.9x faster than sequential"
}
```

---

### **Test 3.2: Fallback Logic**
**File:** `03_orchestration/test_fallback_logic.py`

**Test Cases:**
```python
TEST_CASES = [
    {
        "name": "drug_api_fails",
        "mock": "force drug API to return 500",
        "verify": [
            "unified care still returns response",
            "food recommendations present",
            "error logged in provenance",
            "integrated confidence uses fallback (0.5)"
        ]
    },
    {
        "name": "food_api_fails",
        "mock": "force food API to return 500",
        "verify": [
            "unified care still returns response",
            "drug recommendations present",
            "uses FALLBACK_FOOD_TARGETS",
            "fallback foods are vitamin_d, curcumin, omega3"
        ]
    }
]
```

---

## 🎨 PHASE 4: FRONTEND INTEGRATION (15 mins)

### **Test 4.1: Patient Context Updates**
**File:** `04_frontend_integration/test_patient_context.py`

**Test:**
- Use Playwright/Selenium to interact with real frontend
- Click BRCA1/2 checkbox
- Capture state change
- Verify both biomarkers toggle

**Not Python - use browser automation:**
```javascript
// test_patient_context.spec.js
test('BRCA1/2 checkbox toggles both biomarkers', async ({ page }) => {
  await page.goto('http://localhost:3000/food-validator');
  
  // Get initial state
  const initialState = await page.evaluate(() => ({
    brca1: window.__PATIENT_CONTEXT__.biomarkers.brca1_mutant,
    brca2: window.__PATIENT_CONTEXT__.biomarkers.brca2_mutant
  }));
  
  // Click checkbox
  await page.click('[data-testid="brca-checkbox"]');
  
  // Get updated state
  const updatedState = await page.evaluate(() => ({
    brca1: window.__PATIENT_CONTEXT__.biomarkers.brca1_mutant,
    brca2: window.__PATIENT_CONTEXT__.biomarkers.brca2_mutant
  }));
  
  // Verify both changed
  console.log('Before:', initialState);
  console.log('After:', updatedState);
  // Save to outputs/brca_toggle_test.json
});
```

---

## 📋 TEST EXECUTION SCRIPT

### **Master Test Runner:**
**File:** `run_all_real_tests.sh`

```bash
#!/bin/bash

echo "⚔️ REAL TEST EXECUTION - NO MOCKS"
echo "=================================="

# 1. Start backend
echo "Starting backend..."
cd oncology-coPilot/oncology-backend-minimal
source venv/bin/activate
uvicorn api.main:app --host 0.0.0.0 --port 8000 &
BACKEND_PID=$!
sleep 5

# 2. Run Phase 1: Backend Services
echo "Phase 1: Backend Services..."
python .cursor/ayesha/hypothesis_validator/real_tests/01_backend_services/test_food_validator.py
python .cursor/ayesha/hypothesis_validator/real_tests/01_backend_services/test_unified_care.py
python .cursor/ayesha/hypothesis_validator/real_tests/01_backend_services/test_sae_features.py

# 3. Run Phase 2: Data Extraction
echo "Phase 2: Data Extraction..."
python .cursor/ayesha/hypothesis_validator/real_tests/02_data_extraction/test_chembl_targets.py
python .cursor/ayesha/hypothesis_validator/real_tests/02_data_extraction/test_pubmed_mining.py

# 4. Run Phase 3: Orchestration
echo "Phase 3: Orchestration..."
python .cursor/ayesha/hypothesis_validator/real_tests/03_orchestration/test_parallel_calls.py
python .cursor/ayesha/hypothesis_validator/real_tests/03_orchestration/test_fallback_logic.py

# 5. Collect outputs
echo "Collecting outputs..."
ls -lh .cursor/ayesha/hypothesis_validator/real_tests/outputs/

# 6. Generate summary report
python .cursor/ayesha/hypothesis_validator/real_tests/generate_summary.py

# 7. Cleanup
kill $BACKEND_PID

echo "✅ All tests complete. See outputs/ directory."
```

---

## 📊 OUTPUT FORMAT

### **Each test produces:**
```json
{
  "test_name": "vitamin_d_line3",
  "test_file": "01_backend_services/test_food_validator.py",
  "timestamp": "2024-12-03T15:30:00Z",
  "duration_ms": 2341,
  "status": "success",
  "input": { /* actual input sent */ },
  "output": { /* actual response received */ },
  "verification": {
    "structure_checks": {
      "has_status": true,
      "has_sae_features": true,
      "sae_feature_count": 3
    },
    "data_quality": {
      "targets_extracted": 3,
      "papers_found": 8,
      "confidence_in_range": true
    }
  },
  "issues": []  // Empty if no problems
}
```

---

## 🎯 SUMMARY REPORT

**After all tests, generate:**
`TEST_SUMMARY_REPORT.md`

```markdown
# Real Test Execution Summary

## Overview
- Tests Run: 15
- Duration: 45 minutes
- Backend Uptime: 100%
- API Failures: 0

## Phase 1: Backend Services (8 tests)
✅ Food Validator API (3 cases)
✅ Unified Care API (2 cases)
✅ SAE Features (2 cases)
✅ Evidence Mining (1 case)

## Phase 2: Data Extraction (4 tests)
✅ ChEMBL Targets (2 cases)
✅ PubMed Mining (2 cases)

## Phase 3: Orchestration (3 tests)
✅ Parallel Execution (verified 1.9x speedup)
✅ Fallback Logic (2 cases)

## Key Findings:
- Average API response time: 2.3s
- Parallel execution confirmed (drug + food APIs < 0.5s apart)
- All SAE computations null-safe
- PubMed mining successful (8-12 papers per query)
- Fallback logic works (no crashes when APIs fail)

## Output Files: 15 JSON files in outputs/
```

---

## ⏱️ ESTIMATED TIME

**Total:** 1.5-2 hours

- **Setup:** 10 mins (start backend, create directories)
- **Phase 1:** 30 mins (backend service tests)
- **Phase 2:** 20 mins (data extraction tests)
- **Phase 3:** 15 mins (orchestration tests)
- **Phase 4:** 15 mins (frontend integration)
- **Reporting:** 10 mins (generate summary)

---

## ✅ ACCEPTANCE CRITERIA

**We consider tests successful if:**
1. ✅ All APIs return valid JSON (structure verified)
2. ✅ No crashes or 500 errors on valid inputs
3. ✅ SAE features computed correctly
4. ✅ Parallel execution confirmed (<0.5s timestamp diff)
5. ✅ Fallback logic works (graceful degradation)
6. ✅ All outputs saved to files for Commander review

**NOT testing:**
- ❌ Clinical accuracy of recommendations
- ❌ Scientific validity of mechanisms
- ❌ Performance under load (that's separate)

---

## 🎯 WHAT COMMANDER GETS

**After execution:**
1. **15 JSON files** with real captured input/output
2. **Summary report** with key findings
3. **Verification** that system works end-to-end
4. **Proof** of parallel execution
5. **Evidence** of graceful degradation

**No mocks. No hardcoded values. Real data only.** ⚔️

---

**Ready to execute, Commander?** 

**Your orders:** 
1. **"Fire in the hole"** - Start execution now
2. **"Review first"** - Want to adjust plan before running
3. **"Skip for now"** - Move to Co-Pilot integration instead

— Zo






