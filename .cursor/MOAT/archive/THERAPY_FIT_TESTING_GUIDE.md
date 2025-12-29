# 🧪 Therapy Fit Testing Guide

**Date:** January 2025  
**Status:** ✅ **READY FOR TESTING**  
**Purpose:** Guide for testing and generating results from THERAPY_FIT_TEST_CASES.md

---

## 🎯 Quick Start

### Prerequisites
1. **Backend server running:**
   ```bash
   cd oncology-coPilot/oncology-backend-minimal
   uvicorn api.main:app --reload
   ```

2. **Test script ready:**
   - `scripts/generate_therapy_fit_results.py` - Main test runner

---

## 🚀 Running Tests

### Option 1: Run Single Test Case (Testing Mode)
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/generate_therapy_fit_results.py --test-case AYESHA-001
```

### Option 2: Run All Test Cases (Testing Mode)
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/generate_therapy_fit_results.py --all
```

### Option 3: Run Demo (Presentation Mode)
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/demo_therapy_fit.py --test-case AYESHA-001
```

### Option 4: Run All Demos (Interactive)
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/demo_therapy_fit.py --all --interactive
```

**Difference:**
- `generate_therapy_fit_results.py` - Testing/validation mode (detailed comparison, JSON output)
- `demo_therapy_fit.py` - Demo/presentation mode (formatted output, visual display)

---

## 📊 Available Test Cases

| Patient ID | Disease | Mutations | Expected Top Drug |
|------------|---------|-----------|-------------------|
| **AYESHA-001** | Ovarian (HGSOC) | MBD4+TP53 | Olaparib |
| **MM-001** | Multiple Myeloma | KRAS G12D | Trametinib |
| **MEL-001** | Melanoma | BRAF V600E | Vemurafenib |
| **OV-002** | Ovarian | BRCA1 truncation | Olaparib |
| **MM-002** | Multiple Myeloma | 5 MAPK variants | Trametinib |

---

## 📋 What Gets Tested

For each test case, the script:
1. ✅ Calls `/api/efficacy/predict` with test case mutations
2. ✅ Extracts top 3 drugs with scores
3. ✅ Extracts pathway alignment scores
4. ✅ Compares actual vs. expected results:
   - Top drug match
   - Confidence score range
   - Evidence tier
   - Pathway alignment scores
5. ✅ Generates JSON output with full results

---

## 📄 Output Format

Results are saved to `therapy_fit_generated_results.json`:

```json
{
  "timestamp": "2025-01-XX...",
  "api_base_url": "http://127.0.0.1:8000",
  "test_cases_run": 5,
  "successful": 5,
  "failed": 0,
  "results": [
    {
      "test_case": "Ayesha (MBD4+TP53 HGSOC)",
      "patient_id": "AYESHA-001",
      "disease": "ovarian_cancer",
      "top_drugs": [
        {
          "rank": 1,
          "name": "Olaparib",
          "efficacy_score": 0.82,
          "confidence": 0.85,
          "evidence_tier": "supported",
          "badges": ["DDR", "HRD", "PARP"],
          "has_insights": true
        }
      ],
      "pathway_scores": {
        "DDR": 0.78,
        "RAS/MAPK": 0.12
      },
      "comparison": {
        "top_drug_match": true,
        "confidence_in_range": true,
        "evidence_tier_match": true,
        "pathway_alignment": {
          "DDR": {
            "passed": true,
            "actual": 0.78,
            "expected": "[0.75, 0.85]"
          }
        }
      }
    }
  ]
}
```

---

## 🔍 Interpreting Results

### ✅ Pass Criteria
- **Top Drug Match:** Actual top drug contains expected drug name
- **Confidence Range:** Actual confidence within expected range
- **Evidence Tier:** Actual tier matches expected tier
- **Pathway Alignment:** Actual pathway score within expected range

### ⚠️ Common Issues
1. **API not running:** Make sure backend server is running on port 8000
2. **Timeout:** Evo2 scoring can take 2-5 minutes per request
3. **Drug name mismatch:** Drug names may vary (e.g., "Olaparib" vs "olaparib")
4. **Pathway name mismatch:** Pathway names may differ (e.g., "DDR" vs "ddr")

---

## 📝 Next Steps After Testing

1. **Review Results:**
   - Check `therapy_fit_generated_results.json`
   - Compare actual vs. expected
   - Note any discrepancies

2. **Update Documentation:**
   - Update `THERAPY_FIT_TEST_CASES.md` with actual results if needed
   - Document any discrepancies found

3. **Iterate:**
   - Fix any issues found
   - Re-run tests
   - Validate fixes

---

## 🛠️ Troubleshooting

### Issue: "Connection refused"
**Solution:** Make sure backend server is running:
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload
```

### Issue: "Timeout"
**Solution:** Evo2 scoring takes time. Increase timeout in script or wait longer.

### Issue: "Test case not found"
**Solution:** Check patient ID spelling. Available IDs:
- AYESHA-001
- MM-001
- MEL-001
- OV-002
- MM-002

---

## 📚 Related Files

- **Test Cases:** `.cursor/MOAT/THERAPY_FIT_TEST_CASES.md`
- **Test Script:** `scripts/generate_therapy_fit_results.py` (testing/validation)
- **Demo Script:** `scripts/demo_therapy_fit.py` (presentation/demo)
- **Results:** `scripts/therapy_fit_generated_results.json` (generated)

---

*Created: January 2025*  
*Status: ✅ **READY FOR USE***

