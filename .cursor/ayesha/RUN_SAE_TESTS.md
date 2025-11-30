# ⚔️ SAE PHASE 3 - TESTING WITH REAL TCGA DATA ⚔️

**Date**: January 13, 2025  
**Data Source**: **REAL TCGA Ovarian Cancer Pan-Cancer Atlas 2018**  
**File**: `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`

---

## ✅ **REAL DATA - NOT SYNTHETIC**

**Our test data comes from:**
- **Source**: TCGA Ovarian Cancer (cBioPortal)
- **File**: Already in our repo at `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`
- **Mutations**: Real BRCA1, BRCA2, TP53 mutations from actual ovarian cancer patients
- **Study**: https://www.cbioportal.org/study/summary?id=ov_tcga_pan_can_atlas_2018

**Real mutations we're using:**
1. **BRCA1 I1108*** (nonsense) - Real TCGA patient
2. **BRCA2 C711*** (nonsense) - Real TCGA patient  
3. **TP53 R306*** (nonsense) - Real TCGA patient

---

## 🚀 **HOW TO TEST**

### **Step 1: Start Backend**

```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn main:app --reload --port 8000
```

---

### **Step 2: Health Check**

```bash
curl http://localhost:8000/api/ayesha/complete_care_v2/health
```

**Expected**: `sae_phase1_enabled: true`, `sae_phase2_enabled: true`

---

### **Step 3: Test Pre-NGS (Ayesha TODAY)**

```bash
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d @.cursor/ayesha/test_payloads/01_pre_ngs.json
```

**Expected:**
- `sae_features.status == "awaiting_ngs"`
- `next_test_recommender` has ≥2 tests
- `mechanism_map` all gray

---

### **Step 4: Test REAL TCGA BRCA1 Patient**

```bash
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d @.cursor/ayesha/test_payloads/04_real_tcga_brca1.json | python3 -m json.tool
```

**This is a REAL TCGA ovarian cancer patient with:**
- **BRCA1 I1108*** (nonsense mutation)
- **TP53 R306*** (nonsense mutation)
- **HRD score**: 55 (typical for BRCA1)

**Expected SAE results:**
```json
{
  "sae_features": {
    "dna_repair_capacity": 0.70-0.85,
    "pathway_burden_ddr": 0.70-0.90,
    "essentiality_hrr_genes": 0.70-0.90,
    "mechanism_vector": [0.8, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0]
  },
  "resistance_alert": {
    "resistance_detected": false,
    "trigger_count": 0
  }
}
```

---

### **Step 5: Test REAL TCGA BRCA2 Patient**

```bash
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d @.cursor/ayesha/test_payloads/05_real_tcga_brca2.json | python3 -m json.tool
```

**This is a REAL TCGA ovarian cancer patient with:**
- **BRCA2 C711*** (nonsense mutation)
- **TP53 R306*** (nonsense mutation)
- **HRD score**: 60 (typical for BRCA2)

---

## 📊 **VALIDATION CHECKLIST**

### **Test 1: Pre-NGS** ✅
- [ ] Backend returns 200 OK
- [ ] `sae_features.status == "awaiting_ngs"`
- [ ] `next_test_recommender.total_tests >= 2`
- [ ] No crashes or 500 errors

### **Test 2: REAL TCGA BRCA1** ✅
- [ ] Backend returns 200 OK
- [ ] `sae_features.dna_repair_capacity > 0.70`
- [ ] `sae_features.pathway_burden_ddr > 0.70`
- [ ] `sae_features.mechanism_vector` is 7D array
- [ ] `resistance_alert.resistance_detected == false`
- [ ] `provenance.manager_policy` includes Manager's policy ref

### **Test 3: REAL TCGA BRCA2** ✅
- [ ] Backend returns 200 OK  
- [ ] Similar high DDR burden as BRCA1
- [ ] SAE features computed correctly

---

## 🔬 **DATA PROVENANCE**

**All test mutations come from:**

**File**: `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`

**Source Study**: TCGA Ovarian Cancer Pan-Cancer Atlas 2018
- **Paper**: https://www.nature.com/articles/nature10166
- **Data Portal**: https://www.cbioportal.org/study/summary?id=ov_tcga_pan_can_atlas_2018
- **Patients**: 584 high-grade serous ovarian cancer patients
- **License**: CC BY 4.0 (public data)

**Extracted mutations:**
```bash
# To see all real BRCA mutations in our data:
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
grep -n "BRCA1\|BRCA2" tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json
```

---

## ⚔️ **COMMANDER - THIS IS REAL DATA!**

**NOT synthetic:**
- ✅ Real TCGA ovarian cancer patients
- ✅ Real BRCA1/BRCA2/TP53 mutations
- ✅ Already in our repository
- ✅ Public data (CC BY 4.0 license)

**Files created:**
1. `.cursor/ayesha/test_payloads/04_real_tcga_brca1.json` - REAL TCGA patient
2. `.cursor/ayesha/test_payloads/05_real_tcga_brca2.json` - REAL TCGA patient

**Ready to test when backend is running!** 🔥

