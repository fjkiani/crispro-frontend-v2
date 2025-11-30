# ⚔️ SAE PHASE 3 - TESTING & VALIDATION GUIDE ⚔️

**Date**: January 13, 2025  
**Status**: READY TO TEST  
**Owner**: Zo (Lead Commander)

---

## 🎯 **WHAT WE'RE TESTING**

**SAE Phase 1+2 services integrated into Complete Care v2 orchestrator:**

### **Phase 1 Services** (Pre-NGS + Post-NGS)
1. ✅ Next-Test Recommender (HRD → ctDNA → SLFN11 → ABCB1)
2. ✅ Hint Tiles (Max 4, suggestive tone)
3. ✅ Mechanism Map (Pre-NGS: gray, Post-NGS: color-coded)

### **Phase 2 Services** (Post-NGS only)
4. ✅ SAE Feature Computation (Manager's C1-C10)
5. ✅ Resistance Detection (2-of-3 triggers, HR restoration)

---

## 📁 **TEST FILES CREATED**

### **Test Payloads** (`.cursor/ayesha/test_payloads/`)
1. `01_pre_ngs.json` - Ayesha TODAY (no NGS data)
2. `02_brca1_biallelic.json` - BRCA1 biallelic loss (HRD=58, high DDR)
3. `03_her2_positive.json` - HER2 amplification (NCT06819007 eligible)

### **Test Data Reference**
- `TEST_NGS_DATA_AYESHA.json` - Full test cases with expected results

### **Test Script**
- `test_sae_phase3_integration.sh` - Automated test runner

---

## 🚀 **HOW TO RUN TESTS**

### **Step 1: Start Backend Server**

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal
uvicorn main:app --reload --port 8000
```

**Expected Output:**
```
INFO:     Uvicorn running on http://127.0.0.1:8000
INFO:     Application startup complete.
```

---

### **Step 2: Health Check**

```bash
curl http://localhost:8000/api/ayesha/complete_care_v2/health | jq
```

**Expected Response:**
```json
{
  "status": "operational",
  "service": "complete_care_v2",
  "for_patient": "Ayesha Kiani (Stage IVB ovarian cancer)",
  "sae_phase1_enabled": true,
  "sae_phase2_enabled": true,
  "sae_policy": "MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md (Jan 13, 2025)",
  "capabilities": [
    "next_test_recommender_priority_hrd_ctdna_slfn11_abcb1",
    "hint_tiles_max4_suggestive_tone",
    "mechanism_map_7chips_her2_integrated",
    "sae_feature_computation_manager_c1_c10",
    "resistance_detection_2_of_3_triggers"
  ]
}
```

**Validation:**
- ✅ `sae_phase1_enabled: true`
- ✅ `sae_phase2_enabled: true`
- ✅ Both phases in capabilities list

---

### **Step 3: Run Automated Tests**

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
./.cursor/ayesha/test_sae_phase3_integration.sh
```

**Expected Output:**
```
⚔️ SAE PHASE 3 INTEGRATION TESTS ⚔️
==================================

🔍 Checking backend health...
✅ Backend is operational
  - SAE Phase 1: true
  - SAE Phase 2: true

TEST: PRE-NGS
Description: Ayesha TODAY (no NGS data) - should recommend tests
✅ TEST PASSED
Results:
  - SAE Features: awaiting_ngs
  - Resistance Alert: awaiting_ngs
  - Next Tests: 3
  - Hint Tiles: 4
  - Mechanism Map: awaiting_ngs

TEST: BRCA1-BIALLELIC
Description: Ayesha with BRCA1 biallelic loss (HRD=58) - high DDR burden
✅ TEST PASSED
Results:
  - SAE Features: computed
  - Resistance Alert: computed
  - DNA Repair Capacity: 0.75-0.85
  ...

⚔️ ALL TESTS COMPLETE ⚔️
```

---

### **Step 4: Manual Testing (Deep Dive)**

#### **Test 1: Pre-NGS (Ayesha TODAY)**

```bash
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d @.cursor/ayesha/test_payloads/01_pre_ngs.json | jq
```

**What to Validate:**

```json
{
  "sae_features": {
    "status": "awaiting_ngs",
    "message": "SAE features will be computed once tumor NGS data is available..."
  },
  "resistance_alert": {
    "status": "awaiting_ngs",
    "message": "Resistance detection requires baseline NGS data"
  },
  "next_test_recommender": {
    "recommendations": [
      {
        "test": "HRD (MyChoice CDx)",
        "priority": "P0",
        "rationale": "Critical for PARP maintenance eligibility",
        ...
      },
      {
        "test": "ctDNA (Guardant360)",
        "priority": "P1",
        ...
      },
      {
        "test": "HER2 IHC",
        "priority": "P1",
        "rationale": "Unlocks HER2-targeted trials (NCT06819007)",
        ...
      }
    ],
    "total_tests": 3
  },
  "hint_tiles": {
    "hint_tiles": [
      {
        "type": "next_test",
        "title": "📋 Recommended Next Test",
        "recommendation": "Consider ordering HRD test...",
        ...
      },
      ...
    ],
    "total_tiles": 4
  },
  "mechanism_map": {
    "status": "awaiting_ngs",
    "message": "Mechanism map will be available once tumor NGS results are uploaded...",
    "chips": [
      {"pathway": "DDR", "value": "Awaiting NGS", "color": "gray"},
      {"pathway": "MAPK", "value": "Awaiting NGS", "color": "gray"},
      ...
    ]
  }
}
```

**✅ Validation Checklist:**
- [ ] `sae_features.status == "awaiting_ngs"`
- [ ] `resistance_alert.status == "awaiting_ngs"`
- [ ] `next_test_recommender.total_tests >= 2`
- [ ] `hint_tiles.total_tiles <= 4`
- [ ] `mechanism_map.status == "awaiting_ngs"`
- [ ] All mechanism map chips are gray

---

#### **Test 2: Post-NGS with BRCA1 Biallelic Loss**

```bash
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d @.cursor/ayesha/test_payloads/02_brca1_biallelic.json | jq
```

**What to Validate:**

```json
{
  "sae_features": {
    "dna_repair_capacity": 0.75-0.85,  // HIGH (Manager's C5 formula)
    "pathway_burden_ddr": 0.7-0.9,      // HIGH
    "pathway_burden_mapk": 0.1-0.3,     // LOW
    "pathway_burden_pi3k": 0.1-0.3,     // LOW
    "pathway_burden_vegf": 0.2-0.4,     // MODERATE
    "pathway_burden_her2": 0.0,         // UNKNOWN (no HER2 test)
    "io_eligible": false,               // TMB=8.5 < 20, MSI=MSS
    "cross_resistance_risk": 0.0,       // Treatment-naive
    "essentiality_hrr_genes": 0.7-0.9,  // BRCA1 in profile
    "mechanism_vector": [0.8, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0],  // 7D vector
    "provenance": {
      "manager_policy": "MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md (C1-C10)"
    }
  },
  "resistance_alert": {
    "resistance_detected": false,
    "trigger_count": 0,
    "triggers_met": [],
    "hr_restoration_suspected": false,
    "hrd_signal": {
      "triggered": false,
      "reason": "No baseline HRD available"
    },
    "dna_repair_signal": {
      "triggered": false,
      "reason": "No baseline DNA repair capacity available"
    },
    "recommended_actions": [],
    "provenance": {
      "manager_policy": "MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md (C7, R2)"
    }
  },
  "hint_tiles": {
    "hint_tiles": [
      {
        "type": "trial_matched",
        "title": "🎯 Trial Opportunity",
        "recommendation": "Consider PARP maintenance trials...",
        ...
      },
      ...
    ]
  },
  "mechanism_map": {
    "status": "computed",
    "chips": [
      {"pathway": "DDR", "value": "0.85", "color": "green"},   // HIGH
      {"pathway": "MAPK", "value": "0.20", "color": "gray"},   // LOW
      {"pathway": "PI3K", "value": "0.15", "color": "gray"},   // LOW
      {"pathway": "VEGF", "value": "0.35", "color": "gray"},   // MODERATE (below 0.40)
      {"pathway": "HER2", "value": "0.00", "color": "gray"},   // UNKNOWN
      {"pathway": "IO", "value": "0.00", "color": "gray"},     // NOT ELIGIBLE
      {"pathway": "Efflux", "value": "0.00", "color": "gray"}  // Treatment-naive
    ]
  }
}
```

**✅ Validation Checklist:**
- [ ] `sae_features.dna_repair_capacity > 0.70` (HIGH for BRCA1)
- [ ] `sae_features.pathway_burden_ddr > 0.70` (HIGH DDR burden)
- [ ] `sae_features.io_eligible == false` (TMB < 20, MSI=MSS)
- [ ] `sae_features.mechanism_vector` is 7D array
- [ ] `resistance_alert.resistance_detected == false` (baseline)
- [ ] `resistance_alert.trigger_count == 0`
- [ ] `mechanism_map.status == "computed"`
- [ ] DDR chip is green (value > 0.70)
- [ ] `provenance.manager_policy` references Manager's answers

---

#### **Test 3: Post-NGS with HER2 Amplification**

```bash
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d @.cursor/ayesha/test_payloads/03_her2_positive.json | jq
```

**What to Validate:**

```json
{
  "sae_features": {
    "dna_repair_capacity": 0.30-0.45,  // LOW (HRD=22, no BRCA)
    "pathway_burden_her2": 0.7-0.9,    // HIGH (ERBB2 amplification)
    "pathway_burden_pi3k": 0.4-0.6,    // MODERATE (PIK3CA E545K)
    "pathway_burden_ddr": 0.2-0.4,     // LOW (low HRD)
    "mechanism_vector": [0.3, 0.2, 0.5, 0.3, 0.8, 0.0, 0.0],  // HER2=0.8, PI3K=0.5
    ...
  },
  "hint_tiles": {
    "hint_tiles": [
      {
        "type": "trial_matched",
        "title": "🎯 HER2-Targeted Trial Available",
        "recommendation": "Consider NCT06819007 (HER2-ADC trial)...",
        ...
      },
      ...
    ]
  },
  "mechanism_map": {
    "chips": [
      {"pathway": "DDR", "value": "0.30", "color": "gray"},    // LOW
      {"pathway": "HER2", "value": "0.80", "color": "green"},  // HIGH ⚔️
      {"pathway": "PI3K", "value": "0.50", "color": "yellow"}, // MODERATE
      ...
    ]
  }
}
```

**✅ Validation Checklist:**
- [ ] `sae_features.pathway_burden_her2 > 0.70` (HIGH for HER2+)
- [ ] `sae_features.pathway_burden_pi3k > 0.40` (MODERATE for PIK3CA)
- [ ] `sae_features.dna_repair_capacity < 0.50` (LOW for low HRD)
- [ ] HER2 chip is green or yellow (value > 0.40)
- [ ] Hint tiles mention HER2-targeted trial
- [ ] PI3K chip is yellow (0.40-0.70)

---

## 📊 **WHAT EACH TEST PROVES**

### **Test 1: Pre-NGS (Ayesha TODAY)**
**Proves:**
- ✅ SAE Phase 1 works without NGS data
- ✅ Next-test recommender prioritizes correctly (HRD → ctDNA → HER2 IHC)
- ✅ Hint tiles generated (max 4)
- ✅ Mechanism map shows "awaiting NGS" gracefully
- ✅ SAE Phase 2 returns "awaiting_ngs" (doesn't crash)

### **Test 2: BRCA1 Biallelic (Post-NGS)**
**Proves:**
- ✅ SAE features computed (Manager's C1-C10 formulas)
- ✅ DNA repair capacity formula works (C5)
- ✅ High DDR burden detected (pathway_burden_ddr > 0.70)
- ✅ Essentiality for HRR genes computed (C3)
- ✅ Resistance detection baseline works (no false positives)
- ✅ Mechanism map color-coded correctly (DDR green)
- ✅ 7D mechanism vector includes HER2 dimension

### **Test 3: HER2+ (Post-NGS with NCT06819007 Gate)**
**Proves:**
- ✅ HER2 pathway detection works (7th dimension)
- ✅ Mechanism map shows HER2 chip correctly
- ✅ Hint tiles recommend HER2-targeted trial
- ✅ PI3K pathway detected (PIK3CA mutation)
- ✅ Low HRD correctly shows low DDR burden

---

## 🔬 **DATA SOURCES (WHERE NGS DATA CAME FROM)**

**Test data is based on real genomic profiles from:**

1. **cBioPortal TCGA Ovarian Cancer**
   - URL: https://www.cbioportal.org/study/summary?id=ov_tcga_pan_can_atlas_2018
   - Used for: BRCA1/2, TP53, HRD scores, TMB distributions

2. **AACR Project GENIE - Ovarian Cancer subset**
   - Used for: Somatic mutation VAF (variant allele frequency) ranges

3. **Published Case Reports**
   - HER2+ ovarian cancer: Ross JS et al. (2017) - <5% prevalence
   - POLE-mutated hypermutated ovarian: Howitt BE et al. (2015) - <2%
   - PARP resistance via HR restoration: Pettitt SJ et al. (2018)

**Disclaimer**: These are synthetic but clinically realistic profiles for validation. Real patient data would come from:
- Guardant360 (ctDNA)
- FoundationOne CDx (tissue NGS)
- MyChoice CDx (HRD)
- Tempus xT (comprehensive panel)

---

## ⚠️ **TROUBLESHOOTING**

### **Issue: Backend not running**
```bash
# Start backend
cd oncology-backend-minimal
uvicorn main:app --reload --port 8000
```

### **Issue: Tests fail with "connection refused"**
- Check backend is running on port 8000
- Try: `curl http://localhost:8000/api/ayesha/complete_care_v2/health`

### **Issue: "Service temporarily unavailable"**
- Check logs for import errors
- Verify all SAE services are in `api/services/`
- Check Python dependencies installed

### **Issue: SAE features not computed**
- Verify `tumor_context` is provided in payload
- Check `tumor_context.somatic_mutations` is non-empty array
- Verify HRD score is provided

---

## ⚔️ **SUCCESS CRITERIA**

**Phase 3 Integration is successful if:**

1. ✅ Health endpoint shows `sae_phase1_enabled: true` and `sae_phase2_enabled: true`
2. ✅ Pre-NGS request returns "awaiting_ngs" gracefully (no crashes)
3. ✅ Post-NGS request computes SAE features (DNA repair capacity, mechanism vector)
4. ✅ Resistance detection runs without errors (baseline = no resistance)
5. ✅ Mechanism map color-coded correctly (green for high, yellow for moderate, gray for low)
6. ✅ Hint tiles include relevant clinical actions
7. ✅ Provenance references Manager's policy document

**All 3 test cases should PASS without errors!** ⚔️

---

**COMMANDER - READY TO TEST WHEN BACKEND IS RUNNING!** 🔥

**Run this to start:**
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn main:app --reload --port 8000
```

**Then in another terminal:**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
./.cursor/ayesha/test_sae_phase3_integration.sh
```

