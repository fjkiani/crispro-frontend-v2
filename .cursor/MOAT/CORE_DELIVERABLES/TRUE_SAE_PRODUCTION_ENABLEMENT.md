# TRUE SAE Production Enablement Guide

**Date:** January 28, 2025  
**Status:** ✅ **READY FOR ENABLEMENT**  
**Purpose:** Step-by-step guide to enable TRUE SAE in production  
**Location:** `.cursor/MOAT/CORE_DELIVERABLES/TRUE_SAE_PRODUCTION_ENABLEMENT.md`

---

## 🎯 Quick Start

**Enable TRUE SAE in 3 steps:**

1. **Set Environment Variable:**
   ```bash
   export ENABLE_TRUE_SAE_PATHWAYS=true
   ```

2. **Restart Backend Server:**
   ```bash
   # Restart your API server to load the new environment variable
   ```

3. **Test with MBD4+TP53 Case:**
   - Navigate to Clinical Dossier
   - Verify "TRUE SAE" badge appears
   - Verify DDR_bin gauge shows ~0.88

---

## 📋 Detailed Enablement Steps

### **Step 1: Configure Environment Variable**

**Option 1: .env File (Recommended)**
```bash
# Add to .env file in oncology-backend-minimal directory
ENABLE_TRUE_SAE_PATHWAYS=true
```

**Option 2: Environment Variable**
```bash
export ENABLE_TRUE_SAE_PATHWAYS=true
```

**Option 3: Docker/Container**
```yaml
# Add to docker-compose.yml
environment:
  - ENABLE_TRUE_SAE_PATHWAYS=true
```

**Option 4: System Environment**
```bash
# Add to ~/.bashrc or ~/.zshrc
export ENABLE_TRUE_SAE_PATHWAYS=true
```

### **Step 2: Verify Configuration**

**Check Flag is Loaded:**
```python
# In Python shell or test script
from api.config import get_feature_flags
flags = get_feature_flags()
print(f"TRUE SAE Pathways Enabled: {flags.get('enable_true_sae_pathways', False)}")
```

**Expected Output:**
```
TRUE SAE Pathways Enabled: True
```

**Check Configuration File:**
```bash
# Verify flag is set in config.py
grep -n "ENABLE_TRUE_SAE_PATHWAYS" api/config.py
```

**Expected Output:**
```
57:ENABLE_TRUE_SAE_PATHWAYS = os.getenv("ENABLE_TRUE_SAE_PATHWAYS", "false").lower() in ("true", "1", "yes")
```

### **Step 3: Restart Backend Server**

**Restart API Server:**
```bash
# Stop current server (Ctrl+C or kill process)
# Then restart
cd oncology-coPilot/oncology-backend-minimal
python -m uvicorn api.main:app --reload
```

**Verify Server Started:**
```bash
# Check logs for flag status
# Should see: "DEBUG: ENABLE_TRUE_SAE resolved to: True" (if debug logging enabled)
```

### **Step 4: Test with MBD4+TP53 Case**

**Test Patient Profile:**
```json
{
  "mutations": [
    {"gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2", "type": "germline"},
    {"gene": "TP53", "hgvs_p": "p.Arg175His", "type": "somatic"}
  ],
  "disease": "ovarian_cancer_hgsoc",
  "stage": "IVB",
  "tmb": 8.5,
  "msi_status": "MSS"
}
```

**Expected API Response:**
```json
{
  "sae_features": {
    "provenance": {
      "sae": "true_sae",
      "sae_version": "v1",
      "mapping_version": "v1",
      "sae_diagnostics": {
        "ddr_bin_score": 0.88,
        "ddr_sae_score": 0.88
      }
    }
  },
  "trials": [
    {
      "nct_id": "NCT04284969",
      "sae_source": "true_sae",
      "ddr_bin_score": 0.88,
      "mechanism_fit_score": 0.989,
      "combined_score": 0.892
    }
  ]
}
```

**Expected Frontend Display:**
- ✅ "TRUE SAE" badge in trial matching section
- ✅ DDR_bin gauge showing ~0.88 in pathway disruption section
- ✅ Mechanism alignment showing DDR_bin scores for DDR pathway
- ✅ Tooltips explaining TRUE SAE vs PROXY SAE

---

## 🔍 Verification Checklist

### **Backend Verification:**
- [ ] `ENABLE_TRUE_SAE_PATHWAYS=true` set in environment
- [ ] Backend server restarted
- [ ] Flag loaded correctly (`get_feature_flags()` returns `True`)
- [ ] API response includes `sae_features.provenance.sae = "true_sae"`
- [ ] API response includes `sae_features.provenance.sae_diagnostics.ddr_bin_score`
- [ ] Trial responses include `sae_source = "true_sae"`
- [ ] Trial responses include `ddr_bin_score ≈ 0.88` (for MBD4+TP53)

### **Frontend Verification:**
- [ ] "TRUE SAE" badge displays in trial cards
- [ ] DDR_bin gauge displays in pathway disruption section
- [ ] Mechanism alignment shows DDR_bin scores for DDR pathway
- [ ] Tooltips explain TRUE SAE vs PROXY SAE difference
- [ ] All components render without errors

---

## 🐛 Troubleshooting

### **Issue: Flag Not Loading**

**Symptom:** `get_feature_flags()` returns `False` even after setting environment variable

**Solution:**
1. Verify environment variable is set: `echo $ENABLE_TRUE_SAE_PATHWAYS`
2. Restart backend server (environment variables loaded at startup)
3. Check `.env` file is in correct location (same directory as `api/config.py`)
4. Verify `.env` file syntax: `ENABLE_TRUE_SAE_PATHWAYS=true` (no quotes, no spaces)

### **Issue: TRUE SAE Not Used**

**Symptom:** `sae_features.provenance.sae = "proxy"` even with flag enabled

**Solution:**
1. Verify TRUE SAE features are provided to `compute_sae_features()`:
   - Check `sae_features` parameter is not `None`
   - Verify Evo2 activations are available
   - Check `ENABLE_TRUE_SAE` flag is also enabled (if using TRUE SAE feature extraction)
2. Verify Feature→Pathway Mapping exists:
   - Check `api/resources/sae_feature_mapping.true_sae_diamonds.v1.json` exists
   - Verify mapping file is valid JSON
3. Check logs for errors in `_compute_sae_diagnostics()`

### **Issue: Frontend Not Displaying TRUE SAE**

**Symptom:** No "TRUE SAE" badge or DDR_bin gauge visible

**Solution:**
1. Verify API response includes `sae_source` and `ddr_bin_score`:
   - Check browser network tab for API response
   - Verify `trial.sae_source = "true_sae"`
   - Verify `trial.ddr_bin_score` is present
2. Verify frontend components are updated:
   - Check `SAESourceIndicator.jsx` exists
   - Check `DDRBinGauge.jsx` exists
   - Verify components are imported correctly
3. Check browser console for errors

---

## 📊 Expected Behavior

### **When TRUE SAE is Enabled:**
- ✅ `sae_features.provenance.sae = "true_sae"`
- ✅ `sae_features.provenance.sae_diagnostics.ddr_bin_score` computed
- ✅ Trial responses include `sae_source = "true_sae"`
- ✅ Trial responses include `ddr_bin_score` (0.0-1.0)
- ✅ Frontend displays "TRUE SAE" badge
- ✅ Frontend displays DDR_bin gauge

### **When TRUE SAE is Disabled (Default):**
- ✅ `sae_features.provenance.sae = "proxy"`
- ✅ `sae_features.provenance.sae_diagnostics.ddr_bin_score = None` (or not computed)
- ✅ Trial responses include `sae_source = "proxy"` (or not included)
- ✅ Frontend displays "PROXY SAE" badge (or no badge)
- ✅ Frontend does not display DDR_bin gauge

---

## 🔗 Related Documents

- **Implementation Details:** [DELIVERABLES_1.5_AND_2_TRUE_SAE_INTEGRATION.md](DELIVERABLES_1.5_AND_2_TRUE_SAE_INTEGRATION.md)
- **Impact Analysis:** [TRUE_SAE_TRIAL_MATCHING_IMPACT.md](TRUE_SAE_TRIAL_MATCHING_IMPACT.md)
- **Test Report:** [DELIVERABLE_1_5_AND_2_TEST_REPORT.md](../SAE_INTELLIGENCE/DELIVERABLE_1_5_AND_2_TEST_REPORT.md)
- **Configuration:** `api/config.py` (line 57)

---

## ✅ Success Criteria

**Enablement Complete When:**
1. ✅ `ENABLE_TRUE_SAE_PATHWAYS=true` set in environment
2. ✅ Backend server restarted
3. ✅ API response includes `sae_source = "true_sae"` for MBD4+TP53 case
4. ✅ API response includes `ddr_bin_score ≈ 0.88` for MBD4+TP53 case
5. ✅ Frontend displays "TRUE SAE" badge
6. ✅ Frontend displays DDR_bin gauge showing ~0.88

---

*Guide Created: January 28, 2025*  
*Status: ✅ READY FOR ENABLEMENT*


