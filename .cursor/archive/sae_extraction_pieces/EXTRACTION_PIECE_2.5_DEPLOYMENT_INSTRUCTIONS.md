# 📖 EXTRACTION PIECE 2.5: Deployment Instructions
**Source**: Lines 10370-10800 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents the complete deployment instructions for Evo2 and SAE services on Modal, including CLI setup, deployment commands, environment variable configuration, and health check procedures.

---

## 🔍 KEY FINDINGS

### **1. Modal CLI & Auth Setup**

**One-time setup per machine:**

```bash
# From repo root
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

# Install Modal (if not installed)
pip3 install modal

# Login (will open browser)
modal token new
```

**Key Point**: Authentication is required before deployment.

---

### **2. Deploy Evo2 Service (evo-service)**

**Deployment Command:**

```bash
cd src/services/evo_service
modal deploy main.py
```

**After Deployment:**
- Modal prints URL like: `https://<something>--evo-service-fastapi-app.modal.run`
- Set environment variables:

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

export EVO_URL_1B="https://<your-evo-url>.modal.run"
export ENABLE_EVO2_SAE=1
```

**Note**: Add these to `.env` file if using one for backend.

**Key Point**: Evo2 service must be deployed first (SAE depends on it).

---

### **3. Deploy SAE Service (sae-service)**

**Deployment Command:**

```bash
cd src/services/sae_service
modal deploy main.py
```

**After Deployment:**
- Modal prints URL like: `https://<something>--sae-service-fastapi-app.modal.run`
- Set environment variables:

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

export SAE_SERVICE_URL="https://<your-sae-url>.modal.run"
export ENABLE_TRUE_SAE=1
```

**Key Point**: SAE service requires Evo2 service to be deployed first.

---

### **4. Health Checks**

**Evo2 Activations Health Check:**

```bash
# Start backend
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload

# In another terminal, test activations endpoint:
curl -X POST http://127.0.0.1:8000/api/evo/score_variant_with_activations \
  -H 'Content-Type: application/json' \
  -d '{
    "assembly": "GRCh38",
    "chrom": "17",
    "pos": 43044295,
    "ref": "T",
    "alt": "G",
    "window": 8192,
    "model_id": "evo2_1b"
  }' | jq 'keys'
```

**Expected Response**: Should see keys like `delta_score`, `layer_26_alt_activations`, `provenance`.

**SAE Health Check:**

```bash
# Direct Modal endpoint (if implemented)
curl "$SAE_SERVICE_URL/health" || echo "Health endpoint not implemented yet"

# Via backend RUO endpoint:
curl -X POST http://127.0.0.1:8000/api/sae/extract_features \
  -H 'Content-Type: application/json' \
  -d '{"activations": [[[0.0]*4096]]}' | jq 'keys'
```

**Key Point**: Health checks verify services are operational before extraction.

---

### **5. Post-Deployment: Extract Real SAE Features**

**Once Both Services Deployed:**

```bash
# Run cohort extraction script
cd scripts/sae
python extract_sae_features_cohort.py \
  --input data/tcga_ov_platinum_cohort.json \
  --output data/sae_features_tcga_ov_platinum.json \
  --checkpoint data/sae_extraction_checkpoint.json
```

**Key Point**: Real extraction can begin once both services are live.

---

## 📊 KEY INSIGHTS

### **Deployment Sequence**

1. **Evo2 First**: SAE depends on Evo2 activations
2. **SAE Second**: Requires Evo2 URL to be set
3. **Health Checks**: Verify both services before extraction
4. **Environment Variables**: Must be set correctly for backend to use services

### **Dependencies**

- **Evo2 Service**: Provides layer 26 activations
- **SAE Service**: Decodes activations → 32K features
- **Backend**: Orchestrates calls to both services
- **Environment Variables**: Control which services are enabled

### **Resource Requirements**

- **Evo2 Service**: H100 GPU (for model inference)
- **SAE Service**: H100 GPU (for feature decoding)
- **Timeline**: 2-4 hours for deployment (as stated in blocker)

---

## 🔗 CONTEXT & CONNECTIONS

- **Enables**: Real SAE feature extraction (Sprint 1 completion)
- **Depends on**: Modal account setup and authentication
- **Required for**: Biomarker discovery with real data
- **Key Insight**: Deployment is the blocker preventing real extraction

---

## 📝 NOTES

- Deployment instructions are step-by-step and executable
- Health checks verify services before use
- Environment variables must be set correctly
- Post-deployment extraction can begin immediately
- All commands are bash-compatible

---

## 🎯 QUESTIONS RESOLVED

- ✅ How to deploy Evo2? → `cd src/services/evo_service && modal deploy main.py`
- ✅ How to deploy SAE? → `cd src/services/sae_service && modal deploy main.py`
- ✅ What environment variables? → `EVO_URL_1B`, `ENABLE_EVO2_SAE`, `SAE_SERVICE_URL`, `ENABLE_TRUE_SAE`
- ✅ How to verify? → Health check endpoints and test API calls
- ✅ What's the sequence? → Evo2 first, then SAE, then health checks, then extraction

