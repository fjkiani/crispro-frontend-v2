# SAE Redeploy and Re-extraction Plan

**Date**: January 20, 2025  
**Status**: Ready to Execute  
**Reason**: Backend routing fix ensures 7B model is used correctly

---

## 🔧 Changes Made

1. **Backend Routing Fixed**: `api/routers/evo.py` now respects `model_id="evo2_7b"` instead of forcing 1B
2. **Centralized Default**: `api/config.py` now has `DEFAULT_EVO_MODEL="evo2_7b"` as single source of truth
3. **SAE Service**: Already configured to use `evo2_7b_base` internally (line 155 of `src/services/sae_service/main.py`)

---

## 📋 Step 1: Redeploy SAE Service

**Purpose**: Ensure SAE service is using latest code and correctly configured for 7B

```bash
# Navigate to SAE service directory
cd src/services/sae_service

# Deploy to Modal
modal deploy main.py

# Note the deployed URL (should be similar to existing)
# Expected: https://crispro--sae-service-saeservice-api.modal.run
```

**Verification**:
```bash
# Health check
curl https://crispro--sae-service-saeservice-api.modal.run/health

# Should return: {"status": "healthy", ...}
```

---

## 📋 Step 2: Verify Backend Configuration

**Ensure environment variables are set** (in `.env` or exported):

```bash
# From project root
cd oncology-coPilot/oncology-backend-minimal

# Check .env file or export:
export SAE_SERVICE_URL="https://crispro--sae-service-saeservice-api.modal.run"
export EVO_URL_7B="https://crispro--evo-service-evoservice7b-api-7b.modal.run"
export ENABLE_EVO2_SAE=1
export ENABLE_TRUE_SAE=1
export DEFAULT_EVO_MODEL="evo2_7b"  # Optional, already defaults to 7B in config.py
```

**Restart backend** (if running):
```bash
# Kill existing backend
lsof -ti:8000 | xargs -r kill -9

# Start backend with new config
cd oncology-coPilot/oncology-backend-minimal
source venv/bin/activate  # If using venv
PYTHONPATH=. uvicorn api.main:app --host 127.0.0.1 --port 8000 --reload
```

---

## 📋 Step 3: Test with Small Batch

**Before full re-extraction, test with 1-2 patients**:

```bash
# From project root
cd scripts/sae

# Run small batch test (1 patient)
python extract_sae_features_cohort.py \
  --max-patients 1 \
  --max-total-variants 50 \
  --model-id evo2_7b \
  --clear-checkpoint

# Verify output shows:
# - model_id: "evo2_7b" in provenance
# - d_in: 4096 (7B dimension)
# - trained weights loaded
```

**Expected Output**:
```json
{
  "provenance": {
    "model_id": "evo2_7b",
    "d_in": 4096,
    "model": "Goodfire/Evo-2-Layer-26-Mixed (trained weights)",
    ...
  }
}
```

---

## 📋 Step 4: Full Cohort Re-extraction

**After small batch test succeeds**:

```bash
# Clear previous checkpoint (if exists)
rm -f data/validation/sae_cohort/extraction_checkpoint.json

# Run full extraction
python extract_sae_features_cohort.py \
  --max-patients 66 \
  --max-total-variants 3000 \
  --model-id evo2_7b

# Monitor progress
# Expected time: ~10-20 minutes
# Expected cost: Monitor Modal dashboard
```

**Output File**:
- `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

**Verification**:
```bash
# Check output file
python -c "
import json
with open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json') as f:
    data = json.load(f)
    print(f'Patients: {len(data[\"patients\"])}')
    print(f'Model ID: {data[\"provenance\"][\"model_id\"]}')
    print(f'Dimensions: d_in={data[\"provenance\"][\"d_in\"]}, d_hidden={data[\"provenance\"][\"d_hidden\"]}')
    print(f'Weights: {data[\"provenance\"][\"model\"]}')
"
```

---

## ✅ Success Criteria

1. **SAE Service**: Deployed and health check passes
2. **Backend**: Routes to 7B when `model_id="evo2_7b"` is specified
3. **Small Batch Test**: Completes successfully with 7B model
4. **Full Extraction**: All patients processed with 7B model and trained weights
5. **Output Provenance**: Shows `model_id: "evo2_7b"`, `d_in: 4096`, trained weights loaded

---

## ⚠️ Cost Control

- **Small batch test first** (1 patient, ~50 variants)
- **Monitor Modal dashboard** during extraction
- **Circuit breakers active** (max 50 errors per 100 requests)
- **Checkpoint system** allows resuming if interrupted

---

## 🔍 Troubleshooting

**If SAE service health check fails**:
- Wait 30-60 seconds (cold start)
- Check Modal dashboard for service status
- Verify Modal token is valid: `modal token list`

**If backend still routes to 1B**:
- Verify `DEFAULT_EVO_MODEL` in `api/config.py` is `"evo2_7b"`
- Check `FORCE_MODEL_ID` is empty (not set in env)
- Restart backend after config changes

**If extraction fails**:
- Check backend logs: `tail -f logs/backend.log` (if logging enabled)
- Verify `ENABLE_EVO2_SAE=1` and `ENABLE_TRUE_SAE=1`
- Check Modal service URLs are correct



