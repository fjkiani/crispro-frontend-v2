# Next Steps After Backend Restart

**Date**: January 20, 2025  
**Status**: ✅ Code Complete, ⏸️ Awaiting Backend Restart

---

## ✅ COMPLETED

1. **Code Changes**: evo2_7b migration complete (3 files, 9 changes)
2. **Environment Variables**: Added to `.env` file:
   - `ENABLE_TRUE_SAE=1`
   - `ENABLE_EVO2_SAE=1`
   - `SAE_SERVICE_URL=https://crispro--sae-service-saeservice-api.modal.run`
   - `EVO_URL_7B=https://crispro--evo-service-evoservice7b-api-7b.modal.run`

---

## ⏸️ REQUIRED: Backend Restart

**The backend needs to be restarted to pick up the new `.env` variables.**

**Location**: `.env` file in `oncology-coPilot/oncology-backend-minimal/`

**After restart, test with**:
```bash
curl -X POST http://localhost:8000/api/sae/extract_features \
  -H "Content-Type: application/json" \
  -d '{
    "chrom": "17",
    "pos": 43044295,
    "ref": "T",
    "alt": "G",
    "assembly": "GRCh38",
    "window": 8192,
    "model_id": "evo2_7b"
  }' | jq '.provenance.model'
```

**Expected**: Should show "trained weights" (not "random init")

---

## 🧪 TESTING PLAN (After Backend Restart)

### Step 1: Single Variant Test
```bash
python3 -c "
import asyncio
import httpx

async def test():
    async with httpx.AsyncClient(timeout=180.0) as client:
        response = await client.post(
            'http://localhost:8000/api/sae/extract_features',
            json={
                'chrom': '17',
                'pos': 43044295,
                'ref': 'T',
                'alt': 'G',
                'assembly': 'GRCh38',
                'window': 8192,
                'model_id': 'evo2_7b'
            }
        )
        if response.status_code == 200:
            result = response.json()
            prov = result.get('provenance', {})
            print(f'Model: {prov.get(\"model\")}')
            print(f'd_in: {prov.get(\"d_in\")}')
            if 'trained weights' in prov.get('model', '').lower():
                print('✅ TRAINED WEIGHTS LOADED!')
            else:
                print('⚠️  Still using random weights')
        else:
            print(f'Error: {response.text}')

asyncio.run(test())
"
```

### Step 2: Small Batch Test (1 patient)
```bash
export ENABLE_SAE_COHORT_RUN=1
export MAX_PATIENTS=1
export MAX_TOTAL_VARIANTS=10

python3 scripts/sae/extract_sae_features_cohort.py
```

**Verify**:
- Provenance shows "trained weights"
- `d_in` is 4096 (not 1920)
- Features extracted successfully

### Step 3: Full Cohort Re-Extraction (After Step 2 Passes)
```bash
export MAX_PATIENTS=66  # Or remove limit
export MAX_TOTAL_VARIANTS=3000

python3 scripts/sae/extract_sae_features_cohort.py
```

---

## 📊 EXPECTED RESULTS

**Before (Random Weights)**:
- Model: "Goodfire/Evo-2-Layer-26-Mixed (random init for RUO)"
- d_in: 1920
- Features: May be noise

**After (Trained Weights)**:
- Model: "Goodfire/Evo-2-Layer-26-Mixed (trained weights)"
- d_in: 4096
- Features: Biologically meaningful

---

**Status**: ✅ **READY** - Backend restart required, then proceed with testing

