# SAE EXTRACTION AUDIT: TECHNICAL RETROSPECTIVE

**Status**: ✅ COMPREHENSIVE A-Z EXTRACTION
**Objective**: Complete audit of SAE extraction process (TCGA-OV).

## 1. EXECUTIVE SUMMARY

**Extraction Target**: SAE features from Evo2 layer 26 activations (TCGA-OV).
**Results**:
- ✅ 149 patients extracted (TCGA-OV)
- ✅ ~2,897 variants processed
- ✅ ~33 hours runtime (2 min/variant)
- ✅ >95% success rate

**Critical Lessons**:
1. **Mock data first** - Verify pipeline before expensive extraction.
2. **Circuit breaker** - Stop at 30% error rate.
3. **Assembly mismatch** - GRCh37 vs GRCh38 causes "Reference allele mismatch".
4. **Feature index bug** - Aggregate BEFORE top-k, not after.

---

## 2. EXTRACTION SCRIPTS

### Primary: `extract_sae_features_cohort.py`
**Location**: `scripts/sae/extract_sae_features_cohort.py`
**Features**:
- **Checkpoint/Resume**: Saves every 10 patients.
- **Circuit Breaker**: Stops at 30% error rate.
- **Cost Control**: `MAX_PATIENTS`, `MAX_TOTAL_VARIANTS` caps.
- **Assembly Handling**: Hardcoded `GRCh37` for TCGA data.

### Alternative (Robust): `extract_true_sae_cohort_from_cbioportal.py`
**Location**: `oncology-coPilot/oncology-backend-minimal/sae_validation/scripts/extract_true_sae_cohort_from_cbioportal.py`
**Features**:
- **Preflight Validation**: Auto-detects assembly (GRCh37 vs 38).
- **Caching**: SHA256 hash of variant params.
- **Budget Guards**: `--budget_seconds` hard stop.
- **Atomic Writes**: `.tmp` -> rename to prevent corruption.

---

## 3. CRITICAL BUGS FIXED

### BUG 1: Feature Index Bug (Critical)
**Problem**: Flattened 3D tensor `[1, 8193, 32768]` -> `[268M]` before top-k. Resulted in invalid indices > 32767.
**Fix**: Aggregate `mean(dim=1)` -> `[32768]` -> then `topk`.
**Impact**: All previous analysis invalidated.

### BUG 2: Modal Payload Size
**Problem**: Returning full 32K float vector (1-2GB JSON).
**Fix**: Return only `top_features` (Top-64 sparse).

### BUG 3: Assembly Mismatch
**Problem**: TCGA is GRCh37, Evo2 defaults GRCh38. Caused 100% failure for some patients.
**Fix**: Explicit assembly specification and preflight detection.

---

## 4. ERROR PATTERNS

1. **"Reference allele mismatch"**: Wrong assembly (GRCh37 vs 38).
2. **Ensembl 400 errors**: Invalid chromosome positions.
3. **HTTP 403**: Feature flags not set (`ENABLE_EVO2_SAE=1`).

---

## 5. PIPELINE & BEST PRACTICES

### Pre-Extraction
- Verify feature flags.
- Run preflight check (Script 2).
- Set budget limits.

### During Extraction
- Use `--resume`.
- Monitor `extract_errors.jsonl`.
- Concurrency: 4-8 threads.

### Post-Extraction
- Validate JSON format.
- Run unit tests on extracted features.
- Clean temp files.

---

## 6. FILE STRUCTURE

**Output**: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`
**Format**:
```json
{
  "cohort": "TCGA-OV",
  "patients": [
    {
      "patient_id": "TCGA-XX-XXXX",
      "variants": [
        {
          "top_features": [{"index": 1234, "value": 0.85}, ...]
        }
      ],
      "aggregated_features": {
        "mean_features": [...]
      }
    }
  ]
}
```
