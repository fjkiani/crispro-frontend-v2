# 📖 EXTRACTION PIECE 3.4: Circuit Breaker and Error Handling
**Source**: Lines 31000-31420 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents the circuit breaker mechanism implemented in the cohort extraction script, how it prevents credit burn, error handling patterns, and the actual circuit breaker trigger event with analysis.

---

## 🔍 KEY FINDINGS

### **Circuit Breaker Mechanism**

**Purpose**: Prevent excessive credit burn when error rates are high during cohort extraction.

**Implementation Location**: `scripts/sae/extract_sae_features_cohort.py`

**Configuration:**
```python
MAX_ERROR_RATE = 0.30  # Stop if >30% of calls fail
MIN_CALLS_BEFORE_CHECK = 20  # Wait for at least 20 calls before checking error rate
```

**Logic:**
```python
# Circuit breaker: check error rate
if total_variants_processed + total_variants_failed >= MIN_CALLS_BEFORE_CHECK:
    error_rate = total_variants_failed / (total_variants_processed + total_variants_failed)
    if error_rate > MAX_ERROR_RATE:
        logger.error(f"🚨 Circuit breaker triggered! Error rate: {error_rate:.1%} (>{MAX_ERROR_RATE:.1%})")
        logger.error(f"   Variants processed: {total_variants_processed}, failed: {total_variants_failed}")
        logger.error("   Stopping extraction to prevent credit burn. Check configuration and logs.")
        break
```

**Key Features:**
- **Threshold**: 30% error rate triggers circuit breaker
- **Minimum Calls**: Requires at least 20 calls before checking (prevents false positives)
- **Automatic Stop**: Breaks extraction loop when threshold exceeded
- **Checkpoint Save**: Saves progress before stopping

---

### **Actual Circuit Breaker Event**

**Trigger Details:**
```
🚨 Circuit breaker triggered! Error rate: 100.0% (>30.0%)
   Variants processed: 0, failed: 50
   Stopping extraction to prevent credit burn. Check configuration and logs.
```

**Context:**
- **Patient**: TCGA-13-0889
- **Error Pattern**: All 50 variants failed with "Reference allele mismatch"
- **Error Type**: HTTP 400 from SAE service
- **Root Cause**: Variant data doesn't match Ensembl reference sequence

**Error Messages (Sample):**
```
⚠️  HTTP 400 for 6:70049259 ->T: {"detail":"SAE service error: {\"detail\":\"Reference allele mismatch\"}"}
⚠️  HTTP 400 for 3:187449580 ->C: {"detail":"SAE service error: {\"detail\":\"Reference allele mismatch\"}"}
⚠️  HTTP 400 for 13:32972497 ->C: {"detail":"SAE service error: {\"detail\":\"Reference allele mismatch\"}"}
```

**Result:**
- Circuit breaker correctly stopped extraction
- 30 patients successfully processed before trigger
- 1 patient failed (TCGA-13-0889)
- Checkpoint saved with failed patient marked
- Extraction can resume from checkpoint

---

### **Error Handling Patterns**

#### **1. Retry Logic**

**Configuration:**
```python
MAX_RETRIES = 3
RETRY_DELAY = 5.0
REQUEST_TIMEOUT = 180.0  # 3 minutes per request
```

**Implementation:**
- Retries failed requests up to 3 times
- 5-second delay between retries
- Handles timeout exceptions separately
- Logs warnings for each retry attempt

#### **2. Failed Patient Handling**

**Pattern:**
- Failed patients added to `failed_patients` list
- Checkpoint tracks failed patients
- Extraction continues with next patient
- Failed patients skipped on resume

**Checkpoint Structure:**
```python
checkpoint = {
    "completed_patients": [...],
    "failed_patients": [...],
    "last_processed_index": ...
}
```

#### **3. Variant-Level Error Handling**

**Per-Variant Errors:**
- HTTP 400 (bad request): Logged as warning, retried
- HTTP 403 (forbidden): Feature flags not enabled, stops
- HTTP 503 (service unavailable): Service not configured, stops
- Timeout: Retried with exponential backoff
- Other exceptions: Logged and skipped

---

### **Cost Control Mechanisms**

**Multiple Layers:**

1. **Feature Flags:**
   - `ENABLE_SAE_COHORT_RUN`: Explicit flag required to run
   - `ENABLE_EVO2_SAE`: Must be enabled
   - `ENABLE_TRUE_SAE`: Must be enabled

2. **Hard Limits:**
   - `MAX_PATIENTS`: Default 50, configurable via env var
   - `MAX_TOTAL_VARIANTS`: Default 2500, configurable via env var
   - `MAX_VARIANTS_PER_PATIENT`: Hard limit of 50 variants per patient

3. **Circuit Breaker:**
   - Stops if error rate >30% after 20+ calls
   - Prevents runaway costs from systematic failures

4. **Checkpointing:**
   - Saves progress every 10 patients
   - Allows resume from last checkpoint
   - Prevents re-processing completed patients

---

### **Circuit Breaker Analysis**

**What Happened:**
- Patient TCGA-13-0889 had 50 variants
- All 50 variants failed with "Reference allele mismatch"
- Error rate: 100% (50 failed / 50 attempted)
- Circuit breaker triggered correctly
- Extraction stopped to prevent further credit burn

**Why It Worked:**
- Circuit breaker correctly identified systematic failure
- Stopped before processing more patients with same issue
- Saved checkpoint with 30 successful patients
- Failed patient marked for investigation

**What Was Saved:**
- ✅ 30 patients successfully extracted
- ✅ 1,000+ variants with SAE features
- ✅ Checkpoint file for resume
- ✅ Failed patient identified for investigation

**Next Steps:**
- Investigate "Reference allele mismatch" errors
- Fix variant data or reference sequence handling
- Resume extraction from checkpoint (skips failed patient)

---

## 📊 KEY INSIGHTS

### **Circuit Breaker Design**

1. **Threshold Selection**: 30% error rate is reasonable - allows some failures but stops systematic issues
2. **Minimum Calls**: 20 calls prevents false positives from small sample sizes
3. **Automatic Stop**: Prevents manual intervention during runaway failures
4. **Checkpoint Save**: Ensures progress isn't lost when circuit breaker triggers

### **Error Patterns**

1. **Reference Allele Mismatch**: Common error when variant data doesn't match Ensembl reference
2. **Systematic Failures**: One patient with all variants failing triggers circuit breaker
3. **Retry Logic**: Handles transient errors but not systematic data issues
4. **Graceful Degradation**: Continues with successful patients, marks failures

### **Cost Protection**

1. **Multiple Safeguards**: Feature flags + limits + circuit breaker
2. **Early Detection**: Circuit breaker triggers before excessive costs
3. **Progress Preservation**: Checkpointing ensures work isn't lost
4. **Resume Capability**: Can restart from checkpoint after fixing issues

---

## 🔗 CONTEXT & CONNECTIONS

- **Builds on**: Bug fixes (Piece 3.3), Real data extraction (Piece 3.2)
- **Prevents**: Excessive credit burn from systematic failures
- **Enables**: Safe cohort extraction with automatic failure detection
- **Key Insight**: Circuit breaker is a protective mechanism, not a bug

---

## 📝 NOTES

- Circuit breaker worked as designed
- "Reference allele mismatch" is a data quality issue, not a code bug
- Failed patient can be skipped on resume
- 30 patients successfully extracted before trigger
- Checkpoint mechanism allows safe resume

---

## 🎯 QUESTIONS RESOLVED

- ✅ What triggers circuit breaker? → Error rate >30% after 20+ calls
- ✅ Why did it trigger? → Patient TCGA-13-0889 had 100% failure rate (50/50 variants failed)
- ✅ What was saved? → 30 patients successfully extracted, checkpoint saved
- ✅ Is this a bug? → No, circuit breaker worked correctly to prevent credit burn
- ✅ How to proceed? → Investigate reference allele mismatch, resume from checkpoint

