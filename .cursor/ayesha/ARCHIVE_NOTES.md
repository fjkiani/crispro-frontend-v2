# 📦 Archive Notes

**Date:** January 28, 2025  
**Purpose:** Document archived/superseded files

---

## Synthetic Lethality Documentation

### ✅ Active Document (Single Source of Truth)
- **`SYNTHETIC_LETHALITY_COMPLETE.md`** - Complete consolidated documentation
  - Contains: V1 features, V2 enhancements, benchmark results, all context

### 📦 Archived Documents (Superseded)

These documents are kept for historical reference but are superseded by `SYNTHETIC_LETHALITY_COMPLETE.md`:

1. **`SYNTHETIC_LETHALITY_FRONTEND_PLAN.md`** - V1 implementation plan
   - Status: All features implemented
   - Superseded by: Complete doc V1 section
   - **Action:** Archived (reference only)

2. **`SYNTHETIC_LETHALITY_V2_IMPROVEMENT_PLAN.md`** - V2 enhancement plan
   - Status: All HIGH and MEDIUM priority tasks implemented
   - Superseded by: Complete doc V2 section
   - **Action:** Archived (reference only)

3. **`SYNTHETIC_LETHALITY_V2_IMPLEMENTATION_SUMMARY.md`** - V2 implementation summary
   - Status: Historical record of V2 completion
   - Superseded by: Complete doc V2 section
   - **Action:** Archived (reference only)

4. **`SYNTHETIC_LETHALITY_V2_VERIFICATION.md`** - V2 verification checklist
   - Status: All items verified complete
   - Superseded by: Complete doc implementation status
   - **Action:** Archived (reference only)

5. **`SYNTHETIC_LETHALITY_BENCHMARK_PLAN.md`** - Original benchmark plan (V1)
   - Status: Superseded by V2 plan
   - Superseded by: `SYNTHETIC_LETHALITY_BENCHMARK_PLAN_V2.md`
   - **Action:** Archived (reference only)

6. **`SYNTHETIC_LETHALITY_BENCHMARK_PLAN_V2.md`** - Detailed benchmark implementation plan
   - Status: Superseded by actual implementation and results
   - **Note:** Contains detailed technical plan, but actual results documented in Complete doc
   - **Action:** Archived (reference only)

---

## Benchmark Implementation Files

### ✅ Active Files (In Use)
- `scripts/benchmark_sl/benchmark_efficacy.py` - CORRECT benchmark (uses Evo2)
- `scripts/benchmark_sl/test_cases_pilot.json` - 10 test cases
- `scripts/benchmark_sl/ISSUES_FOUND.md` - Documents critical issues
- `scripts/benchmark_sl/README.md` - Usage guide

### 📦 Reference Files
- `scripts/benchmark_sl/benchmark_synthetic_lethality.py` - Original (tests rules only)
- `scripts/benchmark_sl/create_pilot_dataset.py` - Dataset creation
- `scripts/benchmark_sl/download_depmap.py` - DepMap data processing

---

## Usage

- **For current status:** See `SYNTHETIC_LETHALITY_COMPLETE.md`
- **For benchmark issues:** See `scripts/benchmark_sl/ISSUES_FOUND.md`
- **For historical reference:** See archived documents above

---

## Summary of Changes

**What Was Fixed:**
1. ✅ Created corrected benchmark using `/api/efficacy/predict` (actually uses Evo2)
2. ✅ Identified GUIDANCE_FAST bypass issue
3. ✅ Documented that original benchmark tested rules, not ML
4. ✅ Validated Evo2 is working (100% usage rate)
5. ✅ Established realistic baseline (50% drug match)

**What Was Archived:**
- All redundant documentation files
- Original benchmark plan (superseded by actual implementation)
- V1/V2 implementation plans (consolidated into Complete doc)

---

**Last Updated:** January 28, 2025
