# 📖 EXTRACTION PIECE 2.1: Phase 1 SAE Implementation Complete
**Source**: Lines 2644-3050 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents the completion of Phase 1 SAE implementation, including what was built, guardrails, status, and next steps. It also includes the learning from Evo2 notebook review and the high-level plan for SAE integration.

---

## 🔍 KEY FINDINGS

### **What Was Learned (New Since Last Iteration)**

From Evo2 notebook review (lines 2644-2689):

1. **Evo2 Notebook Pattern**
   - Official Evo2 notebook shows how to extract activations
   - Pattern: `ObservableEvo2` wrapper + `BatchTopKTiedSAE`
   - SAE weights from Hugging Face: `Goodfire/Evo-2-Layer-26-Mixed`

2. **Current Codebase Status**
   - Evo2 service exists but doesn't expose activations endpoint
   - No SAE service exists yet
   - Proxy SAE features computed in `sae_feature_service.py`

3. **Gap Analysis**
   - Need: Evo2 activations endpoint
   - Need: SAE Modal service (following notebook pattern)
   - Need: SAE router in backend
   - Need: Integration into `sae_feature_service.py` (optional, diagnostics)

---

### **How This Will Be Applied (Near-Term, Non-Breaking)**

**Phase 1 Approach:**
- Build infrastructure WITHOUT changing production behavior
- Add new endpoints behind feature flags
- Compute SAE features as diagnostics only
- Full provenance tracking
- No changes to `/api/efficacy/predict` or WIWFM scoring

**Guardrails:**
- Feature flags: `ENABLE_EVO2_SAE`, `ENABLE_TRUE_SAE` (default: disabled)
- Diagnostics only: No scoring changes
- Proxy SAE remains default
- All true SAE features are optional diagnostics

---

### **High-Level Plan**

#### **A. Guardrails / Scope**
- No changes to `/api/efficacy/predict` or WIWFM scoring
- Proxy SAE remains default
- True SAE gated behind TWO flags
- All true SAE features are diagnostics-only (Phase 1)
- Full provenance tracking

#### **B. Phase 1 – Minimal End-to-End SAE Extraction (Foundation)**
1. **Evo2 Activations Endpoint** (Modal)
   - Add `POST /score_variant_with_activations` to `evo_service`
   - Returns layer 26 activations [batch, seq_len, 4096]
   - Gated by `ENABLE_EVO2_SAE` flag

2. **SAE Modal Service**
   - File: `src/services/sae_service/main.py`
   - `ObservableEvo2` + `BatchTopKTiedSAE` (from notebook)
   - Downloads HF weights: `Goodfire/Evo-2-Layer-26-Mixed`
   - Endpoint: `POST /extract_features`
   - Input: variant (chr, pos, ref, alt) OR activations tensor
   - Output: top-k SAE features (indices + values)

3. **SAE Router** (Backend)
   - File: `api/routers/sae.py`
   - Endpoint: `POST /api/sae/extract_features`
   - Gated by `ENABLE_TRUE_SAE` flag
   - Calls Modal SAE service

4. **SAE Client Service**
   - File: `api/services/sae_model_service.py`
   - Methods: `extract_features_from_activations()`, `extract_features_from_variant()`, `health_check()`

5. **SAE Feature Service Integration** (Optional)
   - Extend `compute_sae_features(sae_features: Optional[Dict] = None)`
   - If `sae_features` provided (from true SAE), compute diagnostics
   - Provenance: `sae = "proxy" | "proxy+true"`

#### **C. Phase 2 – Cohort SAE Extraction & Correlation (Later)**
- Extract SAE features for labeled patient cohorts
- Run correlation analysis: which of 32K features correlate with outcomes?
- Build `scripts/extract_sae_features_cohort.py`
- Build `api/services/biomarker_correlation_service.py`

#### **D. Phase 3 – Carefully Using SAE in Product (After Validation)**
- Only after validation shows clear benefit
- Add SAE-based lifts to confidence scores
- Integrate into S/P/E pipeline

---

### **What Was Built (Phase 1 Complete)**

1. **Evo2 Activations Endpoint** ✅
   - Modal: `POST /score_variant_with_activations` in `src/services/evo_service/main.py`
   - Router: `POST /api/evo/score_variant_with_activations` (gated by `ENABLE_EVO2_SAE`)

2. **SAE Modal Service** ✅
   - File: `src/services/sae_service/main.py`
   - `ObservableEvo2` + `BatchTopKTiedSAE` (from official Evo2 notebook)
   - Downloads HF weights: `Goodfire/Evo-2-Layer-26-Mixed`
   - Endpoint: `POST /extract_features`

3. **SAE Router** ✅
   - File: `oncology-coPilot/oncology-backend-minimal/api/routers/sae.py`
   - Endpoint: `POST /api/sae/extract_features` (gated by `ENABLE_TRUE_SAE`)
   - Registered in `main.py`

4. **SAE Client Service** ✅
   - File: `oncology-coPilot/oncology-backend-minimal/api/services/sae_model_service.py`
   - Methods: `extract_features_from_activations()`, `extract_features_from_variant()`, `health_check()`

5. **SAE Feature Service Integration** ✅
   - Extended `compute_sae_features(sae_features: Optional[Dict] = None)`
   - Computes diagnostics: `ddr_sae_score`, `io_sae_score`, `mapk_sae_score`
   - **Diagnostics only** — no scoring changes
   - Provenance: `sae = "proxy" | "proxy+true"`

6. **Feature Flags** ✅
   - `ENABLE_EVO2_SAE` — Evo2 activations endpoint (default: disabled)
   - `ENABLE_TRUE_SAE` — True SAE features (default: disabled)

7. **Provenance Tracking** ✅
   - All SAE paths include explicit provenance
   - Clear markers: `sae_ready`, `sae_version`, `method`, `source`

---

### **Guardrails**

✅ No changes to `/api/efficacy/predict` or WIWFM scoring  
✅ Proxy SAE remains default  
✅ True SAE gated behind TWO flags  
✅ All true SAE features are diagnostics-only (Phase 1)  
✅ Full provenance tracking  

---

### **Status**

**Phase 1: Complete and ready for deployment**

**Next:** Phase 2 requires manager approval before starting cohort analysis and correlation studies.

---

## 📊 KEY INSIGHTS

### **Implementation Strategy**

1. **Non-Breaking First**: Phase 1 adds infrastructure without changing production behavior
2. **Feature Flags**: Everything gated behind flags for safety
3. **Diagnostics Only**: True SAE features don't affect scoring in Phase 1
4. **Provenance**: Full tracking of SAE source and method
5. **Modular**: Each component can be tested independently

### **Architecture Decisions**

1. **Modal Services**: Evo2 and SAE services deployed on Modal cloud
2. **Backend Router**: Thin wrapper around Modal services
3. **Client Service**: Abstraction layer for calling SAE
4. **Integration Point**: Optional integration into existing `sae_feature_service.py`

---

## 🔗 CONTEXT & CONNECTIONS

- **Builds on**: All previous pieces (understanding SAE, manager policy, Evo2 relationship)
- **Implements**: The infrastructure needed for real SAE features
- **Leads to**: Phase 2 (cohort extraction) and Phase 3 (integration)
- **Key Insight**: Phase 1 is foundation-building without production risk

---

## 📝 NOTES

- Phase 1 is intentionally conservative (diagnostics only)
- All infrastructure is built and ready
- Next step requires manager approval for Phase 2
- The implementation follows the official Evo2 notebook pattern
- Feature flags provide safety and control

---

## 🎯 QUESTIONS RESOLVED

- ✅ How do we extract real SAE features? → Modal services + Evo2 activations endpoint
- ✅ What was built in Phase 1? → Complete infrastructure (Evo2 activations + SAE service)
- ✅ Are production systems affected? → No, everything is gated and diagnostics-only
- ✅ What's next? → Phase 2 requires manager approval

