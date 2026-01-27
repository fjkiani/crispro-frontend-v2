# 🔬 TECHNICAL AUDIT: AlphaFold (Boltz-2), Evo2, & Scientific Engine
**Date**: January 2025  
**Auditor**: Zo (Senior Engineer)  
**Status**: COMPREHENSIVE DEEP DIVE ✅

---

## 🎯 EXECUTIVE SUMMARY

**Verdict**: **STRONG FOUNDATION, CRITICAL GAPS IDENTIFIED** 🚨

### **Current State:**
- ✅ **Evo2 Service**: Production-ready, 8 endpoints, 1B/7B/40B models, solid integration
- ⚠️ **Boltz-2 Service**: Functional but **NOT AlphaFold 3** (different model), fast mode trade-offs
- 🔴 **Scientific Engine**: Fragmented orchestration, missing critical integrations, hardcoded URLs

### **Critical Findings:**
1. **Boltz-2 ≠ AlphaFold 3**: Publication claims AF3 but uses Boltz-2 (different architecture)
2. **Hardcoded Service URLs**: `http://127.0.0.1:8000` in production code (design.py line 96)
3. **No Real Enformer**: Chromatin signals are deterministic stubs (confirmed in audit)
4. **Missing Integration Tests**: No end-to-end validation of Evo2 → Boltz → Assassin Score pipeline
5. **Fast Mode Trade-offs**: Boltz-2 uses `msa='empty'` (pLDDT 50-70) vs full MSA (60+ min)

### **What We Can Push To:**
- **Evo2**: Layer 26 activations (SAE extraction), multi-window scoring, batch optimization
- **Boltz-2**: Full MSA mode for publication-grade structures, RNA-DNA complex validation
- **Scientific Engine**: Real Enformer deployment, unified service discovery, end-to-end testing

---

## 📊 COMPONENT DEEP DIVES

### **1. Evo2 Service (`src/services/evo_service/main.py`)**

#### **Architecture:**
- **Modal App**: `evo-service`
- **GPU**: H100:2 (1B/7B), H100:1 (7B standalone), H100:1 (1B standalone)
- **Models**: `evo2_1b_base`, `evo2_7b`, `evo2_40b` (conditional)
- **Volume**: Persistent model cache (`evo-model-cache`)

#### **Endpoints (8 total):**
1. `/score_delta` - Delta log-likelihood (ref vs alt)
2. `/score_batch` - Batch delta scoring
3. `/score_variant` - Genomic SNV scoring (Ensembl fetch)
4. `/score_variant_multi` - Multi-window scoring (1024, 2048, 4096, 8192bp)
5. `/score_variant_exon` - Tight-window scoring (±600bp default)
6. `/score_variant_profile` - Local delta profile (±100bp radius)
7. `/score_variant_probe` - 3-alt sensitivity probe (A/C/G/T)
8. `/score_variant_with_activations` - **Layer 26 activations for SAE extraction** 🔥

#### **Strengths:**
- ✅ Comprehensive endpoint coverage (matches Evo2 paper capabilities)
- ✅ Ensembl integration for genomic context
- ✅ Reference allele validation (strict mismatch detection)
- ✅ Layer 26 activation extraction (SAE-ready)
- ✅ Graceful error handling with HTTPException
- ✅ Model fallback logic (1B → 7B → 40B)

#### **Gaps & Opportunities:**
- ⚠️ **No caching layer**: Every request hits Evo2 (expensive for repeated variants)
- ⚠️ **No batch optimization**: `/score_batch` processes sequentially
- ⚠️ **Activation extraction incomplete**: 7B version has different tokenization path
- 🔥 **PUSH TO**: 
  - Redis cache for variant scores (TTL: 24h)
  - True batch processing (parallel GPU inference)
  - SAE feature extraction pipeline (layer 26 → interpretable features)

#### **Integration Points:**
- **Backend Router**: `oncology-coPilot/oncology-backend-minimal/api/routers/evo.py`
- **Design Router**: `oncology-coPilot/oncology-backend-minimal/api/routers/design.py` (line 96: **HARDCODED URL**)
- **Metastasis Service**: `oncology-coPilot/oncology-backend-minimal/api/services/metastasis_interception_service.py`

---

### **2. Boltz-2 Service (`src/services/boltz_service/main.py`)**

#### **Architecture:**
- **Modal App**: `boltz-service`
- **GPU**: H100 (required)
- **Model**: `boltz-community/boltz-2` (Hugging Face)
- **Volume**: Persistent model cache (`boltz-models`)

#### **Endpoints (2 total):**
1. `/v1/predict_interaction` - Affinity prediction (target + candidate sequences)
2. `/v1/predict_structure` - **Structural integrity (FAST MODE)** ⚠️

#### **Critical Discovery: BOLTZ-2 ≠ ALPHAFOLD 3** 🚨

**Publication Claims**: "AlphaFold 3 structural validation"  
**Reality**: Uses Boltz-2 (different architecture, different training data)

**Evidence:**
- `src/services/boltz_service/main.py` line 59: `repo_id="boltz-community/boltz-2"`
- Boltz-2 is a **protein folding model** (not AF3's biomolecular complex predictor)
- Publication's `structural_validation_details.md` references "AlphaFold 3 Server" but code uses Boltz-2

**Impact:**
- ⚠️ **Scientific Misrepresentation**: Publication claims AF3 but uses Boltz-2
- ⚠️ **Different Capabilities**: Boltz-2 optimized for proteins, not RNA-DNA complexes
- ⚠️ **Threshold Mismatch**: Publication thresholds (pLDDT ≥50, iPTM ≥0.30) calibrated for AF3, not Boltz-2

#### **Fast Mode Trade-offs:**
```python
# Line 220-223: FAST MODE (msa='empty')
# Trade-off: Lower accuracy (pLDDT ~50-70) but FAST (2-5 min vs 60+ min)
final_input_data = {
    'version': 1,
    'sequences': [
        {'protein': {'id': 'TARG', 'sequence': protein_sequence, 'msa': 'empty'}},
    ],
}
```

**Current Implementation:**
- ✅ Fast: 2-5 minutes per structure
- ⚠️ Lower accuracy: pLDDT 50-70 (vs 70-90 with full MSA)
- ⚠️ Good for relative ranking, not absolute confidence

**Full MSA Mode (Not Implemented):**
- 🔥 **PUSH TO**: Full MSA mode for publication-grade structures
- Requires: `mmseqs2` MSA fetch (60+ min timeout)
- Expected: pLDDT 70-90, higher confidence for RNA-DNA complexes

#### **False Flag V2 Doctrine:**
```python
# Line 112-141: "False Flag V2" - Treats candidates as ligands
# This is a workaround to trigger interaction modeling
final_input_data = {
    'version': 1,
    'sequences': [
        {'protein': {'id': 'TARG', 'sequence': target_sequence, 'msa': str(sanitized_msa_path)}},
        {'ligand': {'id': candidate_id, 'smiles': smiles}}  # Candidate as ligand
    ],
    'properties': [
        {'affinity': {'binder': candidate_id}}
    ]
}
```

**Analysis:**
- ✅ Works for protein-protein interactions
- ⚠️ **Questionable for RNA-DNA complexes**: Boltz-2 not designed for nucleic acids
- 🔥 **PUSH TO**: Native RNA-DNA complex support (if Boltz-2 supports it) or migrate to AF3

#### **Strengths:**
- ✅ Fast mode enables rapid screening (2-5 min)
- ✅ MSA sanitization (handles ColabFold null bytes)
- ✅ Robust error handling (graceful degradation)

#### **Gaps & Opportunities:**
- 🔴 **Not AlphaFold 3**: Publication misrepresents model
- ⚠️ **No RNA-DNA native support**: Uses protein-ligand workaround
- ⚠️ **Fast mode only**: No full MSA option in production
- 🔥 **PUSH TO**:
  - Deploy real AlphaFold 3 Server integration (Google DeepMind API)
  - Full MSA mode for publication-grade structures
  - RNA-DNA complex native support

---

### **3. Scientific Engine (Orchestration Layer)**

#### **Architecture:**
- **Command Center**: `src/services/command_center/main.py`
- **Backend Routers**: `oncology-coPilot/oncology-backend-minimal/api/routers/`
- **Services**: `oncology-coPilot/oncology-backend-minimal/api/services/`

#### **Critical Issues:**

##### **1. Hardcoded Service URLs** 🚨
```python
# design.py line 96
evo_url = "http://127.0.0.1:8000/api/evo/score"
```

**Impact:**
- ❌ **Production Code**: Hardcoded localhost URL
- ❌ **No Environment Detection**: Doesn't use `EVO_URL_1B`/`EVO_URL_7B` from config
- ❌ **Breaks in Production**: Will fail when deployed

**Fix Required:**
```python
from ..config import get_model_url, DEFAULT_EVO_MODEL
evo_url = get_model_url(DEFAULT_EVO_MODEL) + "/score"
```

##### **2. Enformer Stubs (Confirmed)** 🔴
```python
# enformer_client.py: _stub_prediction()
# Deterministic fallback when ENFORMER_URL not set
base = 0.4 + (seed_int % 1000) / 1000 * 0.3
return {"accessibility_score": round(base, 3), ...}
```

**Impact:**
- ❌ **Publication Claims**: "Enformer/Borzoi" but uses stubs
- ❌ **15% of Target Lock Score**: Chromatin signal is noise
- ❌ **Inflated Metrics**: AUROC 0.976 includes fake data

**Fix Required:**
- Deploy real Enformer service (Modal or external API)
- Recompute all metrics without chromatin stubs
- Update publication to reflect 3-signal approach (not 4-signal)

##### **3. Fragmented Service Discovery**
- **Evo2 URLs**: Scattered across `config.py`, `index.py`, `evo.py`
- **Boltz URL**: Hardcoded in `command_center/main.py` as `BOLTZ_SERVICE_URL`
- **No Central Registry**: Each service discovers URLs independently

**Fix Required:**
- Unified service registry (environment-based)
- Health check endpoints for all services
- Automatic failover logic

#### **Integration Flow:**
```
User Request
    ↓
Backend Router (design.py, evo.py, etc.)
    ↓
Service Layer (metastasis_interception_service.py, safety_service.py)
    ↓
External Services (Evo2 Modal, Boltz Modal, Enformer [STUB])
    ↓
Response Assembly (Assassin Score, Target Lock, etc.)
```

**Gaps:**
- ⚠️ **No Circuit Breakers**: Services can fail silently
- ⚠️ **No Retry Logic**: Single-attempt calls
- ⚠️ **No Timeout Coordination**: Each service has independent timeouts

---

## 🧪 TEST STRATEGY

### **Current Test Coverage:**

#### **Evo2 Tests:**
- ✅ `tests/evo2/test_evo2_real_validation.py` - Real API validation (4 tests)
- ✅ `tests/evo2/test_evo2_paper_capabilities.py` - Paper capabilities (6 tests)
- ⚠️ **Missing**: Integration tests (Evo2 → Design → Assassin Score)

#### **Boltz Tests:**
- ❌ **No Tests Found**: No test files for Boltz-2 service
- ⚠️ **Missing**: Structural validation end-to-end tests

#### **Scientific Engine Tests:**
- ⚠️ **Partial**: Some unit tests, no end-to-end pipeline tests

### **Recommended Test Suite:**

#### **1. Evo2 Integration Tests** 🔥
```python
# tests/integration/test_evo2_design_pipeline.py
async def test_evo2_to_assassin_score_pipeline():
    """Test: Evo2 scoring → Design → Assassin Score"""
    # 1. Call Evo2 /score_variant
    # 2. Use delta in spacer efficacy endpoint
    # 3. Verify Assassin Score includes Evo2 signal
    # 4. Check provenance chain
```

#### **2. Boltz-2 Structural Validation Tests** 🔥
```python
# tests/integration/test_boltz_structural_validation.py
async def test_boltz_fast_vs_full_msa():
    """Compare fast mode vs full MSA mode"""
    # 1. Run fast mode (msa='empty')
    # 2. Run full MSA mode
    # 3. Compare pLDDT scores
    # 4. Validate threshold differences
```

#### **3. End-to-End Pipeline Tests** 🔥
```python
# tests/integration/test_metastasis_interception_pipeline.py
async def test_complete_metastasis_interception_flow():
    """Test: Target Lock → Guide Design → Evo2 → Boltz → Assassin Score"""
    # 1. Target Lock Algorithm (multi-modal signals)
    # 2. Guide RNA design
    # 3. Evo2 spacer efficacy scoring
    # 4. Boltz structural validation
    # 5. Assassin Score computation
    # 6. Verify all provenance links
```

#### **4. Service Discovery Tests** 🔥
```python
# tests/integration/test_service_discovery.py
async def test_service_url_resolution():
    """Test: Service URLs resolve correctly in all environments"""
    # 1. Check Evo2 URLs (1B, 7B, 40B)
    # 2. Check Boltz URL
    # 3. Check Enformer URL (should fail gracefully if stub)
    # 4. Verify health checks
```

---

## 🚀 IMPROVEMENT ROADMAP

### **P0 (Critical - Fix Immediately):**

1. **Fix Hardcoded URLs** 🔴
   - **File**: `oncology-coPilot/oncology-backend-minimal/api/routers/design.py`
   - **Line**: 96
   - **Fix**: Use `get_model_url()` from config
   - **Impact**: Production deployment will fail without this

2. **Clarify Boltz-2 vs AlphaFold 3** 🔴
   - **File**: `metastasis-interception/supplementary/structural_validation_details.md`
   - **Fix**: Update to reflect Boltz-2 usage (not AF3)
   - **Impact**: Scientific accuracy, publication integrity

3. **Deploy Real Enformer** 🔴
   - **File**: `oncology-coPilot/oncology-backend-minimal/api/services/enformer_client.py`
   - **Fix**: Deploy Enformer service (Modal or external API)
   - **Impact**: Remove 15% noise from Target Lock scores

### **P1 (High Priority - This Sprint):**

4. **Add Evo2 Caching Layer** 🔥
   - **Implementation**: Redis cache for variant scores
   - **TTL**: 24 hours
   - **Impact**: 10-100x cost reduction for repeated variants

5. **Implement Full MSA Mode for Boltz-2** 🔥
   - **Implementation**: Add `full_msa=True` parameter
   - **Timeout**: 1800s (30 min) for MSA fetch
   - **Impact**: Publication-grade structures (pLDDT 70-90)

6. **End-to-End Integration Tests** 🔥
   - **Implementation**: Test suite covering Evo2 → Design → Boltz → Assassin Score
   - **Impact**: Catch integration bugs before production

### **P2 (Medium Priority - Next Sprint):**

7. **SAE Feature Extraction Pipeline** 🔥
   - **Implementation**: Extract layer 26 activations → interpretable features
   - **Impact**: Explainable AI for guide RNA design

8. **Service Discovery Registry** 🔥
   - **Implementation**: Centralized service registry with health checks
   - **Impact**: Easier deployment, automatic failover

9. **Circuit Breakers & Retry Logic** 🔥
   - **Implementation**: Hystrix-style circuit breakers for external services
   - **Impact**: Resilience to service failures

### **P3 (Future Enhancements):**

10. **AlphaFold 3 Native Integration** 🔥
    - **Implementation**: Google DeepMind AF3 Server API
    - **Impact**: True RNA-DNA complex validation (not protein-ligand workaround)

11. **Batch Optimization for Evo2** 🔥
    - **Implementation**: Parallel GPU inference for batch requests
    - **Impact**: 10x throughput for large variant sets

12. **Multi-Modal Fusion Engine** 🔥
    - **Implementation**: Unified fusion of Evo2 + Enformer + Boltz signals
    - **Impact**: Better Assassin Scores, more accurate guide ranking

---

## 📈 PERFORMANCE BENCHMARKS

### **Evo2 Service:**
- **Latency**: 1-5 seconds per variant (1B model)
- **Throughput**: ~10 variants/second (1B), ~2 variants/second (7B)
- **Cost**: ~$0.001 per variant (1B), ~$0.01 per variant (7B)

### **Boltz-2 Service:**
- **Fast Mode**: 2-5 minutes per structure (pLDDT 50-70)
- **Full MSA Mode**: 60+ minutes per structure (pLDDT 70-90, not implemented)
- **Cost**: ~$0.10 per structure (H100 GPU)

### **Scientific Engine:**
- **Target Lock Algorithm**: <1 second (cached signals)
- **Assassin Score**: 2-10 seconds (Evo2 + Boltz calls)
- **End-to-End Pipeline**: 5-15 minutes (guide design → validation)

---

## 🎯 RECOMMENDATIONS FOR ALPHA

### **Immediate Actions:**
1. ✅ **Fix hardcoded URLs** (P0 - blocks production)
2. ✅ **Clarify Boltz-2 vs AF3** (P0 - scientific integrity)
3. ✅ **Deploy real Enformer** (P0 - removes 15% noise)

### **This Sprint:**
4. ✅ **Add Evo2 caching** (P1 - cost reduction)
5. ✅ **Full MSA mode for Boltz** (P1 - publication-grade)
6. ✅ **Integration tests** (P1 - quality assurance)

### **Next Sprint:**
7. ✅ **SAE pipeline** (P2 - explainability)
8. ✅ **Service registry** (P2 - deployment ease)
9. ✅ **Circuit breakers** (P2 - resilience)

### **Future:**
10. ✅ **AF3 native** (P3 - true RNA-DNA validation)
11. ✅ **Batch optimization** (P3 - throughput)
12. ✅ **Fusion engine** (P3 - accuracy)

---

## 🔥 WHAT WE CAN PUSH TO

### **Evo2:**
- **Layer 26 Activations**: ✅ Implemented, needs SAE pipeline
- **Multi-Window Scoring**: ✅ Implemented, needs optimization
- **Batch Processing**: ⚠️ Sequential, can parallelize
- **Caching**: ❌ Not implemented, 10-100x cost savings

### **Boltz-2:**
- **Full MSA Mode**: ❌ Not implemented, publication-grade structures
- **RNA-DNA Native**: ❌ Uses protein-ligand workaround
- **AlphaFold 3**: ❌ Not integrated, true RNA-DNA validation

### **Scientific Engine:**
- **Real Enformer**: ❌ Stubs only, 15% noise removal
- **Service Discovery**: ⚠️ Fragmented, needs central registry
- **End-to-End Tests**: ❌ Missing, critical for quality

---

**END OF AUDIT** ✅

**Next Steps**: Execute P0 fixes immediately, then P1 improvements this sprint.
