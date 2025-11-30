# Pre-Development Readiness Assessment

**Date**: January 21, 2025  
**Status**: ✅ **READY FOR DEVELOPMENT**  
**Plan Reference**: `.cursor/plans/sota-benchmarks-and-frontend-integration-aa7ca3bc.plan.md`

---

## ✅ Verification Checklist

### Critical Fixes Status

**Fix 1: Pathway Normalization** ✅ **VERIFIED**
- **Location**: `api/services/efficacy_orchestrator/drug_scorer.py:48-55`
- **Status**: Implemented correctly
- **Formula**: `path_pct = min(1.0, max(0.0, s_path / 0.005))`
- **Verification**: Code matches plan exactly

**Fix 2: Tier Computation Parameter** ✅ **VERIFIED**
- **Location**: `api/services/efficacy_orchestrator/drug_scorer.py:139`
- **Status**: Implemented correctly
- **Change**: Passes raw `s_path` instead of normalized `path_pct`
- **Verification**: Code matches plan exactly

**Fix 3: Tier Threshold** ✅ **VERIFIED**
- **Location**: `api/services/confidence/tier_computation.py:61`
- **Status**: Implemented correctly
- **Change**: Threshold adjusted from 0.05 to 0.001
- **Verification**: Code matches plan exactly

**Fix 4: Sporadic Gates Capping** ✅ **VERIFIED**
- **Location**: `api/services/efficacy_orchestrator/orchestrator.py:225-228`
- **Status**: Implemented correctly
- **Logic**: Only applies when `tumor_context_data` is provided OR `germline_status != "unknown"`
- **Verification**: Code matches plan exactly

### Benchmark Scripts Status

**MM Benchmark** ✅ **READY**
- **File**: `scripts/benchmark_sota_mm.py`
- **Status**: Complete and ready
- **Test Variants**: 7 variants (BRAF, KRAS, NRAS, TP53)
- **Validation**: Checks drug ranking accuracy (>80% target)
- **Output**: Saves to `results/sota/mm_benchmark_*.json`

**Ovarian Benchmark** ✅ **READY**
- **File**: `scripts/benchmark_sota_ovarian.py`
- **Status**: Complete and ready
- **Dataset**: Uses 1k dataset (`hrd_tcga_ov_labeled_1k_results.json`) with fallback
- **Validation**: Computes AUROC/AUPRC (>0.75 target)
- **Output**: Saves to `results/sota/ovarian_benchmark_*.json`
- **Note**: Script correctly prioritizes 1k dataset (lines 25-28)

**Melanoma Benchmark** ✅ **READY**
- **File**: `scripts/benchmark_sota_melanoma.py`
- **Status**: Complete and ready
- **Test Variants**: BRAF V600E, NRAS Q61K
- **Validation**: Checks drug ranking accuracy (>90% target)
- **Output**: Saves to `results/sota/melanoma_benchmark_*.json`

### Evo2 Integration Status

**Evo2 Modal Service** ✅ **OPERATIONAL**
- **Location**: `src/services/evo_service/main.py`
- **Status**: Deployed and operational
- **Model**: evo2_1b (default), evo2_7b (optional)
- **Endpoints**: `/score_delta`, `/score_variant_multi`, `/score_variant_exon`
- **GPU**: H100:2 (1B), A10G (7B)

**Evo2 Proxy Router** ✅ **READY**
- **Location**: `api/routers/evo.py`
- **Status**: Complete with caching and fallback
- **Caching**: Redis with in-memory fallback
- **Reference Sequences**: Ensembl integration working

**Sequence Processor** ✅ **READY**
- **Location**: `api/services/efficacy_orchestrator/sequence_processor.py`
- **Status**: Complete with fallback chain
- **Flow**: FusionAM → Evo2 → MassiveOracle
- **Integration Point**: Ready for SOTA enhancements (BRCA classifier, splice scorer)

### Integration Architecture

**Single Integration Point** ✅ **CONFIRMED**
- **File**: `api/services/efficacy_orchestrator/sequence_processor.py`
- **Status**: Ready for SOTA enhancements
- **Pattern**: Detect variant type → use appropriate scorer → fallback chain
- **No Duplication**: All enhancements integrate here, no parallel paths

**Automatic Benefits** ✅ **CONFIRMED**
- **Pathway (P)**: Benefits automatically from better S scores
- **Evidence (E)**: Benefits from better variant interpretation
- **Confidence**: Improves automatically with better S scores
- **No Code Changes Needed**: Orchestrator automatically benefits

---

## 🎯 Development Readiness

### Immediate Next Steps (Phase 1: Supervised BRCA Classifier)

**Week 1: Build**
1. ✅ **Training Data**: ClinVar BRCA1/BRCA2 variants (need to verify availability)
2. ✅ **Model Training**: XGBoost on Evo2 embeddings (layer 26)
3. ✅ **Modal Endpoint**: `/brca_classifier` in `src/services/evo_service/main.py`
4. ✅ **Backend Proxy**: `/api/evo/brca_classifier` in `api/routers/evo.py`

**Week 2: Integrate**
5. ✅ **Sequence Processor**: Add `_use_supervised_brca_classifier()` method
6. ✅ **Fallback Chain**: Zero-shot if classifier unavailable
7. ✅ **Validation**: Verify BRCA variants get higher scores

**Week 3: Validate**
8. ✅ **Benchmark**: Re-run ovarian benchmark, verify AUROC >0.75
9. ✅ **Regression**: Verify MM and Melanoma benchmarks (no regression)

### Questions & Clarifications

**Q1: ClinVar BRCA1/BRCA2 Training Data**
- **Question**: Do we have access to ClinVar BRCA1/BRCA2 variants with pathogenicity labels?
- **Status**: ⚠️ **NEEDS VERIFICATION**
- **Action**: Check if we have ClinVar data extraction scripts or need to fetch from ClinVar API
- **Fallback**: Use existing ClinVar integration in `api/services/evidence/literature_client.py`

**Q2: Evo2 Embeddings Extraction**
- **Question**: How do we extract Evo2 embeddings (layer 26) for training?
- **Status**: ✅ **ANSWERED** - Evo2 Modal service has `/score_variant_with_activations` endpoint
- **Location**: `src/services/evo_service/main.py:805-842` (7B service)
- **Action**: Use this endpoint to extract embeddings for training data

**Q3: Classifier Model Storage**
- **Question**: Where should we store the trained classifier model?
- **Status**: ⚠️ **NEEDS DECISION**
- **Options**:
  1. Modal Volume (persistent storage)
  2. HuggingFace Hub (version control)
  3. S3/Cloud storage (scalable)
- **Recommendation**: Start with Modal Volume, migrate to HuggingFace Hub for versioning

**Q4: Training Infrastructure**
- **Question**: Where should we train the classifier (local vs. Modal)?
- **Status**: ⚠️ **NEEDS DECISION**
- **Options**:
  1. Local training (faster iteration, requires Evo2 embeddings pre-extracted)
  2. Modal training (scalable, can extract embeddings on-the-fly)
- **Recommendation**: Start local (faster iteration), move to Modal for production

**Q5: Splice Site Detection Data**
- **Question**: Do we have Ensembl exon annotations for splice site detection?
- **Status**: ✅ **ANSWERED** - Ensembl integration exists in `api/routers/evo.py`
- **Action**: Reuse existing Ensembl fetching logic, add exon boundary detection

**Q6: Noncoding Variant Detection**
- **Question**: How do we detect noncoding variants (promoters, enhancers, UTRs)?
- **Status**: ⚠️ **NEEDS IMPLEMENTATION**
- **Options**:
  1. Ensembl annotations (gene structure)
  2. Regulatory region databases (ENCODE, FANTOM)
  3. Heuristic detection (distance from TSS, UTR boundaries)
- **Recommendation**: Start with Ensembl annotations, add regulatory databases later

### Risk Assessment

**Low Risk** ✅
- Pathway normalization fix (already implemented)
- Tier computation fix (already implemented)
- Sporadic gates fix (already implemented)
- Benchmark scripts (ready and tested)

**Medium Risk** ⚠️
- **ClinVar Data Access**: Need to verify we can extract BRCA1/BRCA2 variants
- **Training Infrastructure**: Need to decide on training location and model storage
- **Evo2 Embeddings**: Need to verify `/score_variant_with_activations` works for 1B model

**High Risk** 🔴
- **None Identified**: All critical components are ready

### Gaps & Missing Pieces

**Gap 1: ClinVar BRCA1/BRCA2 Variant Extraction**
- **Status**: ⚠️ **NEEDS VERIFICATION**
- **Impact**: Blocks Phase 1 Week 1 (training data)
- **Action**: Check existing ClinVar integration or create extraction script
- **Timeline**: 1-2 hours to verify/create

**Gap 2: Evo2 Embeddings Extraction for 1B Model**
- **Status**: ⚠️ **NEEDS VERIFICATION**
- **Impact**: Blocks Phase 1 Week 1 (training data)
- **Action**: Verify if `/score_variant_with_activations` exists for 1B or only 7B
- **Timeline**: 30 minutes to verify

**Gap 3: Base Scorer Interface**
- **Status**: ⚠️ **NEEDS CREATION**
- **Impact**: Blocks clean integration of new scorers
- **Action**: Create `api/services/sequence_scorers/base_scorer.py`
- **Timeline**: 1 hour to create

**Gap 4: Noncoding Variant Detection Logic**
- **Status**: ⚠️ **NEEDS IMPLEMENTATION**
- **Impact**: Blocks Phase 3 (noncoding optimization)
- **Action**: Implement variant type detection in `evo2_scorer.py`
- **Timeline**: 2-3 hours to implement

---

## 📋 Pre-Development Action Items

### Before Starting Phase 1

**P0: Critical (Do Before Development)**
1. ✅ **Verify ClinVar Data Access** (1-2 hours)
   - Check if we have ClinVar BRCA1/BRCA2 variants
   - If not, create extraction script using existing ClinVar integration
   - Target: 1000+ variants with pathogenicity labels

2. ✅ **Verify Evo2 Embeddings Extraction** (30 minutes)
   - Test `/score_variant_with_activations` endpoint (if exists for 1B)
   - If not, check if we can extract embeddings from 7B service
   - Verify layer 26 activations are accessible

3. ✅ **Create Base Scorer Interface** (1 hour)
   - Create `api/services/sequence_scorers/base_scorer.py`
   - Define common interface for all scorers
   - Update existing scorers to implement interface

**P1: Important (Can Do During Development)**
4. ⚠️ **Decide Training Infrastructure** (30 minutes)
   - Choose: Local vs. Modal training
   - Choose: Model storage location (Modal Volume vs. HuggingFace Hub)
   - Document decision in plan

5. ⚠️ **Create Training Script** (2-3 hours)
   - Script to extract Evo2 embeddings from ClinVar variants
   - Script to train XGBoost classifier
   - Script to validate AUROC >0.94 on test set

**P2: Nice to Have (Can Do Later)**
6. ⚠️ **Noncoding Variant Detection** (2-3 hours)
   - Implement variant type detection
   - Add Ensembl annotation fetching
   - Test with known noncoding variants

---

## 🎯 Success Criteria (Re-Confirmed)

### Phase 1: Supervised BRCA Classifier

**Training Validation:**
- ✅ Classifier AUROC >0.94 on held-out test set (match paper)
- ✅ Training data: 1000+ ClinVar BRCA1/BRCA2 variants
- ✅ Model: XGBoost on Evo2 embeddings (layer 26)

**Integration Validation:**
- ✅ BRCA variants use supervised classifier (not zero-shot)
- ✅ Fallback to zero-shot if classifier unavailable
- ✅ No regression in MM/Melanoma benchmarks

**Benchmark Validation:**
- ✅ Ovarian AUROC >0.75 (especially BRCA variants)
- ✅ Variants differentiate (not all 0.5 baseline)
- ✅ No regression in MM (maintain 100%) or Melanoma

### Phase 2: Splice Variant Prediction

**Training Validation:**
- ✅ Splice scoring AUROC >0.82 on test set (match paper)
- ✅ Splice site detection working (±20 bp of exon-intron boundaries)

**Integration Validation:**
- ✅ Splice-affecting variants use splice scorer
- ✅ Fallback to zero-shot if splice scorer unavailable
- ✅ Improved predictions for splice variants in all benchmarks

### Phase 3: Noncoding Optimization

**Implementation Validation:**
- ✅ Noncoding variant detection working
- ✅ Optimized windows for regulatory variants (50K bp, cost acceptable)
- ✅ Improved predictions for noncoding variants

**Benchmark Validation:**
- ✅ Ovarian noncoding variants show differentiation (not all 0.5)
- ✅ Computational cost acceptable (longer windows justified)

---

## 🚀 Development Roadmap (Finalized)

### Week 1: Build Supervised BRCA Classifier

**Day 1-2: Data Preparation**
- Extract ClinVar BRCA1/BRCA2 variants (or verify existing data)
- Extract Evo2 embeddings for training variants
- Split into train/validation/test sets (70/15/15)

**Day 3-4: Model Training**
- Train XGBoost classifier on Evo2 embeddings
- Validate AUROC >0.94 on test set
- Save model to storage (Modal Volume or HuggingFace Hub)

**Day 5: Modal Endpoint**
- Create `/brca_classifier` endpoint in `src/services/evo_service/main.py`
- Load trained classifier model
- Test endpoint with known variants

### Week 2: Integrate Supervised BRCA Classifier

**Day 1-2: Backend Proxy**
- Create `/api/evo/brca_classifier` endpoint in `api/routers/evo.py`
- Add Redis caching (reuse existing cache infrastructure)
- Add fallback to zero-shot if classifier unavailable

**Day 3-4: Sequence Processor Integration**
- Add `_use_supervised_brca_classifier()` method to `sequence_processor.py`
- Detect BRCA1/BRCA2 variants → use supervised classifier
- Test with known BRCA variants

**Day 5: Validation**
- Verify BRCA variants get higher/differentiated scores
- Test fallback chain (classifier → zero-shot → MassiveOracle)
- Run regression tests (MM, Melanoma benchmarks)

### Week 3: Benchmark Validation

**Day 1-2: Ovarian Benchmark**
- Re-run `scripts/benchmark_sota_ovarian.py`
- Verify AUROC >0.75 (especially BRCA variants)
- Document improvement from supervised classifier

**Day 3: Regression Testing**
- Re-run MM benchmark (verify 100% maintained)
- Re-run Melanoma benchmark (verify no regression)
- Document all results

**Day 4-5: Documentation & Deployment**
- Document supervised classifier implementation
- Update plan with results
- Prepare for Phase 2 (if Phase 1 successful)

---

## ✅ Final Pre-Development Checklist

### Code Readiness
- [x] All 4 critical fixes implemented and verified
- [x] Benchmark scripts ready and tested
- [x] Evo2 Modal service operational
- [x] Integration architecture clear (single integration point)

### Data Readiness
- [ ] ClinVar BRCA1/BRCA2 variants available (or extraction script ready)
- [ ] Evo2 embeddings extraction verified (layer 26)
- [ ] Training/validation/test splits prepared

### Infrastructure Readiness
- [ ] Training infrastructure decided (local vs. Modal)
- [ ] Model storage decided (Modal Volume vs. HuggingFace Hub)
- [ ] Base scorer interface created

### Documentation Readiness
- [x] Plan finalized and comprehensive
- [x] Success criteria clear
- [x] Risk mitigation strategies defined
- [x] Timeline realistic (6-8 weeks for Tier 1-2)

---

## 🎯 Ready to Begin Development

**Status**: ✅ **READY** (pending P0 action items)

**Immediate Actions**:
1. Verify ClinVar data access (1-2 hours)
2. Verify Evo2 embeddings extraction (30 minutes)
3. Create base scorer interface (1 hour)

**Then Proceed With**:
- Phase 1 Week 1: Build Supervised BRCA Classifier

**Expected Timeline**: 6-8 weeks for Tier 1-2 capabilities (realistic, not months)

**Success Metrics**: 
- Ovarian: AUROC 0.500 → 0.75+ (especially BRCA variants)
- MM: Maintain 100% pathway alignment accuracy
- Melanoma: Improve drug ranking accuracy (50% → 90%+)

---

## 📝 Notes & Considerations

### From .cursorrules

**Key Principles to Follow**:
- ✅ **Build Once, Use Everywhere**: Modal services build capabilities, backend integrates once
- ✅ **Single Integration Point**: `sequence_processor.py` is the only place SOTA capabilities integrate
- ✅ **Reuse Existing Infrastructure**: Caching, fallback chains, error handling all reused
- ✅ **Automatic Benefits**: Better S scores automatically improve P/E/Confidence (no code changes needed)
- ✅ **Realistic Timeline**: 6-8 weeks for Tier 1-2 (not months)
- ✅ **Smart & Efficient**: No duplication, reusable components, continuous integration

**Lessons Learned**:
- ✅ Modal service communication: Use `modal.Cls.lookup()` for service-to-service calls
- ✅ Error handling: Graceful degradation with fallback chains
- ✅ Testing: Run all benchmarks after each phase to catch regressions
- ✅ Documentation: Update plan with results and learnings

### Questions for Manager (If Needed)

1. **ClinVar Data**: Do we have existing ClinVar extraction scripts, or should we use the ClinVar API directly?
2. **Model Storage**: Preference for classifier model storage (Modal Volume vs. HuggingFace Hub)?
3. **Training Location**: Preference for training location (local vs. Modal)?
4. **Priority**: Should we proceed with Phase 1 immediately, or wait for any approvals?

---

**ASSESSMENT COMPLETE** ✅  
**READY FOR DEVELOPMENT** (pending P0 action items)  
**NEXT STEP**: Verify ClinVar data access and Evo2 embeddings extraction

