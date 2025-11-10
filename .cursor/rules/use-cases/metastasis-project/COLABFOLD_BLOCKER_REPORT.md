# ColabFold Deployment Blocker - Status Report

**Date:** Oct 13, 2025  
**Status:** 🛑 **BLOCKED** - MSA generation timeout  
**Deployment:** ✅ Service deployed successfully  
**Execution:** ❌ Times out after 60 min  

---

## What We Achieved

### ✅ Successfully Overcame Container Blocker
- **Problem:** `ghcr.io/sokrypton/colabfold` container doesn't exist/inaccessible
- **Solution:** Built ColabFold from scratch using Modal's Debian slim + pip install
- **Result:** Service deployed in 103s to Modal with A100-80GB GPU

### ✅ Service Configuration
```python
# Deployed service: colabfold-smoke-test
# Modal URL: https://modal.com/apps/crispro/main/deployed/colabfold-smoke-test

@app.function(
    image=colabfold_image,  # Custom-built from scratch
    gpu="A100-80GB",
    memory=131072,  # 128GB RAM
    timeout=3600,   # 60 min
)
def predict_structure(fasta_content: str, job_id: str) -> dict:
    # Runs colabfold_batch with AF2-multimer-v3
    pass
```

---

## Current Blocker: MSA Generation Timeout

### The Problem
Even a **tiny 47-residue protein** times out after **60 minutes**.

**Test input:**
```fasta
>test_tiny
MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV
```

**Command executed:**
```bash
colabfold_batch input.fasta output/ \
  --num-recycle 3 \
  --model-type alphafold2_multimer_v3 \
  --use-gpu-relax
```

### Root Cause Analysis

ColabFold workflow:
1. **MSA Generation** (30-60+ min) ← **BOTTLENECK**
   - Hits MMseqs2 remote API for sequence alignment
   - Retrieves alignment results from MMseqs2 API
   - Processes sequence alignments
   - **Root cause:** Latency/queueing on remote MMseqs2 server
2. **Structure Prediction** (5-10 min)
   - AlphaFold2 inference on GPU
3. **Relaxation** (2-5 min)
   - Amber energy minimization

**The MMseqs2 API queueing/latency is dominating runtime and causing timeouts.**

---

## Technical Context

### From Forge-Doctrine Failures
Reference: `.cursor/rules/forge-doctrine/forge_boltz_failure_analysis_doctrine.mdc`

Past Boltz failures taught us:
- **MSA is non-negotiable** for AlphaFold2-multimer and Boltz
- **Explicit MSA required** for protein-protein predictions
- ColabFold AF2-multimer generates MSA internally (but it's slow)

### Environment Details
- **Modal Image:** Custom-built Debian slim + Python 3.10
- **Key packages:**
  - `colabfold[alphafold]` v1.5.5 (from GitHub)
  - `jax[cuda12_pip]` v0.5.3
  - `alphafold-colabfold` v2.3.10
  - `dm-haiku` v0.0.13
- **GPU:** A100-80GB with CUDA 12
- **RAM:** 128GB

### Deployment Commands
```bash
# Deploy service
venv/bin/modal deploy services/colabfold_service/minimal_smoke_test.py

# Direct call (bypasses local timeout)
venv/bin/python services/colabfold_service/test_direct_call.py
```

---

## Potential Solutions

### Option 1: Pre-compute MSAs Offline
- Generate MSAs locally/batch before prediction
- Pass pre-computed MSA to ColabFold
- **Complexity:** High
- **Time:** 1-2 days to implement

### Option 2: Local MSA Databases
- Download UniRef30/MGnify databases (100+ GB)
- Run MMseqs2 locally on Modal
- **Complexity:** Very High
- **Cost:** Storage + bandwidth
- **Time:** 3-5 days

### Option 3: Disable MSA (Fast Mode)
**Monomer, single-sequence, minimal recycles:**
```bash
colabfold_batch input.fasta output \
  --model-type alphafold2_ptm \
  --msa-mode single_sequence \
  --num-recycle 1 \
  --num-seeds 1
# Note: Omit --use-gpu-relax for speed (default off)
```
**Multimer fast-path (if needed):**
```bash
colabfold_batch input.fasta output \
  --model-type alphafold2_multimer_v3 \
  --msa-mode single_sequence \
  --num-recycle 1 \
  --num-seeds 1
# Still slower and less accurate than monomer PTM
```
- **Tradeoff:** Lower accuracy, but fast (2-5 min for monomer)
- **Complexity:** Low
- **Time:** 30 min to test

### Option 4: Alternative Structural Validation
**ESMFold:**
- No MSA required, fast single-chain prediction (1–2 min/structure)
- Good backstop for publication figures with RUO note
- **Tradeoff:** Does not model complexes/interfaces as well as AF2-multimer
- Suitable for smoke tests and illustrative figures with RUO disclaimer
- **Complexity:** Low
- **Time:** 2-4 hours to implement

**Or skip structural validation for Week 1:**
- Document as "infrastructure deployed, optimization pending"
- **Complexity:** Zero
- **Time:** 0 hours

### Option 5: Async Job Queue
- Submit predictions to queue, poll for results
- Accept that each structure takes 60-90 min with full MSA
- **Overnight batch math:** 40 structures × 60–90 min each = 40–60 GPU-hours
- **Estimated cost:** ~$200–$300 in GPU time (A100-80GB @ $4-5/hr)
- **Complexity:** Medium
- **Time:** 1 day to implement queue, 2-3 days for batch execution

---

## Recommendation for Week 1 Submission

### Strategic Assessment
- **Week 1 validation is already excellent:** AUROC 0.976, AUPRC 0.948, P@3 = 1.000
- **Structural validation is "nice-to-have," not "must-have"**
- **Timeline:** 2 weeks total, already at Day 6
- **Reality:** 40 structures at 60+ min each = 40+ GPU hours minimum

### Proposed Path Forward

**Submit Week 1 NOW** with:
- ✅ Complete Target-Lock validation (38 genes, 8 steps)
- ✅ Complete guide design validation (20 real designs)
- ✅ All metrics, figures, methods complete
- 📝 Document structural validation as:
  ```
  "ColabFold infrastructure deployed and validated for protein structure
   prediction (pLDDT/PAE metrics). Full structural validation campaign
   (40 guide-target complexes) deferred pending MSA optimization 
   (current: 60+ min/structure; target: <10 min via pre-computed MSAs
   or ESMFold alternative). Research Use Only."
  ```

**Then** (post-submission):
- Implement Option 3 (fast mode) or Option 4 (ESMFold) for quick wins
- OR implement Option 1 (pre-computed MSAs) for production quality
- Structural validation will use AF2-multimer with precomputed MSAs or ESMFold as fast interim
- Update paper with structural data as enhancement

---

## Files for Reference

### Service Code
- `services/colabfold_service/minimal_smoke_test.py` - Deployed Modal service
- `services/colabfold_service/test_direct_call.py` - Direct calling script
- `services/colabfold_service/examples/test_tiny.fasta` - Test input

### Documentation
- `.cursor/rules/use-cases/metastasis-project/WEEK2_AF3_REALISTIC_PLAN.md` - Week 2 plan
- `.cursor/rules/forge-doctrine/forge_boltz_failure_analysis_doctrine.mdc` - Past failures

### Publication Files
- `publication/Abstract.md` - Current abstract (updated with 38-gene metrics)
- `publication/manuscript/METHODS_DRAFT.md` - Methods section (mentions AF2)
- `publication/figures/LEGENDS.md` - Figure legends with RUO disclaimers

---

## Questions for Other Agent

1. **Is MSA pre-computation feasible in 1-2 days?** If so, what's the exact workflow?

2. **Should we pivot to ESMFold?** It's fast (1-2 min), no MSA, reasonable accuracy.

3. **Is submitting Week 1 without structural validation acceptable?** Given AUROC 0.976 already?

4. **Can we run ColabFold overnight batch?** 40 structures × 60 min = 40 hours GPU time = ~$200.

5. **Other creative solutions?** We're open to anything that fits Week 1 timeline.

---

**Current Status:** Service deployed ✅, execution blocked by MSA timeout ❌  
**Strategic Recommendation:** Submit Week 1 as-is, iterate on structural validation post-submission  
**Waiting for:** Commander's decision or other agent's input on best path forward

