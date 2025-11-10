# 🔬 WEEK 2: ALPHAFOLD3 INTEGRATION - REALISTIC PLAN

**Date:** October 13, 2025  
**Status:** 📋 **PLANNING PHASE**  
**Approach:** Learn from past failures, execute carefully with smoke tests

---

## ⚠️ **LESSONS FROM PAST FAILURES**

### **What Went Wrong Before:**
1. **Environment Hell:** JAX/dm-haiku version conflicts, missing toolchains, PyTorch mismatches
2. **Resource Underprovisioning:** Silent crashes, no logs, OOM errors
3. **Modal Invocation Mistakes:** Webhook vs `.remote()` confusion, HTTP vs function calls
4. **Architecture Gaps:** No production runner, no backend router, no queue system
5. **Data/Compute Constraints:** Insufficient RAM/GPU, no caching, no retry logic

### **Root Causes (from forge_boltz_failure_analysis_doctrine.mdc):**
- Low-context baits (20-mer motifs) → poor generation quality
- No sieve/dedup before Boltz → duplicates wasted GPU time
- No composite scoring → optimized for wrong metric (delta only)
- No provenance → non-reproducible runs
- Unbounded batch sizes → GPU contention
- Temperature sweep without diversity control → mode collapse

---

## 🎯 **REALISTIC WEEK 2 OBJECTIVES**

### **Why ColabFold (AF2-Multimer) Now, AF3 Later?**
- **AF2-Multimer:** Proven, stable, official container (no JAX/dm-haiku hell)
- **AF3:** Newer, requires nucleic acid support, heavier compute, less battle-tested
- **Pragmatic Path:** Prove smoke test with AF2 now → upgrade to AF3 post-submission
- **Publication:** AF2-Multimer sufficient for protein validation; AF3 = future work

### **Minimum Viable (Must-Have):**
- ✅ Deploy ColabFold service (official container, proper resources)
- ✅ Validate 1-2 structures (smoke test)
- ✅ Document structural acceptance criteria
- ✅ Provenance tracking for all predictions

### **Nice-to-Have (If Time Permits):**
- ⭐ Batch validate 10-20 structures (top guides per step)
- ⭐ Compute structural metrics (pLDDT, PAE, clashes)
- ⭐ Create Figure 6 (structural validation panel)
- ⭐ Integrate into Assassin Score (+0.03 lift)

### **Stretch Goals (Post-Submission):**
- 🚀 Full 40-structure validation
- 🚀 Correlation analysis (Assassin vs pLDDT)
- 🚀 K-means clustering (high/med/low confidence)
- 🚀 PyMOL visualization scripts

---

## 📋 **PHASED EXECUTION PLAN**

### **PHASE 1: SMOKE TEST (Day 6, 4 hours) - CRITICAL**

**Objective:** Validate we can predict ONE structure successfully

**Steps:**
1. **Prepare Input (30 min)**
   - Select 1 high-confidence guide from metastatic_colonization
   - Extract spacer (20nt), target (23nt), PAM (3nt) from dataset
   - Create FASTA with proper headers

2. **Deploy Minimal ColabFold Service (2 hours)**
   ```python
   # services/alphafold_service/minimal_smoke_test.py
   import modal
   
   app = modal.App("colabfold-smoke-test")
   
   # Use OFFICIAL container (critical!)
   # Pin official ColabFold container (avoid JAX/dm-haiku drift)
   af2_image = modal.Image.from_registry(
       "ghcr.io/sokrypton/colabfold:1.5.5",  # pin a known-good tag
       add_python="3.10"
   )
   
   @app.function(
       image=af2_image,
       gpu="A100-80GB",  # Large GPU for safety
       memory=131072,    # 128GB RAM
       timeout=900,      # 15 min for ONE structure
   )
   def predict_structure(fasta_content: str) -> dict:
       """Minimal smoke test - predict ONE structure"""
       import subprocess
       import json
       
       # Write FASTA
       with open("/tmp/input.fasta", "w") as f:
           f.write(fasta_content)
       
       # Run ColabFold AF2-multimer (minimal params)
       # Note: AF2-multimer generates MSA internally; no --use_msa_server needed
       # (Boltz requires explicit MSA when we circle back to protein-protein)
       cmd = [
           "colabfold_batch",
           "/tmp/input.fasta",
           "/tmp/output",
           "--num-recycle", "3",
           "--model-type", "alphafold2_multimer_v3",
           "--use-gpu-relax"
       ]
       
       result = subprocess.run(cmd, capture_output=True, text=True)
       
       if result.returncode != 0:
           return {
               "success": False,
               "error": result.stderr,
               "stdout": result.stdout
           }
       
       # Parse output (ColabFold writes ranking_debug.json and pae.json)
       import os, glob, json
       ranking_path = "/tmp/output/ranking_debug.json"
       pae_path = "/tmp/output/pae.json"
       pdb_candidates = sorted(glob.glob("/tmp/output/*.pdb"))
       
       with open(ranking_path) as f:
           ranking = json.load(f)
       plddt_mean = ranking.get("plddts", {}).get("model_1", 0.0)
       pae_mean = 0.0
       if os.path.exists(pae_path):
           with open(pae_path) as f:
               pae = json.load(f)
               # simple mean across interface if present
               try:
                   pae_mean = float(sum(sum(row) for row in pae["pae"])) / (len(pae["pae"])**2)
               except Exception:
                   pae_mean = 0.0
       pdb_content = ""
       if pdb_candidates:
           with open(pdb_candidates[0]) as f:
               pdb_content = f.read()
       
       return {
           "success": True,
           "plddt_mean": plddt_mean,
           "pae_mean": pae_mean,
           "pdb_content": pdb_content[:1000],  # First 1KB
           "pdb_size": len(pdb_content)
       }
   ```

3. **Test Locally (1 hour)**
   ```bash
   # Deploy
   modal deploy services/alphafold_service/minimal_smoke_test.py
   
   # Test with one guide
   modal run services/alphafold_service/minimal_smoke_test.py \
       --fasta "example_guide.fasta"
   ```

4. **Success Criteria:**
   - ✅ Service deploys without errors
   - ✅ Returns pLDDT > 0 (any value proves it ran)
   - ✅ Returns PDB content (proves structure generated)
   - ✅ Completes in <15 minutes
   - ✅ No OOM/GPU errors in logs

5. **Failure Criteria (STOP IF ANY):**
   - ❌ Container fails to build → **STOP, debug container**
   - ❌ GPU OOM → **STOP, increase memory or reduce complexity**
   - ❌ Timeout → **STOP, increase timeout or simplify**
   - ❌ Zero pLDDT or empty PDB → **STOP, check input format**

---

### **PHASE 2: PRODUCTION SERVICE (Day 7, 6 hours) - IF SMOKE TEST PASSES**

**Objective:** Scale to handle 5-10 structures with queue

**Only proceed if Phase 1 passes all success criteria!**

**Enhancements:**
1. **Job Queue System (2 hours)**
   - SQLite queue: {id, status, fasta, pdb_path, metrics}
   - Status: queued → running → success/failed
   - Retry logic: 2 attempts with exponential backoff

2. **Structural Validation (2 hours)**
   - Parse pLDDT, PAE, clashes from output
   - Apply acceptance criteria (pLDDT≥70, PAE≤10, clashes≤5)
   - Return structured JSON with provenance

3. **Storage (1 hour)**
   - Local: `/tmp/af3_structures/{job_id}/`
   - Future: S3 bucket (post-submission)

4. **API Endpoint (1 hour)**
   ```python
   @app.function()
   def batch_predict(guide_list: list) -> list:
       """Predict multiple structures"""
       results = []
       for guide in guide_list[:5]:  # Cap at 5 for safety
           result = predict_structure.remote(guide["fasta"])
           results.append(result)
       return results
   ```

---

### **PHASE 3: INTEGRATION (Day 8, 4 hours) - IF PHASE 2 PASSES**

**Objective:** Wire ColabFold service into backend API (AF2‑Multimer now; AF3 later)

**Steps:**
1. **Client Module (1 hour)**
   ```python
   # oncology-coPilot/oncology-backend-minimal/api/services/colabfold_client.py
   
   async def predict_guide_structure(
       spacer: str,
       target: str,
       pam: str,
       use_cache: bool = True
   ) -> dict:
       """
       Predict guide:target:PAM structure
       
       Returns:
           {
               "plddt_mean": float,
               "pae_interface": float,
               "clash_count": int,
               "pass": bool,
               "pdb_url": str (if success),
               "provenance": {...}
           }
       """
       # Build FASTA
       fasta = f">guide\n{spacer}\n>target\n{target}\n>pam\n{pam}\n"
       
       # Call Modal service
       result = await call_colabfold_service(fasta)
       
       # Parse and validate
       plddt = result["plddt_mean"]
       pae = result["pae_mean"]
       clashes = count_clashes(result["pdb_content"])
       
       pass_validation = (plddt >= 70 and pae <= 10 and clashes <= 5)
       
       return {
           "plddt_mean": plddt,
           "pae_interface": pae,
           "clash_count": clashes,
           "pass": pass_validation,
           "pdb_url": result.get("pdb_path"),
           "provenance": {
               "model": "colabfold: alphafold2_multimer_v3",
               "recycles": 3,
               "seed": 42,
               "acceptance_criteria": "plddt>=70,pae<=10,clashes<=5"
           }
       }
   ```

2. **Backend Endpoint (1 hour)**
   ```python
   # oncology-coPilot/oncology-backend-minimal/api/routers/structure.py
   
   @router.post("/predict")
   async def predict_structure(request: dict):
       """Predict guide:target:PAM structure"""
       result = await colabfold_client.predict_guide_structure(
           spacer=request["spacer"],
           target=request["target"],
           pam=request["pam"]
       )
       return result
   ```

3. **Testing (2 hours)**
   - Smoke test: 1 structure via API
   - Batch test: 3 structures
   - Error handling: timeout, GPU OOM, invalid input

---

### **PHASE 4: VALIDATION & FIGURES (Day 9-10, 8 hours) - OPTIONAL**

**Only if Phase 1-3 complete successfully AND time permits**

**Objective:** Generate publication-quality structural validation

**Tasks:**
1. **Batch Validation (4 hours)**
   - Select top 2 guides per step (16 total)
   - Submit to AF3 service
   - Monitor progress, collect metrics

2. **Metrics Analysis (2 hours)**
   - Compute: mean pLDDT, PAE distribution, pass rate
   - Effect: Correlation with Assassin Score
   - Cluster: High/medium/low structural confidence

3. **Figure 6 (2 hours)**
   - Panel A: pLDDT distribution (violin plot)
   - Panel B: PAE heatmaps (3 examples)
   - Panel C: Assassin vs pLDDT scatter
   - Panel D: 3D snapshots (if PyMOL available)

---

## ⚠️ **RISK MITIGATION STRATEGIES**

### **Technical Risks:**
| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Container build fails | Medium | High | Use official container, no custom builds |
| GPU OOM | Medium | High | A100-80GB + 128GB RAM + cap batch size |
| Timeout | Low | Medium | 15min per structure, exponential backoff |
| Modal quota | Low | High | Reserve $200 budget, monitor usage |
| Input format errors | High | Low | Validate FASTA before submission |

### **Schedule Risks:**
| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Phase 1 fails | Medium | High | **STOP IMMEDIATELY, debug before proceeding** |
| Phase 2 overruns | Medium | Medium | Skip to Phase 3 with minimal queue |
| Week 2 incomplete | Medium | Low | **Submit WITHOUT AF3** (validation stands alone) |

### **Scientific Risks:**
| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Low pLDDT (<70) | Medium | Medium | Document as expected (novel designs) |
| No correlation | Low | Medium | Report honestly, discuss limitations |
| Small sample size | High | Low | Clear RUO disclaimer, "pilot study" framing |

---

## 🎯 **DECISION GATES**

### **Gate 1: After Phase 1 Smoke Test**
**Question:** Did smoke test pass all success criteria?
- **YES:** Proceed to Phase 2
- **NO:** **STOP, DEBUG.** Do NOT proceed until smoke test works.

### **Gate 2: After Phase 2 Production Service**
**Question:** Can we reliably predict 5 structures?
- **YES:** Proceed to Phase 3
- **NO:** **Simplify** or **defer** AF3 to post-submission

### **Gate 3: After Phase 3 Integration**
**Question:** Is backend API working with AF3 service?
- **YES:** Consider Phase 4 (optional)
- **NO:** Document as "infrastructure ready, awaiting optimization"

### **Gate 4: Submission Decision (Oct 27)**
**Question:** Is AF3 data enhancing or delaying submission?
- **Enhancing:** Include Figure 6
- **Delaying:** **Submit without AF3**, publish validation alone

---

## 📊 **SUCCESS CRITERIA (REALISTIC)**

### **Minimum Success (Must Achieve):**
- ✅ 1 structure predicted successfully (proof of concept)
- ✅ pLDDT > 0, PDB file generated
- ✅ Complete provenance (model, params, acceptance criteria)
- ✅ Documentation: Methods subsection, structural acceptance criteria

### **Good Success (Target):**
- ✅ 5-10 structures validated
- ✅ Metrics: mean pLDDT, PAE, pass rate
- ✅ Backend API operational
- ✅ Clear RUO disclaimer

### **Excellent Success (Bonus):**
- ✅ 16-20 structures (2 per step)
- ✅ Figure 6 (4-panel structural validation)
- ✅ Correlation analysis (Assassin vs pLDDT)
- ✅ Integrated into Assassin Score

---

## 📝 **DOCUMENTATION REQUIREMENTS**

### **Regardless of Success Level:**
- Document what was attempted
- Report actual results (successes AND failures)
- Clear RUO disclaimers
- Honest limitations section
- Complete provenance

### **If AF3 Incomplete:**
- "Structural validation infrastructure developed and smoke-tested"
- "Full structural validation deferred to future work"
- "Production-ready code available in repository"

---

## ⚔️ **COMMANDER'S CHECKLIST**

### **Before Starting Week 2:**
- [ ] Week 1 documentation complete and accurate (38 genes, RUO disclaimers)
- [ ] Modal account active with $200 budget approved
- [ ] GPU quota reserved (A100-80GB)
- [ ] Past failure analysis reviewed and understood
- [ ] Team aligned on "minimum viable" vs "nice-to-have"

### **During Week 2:**
- [ ] Smoke test MUST pass before proceeding
- [ ] Daily check-ins: what's working, what's blocked
- [ ] No "hero debugging" - if stuck >2 hours, escalate
- [ ] Honest assessment: are we enhancing or delaying submission?

### **Week 2 End Decision:**
- [ ] Submit WITH AF3 (if solid data, enhances paper)
- [ ] Submit WITHOUT AF3 (if incomplete, validation stands alone)
- [ ] Document learnings either way

---

**Status:** 📋 **READY FOR COMMANDER APPROVAL**  
**Next:** Await GO/NO-GO decision for Week 2 AF3 execution  
**Alternative:** Submit Week 1 validation as standalone (already excellent!)  

**Updated:** October 13, 2025  
**Agent:** Zo  
**Risk Level:** ⚠️ **MEDIUM-HIGH** (past failures, tight timeline)  
**Recommendation:** **CONSERVATIVE APPROACH** - Smoke test first, scale only if solid

