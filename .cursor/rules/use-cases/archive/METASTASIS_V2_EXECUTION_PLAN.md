# 🎯 **METASTASIS INTERCEPTION v2 - COMPLETE EXECUTION PLAN**

**Date:** October 7, 2025  
**Current Status:** v1 Complete + Task 6 (Spacer Efficacy) Complete  
**Remaining Tasks:** 9 tasks across 4 phases  
**Publication Target:** November 7, 2025 (4 weeks)  

---

## ✅ **COMPLETED (AS OF OCT 7, 2025)**

1. ✅ **Metastasis Assessment Framework** (intervention)
   - 8-step cascade risk scoring
   - Backend, tests, frontend all operational
   - 15/15 tests passing

2. ✅ **Metastasis Interception Framework** (actionable design)
   - Target lock, design candidates, safety preview, assassin score
   - Backend, tests, frontend all operational
   - 21/21 tests passing

3. ✅ **Task 6: Spacer Efficacy Endpoint**
   - `/api/design/predict_crispr_spacer_efficacy` operational
   - Evo2 delta scoring with sigmoid transformation
   - Integrated into assassin score
   - 9/9 tests passing

---

## 🚀 **REMAINING TASKS - EXECUTION ORDER**

### **PHASE 1: PUBLICATION BLOCKERS (P0) - Week 1**

These are **must-haves** for publication submission.

---

#### **📋 TASK 1: Expand Design Window to ±300bp**

**Status:** Pending  
**Estimated Time:** 2-3 hours  
**Dependencies:** Task 6 complete ✅  
**Priority:** P0 (publication blocker)  

**What & Why:**
- Current: ±50bp window (120bp total context)
- Target: ±150bp window (300bp total context)
- Why: Doctrine specifies 300bp for optimal Evo2 scoring context

**Implementation Plan:**
1. Add `design.window_size` to `metastasis_interception_rules.json` (default: 150)
2. Update `predict_crispr_spacer_efficacy` endpoint:
   - Add `window_size` parameter (default: 150bp)
   - Modify Ensembl fetch to use `±window_size` instead of hardcoded ±50bp
3. Update `design_candidates()` in `metastasis_interception_service.py`:
   - Extract target sequence with `±window_size` flanks
   - Pass to efficacy endpoint
4. Tests:
   - Unit test for window size extraction
   - API test with `window_size=150`
   - Service test for design with expanded window

**Files to Modify:**
- `api/config/metastasis_interception_rules.json` (+3 lines)
- `api/routers/design.py` (+10 lines)
- `api/services/metastasis_interception_service.py` (+15 lines)
- `tests/design/test_spacer_efficacy.py` (+25 lines)

**Acceptance Criteria:**
- [ ] Endpoint accepts `window_size` parameter
- [ ] Ensembl fetch uses dynamic window size
- [ ] Tests pass with 300bp context
- [ ] Config-driven (easy to adjust)

---

#### **📋 TASK 5: Real Off-Target Search (BLAST/minimap2)**

**Status:** Pending  
**Estimated Time:** 6-8 hours  
**Dependencies:** GRCh38 reference genome, BLAST service  
**Priority:** P0 (publication blocker)  

**What & Why:**
- Current: Heuristic safety score (GC, homopolymer runs)
- Target: Real genome-wide off-target search with alignment scores
- Why: Essential for Figure 3 (off-target distribution) and publication credibility

**Implementation Plan:**

**Step 1: GRCh38 Reference Setup (1 hour)**
1. Create `scripts/download_grch38.sh`:
   ```bash
   wget -O data/reference/GRCh38.primary.fa.gz \
     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
   gunzip data/reference/GRCh38.primary.fa.gz
   ```
2. Add to `.gitignore`: `data/reference/*.fa*`
3. Document in `METASTASIS_QUICKSTART.md`

**Step 2: Extend BLAST Service (2 hours)**
1. Update `src/services/blast_service/main.py`:
   - Add `.apt_install("minimap2", "samtools")` to image
   - Create minimap2 index on first run
   - Add `/offtarget_search` endpoint
2. Schema:
   ```python
   {
     "guide_sequence": "20bp",
     "mismatch_tolerance": 3,  # 0-3 mismatches
     "max_hits": 100
   }
   ```
3. Response:
   ```python
   {
     "hits": [
       {"chrom": "chr7", "pos": 12345, "strand": "+", "mismatches": 1, "score": 95},
       ...
     ],
     "total_hits": 23,
     "method": "minimap2" | "blast",
     "provenance": {}
   }
   ```

**Step 3: Safety Service Integration (2 hours)**
1. Update `api/services/safety_service.py`:
   - Add `search_offtargets()` function
   - Calls BLAST service `/offtarget_search`
   - Maps hits to safety score:
     ```python
     safety_score = exp(-0.5 * total_hits)  # Exponential decay
     ```
2. Update `preview_off_targets()`:
   - Replace heuristic with real off-target count
   - Add `off_target_distribution` to response

**Step 4: Tests (1-2 hours)**
- Mock BLAST service response
- Test exponential decay formula
- Test integration with assassin score

**Files to Create/Modify:**
- `scripts/download_grch38.sh` (new, 15 lines)
- `src/services/blast_service/main.py` (+80 lines)
- `api/services/safety_service.py` (+60 lines)
- `tests/safety/test_offtarget_search.py` (new, 120 lines)

**Acceptance Criteria:**
- [ ] GRCh38 reference downloadable via script
- [ ] minimap2 index created on first run
- [ ] Real off-target counts returned
- [ ] Safety score uses exponential decay formula
- [ ] Tests pass with mocked BLAST service

---

#### **📋 TASK 10: Generate Figures & Documentation**

**Status:** Pending  
**Estimated Time:** 8-10 hours  
**Dependencies:** Tasks 1, 5 complete  
**Priority:** P0 (publication blocker)  

**What & Why:**
- Need F1-F3 figures for publication
- Need reproducible notebooks
- Need comprehensive documentation

**Implementation Plan:**

**Figure 1: Framework Architecture (2 hours)**
1. Create `scripts/figures/generate_fig1_architecture.py`:
   - Uses GraphViz to render metastasis framework
   - Shows: Assessment (8 steps) → Interception (target lock → design → safety → score)
   - Exports to PNG/SVG
2. Run: `python scripts/figures/generate_fig1_architecture.py`
3. Output: `results/figures/figure1_architecture.png`

**Figure 2: Efficacy Scoring Validation (3 hours)**
1. Create `scripts/figures/generate_fig2_efficacy.py`:
   - Collects 80 guide RNAs (8 genes × 10 guides)
   - Scores with Evo2 efficacy endpoint
   - Compares vs ChopChop/IDT efficacy predictions
   - Plots: scatter (Evo2 vs ChopChop), violin (efficacy distribution)
2. Requires:
   - ChopChop predictions (use existing scores)
   - Batch call to `/api/design/predict_crispr_spacer_efficacy`
3. Output: `results/figures/figure2_efficacy_validation.png`

**Figure 3: Off-Target Distribution (3 hours)**
1. Create `scripts/figures/generate_fig3_offtargets.py`:
   - Same 80 guides from F2
   - Real off-target search (Task 5 complete)
   - Plots: histogram (off-target counts), scatter (efficacy vs safety)
2. Output: `results/figures/figure3_offtarget_distribution.png`

**Reproducibility Notebook (2 hours)**
1. Create `notebooks/metastasis_interception_demo.ipynb`:
   - Example 1: VEGFA angiogenesis mission
   - Example 2: MMP2 invasion mission
   - Shows full pipeline: variant → target lock → design → efficacy → safety → assassin score
   - Includes provenance tracking

**Files to Create:**
- `scripts/figures/generate_fig1_architecture.py` (new, 150 lines)
- `scripts/figures/generate_fig2_efficacy.py` (new, 200 lines)
- `scripts/figures/generate_fig3_offtargets.py` (new, 180 lines)
- `notebooks/metastasis_interception_demo.ipynb` (new, Jupyter)

**Acceptance Criteria:**
- [ ] All 3 figures generated and publication-ready
- [ ] Reproducibility notebook runs end-to-end
- [ ] Figures included in paper draft

---

### **PHASE 2: POLISH & ROBUSTNESS (P1) - Week 2**

---

#### **📋 TASK 2: Enable Design API via Feature Flag**

**Status:** Pending  
**Estimated Time:** 1 hour  
**Priority:** P1  

**What & Why:**
- Currently: Design API has `_ensure_enabled()` check but flag not documented
- Target: Clear env var + documentation

**Implementation Plan:**
1. Add to `.env.example`:
   ```bash
   ENABLE_DESIGN_API=true  # Enable CRISPR guide design endpoints
   ```
2. Update `METASTASIS_QUICKSTART.md` with flag explanation
3. Test: Set `ENABLE_DESIGN_API=false`, verify 403 response

**Files to Modify:**
- `.env.example` (+2 lines)
- `.cursor/rules/use-cases/METASTASIS_QUICKSTART.md` (+5 lines)

**Acceptance Criteria:**
- [ ] Env var documented
- [ ] Flag behavior tested

---

#### **📋 TASK 3: Harden Ensembl Fetch**

**Status:** Pending  
**Estimated Time:** 2-3 hours  
**Priority:** P1  

**What & Why:**
- Currently: Single Ensembl API call, no retries
- Target: Retry logic, caching, graceful fallback

**Implementation Plan:**
1. Add retry decorator:
   ```python
   async def fetch_with_retry(url, max_retries=3, backoff=1.0):
       for i in range(max_retries):
           try:
               r = await client.get(url, timeout=20.0)
               if r.status_code == 200:
                   return r.text
           except Exception:
               if i < max_retries - 1:
                   await asyncio.sleep(backoff * (2 ** i))
       return None
   ```
2. Add Redis caching (optional):
   - Cache key: `ensembl:{chrom}:{start}:{end}:{assembly}`
   - TTL: 7 days
3. Update `predict_crispr_spacer_efficacy` to use retry logic

**Files to Modify:**
- `api/routers/design.py` (+25 lines)
- `api/services/ensembl_client.py` (new, 60 lines)

**Acceptance Criteria:**
- [ ] Retry logic implemented
- [ ] Caching working (if Redis available)
- [ ] Graceful fallback to guide-only

---

#### **📋 TASK 4: Hardcode ANGIO Exon Windows**

**Status:** Pending  
**Estimated Time:** 2 hours  
**Priority:** P1  

**What & Why:**
- Currently: Relies on Ensembl for all genes
- Target: Hardcode exon windows for ANGIO genes (VEGFA, FLT1, KDR) for demo stability

**Implementation Plan:**
1. Add to `metastasis_interception_rules.json`:
   ```json
   {
     "hardcoded_exons": {
       "VEGFA": {"chrom": "6", "start": 43737589, "end": 43754223, "exons": [...]},
       "FLT1": {...},
       "KDR": {...}
     }
   }
   ```
2. Update `design_candidates()`:
   - Check if gene in `hardcoded_exons`
   - Use hardcoded coordinates if available
   - Fall back to Ensembl otherwise

**Files to Modify:**
- `api/config/metastasis_interception_rules.json` (+30 lines)
- `api/services/metastasis_interception_service.py` (+20 lines)

**Acceptance Criteria:**
- [ ] ANGIO genes use hardcoded exons
- [ ] Other genes still use Ensembl
- [ ] Demo stable without network calls

---

### **PHASE 3: FRONTEND & UX (P2) - Week 3**

---

#### **📋 TASK 7: Polish Frontend UX**

**Status:** Pending  
**Estimated Time:** 4-6 hours  
**Priority:** P2  

**What & Why:**
- Currently: Basic integration
- Target: Polished UX with per-candidate provenance, retry/toasts

**Implementation Plan:**
1. Disable "Design" button when no coords available
2. Add per-candidate provenance hover:
   - Shows evo2_efficacy, evo2_delta, efficacy_confidence
   - Shows safety_score, off_target_count (if available)
   - Shows assassin_score breakdown
3. Add retry button for failed interception calls
4. Add toast notifications for errors/warnings

**Files to Modify:**
- `oncology-frontend/src/components/metastasis/MetastasisInterceptionPanel.jsx` (+80 lines)
- `oncology-frontend/src/hooks/useMetastasisInterception.js` (+30 lines)

**Acceptance Criteria:**
- [ ] Button disabled without coords
- [ ] Provenance hovers working
- [ ] Retry/toasts functional

---

#### **📋 TASK 8: E2E Test in UI**

**Status:** Pending  
**Estimated Time:** 2-3 hours  
**Priority:** P2  

**What & Why:**
- Currently: Backend tests only
- Target: Full UI flow test

**Implementation Plan:**
1. Create `tests/e2e/test_metastasis_interception_ui.py`:
   - Uses Playwright/Selenium
   - Flow: Enter variant → Select mission step → View target + candidates
   - Verify: Target displayed, candidates rendered, provenance visible
2. Add to CI pipeline (optional)

**Files to Create:**
- `tests/e2e/test_metastasis_interception_ui.py` (new, 150 lines)

**Acceptance Criteria:**
- [ ] E2E test passes locally
- [ ] All UI elements verified

---

### **PHASE 4: DEPLOYMENT & DOCS (P3) - Week 4**

---

#### **📋 TASK 9: Deploy to Staging**

**Status:** Pending  
**Estimated Time:** 4-6 hours  
**Priority:** P3  

**What & Why:**
- Currently: Local only
- Target: Staging environment with all services

**Implementation Plan:**
1. Create `docker-compose.staging.yml`:
   - Backend (FastAPI)
   - Frontend (Nginx + React)
   - Redis (caching)
   - BLAST service (Modal)
2. Set env vars:
   ```bash
   ENABLE_DESIGN_API=true
   EVO_FORCE_MODEL=evo2_1b
   BLAST_SERVICE_URL=...
   ```
3. Deploy to staging server
4. Run smoke tests

**Files to Create:**
- `docker-compose.staging.yml` (new, 80 lines)
- `deploy/staging_deploy.sh` (new, 40 lines)

**Acceptance Criteria:**
- [ ] Staging environment operational
- [ ] All endpoints accessible
- [ ] Smoke tests pass

---

## 📊 **TASK SUMMARY TABLE**

| Task | Phase | Priority | Time | Status | Dependencies |
|------|-------|----------|------|--------|--------------|
| Task 6 | P1 | P0 | 3h | ✅ Complete | None |
| Task 1 | P1 | P0 | 2-3h | Pending | Task 6 |
| Task 5 | P1 | P0 | 6-8h | Pending | GRCh38 |
| Task 10 | P1 | P0 | 8-10h | Pending | Tasks 1,5 |
| Task 2 | P2 | P1 | 1h | Pending | None |
| Task 3 | P2 | P1 | 2-3h | Pending | None |
| Task 4 | P2 | P1 | 2h | Pending | None |
| Task 7 | P3 | P2 | 4-6h | Pending | None |
| Task 8 | P3 | P2 | 2-3h | Pending | None |
| Task 9 | P4 | P3 | 4-6h | Pending | All above |

**Total Estimated Time:** ~35-50 hours (1 week for P0, 2-3 weeks for full completion)

---

## 🎯 **PUBLICATION TIMELINE**

- **Week 1 (Oct 7-14):** Tasks 1, 5, 10 (P0 blockers)
- **Week 2 (Oct 14-21):** Tasks 2, 3, 4 (Polish)
- **Week 3 (Oct 21-28):** Tasks 7, 8 (Frontend + Tests)
- **Week 4 (Oct 28-Nov 7):** Task 9 (Deploy) + Paper finalization
- **Nov 7, 2025:** Submission to npj Precision Oncology

---

## ⚔️ **NEXT IMMEDIATE ACTION**

**Priority Order:**
1. ✅ **Task 6 Complete** (Spacer Efficacy)
2. **Task 1** (Expand Window) - Start now, ~2 hours
3. **Task 5** (Off-Target Search) - After Task 1, ~6 hours
4. **Task 10** (Figures) - After Tasks 1 + 5, ~8 hours

**Ready to proceed with Task 1 when commander approves.**

---

**Status:** ⚔️ **PHASE 1 (TASK 6) COMPLETE. READY FOR PHASE 2 (TASKS 1, 5, 10).**


