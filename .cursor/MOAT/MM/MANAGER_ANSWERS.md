# Manager Questions - Zo's Answers for Plumber

**Date:** January 29, 2025  
**Answered by:** Zo  
**Context:** Comprehensive answers based on TRUE SAE validation results, code audit, and validated artifacts

---

## 🔴 CRITICAL STRATEGIC QUESTIONS

### **Question 1: TRUE SAE Validation Strategy**

**If MM TRUE SAE validation fails to reach AUROC ≥ 0.70, what is the fallback strategy?**

**Answer: B (primary) + A (fallback)**

**Reasoning:**
1. **Primary (B):** Iterate on diamond mining methodology. OV found 9 features → DDR_bin with AUROC 0.783. MM may need different feature selection (UPR_bin, Proteasome_bin) rather than DDR_bin. The methodology works—we just need to find MM-specific diamond features.

2. **Fallback (A):** If iteration fails after reasonable effort (2-3 attempts), ship with Proxy SAE (gene-level). DIS3 (RR=2.08, p=0.0145) and TP53 (RR=1.90) are already validated on MMRF GDC. This is production-ready.

**Why NOT C (cross-disease transfer):** OV and MM have different resistance mechanisms. OV = DDR/platinum resistance. MM = proteasome/IMiD/anti-CD38 resistance. Cross-disease transfer is biologically unsound.

**Why NOT D (delay):** Proxy SAE is validated and clinically useful. Delaying for TRUE SAE is unnecessary when we have a working production solution.

---

### **Question 2: Priority Trade-off**

**If we can only complete ONE work item before a critical deadline, which do we choose?**

**Answer: B (PSMB5/CRBN) - Work Item 2**

**Reasoning:**
1. **Work Item 2 (2-3 hours):** Adds immediate, clinically relevant functionality. PSMB5 → PI resistance and CRBN → IMiD resistance are the #1 clinical questions for MM oncologists. Even literature-based, this is useful.

2. **Work Item 1 (2-3 days):** Higher long-term value but unproven for MM. We don't know if MM TRUE SAE will achieve AUROC ≥ 0.70. Risk of spending 2-3 days and failing.

**Trade-off:** Immediate utility beats speculative long-term value under time pressure. We can always add TRUE SAE later as enhancement.

**Evidence:** OV TRUE SAE took months of iteration to validate. MM may have similar learning curve.

---

### **Question 3: Validation Data Dependency**

**If MMRF CoMMpass access is delayed and we only have cBioPortal data (500 patients)?**

**Answer: D (Hybrid validation) → A (proceed with smaller cohort)**

**Reasoning:**
1. **D first:** Combine cBioPortal (500 patients) + literature for hybrid validation. Literature gives us expected RR values (PSMB5→PI, CRBN→IMiD from published studies). cBioPortal gives us actual patient data for power calculation.

2. **A if D insufficient:** Proceed with smaller cohort. N=500 is still sufficient for detecting RR≥2.0 effects (power ≈80% for 10% mutation rate). Our validated DIS3 analysis used N=219 and achieved p=0.0145.

**Why NOT B (wait):** 1-2 week delay is significant for a 3-week project. Time-to-market matters.

**Why NOT C (literature only):** We need SOME cohort validation to be credible. Literature-only is acceptable for rare mutations (PSMB5, CRBN) but not for the entire system.

---

## 🟡 TECHNICAL DEEP DIVE QUESTIONS

### **Question 4: Diamond Feature Selection**

**If MM diamond mining finds features mapping to multiple pathways, how do we handle it?**

**Answer: D (multi-pathway is a feature, not a problem) + A (prioritize by effect size)**

**Reasoning:**
1. **MM biology is multi-pathway.** Unlike OV (DDR-dominant), MM has distinct resistance mechanisms per drug class:
   - Proteasome inhibitors → Proteasome/UPR pathway
   - IMiDs → Cereblon pathway
   - Anti-CD38 → CD38 expression pathway

2. **Prioritize by effect size (A):** Use |Cohen's d| ≥ 0.5 threshold (same as OV). Features with largest effect sizes go into drug-class-specific bins.

3. **Create drug-class bins:**
   - `Proteasome_bin`: Features from diamond mining on PI-resistant vs PI-sensitive
   - `Cereblon_bin`: Features from diamond mining on IMiD-resistant vs IMiD-sensitive

**Validation:** Each bin should have pathway coherence like OV's 9/9 DDR. If not, that bin may lack signal.

---

### **Question 5: PSMB5/CRBN Mutation Rarity**

**If validation shows "INSUFFICIENT_DATA" for PSMB5/CRBN (too few mutations)?**

**Answer: D (mark as LITERATURE_ONLY) + C (combine with pathway-level signals)**

**Reasoning:**
1. **Reality check:** PSMB5 and CRBN resistance mutations are genuinely rare (n=2-3). We cannot fake statistical power.

2. **D: Mark as LITERATURE_ONLY.** Be transparent. Users see: "PSMB5 mutation detected → PI resistance predicted (RR=3.5, LITERATURE_ONLY)". This is clinically useful even without cohort validation.

3. **C: Combine with pathway-level.** If PSMB5 is mutated, also compute `proteasome_pathway_burden` (PSMB5 + PSMB8 + NFE2L2 + XBP1). Pathway-level has more mutations → more power.

**Why NOT A (accept unvalidated):** We must be transparent about evidence quality.

**Why NOT B (broader gene-level):** PSMB5 gene-level is what we're already doing. Specific variants (p.Ala49Thr) are literature-based.

---

### **Question 6: TRUE SAE vs. Proxy SAE Performance**

**If MM TRUE SAE achieves AUROC 0.72 but Proxy SAE achieves 0.68, is 0.04 improvement worth it?**

**Answer: Yes, conditionally**

**Reasoning:**
1. **0.04 AUROC improvement IS clinically meaningful** at this range. AUROC 0.68→0.72 moves from "fair" to "good" discrimination.

2. **But context matters:**
   - If TRUE SAE provides better interpretability (DDR_bin-like pathway mapping) → **YES, worth it**
   - If TRUE SAE is a black box (no pathway coherence) → **NO, stick with Proxy SAE**

3. **Cost-benefit:**
   - 2-3 days development: Acceptable for a 3-week project (20-30% of time)
   - $20-30 Modal cost: Trivial
   - Added complexity: Mitigated by keeping Proxy SAE as fallback
   - External dependency: Modal is already proven reliable (used for OV)

**Decision rule:** If TRUE SAE achieves AUROC ≥ 0.70 AND has pathway coherence (like OV's DDR_bin), use it. Otherwise, stick with Proxy SAE.

---

## 🟢 ARCHITECTURAL & DESIGN QUESTIONS

### **Question 7: Pathway Service Integration**

**How do we integrate pathway-level signals (Work Item 3) with TRUE SAE diamond features (Work Item 1)?**

**Answer: C (use pathways to validate TRUE SAE coherence) → A (use pathways as feature bins)**

**Reasoning:**
1. **C first:** When diamond mining finds features, map them to pathways. Check if features cluster by pathway (like OV's 9/9 DDR coherence). This validates that TRUE SAE is capturing biology, not noise.

2. **A second:** Once validated, use pathways as bins. Create `Proteasome_bin`, `Cereblon_bin`, etc. from diamond features that map to those pathways.

3. **Integration architecture:**
   ```
   TRUE SAE features → Pathway mapping → Pathway bins → Mechanism vector
                                                      ↓
   Proxy SAE (gene-level) → Pathway burden ──────────────→ Combined score
   ```

**Not B (independent):** Keeping them separate wastes synergy.

**Not D (fallback only):** They complement each other—TRUE SAE for variant-level, pathways for gene-level.

---

### **Question 8: Validation Test Failure Handling**

**If only 2/5 tests pass, which tests are acceptable to fail?**

**Answer priority (most acceptable to fail → least acceptable):**

1. **Test 1 (PSMB5) - ACCEPTABLE TO FAIL** ⚠️
   - Rare mutations (n=2-3), INSUFFICIENT_DATA is expected
   - Mark as LITERATURE_ONLY, still clinically useful

2. **Test 2 (CRBN) - ACCEPTABLE TO FAIL** ⚠️
   - Same as PSMB5—rare mutations
   - Mark as LITERATURE_ONLY

3. **Test 3 (del(17p)) - ACCEPTABLE TO FAIL** ⚠️
   - Cytogenetics data may be missing from MMRF GDC extract
   - Literature HR=2.5 is well-established (IMWG guidelines)
   - Mark as LITERATURE_ONLY

4. **Test 5 (TRUE SAE vs. PROXY) - ACCEPTABLE TO FAIL** ⚠️
   - If TRUE SAE doesn't beat Proxy SAE, we still have Proxy SAE (production-ready)
   - This is a "nice to have" comparison, not a launch blocker

5. **Test 4 (DIS3+TP53) - MUST PASS** 🔴
   - Already validated (DIS3 p=0.0145)
   - If this fails, something is broken in the system
   - This is our validated production baseline

**Minimum viable:** Test 4 MUST pass. Tests 1-3 can be LITERATURE_ONLY. Test 5 is optional.

---

### **Question 9: Frontend Integration Timing**

**If Work Item 1 (TRUE SAE) is delayed but Work Items 2-3 complete?**

**Answer: C (build frontend with TRUE SAE placeholder)**

**Reasoning:**
1. **Build frontend now.** Users need to see resistance predictions. Don't block UI on TRUE SAE.

2. **Design for future TRUE SAE:** Include placeholder sections:
   - "Mechanism Breakdown" → Show pathway bins (empty or Proxy-based if no TRUE SAE)
   - "Feature Insights" → Hidden section, enable when TRUE SAE ready

3. **Progressive enhancement:** Ship Proxy SAE first → Add TRUE SAE overlay later → Users see improvement without breaking change.

**Why NOT A (Proxy only):** Limits future flexibility.

**Why NOT B (wait):** Delays user value.

**Why NOT D (two versions):** Maintenance nightmare.

---

## 🔵 RISK & MITIGATION QUESTIONS

### **Question 10: Modal SAE Service Failure**

**If Modal SAE service is unavailable during MM extraction?**

**Answer: B (fall back to Proxy SAE) → D (build local pipeline for long-term)**

**Reasoning:**
1. **Immediate (B):** Fall back to Proxy SAE (gene-level). This is production-ready and doesn't depend on Modal. Ship with Proxy SAE if Modal is down.

2. **Long-term (D):** Build local SAE extraction pipeline. This removes external dependency. Priority: P2 (after launch).

**Why NOT A (cross-disease):** Biologically unsound (OV DDR ≠ MM proteasome).

**Why NOT C (delay):** Unacceptable for launch timeline.

---

### **Question 11: MMRF Data Access Failure**

**If both MMRF and cBioPortal access fail?**

**Answer: A (use existing MMRF GDC subset N=219)**

**Reasoning:**
1. **We already have validated results on MMRF GDC (N=219):**
   - DIS3: RR=2.08, p=0.0145 ✅
   - TP53: RR=1.90 ✅
   - This is sufficient for launch.

2. **For TRUE SAE:** We need TRUE SAE features extracted, which requires mutation data. MMRF GDC has mutations for 219 patients. Use this subset.

3. **Degrade gracefully:**
   - Proxy SAE: Validated (DIS3, TP53)
   - TRUE SAE: Extract on N=219 (may have lower power, but better than nothing)

**Why NOT B (literature only):** We have real data—use it.

**Why NOT C (delay):** We have MMRF GDC already.

**Why NOT D (synthetic):** Synthetic data is for pipeline testing, not validation.

---

### **Question 12: Performance Regression**

**If TRUE SAE integration causes performance regression (lower AUROC than Proxy SAE)?**

**Answer: C (use ensemble) → B (debug if significant)**

**Reasoning:**
1. **C first:** Try ensemble (TRUE SAE + Proxy SAE combined). Often, different signals are complementary. Ensemble may outperform both individually.

2. **B if ensemble fails:** Debug TRUE SAE integration. Common issues:
   - Feature selection leakage (cross-validation bug)
   - Wrong aggregation method (mean vs. max)
   - Pathway mapping errors

3. **A if all else fails:** Revert to Proxy SAE only. Don't ship a regression.

**Why NOT D (accept for interpretability):** Performance regression is unacceptable for clinical tool. Interpretability doesn't justify worse predictions.

---

## 🟣 BUSINESS & PRODUCT QUESTIONS

### **Question 13: Launch Criteria**

**What is the minimum viable launch criteria?**

**Answer: C (validation framework runs, 2/5 tests pass) + D (backend only acceptable)**

**Reasoning:**
1. **C is acceptable IF:**
   - Test 4 (DIS3+TP53) passes ✅ (our validated baseline)
   - Test 1-3 failures are due to rare mutations/missing data (expected)
   - We document failed tests as LITERATURE_ONLY

2. **D is acceptable:**
   - Backend-only launch is viable for API consumers
   - Frontend can follow in Phase 2
   - Many users (developers, researchers) prefer API access

3. **A (AUROC 0.68) is NOT acceptable:**
   - Below 0.70 threshold
   - Would need to ship as "experimental" or "research only"

4. **B (PSMB5/CRBN INSUFFICIENT_DATA) IS acceptable:**
   - Mark as LITERATURE_ONLY
   - Still clinically useful

---

### **Question 14: Competitive Advantage**

**If competitors launch MM resistance with Proxy SAE before we complete TRUE SAE?**

**Answer: C (launch with both: Proxy SAE now, TRUE SAE later)**

**Reasoning:**
1. **Don't delay launch.** First-mover advantage matters. Competitors with Proxy SAE are at parity—we can match immediately.

2. **TRUE SAE as differentiator for V2.** Once TRUE SAE is validated, release as upgrade. This becomes marketing story: "Now with TRUE SAE—25% better discrimination."

3. **Other differentiators matter:**
   - Pathway service (6 MM-specific pathways)
   - Drug-class specific predictions (PI, IMiD, anti-CD38, BCMA)
   - Frontend UX (resistance panel with next-line options)

**Why NOT A (abandon TRUE SAE):** Loses long-term advantage.

**Why NOT B (delay):** Loses market position.

---

### **Question 15: Resource Allocation**

**Option A (Work Items 1+2) vs. Option B (Work Items 2+3+5)?**

**Answer: Option B (Work Items 2+3+5 - Full product)**

**Reasoning:**
1. **Option B delivers a complete product:**
   - PSMB5/CRBN mutations (Work Item 2): Core clinical functionality
   - Pathway service (Work Item 3): Mechanism-level insights
   - Frontend (Work Item 5): User-facing value

2. **Option A is incomplete:**
   - TRUE SAE (Work Item 1): May not validate for MM
   - PSMB5/CRBN (Work Item 2): Good, but no UI to show it
   - No frontend = no user value

3. **User value > technical excellence.** A complete product with Proxy SAE beats an incomplete product with TRUE SAE.

4. **TRUE SAE can be added later** as enhancement without breaking existing functionality.

---

## 🔴 AUDIT DISCREPANCY QUESTIONS

### **Question 16: Mission Claims vs. Reality**

**PSMB5/CRBN "defined" but not implemented. Evo2 "working" but not integrated.**

**Answer: B (trust code audit) + D (update mission document)**

**Reasoning:**
1. **Code is source of truth.** The audit found:
   - PSMB5/CRBN: ❌ NOT in `resistance_prophet_service.py`
   - Evo2: ✅ Service exists, ❌ NOT called from `predict_mm_resistance()`

2. **Mission document is aspirational.** "Defined" likely means "defined in doctrine/rules" not "implemented in code."

3. **Action:** Update mission document to reflect actual state:
   - PSMB5/CRBN: "❌ Not implemented (Work Item 2)"
   - Evo2: "⚠️ Service exists, not integrated (P2)"

---

### **Question 17: TRUE SAE Validation Contradiction**

**AUROC 0.783 (good) vs. 0 FDR-significant features (bad) vs. DDR_bin p=0.0020 (significant)?**

**Answer: A (methodology difference) + C (pathway aggregation)**

**Reasoning:**
1. **These answer DIFFERENT questions:**

   | Analysis | Question | Answer |
   |----------|----------|--------|
   | FDR biomarker scan | "Does any SINGLE feature survive multiple testing?" | No (0 features) |
   | AUROC baseline | "Can a MULTIVARIATE feature set separate labels?" | Yes (0.783) |
   | DDR_bin p=0.0020 | "Does a PATHWAY BIN separate labels?" | Yes |

2. **Reconciliation:**
   - Individual features: Noise (underpowered, 24 resistant vs. 125 sensitive)
   - Combined features (29 total): Signal (AUROC 0.783)
   - Aggregated to pathway: Strong signal (DDR_bin p=0.0020)

3. **This is GOOD NEWS:** Pathway-level aggregation is biologically meaningful and overcomes individual feature noise.

**Updated SAE_READINESS_STATUS.md now explains this explicitly.**

---

### **Question 18: Proxy SAE Production Ready vs. TRUE SAE Requirements**

**If Proxy SAE is production-ready, why does mission require TRUE SAE?**

**Answer: A (ship Proxy SAE now, add TRUE SAE later)**

**Reasoning:**
1. **Proxy SAE is validated and sufficient for launch:**
   - DIS3: RR=2.08, p=0.0145 ✅
   - TP53: RR=1.90 ✅
   - Gene-level is clinically interpretable

2. **TRUE SAE is competitive advantage, not launch requirement:**
   - OV validation proves methodology works
   - MM validation is additional effort (2-3 days)
   - Can be added as V2 enhancement

3. **Mission requirement is aspirational.** Reframe as:
   - V1: Proxy SAE (production-ready)
   - V2: TRUE SAE (when validated)

---

### **Question 19: Data Source Discrepancy**

**MMRF GDC (N=219) vs. full MMRF CoMMpass (N=1,154)?**

**Answer: A (accept MMRF GDC as sufficient for Proxy SAE) + B (require full cohort for TRUE SAE)**

**Reasoning:**
1. **For Proxy SAE (gene-level):** N=219 is sufficient.
   - DIS3 achieved p=0.0145 on N=219
   - Power is adequate for RR≥2.0 effects

2. **For TRUE SAE (32K features):** N=1,154 is better.
   - More patients = more resistant cases = better power
   - TRUE SAE feature selection benefits from larger cohort

3. **Practical approach:**
   - Launch V1 with Proxy SAE (validated on N=219)
   - Run TRUE SAE extraction on full cohort when access available
   - Add TRUE SAE in V2

---

### **Question 20: Architecture: Dedicated Service vs. Prophet Service**

**Create `MMResistanceService` (dedicated) vs. keep logic in `ResistanceProphetService` (disease-agnostic)?**

**Answer: C (refactor to hybrid)**

**Reasoning:**
1. **Current architecture is acceptable:** `ResistanceProphetService` already has `predict_mm_resistance()` method. This is disease-agnostic design.

2. **Create MM-specific helpers:**
   ```python
   # In resistance_prophet_service.py
   def predict_mm_resistance(self, ...):
       # MM-specific logic
       return self._mm_helper.compute_mm_signals(...)
   
   # In mm_resistance_helpers.py (new file)
   class MMResistanceHelper:
       def compute_mm_signals(self, ...):
           # PSMB5, CRBN, cytogenetics, etc.
   ```

3. **Why NOT A (dedicated service):** Duplicates routing/orchestration logic.

4. **Why NOT B (keep in prophet):** MM-specific logic will grow. Helpers keep prophet service clean.

5. **Why NOT D (accept mismatch):** Architecture debt compounds. Fix it now.

---

### **Question 21: Validation Standards for Rare Mutations**

**PSMB5/CRBN too rare to validate (n=2-3). Mission requires RR thresholds.**

**Answer: D (mark as LITERATURE_ONLY) + B (aggregate to pathway)**

**Same as Question 5.** This is the same question rephrased.

**Key insight:** Rare mutations cannot be statistically validated. Accept this limitation and be transparent.

---

### **Question 22: Evo2 Integration Status**

**Mission says "working" but audit shows "not integrated".**

**Answer: B (trust audit - not integrated, needs implementation) + D (defer to P2)**

**Reasoning:**
1. **Audit is correct:** Evo2 service exists (`/api/evo/score_variant_multi`) but `predict_mm_resistance()` does NOT call it.

2. **Evo2 integration is P2:**
   - Requires: Call Evo2 API for each variant, compute delta scores
   - Mission target: r ≥ 0.3 correlation with response
   - Not critical for launch (Proxy SAE works without Evo2)

3. **Action:** Mark as P2, implement after TRUE SAE validation.

---

### **Question 23: Cytogenetics Validation Gap**

**Cytogenetics data missing from MMRF GDC extract. del(17p) HR=2.5 is literature-based.**

**Answer: A (accept literature-based HR values)**

**Reasoning:**
1. **del(17p) HR=2.5 is IMWG guideline.** This is the gold standard for MM risk stratification. Not controversial.

2. **FISH data is rarely in genomic databases.** Cytogenetics is clinical lab data, not WGS data.

3. **Action:**
   - Mark cytogenetics as "LITERATURE_ONLY"
   - HR values are clinically accepted
   - No need to delay for FISH validation

---

### **Question 24: Treatment Line Context**

**Mission says both "implemented" and "missing".**

**Answer: D (mission document has inconsistent status)**

**Reasoning:**
1. **Audit shows:** Treatment line logic exists in `predict_mm_resistance()` (multipliers for 1L/2L/3L).

2. **Mission inconsistency:** Line 34 says "missing" but line 57 says "implemented."

3. **Action:** Update mission document. Treatment line IS implemented.

---

### **Question 25: Frontend Integration Status**

**Frontend "exists" but resistance panel is missing.**

**Answer: D (needs enhancement - add resistance panel to existing page)**

**Reasoning:**
1. **Audit shows:**
   - `MyelomaDigitalTwin.jsx` ✅ EXISTS
   - `MyelomaResponseDisplay.jsx` ✅ EXISTS
   - `MMResistancePanel.jsx` ❌ NOT FOUND

2. **Resolution:** Create `MMResistancePanel.jsx` and add to `MyelomaDigitalTwin.jsx`.

3. **This is Work Item 5** (1-2 days effort).

---

### **Question 26: Pathway Service vs. TRUE SAE Diamond Features**

**Same as Question 7.** See answer above.

---

### **Question 27: Validation Test 4: RAS/MAPK Signal**

**KRAS shows NO SIGNAL (RR=0.93, p=0.87). Mission requires RAS/MAPK → Treatment Line Impact.**

**Answer: D (accept NO SIGNAL as valid result)**

**Reasoning:**
1. **Negative validation is still validation.** We tested RAS/MAPK on MMRF cohort. Result: no signal for mortality.

2. **This is biologically plausible:** RAS/MAPK in MM is not as strong a driver as in solid tumors.

3. **Action:**
   - Document: "RAS/MAPK: No mortality signal detected (RR=0.93, p=0.87)"
   - Consider: Re-analyze by treatment line (maybe signal exists in later lines)
   - Consider: Use different endpoint (PFS, response rate)

4. **Not a launch blocker.** Focus on validated signals (DIS3, TP53).

---

### **Question 28: MMRF GDC vs. Full CoMMpass**

**Same as Question 19.** See answer above.

---

### **Question 29: Mission Document Accuracy**

**Mission has multiple contradictions.**

**Answer: A (trust code audit) + B (update mission document)**

**Action for Plumber:**
1. **Create `MISSION_STATUS_RECONCILIATION.md`** that maps:
   - Mission claim → Audit finding → Action needed
   
2. **Update mission document** with accurate status.

3. **Source of truth:** Code audit > Mission document

---

### **Question 30: Clinical Validation Requirements**

**MMRF GDC is research cohort. Mission requires "clinical validation."**

**Answer: D (define "clinical validation" as "validated on clinical cohort")**

**Reasoning:**
1. **MMRF CoMMpass IS clinical data:** Real patients, real treatments, real outcomes. This is not synthetic or simulated.

2. **"Clinical validation" ≠ "Clinical trial."** Clinical trial validation is FDA-level requirement. For research tool, cohort validation is sufficient.

3. **Our validation is appropriate:**
   - Real patient mutations (from WGS)
   - Real outcomes (vital status, treatment response)
   - Published cohort (MMRF CoMMpass is peer-reviewed)

4. **Mark as "RUO" (Research Use Only)** if regulatory concerns. Clinical trial validation is future work.

---

## 📊 SUMMARY: KEY DECISIONS FOR PLUMBER

### **Launch Strategy:**
1. **V1 (Launch):** Proxy SAE (validated) + PSMB5/CRBN (literature-based) + Pathway Service + Frontend
2. **V2 (Enhancement):** TRUE SAE (when validated on MM) + Evo2 integration

### **Acceptable Failures:**
- Tests 1-3: LITERATURE_ONLY acceptable
- Test 5: TRUE SAE not required for launch
- Test 4 (DIS3+TP53): MUST PASS

### **Source of Truth:**
- Code audit > Mission document
- Update mission document to reflect reality

### **Risk Mitigation:**
- Modal failure → Proxy SAE fallback
- MMRF access failure → Use existing MMRF GDC (N=219)
- TRUE SAE validation failure → Ship Proxy SAE only

---

**Last Updated:** January 29, 2025  
**Author:** Zo  
**For:** Plumber (implementation guidance)
