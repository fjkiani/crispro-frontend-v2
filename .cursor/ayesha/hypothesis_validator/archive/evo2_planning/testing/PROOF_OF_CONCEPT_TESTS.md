# 🧪 PROOF OF CONCEPT: 10 CRITICAL TEST CASES

**Mission:** Validate the Evo2 Food Validator approach BEFORE building all 6 phases.

**Strategy:** Test core hypotheses with minimal implementation to confirm feasibility.

**Audience:** Lower-level agents executing this test plan

---

## **⚔️ TEST PHILOSOPHY**

### **What We're Testing:**
1. ✅ Can Evo2 variant scoring proxy work for compounds? (CRITICAL)
2. ✅ Can we fetch gene sequences reliably?
3. ✅ Can we extract targets dynamically?
4. ✅ Do P/E/SAE scores differentiate compounds?
5. ✅ Can we integrate with existing S/P/E framework?

### **What We're NOT Testing:**
- ❌ Full UI implementation
- ❌ Production error handling
- ❌ Complete SAE features
- ❌ Full treatment line logic

---

## **🎯 CRITICAL GO/NO-GO DECISIONS**

### **Decision Gate 1: After Test 1+2 (Infrastructure)**
- ✅ **GO:** Evo2 API + Ensembl work → Proceed to Test 4+6
- ❌ **NO-GO:** API format unknown or broken → Fix infrastructure first

### **Decision Gate 2: After Test 6 (Primary Evo2 Path)**
- ✅ **GO:** Variant proxy yields meaningful deltas → Use S in MVP
- ⚠️ **PARTIAL:** Variant proxy inconclusive → Skip S (neutral 0.5), rely on P/E/SAE
- ❌ **NO-GO:** Everything fails → Reconsider approach

**⚠️ STRATEGIC QUESTION (Pending Manager Review):**
- Should agents run Test 6 BEFORE building anything (test-first approach)?
- OR build P/E/SAE MVP first, then add Evo2 toggle later (build-first approach)?
- See `EXECUTION_DECISIONS.md` Q10 for full details.

### **Decision Gate 3: After Test 8+9 (Integration)**
- ✅ **GO:** End-to-end + multi-compound work → Proceed with full build
- ❌ **NO-GO:** Integration broken → Debug components

### **Decision Gate 4: After Test 10 (Validation)**
- ✅ **GO:** Scores correlate with biology → Approach validated
- ⚠️ **PARTIAL:** Mixed results → Tune thresholds/weights
- ❌ **NO-GO:** Random scores → Fundamental rethinking needed

**⚠️ STRATEGIC QUESTION (Pending Manager Review):**
- What exactly does "Evo2 works" mean?
  - Technical success (deltas > 0.2, stable)?
  - Biological correlation (Test 10 passes)?
  - Business value (meaningful differentiation)?
- See `EXECUTION_DECISIONS.md` Q8 for full details.

---

## **📋 TEST CASE 1: EVO2 API FORMAT DISCOVERY**

### **Purpose:**
Understand what `/api/evo/score` returns so we can parse it correctly.

### **Execution:**
```bash
cd .cursor/ayesha/hypothesis_validator/evo2_food_validator_doctrine/testing
python3 test_01_evo2_api_format.py
```

### **What This Tests:**
- Endpoint connectivity
- Response format (score vs likelihood)
- Score range (critical for normalization)
- Response time

### **Success → Proceed to Test 2**
### **Failure → BLOCKING - Fix Evo2 connection**

---

## **📋 TEST CASE 2: ENSEMBL GENE SEQUENCE FETCHING**

### **Purpose:**
Verify we can fetch gene sequences for all target genes.

### **Execution:**
```bash
python3 test_02_ensembl_gene_fetch.py
```

### **What This Tests:**
- Ensembl REST API connectivity
- Gene ID lookup (VDR, TP53, NFKB1)
- Sequence retrieval
- Response time per gene

### **Success → Proceed to Test 4**
### **Failure → Need fallback (local FASTA or genomic_intel endpoint)**

---

## **📋 TEST CASE 4: KNOWLEDGE BASE TARGET EXTRACTION**

### **Purpose:**
Verify `food_targets.json` exists and contains valid targets.

### **Execution:**
```bash
# Quick check
cat ../../../../data/food_targets.json | jq '.compounds[] | {compound: .compound, targets: .targets}'

# OR run test script
python3 test_04_knowledge_base.py
```

### **What This Tests:**
- File exists and is valid JSON
- Vitamin D has targets (VDR, TP53, etc.)
- At least 3 compounds have targets
- Target names are valid gene symbols

### **Success → Proceed to Test 6**
### **Failure → Create food_targets.json first**

---

## **📋 TEST CASE 6 (PRIMARY): VARIANT SCORING PROXY** ⭐

### **Purpose:**
**CRITICAL TEST:** Can we use variant scoring as a proxy for compound modulation?

### **Approach:**
Score a synthetic variant at a target gene's promoter/functional region as a proxy for compound effect.

### **Execution:**
```bash
python3 test_06_variant_scoring_proxy.py
```

### **What This Tests:**
1. **VDR Promoter Variant (Vitamin D proxy):**
   - Create variant in VDR promoter region (TSS - 500bp)
   - Score: Does variant impact expression?
   - Interpretation: Delta represents "how much VDR can be modulated"

2. **NFKB1 Promoter Variant (Curcumin proxy):**
   - Create variant in NFKB1 promoter region
   - Score: Does variant impact NF-κB pathway?
   - Interpretation: Delta represents "modulation potential"

3. **Comparison:**
   - High-impact variants → HIGH plausibility
   - Low-impact variants → LOW plausibility

### **Test Implementation:**
```python
# Example: VDR promoter variant as Vitamin D proxy

# Step 1: Get VDR TSS (transcription start site) from Ensembl
vdr_gene_info = fetch_ensembl_gene_info("VDR")
vdr_tss = vdr_gene_info['start']  # GRCh38 coordinate

# Step 2: Create synthetic variant at promoter (TSS - 500bp)
promoter_variant = {
    "assembly": "GRCh38",
    "chrom": "12",
    "pos": vdr_tss - 500,  # Upstream promoter
    "ref": "A",  # Get from reference
    "alt": "G",  # Synthetic mutation
    "model_id": "evo2_1b"
}

# Step 3: Score variant
response = await client.post(
    "http://127.0.0.1:8000/api/evo/score_variant_multi",
    json=promoter_variant
)

# Step 4: Extract delta (min_delta or other score)
delta = abs(response.json().get('min_delta', 0))

# Step 5: Classify
plausibility = "HIGH" if delta > 0.5 else "MODERATE" if delta > 0.2 else "LOW"

print(f"VDR promoter variant delta: {delta} ({plausibility})")
```

### **Success Criteria:**
- [ ] Variant scoring returns non-zero delta for promoter variants
- [ ] Different genes yield different deltas (not all same)
- [ ] Deltas are stable (repeat test shows similar value)
- [ ] Promoter variants show higher impact than random sites

### **GO/NO-GO:**
- ✅ **GO:** Delta > 0.2 for promoter variants, differentiation exists → Use this for S
- ⚠️ **PARTIAL:** Deltas too small or noisy → Set S=neutral (0.5), rely on P/E/SAE
- ❌ **NO-GO:** All deltas = 0 or random → Approach doesn't work

---

## **📋 TEST CASE 8: END-TO-END INTEGRATION (VITAMIN D)**

### **Purpose:**
Run complete workflow for one compound (Vitamin D).

### **Execution:**
```bash
python3 test_08_end_to_end_minimal.py
```

### **What This Tests:**
```
[1] Target Extraction
    → Input: "Vitamin D"
    → Output: ["VDR", "TP53"]
    → Source: food_targets.json

[2] Sequence Fetching
    → Fetch VDR sequence (Ensembl)
    → Fetch TP53 sequence (Ensembl)
    → Cache results

[3] Evo2 Scoring (if Test 6 passed)
    → Score VDR promoter variant (proxy for Vitamin D modulation)
    → Score TP53 promoter variant
    → Aggregate deltas → S = 0.45

[4] Pathway Alignment (P)
    → VDR targets → DNA repair pathway
    → TP53 → DNA repair pathway
    → Disease context: HRD+ (DNA repair deficient)
    → Alignment: HIGH (0.85)
    → P = 0.85

[5] Evidence Grade (E)
    → LLM search: "Vitamin D AND ovarian cancer"
    → Found: 15 papers
    → Grade: MODERATE
    → E = 0.60

[6] SAE Features
    → Input: HRD+, L3 post-platinum
    → Vitamin D rule: HIGH for HRD+
    → line_appropriateness = 0.9
    → cross_resistance = 0.0
    → sequencing_fitness = 0.85

[7] S/P/E Aggregation
    → overall_score = 0.4×0.45 + 0.3×0.85 + 0.3×0.60
    → overall_score = 0.18 + 0.255 + 0.18 = 0.615

[8] Confidence Modulation
    → base = (0.45+0.85+0.60)/3 = 0.63
    → evo2_boost = +0.05 (no HIGH targets)
    → sae_boost = (0.9+0.85)×0.05 = +0.0875
    → biomarker_boost = +0.05 (HRD+ → DNA repair match)
    → final = min(0.63+0.05+0.0875+0.05, 0.95) = 0.8175

[9] Verdict
    → score ≥ 0.65? NO (0.615)
    → score ≥ 0.45 AND confidence ≥ 0.50? YES (0.615, 0.82)
    → Verdict: WEAK_SUPPORT

[10] Recommendation
    → Dosage: "400-2000 IU daily"
    → Bioavailability: "GOOD"
    → Timing: "Best as adjunct to L3+ therapy"
    → Safety: "GOOD - minimal interactions"
```

### **Expected Output:**
```json
{
  "status": "SUCCESS",
  "compound": "Vitamin D",
  "verdict": "WEAK_SUPPORT",
  "overall_score": 0.615,
  "confidence": 0.82,
  "spe_breakdown": {
    "sequence": 0.45,
    "pathway": 0.85,
    "evidence": 0.60
  },
  "sae_features": {
    "line_appropriateness": 0.9,
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.85
  }
}
```

### **Success Criteria:**
- [ ] All 10 steps execute without errors
- [ ] Overall score is computed (0-1 range)
- [ ] Verdict is classified correctly
- [ ] Confidence modulation applied correctly
- [ ] Total execution time < 60s

### **GO/NO-GO:**
- ✅ **GO:** End-to-end works → Proceed with full build (Phases 1-6)
- ❌ **NO-GO:** Any step fails → Debug that component

---

## **📋 TEST CASE 9: MULTI-COMPOUND DIFFERENTIATION**

### **Purpose:**
Verify that different compounds get different scores (not all identical).

### **Execution:**
```bash
python3 test_09_multi_compound.py
```

### **What This Tests:**
Run 5 compounds through complete workflow:

| Compound | Targets | Expected P/E/SAE Behavior |
|----------|---------|---------------------------|
| **Vitamin D** | VDR, TP53 | P=HIGH (DNA repair), E=MODERATE, SAE=HIGH (HRD+) |
| **Curcumin** | NFKB1, PTGS2 | P=MODERATE (inflammation), E=MODERATE, SAE=MODERATE |
| **NAC** | Nrf2, GSR | P=MODERATE (oxidative stress), E=WEAK, SAE=HIGH (post-platinum) |
| **Omega-3** | PPARA, NFKB1 | P=MODERATE (inflammation), E=MODERATE, SAE=HIGH (post-chemo) |
| **Green Tea** | EGCG, TP53 | P=LOW (weak pathway match), E=WEAK, SAE=LOW |

### **Expected Ranking (by overall_score):**
1. NAC (HIGH SAE for post-platinum)
2. Vitamin D (HIGH P+SAE for HRD+)
3. Omega-3 (HIGH SAE for post-chemo)
4. Curcumin (MODERATE across board)
5. Green Tea (LOW - weak evidence)

### **Success Criteria:**
- [ ] All 5 compounds process successfully
- [ ] Scores are different (not all identical)
- [ ] Ranking makes biological sense
- [ ] NAC scores highest for post-platinum context
- [ ] Green Tea scores lowest (weakest evidence)

### **GO/NO-GO:**
- ✅ **GO:** Compounds differentiate based on P/E/SAE → Approach works
- ⚠️ **PARTIAL:** All scores similar but S=neutral by design → Expected (P/E/SAE still provide some differentiation)
- ❌ **NO-GO:** All scores identical despite different P/E/SAE → Logic broken

---

## **📋 TEST CASE 10: BIOLOGICAL SANITY CHECK**

### **Purpose:**
Verify scores correlate with known biological mechanisms (not random).

### **Execution:**
```bash
python3 test_10_sanity_check.py
```

### **What This Tests:**

**Test A: NAC Post-Platinum (KNOWN GOOD)**
- Context: Ovarian cancer, HRD+, post-carboplatin (oxidative stress)
- Expected: HIGH/MODERATE (SAE should be HIGH for post-platinum)
- Mechanism: NAC reduces oxidative stress from platinum

**Test B: Random Placebo (KNOWN BAD)**
- Context: Ovarian cancer, HRD+
- Compound: "Supplement X" (unknown targets)
- Expected: LOW (no targets, no pathways, no evidence)

**Test C: Curcumin for Inflammation (KNOWN MODERATE)**
- Context: Chronic inflammation, NF-κB active
- Expected: MODERATE (P=HIGH for inflammation, E=MODERATE, but bioavailability poor)

### **Success Criteria:**
- [ ] NAC scores HIGH/MODERATE for post-platinum
- [ ] Random placebo scores LOW
- [ ] Curcumin scores MODERATE for inflammation context
- [ ] Scores align with biological literature
- [ ] False positives minimized (bad compounds don't score HIGH)

### **GO/NO-GO:**
- ✅ **GO:** Scores correlate with biology → Approach validated
- ❌ **NO-GO:** Scores are random or inverse of biology → Fundamental problem

---

## **⚔️ SIMPLIFIED EXECUTION PLAN FOR AGENTS**

### **Phase 1: Infrastructure (30 min)**
Run these tests to confirm infrastructure works:

```bash
cd .cursor/ayesha/hypothesis_validator/evo2_food_validator_doctrine/testing

# Test 1: Evo2 API format
python3 test_01_evo2_api_format.py
# Expected: ✅ Gets response, identifies score format

# Test 2: Ensembl gene fetching
python3 test_02_ensembl_gene_fetch.py
# Expected: ✅ Fetches VDR, TP53, NFKB1 sequences

# Test 4: Knowledge base
cat ../../../../data/food_targets.json | jq '.'
# Expected: ✅ File exists, has Vitamin D/Curcumin/NAC/etc.
```

**Decision Gate 1:**
- ✅ All 3 pass → Proceed to Phase 2
- ❌ Any fail → Fix infrastructure, do not proceed

---

### **Phase 2: Primary Evo2 Path (30 min)**

**CRITICAL:** Test if variant scoring proxy works for compounds.

```bash
# Test 6: Variant scoring proxy (PRIMARY APPROACH)
python3 test_06_variant_scoring_proxy.py
# This will:
# 1. Fetch VDR gene info (TSS coordinate)
# 2. Create synthetic variant at VDR promoter (TSS - 500bp)
# 3. Score variant with /api/evo/score_variant_multi
# 4. Extract delta
# 5. Check if delta > 0.2 (meaningful impact)
```

**Decision Gate 2:**
- ✅ **Delta > 0.2:** Variant proxy works → Use S in MVP
- ⚠️ **Delta < 0.2 or noisy:** Skip S (set to neutral 0.5) → Rely on P/E/SAE only
- ❌ **Delta = 0 or errors:** Approach broken → Stop

---

### **Phase 3: Integration (1 hour)**

**IF Test 6 passed OR we're proceeding with neutral S:**

```bash
# Test 8: End-to-end minimal (Vitamin D)
python3 test_08_end_to_end_minimal.py
# This will:
# 1. Extract targets (Vitamin D → VDR, TP53)
# 2. Fetch sequences
# 3. Score S (variant proxy OR neutral if skipped)
# 4. Compute P (pathway alignment)
# 5. Compute E (literature grade - mocked for now)
# 6. Compute SAE (supplement rules)
# 7. Aggregate S/P/E → overall_score
# 8. Modulate confidence
# 9. Classify verdict
# Expected: overall_score ≈ 0.6-0.7, confidence ≈ 0.75-0.85, verdict = WEAK_SUPPORT

# Test 9: Multi-compound
python3 test_09_multi_compound.py
# This will:
# 1. Run Test 8 for 5 compounds
# 2. Rank by overall_score
# 3. Verify differentiation exists
# Expected: NAC > Vitamin D > Omega-3 > Curcumin > Green Tea
```

**Decision Gate 3:**
- ✅ Tests pass → Proceed to Phase 4 (validation)
- ❌ Any fail → Debug specific component

---

### **Phase 4: Validation (30 min)**

**Confirm scores make biological sense:**

```bash
# Test 10: Sanity check
python3 test_10_sanity_check.py
# This will:
# 1. Test NAC post-platinum (should score HIGH/MODERATE)
# 2. Test random placebo (should score LOW)
# 3. Test Curcumin for inflammation (should score MODERATE)
# 4. Validate scores align with biology
```

**Decision Gate 4:**
- ✅ Scores correlate with biology → **APPROACH VALIDATED - PROCEED WITH FULL BUILD**
- ❌ Scores random → **NO-GO - RECONSIDER APPROACH**

---

## **📊 MINIMUM VIABLE SUCCESS CRITERIA**

To proceed with full build, we need:

### **MUST HAVE:**
- ✅ Test 1: Evo2 API format known
- ✅ Test 2: Gene sequences fetchable
- ✅ Test 4: Knowledge base has targets
- ✅ Test 6: Variant proxy works OR we accept neutral S
- ✅ Test 8: End-to-end integration works
- ✅ Test 9: Compounds differentiate (via P/E/SAE at minimum)
- ✅ Test 10: Scores correlate with biology

### **NICE TO HAVE:**
- ⚠️ Test 5: ChEMBL extraction (can add later)
- ⚠️ Test 3: Baseline vs intervention (documentation only)
- ⚠️ Test 7: Generation approach (exploratory)

---

## **🎯 FILES PROVIDED FOR AGENTS**

**Test Scripts (ready to run):**
1. ✅ `test_01_evo2_api_format.py` - Discovers Evo2 response format
2. ✅ `test_02_ensembl_gene_fetch.py` - Tests gene sequence fetching
3. ✅ `test_03_baseline_vs_intervention.py` - Documents why sequence-only fails
4. ⚠️ `test_04_knowledge_base.py` - TO BE CREATED
5. ⚠️ `test_06_variant_scoring_proxy.py` - TO BE CREATED (CRITICAL)
6. ⚠️ `test_07_generate_alternative.py` - TO BE CREATED (optional)
7. ⚠️ `test_08_end_to_end_minimal.py` - TO BE CREATED
8. ⚠️ `test_09_multi_compound.py` - TO BE CREATED
9. ⚠️ `test_10_sanity_check.py` - TO BE CREATED

**Runner Script:**
- ✅ `run_all_tests.sh` - Executes all tests in order

---

## **⚔️ WHAT TO BUILD NEXT (FOR AGENTS)**

### **Immediate (TONIGHT):**
1. Create missing test files (test_04, test_06, test_08, test_09, test_10)
2. Run Phase 1 tests (1, 2, 4) - Infrastructure validation
3. Run Phase 2 tests (6) - Primary Evo2 path validation
4. **GO/NO-GO decision** based on Test 6 results

### **If Test 6 Passes (Variant Proxy Works):**
- Build full MVP with S/P/E/SAE (all 6 phases)
- Total time: 7-8 hours

### **If Test 6 Fails (Variant Proxy Doesn't Work):**
- Build P/E/SAE-only MVP (skip Phase 1, neutral S)
- Total time: 4-5 hours
- Add Evo2 later as experimental feature

**⚠️ STRATEGIC QUESTION (Pending Manager Review):**
- Manager recommends building P/E/SAE first (guaranteed 90% confidence), then adding Evo2 toggle (60% confidence, questionable but novel).
- Should agents follow this conservative approach, or stick with test-first strategy?
- See `EXECUTION_DECISIONS.md` Q7, Q10, Q13 for full details.

---

## **📝 AGENT INSTRUCTIONS**

**Step 1:** Run infrastructure tests (1, 2, 4)
**Step 2:** Create and run Test 6 (variant proxy - CRITICAL)
**Step 3:** Based on Test 6 results, create Test 8+9+10
**Step 4:** Run integration/validation tests
**Step 5:** Report GO/NO-GO decision with evidence

**DO NOT:**
- ❌ Build all 6 phases before testing
- ❌ Assume variant proxy will work without testing
- ❌ Proceed to full build if tests fail

**DO:**
- ✅ Test infrastructure first (30 min)
- ✅ Test primary approach (Test 6 - 30 min)
- ✅ Make GO/NO-GO decision based on results
- ✅ Build MVP (P/E/SAE only if Evo2 fails)

---

**Agent Checklist:**
- [ ] Created all missing test files
- [ ] Ran Phase 1 tests (infrastructure)
- [ ] Ran Phase 2 tests (variant proxy)
- [ ] Made GO/NO-GO decision
- [ ] Reported results to commander

**Status: ⚔️ READY FOR AGENT EXECUTION** 🎯
