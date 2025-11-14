# 🎯 EVO2 REAL INTEGRATION PLAN - NO MORE MOCKS

**Date:** November 6, 2025  
**Purpose:** Explain REAL Evo2 integration based on lessons from shark cartilage failures  
**Commander's Question:** "What are we sending to Evo2? How are we preventing junk DNA?"

---

## 🚨 **CRITICAL LESSONS FROM SHARK CARTILAGE FAILURES**

### **Failure #1: Evo2's Dual Nature**
- **As Scorer (Sniper)**: ✅ Precision variant effect prediction
- **As Generator (Poet)**: ⚠️ Needs rich context, NOT short prompts

### **Failure #2: The "Pathological Attractor"**
- **Problem**: Low-complexity, repetitive sequences force Evo2 into failure states
- **Result**: Generates JUNK DNA (poly-A, poly-AT, meaningless repeats)
- **Solution**: ALWAYS use biologically rich context

### **Failure #3: Wrong Language**
- **Problem**: Tried using English prompts
- **Reality**: Evo2 trained on DNA language (ATCG), not English
- **Solution**: Use DNA sequences as input

---

## ✅ **HOW EVO2 ACTUALLY WORKS (SCORING)**

### **For Variant Impact Prediction (What We SHOULD Use):**

**Input Format:**
```python
{
    "sequence": "ATGCGATCG...",  # Reference sequence (DNA)
    "model_id": "evo2_1b",
    "window_size": 4096,  # Context window
    "variants": [
        {
            "pos": 100,
            "ref": "A",
            "alt": "T"
        }
    ]
}
```

**What Evo2 Does:**
1. Scores reference sequence → `likelihood_ref`
2. Scores variant sequence → `likelihood_alt`
3. Returns `delta = likelihood_alt - likelihood_ref`

**Interpretation:**
- `delta < 0`: Variant **decreases** sequence likelihood (likely disruptive)
- `delta > 0`: Variant **increases** sequence likelihood (likely stabilizing)
- `|delta|` large: High confidence prediction

---

## ❌ **WHAT WE'RE CURRENTLY DOING (WRONG)**

**File:** `therapeutic_optimizer.py` Line 388

```python
# 4. Evo2 delta score (STUB - would call actual Evo2)
# TODO: Implement real Evo2 delta scoring
evo2_score = 0.75  # Mock score  <-- HARDCODED! NO API CALL!
```

**Problems:**
1. ❌ Not calling Evo2 API at all
2. ❌ Can't distinguish known vs random sequences
3. ❌ Can't prevent junk DNA generation
4. ❌ Ignoring `target_sequence` parameter
5. ❌ No biological context used

---

## ✅ **CORRECT IMPLEMENTATION - REAL EVO2 CALLS**

### **Step 1: Score Candidate Sequence**

```python
async def _score_candidate_REAL(
    self,
    sequence: str,
    target_gene: str,
    target_sequence: str = None,
    therapeutic_type: str = "protein"
) -> Tuple[float, Dict[str, any]]:
    """
    Score therapeutic candidate using REAL Evo2 API.
    
    Prevents junk DNA by:
    1. Safety validation FIRST (viral, GC, homopolymer)
    2. Biological context (target sequence)
    3. Real Evo2 delta scoring
    """
    
    metrics = {}
    
    # STEP 1: Safety validation (prevent junk DNA)
    safety_result = self.safety_validator.validate_sequence(sequence)
    if not safety_result.is_safe:
        return 0.0, {
            "error": "Safety blocked",
            "reason": safety_result.reason,
            "evo2_delta": 0.0
        }
    
    # STEP 2: Build biological context
    # For protein therapeutics, need DNA context
    if therapeutic_type == "protein":
        # Translate protein → DNA (use standard codon table)
        from Bio.Seq import Seq
        dna_sequence = str(Seq(sequence).back_translate())
    else:
        dna_sequence = sequence
    
    # STEP 3: Call REAL Evo2 API
    try:
        # Get Evo2 proxy service
        evo2_client = await self._get_evo2_client()
        
        # Score reference (target sequence)
        if target_sequence:
            ref_result = await evo2_client.score_sequence(
                sequence=target_sequence,
                model_id="evo2_1b"  # CRITICAL: Use 1B for cost
            )
            ref_likelihood = ref_result["likelihood"]
        else:
            ref_likelihood = 0.0  # No baseline
        
        # Score candidate (therapeutic sequence)
        alt_result = await evo2_client.score_sequence(
            sequence=dna_sequence,
            model_id="evo2_1b"  # CRITICAL: Use 1B for cost
        )
        alt_likelihood = alt_result["likelihood"]
        
        # Calculate delta
        evo2_delta = alt_likelihood - ref_likelihood
        
        metrics["evo2_delta"] = evo2_delta
        metrics["ref_likelihood"] = ref_likelihood
        metrics["alt_likelihood"] = alt_likelihood
        metrics["api_called"] = "REAL"  # Prove it's not mocked
        
    except Exception as e:
        self.logger.error(f"Evo2 API call failed: {e}")
        # Graceful degradation
        evo2_delta = 0.0
        metrics["evo2_delta"] = 0.0
        metrics["error"] = str(e)
        metrics["api_called"] = "FAILED"
    
    # STEP 4: Other heuristics (GC, homopolymer)
    gc_content = self._calculate_gc_content(sequence)
    metrics["gc_content"] = gc_content
    
    max_homopolymer = self._find_max_homopolymer(sequence)
    metrics["max_homopolymer"] = max_homopolymer
    
    # STEP 5: Composite score
    # Weight Evo2 heavily (80%) since it's the biological signal
    gc_score = 1.0 if 0.40 <= gc_content <= 0.60 else 0.5
    homopoly_score = 1.0 if max_homopolymer <= 4 else 0.5
    
    # Normalize Evo2 delta to 0-1 range
    # Typical deltas: -10 to +10
    evo2_normalized = (evo2_delta + 10) / 20  # Maps [-10,10] → [0,1]
    evo2_normalized = max(0.0, min(1.0, evo2_normalized))  # Clamp
    
    composite_score = (
        0.80 * evo2_normalized +  # HEAVY weight on Evo2
        0.10 * gc_score +
        0.10 * homopoly_score
    )
    
    metrics["composite_score"] = composite_score
    
    return composite_score, metrics
```

---

## 🎯 **HOW THIS PREVENTS JUNK DNA**

### **1. Safety Validation First**
- Blocks poly-A, poly-AT, low-complexity BEFORE Evo2 call
- Checks viral content, GC extremes, homopolymers
- **Result**: Only biologically plausible sequences reach Evo2

### **2. Biological Context**
- Uses `target_sequence` to establish baseline
- Compares candidate vs target (delta scoring)
- **Result**: Evo2 has context to judge biological relevance

### **3. Real Evo2 Delta Scoring**
- Actually calls Evo2 API (not hardcoded 0.75)
- Returns REAL likelihood differences
- **Result**: Can distinguish known therapeutic vs junk

### **4. Cost Control**
- ALWAYS use `evo2_1b` (1B parameter, not 7B/40B)
- Cost: ~$0.0001 per sequence (~500 tokens)
- **Result**: 100 tests = $0.01 (acceptable)

---

## 🧪 **TEST VALIDATION STRATEGY**

### **Test 1: Known Anti-VEGF vs Random**
**Expected:**
- Known anti-VEGF: Evo2 delta ~ -5 to -10 (disrupts VEGF)
- Random sequence: Evo2 delta ~ 0 (no biological signal)
- **Validation**: Known should score LOWER delta (more disruptive)

### **Test 2: Vitamin D Agonist vs Antagonist**
**Expected:**
- Agonist: Evo2 delta ~ +3 to +5 (stabilizes VDR)
- Antagonist: Evo2 delta ~ -3 to -5 (disrupts VDR)
- **Validation**: Deltas should be OPPOSITE signs

### **Test 3: Curcumin Multi-Target**
**Expected:**
- NFKB1: Evo2 delta ~ -4
- PTGS2: Evo2 delta ~ -6
- AKT1: Evo2 delta ~ -3
- **Validation**: Deltas should VARY (not all 0.75!)

### **Test 4: Nonsense Control**
**Expected:**
- Poly-A: **BLOCKED by safety** (never reaches Evo2)
- Low-complexity: Evo2 delta ~ 0 (no signal)
- Biological: Evo2 delta ~ -5
- **Validation**: Junk scores LOWER than biological

---

## 📊 **ACCEPTANCE CRITERIA**

**Evo2 Integration is REAL when:**
1. ✅ API calls return DIFFERENT scores for different sequences
2. ✅ Known therapeutic scores BETTER than random
3. ✅ Agonist vs antagonist scores have OPPOSITE signs
4. ✅ Multi-target shows VARYING deltas (not hardcoded)
5. ✅ Junk DNA blocked OR scores near zero
6. ✅ Cost < $0.001 per test (1B model confirmed)
7. ✅ Provenance shows `"api_called": "REAL"`

**If ANY test shows all scores = 0.75 → STILL MOCKED!**

---

## 🚀 **NEXT STEPS**

1. **REPLACE** `therapeutic_optimizer.py` line 388 with REAL Evo2 call
2. **ADD** `_get_evo2_client()` method to get API client
3. **UPDATE** tests to check for `"api_called": "REAL"` in metrics
4. **RUN** all 4 test cases and verify DIFFERENT scores
5. **VALIDATE** cost is ~$0.0004 (confirms 1B model used)

---

## ⚔️ **COMMANDER'S QUESTIONS ANSWERED**

**Q: What are we sending to Evo2?**
**A:** DNA sequences (not English) with biological context (target gene sequence)

**Q: How does Evo2 work?**
**A:** Scores sequence likelihood, compares reference vs candidate, returns delta

**Q: How are we preventing junk DNA?**
**A:** 
1. Safety validation FIRST (blocks poly-A, low-complexity)
2. Biological context (target sequence establishes baseline)
3. Real delta scoring (junk scores near zero vs biological scores high)

**Q: How will we not generate junk vs actual sequence?**
**A:** 
- Generation is SEPARATE from scoring (we're only testing scoring now)
- When we DO generate, we'll use rich prompts (from shark cartilage lessons)
- Safety validator blocks junk before it reaches Evo2

---

**STATUS:** ⚔️ **READY TO IMPLEMENT REAL EVO2 INTEGRATION**





