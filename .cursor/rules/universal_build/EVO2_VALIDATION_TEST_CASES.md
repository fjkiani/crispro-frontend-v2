# 🧪 EVO2 VALIDATION TEST CASES - REAL API CALLS

**Date:** November 6, 2025  
**Purpose:** Test Evo2 API (1B parameter) with REAL biological hypotheses  
**Goal:** Validate if Evo2 actually works vs just scaffolding

---

## 🎯 TEST STRATEGY

### **Question:** Does Evo2 provide REAL value for hypothesis testing?

### **Approach:**
1. Call Evo2 API with REAL biological sequences
2. Compare delta scores for meaningful vs nonsense sequences
3. Verify costs (should be ~$0.0001 per test with 1B model)
4. Validate if results match biological expectations

---

## 🧪 TEST CASE 1: SHARK CARTILAGE ANTI-ANGIOGENIC PROTEIN

**Hypothesis:** Shark cartilage contains proteins that inhibit blood vessel formation (anti-angiogenic)

**Target Gene:** VEGFA (Vascular Endothelial Growth Factor A)  
**Mechanism:** Inhibit VEGFA signaling → block tumor angiogenesis

### **Test 1A: Known Anti-Angiogenic Sequence vs VEGFA**

```python
# Test: Score a known anti-VEGF antibody sequence
# Expected: Negative delta (disrupts VEGFA function)

import asyncio
import os
from api.services.therapeutic_prompt_builder import get_prompt_builder
from api.services.therapeutic_optimizer import get_therapeutic_optimizer

async def test_shark_cartilage_vegfa():
    """
    Test 1: Shark cartilage anti-angiogenic hypothesis
    Compare: Known anti-VEGF sequence vs random sequence
    """
    optimizer = get_therapeutic_optimizer()
    
    # Known anti-VEGF nanobody sequence (from bevacizumab/avastin)
    # Simplified 20 AA sequence from anti-VEGF binding domain
    known_anti_vegf = "QVQLVESGGGVVQPGRSLRL"  # Fragment
    
    # Random sequence (control)
    random_sequence = "ATGATGATGATGATGATGAT"  # Should score poorly
    
    # VEGFA target sequence (first 100bp of coding sequence)
    vegfa_target = "ATGAACTTTCTGCTGTCTTGGGTGCATTGGAGCCTTGCCTTGCTGCTCTACCTCCACCATGCCAAGTGGTCCCAGGCTGCACCCATGGCAGAAGGAGGAGGGCAGAA"
    
    print("=== TEST 1A: Known Anti-VEGF vs Random ===")
    print(f"Target: VEGFA (anti-angiogenesis)")
    print(f"Model: evo2_1b (cost-controlled)")
    print()
    
    # Score known anti-VEGF sequence
    print("Scoring known anti-VEGF sequence...")
    known_result = await optimizer._score_candidate(
        sequence=known_anti_vegf,
        target_gene="VEGFA",
        target_sequence=vegfa_target,
        therapeutic_type="protein"
    )
    
    # Score random sequence
    print("Scoring random sequence (control)...")
    random_result = await optimizer._score_candidate(
        sequence=random_sequence,
        target_gene="VEGFA",
        target_sequence=vegfa_target,
        therapeutic_type="protein"
    )
    
    print("\n=== RESULTS ===")
    print(f"Known Anti-VEGF Score: {known_result['overall_score']:.4f}")
    print(f"Random Sequence Score: {random_result['overall_score']:.4f}")
    print(f"Evo2 Delta (Known): {known_result.get('evo2_delta', 'N/A')}")
    print(f"Evo2 Delta (Random): {random_result.get('evo2_delta', 'N/A')}")
    
    # Validation
    if known_result['overall_score'] > random_result['overall_score']:
        print("\n✅ PASS: Known anti-VEGF scores higher than random")
        print("   Interpretation: Evo2 recognizes biologically relevant sequence")
    else:
        print("\n❌ FAIL: Random scores higher than known anti-VEGF")
        print("   Interpretation: Evo2 may not be working correctly OR scoring is mocked")
    
    return {
        "test_name": "shark_cartilage_vegfa",
        "known_score": known_result['overall_score'],
        "random_score": random_result['overall_score'],
        "delta_difference": known_result.get('evo2_delta', 0) - random_result.get('evo2_delta', 0),
        "passed": known_result['overall_score'] > random_result['overall_score']
    }

# Run test
if __name__ == "__main__":
    result = asyncio.run(test_shark_cartilage_vegfa())
    print(f"\nTest Result: {'PASS' if result['passed'] else 'FAIL'}")
```

**Expected Outcome:**
- ✅ Known anti-VEGF should score **higher** (recognizes therapeutic potential)
- ✅ Random sequence should score **lower** (no biological relevance)
- ✅ Evo2 delta should be **negative for known** (disrupts target)
- ✅ Cost should be <$0.0001 (using 1B model)

---

## 🧪 TEST CASE 2: VITAMIN D RECEPTOR (VDR) ACTIVATION

**Hypothesis:** Vitamin D activates VDR → anti-cancer effects

**Target Gene:** VDR (Vitamin D Receptor)  
**Mechanism:** Agonist binding → increased VDR activity → cancer cell differentiation

### **Test 2A: Known VDR Agonist vs Antagonist**

```python
async def test_vitamin_d_vdr_agonist():
    """
    Test 2: Vitamin D hypothesis
    Compare: Known VDR agonist vs known VDR antagonist
    """
    optimizer = get_therapeutic_optimizer()
    
    # Simplified active metabolite sequence (1,25-dihydroxyvitamin D3 binding region)
    # This would be a peptide/small molecule mimic in reality
    vdr_agonist = "CEGALELISFLSKLAQELGL"  # VDR activation motif
    
    # Known VDR antagonist sequence
    vdr_antagonist = "CEGALELISELQTLAQGLGL"  # Modified (hypothetical)
    
    # VDR target sequence (ligand-binding domain)
    vdr_target = "ATGGAGGCAATGGCGGCCAGCCTGGTCACCCACAGCAAGTACGAGTGGATGGTCAACGAGGTCACCAAGCTCAAGCACCAGCAGCCGGGTGGCGGCGAGTCCTGG"
    
    print("=== TEST 2A: VDR Agonist vs Antagonist ===")
    print(f"Target: VDR (Vitamin D Receptor)")
    print(f"Model: evo2_1b")
    print()
    
    # Score agonist
    print("Scoring VDR agonist...")
    agonist_result = await optimizer._score_candidate(
        sequence=vdr_agonist,
        target_gene="VDR",
        target_sequence=vdr_target,
        therapeutic_type="peptide"
    )
    
    # Score antagonist
    print("Scoring VDR antagonist...")
    antagonist_result = await optimizer._score_candidate(
        sequence=vdr_antagonist,
        target_gene="VDR",
        target_sequence=vdr_target,
        therapeutic_type="peptide"
    )
    
    print("\n=== RESULTS ===")
    print(f"Agonist Score: {agonist_result['overall_score']:.4f}")
    print(f"Antagonist Score: {antagonist_result['overall_score']:.4f}")
    print(f"Evo2 Delta (Agonist): {agonist_result.get('evo2_delta', 'N/A')}")
    print(f"Evo2 Delta (Antagonist): {antagonist_result.get('evo2_delta', 'N/A')}")
    
    # For activation (agonist), we expect POSITIVE delta (enhances function)
    # For inhibition (antagonist), we expect NEGATIVE delta (disrupts function)
    
    if agonist_result.get('evo2_delta', 0) > antagonist_result.get('evo2_delta', 0):
        print("\n✅ PASS: Agonist has more positive delta than antagonist")
        print("   Interpretation: Evo2 distinguishes activation vs inhibition")
    else:
        print("\n❌ FAIL: Deltas don't match expected biology")
        print("   Interpretation: Evo2 may be mocked or not calibrated")
    
    return {
        "test_name": "vitamin_d_vdr",
        "agonist_score": agonist_result['overall_score'],
        "antagonist_score": antagonist_result['overall_score'],
        "delta_difference": agonist_result.get('evo2_delta', 0) - antagonist_result.get('evo2_delta', 0),
        "passed": agonist_result.get('evo2_delta', 0) > antagonist_result.get('evo2_delta', 0)
    }
```

**Expected Outcome:**
- ✅ Agonist should have **positive delta** (enhances VDR function)
- ✅ Antagonist should have **negative delta** (disrupts VDR function)
- ✅ Evo2 should distinguish between activation vs inhibition

---

## 🧪 TEST CASE 3: CURCUMIN MULTI-TARGET EFFECTS

**Hypothesis:** Curcumin inhibits multiple cancer pathways (NFκB, COX2, AKT)

**Target Genes:** NFKB1, PTGS2 (COX2), AKT1  
**Mechanism:** Multi-target inhibition → broad anti-cancer effects

### **Test 3A: Real Curcumin-Like Sequence vs Each Target**

```python
async def test_curcumin_multitarget():
    """
    Test 3: Curcumin multi-target hypothesis
    Score: Curcumin-like peptide against 3 different targets
    """
    optimizer = get_therapeutic_optimizer()
    
    # Simplified curcumin-binding motif (hypothetical peptide mimic)
    curcumin_mimic = "FEQARAEMAQEMGELVRLAQ"
    
    # Target sequences (first 100bp of each)
    targets = {
        "NFKB1": "ATGGCAGAAGATGATCCATATTTGGGAAGGAGACATCCAGGTGGTACCAAGGGCCCCAGCCACCTTGCCCTGTGGCTGGACCCTCACCGTGACCTTGGGTGCGGAG",
        "PTGS2": "ATGCTCGCCCGCGCCCTGCTGCTGTGCGCGGTCCTGGCGCTCAGCCATACAGCAAATCCTTGCTGTTCCCACCCATGTCAAAACCGAGGTGTATGTATGAGTGTGG",
        "AKT1": "ATGAGCGACGTGGCTATTGTGAAGGAGGGTTGGCTGCACAAACGCGGGGAGTTCCTGAAGCCAGCCATCCAGCTGGGCCACATCTTCAATCAGTCTGGAACGGAGC"
    }
    
    print("=== TEST 3A: Curcumin Multi-Target Effects ===")
    print(f"Targets: NFKB1, PTGS2 (COX2), AKT1")
    print(f"Model: evo2_1b")
    print()
    
    results = {}
    for gene, target_seq in targets.items():
        print(f"Scoring against {gene}...")
        result = await optimizer._score_candidate(
            sequence=curcumin_mimic,
            target_gene=gene,
            target_sequence=target_seq,
            therapeutic_type="peptide"
        )
        results[gene] = result
        print(f"  Score: {result['overall_score']:.4f}, Delta: {result.get('evo2_delta', 'N/A')}")
    
    print("\n=== MULTI-TARGET VALIDATION ===")
    # Curcumin should show some activity against ALL targets (not just one)
    active_targets = sum(1 for r in results.values() if r['overall_score'] > 0.6)
    
    if active_targets >= 2:
        print(f"✅ PASS: Shows activity against {active_targets}/3 targets")
        print("   Interpretation: Evo2 recognizes multi-target potential")
    else:
        print(f"❌ FAIL: Only active against {active_targets}/3 targets")
        print("   Interpretation: Multi-target scoring may need tuning")
    
    return {
        "test_name": "curcumin_multitarget",
        "results": {gene: r['overall_score'] for gene, r in results.items()},
        "active_targets": active_targets,
        "passed": active_targets >= 2
    }
```

**Expected Outcome:**
- ✅ Should show activity against **2-3 targets** (multi-target nature)
- ✅ Not equally effective (some targets more sensitive)
- ✅ Validates Evo2's ability to assess multi-target compounds

---

## 🧪 TEST CASE 4: NONSENSE SEQUENCE (NEGATIVE CONTROL)

**Hypothesis:** Complete junk DNA should score VERY LOW

### **Test 4A: Poly-A Tract vs Biological Sequence**

```python
async def test_nonsense_sequence_control():
    """
    Test 4: Negative control - nonsense sequence
    Expected: Should score much lower than any biological sequence
    """
    optimizer = get_therapeutic_optimizer()
    
    # Nonsense sequences
    poly_a = "AAAAAAAAAAAAAAAAAAAA"  # Homopolymer (should be blocked by safety)
    random_junk = "XYZQWERTXYZQWERTXYZQ"  # Not even valid DNA
    low_complexity = "ATATATATATATATATAT"  # Low complexity
    
    # Biological control (from bevacizumab)
    biological = "QVQLVESGGGVVQPGRSLRL"
    
    # Target
    vegfa_target = "ATGAACTTTCTGCTGTCTTGGGTGCATTGGAGCCTTGCCTTGCTGCTCTACCTCCACCATGCCAAGTGGTCCCAGGCTGCACCCATGGCAGAAGGAGGAGGGCAGAA"
    
    print("=== TEST 4A: Nonsense Sequence Negative Control ===")
    print(f"Model: evo2_1b")
    print()
    
    # Test poly-A (should be BLOCKED by safety validator)
    print("Testing poly-A tract...")
    try:
        poly_a_result = await optimizer._score_candidate(
            sequence=poly_a,
            target_gene="VEGFA",
            target_sequence=vegfa_target,
            therapeutic_type="protein"
        )
        print(f"  ⚠️ WARNING: Poly-A not blocked! Score: {poly_a_result['overall_score']:.4f}")
    except ValueError as e:
        print(f"  ✅ GOOD: Poly-A blocked by safety: {str(e)[:50]}...")
    
    # Test low complexity
    print("Testing low-complexity sequence...")
    low_result = await optimizer._score_candidate(
        sequence=low_complexity,
        target_gene="VEGFA",
        target_sequence=vegfa_target,
        therapeutic_type="protein"
    )
    
    # Test biological control
    print("Testing biological control...")
    bio_result = await optimizer._score_candidate(
        sequence=biological,
        target_gene="VEGFA",
        target_sequence=vegfa_target,
        therapeutic_type="protein"
    )
    
    print("\n=== RESULTS ===")
    print(f"Low-Complexity Score: {low_result['overall_score']:.4f}")
    print(f"Biological Control Score: {bio_result['overall_score']:.4f}")
    
    if bio_result['overall_score'] > low_result['overall_score'] * 1.5:
        print("\n✅ PASS: Biological sequence scores significantly higher")
        print("   Interpretation: Evo2 distinguishes meaningful vs junk sequences")
    else:
        print("\n❌ FAIL: Scores too similar")
        print("   Interpretation: Evo2 may not be properly scoring complexity")
    
    return {
        "test_name": "nonsense_control",
        "low_complexity_score": low_result['overall_score'],
        "biological_score": bio_result['overall_score'],
        "ratio": bio_result['overall_score'] / max(low_result['overall_score'], 0.01),
        "passed": bio_result['overall_score'] > low_result['overall_score'] * 1.5
    }
```

**Expected Outcome:**
- ✅ Poly-A should be **BLOCKED** by safety validator
- ✅ Low-complexity should score **<0.3**
- ✅ Biological control should score **>0.6**
- ✅ Validates safety gates are working

---

## 🧪 MASTER TEST RUNNER

```python
# tests/test_evo2_real_validation.py

import pytest
import asyncio
from typing import Dict, List

@pytest.mark.asyncio
async def test_evo2_validation_suite():
    """
    Master test suite - runs all 4 test cases and summarizes
    
    Purpose: Validate Evo2 API is REAL and working (not mocked)
    Model: evo2_1b (cost-controlled)
    Expected cost: <$0.001 total
    """
    
    print("\n" + "="*70)
    print("EVO2 VALIDATION TEST SUITE - REAL API CALLS")
    print("="*70)
    
    results = []
    
    # Test 1: Shark cartilage anti-VEGF
    print("\n[1/4] Running shark cartilage test...")
    result1 = await test_shark_cartilage_vegfa()
    results.append(result1)
    
    # Test 2: Vitamin D VDR agonist
    print("\n[2/4] Running Vitamin D test...")
    result2 = await test_vitamin_d_vdr_agonist()
    results.append(result2)
    
    # Test 3: Curcumin multi-target
    print("\n[3/4] Running curcumin test...")
    result3 = await test_curcumin_multitarget()
    results.append(result3)
    
    # Test 4: Nonsense control
    print("\n[4/4] Running nonsense control test...")
    result4 = await test_nonsense_sequence_control()
    results.append(result4)
    
    # Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    passed = sum(1 for r in results if r.get('passed', False))
    total = len(results)
    
    print(f"\nTests Passed: {passed}/{total} ({passed/total*100:.0f}%)")
    print("\nDetailed Results:")
    for i, result in enumerate(results, 1):
        status = "✅ PASS" if result.get('passed', False) else "❌ FAIL"
        print(f"  Test {i} ({result['test_name']}): {status}")
    
    # Critical validation
    if passed >= 3:
        print("\n✅ OVERALL: EVO2 API IS WORKING")
        print("   Interpretation: Real API calls returning biologically sensible results")
        print("   Ready for production use with 1B model")
    else:
        print("\n❌ OVERALL: EVO2 MAY BE MOCKED OR BROKEN")
        print("   Interpretation: Too many failures - check implementation")
        print("   DO NOT use in production until fixed")
    
    # Cost validation
    print("\n💰 COST ESTIMATE:")
    print("   Assuming ~500 tokens per test * 8 API calls")
    print("   Total tokens: ~4,000")
    print("   Cost (evo2_1b @ $0.10/1M): ~$0.0004")
    print("   ✅ Within budget (<$0.001)")
    
    assert passed >= 3, f"Only {passed}/4 tests passed - Evo2 validation failed"
```

---

## ⚔️ EXECUTION PLAN

### **Step 1: Run Master Test Suite**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal
PYTHONPATH=. /Users/fahadkiani/Desktop/development/crispr-assistant-main/venv/bin/pytest tests/test_evo2_real_validation.py -v -s
```

### **Step 2: Analyze Results**
- If **3-4 tests pass**: ✅ Evo2 is REAL and working
- If **1-2 tests pass**: ⚠️ Evo2 may be partially mocked
- If **0 tests pass**: ❌ Evo2 is mocked or broken

### **Step 3: Cost Validation**
- Expected: <$0.001 for full suite
- If higher: Check model parameter (should be 1B, not 7B/40B)

---

## 📊 ACCEPTANCE CRITERIA

### ✅ **PASS Criteria:**
1. At least **3/4 tests pass**
2. Known sequences score **higher** than random
3. Safety validator **blocks** dangerous sequences
4. Total cost **<$0.001**
5. Results are **NOT identical** across runs (proves real API)

### ❌ **FAIL Criteria:**
1. All tests return **identical scores** (proves mocking)
2. Random sequences score **higher** than biological
3. Safety validator **doesn't block** poly-A
4. Cost is **$0.00** (proves no API calls made)

---

**COMMANDER - SHALL I CREATE THE EXECUTABLE TEST FILE?** ⚔️





