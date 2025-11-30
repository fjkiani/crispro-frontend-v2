# Mechanism-Based Trial Matching: Focused Mission Plan

**Date**: January 28, 2025  
**Mission**: Deliver mechanism-based trial matching for MBD4+TP53 HGSOC patients  
**Scope**: Strip away noise, focus on core deliverable

---

## The One Thing We're Building

**Mechanism-Based Trial Matching** - matching patient tumor pathways to trial drug mechanisms.

### What We CAN Say (Our Value Proposition)

> "Based on MBD4+TP53 pathway burden (DDR: 0.88), these trials have high mechanism alignment:
> - NCT05678901: PARP + ATR inhibitor (mechanism fit: 0.92)
> - NCT04567890: Basket trial DDR-deficient (mechanism fit: 0.87)
> 
> This suggests these trials target your tumor's vulnerabilities."

### What We CANNOT Say (Not Our Scope)

> "You have 85% probability of responding to this trial."
> "NCT05678901 will extend your PFS by 6 months."

**This is mechanism alignment, not outcome prediction. Accept this constraint.**

---

## Current State (What's Working)

### ✅ Pathway Computation (P)
- Gene→pathway mapping: `pathway/aggregation.py`
- 7D mechanism vector: `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
- MBD4 mapped to DDR pathway ✅
- Code: `api/services/pathway/drug_mapping.py`

### ✅ Mechanism Fit Ranker
- Formula: `combined_score = 0.7×eligibility + 0.3×mechanism_fit`
- Cosine similarity between patient vector and trial MoA vector
- 47 trials tagged with MoA vectors
- Code: `api/services/mechanism_fit_ranker.py`

### ✅ Trial Search
- Autonomous trial agent: `api/services/autonomous_trial_agent.py`
- AstraDB semantic search
- Endpoint: `POST /api/trials/agent/search`

### ✅ Pathway→Mechanism Vector Conversion
- Function: `convert_pathway_scores_to_mechanism_vector()`
- File: `api/services/pathway_to_mechanism_vector.py`
- Maps: `{ddr, ras_mapk, pi3k, vegf, tp53}` → 7D vector

---

## The Gap (What's Missing)

### Gap 1: Mechanism Fit NOT Wired to Trial Response

**Current**: Trial search returns trials, mechanism fit computed separately.  
**Needed**: Single endpoint that returns ranked trials WITH mechanism fit scores.

**Fix**: Wire `MechanismFitRanker.rank_trials()` into trial search response.

### Gap 2: Advanced Query Generation

**Current**: Autonomous agent generates 3 basic queries.  
**Needed**: 5-10 queries including basket trials, DDR deficiency, rare mutations.

**Fix**: Add query templates to `autonomous_trial_agent.py`:
```python
QUERY_TEMPLATES = {
    'basket_trial': "{condition} basket trial tumor agnostic",
    'dna_repair': "{condition} DNA repair deficiency",
    'parp_inhibitor': "{condition} PARP inhibitor",
    'rare_mutation': "{gene} mutation rare disease registry"
}
```

### Gap 3: No Mechanism Fit in Trial Response

**Current**: Trial results don't include mechanism fit scores.  
**Needed**: Each trial has `mechanism_fit_score`, `combined_score`, `mechanism_alignment`.

**Fix**: Enhance trial response schema:
```python
class TrialResult:
    nct_id: str
    title: str
    eligibility_score: float
    mechanism_fit_score: float  # NEW
    combined_score: float        # NEW (0.7×eligibility + 0.3×mechanism_fit)
    mechanism_alignment: Dict[str, float]  # NEW (per-pathway breakdown)
```

---

## Implementation Plan (3 Days)

### Day 1: Wire Mechanism Fit to Trial Search

**Task 1.1**: Modify trial search to accept mechanism vector
```python
# In /api/trials/agent/search
request.mechanism_vector = [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
```

**Task 1.2**: Call MechanismFitRanker after search
```python
from api.services.mechanism_fit_ranker import MechanismFitRanker

ranker = MechanismFitRanker(alpha=0.7, beta=0.3)
ranked_trials = ranker.rank_trials(
    trials=search_results,
    sae_mechanism_vector=mechanism_vector,
    min_eligibility=0.60,
    min_mechanism_fit=0.50
)
```

**Task 1.3**: Add mechanism fit to response
```python
for trial in ranked_trials:
    trial['mechanism_fit_score'] = ranker.compute_mechanism_fit(...)
    trial['combined_score'] = 0.7 * eligibility + 0.3 * mechanism_fit
    trial['mechanism_alignment'] = {"DDR": 0.92, "MAPK": 0.10, ...}
```

**Files to Modify**:
- `api/services/autonomous_trial_agent.py`
- `api/routers/trials_agent.py`

### Day 2: Enhance Query Generation

**Task 2.1**: Add query templates for MBD4/DDR/basket trials
```python
# In autonomous_trial_agent.py
def generate_search_queries(self, patient_context):
    queries = []
    
    # Existing basic queries
    queries.append(f"{disease} clinical trial")
    
    # NEW: DDR/DNA repair queries
    if self._has_ddr_mutations(patient_context):
        queries.append(f"{disease} DNA repair deficiency syndrome")
        queries.append(f"{disease} PARP inhibitor")
        queries.append(f"{disease} basket trial tumor agnostic")
    
    # NEW: Rare mutation queries
    for gene in patient_context.get('rare_genes', []):
        queries.append(f"{gene} mutation rare disease registry")
    
    return queries[:10]  # Max 10 queries
```

**Task 2.2**: Detect DDR pathway mutations
```python
DDR_GENES = {'BRCA1', 'BRCA2', 'MBD4', 'ATM', 'ATR', 'CHEK2', 'RAD51', 'PALB2'}

def _has_ddr_mutations(self, context):
    genes = {m['gene'] for m in context.get('mutations', [])}
    return bool(genes & DDR_GENES)
```

**Files to Modify**:
- `api/services/autonomous_trial_agent.py`

### Day 3: Test with MBD4+TP53 Case

**Task 3.1**: Run MBD4+TP53 trial search
```python
request = {
    "mutations": [
        {"gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2"},
        {"gene": "TP53", "hgvs_p": "p.Arg175His"}
    ],
    "disease": "ovarian_cancer",
    "mechanism_vector": [0.88, 0.12, 0.15, 0.10, 0.05, 0.0, 0.0]  # DDR high
}
# POST /api/trials/agent/search
```

**Task 3.2**: Verify mechanism fit in response
```python
# Expected: DDR-targeting trials rank higher
assert trials[0]['mechanism_fit_score'] > 0.80  # PARP trial
assert trials[0]['mechanism_alignment']['DDR'] > 0.80
```

**Task 3.3**: Document results
- Create `MBD4_TRIAL_MATCHING_RESULTS.md`
- Show: queries used, trials found, mechanism fit scores
- Highlight: DDR trials ranked higher due to mechanism alignment

---

## Success Criteria

| Metric | Target | Measurement |
|--------|--------|-------------|
| MBD4 trial search works | ✅ | Returns DDR/PARP trials |
| Mechanism fit in response | ✅ | Each trial has `mechanism_fit_score` |
| DDR trials rank higher | ✅ | Top trials have DDR alignment >0.70 |
| Combined score correct | ✅ | 0.7×eligibility + 0.3×mechanism_fit |
| Response time | <10s | For <100 trials |

---

## What We're NOT Doing (Scope Exclusion)

❌ **NOT outcome prediction** - No PFS/OS correlation  
❌ **NOT efficacy validation** - Our benchmark shows r=0.037 with PFS  
❌ **NOT TRUE SAE** - Proxy SAE works, TRUE SAE is enhancement not requirement  
❌ **NOT resistance detection** - Out of scope for trial matching  
❌ **NOT drug efficacy modulation** - SAE doesn't change drug confidence  

**We focus ONLY on mechanism-based trial matching.**

---

## Key Files

| File | Purpose | Status |
|------|---------|--------|
| `mechanism_fit_ranker.py` | Core ranking logic | ✅ Ready |
| `autonomous_trial_agent.py` | Query generation | 🔄 Needs enhancement |
| `trials_agent.py` | Trial search endpoint | 🔄 Needs mechanism fit |
| `pathway_to_mechanism_vector.py` | Pathway→vector conversion | ✅ Ready |
| `drug_mapping.py` | Gene→pathway mapping | ✅ Ready |

---

## The Deliverable

**For MBD4+TP53 Patient**:

```json
{
  "query": "MBD4+TP53 ovarian cancer trial matching",
  "patient_mechanism_vector": {
    "DDR": 0.88,
    "MAPK": 0.12,
    "PI3K": 0.15,
    "VEGF": 0.10,
    "HER2": 0.05,
    "IO": 0.0,
    "Efflux": 0.0
  },
  "top_trials": [
    {
      "nct_id": "NCT05678901",
      "title": "PARP + ATR Inhibitor in DDR-Deficient Ovarian Cancer",
      "phase": "Phase 2",
      "eligibility_score": 0.85,
      "mechanism_fit_score": 0.92,
      "combined_score": 0.87,
      "mechanism_alignment": {
        "DDR": 0.95,
        "MAPK": 0.10,
        "PI3K": 0.05
      },
      "why_matched": "High DDR pathway alignment (0.95) with PARP+ATR mechanism"
    }
  ],
  "interpretation": "Trials targeting DDR pathway rank highest due to MBD4+TP53 DDR burden (0.88). This is mechanism alignment, not outcome prediction."
}
```

---

## Bottom Line

**Mission**: Mechanism-based trial matching for MBD4+TP53  
**Deliverable**: Ranked trials with mechanism fit scores  
**Timeline**: 3 days  
**Success**: DDR-targeting trials rank higher for DDR-high patients  

**No noise. No scope creep. Just mechanism-based trial matching.**

