# 🎯 TOXICITY MOAT IMPLEMENTATION - AGENT JR TASKS

**Status**: Ready for Implementation  
**Priority**: P0 (MOAT - Competitive Advantage)  
**Connection**: Implements Section 4 (Toxicity & Pharmacogenomics) + Section 7 (Nutraceutical Synergy) from ADVANCED_CARE_PLAN_EXPLAINED.md

---

## 🔍 CURRENT STATE VERIFICATION

### What EXISTS (Verified):

| Component | Location | Status |
|-----------|----------|--------|
| `compute_pathway_overlap()` | `toxicity_pathway_mappings.py:206` | ✅ Works |
| `is_pharmacogene()` | `toxicity_pathway_mappings.py:177` | ✅ Works |
| `MOA_TO_TOXIC_PATHWAYS` | `toxicity_pathway_mappings.py:73` | ✅ 9 MoAs mapped |
| `ToxicityRiskResponse` | `safety.py:66` | ✅ Exists (missing `mitigating_foods`) |
| `compute_toxicity_risk()` | `safety_service.py:32` | ✅ Works (doesn't call `get_mitigating_foods`) |

### What DOES NOT EXIST (Must Create):

| Component | Location | Status |
|-----------|----------|--------|
| `get_mitigating_foods()` | `toxicity_pathway_mappings.py` | ❌ **MISSING** |
| `DRUG_TO_MOA` mapping | `toxicity_pathway_mappings.py` | ❌ **MISSING** |
| `get_drug_moa()` | `toxicity_pathway_mappings.py` | ❌ **MISSING** |
| `mitigating_foods` field | `ToxicityRiskResponse` schema | ❌ **MISSING** |
| Call to `get_mitigating_foods()` | `safety_service.py` | ❌ **MISSING** |
| `toxicity_mitigation` field | `hypothesis_validator.py` response | ❌ **MISSING** |
| Toxicity badge display | `FoodRankingPanel.jsx` | ❌ **MISSING** |

---

## 📋 TASK LIST FOR AGENT JR

### Pre-Flight Check (Run First)

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 -c "
print('=== PRE-FLIGHT CHECK ===')

# Test 1: compute_pathway_overlap works
from api.services.toxicity_pathway_mappings import compute_pathway_overlap
overlap = compute_pathway_overlap(['BRCA1'], 'platinum_agent')
print(f'✅ compute_pathway_overlap works: {overlap}')

# Test 2: Check if get_mitigating_foods exists
try:
    from api.services.toxicity_pathway_mappings import get_mitigating_foods
    print('✅ get_mitigating_foods EXISTS')
except ImportError:
    print('❌ get_mitigating_foods MISSING - Task 1 will create it')

# Test 3: Check if mitigating_foods field exists in schema
from api.schemas.safety import ToxicityRiskResponse
if 'mitigating_foods' in ToxicityRiskResponse.__fields__:
    print('✅ mitigating_foods field EXISTS in schema')
else:
    print('❌ mitigating_foods field MISSING - Task 2 will add it')

print('\\n=== PRE-FLIGHT COMPLETE ===')
"
```

---

## TASK 1: Add `get_mitigating_foods()` Function

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/toxicity_pathway_mappings.py`

**Action**: Add this code at the END of the file (after line 234):

```python
# ============================================================================
# MITIGATING FOODS (THE MOAT)
# ============================================================================

def get_mitigating_foods(pathway_overlap: Dict[str, float]) -> List[Dict[str, Any]]:
    """
    Map toxicity pathway overlap to mitigating foods.
    
    THE MOAT - connecting toxicity assessment to food recommendations.
    
    This implements Section 4 + Section 7 from ADVANCED_CARE_PLAN_EXPLAINED.md:
    - Section 4: Toxicity & Pharmacogenomics - detects drug toxicity risks
    - Section 7: Nutraceutical Synergy/Antagonism - recommends timing and foods
    
    Args:
        pathway_overlap: Dict from compute_pathway_overlap() 
                        e.g., {"dna_repair": 1.0, "inflammation": 0.0, "cardiometabolic": 0.0}
    
    Returns:
        List of mitigating food recommendations with timing guidance
    """
    from typing import Any  # Ensure Any is imported
    recommendations = []
    
    # DNA REPAIR pathway - mitigating foods
    # From ADVANCED_CARE_PLAN Section 7: "Vitamin D for HRD context repletion (DNA repair support)"
    if pathway_overlap.get("dna_repair", 0) > 0.3:
        recommendations.extend([
            {
                "compound": "NAC (N-Acetyl Cysteine)",
                "dose": "600mg twice daily",
                "timing": "post-chemo (not during infusion)",
                "mechanism": "Glutathione precursor, supports DNA repair enzymes",
                "evidence_tier": "MODERATE",
                "pathway": "dna_repair",
                "care_plan_ref": "Section 7: Antioxidants - After chemo: Use - helps recovery"
            },
            {
                "compound": "Vitamin D3",
                "dose": "5000 IU daily",
                "timing": "continuous, with fatty meal",
                "mechanism": "Modulates DNA repair gene expression (VDR-mediated)",
                "evidence_tier": "MODERATE",
                "pathway": "dna_repair",
                "care_plan_ref": "Section 7: Vitamin D - HRD context repletion (DNA repair support)"
            },
            {
                "compound": "Folate (5-MTHF)",
                "dose": "400-800mcg daily",
                "timing": "continuous",
                "mechanism": "DNA synthesis and repair cofactor",
                "evidence_tier": "MODERATE",
                "pathway": "dna_repair",
                "care_plan_ref": "DNA repair support - continuous supplementation"
            }
        ])
    
    # INFLAMMATION pathway - mitigating foods
    # From ADVANCED_CARE_PLAN Section 7: "Omega-3 - Post-chemo inflammation control"
    if pathway_overlap.get("inflammation", 0) > 0.3:
        recommendations.extend([
            {
                "compound": "Omega-3 (EPA+DHA)",
                "dose": "2-3g combined EPA+DHA daily",
                "timing": "post-infusion (anti-inflammatory)",
                "mechanism": "Resolvin precursor, inhibits NF-κB, reduces IL-6",
                "evidence_tier": "SUPPORTED",
                "pathway": "inflammation",
                "care_plan_ref": "Section 7: Omega-3 - Post-chemo inflammation control"
            },
            {
                "compound": "Curcumin (Turmeric Extract)",
                "dose": "500-1000mg daily (with piperine)",
                "timing": "between meals, post-chemo",
                "mechanism": "NF-κB inhibitor, COX-2 inhibitor, reduces cytokine storm",
                "evidence_tier": "MODERATE",
                "pathway": "inflammation",
                "care_plan_ref": "Anti-inflammatory - for checkpoint inhibitor iRAEs"
            },
            {
                "compound": "EGCG (Green Tea Extract)",
                "dose": "400-800mg daily",
                "timing": "between meals (not with iron supplements)",
                "mechanism": "Anti-inflammatory, STAT3 inhibitor",
                "evidence_tier": "MODERATE",
                "pathway": "inflammation",
                "care_plan_ref": "STAT3 pathway modulation"
            }
        ])
    
    # CARDIOMETABOLIC pathway - mitigating foods  
    # From ADVANCED_CARE_PLAN Section 4: Anthracycline cardiotoxicity prevention
    if pathway_overlap.get("cardiometabolic", 0) > 0.3:
        recommendations.extend([
            {
                "compound": "CoQ10 (Ubiquinol)",
                "dose": "200-400mg daily",
                "timing": "with fatty meal, continuous during anthracycline",
                "mechanism": "Mitochondrial support, cardioprotective against doxorubicin",
                "evidence_tier": "SUPPORTED",
                "pathway": "cardiometabolic",
                "care_plan_ref": "Section 4: Anthracycline cardiotoxicity - CoQ10 for protection"
            },
            {
                "compound": "L-Carnitine",
                "dose": "1000-2000mg daily",
                "timing": "morning, with breakfast",
                "mechanism": "Fatty acid transport, cardiac energy metabolism",
                "evidence_tier": "MODERATE",
                "pathway": "cardiometabolic",
                "care_plan_ref": "Cardiac metabolism support"
            },
            {
                "compound": "Magnesium Glycinate",
                "dose": "400mg daily",
                "timing": "evening (aids sleep, cardiac rhythm)",
                "mechanism": "Cardiac rhythm support, ATP synthesis, QT stabilization",
                "evidence_tier": "MODERATE",
                "pathway": "cardiometabolic",
                "care_plan_ref": "QT prolongation risk mitigation"
            }
        ])
    
    return recommendations
```

**Also update imports at top of file** (line 13):

Change:
```python
from typing import Dict, List, Set
```

To:
```python
from typing import Dict, List, Set, Any
```

**Test Command**:

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 -c "
from api.services.toxicity_pathway_mappings import get_mitigating_foods, compute_pathway_overlap

# Test 1: DNA repair pathway (platinum + BRCA)
overlap1 = compute_pathway_overlap(['BRCA1'], 'platinum_agent')
foods1 = get_mitigating_foods(overlap1)
print('Test 1 - BRCA1 + platinum:')
print(f'  Pathway overlap: {overlap1}')
print(f'  Mitigating foods: {len(foods1)}')
for f in foods1:
    print(f'    - {f[\"compound\"]}: {f[\"timing\"]}')

# Test 2: Inflammation pathway (checkpoint inhibitor)
overlap2 = compute_pathway_overlap(['TNF', 'IL6'], 'checkpoint_inhibitor')
foods2 = get_mitigating_foods(overlap2)
print(f'\\nTest 2 - TNF/IL6 + checkpoint_inhibitor:')
print(f'  Mitigating foods: {len(foods2)}')
for f in foods2:
    print(f'    - {f[\"compound\"]}')

# Test 3: Cardiometabolic pathway (anthracycline)
overlap3 = compute_pathway_overlap(['KCNQ1', 'RYR2'], 'anthracycline')
foods3 = get_mitigating_foods(overlap3)
print(f'\\nTest 3 - KCNQ1/RYR2 + anthracycline:')
print(f'  Mitigating foods: {len(foods3)}')
for f in foods3:
    print(f'    - {f[\"compound\"]}')

print('\\n✅ TASK 1 COMPLETE!')
"
```

**Expected Output**:
```
Test 1 - BRCA1 + platinum:
  Pathway overlap: {'dna_repair': 1.0, 'inflammation': 0.0, 'cardiometabolic': 0.0}
  Mitigating foods: 3
    - NAC (N-Acetyl Cysteine): post-chemo (not during infusion)
    - Vitamin D3: continuous, with fatty meal
    - Folate (5-MTHF): continuous

Test 2 - TNF/IL6 + checkpoint_inhibitor:
  Mitigating foods: 3
    - Omega-3 (EPA+DHA)
    - Curcumin (Turmeric Extract)
    - EGCG (Green Tea Extract)

Test 3 - KCNQ1/RYR2 + anthracycline:
  Mitigating foods: 3
    - CoQ10 (Ubiquinol)
    - L-Carnitine
    - Magnesium Glycinate

✅ TASK 1 COMPLETE!
```

---

## TASK 2: Add `mitigating_foods` Field to Schema

**File**: `oncology-coPilot/oncology-backend-minimal/api/schemas/safety.py`

**Action**: Update `ToxicityRiskResponse` class (around line 66-74):

Change:
```python
class ToxicityRiskResponse(BaseModel):
    """Response from toxicity risk assessment."""
    risk_score: float = Field(..., description="Overall toxicity risk (0-1, higher = more risk)")
    confidence: float = Field(..., description="Confidence in assessment (0-1)")
    reason: str = Field(..., description="Plain-English summary of risk")
    factors: List[ToxicityFactor] = Field(default_factory=list, description="Individual risk factors")
    evidence: Dict[str, Any] = Field(default_factory=dict, description="Citations and badges")
    provenance: Dict[str, Any] = Field(default_factory=dict, description="Run ID, profile, methods, cache")
```

To:
```python
class ToxicityRiskResponse(BaseModel):
    """Response from toxicity risk assessment."""
    risk_score: float = Field(..., description="Overall toxicity risk (0-1, higher = more risk)")
    confidence: float = Field(..., description="Confidence in assessment (0-1)")
    reason: str = Field(..., description="Plain-English summary of risk")
    factors: List[ToxicityFactor] = Field(default_factory=list, description="Individual risk factors")
    mitigating_foods: List[Dict[str, Any]] = Field(
        default_factory=list, 
        description="Foods that mitigate detected toxicity pathways (THE MOAT)"
    )
    evidence: Dict[str, Any] = Field(default_factory=dict, description="Citations and badges")
    provenance: Dict[str, Any] = Field(default_factory=dict, description="Run ID, profile, methods, cache")
```

**Test Command**: No direct test - schema validation happens at runtime with Task 3.

---

## TASK 3: Call `get_mitigating_foods()` in Safety Service

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/safety_service.py`

**Action 1**: Update imports (line 20-23):

Change:
```python
from api.services.toxicity_pathway_mappings import (
    is_pharmacogene, get_pharmacogene_risk_weight,
    compute_pathway_overlap, get_moa_toxicity_weights
)
```

To:
```python
from api.services.toxicity_pathway_mappings import (
    is_pharmacogene, get_pharmacogene_risk_weight,
    compute_pathway_overlap, get_moa_toxicity_weights,
    get_mitigating_foods  # NEW - THE MOAT
)
```

**Action 2**: Store pathway_overlaps for later use. Find line 60-72 and update:

Change:
```python
        # Factor 2: MoA → Toxicity pathway overlap (P signal)
        pathway_weight = 0.0
        if request.candidate.moa and germline_genes:
            pathway_overlaps = compute_pathway_overlap(germline_genes, request.candidate.moa)
```

To:
```python
        # Factor 2: MoA → Toxicity pathway overlap (P signal)
        pathway_weight = 0.0
        pathway_overlaps = {}  # Store for mitigating foods calculation
        if request.candidate.moa and germline_genes:
            pathway_overlaps = compute_pathway_overlap(germline_genes, request.candidate.moa)
```

**Action 3**: Add mitigating foods calculation before return statement (around line 130):

Add this BEFORE the return statement:
```python
        # THE MOAT: Get mitigating foods based on detected toxicity pathways
        mitigating_foods = get_mitigating_foods(pathway_overlaps) if pathway_overlaps else []
```

**Action 4**: Update return statement to include mitigating_foods:

Change:
```python
        return ToxicityRiskResponse(
            risk_score=risk_score,
            confidence=confidence,
            reason=reason,
            factors=factors,
            evidence=evidence,
            provenance=provenance
        )
```

To:
```python
        return ToxicityRiskResponse(
            risk_score=risk_score,
            confidence=confidence,
            reason=reason,
            factors=factors,
            mitigating_foods=mitigating_foods,  # THE MOAT
            evidence=evidence,
            provenance=provenance
        )
```

**Test Command**:

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 -c "
import asyncio
from api.services.safety_service import get_safety_service
from api.schemas.safety import (
    ToxicityRiskRequest, PatientContext, GermlineVariant,
    TherapeuticCandidate, ClinicalContext
)

async def test():
    service = get_safety_service()
    
    # Test: BRCA1 + platinum_agent
    request = ToxicityRiskRequest(
        patient=PatientContext(
            germlineVariants=[
                GermlineVariant(chrom='17', pos=41276045, ref='A', alt='G', gene='BRCA1')
            ]
        ),
        candidate=TherapeuticCandidate(type='drug', moa='platinum_agent'),
        context=ClinicalContext(disease='ovarian_cancer'),
        options={'profile': 'baseline'}
    )
    
    result = await service.compute_toxicity_risk(request)
    
    print('=== TOXICITY RISK RESPONSE ===')
    print(f'Risk Score: {result.risk_score}')
    print(f'Confidence: {result.confidence}')
    print(f'Reason: {result.reason}')
    print(f'\\nMitigating Foods: {len(result.mitigating_foods)}')
    for food in result.mitigating_foods:
        print(f'  - {food[\"compound\"]}')
        print(f'    Dose: {food[\"dose\"]}')
        print(f'    Timing: {food[\"timing\"]}')
        print(f'    Mechanism: {food[\"mechanism\"]}')
    
    if len(result.mitigating_foods) >= 3:
        print('\\n✅ TASK 3 COMPLETE - Mitigating foods returned!')
    else:
        print('\\n❌ TASK 3 FAILED - Expected 3+ mitigating foods')

asyncio.run(test())
"
```

**Expected Output**:
```
=== TOXICITY RISK RESPONSE ===
Risk Score: 1.0
Confidence: 0.48
Reason: MoA overlaps toxicity pathways with germline variants

Mitigating Foods: 3
  - NAC (N-Acetyl Cysteine)
    Dose: 600mg twice daily
    Timing: post-chemo (not during infusion)
    Mechanism: Glutathione precursor, supports DNA repair enzymes
  - Vitamin D3
    Dose: 5000 IU daily
    Timing: continuous, with fatty meal
    Mechanism: Modulates DNA repair gene expression (VDR-mediated)
  - Folate (5-MTHF)
    Dose: 400-800mcg daily
    Timing: continuous
    Mechanism: DNA synthesis and repair cofactor

✅ TASK 3 COMPLETE - Mitigating foods returned!
```

---

## TASK 4: Add Drug → MoA Lookup Helper

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/toxicity_pathway_mappings.py`

**Action**: Add this code AFTER `get_mitigating_foods()` function:

```python
# ============================================================================
# DRUG TO MOA MAPPING
# ============================================================================

DRUG_TO_MOA: Dict[str, str] = {
    # Platinum agents
    "carboplatin": "platinum_agent",
    "cisplatin": "platinum_agent",
    "oxaliplatin": "platinum_agent",
    
    # Anthracyclines (HIGH cardiotoxicity - Section 4)
    "doxorubicin": "anthracycline",
    "adriamycin": "anthracycline",
    "epirubicin": "anthracycline",
    "daunorubicin": "anthracycline",
    
    # PARP inhibitors
    "olaparib": "PARP_inhibitor",
    "niraparib": "PARP_inhibitor",
    "rucaparib": "PARP_inhibitor",
    "talazoparib": "PARP_inhibitor",
    
    # Checkpoint inhibitors (inflammation - iRAEs)
    "pembrolizumab": "checkpoint_inhibitor",
    "nivolumab": "checkpoint_inhibitor",
    "atezolizumab": "checkpoint_inhibitor",
    "ipilimumab": "checkpoint_inhibitor",
    "durvalumab": "checkpoint_inhibitor",
    "avelumab": "checkpoint_inhibitor",
    
    # Alkylating agents
    "cyclophosphamide": "alkylating_agent",
    "temozolomide": "alkylating_agent",
    "ifosfamide": "alkylating_agent",
    "bendamustine": "alkylating_agent",
    
    # BRAF/MEK inhibitors
    "vemurafenib": "BRAF_inhibitor",
    "dabrafenib": "BRAF_inhibitor",
    "encorafenib": "BRAF_inhibitor",
    "trametinib": "MEK_inhibitor",
    "cobimetinib": "MEK_inhibitor",
    "binimetinib": "MEK_inhibitor",
    
    # Proteasome inhibitors
    "bortezomib": "proteasome_inhibitor",
    "carfilzomib": "proteasome_inhibitor",
    "ixazomib": "proteasome_inhibitor",
    
    # IMiDs
    "lenalidomide": "immunomodulatory",
    "pomalidomide": "immunomodulatory",
    "thalidomide": "immunomodulatory",
}


def get_drug_moa(drug_name: str) -> str:
    """
    Get mechanism of action for a drug.
    
    Args:
        drug_name: Drug name (case-insensitive)
    
    Returns:
        MoA string or "unknown" if not found
    """
    return DRUG_TO_MOA.get(drug_name.lower().strip(), "unknown")
```

**Test Command**:

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 -c "
from api.services.toxicity_pathway_mappings import get_drug_moa

tests = [
    ('Carboplatin', 'platinum_agent'),
    ('DOXORUBICIN', 'anthracycline'),
    ('olaparib', 'PARP_inhibitor'),
    ('Pembrolizumab', 'checkpoint_inhibitor'),
    ('vemurafenib', 'BRAF_inhibitor'),
    ('unknown_drug', 'unknown'),
]

print('=== DRUG → MOA LOOKUP ===')
for drug, expected in tests:
    result = get_drug_moa(drug)
    status = '✅' if result == expected else '❌'
    print(f'{status} {drug} → {result}')

print('\\n✅ TASK 4 COMPLETE!')
"
```

---

## TASK 5: Create Integration Test

**File**: Create `oncology-coPilot/oncology-backend-minimal/test_moat_integration.py`

```python
"""
MOAT Integration Test - Validates the complete toxicity → food mitigation flow.

This test validates the implementation of:
- ADVANCED_CARE_PLAN Section 4: Toxicity & Pharmacogenomics
- ADVANCED_CARE_PLAN Section 7: Nutraceutical Synergy/Antagonism
"""

import asyncio
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent))


async def test_complete_moat():
    print("=" * 70)
    print("MOAT INTEGRATION TEST")
    print("Testing: Toxicity → Mitigating Foods Flow")
    print("=" * 70)
    
    # =========================================================================
    # TEST 1: Core Functions Exist
    # =========================================================================
    print("\n[TEST 1] Core Functions Exist")
    print("-" * 50)
    
    from api.services.toxicity_pathway_mappings import (
        compute_pathway_overlap,
        get_mitigating_foods,
        get_drug_moa
    )
    print("   ✅ compute_pathway_overlap imported")
    print("   ✅ get_mitigating_foods imported")
    print("   ✅ get_drug_moa imported")
    
    # =========================================================================
    # TEST 2: Pathway → Foods Mapping
    # =========================================================================
    print("\n[TEST 2] Pathway → Foods Mapping")
    print("-" * 50)
    
    # DNA repair
    dna_foods = get_mitigating_foods({"dna_repair": 1.0, "inflammation": 0.0, "cardiometabolic": 0.0})
    assert len(dna_foods) == 3, f"Expected 3 DNA repair foods, got {len(dna_foods)}"
    print(f"   ✅ DNA repair → {len(dna_foods)} foods: {[f['compound'] for f in dna_foods]}")
    
    # Inflammation
    inf_foods = get_mitigating_foods({"dna_repair": 0.0, "inflammation": 1.0, "cardiometabolic": 0.0})
    assert len(inf_foods) == 3, f"Expected 3 inflammation foods, got {len(inf_foods)}"
    print(f"   ✅ Inflammation → {len(inf_foods)} foods: {[f['compound'] for f in inf_foods]}")
    
    # Cardiometabolic
    cardio_foods = get_mitigating_foods({"dna_repair": 0.0, "inflammation": 0.0, "cardiometabolic": 1.0})
    assert len(cardio_foods) == 3, f"Expected 3 cardio foods, got {len(cardio_foods)}"
    print(f"   ✅ Cardiometabolic → {len(cardio_foods)} foods: {[f['compound'] for f in cardio_foods]}")
    
    # =========================================================================
    # TEST 3: Drug → MoA Lookup
    # =========================================================================
    print("\n[TEST 3] Drug → MoA Lookup")
    print("-" * 50)
    
    assert get_drug_moa("carboplatin") == "platinum_agent"
    assert get_drug_moa("doxorubicin") == "anthracycline"
    assert get_drug_moa("pembrolizumab") == "checkpoint_inhibitor"
    assert get_drug_moa("olaparib") == "PARP_inhibitor"
    assert get_drug_moa("unknown") == "unknown"
    
    print("   ✅ carboplatin → platinum_agent")
    print("   ✅ doxorubicin → anthracycline")
    print("   ✅ pembrolizumab → checkpoint_inhibitor")
    print("   ✅ olaparib → PARP_inhibitor")
    print("   ✅ unknown → unknown (fallback works)")
    
    # =========================================================================
    # TEST 4: Full API Response Includes mitigating_foods
    # =========================================================================
    print("\n[TEST 4] Safety API Returns mitigating_foods")
    print("-" * 50)
    
    from api.services.safety_service import get_safety_service
    from api.schemas.safety import (
        ToxicityRiskRequest, PatientContext, GermlineVariant,
        TherapeuticCandidate, ClinicalContext
    )
    
    service = get_safety_service()
    
    # Scenario: BRCA1 patient on carboplatin (platinum_agent)
    request = ToxicityRiskRequest(
        patient=PatientContext(
            germlineVariants=[
                GermlineVariant(chrom="17", pos=41276045, ref="A", alt="G", gene="BRCA1")
            ]
        ),
        candidate=TherapeuticCandidate(type="drug", moa="platinum_agent"),
        context=ClinicalContext(disease="ovarian_cancer"),
        options={"profile": "baseline"}
    )
    
    result = await service.compute_toxicity_risk(request)
    
    assert hasattr(result, 'mitigating_foods'), "mitigating_foods field missing from response"
    assert len(result.mitigating_foods) >= 3, f"Expected 3+ mitigating foods, got {len(result.mitigating_foods)}"
    
    print(f"   ✅ Risk Score: {result.risk_score}")
    print(f"   ✅ mitigating_foods field present: {len(result.mitigating_foods)} foods")
    for food in result.mitigating_foods:
        print(f"      - {food['compound']}: {food['timing']}")
    
    # =========================================================================
    # TEST 5: Clinical Scenarios
    # =========================================================================
    print("\n[TEST 5] Clinical Scenarios")
    print("-" * 50)
    
    # Scenario A: Anthracycline patient (cardiotoxicity)
    request_anthracycline = ToxicityRiskRequest(
        patient=PatientContext(
            germlineVariants=[
                GermlineVariant(chrom="7", pos=100000, ref="C", alt="T", gene="KCNQ1")
            ]
        ),
        candidate=TherapeuticCandidate(type="drug", moa="anthracycline"),
        context=ClinicalContext(disease="breast_cancer"),
        options={"profile": "baseline"}
    )
    
    result_anthracycline = await service.compute_toxicity_risk(request_anthracycline)
    cardio_recommended = [f for f in result_anthracycline.mitigating_foods if f["pathway"] == "cardiometabolic"]
    
    print(f"   Scenario A: Anthracycline + KCNQ1 variant")
    print(f"   ✅ Risk Score: {result_anthracycline.risk_score}")
    print(f"   ✅ Cardiometabolic foods: {len(cardio_recommended)}")
    assert any("CoQ10" in f["compound"] for f in cardio_recommended), "CoQ10 should be recommended for cardiotoxicity"
    print(f"   ✅ CoQ10 recommended for cardioprotection")
    
    # Scenario B: Checkpoint inhibitor patient (inflammation/iRAEs)
    request_io = ToxicityRiskRequest(
        patient=PatientContext(
            germlineVariants=[
                GermlineVariant(chrom="6", pos=31575000, ref="G", alt="A", gene="TNF")
            ]
        ),
        candidate=TherapeuticCandidate(type="drug", moa="checkpoint_inhibitor"),
        context=ClinicalContext(disease="melanoma"),
        options={"profile": "baseline"}
    )
    
    result_io = await service.compute_toxicity_risk(request_io)
    inflammation_recommended = [f for f in result_io.mitigating_foods if f["pathway"] == "inflammation"]
    
    print(f"\n   Scenario B: Checkpoint Inhibitor + TNF variant")
    print(f"   ✅ Risk Score: {result_io.risk_score}")
    print(f"   ✅ Inflammation foods: {len(inflammation_recommended)}")
    assert any("Omega-3" in f["compound"] for f in inflammation_recommended), "Omega-3 should be recommended for inflammation"
    print(f"   ✅ Omega-3 recommended for iRAE prevention")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 70)
    print("✅ ALL MOAT INTEGRATION TESTS PASSED!")
    print("=" * 70)
    print("""
MOAT CAPABILITY VALIDATED:
- Toxicity detection → Pathway overlap calculation
- Pathway overlap → Mitigating foods recommendation
- Drug name → MoA lookup
- Full API response includes mitigating_foods
- Clinical scenarios (cardiotoxicity, iRAEs) work correctly

NEXT STEPS:
- Task 5: Add toxicity_mitigation to food validation response
- Task 6: Update FoodRankingPanel.jsx to display mitigation badge
""")


if __name__ == "__main__":
    asyncio.run(test_complete_moat())
```

**Run Test**:

```bash
cd oncology-coPilot/oncology-backend-minimal && python3 test_moat_integration.py
```

---

## 📊 TASK CHECKLIST

### Phase 1: Backend MOAT (COMPLETE ✅)

| # | Task | File | Status |
|---|------|------|--------|
| 1 | Add `get_mitigating_foods()` | `toxicity_pathway_mappings.py` | ✅ |
| 2 | Add `mitigating_foods` field | `safety.py` schema | ✅ |
| 3 | Call `get_mitigating_foods()` | `safety_service.py` | ✅ |
| 4 | Add `get_drug_moa()` | `toxicity_pathway_mappings.py` | ✅ |
| 5 | Create integration test | `test_moat_integration.py` | ✅ |

**Execution Order**: 1 → 2 → 3 → 4 → 5 (test)

### Phase 2: Food Validation Integration (COMPLETE ✅)

| # | Task | File | Status |
|---|------|------|--------|
| 6 | Add `toxicity_mitigation` to food validation response | `hypothesis_validator.py` | ✅ |
| 7 | Update FoodRankingPanel.jsx to display toxicity badge | `FoodRankingPanel.jsx` | ✅ |
| 8 | Create Phase 2 integration test | `test_toxicity_food_integration.py` | ✅ |

**Execution Order**: 6 → 7 → 8 (test)

---

## 🔗 CONNECTION TO ADVANCED CARE PLAN

This implementation directly addresses:

| Advanced Care Plan Section | What We're Implementing |
|---------------------------|------------------------|
| **Section 4: Toxicity & Pharmacogenomics** | Detecting toxicity pathways from germline + drug MoA |
| **Section 7: Nutraceutical Synergy/Antagonism** | Timing guide + mitigating foods |

**The MOAT**: No competitor has integrated toxicity detection with food recommendations. This is "what foods mitigate YOUR drug's toxicity for YOUR germline profile."

---

## 🚀 AGENT JR: START HERE

1. **Run pre-flight check** (verify current state)
2. **Execute Task 1** (add `get_mitigating_foods()`)
3. **Run Task 1 test** (verify it works)
4. **Execute Task 2** (add schema field)
5. **Execute Task 3** (wire up in safety_service)
6. **Run Task 3 test** (verify API returns foods)
7. **Execute Task 4** (add drug→MoA helper)
8. **Run Task 4 test** (verify lookups work)
9. **Execute Task 5** (create integration test)
10. **Run integration test** (verify everything works)

---

---

## PHASE 2: FOOD VALIDATION INTEGRATION (COMPLETE ✅)

### TASK 6: Add `toxicity_mitigation` to Food Validation Response

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

**Implementation**: Added toxicity mitigation check in `validate_food_dynamic` endpoint (after dietician recommendations, before response assembly).

**Logic**:
1. Extract genes from `disease_context.mutations` (conservative: treat as potential germline)
2. For each medication in `patient_medications`:
   - Get drug MoA using `get_drug_moa()`
   - Compute pathway overlap using `compute_pathway_overlap()`
   - Get mitigating foods using `get_mitigating_foods()`
   - Check if current compound matches any mitigating food (flexible matching)
3. Add `toxicity_mitigation` field to response if match found

**Response Field**:
```python
"toxicity_mitigation": {
    "mitigates": True,
    "target_drug": "carboplatin",
    "target_moa": "platinum_agent",
    "pathway": "dna_repair",
    "mechanism": "Glutathione precursor, supports DNA repair enzymes",
    "timing": "post-chemo (not during infusion)",
    "evidence_tier": "MODERATE",
    "dose": "600mg twice daily",
    "care_plan_ref": "Section 7: Antioxidants - After chemo: Use - helps recovery"
}
```

### TASK 7: Update FoodRankingPanel.jsx

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/FoodRankingPanel.jsx`

**Changes**:
1. Added `SecurityIcon` import
2. Added toxicity mitigation badge display (before biomarker matches)
3. Updated PropTypes to include `toxicity_mitigation` field

**Display**: Green chip with security icon showing "Mitigates {drug} {pathway} stress"

### TASK 8: Create Phase 2 Integration Test

**File**: `oncology-coPilot/oncology-backend-minimal/test_toxicity_food_integration.py`

**Tests**:
- NAC validation for BRCA1 + carboplatin → should show toxicity_mitigation
- Vitamin D validation for BRCA1 + carboplatin → should show toxicity_mitigation

**Run Test**:
```bash
cd oncology-coPilot/oncology-backend-minimal
python3 test_toxicity_food_integration.py
```

**Note**: Requires API server running (`uvicorn api.main:app --reload`)

---

**Created**: 2025-01-XX  
**Status**: ✅ Phase 1 + Phase 2 COMPLETE - Ready for Phase 3 (LLM Enhancement)

---

---

## PHASE 3: LLM ENHANCEMENT FOR TOXICITY MOAT 🆕

**Goal**: Use LLM (`src/tools/llm_api.py`) to generate personalized, context-aware rationales for toxicity mitigation.

**Current State**:
- ✅ `get_mitigating_foods()` returns structured food recommendations
- ✅ `toxicity_mitigation` field in food validation response
- ❌ Rationales are static/templated
- ❌ No patient-specific synthesis
- ❌ No plain-English explanations for patients

**LLM API Available** (`src/tools/llm_api.py`):
- Providers: Gemini
- Function: `query_llm(prompt, provider="gemini")`
- Supports: Image input, conversation history

---

### TASK 9: Create LLM Toxicity Rationale Service

**File**: Create `oncology-coPilot/oncology-backend-minimal/api/services/llm_toxicity_service.py`

**Purpose**: Generate LLM-enhanced rationales for toxicity mitigation

```python
"""
LLM-Enhanced Toxicity Rationale Service

Uses LLM to generate personalized, context-aware explanations for:
1. Why a food mitigates drug toxicity
2. Patient-friendly summaries
3. Timing/dosage rationale based on treatment context
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Any, Optional

# Add project root to path for LLM API access
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    from src.tools.llm_api import query_llm
    LLM_AVAILABLE = True
except ImportError:
    LLM_AVAILABLE = False
    print("Warning: LLM API not available. Falling back to static rationales.")


TOXICITY_RATIONALE_PROMPT = """You are a precision oncology nutritionist. Generate a concise, evidence-based explanation for why this food/supplement helps mitigate drug-induced toxicity.

PATIENT CONTEXT:
- Cancer Type: {cancer_type}
- Current Drug: {drug_name} (Mechanism: {drug_moa})
- Germline Variants: {germline_genes}
- Toxicity Pathway: {toxicity_pathway}
- Treatment Phase: {treatment_phase}

FOOD/SUPPLEMENT:
- Compound: {compound}
- Mechanism: {base_mechanism}
- Recommended Timing: {timing}
- Dose: {dose}

INSTRUCTIONS:
1. Explain why this compound specifically helps with {drug_name} toxicity
2. Mention the pathway connection ({toxicity_pathway})
3. Explain optimal timing around chemotherapy
4. Keep to 2-3 sentences, patient-friendly language
5. Do NOT claim it "cures" or "treats" - use "may help support", "can help protect"

Generate a personalized rationale:"""


PATIENT_SUMMARY_PROMPT = """You are explaining toxicity mitigation to a cancer patient in simple terms.

The patient is on {drug_name} chemotherapy.
We've identified that {compound} may help protect against {toxicity_type}.

Write a 2-sentence patient-friendly explanation that:
1. Explains what the drug does that causes this side effect
2. Explains how the supplement may help
3. Uses simple language (8th grade reading level)
4. Avoids medical jargon

Patient explanation:"""


async def generate_toxicity_rationale(
    compound: str,
    drug_name: str,
    drug_moa: str,
    toxicity_pathway: str,
    germline_genes: List[str],
    cancer_type: str = "cancer",
    treatment_phase: str = "active treatment",
    base_mechanism: str = "",
    timing: str = "",
    dose: str = "",
    provider: str = "gemini"
) -> Dict[str, Any]:
    """
    Generate LLM-enhanced rationale for toxicity mitigation.
    
    Returns:
        Dict with:
        - rationale: Personalized mechanism explanation
        - patient_summary: Patient-friendly explanation
        - confidence: Confidence in the recommendation
        - llm_enhanced: Whether LLM was used
    """
    
    result = {
        "rationale": base_mechanism,  # Fallback
        "patient_summary": f"{compound} may help support your body during {drug_name} treatment.",
        "confidence": 0.6,
        "llm_enhanced": False
    }
    
    if not LLM_AVAILABLE:
        return result
    
    try:
        # Generate personalized rationale
        rationale_prompt = TOXICITY_RATIONALE_PROMPT.format(
            cancer_type=cancer_type,
            drug_name=drug_name,
            drug_moa=drug_moa,
            toxicity_pathway=toxicity_pathway,
            germline_genes=", ".join(germline_genes) if germline_genes else "none identified",
            treatment_phase=treatment_phase,
            compound=compound,
            base_mechanism=base_mechanism,
            timing=timing,
            dose=dose
        )
        
        rationale = query_llm(rationale_prompt, provider=provider)
        
        if rationale and not rationale.startswith("Error"):
            result["rationale"] = rationale.strip()
            result["llm_enhanced"] = True
            result["confidence"] = 0.75
        
        # Generate patient summary
        toxicity_type_map = {
            "dna_repair": "DNA damage",
            "inflammation": "inflammation and immune reactions",
            "cardiometabolic": "heart stress"
        }
        toxicity_type = toxicity_type_map.get(toxicity_pathway, "side effects")
        
        patient_prompt = PATIENT_SUMMARY_PROMPT.format(
            drug_name=drug_name,
            compound=compound,
            toxicity_type=toxicity_type
        )
        
        patient_summary = query_llm(patient_prompt, provider=provider)
        
        if patient_summary and not patient_summary.startswith("Error"):
            result["patient_summary"] = patient_summary.strip()
    
    except Exception as e:
        print(f"LLM rationale generation failed: {e}")
        # Keep fallback values
    
    return result


async def generate_mitigation_dossier(
    patient_context: Dict[str, Any],
    medications: List[str],
    mitigating_foods: List[Dict[str, Any]],
    provider: str = "gemini"
) -> Dict[str, Any]:
    """
    Generate a complete toxicity mitigation dossier using LLM.
    
    This creates a structured document similar to MBD4_TP53_CLINICAL_DOSSIER.md
    but focused on toxicity mitigation and food recommendations.
    
    Returns:
        Dict with dossier sections:
        - executive_summary
        - toxicity_assessment
        - food_recommendations (with LLM rationales)
        - timing_protocol
        - monitoring_recommendations
    """
    
    germline_genes = patient_context.get("germline_genes", [])
    cancer_type = patient_context.get("cancer_type", "cancer")
    treatment_line = patient_context.get("treatment_line", "active treatment")
    
    dossier = {
        "executive_summary": "",
        "toxicity_assessment": {},
        "food_recommendations": [],
        "timing_protocol": "",
        "monitoring_recommendations": "",
        "llm_enhanced": False
    }
    
    # Enhance each food recommendation with LLM rationale
    enhanced_foods = []
    for food in mitigating_foods:
        for drug in medications:
            from api.services.toxicity_pathway_mappings import get_drug_moa
            drug_moa = get_drug_moa(drug)
            
            if drug_moa == "unknown":
                continue
            
            enhanced = await generate_toxicity_rationale(
                compound=food.get("compound", ""),
                drug_name=drug,
                drug_moa=drug_moa,
                toxicity_pathway=food.get("pathway", ""),
                germline_genes=germline_genes,
                cancer_type=cancer_type,
                treatment_phase=treatment_line,
                base_mechanism=food.get("mechanism", ""),
                timing=food.get("timing", ""),
                dose=food.get("dose", ""),
                provider=provider
            )
            
            enhanced_food = {**food, **enhanced}
            enhanced_foods.append(enhanced_food)
    
    dossier["food_recommendations"] = enhanced_foods
    
    # Generate executive summary if LLM available
    if LLM_AVAILABLE and enhanced_foods:
        try:
            summary_prompt = f"""Summarize the following toxicity mitigation plan in 3-4 sentences:

Patient: {cancer_type}, on {', '.join(medications)}
Germline variants: {', '.join(germline_genes) if germline_genes else 'none'}
Recommended foods: {', '.join([f['compound'] for f in enhanced_foods[:3]])}

Focus on: What toxicities are being addressed and how the foods help."""

            summary = query_llm(summary_prompt, provider=provider)
            if summary and not summary.startswith("Error"):
                dossier["executive_summary"] = summary.strip()
                dossier["llm_enhanced"] = True
        except Exception as e:
            print(f"Executive summary generation failed: {e}")
    
    return dossier


# Singleton service
_llm_toxicity_service = None

def get_llm_toxicity_service():
    global _llm_toxicity_service
    if _llm_toxicity_service is None:
        _llm_toxicity_service = {
            "generate_rationale": generate_toxicity_rationale,
            "generate_dossier": generate_mitigation_dossier,
            "available": LLM_AVAILABLE
        }
    return _llm_toxicity_service
```

**Test Command**:

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 -c "
import asyncio
import sys
from pathlib import Path

# Ensure .env is loaded
try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).parent.parent.parent.parent / '.env')
except: pass

from api.services.llm_toxicity_service import generate_toxicity_rationale, LLM_AVAILABLE

async def test():
    print(f'LLM Available: {LLM_AVAILABLE}')
    
    result = await generate_toxicity_rationale(
        compound='NAC (N-Acetyl Cysteine)',
        drug_name='carboplatin',
        drug_moa='platinum_agent',
        toxicity_pathway='dna_repair',
        germline_genes=['BRCA1'],
        cancer_type='ovarian cancer',
        treatment_phase='first-line chemotherapy',
        base_mechanism='Glutathione precursor, supports DNA repair enzymes',
        timing='post-chemo (not during infusion)',
        dose='600mg twice daily',
        provider='gemini'
    )
    
    print(f'\\nRationale (LLM Enhanced: {result[\"llm_enhanced\"]}):')
    print(result['rationale'])
    print(f'\\nPatient Summary:')
    print(result['patient_summary'])

asyncio.run(test())
"
```

---

### TASK 10: Integrate LLM Rationales into Food Validation

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

**Action**: After finding `toxicity_mitigation`, optionally enhance with LLM rationale

Add import at top:
```python
from api.services.llm_toxicity_service import get_llm_toxicity_service
```

Add after the `toxicity_mitigation` block (around line ~950):
```python
        # === LLM ENHANCEMENT (Optional) ===
        if toxicity_mitigation and request.enable_llm_enhancement:
            try:
                llm_service = get_llm_toxicity_service()
                if llm_service["available"]:
                    germline_genes = []
                    if disease_context:
                        for mut in disease_context.get("mutations", []):
                            if isinstance(mut, dict) and mut.get("gene"):
                                germline_genes.append(mut["gene"])
                    
                    enhanced = await llm_service["generate_rationale"](
                        compound=compound,
                        drug_name=toxicity_mitigation["target_drug"],
                        drug_moa=toxicity_mitigation["target_moa"],
                        toxicity_pathway=toxicity_mitigation["pathway"],
                        germline_genes=germline_genes,
                        cancer_type=disease_context.get("disease", "cancer") if disease_context else "cancer",
                        base_mechanism=toxicity_mitigation["mechanism"],
                        timing=toxicity_mitigation["timing"],
                        dose=toxicity_mitigation.get("dose", ""),
                        provider="gemini"
                    )
                    
                    toxicity_mitigation["llm_rationale"] = enhanced.get("rationale")
                    toxicity_mitigation["patient_summary"] = enhanced.get("patient_summary")
                    toxicity_mitigation["llm_enhanced"] = enhanced.get("llm_enhanced", False)
            except Exception as e:
                logger.warning(f"LLM enhancement failed: {e}")
```

**Request Schema Update**: Add `enable_llm_enhancement: bool = False` to request model

---

### TASK 11: Create LLM Enhancement Test

**File**: Create `oncology-coPilot/oncology-backend-minimal/test_llm_toxicity.py`

```python
"""
Test LLM Enhancement for Toxicity Mitigation
"""

import asyncio
import sys
from pathlib import Path

# Load .env
try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).parent.parent.parent / '.env')
except:
    pass

sys.path.insert(0, str(Path(__file__).parent))


async def test_llm_enhancement():
    print("=" * 60)
    print("LLM TOXICITY ENHANCEMENT TEST")
    print("=" * 60)
    
    from api.services.llm_toxicity_service import (
        generate_toxicity_rationale, 
        generate_mitigation_dossier,
        LLM_AVAILABLE
    )
    
    print(f"\n[CHECK] LLM Available: {LLM_AVAILABLE}")
    
    if not LLM_AVAILABLE:
        print("⚠️  LLM not available - using static rationales")
        print("    To enable: Set GEMINI_API_KEY in .env")
    
    # Test 1: Single rationale generation
    print("\n[TEST 1] Single Rationale Generation")
    print("-" * 40)
    
    result = await generate_toxicity_rationale(
        compound="NAC (N-Acetyl Cysteine)",
        drug_name="carboplatin",
        drug_moa="platinum_agent",
        toxicity_pathway="dna_repair",
        germline_genes=["BRCA1"],
        cancer_type="ovarian cancer (HGSOC)",
        treatment_phase="first-line chemotherapy",
        base_mechanism="Glutathione precursor, supports DNA repair enzymes",
        timing="post-chemo (not during infusion)",
        dose="600mg twice daily"
    )
    
    print(f"   LLM Enhanced: {result['llm_enhanced']}")
    print(f"   Confidence: {result['confidence']}")
    print(f"\n   Rationale:")
    print(f"   {result['rationale']}")
    print(f"\n   Patient Summary:")
    print(f"   {result['patient_summary']}")
    
    # Test 2: Cardiotoxicity scenario
    print("\n[TEST 2] Cardiotoxicity Mitigation")
    print("-" * 40)
    
    result2 = await generate_toxicity_rationale(
        compound="CoQ10 (Ubiquinol)",
        drug_name="doxorubicin",
        drug_moa="anthracycline",
        toxicity_pathway="cardiometabolic",
        germline_genes=["KCNQ1"],
        cancer_type="breast cancer",
        treatment_phase="adjuvant chemotherapy",
        base_mechanism="Mitochondrial support, cardioprotective",
        timing="with fatty meal, continuous during treatment",
        dose="200-400mg daily"
    )
    
    print(f"   LLM Enhanced: {result2['llm_enhanced']}")
    print(f"\n   Rationale: {result2['rationale'][:200]}...")
    
    # Test 3: Full dossier generation
    print("\n[TEST 3] Full Mitigation Dossier")
    print("-" * 40)
    
    from api.services.toxicity_pathway_mappings import get_mitigating_foods
    
    mitigating_foods = get_mitigating_foods({"dna_repair": 1.0})
    
    dossier = await generate_mitigation_dossier(
        patient_context={
            "cancer_type": "ovarian cancer",
            "germline_genes": ["BRCA1"],
            "treatment_line": "first-line"
        },
        medications=["carboplatin"],
        mitigating_foods=mitigating_foods
    )
    
    print(f"   Dossier Generated: {bool(dossier)}")
    print(f"   LLM Enhanced: {dossier.get('llm_enhanced', False)}")
    print(f"   Food Recommendations: {len(dossier.get('food_recommendations', []))}")
    
    if dossier.get("executive_summary"):
        print(f"\n   Executive Summary:")
        print(f"   {dossier['executive_summary'][:300]}...")
    
    print("\n" + "=" * 60)
    print("✅ LLM TOXICITY ENHANCEMENT TEST COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(test_llm_enhancement())
```

**Run Test**:
```bash
cd oncology-coPilot/oncology-backend-minimal && python3 test_llm_toxicity.py
```

---

## Phase 3 Task Checklist

| # | Task | File | Status |
|---|------|------|--------|
| 9 | Create `llm_toxicity_service.py` | `api/services/llm_toxicity_service.py` | ⬜ |
| 10 | Integrate LLM rationales into `validate_food_dynamic` | `hypothesis_validator.py` | ⬜ |
| 11 | Create LLM enhancement test | `test_llm_toxicity.py` | ⬜ |

**Execution Order**: 9 → 11 (test) → 10 → API test with `enable_llm_enhancement=true`

---

---

## PHASE 4: TOXICITY NUTRITION DOSSIER 🆕

**Goal**: Generate comprehensive clinical dossiers for toxicity-aware nutrition, similar to `MBD4_TP53_CLINICAL_DOSSIER.md`

**Connection**: Uses existing dossier infrastructure from `ClinicalDossier/` components

---

### TASK 12: Create Toxicity Nutrition Dossier Generator

**File**: Create `oncology-coPilot/oncology-backend-minimal/api/services/toxicity_dossier_generator.py`

**Purpose**: Generate markdown-formatted clinical dossiers for toxicity mitigation

```python
"""
Toxicity Nutrition Dossier Generator

Generates comprehensive clinical dossiers similar to MBD4_TP53_CLINICAL_DOSSIER.md
but focused on toxicity mitigation and nutritional support.
"""

from datetime import datetime
from typing import Dict, List, Any
import uuid

from api.services.toxicity_pathway_mappings import (
    compute_pathway_overlap,
    get_mitigating_foods,
    get_drug_moa
)
from api.services.llm_toxicity_service import generate_mitigation_dossier


async def generate_toxicity_nutrition_dossier(
    patient_context: Dict[str, Any],
    medications: List[str],
    provider: str = "gemini"
) -> Dict[str, Any]:
    """
    Generate a complete toxicity nutrition dossier.
    
    Args:
        patient_context: Dict with:
            - cancer_type: str
            - germline_genes: List[str]
            - treatment_line: str
            - biomarkers: Dict (optional)
        medications: List of drug names
        provider: LLM provider
    
    Returns:
        Dict with dossier data and markdown content
    """
    
    report_id = f"TOX-NUT-{datetime.utcnow().strftime('%Y-%m-%d-%H%M%S')}-{uuid.uuid4().hex[:6]}"
    
    germline_genes = patient_context.get("germline_genes", [])
    cancer_type = patient_context.get("cancer_type", "cancer")
    treatment_line = patient_context.get("treatment_line", "active treatment")
    
    # Compute toxicity pathways for each medication
    toxicity_assessments = []
    all_mitigating_foods = []
    
    for drug in medications:
        moa = get_drug_moa(drug)
        if moa == "unknown":
            continue
        
        pathway_overlap = compute_pathway_overlap(germline_genes, moa)
        mitigating_foods = get_mitigating_foods(pathway_overlap)
        
        toxicity_assessments.append({
            "drug": drug,
            "moa": moa,
            "pathway_overlap": pathway_overlap,
            "risk_pathways": [k for k, v in pathway_overlap.items() if v > 0.3]
        })
        
        all_mitigating_foods.extend(mitigating_foods)
    
    # Remove duplicates
    seen = set()
    unique_foods = []
    for food in all_mitigating_foods:
        if food["compound"] not in seen:
            seen.add(food["compound"])
            unique_foods.append(food)
    
    # Generate LLM-enhanced dossier
    llm_dossier = await generate_mitigation_dossier(
        patient_context=patient_context,
        medications=medications,
        mitigating_foods=unique_foods,
        provider=provider
    )
    
    # Build structured dossier
    dossier = {
        "report_id": report_id,
        "generated_at": datetime.utcnow().isoformat(),
        "patient_context": {
            "cancer_type": cancer_type,
            "treatment_line": treatment_line,
            "germline_variants": germline_genes,
            "current_medications": medications
        },
        "toxicity_assessment": toxicity_assessments,
        "food_recommendations": llm_dossier.get("food_recommendations", unique_foods),
        "executive_summary": llm_dossier.get("executive_summary", ""),
        "timing_protocol": _generate_timing_protocol(unique_foods, medications),
        "monitoring": _generate_monitoring_recommendations(toxicity_assessments),
        "tier": _calculate_dossier_tier(toxicity_assessments, unique_foods),
        "match_score": _calculate_match_score(toxicity_assessments),
        "llm_enhanced": llm_dossier.get("llm_enhanced", False)
    }
    
    # Generate markdown
    dossier["markdown"] = _generate_markdown(dossier)
    
    return dossier


def _generate_timing_protocol(foods: List[Dict], medications: List[str]) -> str:
    """Generate timing protocol for supplements around chemotherapy."""
    
    protocol_lines = [
        "## Supplement Timing Protocol",
        "",
        "### Pre-Chemotherapy (Day Before)",
        "- Ensure adequate hydration",
        "- Take continuous supplements (Vitamin D, Magnesium) as usual",
        "",
        "### Day of Chemotherapy",
        "- **BEFORE infusion**: Continue continuous supplements",
        "- **DURING infusion**: No supplements (avoid interactions)",
        "- **AFTER infusion (4+ hours)**: Resume post-chemo supplements",
        ""
    ]
    
    # Group by timing
    post_chemo = [f for f in foods if "post" in f.get("timing", "").lower()]
    continuous = [f for f in foods if "continuous" in f.get("timing", "").lower()]
    
    if post_chemo:
        protocol_lines.append("### Post-Chemotherapy Supplements")
        for f in post_chemo:
            protocol_lines.append(f"- **{f['compound']}**: {f['dose']} - {f['timing']}")
        protocol_lines.append("")
    
    if continuous:
        protocol_lines.append("### Continuous Supplements")
        for f in continuous:
            protocol_lines.append(f"- **{f['compound']}**: {f['dose']} - {f['timing']}")
        protocol_lines.append("")
    
    return "\n".join(protocol_lines)


def _generate_monitoring_recommendations(assessments: List[Dict]) -> str:
    """Generate monitoring recommendations based on toxicity pathways."""
    
    recommendations = ["## Monitoring Recommendations", ""]
    
    risk_pathways = set()
    for a in assessments:
        risk_pathways.update(a.get("risk_pathways", []))
    
    if "dna_repair" in risk_pathways:
        recommendations.extend([
            "### DNA Repair Pathway Monitoring",
            "- Complete blood count (CBC) weekly during active treatment",
            "- Renal function (creatinine, BUN) before each cycle",
            "- Audiometry baseline and as needed (platinum agents)",
            ""
        ])
    
    if "cardiometabolic" in risk_pathways:
        recommendations.extend([
            "### Cardiovascular Monitoring",
            "- Baseline echocardiogram before anthracycline",
            "- LVEF assessment every 3 cycles",
            "- ECG monitoring for QT prolongation",
            "- Troponin if symptomatic",
            ""
        ])
    
    if "inflammation" in risk_pathways:
        recommendations.extend([
            "### Immune/Inflammation Monitoring",
            "- Thyroid function (TSH, T4) before each cycle",
            "- Liver function tests (AST, ALT, bilirubin)",
            "- Watch for iRAE symptoms: rash, diarrhea, fatigue",
            ""
        ])
    
    return "\n".join(recommendations)


def _calculate_dossier_tier(assessments: List[Dict], foods: List[Dict]) -> str:
    """Calculate dossier tier (TOP_TIER, GOOD_TIER, or ACCEPTABLE)."""
    
    if not assessments:
        return "ACCEPTABLE"
    
    # High-risk pathways detected + matching foods = TOP_TIER
    risk_pathways = set()
    for a in assessments:
        risk_pathways.update(a.get("risk_pathways", []))
    
    food_pathways = set(f.get("pathway") for f in foods)
    
    overlap = risk_pathways & food_pathways
    
    if len(overlap) >= 2 and len(foods) >= 3:
        return "TOP_TIER"
    elif len(overlap) >= 1 and len(foods) >= 2:
        return "GOOD_TIER"
    else:
        return "ACCEPTABLE"


def _calculate_match_score(assessments: List[Dict]) -> float:
    """Calculate overall match score for the dossier."""
    
    if not assessments:
        return 0.5
    
    # Average of pathway overlaps
    total = 0
    count = 0
    for a in assessments:
        for pathway, score in a.get("pathway_overlap", {}).items():
            total += score
            count += 1
    
    return round(total / max(count, 1), 2)


def _generate_markdown(dossier: Dict) -> str:
    """Generate full markdown dossier."""
    
    md = f"""# Toxicity Nutrition Dossier
## {dossier['patient_context']['cancer_type'].replace('_', ' ').title()}

**Report ID**: {dossier['report_id']}  
**Generated**: {dossier['generated_at']}  
**Tier**: {dossier['tier']}  
**Match Score**: {int(dossier['match_score'] * 100)}%

---

## Executive Summary

{dossier.get('executive_summary', 'Toxicity-aware nutritional support recommendations based on current medications and germline profile.')}

---

## Patient Context

| Parameter | Value |
|-----------|-------|
| Cancer Type | {dossier['patient_context']['cancer_type']} |
| Treatment Line | {dossier['patient_context']['treatment_line']} |
| Current Medications | {', '.join(dossier['patient_context']['current_medications'])} |
| Germline Variants | {', '.join(dossier['patient_context']['germline_variants']) or 'None identified'} |

---

## Toxicity Assessment

"""
    
    for assessment in dossier['toxicity_assessment']:
        md += f"""### {assessment['drug'].title()} ({assessment['moa']})

**Risk Pathways**: {', '.join(assessment['risk_pathways']) or 'Low risk'}

| Pathway | Overlap Score |
|---------|--------------|
"""
        for pathway, score in assessment['pathway_overlap'].items():
            md += f"| {pathway} | {score:.2f} |\n"
        md += "\n"
    
    md += """---

## Recommended Nutritional Support

"""
    
    for i, food in enumerate(dossier['food_recommendations'], 1):
        md += f"""### {i}. {food['compound']}

- **Dose**: {food.get('dose', 'Per healthcare provider')}
- **Timing**: {food.get('timing', 'As directed')}
- **Mechanism**: {food.get('mechanism', '')}
- **Evidence Tier**: {food.get('evidence_tier', 'MODERATE')}
- **Target Pathway**: {food.get('pathway', '')}

"""
        if food.get('rationale'):
            md += f"**Rationale**: {food['rationale']}\n\n"
        if food.get('patient_summary'):
            md += f"**For the Patient**: {food['patient_summary']}\n\n"
    
    md += dossier.get('timing_protocol', '')
    md += "\n\n"
    md += dossier.get('monitoring', '')
    
    md += """

---

## Disclaimer

This dossier provides **mechanism-aligned nutritional support recommendations**. It does NOT predict clinical outcomes or replace professional medical advice. All supplement decisions should be discussed with your oncology team.

---

*Generated by Precision Oncology Platform*
"""
    
    return md
```

---

### TASK 13: Create Dossier API Endpoint

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

**Action**: Add new endpoint for generating toxicity nutrition dossiers

```python
@router.post("/toxicity_nutrition_dossier")
async def generate_toxicity_nutrition_dossier_endpoint(
    request: ToxicityNutritionDossierRequest
):
    """
    Generate a comprehensive toxicity nutrition dossier.
    
    Returns structured dossier data + markdown for display.
    """
    from api.services.toxicity_dossier_generator import generate_toxicity_nutrition_dossier
    
    dossier = await generate_toxicity_nutrition_dossier(
        patient_context={
            "cancer_type": request.cancer_type,
            "germline_genes": request.germline_genes,
            "treatment_line": request.treatment_line,
            "biomarkers": request.biomarkers
        },
        medications=request.medications,
        provider=request.llm_provider
    )
    
    return dossier
```

**Request Schema** (add to schemas):
```python
class ToxicityNutritionDossierRequest(BaseModel):
    cancer_type: str
    germline_genes: List[str] = []
    treatment_line: str = "active treatment"
    medications: List[str]
    biomarkers: Dict[str, Any] = {}
    llm_provider: str = "gemini"
```

---

## Phase 4 Task Checklist

| # | Task | File | Status |
|---|------|------|--------|
| 12 | Create `toxicity_dossier_generator.py` | `api/services/toxicity_dossier_generator.py` | ⬜ |
| 13 | Create dossier API endpoint | `hypothesis_validator.py` | ⬜ |
| 14 | Create dossier generation test | `test_dossier_generator.py` | ⬜ |

---

---

## PHASE 5: FRONTEND INTEGRATION 🆕

**Goal**: Display toxicity nutrition dossiers in existing frontend infrastructure

**Existing Components to Reuse**:
- `DossierSummaryCard.jsx` - Already displays dossier metadata (tier, match_score, title)
- `ClinicalDossierView.jsx` - Container for dossier sections
- `useClinicalDossier.js` - API hook pattern

---

### TASK 15: Create ToxicityNutritionDossierView Component

**File**: Create `oncology-coPilot/oncology-frontend/src/components/ayesha/ToxicityNutritionDossierView.jsx`

**Purpose**: Display toxicity nutrition dossiers using existing design patterns

```jsx
/**
 * ToxicityNutritionDossierView
 * 
 * Displays a comprehensive toxicity nutrition dossier.
 * Reuses patterns from ClinicalDossierView.
 */
import React, { useState, useEffect } from 'react';
import {
  Box,
  Card,
  CardContent,
  Typography,
  Chip,
  LinearProgress,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Alert,
  Button,
  Skeleton
} from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import SecurityIcon from '@mui/icons-material/Security';
import RestaurantIcon from '@mui/icons-material/Restaurant';
import WarningIcon from '@mui/icons-material/Warning';
import LocalHospitalIcon from '@mui/icons-material/LocalHospital';
import DownloadIcon from '@mui/icons-material/Download';

import { apiPost } from '../ClinicalGenomicsCommandCenter/apiUtils';

export default function ToxicityNutritionDossierView({
  cancerType,
  germlineGenes = [],
  medications = [],
  treatmentLine = 'active treatment',
  biomarkers = {}
}) {
  const [dossier, setDossier] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    if (medications.length > 0) {
      fetchDossier();
    }
  }, [cancerType, germlineGenes, medications, treatmentLine]);

  const fetchDossier = async () => {
    setLoading(true);
    setError(null);
    
    try {
      const response = await apiPost('/api/hypothesis/toxicity_nutrition_dossier', {
        cancer_type: cancerType,
        germline_genes: germlineGenes,
        treatment_line: treatmentLine,
        medications: medications,
        biomarkers: biomarkers,
        llm_provider: 'gemini'
      });
      
      setDossier(response);
    } catch (err) {
      setError(err.message || 'Failed to generate dossier');
    } finally {
      setLoading(false);
    }
  };

  const getTierColor = (tier) => {
    switch(tier) {
      case 'TOP_TIER': return 'success';
      case 'GOOD_TIER': return 'info';
      default: return 'default';
    }
  };

  const handleExport = () => {
    if (!dossier?.markdown) return;
    
    const blob = new Blob([dossier.markdown], { type: 'text/markdown' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `toxicity_nutrition_dossier_${dossier.report_id}.md`;
    a.click();
    URL.revokeObjectURL(url);
  };

  if (loading) {
    return (
      <Card sx={{ p: 3 }}>
        <Skeleton variant="text" width="60%" height={40} />
        <Skeleton variant="rectangular" height={200} sx={{ mt: 2 }} />
        <Skeleton variant="rectangular" height={150} sx={{ mt: 2 }} />
      </Card>
    );
  }

  if (error) {
    return (
      <Alert severity="error" action={
        <Button color="inherit" size="small" onClick={fetchDossier}>
          Retry
        </Button>
      }>
        {error}
      </Alert>
    );
  }

  if (!dossier) {
    return (
      <Card sx={{ p: 3 }}>
        <Typography color="text.secondary">
          Select medications to generate toxicity nutrition dossier.
        </Typography>
      </Card>
    );
  }

  return (
    <Card sx={{ p: 3 }}>
      {/* Header */}
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 3 }}>
        <SecurityIcon color="success" fontSize="large" />
        <Box sx={{ flex: 1 }}>
          <Typography variant="h5" fontWeight="bold">
            Toxicity Nutrition Dossier
          </Typography>
          <Typography variant="body2" color="text.secondary">
            Report ID: {dossier.report_id}
          </Typography>
        </Box>
        <Chip 
          label={dossier.tier?.replace('_', ' ')} 
          color={getTierColor(dossier.tier)}
        />
        {dossier.llm_enhanced && (
          <Chip label="LLM Enhanced" variant="outlined" color="secondary" size="small" />
        )}
        <Button 
          startIcon={<DownloadIcon />} 
          variant="outlined" 
          onClick={handleExport}
        >
          Export
        </Button>
      </Box>

      {/* Disclaimer */}
      <Alert severity="warning" sx={{ mb: 3 }} icon={<WarningIcon />}>
        <Typography variant="subtitle2" fontWeight="bold">
          Research Use Only - Not Medical Advice
        </Typography>
        <Typography variant="body2">
          These are mechanism-aligned recommendations. Consult your oncology team before making any dietary changes.
        </Typography>
      </Alert>

      {/* Executive Summary */}
      {dossier.executive_summary && (
        <Card variant="outlined" sx={{ mb: 3, bgcolor: 'success.50' }}>
          <CardContent>
            <Typography variant="h6" gutterBottom>
              Executive Summary
            </Typography>
            <Typography variant="body1">
              {dossier.executive_summary}
            </Typography>
          </CardContent>
        </Card>
      )}

      {/* Match Score */}
      <Box sx={{ mb: 3 }}>
        <Typography variant="body2" color="text.secondary">
          Pathway Match Score
        </Typography>
        <LinearProgress 
          variant="determinate" 
          value={dossier.match_score * 100}
          color="success"
          sx={{ height: 10, borderRadius: 1 }}
        />
        <Typography variant="h6" color="success.main">
          {Math.round(dossier.match_score * 100)}%
        </Typography>
      </Box>

      {/* Toxicity Assessment */}
      <Accordion defaultExpanded>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <LocalHospitalIcon sx={{ mr: 1 }} />
          <Typography variant="h6">Toxicity Assessment</Typography>
        </AccordionSummary>
        <AccordionDetails>
          {dossier.toxicity_assessment?.map((assessment, idx) => (
            <Box key={idx} sx={{ mb: 2 }}>
              <Typography variant="subtitle1" fontWeight="bold">
                {assessment.drug} ({assessment.moa})
              </Typography>
              <Box sx={{ display: 'flex', gap: 1, mt: 1, flexWrap: 'wrap' }}>
                {assessment.risk_pathways?.map(pathway => (
                  <Chip 
                    key={pathway} 
                    label={pathway} 
                    color="warning" 
                    size="small" 
                  />
                ))}
              </Box>
            </Box>
          ))}
        </AccordionDetails>
      </Accordion>

      {/* Food Recommendations */}
      <Accordion defaultExpanded>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <RestaurantIcon sx={{ mr: 1 }} />
          <Typography variant="h6">
            Recommended Nutritional Support ({dossier.food_recommendations?.length || 0})
          </Typography>
        </AccordionSummary>
        <AccordionDetails>
          {dossier.food_recommendations?.map((food, idx) => (
            <Card key={idx} variant="outlined" sx={{ mb: 2, p: 2 }}>
              <Typography variant="subtitle1" fontWeight="bold">
                {idx + 1}. {food.compound}
              </Typography>
              <Box sx={{ display: 'flex', gap: 1, mt: 1, mb: 1 }}>
                <Chip label={food.pathway} color="success" size="small" />
                <Chip label={food.evidence_tier} variant="outlined" size="small" />
              </Box>
              <Typography variant="body2" color="text.secondary">
                <strong>Dose:</strong> {food.dose}
              </Typography>
              <Typography variant="body2" color="text.secondary">
                <strong>Timing:</strong> {food.timing}
              </Typography>
              {food.rationale && (
                <Typography variant="body2" sx={{ mt: 1 }}>
                  {food.rationale}
                </Typography>
              )}
              {food.patient_summary && (
                <Alert severity="info" sx={{ mt: 1 }}>
                  <Typography variant="body2">
                    <strong>For the Patient:</strong> {food.patient_summary}
                  </Typography>
                </Alert>
              )}
            </Card>
          ))}
        </AccordionDetails>
      </Accordion>

      {/* Timing Protocol */}
      {dossier.timing_protocol && (
        <Accordion>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="h6">Timing Protocol</Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Typography 
              component="pre" 
              sx={{ 
                whiteSpace: 'pre-wrap', 
                fontFamily: 'inherit',
                fontSize: '0.9rem'
              }}
            >
              {dossier.timing_protocol}
            </Typography>
          </AccordionDetails>
        </Accordion>
      )}

      {/* Monitoring */}
      {dossier.monitoring && (
        <Accordion>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="h6">Monitoring Recommendations</Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Typography 
              component="pre" 
              sx={{ 
                whiteSpace: 'pre-wrap', 
                fontFamily: 'inherit',
                fontSize: '0.9rem'
              }}
            >
              {dossier.monitoring}
            </Typography>
          </AccordionDetails>
        </Accordion>
      )}
    </Card>
  );
}
```

---

### TASK 16: Integrate with Existing FoodRankingPanel

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/FoodRankingPanel.jsx`

**Action**: Add "Generate Full Dossier" button that opens ToxicityNutritionDossierView

Add after the disclaimer section:
```jsx
{/* Generate Full Dossier Button */}
{foods.some(f => f.toxicity_mitigation?.mitigates) && (
  <Box sx={{ mb: 3 }}>
    <Button
      variant="contained"
      color="success"
      startIcon={<SecurityIcon />}
      onClick={() => setShowDossier(true)}
    >
      Generate Full Toxicity Nutrition Dossier
    </Button>
  </Box>
)}
```

---

## Phase 5 Task Checklist

| # | Task | File | Status |
|---|------|------|--------|
| 15 | Create `ToxicityNutritionDossierView.jsx` | Frontend | ⬜ |
| 16 | Add "Generate Dossier" button to FoodRankingPanel | Frontend | ⬜ |
| 17 | Update routes for dossier view | Frontend | ⬜ |
| 18 | Create integration test | Frontend | ⬜ |

---

---

## 📋 COMPLETE IMPLEMENTATION ROADMAP

| Phase | Goal | Tasks | Effort |
|-------|------|-------|--------|
| **Phase 1** ✅ | Backend MOAT | Tasks 1-5 | COMPLETE |
| **Phase 2** ✅ | Food Validation Integration | Tasks 6-8 | COMPLETE |
| **Phase 3** 🆕 | LLM Enhancement | Tasks 9-11 | 2-3 hours |
| **Phase 4** 🆕 | Dossier Generation | Tasks 12-14 | 2-3 hours |
| **Phase 5** 🆕 | Frontend Integration | Tasks 15-18 | 3-4 hours |

---

## 🚀 AGENT JR: START WITH PHASE 3

**Execute in order**:
1. Task 9: Create `llm_toxicity_service.py`
2. Task 11: Test LLM enhancement (requires GEMINI_API_KEY)
3. Task 10: Integrate into `validate_food_dynamic`
4. Re-run end-to-end tests with `enable_llm_enhancement=true`

**Pre-requisite**: Ensure `GEMINI_API_KEY` is set in `.env`

---

---

## 🧪 COMPREHENSIVE TEST STRATEGY

### FAIL FAST PRINCIPLE

> **Philosophy**: Detect failures at the earliest possible moment. Every phase has pre-flight checks that MUST pass before implementation. If pre-flight fails, STOP and fix before proceeding.

---

## Phase 3 Tests: LLM Enhancement

### Pre-Flight Check (MUST PASS BEFORE Task 9)

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
python3 -c "
print('=== PHASE 3 PRE-FLIGHT CHECK ===')
failures = []

# Test 1: LLM API file exists
import os
llm_path = 'src/tools/llm_api.py'
if os.path.exists(llm_path):
    print('✅ LLM API file exists')
else:
    failures.append('LLM API file missing: src/tools/llm_api.py')
    print('❌ LLM API file MISSING')

# Test 2: Can import query_llm
try:
    import sys
    sys.path.insert(0, '.')
    from src.tools.llm_api import query_llm
    print('✅ query_llm importable')
except ImportError as e:
    failures.append(f'Cannot import query_llm: {e}')
    print(f'❌ Cannot import query_llm: {e}')

# Test 3: GEMINI_API_KEY exists
from dotenv import load_dotenv
load_dotenv()
api_key = os.environ.get('GEMINI_API_KEY', '')
if api_key and api_key != 'your_gemini_api_key_here':
    print(f'✅ GEMINI_API_KEY set ({len(api_key)} chars)')
else:
    failures.append('GEMINI_API_KEY not set or placeholder')
    print('⚠️  GEMINI_API_KEY not set - LLM will fallback to static')

# Test 4: toxicity_pathway_mappings accessible from project root
try:
    sys.path.insert(0, 'oncology-coPilot/oncology-backend-minimal')
    from api.services.toxicity_pathway_mappings import get_mitigating_foods
    print('✅ get_mitigating_foods accessible')
except ImportError as e:
    failures.append(f'Cannot import get_mitigating_foods: {e}')
    print(f'❌ Cannot import get_mitigating_foods: {e}')

print()
if failures:
    print('❌ PRE-FLIGHT FAILED - Fix these issues before proceeding:')
    for f in failures:
        print(f'   - {f}')
    exit(1)
else:
    print('✅ PRE-FLIGHT PASSED - Ready for Phase 3')
"
```

### Edge Cases for Task 9 (llm_toxicity_service.py)

| Edge Case | Input | Expected Behavior | Test |
|-----------|-------|-------------------|------|
| **No API key** | GEMINI_API_KEY missing | Return static rationale, `llm_enhanced=False` | ✅ |
| **API timeout** | Slow response >30s | Fallback to static, log warning | Test with mock |
| **API rate limit** | 429 error | Retry with backoff, then fallback | Test with mock |
| **Empty germline genes** | `germline_genes=[]` | Generate generic rationale (no gene-specific) | ✅ |
| **Unknown drug** | `drug_name="unknown_drug"` | Return static rationale | ✅ |
| **Very long prompt** | 10+ genes, complex context | Truncate prompt, still works | ✅ |
| **Non-English chars** | `cancer_type="卵巢癌"` | Handle gracefully | ✅ |
| **LLM returns error text** | `"Error: ..."` | Detect, fallback to static | ✅ |
| **LLM returns empty** | `""` | Fallback to static | ✅ |
| **LLM returns harmful** | Medical misinformation | Not possible to fully prevent, but prompt constrains output | ⚠️ Monitor |

### Task 9 Test Script

```bash
# Create: oncology-coPilot/oncology-backend-minimal/test_llm_edge_cases.py
cd oncology-coPilot/oncology-backend-minimal
python3 << 'EOF'
"""
LLM Service Edge Case Tests
Tests all edge cases BEFORE integration.
"""

import asyncio
import sys
import os
from pathlib import Path

# Load .env
try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).parent.parent.parent / '.env')
except: pass

sys.path.insert(0, str(Path(__file__).parent))

async def test_edge_cases():
    print("=" * 60)
    print("LLM SERVICE EDGE CASE TESTS")
    print("=" * 60)
    
    from api.services.llm_toxicity_service import (
        generate_toxicity_rationale, 
        LLM_AVAILABLE
    )
    
    print(f"\n[INFO] LLM Available: {LLM_AVAILABLE}")
    
    tests = []
    
    # Test 1: Empty germline genes
    print("\n[TEST 1] Empty germline genes")
    result1 = await generate_toxicity_rationale(
        compound="NAC",
        drug_name="carboplatin",
        drug_moa="platinum_agent",
        toxicity_pathway="dna_repair",
        germline_genes=[],  # EMPTY
        cancer_type="ovarian cancer"
    )
    if result1.get("rationale"):
        print("   ✅ PASS: Returns rationale even without genes")
        tests.append(("Empty genes", True))
    else:
        print("   ❌ FAIL: No rationale returned")
        tests.append(("Empty genes", False))
    
    # Test 2: Unknown drug (should still work)
    print("\n[TEST 2] Unknown drug")
    result2 = await generate_toxicity_rationale(
        compound="NAC",
        drug_name="made_up_drug_xyz",
        drug_moa="unknown",
        toxicity_pathway="dna_repair",
        germline_genes=["BRCA1"],
        cancer_type="cancer"
    )
    if result2.get("rationale"):
        print("   ✅ PASS: Returns rationale for unknown drug")
        tests.append(("Unknown drug", True))
    else:
        print("   ❌ FAIL: No rationale returned")
        tests.append(("Unknown drug", False))
    
    # Test 3: Special characters
    print("\n[TEST 3] Special characters in cancer type")
    result3 = await generate_toxicity_rationale(
        compound="CoQ10 (Ubiquinol)",
        drug_name="doxorubicin",
        drug_moa="anthracycline",
        toxicity_pathway="cardiometabolic",
        germline_genes=["KCNQ1"],
        cancer_type="breast cancer (ER+/HER2-)"
    )
    if result3.get("rationale"):
        print("   ✅ PASS: Handles special chars")
        tests.append(("Special chars", True))
    else:
        print("   ❌ FAIL: No rationale returned")
        tests.append(("Special chars", False))
    
    # Test 4: All pathways (stress test)
    print("\n[TEST 4] All three pathways")
    for pathway in ["dna_repair", "inflammation", "cardiometabolic"]:
        result = await generate_toxicity_rationale(
            compound="Omega-3",
            drug_name="pembrolizumab",
            drug_moa="checkpoint_inhibitor",
            toxicity_pathway=pathway,
            germline_genes=["TNF"],
            cancer_type="melanoma"
        )
        if result.get("rationale"):
            print(f"   ✅ {pathway}: OK")
            tests.append((f"Pathway {pathway}", True))
        else:
            print(f"   ❌ {pathway}: FAILED")
            tests.append((f"Pathway {pathway}", False))
    
    # Summary
    print("\n" + "=" * 60)
    passed = sum(1 for _, r in tests if r)
    total = len(tests)
    print(f"RESULTS: {passed}/{total} tests passed")
    
    if passed == total:
        print("✅ ALL EDGE CASE TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
        for name, result in tests:
            if not result:
                print(f"   FAILED: {name}")
    
    return passed == total

if __name__ == "__main__":
    success = asyncio.run(test_edge_cases())
    exit(0 if success else 1)
EOF
```

### Task 10 Integration Test

```bash
# After Task 10: Test LLM integration in validate_food_dynamic
cd oncology-coPilot/oncology-backend-minimal
curl -X POST http://127.0.0.1:8000/api/hypothesis/validate_food_dynamic \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "NAC",
    "disease_context": {
      "disease": "ovarian_cancer",
      "mutations": [{"gene": "BRCA1", "variant": "pathogenic"}]
    },
    "patient_medications": ["carboplatin"],
    "enable_llm_enhancement": true
  }' | jq '.toxicity_mitigation'
```

**Expected Output**:
```json
{
  "mitigates": true,
  "target_drug": "carboplatin",
  "pathway": "dna_repair",
  "mechanism": "Glutathione precursor...",
  "llm_rationale": "NAC supports DNA repair by boosting glutathione levels, which are depleted by carboplatin's platinum-based mechanism. For BRCA1 carriers, this is particularly important as...",
  "patient_summary": "NAC is a supplement that may help protect your cells during carboplatin treatment...",
  "llm_enhanced": true
}
```

---

## Phase 4 Tests: Dossier Generation

### Pre-Flight Check (MUST PASS BEFORE Task 12)

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 -c "
print('=== PHASE 4 PRE-FLIGHT CHECK ===')
failures = []

# Test 1: Phase 3 complete (LLM service exists)
try:
    from api.services.llm_toxicity_service import generate_toxicity_rationale
    print('✅ llm_toxicity_service exists')
except ImportError as e:
    failures.append(f'Phase 3 not complete: {e}')
    print(f'❌ llm_toxicity_service missing - complete Phase 3 first!')

# Test 2: All toxicity functions work
try:
    from api.services.toxicity_pathway_mappings import (
        compute_pathway_overlap, get_mitigating_foods, get_drug_moa
    )
    overlap = compute_pathway_overlap(['BRCA1'], 'platinum_agent')
    foods = get_mitigating_foods(overlap)
    moa = get_drug_moa('carboplatin')
    assert len(foods) >= 3
    assert moa == 'platinum_agent'
    print('✅ Toxicity functions work')
except Exception as e:
    failures.append(f'Toxicity functions failed: {e}')
    print(f'❌ Toxicity functions failed: {e}')

# Test 3: Safety service returns mitigating_foods
try:
    import asyncio
    from api.services.safety_service import get_safety_service
    from api.schemas.safety import (
        ToxicityRiskRequest, PatientContext, GermlineVariant,
        TherapeuticCandidate, ClinicalContext
    )
    
    async def check():
        service = get_safety_service()
        request = ToxicityRiskRequest(
            patient=PatientContext(germlineVariants=[
                GermlineVariant(chrom='17', pos=41276045, ref='A', alt='G', gene='BRCA1')
            ]),
            candidate=TherapeuticCandidate(type='drug', moa='platinum_agent'),
            context=ClinicalContext(disease='ovarian_cancer'),
            options={'profile': 'baseline'}
        )
        result = await service.compute_toxicity_risk(request)
        return hasattr(result, 'mitigating_foods') and len(result.mitigating_foods) >= 3
    
    if asyncio.run(check()):
        print('✅ Safety service returns mitigating_foods')
    else:
        failures.append('Safety service missing mitigating_foods')
        print('❌ Safety service missing mitigating_foods')
except Exception as e:
    failures.append(f'Safety service check failed: {e}')
    print(f'❌ Safety service check failed: {e}')

print()
if failures:
    print('❌ PRE-FLIGHT FAILED:')
    for f in failures:
        print(f'   - {f}')
    exit(1)
else:
    print('✅ PRE-FLIGHT PASSED - Ready for Phase 4')
"
```

### Edge Cases for Task 12 (toxicity_dossier_generator.py)

| Edge Case | Input | Expected Behavior | Test |
|-----------|-------|-------------------|------|
| **No medications** | `medications=[]` | Return empty dossier, tier="ACCEPTABLE" | ✅ |
| **Unknown medications** | `["magic_drug"]` | Skip unknown, return partial dossier | ✅ |
| **No germline genes** | `germline_genes=[]` | Return dossier with no gene-specific foods | ✅ |
| **Multiple medications** | 3+ drugs | Combine assessments, dedupe foods | ✅ |
| **Mixed known/unknown** | `["carboplatin", "unknown_drug"]` | Process known, skip unknown | ✅ |
| **LLM fails mid-dossier** | API error | Continue with static content | ✅ |
| **Very long gene list** | 20+ genes | Handle gracefully, may truncate | ✅ |

### Task 12 Test Script

```bash
# Create: oncology-coPilot/oncology-backend-minimal/test_dossier_edge_cases.py
cd oncology-coPilot/oncology-backend-minimal
python3 << 'EOF'
"""
Dossier Generator Edge Case Tests
"""

import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

async def test_dossier_edge_cases():
    print("=" * 60)
    print("DOSSIER GENERATOR EDGE CASE TESTS")
    print("=" * 60)
    
    from api.services.toxicity_dossier_generator import generate_toxicity_nutrition_dossier
    
    tests = []
    
    # Test 1: No medications
    print("\n[TEST 1] No medications")
    dossier1 = await generate_toxicity_nutrition_dossier(
        patient_context={"cancer_type": "ovarian_cancer", "germline_genes": ["BRCA1"]},
        medications=[]
    )
    if dossier1.get("tier") == "ACCEPTABLE" and dossier1.get("food_recommendations") == []:
        print("   ✅ PASS: Empty dossier for no medications")
        tests.append(("No medications", True))
    else:
        print(f"   ❌ FAIL: tier={dossier1.get('tier')}, foods={len(dossier1.get('food_recommendations', []))}")
        tests.append(("No medications", False))
    
    # Test 2: Unknown medication
    print("\n[TEST 2] Unknown medication")
    dossier2 = await generate_toxicity_nutrition_dossier(
        patient_context={"cancer_type": "cancer", "germline_genes": []},
        medications=["magic_potion_xyz"]
    )
    if dossier2.get("toxicity_assessment") == []:
        print("   ✅ PASS: Skips unknown medications")
        tests.append(("Unknown medication", True))
    else:
        print(f"   ❌ FAIL: Should have empty assessment")
        tests.append(("Unknown medication", False))
    
    # Test 3: Multiple medications
    print("\n[TEST 3] Multiple medications")
    dossier3 = await generate_toxicity_nutrition_dossier(
        patient_context={"cancer_type": "ovarian_cancer", "germline_genes": ["BRCA1", "TP53"]},
        medications=["carboplatin", "doxorubicin", "pembrolizumab"]
    )
    assessments = len(dossier3.get("toxicity_assessment", []))
    if assessments == 3:
        print(f"   ✅ PASS: {assessments} assessments for 3 drugs")
        tests.append(("Multiple medications", True))
    else:
        print(f"   ❌ FAIL: Expected 3 assessments, got {assessments}")
        tests.append(("Multiple medications", False))
    
    # Test 4: No germline genes
    print("\n[TEST 4] No germline genes")
    dossier4 = await generate_toxicity_nutrition_dossier(
        patient_context={"cancer_type": "breast_cancer", "germline_genes": []},
        medications=["doxorubicin"]
    )
    if dossier4.get("report_id"):
        print("   ✅ PASS: Dossier generated without genes")
        tests.append(("No germline genes", True))
    else:
        print("   ❌ FAIL: No dossier generated")
        tests.append(("No germline genes", False))
    
    # Test 5: Markdown generation
    print("\n[TEST 5] Markdown generation")
    dossier5 = await generate_toxicity_nutrition_dossier(
        patient_context={"cancer_type": "ovarian_cancer", "germline_genes": ["BRCA1"]},
        medications=["carboplatin"]
    )
    markdown = dossier5.get("markdown", "")
    if "# Toxicity Nutrition Dossier" in markdown and "## Recommended Nutritional Support" in markdown:
        print("   ✅ PASS: Markdown has required sections")
        tests.append(("Markdown generation", True))
    else:
        print("   ❌ FAIL: Markdown missing sections")
        tests.append(("Markdown generation", False))
    
    # Summary
    print("\n" + "=" * 60)
    passed = sum(1 for _, r in tests if r)
    total = len(tests)
    print(f"RESULTS: {passed}/{total} tests passed")
    
    if passed == total:
        print("✅ ALL DOSSIER EDGE CASE TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
    
    return passed == total

if __name__ == "__main__":
    success = asyncio.run(test_dossier_edge_cases())
    exit(0 if success else 1)
EOF
```

---

## Phase 5 Tests: Frontend Integration

### Pre-Flight Check (MUST PASS BEFORE Task 15)

```bash
echo "=== PHASE 5 PRE-FLIGHT CHECK ==="

# Test 1: Backend API running
if curl -s http://127.0.0.1:8000/docs > /dev/null; then
    echo "✅ Backend API running"
else
    echo "❌ Backend API not running - start it first"
    exit 1
fi

# Test 2: Dossier endpoint exists
if curl -s -o /dev/null -w "%{http_code}" -X POST http://127.0.0.1:8000/api/hypothesis/toxicity_nutrition_dossier \
   -H "Content-Type: application/json" \
   -d '{"cancer_type": "test", "medications": []}' | grep -q "200\|422"; then
    echo "✅ Dossier endpoint exists"
else
    echo "❌ Dossier endpoint missing - complete Phase 4 first"
    exit 1
fi

# Test 3: Frontend node_modules exist
if [ -d "oncology-coPilot/oncology-frontend/node_modules" ]; then
    echo "✅ Frontend dependencies installed"
else
    echo "❌ Run: cd oncology-coPilot/oncology-frontend && npm install"
    exit 1
fi

# Test 4: FoodRankingPanel exists
if [ -f "oncology-coPilot/oncology-frontend/src/components/ayesha/FoodRankingPanel.jsx" ]; then
    echo "✅ FoodRankingPanel.jsx exists"
else
    echo "❌ FoodRankingPanel.jsx missing"
    exit 1
fi

echo ""
echo "✅ PRE-FLIGHT PASSED - Ready for Phase 5"
```

### Edge Cases for Task 15 (ToxicityNutritionDossierView.jsx)

| Edge Case | Input | Expected Behavior | Test |
|-----------|-------|-------------------|------|
| **No medications prop** | `medications=[]` | Show "Select medications" message | ✅ |
| **API timeout** | Slow backend | Show loading, then error with retry | ✅ |
| **API error** | 500 response | Show error alert with retry button | ✅ |
| **Empty dossier** | No assessments | Show minimal dossier, no crash | ✅ |
| **Missing fields** | Partial response | Handle gracefully with fallbacks | ✅ |
| **Export with special chars** | Unicode in filename | Sanitize filename | ✅ |
| **Very long recommendations** | 20+ foods | Scrollable, performant | ✅ |

### Frontend Test Script

```javascript
// Create: oncology-coPilot/oncology-frontend/src/components/ayesha/__tests__/ToxicityNutritionDossierView.test.jsx
import { render, screen, waitFor, fireEvent } from '@testing-library/react';
import ToxicityNutritionDossierView from '../ToxicityNutritionDossierView';

// Mock apiPost
jest.mock('../../ClinicalGenomicsCommandCenter/apiUtils', () => ({
  apiPost: jest.fn()
}));

import { apiPost } from '../../ClinicalGenomicsCommandCenter/apiUtils';

describe('ToxicityNutritionDossierView', () => {
  it('shows message when no medications provided', () => {
    render(<ToxicityNutritionDossierView medications={[]} />);
    expect(screen.getByText(/Select medications/i)).toBeInTheDocument();
  });

  it('shows loading state during fetch', () => {
    apiPost.mockImplementation(() => new Promise(() => {})); // Never resolves
    render(<ToxicityNutritionDossierView medications={['carboplatin']} cancerType="ovarian" />);
    // Check for skeleton loaders
  });

  it('shows error with retry button on API failure', async () => {
    apiPost.mockRejectedValue(new Error('Network error'));
    render(<ToxicityNutritionDossierView medications={['carboplatin']} cancerType="ovarian" />);
    
    await waitFor(() => {
      expect(screen.getByText(/Failed to generate/i)).toBeInTheDocument();
      expect(screen.getByRole('button', { name: /Retry/i })).toBeInTheDocument();
    });
  });

  it('renders dossier data correctly', async () => {
    apiPost.mockResolvedValue({
      report_id: 'TEST-123',
      tier: 'TOP_TIER',
      match_score: 0.85,
      executive_summary: 'Test summary',
      food_recommendations: [
        { compound: 'NAC', pathway: 'dna_repair', dose: '600mg' }
      ]
    });
    
    render(<ToxicityNutritionDossierView medications={['carboplatin']} cancerType="ovarian" />);
    
    await waitFor(() => {
      expect(screen.getByText('TEST-123')).toBeInTheDocument();
      expect(screen.getByText('TOP TIER')).toBeInTheDocument();
      expect(screen.getByText('NAC')).toBeInTheDocument();
    });
  });
});
```

---

## 🚨 FAIL FAST DETECTION POINTS

### Where Things WILL Fail (Detect Early)

| Failure Point | When Detected | How to Detect | Fix |
|---------------|---------------|---------------|-----|
| **Missing API key** | Phase 3 pre-flight | Check env var | Set GEMINI_API_KEY in .env |
| **LLM import path** | Task 9 test | ImportError | Fix sys.path in llm_toxicity_service |
| **LLM timeout** | Task 9 stress test | 30s+ response | Add timeout handling |
| **Schema mismatch** | Task 13 test | Pydantic validation error | Update schema |
| **Frontend build** | Phase 5 pre-flight | npm error | Fix imports, dependencies |
| **CORS issues** | Frontend API call | Browser console | Add CORS middleware |
| **PropTypes warnings** | Frontend render | Console warnings | Update PropTypes |

### Pre-Commit Checklist

Before each PR/merge:

```bash
#!/bin/bash
# Save as: run_all_tests.sh

echo "=== RUNNING ALL TOXICITY MOAT TESTS ==="

cd oncology-coPilot/oncology-backend-minimal

# Phase 1+2: Core MOAT
echo "[1/5] Phase 1+2: Core MOAT Tests"
bash test_e2e_toxicity_moat.sh || exit 1

# Phase 3: LLM Enhancement
echo "[2/5] Phase 3: LLM Edge Cases"
python3 test_llm_edge_cases.py || exit 1

# Phase 4: Dossier Generation
echo "[3/5] Phase 4: Dossier Edge Cases"
python3 test_dossier_edge_cases.py || exit 1

# Phase 4: API Endpoint
echo "[4/5] Phase 4: Dossier API Test"
curl -s -X POST http://127.0.0.1:8000/api/hypothesis/toxicity_nutrition_dossier \
  -H "Content-Type: application/json" \
  -d '{
    "cancer_type": "ovarian_cancer",
    "germline_genes": ["BRCA1"],
    "medications": ["carboplatin"],
    "treatment_line": "first-line"
  }' | jq -e '.report_id' > /dev/null || exit 1
echo "   ✅ Dossier API returns valid response"

# Phase 5: Frontend
echo "[5/5] Phase 5: Frontend Tests"
cd ../oncology-frontend
npm test -- --watchAll=false || exit 1

echo ""
echo "✅ ALL TESTS PASSED - Safe to merge"
```

---

## 📊 COMPLETE TEST MATRIX

| Phase | Test Type | File | Command | Pass Criteria |
|-------|-----------|------|---------|---------------|
| 1+2 | E2E | `test_e2e_toxicity_moat.sh` | `bash test_e2e_toxicity_moat.sh` | All ✅ |
| 3 | Pre-flight | inline | See above | No failures |
| 3 | Edge cases | `test_llm_edge_cases.py` | `python3 test_llm_edge_cases.py` | All ✅ |
| 3 | Integration | curl | See above | `llm_enhanced: true` |
| 4 | Pre-flight | inline | See above | No failures |
| 4 | Edge cases | `test_dossier_edge_cases.py` | `python3 test_dossier_edge_cases.py` | All ✅ |
| 4 | API | curl | See above | Valid `report_id` |
| 5 | Pre-flight | inline | See above | No failures |
| 5 | Unit | Jest | `npm test` | All ✅ |
| 5 | E2E | Browser | Manual | Renders correctly |

---

## 🎯 AGENT JR EXECUTION ORDER (Updated)

```
PHASE 3: LLM Enhancement
├── 1. Run Phase 3 Pre-Flight Check
│   └── If FAIL → Fix before proceeding
├── 2. Task 9: Create llm_toxicity_service.py
├── 3. Run test_llm_edge_cases.py
│   └── If FAIL → Fix Task 9
├── 4. Task 10: Integrate into hypothesis_validator.py
├── 5. Run integration curl test
│   └── If FAIL → Debug Task 10
└── 6. Mark Phase 3 COMPLETE

PHASE 4: Dossier Generation
├── 1. Run Phase 4 Pre-Flight Check
│   └── If FAIL → Phase 3 incomplete
├── 2. Task 12: Create toxicity_dossier_generator.py
├── 3. Run test_dossier_edge_cases.py
│   └── If FAIL → Fix Task 12
├── 4. Task 13: Create API endpoint
├── 5. Run API curl test
│   └── If FAIL → Debug Task 13
└── 6. Mark Phase 4 COMPLETE

PHASE 5: Frontend Integration
├── 1. Run Phase 5 Pre-Flight Check
│   └── If FAIL → Fix backend/dependencies
├── 2. Task 15: Create ToxicityNutritionDossierView.jsx
├── 3. Run npm test
│   └── If FAIL → Fix component
├── 4. Task 16: Update FoodRankingPanel.jsx
├── 5. Browser test
│   └── If FAIL → Debug
└── 6. Mark Phase 5 COMPLETE
```

---

**Updated**: 2025-01-XX  
**Status**: ✅ Phase 1+2 COMPLETE | 🆕 Phase 3-5 PLANNED (with full test coverage)

