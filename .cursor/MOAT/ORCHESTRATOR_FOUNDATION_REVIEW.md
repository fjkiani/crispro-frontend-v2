# 🏗️ ORCHESTRATOR FOUNDATION - COMPLETE REVIEW

**Purpose:** Review the foundation/scaffolding built by Senior Agent  
**Date:** January 28, 2025  
**For:** JR Agents (especially AGENT_06_NUTRITION)

---

## ✅ WHAT'S BUILT (FOUNDATION)

### 1. Core Orchestrator Infrastructure

**Location:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/`

| File | Purpose | Lines | Status |
|------|---------|-------|--------|
| `orchestrator.py` | Main pipeline coordinator | ~774 | ✅ Complete |
| `state.py` | PatientState dataclass | ~303 | ✅ Complete |
| `state_store.py` | Persistence layer | ~250 | ✅ Complete |
| `message_bus.py` | Inter-agent messaging | ~180 | ✅ Complete |

### 2. Pipeline Architecture

**7-Phase Pipeline:**
```
1. EXTRACTING    → Parse files (VCF, PDF, MAF)
2. ANALYZING     → Biomarker + Resistance + Nutrition (PARALLEL)
3. RANKING       → Drug Efficacy (S/P/E)
4. MATCHING      → Trial Matching
5. PLANNING      → Care Plan Generation
6. MONITORING    → Monitoring Setup
7. COMPLETE      → All done!
```

**Key Pattern:** Agents run in parallel when dependencies allow (line 240-253)

### 3. State Management

**PatientState Structure:**
```python
@dataclass
class PatientState:
    patient_id: str
    disease: Optional[str]
    phase: StatePhase
    
    # Agent outputs (populated as agents complete)
    patient_profile: Optional[Dict]        # From 01_DATA_EXTRACTION
    biomarker_profile: Optional[Dict]      # From 02_BIOMARKER
    resistance_prediction: Optional[Dict]  # From 03_RESISTANCE
    drug_ranking: Optional[List[Dict]]     # From 04_DRUG_EFFICACY
    trial_matches: Optional[List[Dict]]     # From 05_TRIAL_MATCHING
    nutrition_plan: Optional[Dict]         # From 06_NUTRITION ⬅️ MY TARGET
    care_plan: Optional[Dict]               # From 07_CARE_PLAN
    monitoring_config: Optional[Dict]      # From 08_MONITORING
    
    # Execution tracking
    agent_executions: List[AgentExecution]
    alerts: List[Alert]
    history: List[StateChange]
```

### 4. Agent Execution Pattern

**Every agent follows this pattern:**

```python
async def _run_<agent>_agent(self, state: PatientState) -> Dict:
    """Run the <agent> agent."""
    execution = state.start_agent('<agent>')
    
    try:
        # 1. Import service
        from ..<service> import get_<service>
        
        # 2. Build request from state
        mutations = state.mutations
        disease = state.disease
        
        # 3. Call service
        service = get_<service>()
        result = await service.<method>(...)
        
        # 4. Convert to dict if needed
        if is_dataclass(result):
            result = self._dataclass_to_dict(result)
        
        # 5. Mark complete
        execution.complete(result)
        return result
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

### 5. What's Already Wired

| Module | Status | Location in Orchestrator |
|--------|--------|-------------------------|
| **02 Biomarker** | ✅ WIRED | `_run_biomarker_agent()` (line 409) |
| **03 Resistance** | ✅ WIRED | `_run_resistance_agent()` (line 461) |
| **05 Trial Matching** | ✅ WIRED | `_run_trial_matching_agent()` (line 601) |
| **07 Care Plan** | ✅ WIRED | `_run_care_plan_agent()` (line 669) |
| **08 Monitoring** | ✅ WIRED | `_run_monitoring_agent()` (line 704) |

### 6. What Needs Building

| Module | Status | Assigned To |
|--------|--------|-------------|
| **01 Data Extraction** | ⏳ TODO | JR Agent A |
| **04 Drug Efficacy** | ⏳ TODO | JR Agent C |
| **06 Nutrition** | ⏳ TODO | **JR Agent G (ME)** |

---

## 🎯 MY TASK: AGENT_06_NUTRITION

### Current State

**Location:** `orchestrator.py` line 572-590

```python
async def _run_nutrition_agent(self, state: PatientState) -> Dict:
    """Run the nutrition planning agent."""
    execution = state.start_agent('nutrition')
    
    try:
        # TODO: Implement actual nutrition service
        # For now, return placeholder
        result = {
            'recommendations': [],
            'drug_food_interactions': [],
            'timing_rules': []
        }
        
        execution.complete(result)
        return result
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

**Status:** ⏳ SKELETON - Placeholder only

### What I Need to Build

#### 1. Create Nutrition Agent Service

**File:** `api/services/nutrition/nutrition_agent.py`

**Structure:**
```python
"""
Nutrition Agent - AGENT_06_NUTRITION

Provides:
- Food/supplement recommendations
- Drug-food interaction warnings
- Timing rules (when to take supplements)
- Toxicity mitigation foods (Phase 3 LLM work)
"""

from typing import Dict, List, Optional
import logging

from ..toxicity_pathway_mappings import get_mitigating_foods, get_drug_moa
from ..llm_toxicity_service import get_llm_toxicity_service
from ..food_spe_integration import validate_food_dynamic  # Or use service directly

logger = logging.getLogger(__name__)


class NutritionAgent:
    """
    Nutrition planning agent for toxicity-aware food recommendations.
    
    Integrates:
    - Toxicity pathway detection
    - Food validation (SPE framework)
    - LLM-enhanced rationales (Phase 3)
    - Drug-food interactions
    """
    
    def __init__(self):
        self.llm_service = get_llm_toxicity_service()
    
    async def generate_nutrition_plan(
        self,
        mutations: List[Dict],
        disease: str,
        patient_profile: Dict,
        biomarker_profile: Optional[Dict] = None,
        resistance_prediction: Optional[Dict] = None
    ) -> Dict:
        """
        Generate complete nutrition plan for patient.
        
        Args:
            mutations: List of mutation dicts
            disease: Disease type
            patient_profile: Patient context (medications, treatment_line)
            biomarker_profile: TMB, MSI, HRD status
            resistance_prediction: Pathway analyses
        
        Returns:
            Dict with:
            - recommendations: List of food recommendations
            - drug_food_interactions: List of interaction warnings
            - timing_rules: List of timing recommendations
            - toxicity_mitigation: Foods that mitigate drug toxicity
        """
        # Extract context
        medications = patient_profile.get('current_regimen', [])
        if isinstance(medications, str):
            medications = [medications]
        
        germline_genes = [m.get('gene') for m in mutations if m.get('germline', False)]
        
        # Step 1: Get toxicity-mitigating foods
        mitigating_foods = []
        if medications:
            for drug in medications:
                moa = get_drug_moa(drug)
                if moa != "unknown":
                    # Compute pathway overlap (simplified)
                    pathway_overlap = self._compute_pathway_overlap(
                        germline_genes, moa
                    )
                    foods = get_mitigating_foods(pathway_overlap)
                    mitigating_foods.extend(foods)
        
        # Step 2: Validate top foods with SPE framework
        recommendations = []
        for food in mitigating_foods[:10]:  # Top 10
            # Call food validation (or use service directly)
            validated = await self._validate_food(
                compound=food['compound'],
                disease=disease,
                mutations=mutations,
                medications=medications
            )
            if validated:
                recommendations.append(validated)
        
        # Step 3: Check drug-food interactions
        interactions = self._check_interactions(medications, recommendations)
        
        # Step 4: Generate timing rules
        timing_rules = self._generate_timing_rules(recommendations, medications)
        
        return {
            'recommendations': recommendations,
            'drug_food_interactions': interactions,
            'timing_rules': timing_rules,
            'toxicity_mitigation': {
                'foods_count': len(mitigating_foods),
                'pathways_addressed': list(set(f['pathway'] for f in mitigating_foods))
            },
            'provenance': {
                'agent': 'AGENT_06_NUTRITION',
                'method': 'toxicity_aware_nutrition'
            }
        }
    
    def _compute_pathway_overlap(self, genes: List[str], moa: str) -> Dict[str, float]:
        """Compute pathway overlap for toxicity detection."""
        from ..toxicity_pathway_mappings import compute_pathway_overlap
        return compute_pathway_overlap(genes, moa)
    
    async def _validate_food(self, compound: str, disease: str, 
                            mutations: List[Dict], medications: List[str]) -> Optional[Dict]:
        """Validate food using existing food validation service."""
        # Use existing validate_food_dynamic endpoint logic
        # Or call food_spe_integration service directly
        pass
    
    def _check_interactions(self, medications: List[str], 
                           recommendations: List[Dict]) -> List[Dict]:
        """Check for drug-food interactions."""
        # Use existing dietician_recommendations service
        pass
    
    def _generate_timing_rules(self, recommendations: List[Dict],
                          medications: List[str]) -> List[Dict]:
        """Generate timing rules for supplements."""
        # Extract timing from recommendations
        return [
            {
                'compound': r['compound'],
                'timing': r.get('timing', 'continuous'),
                'relative_to': medications[0] if medications else None
            }
            for r in recommendations
        ]


# Singleton pattern
_nutrition_agent = None

def get_nutrition_agent() -> NutritionAgent:
    """Get singleton nutrition agent instance."""
    global _nutrition_agent
    if _nutrition_agent is None:
        _nutrition_agent = NutritionAgent()
    return _nutrition_agent
```

#### 2. Wire into Orchestrator

**Update:** `orchestrator.py` line 572-590

```python
async def _run_nutrition_agent(self, state: PatientState) -> Dict:
    """Run the nutrition planning agent."""
    execution = state.start_agent('nutrition')
    
    try:
        from ..nutrition import get_nutrition_agent
        
        agent = get_nutrition_agent()
        
        result = await agent.generate_nutrition_plan(
            mutations=state.mutations,
            disease=state.disease or 'ovarian',
            patient_profile=state.patient_profile or {},
            biomarker_profile=state.biomarker_profile,
            resistance_prediction=state.resistance_prediction
        )
        
        execution.complete({
            'recommendations_count': len(result.get('recommendations', [])),
            'interactions_count': len(result.get('drug_food_interactions', []))
        })
        return result
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

#### 3. File Structure to Create

```
api/services/nutrition/
├── __init__.py
├── nutrition_agent.py          # Main agent (200+ lines)
├── drug_food_interactions.py   # Interaction checking (100+ lines)
└── tests/
    ├── test_nutrition_agent.py
    └── test_integration.py
```

---

## 📋 INTEGRATION CHECKLIST

### Phase 1: Create Service
- [ ] Create `api/services/nutrition/` directory
- [ ] Create `nutrition_agent.py` with `NutritionAgent` class
- [ ] Implement `generate_nutrition_plan()` method
- [ ] Integrate with existing services:
  - [ ] `toxicity_pathway_mappings.get_mitigating_foods()`
  - [ ] `llm_toxicity_service.get_llm_toxicity_service()`
  - [ ] Food validation (SPE framework)

### Phase 2: Wire to Orchestrator
- [ ] Update `_run_nutrition_agent()` in `orchestrator.py`
- [ ] Import `get_nutrition_agent()`
- [ ] Call agent with state data
- [ ] Handle errors gracefully

### Phase 3: Testing
- [ ] Unit tests for `NutritionAgent`
- [ ] Integration test with orchestrator
- [ ] Test with real patient data (mutations, medications)

### Phase 4: Documentation
- [ ] Update `AGENT_OWNERSHIP_REGISTRY.md` status
- [ ] Document nutrition plan structure
- [ ] Add to care plan generation

---

## 🔗 EXISTING SERVICES TO USE

### Already Built (Can Use Directly)

| Service | Location | What It Does |
|---------|----------|--------------|
| `toxicity_pathway_mappings.py` | `api/services/` | `get_mitigating_foods()` - Maps pathways to foods |
| `llm_toxicity_service.py` | `api/services/` | LLM-enhanced rationales (my Phase 3 work) |
| `food_spe_integration.py` | `api/services/` | SPE validation framework |
| `dietician_recommendations.py` | `api/services/` | Drug-food interaction checking |
| `safety_service.py` | `api/services/` | Toxicity risk assessment |

### Integration Points

1. **Toxicity Detection:**
   ```python
   from ..toxicity_pathway_mappings import get_mitigating_foods, get_drug_moa, compute_pathway_overlap
   
   moa = get_drug_moa(drug_name)
   pathway_overlap = compute_pathway_overlap(germline_genes, moa)
   foods = get_mitigating_foods(pathway_overlap)
   ```

2. **LLM Enhancement:**
   ```python
   from ..llm_toxicity_service import get_llm_toxicity_service
   
   llm_service = get_llm_toxicity_service()
   if llm_service['available']:
       enhanced = await llm_service['generate_rationale'](...)
   ```

3. **Food Validation:**
   ```python
   # Use existing validate_food_dynamic endpoint logic
   # Or call food_spe_integration service directly
   ```

---

## 📊 EXPECTED OUTPUT STRUCTURE

### Nutrition Plan Dict

```python
{
    'recommendations': [
        {
            'compound': 'NAC (N-Acetyl Cysteine)',
            'dose': '600mg twice daily',
            'timing': 'post-chemo (not during infusion)',
            'mechanism': 'Glutathione precursor, supports DNA repair enzymes',
            'pathway': 'dna_repair',
            'evidence_tier': 'MODERATE',
            'toxicity_mitigation': {
                'mitigates': True,
                'target_drug': 'carboplatin',
                'target_moa': 'platinum_agent',
                'llm_rationale': '...',  # From Phase 3
                'patient_summary': '...'  # From Phase 3
            },
            'efficacy_score': 0.72,  # From SPE framework
            'confidence': 0.68
        }
    ],
    'drug_food_interactions': [
        {
            'drug': 'warfarin',
            'food': 'Vitamin D',
            'severity': 'moderate',
            'message': 'Monitor INR closely'
        }
    ],
    'timing_rules': [
        {
            'compound': 'NAC',
            'timing': 'post-infusion',
            'relative_to': 'carboplatin',
            'hours_after': 4
        }
    ],
    'toxicity_mitigation': {
        'foods_count': 3,
        'pathways_addressed': ['dna_repair', 'inflammation']
    },
    'provenance': {
        'agent': 'AGENT_06_NUTRITION',
        'method': 'toxicity_aware_nutrition'
    }
}
```

---

## 🎯 SUCCESS CRITERIA

1. ✅ `NutritionAgent` class created
2. ✅ Wired into orchestrator's `_run_nutrition_agent()`
3. ✅ Returns proper nutrition plan structure
4. ✅ Integrates with existing toxicity services
5. ✅ Uses Phase 3 LLM enhancement
6. ✅ Unit tests passing
7. ✅ Integration test with orchestrator passing
8. ✅ Status updated in `AGENT_OWNERSHIP_REGISTRY.md`

---

## 📝 NOTES

### Key Insights

1. **Parallel Execution:** Nutrition runs in parallel with Biomarker and Resistance (line 250-253)
2. **State Storage:** Result goes into `state.nutrition_plan` (line 274)
3. **Error Handling:** Follow pattern - `execution.fail()` on error
4. **Provenance:** Track method and agent in output
5. **Existing Services:** Lots of code already exists - reuse it!

### Dependencies

- **Requires:** Module 01 (Data Extraction) for mutations
- **Optional:** Module 02 (Biomarker) for TMB/MSI/HRD context
- **Optional:** Module 03 (Resistance) for pathway analyses

### Integration with Care Plan

The care plan agent (Module 07) will consume `state.nutrition_plan` and include it in the final care plan document.

---

**Last Updated:** January 28, 2025  
**Status:** Ready for Implementation  
**Assigned To:** JR Agent G (AGENT_06_NUTRITION)





