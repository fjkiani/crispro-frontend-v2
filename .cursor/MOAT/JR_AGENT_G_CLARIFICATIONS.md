# 📋 JR Agent G - Implementation Clarifications

**Agent:** JR Agent G (Nutrition Agent)  
**Module:** `AGENT_06_NUTRITION`  
**Date:** January 28, 2025  
**Status:** ✅ CLARIFIED - Ready to Build

---

## ❓ Questions Asked

### Question 1: Scope - Build full agent or MOAT first?
**Answer:** ✅ **Build the full nutrition agent**

**Reasoning:**
- MOAT Tasks (1-7) are **enhancements to existing food validation** system
- Your responsibility is **building `AGENT_06_NUTRITION`** - an orchestrator module
- These are complementary:
  - MOAT = toxicity mitigation features (building blocks)
  - Your Agent = integrates MOAT into orchestrator pipeline

**Evidence:**
- `AGENT_OWNERSHIP_REGISTRY.md:296-315` - You own `AGENT_06_NUTRITION`
- `06_NUTRITION_AGENT.mdc:1-203` - Full module specification
- `orchestrator.py:596-614` - Your integration point (TODO)

---

### Question 2: Call endpoints or use services directly?
**Answer:** ✅ **Use services directly (no HTTP calls)**

**Pattern from existing agents:**

**Biomarker Agent:**
```python
# Direct calculation - no HTTP
tmb = len(mutations) / exome_size
result = {'tmb': {'value': round(tmb, 2), ...}}
```

**Resistance Agent:**
```python
# Import service directly
from ..resistance_prophet_service import get_resistance_prophet_service
prophet = get_resistance_prophet_service()
prediction = await prophet.predict_mm_resistance(...)
```

**Trial Matching Agent:**
```python
# Import agent class
from ..trials import TrialMatchingAgent
agent = TrialMatchingAgent()
result = await agent.match(patient_profile=..., biomarker_profile=...)
```

**Why no HTTP calls:**
- ⚡ Performance (avoid HTTP overhead)
- 🔒 Type safety (IDE support, compile-time checks)
- 🎯 Consistency (all agents follow this pattern)
- 🧪 Testability (no server needed)

---

### Question 3: Should I modify care plan agent?
**Answer:** ❌ **No - care plan auto-consumes your data**

**How it works:**

**You (Nutrition Agent):**
```python
# orchestrator.py:298
state.nutrition_plan = result
```

**Care Plan Agent:**
```python
# orchestrator.py:693-726 & 07_CARE_PLAN_AGENT.mdc:129-131
nutrition_summary = self._summarize_nutrition(state.nutrition_plan)
key_supplements = state.nutrition_plan.supplements[:5] if state.nutrition_plan else []
foods_to_avoid = state.nutrition_plan.foods_to_avoid if state.nutrition_plan else []
```

**Your responsibilities:**
1. ✅ Return dict matching MDC output schema
2. ✅ Wire to orchestrator - populate `state.nutrition_plan`
3. ❌ Don't modify care plan - Senior Agent handles aggregation

---

## 🏗️ Your Build Plan

### What to Build

**1. Core Service**
- **File:** `api/services/nutrition/nutrition_agent.py`
- **Class:** `NutritionAgent`
- **Method:** `async def generate_nutrition_plan(...) -> Dict`

**2. Integrations:**
- `get_mitigating_foods()` - MOAT toxicity mitigation
- `llm_toxicity_service` - LLM rationales (Phase 3 complete)
- `food_spe_integration` - SPE validation
- Drug-food interactions

**3. Orchestrator Wiring**
- **File:** `api/services/orchestrator/orchestrator.py`
- **Method:** `_run_nutrition_agent()` (lines 596-614)
- **Pattern:** Import agent → Build request → Call → Store in `state.nutrition_plan`

**4. Support Files:**
- `api/services/nutrition/__init__.py`
- `api/services/nutrition/drug_food_interactions.py`
- `api/services/nutrition/tests/test_nutrition_agent.py`

### Output Structure (from MDC)

```python
{
    'patient_id': str,
    'treatment': str,
    'supplements': [
        {
            'name': str,
            'dosage': str,
            'timing': str,
            'mechanism': str,
            'evidence_level': str
        }
    ],
    'foods_to_prioritize': [
        {
            'food': str,
            'reason': str,
            'frequency': str,
            'category': str  # "prioritize", "include", "limit"
        }
    ],
    'foods_to_avoid': [
        {
            'food': str,
            'reason': str
        }
    ],
    'drug_food_interactions': [
        {
            'drug': str,
            'food': str,
            'interaction_type': str,  # "avoid", "caution", "timing"
            'mechanism': str,
            'severity': str
        }
    ],
    'timing_rules': {},
    'provenance': {}
}
```

---

## 🎯 Implementation Template

### Nutrition Agent Service

```python
# api/services/nutrition/nutrition_agent.py

from typing import Dict, List, Optional, Any

class NutritionAgent:
    """Generate toxicity-aware nutrition plans."""
    
    def __init__(self):
        """Initialize with required services."""
        from ..toxicity_pathway_mappings import get_mitigating_foods, get_drug_moa
        from ..llm_toxicity_service import get_llm_toxicity_service
        from ..food_spe_integration import FoodSPEIntegrationService
        
        self.get_mitigating_foods = get_mitigating_foods
        self.get_drug_moa = get_drug_moa
        self.llm_service = get_llm_toxicity_service()
        self.spe_service = FoodSPEIntegrationService()
    
    async def generate_nutrition_plan(
        self,
        patient_id: str,
        disease: str,
        mutations: List[Dict],
        current_medications: List[str] = None,
        germline_genes: List[str] = None,
        options: Dict = None
    ) -> Dict:
        """Generate comprehensive nutrition plan."""
        
        # 1. Get toxicity-mitigating foods (MOAT)
        supplements = []
        if current_medications and germline_genes:
            for drug in current_medications:
                moa = self.get_drug_moa(drug)
                if moa != 'unknown':
                    from ..toxicity_pathway_mappings import compute_pathway_overlap
                    pathway_overlap = compute_pathway_overlap(germline_genes, moa)
                    mitigating = self.get_mitigating_foods(pathway_overlap)
                    
                    for food in mitigating:
                        supplements.append({
                            'name': food['compound'],
                            'dosage': food['dose'],
                            'timing': food['timing'],
                            'mechanism': food['mechanism'],
                            'evidence_level': food['evidence_tier'],
                            'target_drug': drug,
                            'target_pathway': food['pathway']
                        })
        
        # 2. Check drug-food interactions
        interactions = await self._check_interactions(current_medications or [])
        
        # 3. Get disease-specific foods
        foods_to_prioritize = self._get_disease_foods(disease)
        
        # 4. Generate timing rules
        timing_rules = self._generate_timing_rules(supplements, current_medications or [])
        
        # 5. Compile foods to avoid
        foods_to_avoid = []
        for interaction in interactions:
            if interaction['severity'] in ['HIGH', 'CRITICAL']:
                foods_to_avoid.append({
                    'food': interaction['food'],
                    'reason': f"Interacts with {interaction['drug']}: {interaction['mechanism']}"
                })
        
        return {
            'patient_id': patient_id,
            'treatment': ', '.join(current_medications) if current_medications else 'unknown',
            'supplements': supplements,
            'foods_to_prioritize': foods_to_prioritize,
            'foods_to_avoid': foods_to_avoid,
            'drug_food_interactions': interactions,
            'timing_rules': timing_rules,
            'provenance': {
                'method': 'toxicity_aware_nutrition',
                'moat_enabled': True,
                'llm_enhanced': self.llm_service['available']
            }
        }
```

### Orchestrator Integration

```python
# api/services/orchestrator/orchestrator.py:596-614

async def _run_nutrition_agent(self, state: PatientState) -> Dict:
    """Run the nutrition planning agent."""
    execution = state.start_agent('nutrition')
    
    try:
        # Import nutrition agent
        from ..nutrition import NutritionAgent
        
        # Extract germline genes
        germline_genes = []
        if state.mutations:
            germline_genes = [m.get('gene') for m in state.mutations if m.get('gene')]
        
        # Extract current medications
        current_medications = []
        if state.patient_profile:
            if isinstance(state.patient_profile, dict):
                current_medications = state.patient_profile.get('current_medications', [])
            else:
                current_medications = getattr(state.patient_profile, 'current_medications', [])
        
        # Generate nutrition plan
        agent = NutritionAgent()
        result = await agent.generate_nutrition_plan(
            patient_id=state.patient_id,
            disease=state.disease or 'unknown',
            mutations=state.mutations or [],
            current_medications=current_medications,
            germline_genes=germline_genes,
            options={'include_llm_rationales': True}
        )
        
        execution.complete(result)
        return result
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

---

## ✅ Success Checklist

- [ ] `NutritionAgent` class created
- [ ] Integrates MOAT (`get_mitigating_foods`)
- [ ] Integrates LLM service (rationales)
- [ ] Returns dict matching MDC schema
- [ ] Wired to orchestrator (replaces TODO)
- [ ] Populates `state.nutrition_plan`
- [ ] Uses direct service imports (no HTTP)
- [ ] Has unit tests
- [ ] Has integration test with orchestrator
- [ ] Registry status updated

---

## 📚 Reference Documents Created

To help other agents with similar questions, I created:

1. **[Agent Implementation Guide](./AGENT_IMPLEMENTATION_GUIDE.md)** (detailed, 500+ lines)
   - Answers all 3 questions in depth
   - Code examples from existing agents
   - Full implementation templates
   - Testing strategies
   - Troubleshooting guide

2. **[Quick Reference](./QUICK_REFERENCE_AGENT_QUESTIONS.md)** (TL;DR version)
   - One-page summary
   - Copy-paste patterns
   - Quick lookup table

3. **Updated Master Index** (added "For New Agents" section)
   - Links to guides at the top
   - TL;DR of key points

4. **Updated Ownership Registry** (added reference section)
   - Links to implementation guides
   - Common questions answered

---

## 🎯 Next Steps

1. ✅ Create `api/services/nutrition/nutrition_agent.py`
2. ✅ Implement `generate_nutrition_plan()` method
3. ✅ Wire to orchestrator (replace TODO at line 596)
4. ✅ Create tests
5. ✅ Update registry status

**Ready to proceed!**

---

**Status:** ✅ ALL QUESTIONS ANSWERED  
**Documentation:** ✅ GUIDES CREATED FOR OTHER AGENTS  
**Next:** Start implementation






