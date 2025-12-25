# 🏗️ Agent Implementation Guide

**Purpose:** Clarify common questions for JR Agents building MOAT modules  
**Created:** January 28, 2025  
**Status:** ✅ ACTIVE

---

## 📖 Overview

This guide answers the **3 most common questions** agents face when implementing their assigned modules:

1. **Scope:** What exactly should I build?
2. **Service Integration:** Should I call API endpoints or use services directly?
3. **Orchestrator Integration:** How do I wire my agent into the pipeline?

---

## 🎯 Question 1: Scope - What Should I Build?

### Understanding Your Assignment

Each agent has **TWO potential types of work**:

#### Type A: Build a New Orchestrator Module
- You're creating a **NEW agent service** that the orchestrator calls
- Example modules: 01 (Extraction), 02 (Biomarker), 03 (Resistance), 04 (Efficacy), 05 (Trials), 06 (Nutrition)
- Reference: `.cursor/MOAT/orchestration/XX_MODULE_NAME.mdc`
- Orchestrator file: `api/services/orchestrator/orchestrator.py`

#### Type B: Enhance Existing Features (MOAT/Standalone)
- You're adding features to **existing endpoints** outside orchestration
- Example: Toxicity MOAT (Tasks 1-7 in Food Validation System)
- Reference: Task lists in `.cursor/plans/` or specific plan files

### How to Identify Your Type

**Check the Ownership Registry:**
```
.cursor/MOAT/AGENT_OWNERSHIP_REGISTRY.md
```

Look for your assignment:

| Pattern | Type | Your Task |
|---------|------|-----------|
| `AGENT_XX_NAME` with MDC reference | **Type A** | Build orchestrator module |
| Task list with specific file paths | **Type B** | Enhance existing feature |
| Both listed | **Both** | Build module AND integrate enhancements |

### Example: JR Agent G (Nutrition)

**Assignment:**
- **Type A:** `AGENT_06_NUTRITION` - Build orchestrator module
- **Type B:** MOAT Toxicity Tasks (if assigned) - Enhance food validation

**Relationship:**
- The MOAT features (toxicity mitigation) are **components** your agent uses
- Your module **integrates** these into the orchestrator pipeline
- Build the full agent - it consumes the MOAT capabilities

**Rule of Thumb:**
> If you have an `XX_MODULE.mdc` file assigned, build the **full orchestrator module**. Standalone tasks are building blocks you integrate.

---

## 🔌 Question 2: Service Integration Pattern

### The Golden Rule

**✅ DO:** Import and use service classes/functions directly  
**❌ DON'T:** Make HTTP calls to API endpoints from within the orchestrator

### Why?

1. **Performance:** Direct calls avoid HTTP overhead
2. **Type Safety:** Direct imports provide better IDE support
3. **Consistency:** Matches the established pattern
4. **Simplicity:** No need for serialization/deserialization

### The Pattern (Code Examples)

#### ✅ CORRECT: Direct Service Import

**Example 1: Biomarker Agent** (Inline Logic)
```python
# orchestrator.py:450-483
async def _run_biomarker_agent(self, state: PatientState) -> Dict:
    """Run the biomarker calculation agent."""
    execution = state.start_agent('biomarker')
    
    try:
        mutations = state.mutations
        exome_size = 38.0  # Mb
        
        # Direct calculation (no HTTP call)
        tmb = len(mutations) / exome_size
        
        # Check MSI genes
        dmmr_genes = {'MLH1', 'MSH2', 'MSH6', 'PMS2'}
        mutated_dmmr = [m['gene'] for m in mutations if m.get('gene') in dmmr_genes]
        msi_status = 'MSI-H' if mutated_dmmr else 'MSS'
        
        result = {
            'tmb': {'value': round(tmb, 2), 'classification': 'TMB-H' if tmb >= 10 else 'TMB-L'},
            'msi': {'status': msi_status, 'dmmr_genes_mutated': mutated_dmmr},
            'io_eligible': tmb >= 10 or msi_status == 'MSI-H'
        }
        
        execution.complete(result)
        return result
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

**Example 2: Resistance Agent** (External Service)
```python
# orchestrator.py:485-565
async def _run_resistance_agent(self, state: PatientState) -> Dict:
    """Run the resistance prediction agent."""
    execution = state.start_agent('resistance')
    
    try:
        # Import service directly (no HTTP call)
        from ..resistance_prophet_service import get_resistance_prophet_service
        from ..resistance_playbook_service import get_resistance_playbook_service
        
        prophet = get_resistance_prophet_service()
        playbook = get_resistance_playbook_service()
        
        # Call service methods
        prediction_obj = await prophet.predict_mm_resistance(
            mutations=state.mutations,
            drug_class=profile.get('current_drug_class'),
            treatment_line=profile.get('treatment_line', 1)
        )
        
        # Convert to dict if needed
        prediction = self._dataclass_to_dict(prediction_obj)
        
        execution.complete(prediction)
        return prediction
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

**Example 3: Trial Matching Agent** (Agent Class)
```python
# orchestrator.py:625-691
async def _run_trial_matching_agent(self, state: PatientState) -> Dict:
    """Run the clinical trial matching agent."""
    # Import agent class directly
    from ..trials import TrialMatchingAgent
    
    # Build request from state
    patient_profile = {
        'patient_id': state.patient_id,
        'disease': state.disease,
        'mutations': state.mutations
    }
    
    biomarker_profile = {
        'tmb': getattr(state.biomarker_profile, 'tmb', {}),
        'msi': getattr(state.biomarker_profile, 'msi', {})
    } if state.biomarker_profile else None
    
    # Call agent method
    agent = TrialMatchingAgent()
    result = await agent.match(
        patient_profile=patient_profile,
        biomarker_profile=biomarker_profile,
        max_results=10
    )
    
    # Convert result to dict
    matches = []
    for match in result.matches:
        matches.append({
            'nct_id': match.nct_id,
            'title': match.title,
            'mechanism_fit_score': match.mechanism_fit_score
        })
    
    return {'matches': matches, 'trials_found': result.trials_found}
```

#### ❌ WRONG: HTTP Call from Orchestrator

```python
# DON'T DO THIS
async def _run_biomarker_agent(self, state: PatientState) -> Dict:
    import httpx
    
    # ❌ Making HTTP call to internal API
    async with httpx.AsyncClient() as client:
        response = await client.post(
            "http://localhost:8000/api/biomarkers/calculate",
            json={'mutations': state.mutations}
        )
        return response.json()
```

**Why this is wrong:**
- Adds HTTP overhead
- Requires server to be running
- Creates circular dependency
- Loses type information
- Makes testing harder

### Your Implementation Checklist

- [ ] Import service/agent class from relative path
- [ ] Extract data from `PatientState`
- [ ] Call service methods directly
- [ ] Convert result to dict if needed (use `_dataclass_to_dict` helper)
- [ ] Store result in `state.your_field`
- [ ] Use `execution.complete(result)` pattern
- [ ] Handle errors with `execution.fail(str(e))`

### Service Import Patterns

**Pattern 1: Service Factory Function**
```python
from ..my_service import get_my_service

service = get_my_service()
result = await service.process(...)
```

**Pattern 2: Agent Class**
```python
from ..my_agent import MyAgent

agent = MyAgent()
result = await agent.analyze(...)
```

**Pattern 3: Utility Functions**
```python
from ..my_utils import compute_something, map_something

value = compute_something(data)
mapped = map_something(value)
```

---

## 🔄 Question 3: Orchestrator Integration

### The Execution Pattern

Every agent follows this **exact pattern** in `orchestrator.py`:

```python
async def _run_YOUR_agent(self, state: PatientState) -> Dict:
    """Run the YOUR_NAME agent."""
    
    # 1. Start execution tracking
    execution = state.start_agent('your_agent_name')
    
    try:
        # 2. Import your service/agent
        from ..your_module import YourAgent
        
        # 3. Build request from PatientState
        request_data = {
            'field1': state.field1,
            'field2': state.field2
        }
        
        # 4. Call your service
        agent = YourAgent()
        result_obj = await agent.process(request_data)
        
        # 5. Convert to dict (if dataclass)
        if is_dataclass(result_obj):
            result = self._dataclass_to_dict(result_obj)
        else:
            result = result_obj
        
        # 6. Mark execution complete
        execution.complete(result)
        
        # 7. Return result
        return result
        
    except Exception as e:
        # 8. Mark execution failed
        execution.fail(str(e))
        raise
```

### Where to Add Your Agent

**File:** `api/services/orchestrator/orchestrator.py`

**Location:** Search for `async def _run_nutrition_agent` (or your agent name)

**Pattern:**
1. Find the TODO method for your agent (lines ~596-614 for nutrition)
2. Replace the placeholder with your implementation
3. Follow the 8-step pattern above

### Wiring to the Pipeline

Your agent is called in **Phase 2** (parallel execution):

```python
# orchestrator.py:274-298
# Phase 2: Parallel Analysis (Biomarker, Resistance, Nutrition)
tasks = []
task_names = []

if 'biomarker' not in skip_agents:
    tasks.append(self._run_biomarker_agent(state))
    task_names.append('biomarker')

if 'resistance' not in skip_agents:
    tasks.append(self._run_resistance_agent(state))
    task_names.append('resistance')

# Nutrition agent (can run in parallel)
if 'nutrition' not in skip_agents:
    tasks.append(self._run_nutrition_agent(state))  # ← Your method
    task_names.append('nutrition')

# Run in parallel
results = await asyncio.gather(*tasks, return_exceptions=True)

# Store results
for name, result in zip(task_names, results):
    if isinstance(result, Exception):
        logger.error(f"{name} agent failed: {result}")
    else:
        if name == 'biomarker':
            state.biomarker_profile = result
        elif name == 'resistance':
            state.resistance_prediction = result
        elif name == 'nutrition':
            state.nutrition_plan = result  # ← Your result stored here
```

**Key Points:**
- Your agent runs **in parallel** with biomarker and resistance
- Result is stored in `state.your_field` (e.g., `state.nutrition_plan`)
- Exceptions are caught and logged (pipeline continues)

### Updating PatientState

**File:** `api/services/orchestrator/state.py`

**Your field should already exist:**
```python
@dataclass
class PatientState:
    """Central state container for patient data and agent outputs."""
    
    # ... other fields ...
    
    # Agent outputs (Module-specific)
    biomarker_profile: Optional[Dict] = None    # From 02_BIOMARKER
    resistance_prediction: Optional[Dict] = None # From 03_RESISTANCE
    nutrition_plan: Optional[Dict] = None        # From 06_NUTRITION ← YOU
    drug_ranking: Optional[Dict] = None          # From 04_EFFICACY
    trial_matches: Optional[List] = None         # From 05_TRIALS
```

**If your field is missing, add it:**
1. Add field to `PatientState` dataclass
2. Update `state_store.py` to persist it
3. Update `to_dict()` method

---

## 🎨 Care Plan Integration

### The Automatic Integration Pattern

**Question:** Should I modify the care plan agent to include my output?

**Answer:** ✅ **NO** - Care plan automatically consumes your data

### How It Works

**Care Plan Agent** reads from `PatientState` directly:

```python
# orchestrator.py:693-726
async def _run_care_plan_agent(self, state: PatientState) -> Dict:
    """Run the care plan generation agent."""
    return {
        'patient_id': state.patient_id,
        'sections': [
            {
                'title': 'Patient Summary',
                'content': {
                    'mutations': len(state.mutations),
                    'biomarkers': state.biomarker_profile  # ← Reads from state
                }
            },
            {
                'title': 'Resistance Assessment',
                'content': state.resistance_prediction  # ← Reads from state
            },
            {
                'title': 'Drug Recommendations',
                'content': state.drug_ranking  # ← Reads from state
            },
            {
                'title': 'Clinical Trial Options',
                'content': state.trial_matches  # ← Reads from state
            },
            {
                'title': 'Monitoring Plan',
                'content': state.monitoring_config  # ← Reads from state
            }
            # Nutrition section will be added by Senior Agent
        ]
    }
```

### Your Responsibility

1. ✅ **Build your agent** with correct output structure
2. ✅ **Wire to orchestrator** - populate `state.your_field`
3. ✅ **Document output schema** in your MDC file
4. ❌ **Don't modify care plan** - Senior Agent handles aggregation

### Output Structure Requirements

**From MDC Spec** (e.g., `06_NUTRITION_AGENT.mdc`):

```python
# Your agent should return this structure
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
    'foods_to_avoid': [...],
    'drug_food_interactions': [...],
    'timing_rules': {},
    'provenance': {}
}
```

**Key Points:**
- Return a **dict** (not a dataclass) - or use `_dataclass_to_dict`
- Include all fields from MDC schema
- Use `provenance` for tracking metadata

---

## 📚 Implementation Examples

### Example 1: Nutrition Agent (Full Implementation Template)

```python
# api/services/nutrition/nutrition_agent.py

from typing import Dict, List, Optional, Any
from dataclasses import dataclass

class NutritionAgent:
    """Generate toxicity-aware nutrition plans."""
    
    def __init__(self):
        """Initialize with required services."""
        # Import dependencies
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
        """
        Generate comprehensive nutrition plan.
        
        Args:
            patient_id: Patient identifier
            disease: Disease type (e.g., 'ovarian_cancer')
            mutations: List of mutation dicts
            current_medications: List of drug names
            germline_genes: List of germline gene symbols
            options: Additional options
        
        Returns:
            NutritionPlan dict with supplements, foods, interactions, timing
        """
        options = options or {}
        current_medications = current_medications or []
        germline_genes = germline_genes or []
        
        # 1. Get toxicity-mitigating foods (MOAT)
        supplements = []
        if current_medications:
            for drug in current_medications:
                moa = self.get_drug_moa(drug)
                if moa != 'unknown' and germline_genes:
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
        interactions = await self._check_interactions(current_medications)
        
        # 3. Get disease-specific food recommendations
        foods_to_prioritize = self._get_disease_foods(disease)
        
        # 4. Generate timing rules
        timing_rules = self._generate_timing_rules(supplements, current_medications)
        
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
    
    async def _check_interactions(self, medications: List[str]) -> List[Dict]:
        """Check drug-food interactions."""
        # Import interaction checker
        from .drug_food_interactions import check_interactions
        return check_interactions(medications)
    
    def _get_disease_foods(self, disease: str) -> List[Dict]:
        """Get disease-specific food recommendations."""
        # Simplified - could use SPE service for validation
        disease_map = {
            'ovarian_cancer': [
                {'food': 'Cruciferous vegetables', 'reason': 'Glucosinolates support detox', 'frequency': 'daily', 'category': 'prioritize'},
                {'food': 'Omega-3 rich fish', 'reason': 'Anti-inflammatory', 'frequency': '2-3x/week', 'category': 'include'}
            ],
            'breast_cancer': [
                {'food': 'Flaxseed', 'reason': 'Lignans modulate estrogen', 'frequency': 'daily', 'category': 'prioritize'}
            ]
        }
        return disease_map.get(disease.lower(), [])
    
    def _generate_timing_rules(self, supplements: List[Dict], medications: List[str]) -> Dict[str, str]:
        """Generate timing rules for supplements."""
        rules = {}
        for supp in supplements:
            name = supp['name']
            timing = supp.get('timing', 'with meals')
            rules[name] = timing
        return rules
```

### Example 2: Wiring to Orchestrator

```python
# api/services/orchestrator/orchestrator.py

async def _run_nutrition_agent(self, state: PatientState) -> Dict:
    """Run the nutrition planning agent."""
    execution = state.start_agent('nutrition')
    
    try:
        # Import nutrition agent
        from ..nutrition import NutritionAgent
        
        # Extract germline genes from mutations
        germline_genes = []
        if state.mutations:
            germline_genes = [m.get('gene') for m in state.mutations if m.get('gene')]
        
        # Extract current medications from profile
        current_medications = []
        if state.patient_profile:
            if isinstance(state.patient_profile, dict):
                current_medications = state.patient_profile.get('current_medications', [])
            else:
                current_medications = getattr(state.patient_profile, 'current_medications', [])
        
        # Build agent and generate plan
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
        # Return minimal fallback
        return {
            'patient_id': state.patient_id,
            'supplements': [],
            'foods_to_prioritize': [],
            'foods_to_avoid': [],
            'drug_food_interactions': [],
            'timing_rules': {},
            'error': str(e)
        }
```

---

## 🧪 Testing Your Agent

### Unit Tests

**File:** `api/services/nutrition/tests/test_nutrition_agent.py`

```python
import pytest
from api.services.nutrition import NutritionAgent

@pytest.mark.asyncio
async def test_generate_nutrition_plan_basic():
    """Test basic nutrition plan generation."""
    agent = NutritionAgent()
    
    result = await agent.generate_nutrition_plan(
        patient_id='PT-TEST',
        disease='ovarian_cancer',
        mutations=[{'gene': 'BRCA1', 'variant': 'V600E'}],
        current_medications=['carboplatin'],
        germline_genes=['BRCA1']
    )
    
    assert result['patient_id'] == 'PT-TEST'
    assert len(result['supplements']) > 0
    assert 'provenance' in result

@pytest.mark.asyncio
async def test_toxicity_mitigation():
    """Test MOAT integration - toxicity mitigation."""
    agent = NutritionAgent()
    
    result = await agent.generate_nutrition_plan(
        patient_id='PT-MOAT',
        disease='ovarian_cancer',
        mutations=[],
        current_medications=['carboplatin'],
        germline_genes=['BRCA1']  # DNA repair pathway
    )
    
    # Should recommend NAC for platinum + BRCA1
    supplement_names = [s['name'] for s in result['supplements']]
    assert any('NAC' in name for name in supplement_names)
```

### Integration Tests

**Test with orchestrator:**

```python
# tests/test_orchestrator_nutrition.py

import pytest
from api.services.orchestrator import get_orchestrator
from api.services.orchestrator.state import PatientState

@pytest.mark.asyncio
async def test_orchestrator_nutrition_integration():
    """Test nutrition agent integration with orchestrator."""
    orchestrator = get_orchestrator()
    
    state = PatientState(
        patient_id='PT-ORCH',
        disease='ovarian_cancer',
        mutations=[{'gene': 'BRCA1'}],
        patient_profile={'current_medications': ['carboplatin']}
    )
    
    # Run nutrition agent
    result = await orchestrator._run_nutrition_agent(state)
    
    assert result is not None
    assert 'supplements' in result
    assert state.nutrition_plan is not None
```

---

## 📋 Checklist for Your Agent

### Phase 1: Setup
- [ ] Read your MDC spec (`.cursor/MOAT/orchestration/XX_MODULE.mdc`)
- [ ] Identify dependencies (what services/data you need)
- [ ] Create directory structure (`api/services/your_module/`)
- [ ] Create `__init__.py` and main agent file

### Phase 2: Core Implementation
- [ ] Create agent class with main method
- [ ] Import required services (direct imports, not HTTP)
- [ ] Implement core logic
- [ ] Return dict matching MDC schema
- [ ] Add error handling

### Phase 3: Orchestrator Integration
- [ ] Find your `_run_your_agent` method in `orchestrator.py`
- [ ] Import your agent class
- [ ] Build request from `PatientState`
- [ ] Call your agent
- [ ] Convert result to dict if needed
- [ ] Use `execution.complete(result)` pattern
- [ ] Verify `state.your_field` is populated

### Phase 4: Testing
- [ ] Write unit tests for agent class
- [ ] Write integration test with orchestrator
- [ ] Test error handling
- [ ] Validate output schema

### Phase 5: Documentation
- [ ] Update `AGENT_OWNERSHIP_REGISTRY.md` status
- [ ] Document any deviations from MDC spec
- [ ] Add inline code comments
- [ ] Update integration documentation

---

## 🔍 Troubleshooting

### Common Issues

**Issue 1: Import errors**
```
ModuleNotFoundError: No module named 'api.services.your_module'
```
**Solution:** Check `__init__.py` exists and exports your class

**Issue 2: State field not persisting**
```
state.your_field is None after running agent
```
**Solution:** Ensure orchestrator assigns result: `state.your_field = result`

**Issue 3: Dataclass serialization error**
```
Object of type YourClass is not JSON serializable
```
**Solution:** Use `self._dataclass_to_dict(result_obj)` helper

**Issue 4: Agent not running**
```
Your agent method never executes
```
**Solution:** Check it's added to parallel tasks in orchestrator

---

## 📚 Reference Documents

| Document | Purpose |
|----------|---------|
| `.cursor/MOAT/AGENT_OWNERSHIP_REGISTRY.md` | Who owns what |
| `.cursor/MOAT/orchestration/XX_MODULE.mdc` | Your module spec |
| `.cursor/MOAT/orchestration/00_MASTER_INDEX.mdc` | Module directory |
| `api/services/orchestrator/orchestrator.py` | Integration point |
| `api/services/orchestrator/state.py` | State structure |

---

## ✅ Success Criteria

Your agent is complete when:

1. ✅ Returns dict matching MDC output schema
2. ✅ Wired to orchestrator (replaces TODO)
3. ✅ Populates `state.your_field` correctly
4. ✅ Uses direct service imports (no HTTP calls)
5. ✅ Has unit tests with >80% coverage
6. ✅ Has integration test with orchestrator
7. ✅ Documentation updated in registry
8. ✅ Follows established pattern (execution tracking, error handling)

---

**Guide Status:** ✅ ACTIVE  
**Last Updated:** January 28, 2025  
**Maintained By:** Senior Agent + JR Agent G (example implementation)

**Questions?** Reference existing agent implementations:
- `_run_biomarker_agent` (lines 450-483)
- `_run_resistance_agent` (lines 485-565)
- `_run_trial_matching_agent` (lines 625-691)

