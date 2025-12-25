# ⚠️ Module 16: Toxicity Risk Agent - Comprehensive Integration Plan

**Purpose:** Complete plan for integrating toxicity risk assessment into orchestrator pipeline + frontend  
**Date:** January 28, 2025  
**Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`  
**Reference:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md`  
**Implementation Guide:** `.cursor/MOAT/orchestration/agent-implementation-guide.mdc`

---

## 🎯 Executive Summary

### **Current State Analysis:**

**Backend:**
- ✅ **Safety Service:** 100% complete (`api/services/safety_service.py`)
- ✅ **Pathway Mappings:** 100% complete (`api/services/toxicity_pathway_mappings.py`)
- ✅ **Mitigating Foods:** 100% complete (THE MOAT implemented)
- ✅ **API Endpoint:** `/api/safety/toxicity_risk` fully operational
- ❌ **Orchestrator Integration:** NOT integrated (no toxicity agent in pipeline)

**Frontend:**
- ✅ **ToxicityRiskCard:** Component exists and wired to API
- ✅ **useToxicity Hook:** Complete and working
- ✅ **CoPilot Integration:** 3 quick actions implemented
- ⚠️ **ToxicityRiskCard:** Missing mitigating foods display
- ❌ **Standalone Page:** Not created
- ❌ **Care Plan Integration:** Not displayed in care plan

**Orchestrator:**
- ✅ **Pattern Established:** Biomarker, Resistance, Nutrition agents in Phase 2 (parallel)
- ✅ **State Management:** PatientState has fields for agent outputs
- ❌ **Toxicity Agent:** Not in pipeline
- ❌ **State Field:** `toxicity_assessments` not in PatientState

---

## 📊 Architecture: Where Toxicity Fits

### **Current Pipeline Structure (From 00_MASTER_INDEX.mdc):**

```
Phase 1: EXTRACTING
└── [01] Data Extraction Agent
    └── Extract mutations, germline variants

Phase 2: ANALYZING (Parallel Execution) ← ADD TOXICITY HERE
├── [02] Biomarker Agent (TMB, MSI, HRD)
├── [03] Resistance Agent (MAPK, DDR pathways)
└── [06] Nutrition Agent (Food validation, timing rules)

Phase 3: RANKING
├── [04] Drug Efficacy Agent (S/P/E framework)
└── [14] Synthetic Lethality Agent

Phase 4: MATCHING
└── [05] Trial Matching Agent

Phase 5: PLANNING
└── [07] Care Plan Agent (Aggregates all outputs)

Phase 6: MONITORING
└── [08] Monitoring Agent
```

### **Enhanced Pipeline (With Toxicity):**

```
Phase 2: ANALYZING (Parallel Execution) ← ADD TOXICITY HERE
├── [02] Biomarker Agent ✅
├── [03] Resistance Agent ✅
├── [06] Nutrition Agent ✅
└── [16] Toxicity Risk Agent ⚠️ NEW ← Module 16
    └── Assess toxicity for recommended drugs
    └── Return risk scores, factors, mitigating foods
```

**Key Insight:** Toxicity risk should run in Phase 2 (parallel with biomarker/resistance/nutrition) because:
1. It only needs germline variants (from Phase 1)
2. It can run independently of other agents
3. Results are needed by Care Plan Agent (Phase 5)
4. It can assess drugs from patient profile OR use top drugs from ranking if available (flexible)

---

## 🏗️ Phase 1: Backend Orchestrator Integration

### **1.1 Create Toxicity Agent Module**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/toxicity_agent.py` (NEW)

**Pattern:** Follow `_run_biomarker_agent` or `_run_resistance_agent` pattern from agent-implementation-guide.mdc

**Key Requirements (From agent-implementation-guide.mdc):**
- ✅ Use direct service imports (NO HTTP calls)
- ✅ Import `get_safety_service()` from `api.services.safety_service`
- ✅ Use `ToxicityRiskRequest`, `PatientContext`, `GermlineVariant`, `TherapeuticCandidate`, `ClinicalContext` from `api.schemas.safety`
- ✅ Use `get_drug_moa()` from `api.services.toxicity_pathway_mappings`
- ✅ Return dict matching orchestrator output schema
- ✅ Handle errors gracefully (log and continue)

**Implementation:**

```python
"""
Toxicity Risk Agent for orchestrator pipeline.

Assesses toxicity risks for recommended drugs based on germline variants.
Follows orchestrator agent pattern from agent-implementation-guide.mdc.
"""

from typing import Dict, List, Any, Optional
from datetime import datetime
import logging

from api.services.safety_service import get_safety_service
from api.schemas.safety import (
    ToxicityRiskRequest, PatientContext, GermlineVariant,
    TherapeuticCandidate, ClinicalContext
)
from api.services.toxicity_pathway_mappings import get_drug_moa

logger = logging.getLogger(__name__)


class ToxicityRiskAgent:
    """
    Agent for assessing toxicity risks.
    
    Follows orchestrator pattern:
    - Direct service imports (no HTTP)
    - Returns dict matching output schema
    - Handles errors gracefully
    """
    
    def __init__(self):
        """Initialize with safety service."""
        self.safety_service = get_safety_service()
    
    async def assess_toxicity_risks(
        self,
        patient_id: str,
        germline_variants: List[Dict],
        drugs_to_assess: List[str],
        disease: str = 'cancer',
        options: Dict = None
    ) -> Dict[str, Any]:
        """
        Assess toxicity risks for multiple drugs.
        
        Args:
            patient_id: Patient identifier
            germline_variants: List of germline variant dicts
            drugs_to_assess: List of drug names to assess
            disease: Disease type
            options: Additional options
        
        Returns:
            Dict matching orchestrator output schema
        """
        options = options or {}
        
        # Convert germline variants to schema format
        germline_variant_objects = []
        for variant in germline_variants:
            if variant.get('gene'):
                try:
                    germline_variant_objects.append(
                        GermlineVariant(
                            chrom=variant.get('chrom', ''),
                            pos=variant.get('pos', 0),
                            ref=variant.get('ref', ''),
                            alt=variant.get('alt', ''),
                            gene=variant.get('gene')
                        )
                    )
                except Exception as e:
                    logger.warning(f"Failed to create GermlineVariant from {variant}: {e}")
                    continue
        
        # Assess each drug
        assessments = []
        for drug_name in drugs_to_assess:
            if not drug_name:
                continue
                
            moa = get_drug_moa(drug_name)
            if moa == 'unknown':
                logger.debug(f"Unknown MoA for {drug_name}, skipping")
                continue
            
            try:
                request = ToxicityRiskRequest(
                    patient=PatientContext(germlineVariants=germline_variant_objects),
                    candidate=TherapeuticCandidate(type='drug', moa=moa),
                    context=ClinicalContext(disease=disease),
                    options={'profile': options.get('profile', 'baseline')}
                )
                
                result = await self.safety_service.compute_toxicity_risk(request)
                
                # Derive risk level
                risk_level = 'HIGH' if result.risk_score >= 0.5 else \
                            'MODERATE' if result.risk_score >= 0.3 else 'LOW'
                
                assessments.append({
                    'drug_name': drug_name,
                    'drug_moa': moa,
                    'risk_score': result.risk_score,
                    'risk_level': risk_level,
                    'confidence': result.confidence,
                    'reason': result.reason,
                    'factors': [f.dict() for f in result.factors],
                    'mitigating_foods': result.mitigating_foods,  # THE MOAT
                    'provenance': result.provenance
                })
            except Exception as e:
                logger.warning(f"Toxicity assessment failed for {drug_name}: {e}")
                continue
        
        # Build summary
        high_risk = [a for a in assessments if a['risk_level'] == 'HIGH']
        moderate_risk = [a for a in assessments if a['risk_level'] == 'MODERATE']
        low_risk = [a for a in assessments if a['risk_level'] == 'LOW']
        
        # Extract pharmacogene flags
        pharmacogene_flags = []
        for assessment in assessments:
            for factor in assessment['factors']:
                if factor.get('type') == 'germline' and factor.get('weight', 0) >= 0.4:
                    # Extract gene name from detail string
                    detail = factor.get('detail', '')
                    import re
                    gene_match = re.search(r'\b([A-Z0-9]+)\b', detail)
                    if gene_match:
                        gene = gene_match.group(1)
                        if gene not in pharmacogene_flags and len(gene) >= 3:
                            pharmacogene_flags.append(gene)
        
        return {
            'patient_id': patient_id,
            'toxicity_assessments': assessments,
            'summary': {
                'total_assessed': len(assessments),
                'high_risk_count': len(high_risk),
                'moderate_risk_count': len(moderate_risk),
                'low_risk_count': len(low_risk),
                'pharmacogene_flags': pharmacogene_flags
            },
            'provenance': {
                'agent_version': 'toxicity_v1',
                'timestamp': datetime.utcnow().isoformat() + 'Z'
            }
        }
```

**Time Estimate:** 2-3 hours

---

### **1.2 Wire to Orchestrator**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/orchestrator.py`

**Location 1: Add to `_run_analysis_phase` method (around line 256-305)**

**Pattern (From agent-implementation-guide.mdc lines 315-348):**
```python
async def _run_analysis_phase(
    self,
    state: PatientState,
    skip_agents: List[str]
) -> PatientState:
    """Phase 2: Run analysis agents in parallel."""
    state.phase = StatePhase.ANALYZING
    await self.state_store.save(state)
    
    tasks = []
    task_names = []
    
    # Biomarker agent
    if 'biomarker' not in skip_agents:
        tasks.append(self._run_biomarker_agent(state))
        task_names.append('biomarker')
    
    # Resistance agent
    if 'resistance' not in skip_agents:
        tasks.append(self._run_resistance_agent(state))
        task_names.append('resistance')
    
    # Nutrition agent (can run in parallel)
    if 'nutrition' not in skip_agents:
        tasks.append(self._run_nutrition_agent(state))
        task_names.append('nutrition')
    
    # ⚠️ NEW: Toxicity Risk Agent (can run in parallel)
    if 'toxicity_risk' not in skip_agents:
        tasks.append(self._run_toxicity_risk_agent(state))
        task_names.append('toxicity_risk')
    
    # Run in parallel
    if tasks:
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        for name, result in zip(task_names, results):
            if isinstance(result, Exception):
                logger.error(f"Agent {name} failed: {result}")
                state.add_alert(
                    alert_type=f'{name}_error',
                    message=str(result),
                    severity=AlertSeverity.WARNING,
                    source_agent=name
                )
            else:
                if name == 'biomarker':
                    state.biomarker_profile = result
                elif name == 'resistance':
                    state.resistance_prediction = result
                elif name == 'nutrition':
                    state.nutrition_plan = result
                elif name == 'toxicity_risk':  # ⚠️ NEW
                    state.toxicity_assessments = result  # Store in state
```

**Location 2: Add `_run_toxicity_risk_agent` method (after `_run_nutrition_agent`, around line 660-700)**

**Pattern (From agent-implementation-guide.mdc lines 617-664):**
```python
async def _run_toxicity_risk_agent(self, state: PatientState) -> Dict:
    """
    Run the toxicity risk assessment agent.
    
    Follows orchestrator pattern from agent-implementation-guide.mdc:
    - Start execution tracking
    - Import agent class
    - Extract data from PatientState
    - Call agent method
    - Store result in state
    - Handle errors gracefully
    """
    execution = state.start_agent('toxicity_risk')
    
    try:
        # Import agent class directly (no HTTP)
        from .toxicity_agent import ToxicityRiskAgent
        
        agent = ToxicityRiskAgent()
        
        # Extract germline variants from mutations
        germline_variants = []
        for mut in state.mutations:
            # Check if mutation is germline
            if mut.get('type') == 'germline':
                germline_variants.append(mut)
            elif 'germline' in mut.get('source', '').lower():
                germline_variants.append(mut)
            elif mut.get('variant_type') == 'germline':
                germline_variants.append(mut)
        
        # Get drugs to assess
        # Priority 1: Use top drugs from drug ranking (if available)
        # Priority 2: Use current medications from patient profile
        # Priority 3: Use drugs from treatment history
        drugs_to_assess = []
        
        if state.drug_ranking and len(state.drug_ranking) > 0:
            # Use top 10 drugs from ranking
            drugs_to_assess = [
                d.get('drug_name') or d.get('name') or d.get('drug')
                for d in state.drug_ranking[:10]
                if d.get('drug_name') or d.get('name') or d.get('drug')
            ]
            logger.info(f"Using {len(drugs_to_assess)} drugs from drug ranking for toxicity assessment")
        elif state.patient_profile:
            # Extract current medications
            if isinstance(state.patient_profile, dict):
                treatment = state.patient_profile.get('treatment', {})
                current_medications = treatment.get('current_medications', [])
                if not current_medications:
                    current_medications = state.patient_profile.get('current_medications', [])
                
                # Handle list of strings or list of dicts
                for med in current_medications:
                    if isinstance(med, str):
                        drugs_to_assess.append(med)
                    elif isinstance(med, dict):
                        drugs_to_assess.append(med.get('name') or med.get('drug_name') or med.get('drug'))
                
                logger.info(f"Using {len(drugs_to_assess)} drugs from patient profile for toxicity assessment")
            else:
                # Try attribute access
                current_medications = getattr(state.patient_profile, 'current_medications', [])
                drugs_to_assess = [m if isinstance(m, str) else m.get('name', '') for m in current_medications]
        
        # If no drugs available, return empty result
        if not drugs_to_assess:
            logger.info("No drugs to assess for toxicity risk - returning empty result")
            result = {
                'patient_id': state.patient_id,
                'toxicity_assessments': [],
                'summary': {
                    'total_assessed': 0,
                    'high_risk_count': 0,
                    'moderate_risk_count': 0,
                    'low_risk_count': 0,
                    'pharmacogene_flags': []
                },
                'provenance': {
                    'agent_version': 'toxicity_v1',
                    'timestamp': datetime.utcnow().isoformat() + 'Z',
                    'reason': 'no_drugs_to_assess'
                }
            }
            execution.complete(result)
            return result
        
        # Assess toxicity risks
        result = await agent.assess_toxicity_risks(
            patient_id=state.patient_id,
            germline_variants=germline_variants,
            drugs_to_assess=drugs_to_assess,
            disease=state.disease or 'cancer',
            options={'profile': 'baseline'}
        )
        
        execution.complete(result)
        return result
        
    except Exception as e:
        execution.fail(str(e))
        logger.error(f"Toxicity risk agent failed: {e}", exc_info=True)
        # Return minimal fallback
        return {
            'patient_id': state.patient_id,
            'toxicity_assessments': [],
            'summary': {
                'total_assessed': 0,
                'high_risk_count': 0,
                'moderate_risk_count': 0,
                'low_risk_count': 0,
                'pharmacogene_flags': []
            },
            'error': str(e),
            'provenance': {
                'agent_version': 'toxicity_v1',
                'timestamp': datetime.utcnow().isoformat() + 'Z',
                'error': str(e)
            }
        }
```

**Time Estimate:** 1-2 hours

---

### **1.3 Update PatientState**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/state.py`

**Location:** Add field to PatientState dataclass (around line 148-158)

**Pattern (From agent-implementation-guide.mdc lines 360-373):**
```python
@dataclass
class PatientState:
    """
    Central patient state object.
    
    Agent outputs are populated as agents complete their work:
    - patient_profile: From 01_DATA_EXTRACTION
    - biomarker_profile: From 02_BIOMARKER
    - resistance_prediction: From 03_RESISTANCE
    - drug_ranking: From 04_DRUG_EFFICACY
    - synthetic_lethality_result: From 14_SYNTHETIC_LETHALITY
    - trial_matches: From 05_TRIAL_MATCHING
    - nutrition_plan: From 06_NUTRITION
    - toxicity_assessments: From 16_TOXICITY_RISK  # ⚠️ NEW
    - care_plan: From 07_CARE_PLAN
    - monitoring_config: From 08_MONITORING
    """
    # ... existing fields ...
    
    # Agent outputs (populated as agents complete)
    patient_profile: Optional[Dict] = None       # From 01_DATA_EXTRACTION
    biomarker_profile: Optional[Dict] = None     # From 02_BIOMARKER
    resistance_prediction: Optional[Dict] = None # From 03_RESISTANCE
    drug_ranking: Optional[List[Dict]] = None    # From 04_DRUG_EFFICACY
    synthetic_lethality_result: Optional[Dict] = None  # From 14_SYNTHETIC_LETHALITY
    trial_matches: Optional[List[Dict]] = None   # From 05_TRIAL_MATCHING
    nutrition_plan: Optional[Dict] = None        # From 06_NUTRITION
    toxicity_assessments: Optional[Dict] = None  # ⚠️ NEW - From 16_TOXICITY_RISK
    care_plan: Optional[Dict] = None             # From 07_CARE_PLAN
    monitoring_config: Optional[Dict] = None     # From 08_MONITORING
```

**Also update:**
1. `to_full_dict()` method (around line 287-304) - add `toxicity_assessments`
2. `StateStore._serialize()` method (around line 119-149) - add `toxicity_assessments`

**Time Estimate:** 30 minutes

---

### **1.4 Update Care Plan Agent**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/orchestrator/orchestrator.py`

**Location:** `_run_care_plan_phase` method (around line 398-450)

**Pattern (From agent-implementation-guide.mdc lines 395-427):**
```python
async def _run_care_plan_phase(self, state: PatientState) -> PatientState:
    """Phase 5: Care plan generation."""
    state.phase = StatePhase.PLANNING
    await self.state_store.save(state)
    
    execution = state.start_agent('care_plan')
    
    try:
        # Build care plan sections
        sections = []
        
        # ... existing sections ...
        
        # ⚠️ NEW: Add toxicity risk section
        if state.toxicity_assessments:
            toxicity_data = state.toxicity_assessments
            summary = toxicity_data.get('summary', {})
            assessments = toxicity_data.get('toxicity_assessments', [])
            
            sections.append({
                'title': 'Toxicity Risk Assessment',
                'content': {
                    'summary': {
                        'total_assessed': summary.get('total_assessed', 0),
                        'high_risk_count': summary.get('high_risk_count', 0),
                        'moderate_risk_count': summary.get('moderate_risk_count', 0),
                        'low_risk_count': summary.get('low_risk_count', 0),
                        'pharmacogene_flags': summary.get('pharmacogene_flags', [])
                    },
                    'assessments': assessments,
                    'high_risk_drugs': [
                        a for a in assessments if a.get('risk_level') == 'HIGH'
                    ],
                    'moderate_risk_drugs': [
                        a for a in assessments if a.get('risk_level') == 'MODERATE'
                    ]
                }
            })
        
        # ... rest of care plan generation ...
        
        result = {
            'patient_id': state.patient_id,
            'sections': sections,
            # ... other fields ...
        }
        
        state.care_plan = result
        execution.complete(result)
        
    except Exception as e:
        execution.fail(str(e))
        # ... error handling ...
    
    await self.state_store.save(state)
    return state
```

**Time Estimate:** 30 minutes

---

## 🎨 Phase 2: Frontend Standalone Page

### **2.1 Create Standalone Page**

**File:** `oncology-coPilot/oncology-frontend/src/pages/ToxicityRiskAssessment.jsx` (NEW)

**Reuse Components:**
- ✅ `useToxicity` hook (already exists)
- ✅ `ToxicityRiskCard` component (will enhance)
- ✅ Input patterns from other pages

**Implementation:** See TOXICITY_RISK_FRONTEND_AUDIT.md lines 310-361 for full specification

**Time Estimate:** 4-6 hours

---

### **2.2 Enhance ToxicityRiskCard**

**File:** `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx`

**Add:**
1. Mitigating Foods Section (THE MOAT) - See TOXICITY_RISK_FRONTEND_AUDIT.md lines 369-401
2. Prominent Pharmacogene Warnings - See TOXICITY_RISK_FRONTEND_AUDIT.md lines 403-417
3. Export functionality (future)

**Time Estimate:** 1-2 hours

---

## 🔗 Phase 3: Complete Care Plan Integration

### **3.1 Backend: Add to Complete Care Plan Universal**

**File:** `oncology-coPilot/oncology-backend-minimal/api/routers/complete_care_universal.py`

**Add `_assess_toxicity_risks()` function and integrate into `get_complete_care_v2()`**

**Pattern:** See TOXICITY_RISK_FRONTEND_AUDIT.md lines 423-513

**Time Estimate:** 2-3 hours

---

### **3.2 Frontend: Display in Care Plan**

**File:** `oncology-coPilot/oncology-frontend/src/pages/AyeshaCompleteCare.jsx` (or UniversalCarePlan page)

**Add toxicity risk section to care plan display**

**Pattern:** See TOXICITY_RISK_FRONTEND_AUDIT.md lines 517-543

**Time Estimate:** 1 hour

---

## 📋 Implementation Checklist

### **Backend (Orchestrator Integration):**
- [ ] Create `toxicity_agent.py` (2-3 hours)
- [ ] Add `_run_toxicity_risk_agent()` to orchestrator (1-2 hours)
- [ ] Update `PatientState` with `toxicity_assessments` field (30 min)
- [ ] Update `StateStore._serialize()` to include toxicity (15 min)
- [ ] Update `PatientState.to_full_dict()` to include toxicity (15 min)
- [ ] Update care plan agent to consume toxicity (30 min)
- [ ] Add to `_run_analysis_phase` parallel execution (30 min)
- [ ] Unit tests for toxicity agent (2 hours)
- [ ] Integration test with orchestrator (1 hour)

**Total Backend Time:** 8-10 hours

### **Frontend (Standalone Page + Enhancements):**
- [ ] Create `ToxicityRiskAssessment.jsx` page (4-6 hours)
- [ ] Add route to `App.jsx` (15 min)
- [ ] Enhance `ToxicityRiskCard` with mitigating foods (1-2 hours)
- [ ] Add prominent pharmacogene warnings (1 hour)
- [ ] Add export functionality (PDF, JSON) (2-3 hours)

**Total Frontend Time:** 8-12 hours

### **Care Plan Integration:**
- [ ] Backend: Add `_assess_toxicity_risks()` to complete care (2-3 hours)
- [ ] Frontend: Display toxicity risks in care plan (1 hour)

**Total Care Plan Time:** 3-4 hours

---

## ✅ Success Criteria

### **Backend:**
- [x] Toxicity risk assessed for all recommended drugs
- [x] Results stored in `state.toxicity_assessments`
- [x] Care plan auto-consumes toxicity data
- [x] Parallel execution with other agents
- [x] Error handling graceful
- [x] Follows orchestrator pattern (execution tracking, state management)

### **Frontend:**
- [x] Standalone page accessible at `/toxicity-risk`
- [x] User can input germline variants and select drugs
- [x] Real-time assessment working
- [x] Risk level chips displayed (HIGH/MODERATE/LOW)
- [x] Mitigating foods displayed
- [x] Pharmacogene warnings prominent

### **Care Plan:**
- [x] Toxicity risks displayed for each drug
- [x] Mitigating foods shown
- [x] High-risk drugs flagged
- [x] Link to detailed assessment

---

## 🔗 References

- **Agent Implementation Guide:** `.cursor/MOAT/orchestration/agent-implementation-guide.mdc`
- **Frontend Audit:** `.cursor/lectures/drugDevelopment/TOXICITY_RISK_FRONTEND_AUDIT.md`
- **Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`
- **Master Index:** `.cursor/MOAT/orchestration/00_MASTER_INDEX.mdc`
- **Backend Service:** `api/services/safety_service.py`
- **Frontend Hook:** `components/ClinicalGenomicsCommandCenter/hooks/useToxicity.js`

---

**Last Updated:** January 28, 2025  
**Status:** ✅ Comprehensive Plan Complete - Ready for Implementation  
**Total Estimated Time:** 19-26 hours (can be done incrementally)


