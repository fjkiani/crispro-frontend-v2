# 🔬 Comprehensive MOAT Capabilities Analysis

**Date:** January 2025  
**Purpose:** Complete understanding of all MOAT capabilities, data flows, and frontend/backend contracts before implementation  
**Status:** ✅ **ANALYSIS COMPLETE** - Updated with actual implementation status  
**Last Updated:** January 2025 (aligned with actual codebase)

---

## 🚨 CRITICAL FINDING: Documentation vs Reality

**The master index (`01_CURRENT_STATE.md`) is OUTDATED.** The actual implementation is **significantly more complete** than documented.

### **Actual Agent Status:**

| Agent | Documentation Says | **ACTUAL STATUS** | Location |
|-------|-------------------|-------------------|----------|
| **01 Data Extraction** | ⏳ SKELETON | ✅ **INTEGRATED** | `extraction_agent.py` + `_run_extraction_phase()` |
| **02 Biomarker** | ✅ INTEGRATED | ✅ **INTEGRATED** | `_run_biomarker_agent()` line 453 |
| **03 Resistance** | ✅ VALIDATED | ✅ **INTEGRATED** | `_run_resistance_agent()` line 549 |
| **04 Drug Efficacy** | ⏳ SKELETON | ✅ **INTEGRATED** | `_run_drug_efficacy_agent()` line 743 |
| **05 Trial Matching** | ⏳ SKELETON | ✅ **INTEGRATED** | `_run_trial_matching_agent()` line 923 |
| **06 Nutrition** | ⏳ SKELETON | ✅ **INTEGRATED** | `_run_nutrition_agent()` line 660 |
| **07 Care Plan** | ✅ INTEGRATED | ✅ **INTEGRATED** | `_run_care_plan_agent()` line 991 |
| **08 Monitoring** | ✅ INTEGRATED | ✅ **INTEGRATED** | `_run_monitoring_agent()` line 1137 |
| **09 Trigger System** | ⬜ TODO | ⚠️ **EXISTS BUT NOT IN PIPELINE** | `trigger_engine.py` exists, not called |
| **14 Synthetic Lethality** | Not in index | ✅ **INTEGRATED** | `_run_synthetic_lethality_agent()` line 821 |
| **16 Toxicity Risk** | ⬜ TODO | ⚠️ **SERVICE EXISTS, NOT INTEGRATED** | `safety_service.py` exists, no agent |

### **Key Insights:**

1. **9 out of 10 core agents are FULLY INTEGRATED** (not skeletons)
2. **Data Extraction, Drug Efficacy, Nutrition are COMPLETE** (contrary to documentation)
3. **Synthetic Lethality is integrated** (not mentioned in master index)
4. **Trigger System exists but not wired to pipeline**
5. **Toxicity Risk service exists but needs orchestrator integration**
6. **SAE Features can be extracted but not stored in PatientState**

---

## 📋 Table of Contents

1. [Backend Orchestrator Capabilities](#backend-orchestrator-capabilities)
2. [Frontend Component Dependencies](#frontend-component-dependencies)
3. [API Contracts & Data Formats](#api-contracts--data-formats)
4. [Data Flow Mapping](#data-flow-mapping)
5. [Test Requirements](#test-requirements)
6. [Gap Analysis](#gap-analysis)

---

## 🔧 Backend Orchestrator Capabilities

### **Pipeline Phases & Agents**

#### **Phase 1: Data Extraction** (`EXTRACTING`)
- **Agent:** `DataExtractionAgent`
- **Method:** `_run_extraction_phase()` (line 201-254)
- **Service:** `api/services/extraction/extraction_agent.py` ✅ **INTEGRATED**
- **Input:** File (VCF, PDF, MAF, JSON, TXT) or `patient_profile` dict
- **Output:** 
  - `patient_profile`: Dict with patient data
  - `mutations`: List[Dict] of mutations
  - `disease`: str
- **State Updates:**
  - `state.patient_profile`
  - `state.mutations`
  - `state.disease`
- **Quality Validation:** ✅ Includes `_validate_mutation_quality()` with coverage/VAF thresholds

#### **Phase 2: Parallel Analysis** (`ANALYZING`)
Runs three agents in parallel:

##### **2.1 Biomarker Agent** (`biomarker`)
- **Method:** `_run_biomarker_agent()`
- **Input:** `state.mutations`
- **Output:** `biomarker_profile` Dict:
  ```python
  {
    'tmb': {
      'value': float,  # mutations per Mb
      'classification': 'TMB-H' | 'TMB-M' | 'TMB-L',
      'threshold_used': 10.0,
      'mutation_count': int,
      'total_mutations': int,
      'exome_size_mb': 38.0,
      'method': 'nonsilent_variants_per_mb'
    },
    'msi': {
      'status': 'MSI-H' | 'MSS',
      'dmmr_genes_mutated': List[str],
      'method': 'gene_panel_with_impact',
      'confidence': 'high' | 'moderate'
    },
    'hrd': {
      'status': 'HRD+' | 'HRD-' | 'HRD-uncertain',
      'genes_mutated': List[str],
      'method': 'comprehensive_ddr_panel',
      'high_impact_detected': bool,
      'confidence': 'high' | 'moderate'
    },
    'io_eligible': bool,
    'parp_eligible': bool,
    'provenance': {
      'calculation_method': 'enhanced_gene_panel',
      'timestamp': str,
      'disease': str
    }
  }
  ```
- **State Update:** `state.biomarker_profile`

##### **2.2 Resistance Agent** (`resistance`)
- **Method:** `_run_resistance_agent()`
- **Services Used:**
  - `ResistanceProphetService` (prediction)
  - `ResistancePlaybookService` (next-line options)
- **Input:** 
  - `state.mutations`
  - `state.disease`
  - `state.patient_profile` (for treatment context)
- **Output:** `resistance_prediction` Dict:
  ```python
  {
    'risk_level': 'HIGH' | 'MEDIUM' | 'LOW',
    'probability': float,  # 0.0-1.0
    'confidence': float,  # 0.0-1.0
    'detected_genes': List[Dict],  # Genes with resistance signals
    'next_line_options': {
      'alternatives': List[Dict],
      'regimen_changes': List[Dict],
      'monitoring_changes': Dict
    }
  }
  ```
- **State Update:** `state.resistance_prediction`

##### **2.3 Nutrition Agent** (`nutrition`)
- **Method:** `_run_nutrition_agent()` (line 660-741)
- **Service:** `NutritionAgent` (direct import from `api.services.nutrition`)
- **Input:**
  - `state.mutations`
  - `germline_genes` (extracted from patient_profile or mutations)
  - `current_drugs` (from patient_profile or drug_ranking)
  - `state.disease`
  - `treatment_line`
- **Output:** `nutrition_plan` Dict:
  ```python
  {
    'patient_id': str,
    'treatment': str,
    'supplements': List[Dict],
    'foods_to_prioritize': List[str],
    'foods_to_avoid': List[str],
    'drug_food_interactions': List[Dict],
    'timing_rules': Dict,
    'provenance': Dict
  }
  ```
- **State Update:** `state.nutrition_plan`
- **Error Handling:** ✅ Graceful fallback with placeholder on error

#### **Phase 3: Drug Efficacy Ranking** (`RANKING`)
- **Agent:** `drug_efficacy`
- **Method:** `_run_drug_efficacy_agent()` (line 743-819)
- **Service:** `EfficacyOrchestrator` (direct import, no HTTP)
- **Input:**
  - `state.mutations`
  - `state.disease`
  - Model: `evo2_7b`
  - Options: `{'adaptive': True, 'ensemble': True}`
- **Output:** Dict with:
  ```python
  {
    'ranked_drugs': [
      {
        'drug_name': str,
        'moa': str,
        'efficacy_score': float,
        'confidence': float,
        'evidence_tier': str,
        'badges': List[str],
        'rationale': List[str]
      },
      ...
    ],
    'mechanism_vector': List[float],  # 7D vector [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    'evidence_tier': str,
    'run_signature': str,
    'provenance': Dict
  }
  ```
- **State Updates:**
  - `state.drug_ranking`
  - `state.mechanism_vector`
- **SAE Features:** Can be extracted if `include_sae_features: true` in options (line 357-387 in `efficacy_orchestrator.py`), but NOT stored in `PatientState` currently

#### **Phase 3.5: Synthetic Lethality** (`RANKING`)
- **Agent:** `synthetic_lethality`
- **Method:** `_run_synthetic_lethality_agent()`
- **Service:** `SyntheticLethalityAgent`
- **Input:**
  - `state.mutations`
  - `state.disease`
- **Output:** `synthetic_lethality_result` Dict:
  ```python
  {
    'patient_id': str,
    'disease': str,
    'synthetic_lethality_detected': bool,
    'double_hit_description': str,
    'essentiality_scores': [
      {
        'gene': str,
        'essentiality_score': float,
        'essentiality_level': str,
        'sequence_disruption': str,
        'pathway_impact': str,
        'functional_consequence': str,
        'flags': List[str],
        'confidence': float
      },
      ...
    ],
    'broken_pathways': List[Dict],
    'essential_pathways': List[Dict],
    'recommended_drugs': List[Dict],
    'suggested_therapy': str,
    'explanation': Dict,
    'calculation_time_ms': float,
    'evo2_used': bool,
    'provenance': Dict
  }
  ```
- **State Update:** `state.synthetic_lethality_result`

#### **Phase 4: Trial Matching** (`MATCHING`)
- **Agent:** `trial_matching`
- **Method:** `_run_trial_matching_agent()`
- **Service:** `TrialMatchingAgent`
- **Input:**
  - `patient_profile` Dict
  - `biomarker_profile` Dict
  - `mechanism_vector` List[float]
- **Output:** Dict with:
  ```python
  {
    'matches': [
      {
        'nct_id': str,
        'title': str,
        'brief_summary': str,
        'phase': str,
        'status': str,
        'mechanism_fit_score': float,  # 0.0-1.0
        'eligibility_score': float,     # 0.0-1.0
        'combined_score': float,        # 0.7×eligibility + 0.3×mechanism_fit
        'trial_moa': List[float],       # 7D mechanism vector
        'mechanism_alignment': Dict,    # Per-pathway alignment
        'why_matched': str,
        'url': str,
        'locations': List[str],
        'contact': Dict
      },
      ...
    ],
    'queries_used': List[str],
    'trials_found': int,
    'trials_ranked': int,
    'top_match': Dict,
    'search_time_ms': float,
    'provenance': Dict
  }
  ```
- **State Update:** `state.trial_matches`

#### **Phase 5: Care Plan Generation** (`PLANNING`)
- **Agent:** `care_plan`
- **Method:** `_run_care_plan_agent()`
- **Input:** All previous agent outputs
- **Output:** `care_plan` Dict:
  ```python
  {
    'patient_id': str,
    'disease': str,
    'generated_at': str,  # ISO timestamp
    'executive_summary': {
      'patient_id': str,
      'disease': str,
      'mutation_count': int,
      'top_mutations': List[Dict],
      'biomarker_summary': Dict,
      'top_drug': Dict,
      'trial_count': int,
      'risk_level': str
    },
    'sections': [
      {
        'title': str,
        'order': int,
        'content': Any
      },
      ...
    ],
    'alerts': List[Dict],
    'provenance': {
      'pipeline_version': str,
      'generated_by': str,
      'agents_executed': List[str],
      'execution_time_ms': float,
      'data_quality_flags': List[str]
    }
  }
  ```
- **State Update:** `state.care_plan`

#### **Phase 6: Monitoring Setup** (`MONITORING`)
- **Agent:** `monitoring`
- **Method:** `_run_monitoring_agent()`
- **Input:**
  - `state.resistance_prediction`
  - `state.biomarker_profile`
  - `state.disease`
- **Output:** `monitoring_config` Dict:
  ```python
  {
    'frequency': 'weekly' | 'biweekly' | 'monthly',
    'biomarkers': List[str],  # e.g., ['CA-125', 'ctDNA', 'LDH']
    'primary_biomarker': str,
    'imaging': str,  # e.g., 'CT every 2 months'
    'alerts_enabled': bool,
    'escalation_thresholds': {
      'biomarker_rise_percent': int,
      'new_lesion': bool,
      'symptom_worsening': bool,
      'performance_status_drop': bool
    },
    'ctdna_monitoring': {
      'enabled': bool,
      'frequency': str,
      'target_mutations': List[str]
    },
    'risk_level': str,
    'provenance': Dict
  }
  ```
- **State Update:** `state.monitoring_config`

---

## 🎨 Frontend Component Dependencies

### **UniversalCompleteCare.jsx** - Main Page

#### **Current Data Source:** `/api/complete_care/v2`
- **Request Format:**
  ```javascript
  {
    patient_profile: {
      patient_id: string,
      disease: { type: string, stage: string },
      demographics: { age: number, sex: string, ... },
      tumor_context: {
        somatic_mutations: Array<{gene, hgvs_p, consequence, ...}>,
        hrd_score: number
      },
      ...
    },
    include_trials: boolean,
    include_soc: boolean,
    include_biomarker: boolean,
    include_wiwfm: boolean,
    include_resistance: boolean,
    max_trials: number
  }
  ```

#### **Response Format (complete_care/v2):**
```javascript
{
  // Core components
  biomarker_intelligence: {...},      // Maps to BiomarkerCard
  resistance_prediction: {...},       // Maps to ResistanceCard
  wiwfm: {                            // Maps to DrugRankingCard
    drugs: Array<{...}>,
    evidence_tier: string,
    run_signature: string
  },
  trials: {                           // Maps to TrialMatchesCard
    trials: Array<{...}>
  },
  
  // Additional components
  soc_recommendation: {...},          // Maps to SOCRecommendationCard
  next_test_recommender: {...},       // Maps to NextTestCard
  hint_tiles: {...},                  // Maps to HintTilesPanel
  mechanism_map: {                    // Maps to MechanismChips
    vector: Array<number>  // 7D
  },
  resistance_playbook: {...},         // Maps to ResistancePlaybook
  sae_features: {...},                // Maps to AyeshaSAEFeaturesCard
  resistance_alert: {...},            // Maps to ResistanceAlertBanner
  toxicity_assessments: {...},        // Maps to ToxicityRiskCard
  
  // Metadata
  provenance: {
    orchestrator: string,
    run_id: string,
    generated_at: string
  },
  summary: {
    components_included: Array<string>,
    ngs_status: string,
    confidence_level: string
  }
}
```

#### **Component Mapping:**

| Component | Data Source | Prop Name | Status |
|-----------|-------------|-----------|--------|
| `BiomarkerCard` | `result.biomarker_intelligence` | `biomarkerProfile` | ✅ Working |
| `ResistanceCard` | `result.resistance_prediction` | `resistancePrediction` | ✅ Working |
| `DrugRankingCard` | `getDrugRanking()` (from `result.wiwfm`) | `drugRanking` | ✅ Working |
| `TrialMatchesCard` | `getTrials()` (from `result.trials`) | `trialMatches` | ✅ Working |
| `SOCRecommendationCard` | `result.soc_recommendation` | Spread props | ✅ Working |
| `NextTestCard` | `result.next_test_recommender` | `recommendations` | ✅ Working |
| `HintTilesPanel` | `getHintTiles()` (from `result.hint_tiles`) | `tiles` | ✅ Working |
| `MechanismChips` | `result.mechanism_map` | `mechanism_map` | ✅ Working |
| `ResistancePlaybook` | `result.resistance_playbook` | `resistance_playbook` | ✅ Working |
| `AyeshaSAEFeaturesCard` | `result.sae_features` | `sae_features` | ✅ Working |
| `ResistanceAlertBanner` | `result.resistance_alert` | `resistance_alert` | ✅ Working |
| `ToxicityRiskCard` | `result.toxicity_assessments` | `result` | ✅ Working |
| `SafetyGateCard` | `getPGxDrugs()` (from `result.wiwfm.drugs`) | `drug` | ✅ Working |
| `TrialSafetyGate` | `getTrials()` | `trial` | ✅ Working |

---

## 🔄 API Contracts & Data Formats

### **Orchestrator API: `/api/orchestrate/full`**

#### **Request Format:**
```typescript
{
  patient_id?: string,        // Auto-generated if not provided
  disease: string,            // Required: "ovarian", "myeloma", etc.
  mutations: Array<{
    gene: string,
    hgvs_p?: string,
    hgvs_c?: string,
    consequence?: string,
    chrom?: string,
    pos?: number,
    ref?: string,
    alt?: string,
    zygosity?: string
  }>,
  cytogenetics?: {            // MM-specific
    del_17p?: boolean,
    t4_14?: boolean,
    ...
  },
  treatment_line?: number,   // Default: 1
  prior_therapies?: string[],
  current_regimen?: string,
  current_drug_class?: string,
  clinical_data?: {...},
  run_async?: boolean,        // Default: false
  skip_agents?: string[]     // e.g., ["nutrition", "monitoring"]
}
```

#### **Response Format:**
```typescript
{
  patient_id: string,
  disease: string,
  phase: "initialized" | "extracting" | "analyzing" | "ranking" | "matching" | "planning" | "monitoring" | "complete" | "error",
  progress_percent: number,   // 0-100
  completed_agents: string[],
  created_at: string,         // ISO timestamp
  updated_at: string,         // ISO timestamp
  duration_ms?: number,
  mutation_count: number,
  mechanism_vector: number[], // 7D vector [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
  
  // ✅ Agent outputs (ALL AVAILABLE in orchestrator)
  biomarker_profile?: {
    tmb?: {
      value: number,
      classification: "TMB-H" | "TMB-M" | "TMB-L",
      threshold_used: number,
      mutation_count: number,
      total_mutations: number,
      exome_size_mb: number,
      method: string
    },
    msi?: {
      status: "MSI-H" | "MSS",
      dmmr_genes_mutated: string[],
      method: string,
      confidence: "high" | "moderate"
    },
    hrd?: {
      status: "HRD+" | "HRD-" | "HRD-uncertain",
      genes_mutated: string[],
      method: string,
      high_impact_detected: boolean,
      confidence: "high" | "moderate"
    },
    io_eligible: boolean,
    parp_eligible: boolean
  },
  resistance_prediction?: {
    risk_level: "HIGH" | "MEDIUM" | "LOW",
    probability: number,
    confidence: number,
    detected_genes: Array<{...}>,
    next_line_options?: {
      alternatives?: Array<{...}>,
      regimen_changes?: Array<{...}>,
      monitoring_changes?: {...}
    }
  },
  drug_ranking?: Array<{
    drug_name: string,
    drug_class?: string,
    moa?: string,
    efficacy_score: number,
    tier?: string,
    evidence_tier?: string,
    confidence?: number,
    mechanism?: string,
    rationale?: string[],
    badges?: string[]
  }>,
  synthetic_lethality_result?: {
    synthetic_lethality_detected: boolean,
    double_hit_description: string,
    essentiality_scores: Array<{...}>,
    broken_pathways: Array<{...}>,
    recommended_drugs: Array<{...}>,
    suggested_therapy: string
  },
  trial_matches?: Array<{
    nct_id: string,
    title: string,
    brief_summary?: string,
    phase?: string,
    status?: string,
    mechanism_fit_score?: number,
    eligibility_score?: number,
    combined_score?: number,
    why_matched?: string,
    url?: string,
    locations?: string[],
    contact?: {...}
  }>,
  nutrition_plan?: {
    patient_id: string,
    treatment: string,
    supplements: Array<{...}>,
    foods_to_prioritize: string[],
    foods_to_avoid: string[],
    drug_food_interactions: Array<{...}>,
    timing_rules: {...}
  },
  care_plan?: {
    patient_id: string,
    disease: string,
    generated_at: string,
    executive_summary: {...},
    sections: Array<{
      title: string,
      order: number,
      content: any
    }>,
    alerts: Array<{...}>
  },
  monitoring_config?: {
    frequency: "weekly" | "biweekly" | "monthly",
    biomarkers: string[],
    primary_biomarker: string,
    imaging: string,
    alerts_enabled: boolean,
    escalation_thresholds: {...},
    ctdna_monitoring: {...}
  },
  
  // ❌ NOT IN RESPONSE (but services exist)
  // toxicity_assessments?: {...},  // Service exists, not integrated
  // sae_features?: {...},          // Can be extracted, not stored in PatientState
  
  // Quality & alerts
  data_quality_flags: string[],
  alerts: Array<{
    id: string,
    alert_type: string,
    message: string,
    severity: "critical" | "warning" | "info",
    timestamp: string,
    source_agent: string,
    acknowledged: boolean
  }>
}
```

**Note:** To get full state including all agent outputs, use `GET /api/orchestrate/state/{patient_id}` which returns `PatientState.to_full_dict()` with all fields.

### **Status API: `/api/orchestrate/status/{patient_id}`**

#### **Response Format:**
```typescript
{
  patient_id: string,
  phase: string,
  progress_percent: number,
  current_agent?: string,      // Currently running agent
  completed_agents: string[],
  alerts: Array<{...}>,
  errors: string[],
  status_url: string,
  care_plan_url?: string
}
```

---

## 🔀 Data Flow Mapping

### **Current Flow (complete_care/v2):**
```
Frontend (UniversalCompleteCare.jsx)
  ↓ POST /api/complete_care/v2
Backend (ayesha_orchestrator_v2.py)
  ↓ Calls multiple services in parallel
  - /api/trials/universal
  - /api/soc/recommendation
  - /api/biomarker/intelligence
  - /api/efficacy/predict (WIWFM)
  - /api/care/resistance_playbook
  - /api/resistance/predict
  ↓ Aggregates responses
  ↓ Returns CompleteCareV2Response
Frontend
  ↓ Maps response to component props
  ↓ Renders components
```

### **Target Flow (orchestrator):**
```
Frontend (UniversalCompleteCare.jsx)
  ↓ POST /api/orchestrate/full
Backend (orchestrate.py)
  ↓ Orchestrator.run_full_pipeline()
  ↓ Phase 1: Data Extraction (if file)
  ↓ Phase 2: Parallel Analysis
    - Biomarker Agent
    - Resistance Agent
    - Nutrition Agent
  ↓ Phase 3: Drug Efficacy Ranking
  ↓ Phase 3.5: Synthetic Lethality
  ↓ Phase 4: Trial Matching
  ↓ Phase 5: Care Plan Generation
  ↓ Phase 6: Monitoring Setup
  ↓ Returns OrchestratePipelineResponse
Frontend
  ↓ Maps orchestrator response to component props
  ↓ Renders components
```

### **Data Transformation Required:**

#### **Orchestrator → Legacy Format Mapping:**
```javascript
const mapOrchestratorToLegacy = (orchestratorResponse) => {
  return {
    // Core mappings (all available in orchestrator)
    biomarker_intelligence: orchestratorResponse.biomarker_profile,
    resistance_prediction: orchestratorResponse.resistance_prediction,
    wiwfm: {
      drugs: orchestratorResponse.drug_ranking?.map(d => ({
        name: d.drug_name,
        drug_name: d.drug_name,
        moa: d.mechanism || d.drug_class || d.moa,
        efficacy_score: d.efficacy_score,
        confidence: d.confidence || d.efficacy_score,
        evidence_tier: d.tier || d.evidence_tier || 'insufficient',
        badges: d.badges || [],
        rationale: d.rationale || []
      })) || [],
      evidence_tier: orchestratorResponse.drug_ranking?.[0]?.evidence_tier || 'Supported',
      run_signature: orchestratorResponse.patient_id
    },
    trials: {
      trials: orchestratorResponse.trial_matches?.map(t => ({
        nct_id: t.nct_id,
        title: t.title,
        phase: t.phase,
        status: t.status,
        mechanism_fit_score: t.mechanism_fit_score,
        eligibility_score: t.eligibility_score,
        combined_score: t.combined_score,
        why_matched: t.why_matched,
        url: t.url,
        brief_summary: t.brief_summary,
        locations: t.locations,
        contact: t.contact
      })) || []
    },
    
    // Additional mappings (available in orchestrator)
    mechanism_map: {
      vector: orchestratorResponse.mechanism_vector || [0,0,0,0,0,0,0]
    },
    resistance_playbook: orchestratorResponse.resistance_prediction?.next_line_options || null,
    care_plan: orchestratorResponse.care_plan,
    nutrition_plan: orchestratorResponse.nutrition_plan, // Available but not in legacy format
    synthetic_lethality: orchestratorResponse.synthetic_lethality_result, // Available but not in legacy format
    monitoring_config: orchestratorResponse.monitoring_config, // Available but not in legacy format
    
    // Missing in orchestrator (need to handle gracefully)
    sae_features: orchestratorResponse.sae_features || null, // Not in PatientState currently
    soc_recommendation: null, // Not in orchestrator
    next_test_recommender: null, // Not in orchestrator
    hint_tiles: null, // Not in orchestrator
    toxicity_assessments: orchestratorResponse.toxicity_assessments || null, // Service exists, not integrated
    resistance_alert: null, // Not in orchestrator
    
    // Metadata
    provenance: {
      orchestrator: 'MOAT Orchestrator',
      run_id: orchestratorResponse.patient_id,
      generated_at: orchestratorResponse.updated_at,
      phase: orchestratorResponse.phase,
      completed_agents: orchestratorResponse.completed_agents || []
    },
    summary: {
      components_included: orchestratorResponse.completed_agents || [],
      ngs_status: orchestratorResponse.mutation_count > 0 ? 'available' : 'pending',
      confidence_level: 'moderate-high (70-90%)',
      pipeline_version: '2.0'
    }
  };
};
```

---

## 🧪 Test Requirements

### **Unit Tests**

#### **1. Data Transformation Tests**
- ✅ Test `mapOrchestratorToLegacy()` with full orchestrator response
- ✅ Test with partial responses (missing optional fields)
- ✅ Test with empty arrays/null values
- ✅ Test mechanism_vector mapping (7D vector)
- ✅ Test drug ranking format conversion
- ✅ Test trial matches format conversion

#### **2. Component Prop Mapping Tests**
- ✅ Test `BiomarkerCard` receives correct `biomarkerProfile` format
- ✅ Test `ResistanceCard` receives correct `resistancePrediction` format
- ✅ Test `DrugRankingCard` receives correct `drugRanking` format
- ✅ Test `TrialMatchesCard` receives correct `trialMatches` format
- ✅ Test `ResistancePlaybook` receives correct `resistance_playbook` format
- ✅ Test `AyeshaSAEFeaturesCard` receives correct `sae_features` format

#### **3. API Integration Tests**
- ✅ Test `/api/orchestrate/full` request format
- ✅ Test response parsing
- ✅ Test error handling (404, 500, network errors)
- ✅ Test status polling (`/api/orchestrate/status/{patient_id}`)
- ✅ Test file upload (multipart/form-data)

### **Integration Tests**

#### **1. End-to-End Pipeline Tests**
- ✅ Test complete pipeline execution from frontend
- ✅ Test status polling during execution
- ✅ Test component rendering with real orchestrator data
- ✅ Test error recovery (agent failures)

#### **2. Backward Compatibility Tests**
- ✅ Test legacy endpoint still works (if kept as fallback)
- ✅ Test feature flag toggling between endpoints
- ✅ Test gradual migration path

### **Component Tests**

#### **1. UniversalCompleteCare Tests**
- ✅ Test initial render
- ✅ Test patient profile loading
- ✅ Test plan generation button click
- ✅ Test loading state display
- ✅ Test error state display
- ✅ Test result display with all components
- ✅ Test file upload integration
- ✅ Test status polling integration

#### **2. PipelineStatusCard Tests** (NEW)
- ✅ Test phase display
- ✅ Test progress bar
- ✅ Test alerts display
- ✅ Test agent execution status
- ✅ Test auto-stop when complete

#### **3. PatientUpload Tests** (EXISTS)
- ✅ Test file selection
- ✅ Test file type detection
- ✅ Test upload trigger
- ✅ Test error handling

---

## 🔍 Gap Analysis

### **Agent Status (Actual Implementation vs Documentation):**

**CRITICAL FINDING:** The master index (`01_CURRENT_STATE.md`) is **OUTDATED**. Actual implementation is **100% complete** for core agents.

#### **✅ FULLY INTEGRATED AGENTS:**

1. **Data Extraction Agent** (`01`)
   - **Documentation Says:** ⏳ SKELETON
   - **Actual Status:** ✅ **INTEGRATED**
   - **Location:** `api/services/extraction/extraction_agent.py`
   - **Wired In:** `_run_extraction_phase()` (line 201-254)
   - **Capabilities:** VCF, PDF, MAF, JSON, TXT parsing with quality validation
   - **Output:** `patient_profile`, `mutations`, `disease`

2. **Biomarker Agent** (`02`)
   - **Status:** ✅ **INTEGRATED**
   - **Location:** `_run_biomarker_agent()` (line 453-547)
   - **Output:** TMB, MSI, HRD with enhanced gene panels

3. **Resistance Agent** (`03`)
   - **Status:** ✅ **INTEGRATED & VALIDATED**
   - **Location:** `_run_resistance_agent()` (line 549-626)
   - **Services:** `ResistanceProphetService` + `ResistancePlaybookService`
   - **Output:** Risk level, detected genes, next_line_options

4. **Drug Efficacy Agent** (`04`)
   - **Documentation Says:** ⏳ SKELETON
   - **Actual Status:** ✅ **INTEGRATED**
   - **Location:** `_run_drug_efficacy_agent()` (line 743-819)
   - **Service:** `EfficacyOrchestrator` (direct import, no HTTP)
   - **Output:** Ranked drugs, mechanism_vector (7D), evidence_tier
   - **Note:** SAE features can be extracted if `include_sae_features: true` in options

5. **Trial Matching Agent** (`05`)
   - **Status:** ✅ **INTEGRATED**
   - **Location:** `_run_trial_matching_agent()` (line 923-989)
   - **Service:** `TrialMatchingAgent`
   - **Output:** Matched trials with mechanism_fit_score, eligibility_score, combined_score

6. **Nutrition Agent** (`06`)
   - **Documentation Says:** ⏳ SKELETON
   - **Actual Status:** ✅ **INTEGRATED**
   - **Location:** `_run_nutrition_agent()` (line 660-741)
   - **Service:** `NutritionAgent` (direct import)
   - **Output:** Supplements, foods_to_prioritize, foods_to_avoid, drug_food_interactions

7. **Care Plan Agent** (`07`)
   - **Status:** ✅ **INTEGRATED**
   - **Location:** `_run_care_plan_agent()` (line 991-1135)
   - **Output:** 8-section care plan with executive summary

8. **Monitoring Agent** (`08`)
   - **Status:** ✅ **INTEGRATED**
   - **Location:** `_run_monitoring_agent()` (line 1137-1221)
   - **Output:** Disease-specific monitoring config with ctDNA support

9. **Synthetic Lethality Agent** (`14`)
   - **Status:** ✅ **INTEGRATED**
   - **Location:** `_run_synthetic_lethality_agent()` (line 821-921)
   - **Service:** `SyntheticLethalityAgent`
   - **Output:** SL detection, essentiality scores, broken pathways, recommended drugs

#### **⚠️ PARTIALLY INTEGRATED:**

10. **Trigger System** (`09`)
    - **Status:** ✅ **EXISTS** but ⚠️ **NOT IN PIPELINE**
    - **Location:** `api/services/triggers/trigger_engine.py`
    - **Capabilities:** 8 trigger types, 13 action handlers
    - **Issue:** `TriggerEngine` exists but is NOT called in any pipeline phase
    - **Required:** Add `_run_trigger_system_phase()` to orchestrator
    - **Impact:** Event-driven automation not active

#### **❌ NOT INTEGRATED (But Services Exist):**

11. **Toxicity Risk Agent** (`16`)
    - **Status:** ❌ **NOT IN PIPELINE** (but service exists)
    - **Service Location:** `api/services/safety_service.py` ✅ Complete
    - **Plan:** `.cursor/MOAT/orchestration/16_TOXICITY_RISK_AGENT_INTEGRATION_PLAN.md`
    - **Required:** Create `toxicity_agent.py` and wire to `_run_analysis_phase()`
    - **Impact:** `ToxicityRiskCard` won't display in orchestrator response

12. **SAE Features** (Not a separate agent, but capability)
    - **Status:** ⚠️ **SERVICE EXISTS** but ⚠️ **NOT IN ORCHESTRATOR OUTPUT**
    - **Service Location:** `api/services/sae_service.py` ✅ Complete
    - **Integration Point:** Can be extracted in `EfficacyOrchestrator` if `include_sae_features: true`
    - **Current:** SAE features are computed but not stored in `PatientState`
    - **Required:** Add `sae_features` field to `PatientState` and extract in drug efficacy phase
    - **Impact:** `AyeshaSAEFeaturesCard` will show placeholder

#### **❌ MISSING (No Service/Agent):**

13. **SOC Recommendation**
    - **Status:** ❌ Missing
    - **Required:** Add SOC agent or call existing `/api/soc/recommendation` endpoint
    - **Impact:** `SOCRecommendationCard` won't display

14. **Next Test Recommender**
    - **Status:** ❌ Missing
    - **Required:** Add next test recommender agent
    - **Impact:** `NextTestCard` won't display

15. **Hint Tiles**
    - **Status:** ❌ Missing
    - **Required:** Add hint tiles generation to care plan phase
    - **Impact:** `HintTilesPanel` won't display

16. **PGx Safety Gate**
    - **Status:** ❌ Missing
    - **Required:** Add PGx screening to drug efficacy phase or separate agent
    - **Impact:** `SafetyGateCard` and `TrialSafetyGate` won't display

### **Resistance Playbook Status:**

- **Status:** ⚠️ **PARTIALLY AVAILABLE**
- **Location:** In `resistance_prediction.next_line_options` from `ResistancePlaybookService`
- **Format:** Contains `alternatives`, `regimen_changes`, `monitoring_changes`
- **Missing:** Full playbook format with `risks`, `combo_strategies`, `next_line_switches` arrays
- **Impact:** `ResistancePlaybook` component expects full format - may need data transformation or service enhancement

### **Data Format Mismatches:**

1. **Drug Ranking Format**
   - **Orchestrator:** `drug_ranking` array with `drug_name`, `drug_class`, `efficacy_score`
   - **Legacy:** `wiwfm.drugs` array with `name`, `moa`, `efficacy_score`, `confidence`
   - **Fix:** Transformation function handles this

2. **Trial Matches Format**
   - **Orchestrator:** `trial_matches` array with `nct_id`, `title`, `mechanism_fit_score`
   - **Legacy:** `trials.trials` array with same fields but nested
   - **Fix:** Transformation function handles this

3. **Biomarker Format**
   - **Orchestrator:** `biomarker_profile` with `tmb`, `msi`, `hrd` objects
   - **Legacy:** `biomarker_intelligence` with same structure
   - **Fix:** Direct mapping works

4. **Resistance Format**
   - **Orchestrator:** `resistance_prediction` with `risk_level`, `probability`, `detected_genes`
   - **Legacy:** `resistance_prediction` with similar structure
   - **Fix:** Direct mapping works

---

## 📝 Implementation Checklist

### **Phase 1: Core Integration** 🔴 **CRITICAL**

- [ ] **1.1** Update `UniversalCompleteCare.jsx` to use `/api/orchestrate/full`
- [ ] **1.2** Create `mapOrchestratorToLegacy()` transformation function
- [ ] **1.3** Test data transformation with real orchestrator responses
- [ ] **1.4** Update component prop mappings
- [ ] **1.5** Add error handling for orchestrator API
- [ ] **1.6** Test backward compatibility (if keeping legacy endpoint)
- [ ] **1.7** Verify all 9 integrated agents work correctly (Data Extraction, Biomarker, Resistance, Drug Efficacy, Nutrition, Trial Matching, Synthetic Lethality, Care Plan, Monitoring)

### **Phase 2: Status Polling** 🔴 **CRITICAL**

- [ ] **2.1** Create `usePipelineStatus` hook (DONE)
- [ ] **2.2** Create `PipelineStatusCard` component
- [ ] **2.3** Integrate status polling into `UniversalCompleteCare`
- [ ] **2.4** Test polling stops when complete
- [ ] **2.5** Test error handling during polling

### **Phase 3: File Upload** 🔴 **HIGH**

- [ ] **3.1** Integrate `PatientUpload` component
- [ ] **3.2** Test file upload with orchestrator API
- [ ] **3.3** Test file type detection (VCF, PDF, MAF)
- [ ] **3.4** Test error handling for invalid files

### **Phase 4: Missing Components** 🟡 **MEDIUM**

- [ ] **4.1** Verify `ResistancePlaybook` works with orchestrator data (check `resistance_prediction.next_line_options` format)
- [ ] **4.2** Add SAE features extraction to drug efficacy phase and store in `PatientState`
- [ ] **4.3** Integrate Toxicity Risk Agent (service exists, needs orchestrator wiring)
- [ ] **4.4** Add Trigger System phase to pipeline (TriggerEngine exists but not called)
- [ ] **4.5** Document missing components (SOC, Next Test, Hint Tiles, PGx)

### **Phase 5: Testing** 🟢 **ONGOING**

- [ ] **5.1** Write unit tests for data transformation
- [ ] **5.2** Write integration tests for API calls
- [ ] **5.3** Write component tests for `UniversalCompleteCare`
- [ ] **5.4** Write E2E tests for complete flow
- [ ] **5.5** Test error scenarios

---

## 🎯 Success Criteria

### **Phase 1 Complete When:**
- ✅ `UniversalCompleteCare` uses `/api/orchestrate/full`
- ✅ All existing components render correctly with orchestrator data
- ✅ Data transformation works for all component props
- ✅ Error handling works correctly

### **Phase 2 Complete When:**
- ✅ Status polling shows real-time progress
- ✅ `PipelineStatusCard` displays correctly
- ✅ Polling stops when pipeline completes

### **Phase 3 Complete When:**
- ✅ File upload works with orchestrator
- ✅ File types are detected correctly
- ✅ Upload triggers pipeline execution

### **Overall Success:**
- ✅ All existing functionality preserved
- ✅ No regressions in component rendering
- ✅ All tests pass
- ✅ Documentation updated

---

---

## 📊 Implementation Readiness Summary

### **✅ Ready for Integration (9 Agents):**

1. ✅ **Data Extraction** - Fully integrated, supports VCF/PDF/MAF/JSON/TXT
2. ✅ **Biomarker** - TMB, MSI, HRD calculation with enhanced gene panels
3. ✅ **Resistance** - ResistanceProphet + ResistancePlaybook integrated
4. ✅ **Drug Efficacy** - EfficacyOrchestrator with S/P/E framework
5. ✅ **Trial Matching** - TrialMatchingAgent with mechanism fit scoring
6. ✅ **Nutrition** - NutritionAgent with toxicity-aware recommendations
7. ✅ **Synthetic Lethality** - Full SL/Essentiality analysis
8. ✅ **Care Plan** - 8-section unified care plan generation
9. ✅ **Monitoring** - Disease-specific monitoring configuration

### **⚠️ Needs Integration (2 Services):**

1. ⚠️ **Toxicity Risk** - Service exists (`safety_service.py`), needs orchestrator wiring
2. ⚠️ **SAE Features** - Service exists (`sae_service.py`), needs extraction in drug efficacy phase and storage in PatientState

### **❌ Missing (4 Components):**

1. ❌ **SOC Recommendation** - No service/agent
2. ❌ **Next Test Recommender** - No service/agent
3. ❌ **Hint Tiles** - No service/agent
4. ❌ **PGx Safety Gate** - No service/agent

### **🔧 Needs Pipeline Integration (1 System):**

1. ⚠️ **Trigger System** - `TriggerEngine` exists but not called in any pipeline phase

---

**Next Steps:** Review this analysis, then proceed with implementation following the checklist above.

