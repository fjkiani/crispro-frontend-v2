# 🏆 MOAT Agent Ownership Registry

**Purpose:** Definitive assignment of agent responsibilities  
**Created:** January 28, 2025  
**Status:** ✅ ACTIVE

---

## 🆕 NEW AGENTS - READ FIRST

Before starting your assignment, read these implementation guides:

- **[⚡ Quick Reference](./QUICK_REFERENCE_AGENT_QUESTIONS.md)** - Answers 3 common questions (START HERE)
- **[📖 Agent Implementation Guide](./AGENT_IMPLEMENTATION_GUIDE.md)** - Detailed patterns & examples
- **[🗺️ Master Index](./orchestration/00_MASTER_INDEX.mdc)** - Navigate all modules

**Common Questions:**
1. ❓ What should I build? → Check your agent type below + read your MDC
2. ❓ Call endpoints or use services? → Always use services directly
3. ❓ Modify care plan? → No - it auto-consumes your output

---

## 📊 OWNERSHIP MATRIX

### 🔴 CRITICAL PATH (Must Have)

| Agent ID | Module | Owner | Status | Dependencies |
|----------|--------|-------|--------|--------------|
| `LEAD_ORCHESTRATOR` | 10 - State Management | **Senior Agent** | ✅ COMPLETE | None |
| `AGENT_01_EXTRACTOR` | 01 - Data Extraction | **Auto (JR Agent D)** | ✅ **COMPLETE** | None |
| `AGENT_03_RESISTANCE` | 03 - Resistance Prediction | **Senior Agent** | ✅ VALIDATED | 01, 02 |
| `AGENT_11_API` | 11 - API Contracts | **Senior Agent** | ✅ COMPLETE | All |

### 🟡 HIGH PRIORITY

| Agent ID | Module | Owner | Status | Dependencies |
|----------|--------|-------|--------|--------------|
| `AGENT_02_BIOMARKER` | 02 - Biomarker Calculation | **JR Agent B** | ✅ INTEGRATED | 01 |
| `AGENT_04_EFFICACY` | 04 - Drug Efficacy (S/P/E) | **JR Agent C** | ⏳ PENDING | 01, 02, 03 |
| `AGENT_05_TRIALS` | 05 - Trial Matching | **JR Agent D** | ✅ **COMPLETE** | 01, 02, 04 |
| `AGENT_07_CAREPLAN` | 07 - Care Plan Generation | **Senior Agent** | ✅ INTEGRATED | 01-06 |
| `AGENT_08_MONITOR` | 08 - Continuous Monitoring | **JR Agent E** | ⏳ SKELETON | 07 |
| `AGENT_09_TRIGGER` | 09 - Event Trigger System | **Auto (JR Agent D)** | ✅ **COMPLETE** | All |

### 🟢 MEDIUM PRIORITY

| Agent ID | Module | Owner | Status | Dependencies |
|----------|--------|-------|--------|--------------|
| `AGENT_06_NUTRITION` | 06 - Nutrition Planning | **JR Agent G** | ⏳ SKELETON | 01 |
| `AGENT_12_UI` | 12 - UI Dashboard | **Frontend Agent** | ⬜ TODO | 11 |
| `AGENT_13_SECURITY` | 13 - Security & Compliance | **Security Agent** | ⬜ TODO | All |

---

## 📋 DETAILED AGENT ASSIGNMENTS

### `AGENT_01_EXTRACTOR` - Data Extraction Agent

**Owner:** Auto (JR Agent D)  
**Status:** ✅ **COMPLETE**  
**Priority:** 🔴 CRITICAL

#### What Was Built ✅
- `api/services/extraction/extraction_agent.py` - Main agent class
- `api/services/extraction/parsers/` - VCF, MAF, PDF, JSON, TXT parsers
- `api/services/extraction/models.py` - PatientProfile, Mutation models
- Wired to orchestrator in `_run_extraction_phase()`

#### Responsibilities
- Parse VCF files (standard genomic format)
- Parse MAF files (mutation annotation format)
- Parse NGS PDFs (LLM-based extraction)
- Normalize gene symbols to HGNC
- Validate HGVS notation
- Flag data quality issues

#### Files to Create
```
api/services/extraction/
├── __init__.py
├── extraction_agent.py      # Main agent class
├── parsers/
│   ├── vcf_parser.py        # VCF 4.x parsing
│   ├── maf_parser.py        # MAF parsing
│   ├── pdf_parser.py        # LLM-based PDF extraction
│   └── json_parser.py       # JSON input handling
├── validators/
│   ├── gene_validator.py    # HGNC validation
│   └── hgvs_validator.py    # HGVS notation
└── tests/
    └── test_extraction.py
```

#### Acceptance Criteria
- [ ] Parse VCF with >95% accuracy
- [ ] Parse PDFs with >85% mutation extraction rate
- [ ] All genes normalized to HGNC symbols
- [ ] Processing time <10 seconds per file

#### Interface Contract
```python
class DataExtractionAgent:
    async def extract(
        self,
        file: BinaryIO,
        file_type: str,  # 'vcf', 'maf', 'pdf', 'json'
        metadata: Dict = None
    ) -> PatientProfile:
        """Return structured PatientProfile with mutations."""
```

---

### `AGENT_02_BIOMARKER` - Biomarker Calculation Agent

**Owner:** JR Agent B  
**Status:** ✅ INTEGRATED (basic), needs enhancement  
**Priority:** 🔴 CRITICAL

#### Current State
Basic TMB/MSI/HRD calculation exists in orchestrator. Needs:
- Validated TMB calculation (r=0.933 with TCGA)
- dMMR gene panel MSI inference
- HRD score calculation

#### Responsibilities
- Calculate TMB (mut/Mb) with validated method
- Infer MSI status from dMMR genes
- Calculate/infer HRD status
- Determine IO/PARP eligibility

#### Files to Enhance
```
api/services/biomarkers/
├── __init__.py
├── biomarker_agent.py       # Main agent
├── tmb_calculator.py        # Validated TMB (use existing script)
├── msi_calculator.py        # MSI inference
├── hrd_calculator.py        # HRD scoring
└── tests/
```

#### Acceptance Criteria
- [ ] TMB: Pearson r ≥ 0.90 with TCGA ground truth
- [ ] MSI: Correctly identify dMMR from gene panel
- [ ] HRD: BRCA1/2 → HRD+, HRR genes → HRD-inferred

---

### `AGENT_03_RESISTANCE` - Resistance Prediction Agent

**Owner:** Senior Agent  
**Status:** ✅ VALIDATED  
**Priority:** 🔴 CRITICAL

#### What's Done
- ✅ DIS3 RR=2.08, p=0.0145 (SIGNIFICANT)
- ✅ TP53 RR=1.90, p=0.11 (TREND)
- ✅ NF1 RR=2.10 (OV)
- ✅ KRAS RR=1.97 (OV)
- ✅ ResistancePlaybookService (next-line options)
- ✅ Cytogenetics support (del(17p), t(4;14), 1q)
- ✅ Treatment line context (1.2x/1.4x multipliers)
- ✅ Downstream agent handoffs

#### Files
```
api/services/
├── resistance_prophet_service.py   # ✅ Main prediction
├── resistance_playbook_service.py  # ✅ Next-line options
```

#### DO NOT MODIFY (Validated)
- MM_HIGH_RISK_GENES dictionary
- MM_CYTOGENETICS dictionary
- OV_RESISTANCE_GENES dictionary
- Risk stratification thresholds

---

### `AGENT_04_EFFICACY` - Drug Efficacy Agent (S/P/E)

**Owner:** JR Agent C  
**Status:** ⏳ PENDING  
**Priority:** 🟡 HIGH

#### Responsibilities
- **S** (Sequence): Evo2 variant scoring
- **P** (Pathway): Pathway burden alignment
- **E** (Evidence): PubMed/ClinVar synthesis
- Rank drugs with confidence scores
- Generate 7D mechanism vector

#### Files to Create
```
api/services/efficacy/
├── __init__.py
├── drug_efficacy_agent.py      # Main agent
├── sequence_scorer.py          # Evo2 integration
├── pathway_scorer.py           # Pathway burden
├── evidence_scorer.py          # Literature synthesis
├── mechanism_vector.py         # 7D vector generation
└── tests/
```

#### Interface Contract
```python
class DrugEfficacyAgent:
    async def rank_drugs(
        self,
        patient_profile: Dict,
        biomarker_profile: Dict,
        resistance_prediction: Dict,
        drug_classes: List[str] = None
    ) -> DrugEfficacyResult:
        """Return ranked drugs with S/P/E breakdown."""
```

#### Acceptance Criteria
- [ ] S/P/E scores computed for each drug
- [ ] 7D mechanism vector generated
- [ ] Response time <5 seconds

---

### `AGENT_05_TRIALS` - Trial Matching Agent

**Owner:** JR Agent D  
**Status:** ✅ **COMPLETE**  
**Priority:** 🟡 HIGH

#### What Was Built ✅
- `api/services/trials/trial_matching_agent.py` - Main agent class
- Wired `AutonomousTrialAgent` for query generation and search
- Wired `MechanismFitRanker` for mechanism-based ranking (Manager P4)
- Wired `TrialDataEnricher` for MoA extraction and eligibility
- Integrated with orchestrator in `orchestrator.py`

#### Files Created
```
api/services/trials/
├── __init__.py                    ✅
├── trial_matching_agent.py        ✅
└── README.md                      ✅
```

#### Acceptance Criteria ✅
- [x] Wire existing services to orchestrator
- [x] Return top 10 trials with mechanism fit
- [x] Include eligibility breakdown
- [x] Manager P4 compliance (alpha=0.7, beta=0.3, thresholds)
- [x] Manager P3 compliance (Gemini tags preferred, runtime fallback)

---

### `AGENT_14_SL_ESSENTIALITY` - Synthetic Lethality & Essentiality Agent

**Owner:** AI Agent (Synthetic Lethality Specialist)  
**Status:** ⏳ PENDING  
**Priority:** 🟡 HIGH

#### Responsibilities
- Score gene essentiality using Evo2 foundation model
- Map broken pathways from mutations
- Identify essential backup pathways (synthetic lethality relationships)
- Recommend drugs targeting essential backups
- Generate AI explanations (3 audiences: clinician/patient/researcher)
- Validate with benchmarks (50% drug match, 100% Evo2 usage)

#### Files to Create
```
api/services/synthetic_lethality/
├── __init__.py
├── sl_agent.py                    # Main orchestrating agent
├── essentiality_scorer.py         # Evo2 integration
├── pathway_mapper.py              # Pathway disruption mapping
├── dependency_identifier.py       # Essential backup identification
├── drug_recommender.py            # Drug recommendation engine
├── explanation_generator.py       # AI explanation generation
├── constants.py                   # Pathways, genes, thresholds
├── models.py                      # Data models
└── tests/
    ├── test_essentiality.py
    ├── test_pathways.py
    ├── test_recommendations.py
    └── test_integration.py
```

#### Interface Contract
```python
class SyntheticLethalityAgent:
    async def analyze(
        self,
        request: SyntheticLethalityRequest
    ) -> SyntheticLethalityResult:
        """Return complete synthetic lethality analysis with essentiality scores."""
```

#### Acceptance Criteria
- [ ] Gene essentiality scored with Evo2 integration
- [ ] Broken pathways identified from mutations
- [ ] Essential backup pathways determined via SL relationships
- [ ] Drugs recommended targeting essential pathways
- [ ] AI explanations generated for 3 audiences
- [ ] Integrated with orchestrator
- [ ] Benchmark validated (50% drug match, 100% Evo2 usage)
- [ ] Unit test coverage >80%

#### Reference
- **MDC File:** `.cursor/MOAT/orchestration/14_SYNTHETIC_LETHALITY_ESSENTIALITY_AGENT.mdc` (1,725 lines)
- **Frontend:** Already built (`oncology-frontend/src/components/SyntheticLethality/`)
- **Backend Endpoint:** Existing `/api/guidance/synthetic_lethality` (needs MOAT integration)

---

### `AGENT_06_NUTRITION` - Nutrition Planning Agent

**Owner:** JR Agent G  
**Status:** ⏳ SKELETON  
**Priority:** 🟢 MEDIUM

#### Responsibilities
- Drug-food interaction warnings
- Timing recommendations
- Germline-based dietary considerations

#### Files to Create
```
api/services/nutrition/
├── __init__.py
├── nutrition_agent.py
├── drug_food_interactions.py
└── tests/
```

---

### `AGENT_07_CAREPLAN` - Care Plan Generation Agent

**Owner:** Senior Agent  
**Status:** ✅ INTEGRATED  
**Priority:** 🟡 HIGH

#### What's Done
- Aggregates all agent outputs
- Generates unified care plan document
- Includes sections for each agent output

#### Files
```
api/services/orchestrator/orchestrator.py
# _run_care_plan_agent method
```

---

### `AGENT_08_MONITOR` - Continuous Monitoring Agent

**Owner:** JR Agent E  
**Status:** ⏳ SKELETON  
**Priority:** 🟡 HIGH

#### Responsibilities
- Configure monitoring frequency by risk
- CA-125 kinetics tracking
- ctDNA monitoring alerts
- RECIST response tracking

#### Files to Create
```
api/services/monitoring/
├── __init__.py
├── monitoring_agent.py
├── ca125_tracker.py
├── ctdna_tracker.py
└── tests/
```

---

### `AGENT_09_TRIGGER` - Event Trigger System

**Owner:** Auto (JR Agent D)  
**Status:** ✅ **COMPLETE**  
**Priority:** 🟡 HIGH

#### What Was Built ✅
- `api/services/triggers/trigger_engine.py` - Main trigger engine
- `api/services/triggers/models.py` - TriggerResult, TriggerSeverity
- 8 trigger types implemented
- 13 action handlers implemented
- Wired to orchestrator in `process_event()`

#### Responsibilities
- Event detection (resistance, TMB-H, PD)
- Automated actions (alerts, re-ranking)
- Escalation protocols

#### Files to Create
```
api/services/triggers/
├── __init__.py
├── trigger_engine.py
├── event_handlers/
│   ├── resistance_handler.py
│   ├── tmb_handler.py
│   └── pd_handler.py
└── tests/
```

---

## 📡 INTER-AGENT COMMUNICATION

### Message Types

```python
class MessageType(Enum):
    REQUEST = "request"      # Agent-to-agent request
    RESPONSE = "response"    # Response to request
    EVENT = "event"          # Event notification
    ALERT = "alert"          # Clinical alert
    BROADCAST = "broadcast"  # Broadcast to all
```

### Event Types

| Event | Publisher | Subscribers |
|-------|-----------|-------------|
| `resistance_detected` | AGENT_03 | AGENT_04, AGENT_07, AGENT_08 |
| `tmb_high` | AGENT_02 | AGENT_04, AGENT_05 |
| `hrd_positive` | AGENT_02 | AGENT_04, AGENT_05 |
| `drug_ranked` | AGENT_04 | AGENT_05, AGENT_07 |
| `trial_matched` | AGENT_05 | AGENT_07 |
| `care_plan_ready` | AGENT_07 | AGENT_08 |

---

## 🔄 DEPENDENCY GRAPH

```
AGENT_01_EXTRACTOR (no deps)
     │
     ├──────────────────────────────────────┐
     │                                      │
     ▼                                      ▼
AGENT_02_BIOMARKER                    AGENT_06_NUTRITION
     │
     ▼
AGENT_03_RESISTANCE
     │
     ▼
AGENT_04_EFFICACY
     │
     ├────────────────┐
     ▼                ▼
AGENT_05_TRIALS    (mechanism vector)
     │
     ▼
AGENT_07_CAREPLAN ◄── aggregates all
     │
     ▼
AGENT_08_MONITOR
     │
     ▼
AGENT_09_TRIGGER
```

---

## 📋 TASK PICKUP PROCESS

### For JR Agents

1. **Claim a task** by updating this registry
2. **Read the MDC file** for your module (`.cursor/MOAT/orchestration/XX_MODULE.mdc`)
3. **Follow the interface contract** exactly
4. **Write tests** before implementing
5. **Update status** when complete

### Status Legend

| Status | Meaning |
|--------|---------|
| ⬜ TODO | Not started |
| ⏳ PENDING | Assigned, in progress |
| ⏳ SKELETON | Placeholder exists |
| ✅ INTEGRATED | Wired to orchestrator |
| ✅ VALIDATED | Tested with real data |
| ✅ COMPLETE | Production ready |

---

## 📞 ESCALATION

If blocked, escalate to Senior Agent with:
- Agent ID
- Blocking issue
- Attempted solutions
- Proposed resolution

---

**Registry Status:** ✅ ACTIVE  
**Last Updated:** January 28, 2025

