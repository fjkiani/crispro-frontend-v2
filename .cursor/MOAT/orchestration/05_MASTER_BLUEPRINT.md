# Orchestration System: Master Blueprint

**Source:** Consolidated from `ULTIMATE_MOAT_ORCHESTRATION.mdc`  
**Date:** January 28, 2025  
**Status:** ✅ **REFERENCE** - Complete workflow and architecture  
**Note:** This is a summary. See original file for full details.

---

## 🎯 THE VISION

> **"Upload once. Track forever. Never miss a signal."**

A patient's genomic data enters the system once. From that moment, an agentic swarm takes over - calculating biomarkers, predicting resistance, matching trials, generating nutrition plans, and continuously monitoring for signals that require action. The oncologist receives actionable intelligence, not raw data.

---

## 📊 THE COMPLETE MOAT STACK

### Layer 1: Validated Biomarker Capabilities ✅

| Capability | Validation | Status |
|------------|------------|--------|
| TMB Calculation | r=0.933, 95.4% accuracy | ✅ Production |
| TMB-H Classification | Validated against TCGA | ✅ Production |
| MAPK Resistance | RR=1.97, p<0.05 | ✅ Production |
| MSI Detection | Logic ready | ⏳ Need data |
| HRD Inference | MBD4 + BRCA logic | ✅ Production |

### Layer 2: Clinical Intelligence Capabilities ✅

| Capability | Source | Status |
|------------|--------|--------|
| Toxicity-Nutrition | Drug MoA + Food DB | ✅ Production |
| Trial Matching | ClinicalTrials.gov API | ✅ Production |
| SOC Recommendations | NCCN Guidelines | ✅ Production |
| Drug Efficacy (S/P/E) | Evo2 + Pathway + Evidence | ✅ Production |

### Layer 2.5: VUS Resolution

**Goal:** Convert "unknown variant" → "understood + actionable (or ignorable)" with receipts.

**Canonical endpoint:** `POST /api/vus/identify`

**Two valid resolution paths:**
- **Resolved by prior:** ClinVar decisive (Pathogenic/Likely pathogenic/Benign/Likely benign)
- **Resolved by Evo2 (ML):** ClinVar non-decisive/missing and Evo2 signal is decisive

### Layer 3: Agentic Monitoring Capabilities 🔄

| Capability | Trigger | Status |
|------------|---------|--------|
| CA-125 Kinetics | >25% rise | 🔄 Building |
| ctDNA Monitoring | Variant reappearance | 🔄 Building |
| Treatment Response | RECIST criteria | 🔄 Building |
| Adverse Events | CTCAE grading | 🔄 Building |

---

## 🔄 THE MASTER WORKFLOW

### Phase 0: Data Upload & Ingestion

**User Action:** Uploads NGS PDF, VCF, MAF, or enters mutations manually

**Components:**
- File upload portal
- Manual entry form
- Data validation

### Phase 1: Agentic Data Extraction

**Trigger:** File upload detected

**Actions:**
1. Detect file type (PDF, VCF, MAF, JSON)
2. Parse mutations (gene, variant, VAF, coverage)
3. Extract clinical data (stage, histology, biomarkers)
4. Extract demographics (age, sex, ECOG)
5. Validate data quality (missing fields, coverage thresholds)
6. Flag ambiguities for human review

**Output:** Structured PatientProfile object

### Phase 2: Biomarker Calculation Pipeline

**Trigger:** PatientProfile extracted

**Parallel Calculations:**
- TMB Calculation (r=0.933 vs TCGA)
- MSI Detection (dMMR genes)
- HRD Inference (MBD4 + BRCA logic)
- IO Eligibility Determination

**Output:** BiomarkerProfile

### Phase 3: Resistance Prediction

**Trigger:** BiomarkerProfile calculated

**Analyses:**
- MAPK pathway analysis (RR=1.97, p<0.05)
- DDR pathway analysis
- Platinum sensitivity prediction
- PARP sensitivity prediction

**Output:** ResistancePrediction

### Phase 4: Drug Efficacy Ranking (S/P/E)

**Trigger:** BiomarkerProfile + ResistancePrediction available

**Framework:**
- **S (Sequence):** Evo2 variant impact scoring
- **P (Pathway):** Drug-pathway alignment
- **E (Evidence):** Literature + ClinVar synthesis
- **Formula:** `0.3*S + 0.4*P + 0.3*E + clinvar_prior`

**Output:** Ranked drug list with confidence tiers

### Phase 5: Trial Matching

**Trigger:** Drug ranking + BiomarkerProfile available

**Process:**
- 7D Mechanism Vector construction
- ClinicalTrials.gov API search
- Mechanism fit ranking
- Eligibility scoring

**Output:** Matched trials with scores

### Phase 6: Toxicity-Aware Nutrition

**Trigger:** Drug ranking available

**Process:**
- Drug → Toxicity pathway mapping
- Mitigating food recommendations
- Drug-food interaction checking
- Timing rules (when to take supplements)

**Output:** NutritionPlan

### Phase 7: Unified Care Plan Generation

**Trigger:** All agent outputs available

**Process:**
- Aggregates all agent outputs
- Generates unified care plan document
- Includes: Treatment recommendations, resistance monitoring, trials, nutrition

**Output:** UnifiedCarePlan

### Phase 8: Continuous Monitoring & Alerts

**Trigger:** Care plan generated

**Process:**
- CA-125 kinetics tracking setup
- ctDNA monitoring configuration
- Treatment response (RECIST) tracking
- Risk-based frequency recommendations

**Output:** MonitoringConfig

---

## 🚨 EVENT TRIGGERS & AUTOMATED RESPONSES

### Trigger System Architecture

**Event Detection:**
- Resistance signals (MAPK activation, DDR deficiency)
- Biomarker changes (TMB-H, MSI-H, HRD status)
- Treatment response (PD, PR, SD)
- Adverse events (CTCAE grading)

**Automated Actions:**
- Alerts to oncologist
- Re-ranking drugs
- Re-matching trials
- Escalation protocols

---

## 🏗️ AGENTIC ARCHITECTURE

### Agent Hierarchy

```
ORCHESTRATOR (Coordinator)
├── Data Extraction Agent
├── Biomarker Agent
├── Resistance Agent
├── Drug Efficacy Agent
├── Trial Matching Agent
├── Nutrition Agent
├── Care Plan Agent
├── Monitoring Agent
└── Trigger Engine
```

### Agent Communication Protocol

- **State-based:** All agents read/write to PatientState
- **Event-driven:** Agents trigger on state changes
- **Parallel execution:** Where dependencies allow
- **Sequential execution:** For dependent agents

### State Management

- **Single source of truth:** PatientState object
- **Full audit trail:** All state changes logged
- **History tracking:** For debugging and analysis

---

## 🖥️ USER INTERFACE INTEGRATION

### Dashboard Layout

**Tabs:**
1. **Analysis Tab** - All analysis cards (Biomarker, Resistance, Drug Ranking, Trials, Nutrition, Synthetic Lethality)
2. **Care Plan Tab** - Full care plan viewer
3. **Monitoring Tab** - Monitoring dashboard with alerts

**Components:**
- PatientUpload (file upload)
- Analysis cards (lazy-loaded)
- CarePlanViewer
- MonitoringDashboard
- Real-time status updates

---

## 📡 API ENDPOINTS

### Master Orchestration Endpoint

- `POST /api/orchestrate/full` - Run complete pipeline
- `GET /api/orchestrate/status/{patient_id}` - Get status
- `GET /api/orchestrate/state/{patient_id}` - Get full state
- `POST /api/orchestrate/event` - Process event
- `GET /api/orchestrate/states` - List all states
- `GET /api/orchestrate/health` - Health check

---

## 🎯 THE ULTIMATE MOAT

### What No Competitor Has:

| Capability | Us | Competitors |
|------------|-----|-------------|
| **End-to-end orchestration** | ✅ Upload → Care plan → Monitoring | ❌ Fragmented tools |
| **Agentic intelligence** | ✅ Autonomous agents, triggers, alerts | ❌ Manual workflows |
| **Validated biomarkers** | ✅ TMB r=0.933, MAPK RR=1.97 | ❌ Unvalidated or black-box |
| **Mechanism-based matching** | ✅ 7D vector, mechanism fit | ❌ Keyword matching only |
| **Toxicity-nutrition integration** | ✅ Drug MoA → Food recommendations | ❌ Generic advice |
| **Continuous monitoring** | ✅ Real-time triggers, early alerts | ❌ Point-in-time analysis |
| **Transparent reasoning** | ✅ Full provenance, auditable | ❌ Black-box AI |

---

## 📚 Full Reference

**Original Document:** `.cursor/MOAT/ULTIMATE_MOAT_ORCHESTRATION.mdc` (8,170 lines)

**Key Sections:**
- Complete workflow diagrams
- Detailed agent specifications
- Event trigger definitions
- API endpoint specifications
- Security & compliance measures
- Success metrics
- Implementation roadmap

---

**See Also:**
- [00_MISSION.mdc](00_MISSION.mdc) - Mission and vision
- [01_CURRENT_STATE.md](01_CURRENT_STATE.md) - What's built
- [03_DELIVERABLES_PLAN.md](03_DELIVERABLES_PLAN.md) - Deliverables breakdown
- [04_IMPLEMENTATION_ROADMAP.md](04_IMPLEMENTATION_ROADMAP.md) - Implementation plan


