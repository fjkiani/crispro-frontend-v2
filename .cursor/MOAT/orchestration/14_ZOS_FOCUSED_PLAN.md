# Zo's Focused Plan: Core Orchestrator Integration

**Date:** January 28, 2025  
**Purpose:** Zo's focused workload - core orchestrator integration only  
**Status:** ✅ **FOCUSED PLAN COMPLETE**

---

## 🎯 ZO'S MISSION

**Focus:** Core orchestrator agent integration and infrastructure  
**Delegated:** Frontend work, documentation, testing infrastructure, security implementation, monitoring setup  
**Time:** ~2-3 weeks (down from 6-8 weeks)

---

## ✅ ZO'S CORE DELIVERABLES (9 Total)

### **Week 1: Critical Blockers** (12-16 hours)

#### **Deliverable 1: Data Extraction Agent** ✅ ZO
**Time:** 4-6 hours  
**Priority:** 🔴 CRITICAL (BLOCKING)

**Tasks:**
- VCF parser implementation (PyVCF or cyvcf2)
- PDF parser (LLM-based extraction)
- MAF parser (tab-delimited)
- Clinical data extraction (stage, histology, biomarkers)
- Demographics extraction (age, sex, ECOG)
- Wire to orchestrator (`_run_data_extraction_agent()`)
- Unit tests

**Files:**
- `api/services/orchestrator/agents/data_extraction_agent.py` (NEW)
- `api/services/orchestrator/orchestrator.py` (modify)

**Acceptance:**
- ✅ Can parse VCF/PDF/MAF files
- ✅ Outputs structured PatientProfile
- ✅ Validates data quality
- ✅ End-to-end test passes

---

#### **Deliverable 2: Drug Efficacy Integration** ✅ ZO
**Time:** 8-10 hours  
**Priority:** 🔴 CRITICAL

**Tasks:**
- Wire S/P/E framework to orchestrator
- Import services directly (NOT HTTP calls)
- Connect to `_run_drug_efficacy_agent()`
- Test end-to-end

**Files:**
- `api/services/orchestrator/orchestrator.py` (modify `_run_drug_efficacy_agent()`)
- Import from: `api/services/efficacy_orchestrator/` (direct import)

**Acceptance:**
- ✅ S/P/E framework integrated
- ✅ Direct service imports (no HTTP)
- ✅ Drug ranking output in PatientState
- ✅ Performance: <2 seconds
- ✅ End-to-end test passes

---

### **Week 2: Core Intelligence** (1.5-2 days)

#### **Deliverable 3: Nutrition Integration** ✅ ZO
**Time:** 4-6 hours  
**Priority:** 🟡 HIGH

**Tasks:**
- Wire nutrition services to orchestrator
- Import services directly (NOT HTTP calls)
- Connect to `_run_nutrition_agent()`
- Test end-to-end

**Files:**
- `api/services/orchestrator/orchestrator.py` (modify `_run_nutrition_agent()`)
- Import from: `api/services/nutrition/` or `api/services/toxicity_pathway_mappings` (direct import)

**Acceptance:**
- ✅ Nutrition services integrated
- ✅ Direct service imports (no HTTP)
- ✅ Nutrition plan output in PatientState
- ✅ End-to-end test passes

---

#### **Deliverable 15: Data Validation & Quality** ✅ ZO
**Time:** 1 day  
**Priority:** 🟡 HIGH

**Tasks:**
- Input validation rules
- Data quality checks
- Coverage thresholds
- Mutation validation
- Clinical data validation
- Quality scoring
- Validation error reporting

**Files:**
- `api/services/orchestrator/agents/data_extraction_agent.py` (add validation)
- `api/services/orchestrator/validation.py` (NEW - optional)

**Acceptance:**
- ✅ Input validation rules implemented
- ✅ Data quality checks working
- ✅ Coverage thresholds enforced
- ✅ Quality scoring working
- ✅ Validation errors reported clearly

---

### **Week 3: Core Infrastructure** (2-4 days)

#### **Deliverable 10: Error Handling & Recovery** ✅ ZO
**Time:** 1-2 days  
**Priority:** 🟡 HIGH

**Tasks:**
- Agent failure recovery strategy
- Partial failure handling
- State recovery after crash
- Retry logic and backoff
- Circuit breaker patterns
- Graceful degradation

**Files:**
- `api/services/orchestrator/orchestrator.py` (add error handling)
- `api/services/orchestrator/error_recovery.py` (NEW)

**Acceptance:**
- ✅ Agent failures handled gracefully
- ✅ Partial failures don't break pipeline
- ✅ State recovery works after crash
- ✅ Retry logic prevents transient failures
- ✅ Circuit breakers prevent cascade failures

---

#### **Deliverable 13: State Persistence & Recovery** ✅ ZO
**Time:** 1-2 days  
**Priority:** 🟡 HIGH

**Tasks:**
- State persistence implementation (database)
- Recovery mechanism
- State versioning
- Cleanup/archival strategy

**Files:**
- `api/services/orchestrator/state_store.py` (enhance)
- `api/services/orchestrator/state_persistence.py` (NEW)

**Acceptance:**
- ✅ State persisted to database
- ✅ Recovery works after crash
- ✅ State versioning implemented
- ✅ Cleanup/archival working

---

### **Week 4: Automation** (4-6 hours)

#### **Deliverable 5: Trigger System** ✅ ZO
**Time:** 4-6 hours  
**Priority:** 🟡 HIGH

**Tasks:**
- Event detection logic
- Automated actions
- Escalation protocols
- Integration with monitoring agent

**Files:**
- `api/services/orchestrator/trigger_engine.py` (NEW or enhance existing)
- `api/services/orchestrator/orchestrator.py` (integrate triggers)

**Acceptance:**
- ✅ Event detection working
- ✅ Automated actions triggered
- ✅ Escalation protocols implemented
- ✅ Integration with monitoring agent working

---

### **Week 5-6: Agent Integration** (4-6 days)

#### **Deliverable 7: Access & Advocacy Agent Integration** ⚠️ ZO+PLUMBER
**Time:** 2-3 days (Zo's part)  
**Priority:** 🔴 CRITICAL* (not blocking)

**Zo's Tasks:**
- Wire agent to orchestrator
- Define orchestrator integration points
- Test orchestrator integration
- Provide agent interface documentation

**Plumber Tasks (Parallel):**
- Insurance packet generation logic
- Prior authorization support
- Financial assistance navigation
- Agent implementation details

**Files (Zo):**
- `api/services/orchestrator/orchestrator.py` (add `_run_access_advocacy_agent()`)
- `api/services/orchestrator/agents/access_advocacy_agent.py` (interface/skeleton)

**Acceptance (Zo):**
- ✅ Agent wired to orchestrator
- ✅ Integration points defined
- ✅ Orchestrator integration tested
- ✅ Interface documented for plumber

---

#### **Deliverable 8: Toxicity Risk Agent Integration** ⚠️ ZO+PLUMBER
**Time:** 2-3 days (Zo's part)  
**Priority:** 🔴 CRITICAL* (not blocking)

**Zo's Tasks:**
- Wire agent to orchestrator
- Define orchestrator integration points
- Test orchestrator integration
- Provide agent interface documentation

**Plumber Tasks (Parallel):**
- Germline-based toxicity prediction logic
- Pharmacogene flagging
- Pathway overlap analysis
- Mitigating food recommendations

**Files (Zo):**
- `api/services/orchestrator/orchestrator.py` (add `_run_toxicity_risk_agent()`)
- `api/services/orchestrator/agents/toxicity_risk_agent.py` (interface/skeleton)

**Acceptance (Zo):**
- ✅ Agent wired to orchestrator
- ✅ Integration points defined
- ✅ Orchestrator integration tested
- ✅ Interface documented for plumber

---

## 📊 ZO'S WORKLOAD SUMMARY

### **Time Breakdown:**
- **Week 1:** 12-16 hours (1.5-2 days)
- **Week 2:** 1.5-2 days
- **Week 3:** 2-4 days
- **Week 4:** 4-6 hours (0.5-1 day)
- **Week 5-6:** 4-6 days

**Total:** ~2-3 weeks of focused core work

### **Deliverables:**
- ✅ 7 Pure Core Deliverables (Zo only)
- ⚠️ 2 Shared Deliverables (Zo integration + Plumber implementation)

### **Focus Areas:**
1. Orchestrator agent integration
2. Service wiring patterns
3. Core infrastructure (error handling, state persistence)
4. Critical blocking functionality

---

## 🔧 PLUMBER HANDOFFS

### **What Zo Provides to Plumbers:**

#### **For Frontend (Deliverable 4, 12):**
- ✅ Orchestrator endpoints documented (`/api/orchestrate/full`, `/api/orchestrate/status/{patient_id}`)
- ✅ API contracts defined
- ✅ Test patient data provided
- ✅ Example responses provided

#### **For Agent Implementation (Deliverable 7, 8):**
- ✅ Orchestrator integration points defined
- ✅ Agent interface documented
- ✅ Test harness provided
- ✅ Example agent implementation provided

#### **For Testing (Deliverable 9):**
- ✅ Unit test examples provided
- ✅ Test data requirements documented
- ✅ Mock service requirements defined

#### **For Security (Deliverable 6):**
- ✅ Orchestrator authentication/authorization hooks
- ✅ Audit logging hooks
- ✅ Data encryption support

#### **For Monitoring (Deliverable 14):**
- ✅ Orchestrator metrics/logs exposed
- ✅ Metric definitions provided
- ✅ Logging format documented

---

## ✅ ZO's SUCCESS CRITERIA

### **Week 1 Complete:**
- ✅ Data Extraction Agent working
- ✅ Drug Efficacy integrated
- ✅ End-to-end test: Upload → Extract → Rank → Display

### **Week 2 Complete:**
- ✅ Nutrition integrated
- ✅ Data Validation complete
- ✅ All core agents can run end-to-end

### **Week 3 Complete:**
- ✅ Error Handling complete
- ✅ State Persistence complete
- ✅ Orchestrator infrastructure solid

### **Week 4 Complete:**
- ✅ Trigger System working
- ✅ Event automation functional

### **Week 5-6 Complete:**
- ✅ Access & Advocacy Agent integrated (orchestrator wiring)
- ✅ Toxicity Risk Agent integrated (orchestrator wiring)
- ✅ All core orchestrator work complete

---

## 🎯 ZO'S FOCUSED OUTCOMES

### **By End of Week 6, Zo Will Have:**
- ✅ All blocking agents integrated (Data Extraction, Drug Efficacy, Nutrition)
- ✅ Core infrastructure complete (Error Handling, State Persistence)
- ✅ Automation working (Trigger System)
- ✅ All agent integration points defined and tested
- ✅ Orchestrator ready for plumber work (frontend, security, monitoring, etc.)

### **Plumbers Can Then:**
- ✅ Build frontend on top of orchestrator
- ✅ Implement security on orchestrator hooks
- ✅ Set up monitoring on orchestrator metrics
- ✅ Implement agent logic for Access & Advocacy and Toxicity Risk
- ✅ Complete all non-core work in parallel

---

**See Also:**
- [13_TASK_DELEGATION.md](13_TASK_DELEGATION.md) - Complete delegation breakdown
- [03_DELIVERABLES_PLAN.md](03_DELIVERABLES_PLAN.md) - All 18 deliverables
- [04_IMPLEMENTATION_ROADMAP.md](04_IMPLEMENTATION_ROADMAP.md) - Week-by-week plan

