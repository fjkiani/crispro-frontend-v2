# Orchestration System: Deliverables Plan

**Date:** January 28, 2025  
**Status:** ✅ **PLAN COMPLETE** - Priorities defined, timeline established, tasks delegated  
**Source:** Consolidated from ORCHESTRATION_SCOPE_SYNTHESIS, ORCHESTRATOR_DASHBOARD_VERIFICATION, ULTIMATE_MOAT_ORCHESTRATION  
**Delegation:** See [13_TASK_DELEGATION.md](13_TASK_DELEGATION.md) for core vs plumber task breakdown

---

## 🎯 Executive Summary

**Current State:**
- ✅ Foundation: 100% Complete (Orchestrator, State Management, API Contracts)
- ⏳ Core Agents: 60% Complete (7/14 agents integrated)
- ⏳ Frontend: 20% Complete (OrchestratorDashboard verified, Universal Pages need integration)
- ⏳ Advanced Features: 20% Complete

**Critical Gaps:**
1. **Data Extraction Agent** - BLOCKING (nothing can run without it)
2. **Drug Efficacy Integration** - HIGH PRIORITY (core drug ranking)
3. **Nutrition Integration** - HIGH PRIORITY (MOAT feature)
4. **Frontend Integration** - HIGH PRIORITY (Universal Pages need orchestrator)

**Total Estimated Time:** 6-8 weeks of focused development (expanded from 4-6 weeks to include missing elements)

---

## 🔴 CRITICAL DELIVERABLES (Blocking)

### Deliverable 1: Data Extraction Agent
**Priority:** 🔴 CRITICAL  
**Complexity:** Medium (4-6 hours)  
**Status:** ⏳ SKELETON (20% complete)  
**Dependencies:** None

**What's Needed:**
- VCF parser (use PyVCF or cyvcf2)
- PDF parser (LLM-based extraction)
- MAF parser (tab-delimited)
- Clinical data extraction (stage, histology, biomarkers)
- Demographics extraction (age, sex, ECOG)

**Why Critical:**
- Nothing can run without extracted mutations
- All agents depend on PatientProfile from this agent
- Blocks entire pipeline execution

**Acceptance Criteria:**
- ✅ Can parse VCF files and extract mutations
- ✅ Can parse PDF reports (LLM-based extraction)
- ✅ Can parse MAF files
- ✅ Outputs structured PatientProfile object
- ✅ Validates data quality (missing fields, coverage thresholds)
- ✅ Flags ambiguities for human review

**Files to Create/Modify:**
- `api/services/orchestrator/agents/data_extraction_agent.py` (NEW)
- `api/services/orchestrator/orchestrator.py` (modify `_run_data_extraction_agent()`)

**Estimated Time:** 4-6 hours

---

### Deliverable 2: Drug Efficacy Integration
**Priority:** 🔴 CRITICAL  
**Complexity:** High (8-10 hours)  
**Status:** ⏳ SKELETON (80% - framework exists, needs wiring)  
**Dependencies:** 01 (Data Extraction), 02 (Biomarker), 03 (Resistance)

**What's Needed:**
- Wire S/P/E framework to orchestrator
- Import services directly (NOT HTTP calls)
- Connect to `_run_drug_efficacy_agent()`
- Test end-to-end

**Why Critical:**
- Core drug ranking needed for care plan
- S/P/E framework already proven and validated
- Required for trial matching and care plan generation

**Acceptance Criteria:**
- ✅ S/P/E framework integrated into orchestrator
- ✅ Direct service imports (no HTTP calls)
- ✅ Drug ranking output in PatientState
- ✅ End-to-end test passes
- ✅ Performance: <2 seconds for drug ranking

**Files to Modify:**
- `api/services/orchestrator/orchestrator.py` (modify `_run_drug_efficacy_agent()`)
- Import from: `api/services/efficacy_orchestrator/` (direct import)

**Estimated Time:** 8-10 hours

---

## 🟡 HIGH PRIORITY DELIVERABLES

### Deliverable 3: Nutrition Integration
**Priority:** 🟡 HIGH  
**Complexity:** Low-Medium (4-6 hours)  
**Status:** ⏳ SKELETON (70% - services exist, needs wiring)  
**Dependencies:** 01 (Data Extraction)

**What's Needed:**
- Wire toxicity-aware nutrition services to orchestrator
- Import services directly (NOT HTTP calls)
- Connect to `_run_nutrition_agent()`
- Test end-to-end

**Why High Priority:**
- MOAT feature, already has services built
- Toxicity pathway mappings exist
- Required for care plan generation

**Acceptance Criteria:**
- ✅ Nutrition services integrated into orchestrator
- ✅ Direct service imports (no HTTP calls)
- ✅ Nutrition plan output in PatientState
- ✅ End-to-end test passes

**Files to Modify:**
- `api/services/orchestrator/orchestrator.py` (modify `_run_nutrition_agent()`)
- Import from: `api/services/nutrition/` or `api/services/toxicity_pathway_mappings` (direct import)

**Estimated Time:** 4-6 hours

---

### Deliverable 4: Universal Pages Orchestrator Integration
**Priority:** 🟡 HIGH  
**Complexity:** High (9-13 days total across 4 phases)  
**Status:** ⏳ IN PROGRESS (30% - audit complete, deliverables defined)  
**Dependencies:** None for Phase 1 (can proceed in parallel), Phase 2 requires orchestrator agents (1.1, 1.2, 2.1)  
**Owner:** 🔧 **PLUMBER (Frontend Developer)** - Zo ensures orchestrator endpoints ready

**What's Needed:**
- **4.1 Phase 1:** Complete Missing Components (Resistance Playbook, SAE Features) - 2-3 days
- **4.2 Phase 2:** Orchestrator Integration - 3-4 days (requires orchestrator agents working)
- **4.3 Phase 3:** Testing & Validation - 2-3 days
- **4.4 Phase 4:** Enhancements - 2-3 days

**Note:** This deliverable covers UniversalCompleteCare, UniversalTrialIntelligence, UniversalDossierBrowser, UniversalDossierDetail. AyeshaCompleteCare and AyeshaTrialExplorer are handled separately in Deliverable 12 (Legacy Frontend Migration).

**Why High Priority:**
- User-facing interface, blocks user adoption
- Universal Pages currently use legacy `/api/complete_care/v2`
- Need orchestrator integration for real-time updates

**Acceptance Criteria:**
- ✅ UniversalCompleteCare uses `/api/orchestrate/full`
- ✅ File upload works (VCF/PDF/MAF)
- ✅ Status polling shows real-time progress
- ✅ All agent outputs displayed correctly
- ✅ Error handling for failed agents

**See:** `.cursor/MOAT/UNIVERSAL_PAGES_AUDIT_AND_DELIVERABLES.md` for complete plan

**Estimated Time:** 9-13 days (4 phases)

---

### Deliverable 5: Trigger System
**Priority:** 🟡 HIGH  
**Complexity:** Medium (4-6 hours)  
**Status:** ⬜ TODO (0% - not started)  
**Dependencies:** All agents (runs after pipeline complete)  
**Owner:** ✅ **ZO (CORE)** - Core orchestrator functionality

**What's Needed:**
- Event detection logic (resistance, TMB-H, PD, etc.)
- Automated actions (alerts, re-ranking, trial re-match)
- Escalation protocols
- Integration with monitoring agent

**Why High Priority:**
- Enables continuous monitoring
- Automated response to clinical events
- Critical for "Track forever" vision

**Acceptance Criteria:**
- ✅ Event detection working
- ✅ Automated actions triggered
- ✅ Escalation protocols implemented
- ✅ Integration with monitoring agent

**Estimated Time:** 4-6 hours

---

## 🟢 MEDIUM PRIORITY DELIVERABLES

### Deliverable 6: Security & Compliance
**Priority:** 🟢 MEDIUM  
**Complexity:** Medium (4-6 hours)  
**Status:** ⬜ TODO (0% - not started)  
**Dependencies:** All modules

**What's Needed:**
- Authentication/authorization
- Data encryption
- Audit logging
- HIPAA compliance measures

**Why Medium Priority:**
- Production readiness requirement
- Can be done in parallel with other work
- Not blocking core functionality

**Estimated Time:** 4-6 hours

---

### Deliverable 7: Access & Advocacy Agent (Module 15)
**Priority:** 🔴 CRITICAL (but not blocking)  
**Complexity:** High (1-2 weeks)  
**Status:** ⬜ TODO (0% - not started)  
**Dependencies:** 04 (Drug Efficacy), 07 (Care Plan)

**What's Needed:**
- Insurance packet generation
- Prior authorization support
- Financial assistance navigation
- Integration with OrchestratorDashboard

**Why Critical:**
- Required for complete care plan
- Frontend already has placeholder
- Not blocking other deliverables

**Estimated Time:** 1-2 weeks

---

### Deliverable 8: Toxicity Risk Agent (Module 16)
**Priority:** 🔴 CRITICAL (but not blocking)  
**Complexity:** Medium (4-6 hours)  
**Status:** ⬜ TODO (0% - not started)  
**Dependencies:** 01 (Data Extraction)

**What's Needed:**
- Germline-based toxicity prediction
- Pharmacogene flagging
- Pathway overlap analysis
- Mitigating food recommendations

**Why Critical:**
- Safety-critical feature
- Required for complete care plan
- Not blocking other deliverables

**Estimated Time:** 4-6 hours

---

## 📊 Priority Matrix (All 18 Deliverables)

| # | Deliverable | Priority | Time | Dependencies | Impact | Status |
|---|-------------|----------|------|--------------|--------|--------|
| **1** | Data Extraction | 🔴 CRITICAL | 4-6h | None | BLOCKING | ⏳ SKELETON |
| **2** | Drug Efficacy | 🔴 CRITICAL | 8-10h | 01 (Data Extraction), Biomarker Agent (02), Resistance Agent (03) | HIGH | ⏳ SKELETON |
| **3** | Nutrition | 🟡 HIGH | 4-6h | 01 | HIGH | ⏳ SKELETON |
| **4** | Universal Pages | 🟡 HIGH | 9-13d | None (4.1), 1.1,1.2,2.1 (4.2) | HIGH | ⏳ IN PROGRESS |
| **5** | Trigger System | 🟡 HIGH | 4-6h | All | MEDIUM | ⬜ TODO |
| **6** | Security & Compliance | 🟢 MEDIUM | 1-2d | All | HIGH | ⬜ TODO |
| **7** | Access & Advocacy | 🔴 CRITICAL* | 1-2w | 04, 07 | MEDIUM | ⬜ TODO |
| **8** | Toxicity Risk | 🔴 CRITICAL* | 4-6h | 01 | MEDIUM | ⬜ TODO |
| **9** | Testing Infrastructure | 🟡 HIGH | 2-3d | All | HIGH | ⬜ TODO |
| **10** | Error Handling & Recovery | 🟡 HIGH | 1-2d | All | HIGH | ⬜ TODO |
| **11** | Migration Strategy | 🟡 HIGH | 1d | All | HIGH | ⬜ TODO |
| **12** | Legacy Frontend Migration | 🟡 HIGH | 2-3d | 11, 4.2 | MEDIUM | ⬜ TODO |
| **13** | State Persistence & Recovery | 🟡 HIGH | 1-2d | 10 | HIGH | ⬜ TODO |
| **14** | Monitoring & Observability | 🟢 MEDIUM | 1-2d | All | MEDIUM | ⬜ TODO |
| **15** | Data Validation & Quality | 🟡 HIGH | 1d | 01 | HIGH | ⬜ TODO |
| **16** | API Versioning | 🟢 MEDIUM | 0.5d | 11 | MEDIUM | ⬜ TODO |
| **17** | Documentation Updates | 🟢 MEDIUM | 1-2d | All | MEDIUM | ⬜ TODO |
| **18** | Concurrency & Scalability | 🟢 MEDIUM | 1-2d | All | MEDIUM | ⬜ TODO |

**Legend:**
- 🔴 CRITICAL = Critical for system completeness
- 🔴 CRITICAL* = Critical but not blocking (can be done later)
- 🟡 HIGH = High priority, should be done soon
- 🟢 MEDIUM = Medium priority, can be done later

---

## 🚀 RECOMMENDED IMPLEMENTATION ORDER

### Week 1: Critical Blockers
1. **Deliverable 1: Data Extraction Agent** (4-6 hours)
   - Enables all other agents
   - Unblocks entire pipeline

2. **Deliverable 2: Drug Efficacy Integration** (8-10 hours)
   - Core drug ranking
   - Required for care plan

**Total Week 1:** 12-16 hours (1.5-2 days)

### Week 2: High Priority
3. **Deliverable 3: Nutrition Integration** (4-6 hours)
   - MOAT feature
   - Completes care plan

4. **Deliverable 4: Universal Pages Phase 1** (2-3 days)
   - Complete missing components
   - Unblocks frontend integration

**Total Week 2:** 3-4 days

### Week 3: Integration & Testing
5. **Deliverable 4.2: Universal Pages Phase 2** (3-4 days) - Orchestrator Integration
   - Orchestrator integration
   - File upload capability
   - Status polling

**Total Week 3:** 3-4 days

### Week 4: Automation & Polish
6. **Deliverable 5: Trigger System** (4-6 hours)
   - Event automation
   - Continuous monitoring

7. **Deliverable 4: Universal Pages Phases 3-4** (4-6 days)
   - Testing & Validation
   - Enhancements

**Total Week 4:** 5-7 days

### Week 5: Automation & Polish
9. **Deliverable 5: Trigger System** (4-6 hours)
10. **Deliverable 4.3: Universal Pages Phase 3** (2-3 days) - Testing & Validation
11. **Deliverable 4.4: Universal Pages Phase 4** (2-3 days) - Enhancements
12. **Deliverable 13: State Persistence & Recovery** (1-2 days)

### **PLUMBER WORKLOAD (Parallel with Zo)**

**Can Start Immediately:**
- Deliverable 4: Universal Pages (9-13d) - Frontend Developer
- Deliverable 9: Testing Infrastructure (2-3d) - QA Engineer
- Deliverable 11: Migration Strategy (1d) - Technical Lead
- Deliverable 17: Documentation (1-2d) - Technical Writer

**After Week 1 (Orchestrator Ready):**
- Deliverable 12: Legacy Frontend Migration (2-3d) - Frontend Developer
- Deliverable 14: Monitoring (1-2d) - DevOps Engineer
- Deliverable 18: Concurrency (1-2d) - Backend Engineer

**After Week 2 (Agents Working):**
- Deliverable 6: Security (1-2d) - Security Specialist
- Deliverable 16: API Versioning (0.5d) - API Architect

**After Week 5 (Agents Integrated):**
- Deliverable 7: Access & Advocacy (Implementation, 1-2w) - Agent Developer
- Deliverable 8: Toxicity Risk (Implementation, 4-6h) - Agent Developer

---

## 🎯 Success Criteria

### Phase 1 Complete When:
- ✅ Data Extraction Agent working (VCF/PDF/MAF parsers)
- ✅ Drug Efficacy integrated into orchestrator
- ✅ End-to-end test: Upload → Extract → Rank → Display

### Phase 2 Complete When:
- ✅ Nutrition integrated into orchestrator
- ✅ Universal Pages Phase 1 complete (missing components)
- ✅ All agents can run end-to-end

### Phase 3 Complete When:
- ✅ Universal Pages Phase 2 complete (orchestrator integration)
- ✅ File upload working
- ✅ Status polling working
- ✅ Real-time updates working

### Phase 4 Complete When:
- ✅ Trigger System working
- ✅ Universal Pages Phases 3-4 complete (testing, enhancements)
- ✅ Security & Compliance complete
- ✅ Production-ready system

---

## ⚠️ CRITICAL GAPS & MISSING ELEMENTS

### 1. Testing Strategy (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** HIGH - Need comprehensive testing approach

**What's Missing:**
- Unit test coverage targets (>80% for agents)
- Integration test scenarios
- End-to-end test workflows
- Performance test benchmarks
- Load testing strategy
- Error scenario testing
- User acceptance testing plan

**Add to Deliverables:**
- **Deliverable 9: Testing Infrastructure** (2-3 days)
  - Test framework setup
  - Test data fixtures
  - Mock services
  - CI/CD integration
  - Coverage reporting

### 2. Error Handling & Recovery (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** HIGH - Production reliability

**What's Missing:**
- Agent failure recovery strategy
- Partial failure handling (some agents succeed, others fail)
- State recovery after crash
- Retry logic and backoff
- Circuit breaker patterns
- Graceful degradation

**Add to Deliverables:**
- **Deliverable 10: Error Handling & Recovery** (1-2 days)  
  **Owner:** ✅ **ZO (CORE)** - Core orchestrator infrastructure
  - Agent failure recovery
  - State persistence and recovery
  - Retry logic
  - Circuit breakers
  - Graceful degradation

### 3. Performance Benchmarks (Not Specified)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** MEDIUM - User experience

**What's Missing:**
- Specific performance targets (not just "<2 seconds")
- Load testing requirements
- Concurrent user limits
- Scalability targets
- Resource usage limits

**Add to Deliverables:**
- **Performance Targets:**
  - Data Extraction: <5 seconds for VCF, <30 seconds for PDF
  - Agent execution: <2 seconds per agent (parallel)
  - Complete pipeline: <60 seconds end-to-end
  - API response time (P95): <2 seconds
  - Concurrent patients: 100+ simultaneous
  - State retrieval: <500ms

### 4. Migration Strategy (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** HIGH - Transition from legacy

**What's Missing:**
- How to migrate from `/api/complete_care/v2` to `/api/orchestrate/full`
- Backward compatibility strategy
- Feature flag approach
- Gradual rollout plan
- Rollback procedure

**Add to Deliverables:**
- **Deliverable 11: Migration Strategy** (1 day)  
  **Owner:** 🔧 **PLUMBER (Product/Technical Lead)** - Strategy and planning, Zo implements
  - Feature flags for orchestrator
  - Backward compatibility layer
  - Gradual rollout plan
  - Rollback procedure
  - Migration documentation

### 5. Other Frontend Pages (Not Identified)
**Status:** ⚠️ MISSING  
**Impact:** MEDIUM - Incomplete integration

**What's Missing:**
- AyeshaCompleteCare.jsx - Uses legacy endpoints?
- AyeshaTrialExplorer.jsx - Needs orchestrator integration?
- DoctorDashboard.jsx - Needs orchestrator integration?
- Other pages using legacy endpoints

**Add to Deliverables:**
- **Deliverable 12: Legacy Frontend Migration** (2-3 days)  
  **Owner:** 🔧 **PLUMBER (Frontend Developer)** - All frontend migration work
  - Audit all frontend pages for legacy endpoints
  - Migrate AyeshaCompleteCare to orchestrator
  - Migrate AyeshaTrialExplorer to orchestrator
  - Update DoctorDashboard if needed
  - Remove legacy endpoint dependencies

### 6. State Persistence & Recovery (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** HIGH - Data reliability

**What's Missing:**
- State persistence strategy (database, file system?)
- Recovery after orchestrator crash
- State versioning
- State cleanup/archival
- State migration strategy

**Add to Deliverables:**
- **Deliverable 13: State Persistence & Recovery** (1-2 days)  
  **Owner:** ✅ **ZO (CORE)** - Core orchestrator infrastructure
  - State persistence implementation
  - Recovery mechanism
  - State versioning
  - Cleanup/archival strategy

### 7. Monitoring & Observability (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** MEDIUM - Operational visibility

**What's Missing:**
- Logging strategy
- Metrics collection
- Alerting rules
- Dashboard for operations
- Performance monitoring
- Error tracking

**Add to Deliverables:**
- **Deliverable 14: Monitoring & Observability** (1-2 days)  
  **Owner:** 🔧 **PLUMBER (DevOps/SRE)** - Zo exposes metrics, Plumber sets up monitoring
  - Structured logging
  - Metrics collection (Prometheus/Grafana?)
  - Alerting rules
  - Operations dashboard
  - Error tracking (Sentry?)

### 8. Data Validation & Quality (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** HIGH - Data integrity

**What's Missing:**
- Input validation rules
- Data quality checks
- Coverage thresholds
- Mutation validation
- Clinical data validation
- Quality scoring

**Add to Deliverables:**
- **Deliverable 15: Data Validation & Quality** (1 day)  
  **Owner:** ✅ **ZO (CORE)** - Core data extraction quality
  - Input validation rules
  - Data quality checks
  - Coverage thresholds
  - Quality scoring
  - Validation error reporting

### 9. API Versioning (Not Addressed)
**Status:** ⚠️ MISSING  
**Impact:** MEDIUM - Future compatibility

**What's Missing:**
- API versioning strategy
- Version migration plan
- Deprecation policy
- Breaking change handling

**Add to Deliverables:**
- **Deliverable 16: API Versioning** (0.5 day)  
  **Owner:** 🔧 **PLUMBER (API Architect)** - Strategy and design, Zo implements
  - Versioning strategy
  - Migration plan
  - Deprecation policy

### 10. Documentation Updates (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** MEDIUM - Developer experience

**What's Missing:**
- API documentation updates
- Developer guide updates
- User documentation
- Migration guides
- Troubleshooting guides

**Add to Deliverables:**
- **Deliverable 17: Documentation Updates** (1-2 days)  
  **Owner:** 🔧 **PLUMBER (Technical Writer)** - All documentation writing, Zo reviews
  - API documentation
  - Developer guides
  - User documentation
  - Migration guides
  - Troubleshooting guides

### 11. Security & Privacy (Not Detailed)
**Status:** ⚠️ MISSING DETAIL  
**Impact:** HIGH - Compliance

**What's Missing:**
- PHI handling procedures
- HIPAA compliance measures
- Data encryption at rest
- Data encryption in transit
- Access controls
- Audit logging

**Add to Deliverables:**
- **Enhance Deliverable 6: Security & Compliance** (expand to 1-2 days)
  - PHI handling procedures
  - HIPAA compliance measures
  - Data encryption (at rest and in transit)
  - Access controls
  - Audit logging
  - Security audit

### 12. Concurrency & Scalability (Not Addressed)
**Status:** ⚠️ MISSING  
**Impact:** MEDIUM - System capacity

**What's Missing:**
- Concurrent patient handling
- Resource pooling
- Rate limiting
- Queue management
- Load balancing

**Add to Deliverables:**
- **Deliverable 18: Concurrency & Scalability** (1-2 days)  
  **Owner:** 🔧 **PLUMBER (Backend/Infrastructure Engineer)** - Zo ensures orchestrator supports it, Plumber implements infrastructure
  - Concurrent patient handling
  - Resource pooling
  - Rate limiting
  - Queue management
  - Load testing

---

## 📝 Implementation Notes

### Key Patterns (From Agent Implementation Guide)

**✅ DO: Import services directly**  
**❌ DON'T: Make HTTP calls from orchestrator**

**Example Pattern:**
```python
async def _run_drug_efficacy_agent(self, state: PatientState) -> Dict:
    """Run the drug efficacy agent."""
    execution = state.start_agent('drug_efficacy')
    
    try:
        # 1. Import service DIRECTLY (not HTTP)
        from ..efficacy_orchestrator import DrugEfficacyOrchestrator
        
        # 2. Build request from PatientState
        request = {
            'mutations': state.patient_profile.mutations,
            'biomarker_profile': state.biomarker_profile,
            'resistance_prediction': state.resistance_prediction
        }
        
        # 3. Call service method
        orchestrator = DrugEfficacyOrchestrator()
        result = await orchestrator.predict(request)
        
        # 4. Convert to dict if dataclass
        result_dict = result.to_dict() if hasattr(result, 'to_dict') else result
        
        # 5. execution.complete(result)
        execution.complete(result_dict)
        
        # 6. Return result
        return result_dict
    except Exception as e:
        execution.fail(str(e))
        raise
```

**Reference Examples:**
- `_run_biomarker_agent` (lines 450-483) - Inline logic pattern
- `_run_resistance_agent` (lines 485-565) - External service import pattern
- `_run_trial_matching_agent` (lines 625-691) - Agent class pattern

**📖 Full Guide:** `.cursor/MOAT/orchestration/agent-implementation-guide.mdc`

---

**See Also:**
- [01_CURRENT_STATE.md](01_CURRENT_STATE.md) - What's built
- [02_FRONTEND_STATUS.md](02_FRONTEND_STATUS.md) - Frontend status
- [04_IMPLEMENTATION_ROADMAP.md](04_IMPLEMENTATION_ROADMAP.md) - Detailed roadmap

