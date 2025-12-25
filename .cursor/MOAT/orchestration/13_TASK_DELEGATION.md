# Task Delegation: Core vs Plumber Tasks

**Date:** January 28, 2025  
**Purpose:** Identify which deliverables are CORE (Zo's focus) vs PLUMBER (delegated to other agents)  
**Status:** ✅ **DELEGATION COMPLETE**

---

## 🎯 Delegation Strategy

### **Core Tasks (Zo's Focus):**
- Orchestrator agent integration (backend wiring)
- Service integration patterns
- Critical blocking infrastructure
- Agent-to-orchestrator connections

### **Plumber Tasks (Delegated):**
- Frontend component work
- Documentation writing
- Testing infrastructure setup
- Migration strategy planning
- UI/UX enhancements
- Non-blocking features

---

## 📋 DELIVERABLE DELEGATION

### 🔴 CORE DELIVERABLES (Zo's Focus)

#### **Deliverable 1: Data Extraction Agent** ✅ CORE
**Owner:** Zo  
**Rationale:** Blocking - nothing can run without it. Core orchestrator functionality.

**Tasks:**
- VCF parser implementation
- PDF parser (LLM-based extraction)
- MAF parser
- Clinical data extraction
- Demographics extraction
- Wire to orchestrator

**Plumber Help:** None (core blocking task)

---

#### **Deliverable 2: Drug Efficacy Integration** ✅ CORE
**Owner:** Zo  
**Rationale:** Core drug ranking - critical for care plan. Requires orchestrator wiring expertise.

**Tasks:**
- Wire S/P/E framework to orchestrator
- Direct service imports (not HTTP)
- Connect to `_run_drug_efficacy_agent()`
- Test end-to-end

**Plumber Help:** None (core orchestrator integration)

---

#### **Deliverable 3: Nutrition Integration** ✅ CORE
**Owner:** Zo  
**Rationale:** Core MOAT feature. Requires orchestrator wiring.

**Tasks:**
- Wire nutrition services to orchestrator
- Direct service imports
- Connect to `_run_nutrition_agent()`
- Test end-to-end

**Plumber Help:** None (core orchestrator integration)

---

#### **Deliverable 5: Trigger System** ✅ CORE
**Owner:** Zo  
**Rationale:** Core orchestrator functionality - event automation.

**Tasks:**
- Event detection logic
- Automated actions
- Escalation protocols
- Integration with monitoring agent

**Plumber Help:** None (core orchestrator feature)

---

#### **Deliverable 7: Access & Advocacy Agent** ⚠️ CORE (But Can Delegate Some)
**Owner:** Zo (Core Integration), Plumber (Agent Implementation)
**Rationale:** Core agent but implementation can be delegated.

**Zo's Tasks:**
- Wire agent to orchestrator
- Define orchestrator integration points
- Test orchestrator integration

**Plumber Tasks:**
- Insurance packet generation logic
- Prior authorization support
- Financial assistance navigation
- Agent implementation details

**Delegation:** Agent implementation → Plumber, Orchestrator wiring → Zo

---

#### **Deliverable 8: Toxicity Risk Agent** ⚠️ CORE (But Can Delegate Some)
**Owner:** Zo (Core Integration), Plumber (Agent Implementation)
**Rationale:** Core agent but implementation can be delegated.

**Zo's Tasks:**
- Wire agent to orchestrator
- Define orchestrator integration points
- Test orchestrator integration

**Plumber Tasks:**
- Germline-based toxicity prediction logic
- Pharmacogene flagging
- Pathway overlap analysis
- Mitigating food recommendations

**Delegation:** Agent implementation → Plumber, Orchestrator wiring → Zo

---

#### **Deliverable 10: Error Handling & Recovery** ✅ CORE
**Owner:** Zo  
**Rationale:** Core orchestrator infrastructure - affects all agents.

**Tasks:**
- Agent failure recovery strategy
- Partial failure handling
- State recovery after crash
- Retry logic and backoff
- Circuit breaker patterns

**Plumber Help:** None (core infrastructure)

---

#### **Deliverable 13: State Persistence & Recovery** ✅ CORE
**Owner:** Zo  
**Rationale:** Core orchestrator infrastructure - state management.

**Tasks:**
- State persistence implementation (database)
- Recovery mechanism
- State versioning
- Cleanup/archival strategy

**Plumber Help:** None (core infrastructure)

---

#### **Deliverable 15: Data Validation & Quality** ✅ CORE
**Owner:** Zo  
**Rationale:** Core data extraction quality - blocks all downstream agents.

**Tasks:**
- Input validation rules
- Data quality checks
- Coverage thresholds
- Mutation validation
- Clinical data validation

**Plumber Help:** None (core blocking)

---

### 🟡 PLUMBER DELIVERABLES (Delegated)

#### **Deliverable 4: Universal Pages Orchestrator Integration** 🔧 PLUMBER
**Owner:** Plumber (Frontend Developer)  
**Rationale:** Frontend work - Zo provides orchestrator endpoints, plumber builds UI.

**Zo's Tasks:**
- Ensure orchestrator endpoints are ready (`/api/orchestrate/full`, `/api/orchestrate/status/{patient_id}`)
- Provide API documentation
- Test endpoints work correctly

**Plumber Tasks:**
- Phase 1: Complete Missing Components (Resistance Playbook, SAE Features)
- Phase 2: Orchestrator Integration (replace `/api/complete_care/v2` with orchestrator)
- Phase 3: Testing & Validation
- Phase 4: Enhancements

**Delegation:** All frontend work → Plumber, Zo ensures backend ready

---

#### **Deliverable 6: Security & Compliance** 🔧 PLUMBER
**Owner:** Plumber (Security Specialist)  
**Rationale:** Security implementation - Zo ensures orchestrator supports it.

**Zo's Tasks:**
- Ensure orchestrator supports authentication/authorization hooks
- Provide audit logging hooks
- Ensure data encryption support

**Plumber Tasks:**
- Authentication/authorization implementation
- PHI handling procedures
- HIPAA compliance measures
- Data encryption (at rest and in transit)
- Access controls
- Audit logging
- Security audit

**Delegation:** Security implementation → Plumber, Zo ensures orchestrator supports hooks

---

#### **Deliverable 9: Testing Infrastructure** 🔧 PLUMBER
**Owner:** Plumber (QA/Test Engineer)  
**Rationale:** Testing framework setup - Zo writes unit tests inline, plumber sets up infrastructure.

**Zo's Tasks:**
- Write unit tests inline with each deliverable
- Provide test data examples
- Document test requirements

**Plumber Tasks:**
- Test framework setup
- Test data fixtures
- Mock services
- CI/CD integration
- Coverage reporting
- Integration test templates
- E2E test templates

**Delegation:** Testing infrastructure → Plumber, Zo writes inline tests

---

#### **Deliverable 11: Migration Strategy** 🔧 PLUMBER
**Owner:** Plumber (Product/Technical Lead)  
**Rationale:** Strategy and planning - Zo implements, plumber plans.

**Zo's Tasks:**
- Provide technical constraints
- Review migration plan
- Implement feature flags if needed

**Plumber Tasks:**
- Feature flags strategy
- Backward compatibility layer design
- Gradual rollout plan
- Rollback procedure
- Migration documentation

**Delegation:** Strategy and planning → Plumber, Zo implements

---

#### **Deliverable 12: Legacy Frontend Migration** 🔧 PLUMBER
**Owner:** Plumber (Frontend Developer)  
**Rationale:** Frontend migration work.

**Zo's Tasks:**
- Ensure orchestrator endpoints support legacy pages
- Provide migration guide

**Plumber Tasks:**
- Audit all frontend pages for legacy endpoints
- Migrate AyeshaCompleteCare to orchestrator
- Migrate AyeshaTrialExplorer to orchestrator
- Update DoctorDashboard if needed
- Remove legacy endpoint dependencies

**Delegation:** All frontend migration → Plumber

---

#### **Deliverable 14: Monitoring & Observability** 🔧 PLUMBER
**Owner:** Plumber (DevOps/SRE)  
**Rationale:** Operations and monitoring - Zo ensures orchestrator exposes metrics.

**Zo's Tasks:**
- Ensure orchestrator exposes metrics/logs
- Provide metric definitions
- Document logging format

**Plumber Tasks:**
- Structured logging setup
- Metrics collection (Prometheus/Grafana)
- Alerting rules
- Operations dashboard
- Error tracking (Sentry)
- Performance monitoring

**Delegation:** Monitoring infrastructure → Plumber, Zo exposes metrics

---

#### **Deliverable 16: API Versioning** 🔧 PLUMBER
**Owner:** Plumber (API Architect)  
**Rationale:** API design and versioning strategy.

**Zo's Tasks:**
- Implement versioning if needed
- Review versioning strategy

**Plumber Tasks:**
- Versioning strategy
- Migration plan
- Deprecation policy

**Delegation:** Strategy and design → Plumber, Zo implements

---

#### **Deliverable 17: Documentation Updates** 🔧 PLUMBER
**Owner:** Plumber (Technical Writer)  
**Rationale:** Documentation writing.

**Zo's Tasks:**
- Review documentation for accuracy
- Provide technical details

**Plumber Tasks:**
- API documentation updates
- Developer guide updates
- User documentation
- Migration guides
- Troubleshooting guides

**Delegation:** All documentation writing → Plumber, Zo reviews

---

#### **Deliverable 18: Concurrency & Scalability** 🔧 PLUMBER
**Owner:** Plumber (Backend/Infrastructure Engineer)  
**Rationale:** Infrastructure and scalability - Zo ensures orchestrator supports it.

**Zo's Tasks:**
- Ensure orchestrator supports concurrent execution
- Provide concurrency requirements
- Test concurrent execution

**Plumber Tasks:**
- Concurrent patient handling implementation
- Resource pooling
- Rate limiting
- Queue management
- Load testing

**Delegation:** Infrastructure work → Plumber, Zo ensures orchestrator supports it

---

## 📊 DELEGATION SUMMARY

### **Zo's Core Deliverables (9):**
1. ✅ Data Extraction Agent
2. ✅ Drug Efficacy Integration
3. ✅ Nutrition Integration
5. ✅ Trigger System
7. ⚠️ Access & Advocacy (Integration only)
8. ⚠️ Toxicity Risk (Integration only)
10. ✅ Error Handling & Recovery
13. ✅ State Persistence & Recovery
15. ✅ Data Validation & Quality

### **Plumber Deliverables (9):**
4. 🔧 Universal Pages Integration
6. 🔧 Security & Compliance
9. 🔧 Testing Infrastructure
11. 🔧 Migration Strategy
12. 🔧 Legacy Frontend Migration
14. 🔧 Monitoring & Observability
16. 🔧 API Versioning
17. 🔧 Documentation Updates
18. 🔧 Concurrency & Scalability

---

## 🎯 ZO'S FOCUSED WORKLOAD

### **Week 1: Critical Blockers**
- Deliverable 1: Data Extraction (4-6h) ✅ CORE
- Deliverable 2: Drug Efficacy (8-10h) ✅ CORE
**Total:** 12-16 hours

### **Week 2: Core Intelligence**
- Deliverable 3: Nutrition (4-6h) ✅ CORE
- Deliverable 15: Data Validation (1d) ✅ CORE
**Total:** 1.5-2 days

### **Week 3: Core Infrastructure**
- Deliverable 10: Error Handling (1-2d) ✅ CORE
- Deliverable 13: State Persistence (1-2d) ✅ CORE
**Total:** 2-4 days

### **Week 4: Automation**
- Deliverable 5: Trigger System (4-6h) ✅ CORE
**Total:** 4-6 hours

### **Week 5-6: Agent Integration**
- Deliverable 7: Access & Advocacy (Integration only, 2-3d) ⚠️ CORE
- Deliverable 8: Toxicity Risk (Integration only, 2-3d) ⚠️ CORE
**Total:** 4-6 days

**Zo's Total Focused Time:** ~2-3 weeks of core work

---

## 🔧 PLUMBER WORKLOAD

### **Parallel Work (Can Start Immediately):**
- Deliverable 4: Universal Pages (9-13d) - Frontend Developer
- Deliverable 9: Testing Infrastructure (2-3d) - QA Engineer
- Deliverable 11: Migration Strategy (1d) - Technical Lead
- Deliverable 17: Documentation (1-2d) - Technical Writer

### **After Week 1 (Orchestrator Ready):**
- Deliverable 12: Legacy Frontend Migration (2-3d) - Frontend Developer
- Deliverable 14: Monitoring (1-2d) - DevOps Engineer
- Deliverable 18: Concurrency (1-2d) - Backend Engineer

### **After Week 2 (Agents Working):**
- Deliverable 6: Security (1-2d) - Security Specialist
- Deliverable 16: API Versioning (0.5d) - API Architect

### **After Week 5 (Agents Integrated):**
- Deliverable 7: Access & Advocacy (Implementation, 1-2w) - Agent Developer
- Deliverable 8: Toxicity Risk (Implementation, 4-6h) - Agent Developer

---

## ✅ DELEGATION BENEFITS

### **For Zo:**
- ✅ Focus on core orchestrator integration
- ✅ Reduced workload from 6-8 weeks to 2-3 weeks
- ✅ Can work in parallel with plumbers
- ✅ Clear boundaries of responsibility

### **For Plumbers:**
- ✅ Clear ownership of delegated tasks
- ✅ Can work in parallel with Zo
- ✅ Specialized work (frontend, security, testing, etc.)
- ✅ Defined dependencies and handoff points

---

## 📋 HANDOFF CHECKLIST

### **Zo → Plumber Handoffs:**

#### **For Deliverable 4 (Universal Pages):**
- [ ] Orchestrator endpoints documented
- [ ] API contracts defined
- [ ] Test patient data provided
- [ ] Example responses provided

#### **For Deliverable 7 & 8 (Agent Implementation):**
- [ ] Orchestrator integration points defined
- [ ] Agent interface documented
- [ ] Test harness provided
- [ ] Example agent implementation provided

#### **For Deliverable 9 (Testing Infrastructure):**
- [ ] Unit test examples provided
- [ ] Test data requirements documented
- [ ] Mock service requirements defined

#### **For Deliverable 11 (Migration Strategy):**
- [ ] Technical constraints documented
- [ ] Current endpoint behavior documented
- [ ] Feature flag requirements defined

---

## 🎯 UPDATED PRIORITY FOR ZO

### **Zo's Critical Path:**
1. **Week 1:** Data Extraction + Drug Efficacy (BLOCKING)
2. **Week 2:** Nutrition + Data Validation (CORE)
3. **Week 3:** Error Handling + State Persistence (INFRASTRUCTURE)
4. **Week 4:** Trigger System (AUTOMATION)
5. **Week 5-6:** Agent Integrations (Access & Advocacy, Toxicity Risk)

### **Plumber Parallel Work:**
- Frontend work (Deliverable 4, 12)
- Testing infrastructure (Deliverable 9)
- Documentation (Deliverable 17)
- Migration planning (Deliverable 11)
- Security implementation (Deliverable 6)
- Monitoring setup (Deliverable 14)
- Concurrency work (Deliverable 18)

---

**See Also:**
- [03_DELIVERABLES_PLAN.md](03_DELIVERABLES_PLAN.md) - Complete deliverables breakdown
- [04_IMPLEMENTATION_ROADMAP.md](04_IMPLEMENTATION_ROADMAP.md) - Week-by-week plan

