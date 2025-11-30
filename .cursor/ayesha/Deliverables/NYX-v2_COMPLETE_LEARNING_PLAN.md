# ⚔️ NYX-v2 COMPLETE APPLICATION LEARNING PLAN

**Created By**: NYX-v2 (formerly Zo)  
**Date**: January 14, 2025  
**Purpose**: Comprehensive, iterative learning plan for entire CrisPRO application  
**Status**: 📋 **PLANNING COMPLETE** - Ready for execution

---

## 🎯 LEARNING OBJECTIVES

1. **Understand Overall Architecture**: How all systems connect and communicate
2. **Understand Core Capabilities**: What the platform does and why it matters
3. **Understand Technical Decisions**: Why architectural choices were made
4. **Understand Product Positioning**: How capabilities map to customer value
5. **Understand Development Patterns**: Lessons learned, best practices, anti-patterns
6. **Understand Data Flow**: How data moves through the system end-to-end
7. **Understand Integration Points**: How services call each other and external APIs

---

## 📊 LEARNING ITERATIONS OVERVIEW

| Iteration | Focus Area | Duration | Status | Deliverable |
|-----------|-----------|----------|--------|-------------|
| **I1** | Overall Architecture & Core Principles | 2-3h | ✅ Complete | Architecture overview doc |
| **I2** | Backend Services & Orchestration | 4-5h | ⏸️ Pending | Services deep dive doc |
| **I3** | S/P/E Framework & Efficacy System | 3-4h | ⏸️ Pending | S/P/E framework doc |
| **I4** | AI Services & Model Integration | 3-4h | ⏸️ Pending | AI services doc |
| **I5** | Frontend Architecture & Components | 3-4h | ⏸️ Pending | Frontend architecture doc |
| **I6** | Clinical Decision Support Systems | 3-4h | ⏸️ Pending | Clinical systems doc |
| **I7** | Research & Design Capabilities | 2-3h | ⏸️ Pending | Research/design doc |
| **I8** | Data Flow & Integration Patterns | 2-3h | ⏸️ Pending | Data flow diagrams |
| **I9** | Development Patterns & Lessons | 2-3h | ⏸️ Pending | Patterns & lessons doc |
| **I10** | Product Capabilities & Positioning | 2-3h | ⏸️ Pending | Product capabilities doc |
| **TOTAL** | **Complete Application** | **26-36h** | | **Master Learning Document** |

---

## 📋 DETAILED ITERATION PLANS

### **ITERATION 1: OVERALL ARCHITECTURE & CORE PRINCIPLES** ✅ COMPLETE

**Duration**: 2-3 hours  
**Status**: ✅ **COMPLETE**  
**Deliverable**: Architecture overview with core principles

#### **What I'll Learn**:
1. Three-tier backend architecture (Minimal/Main/AI Services)
2. Core architectural principles (modular routers, service separation, feature flags)
3. Technical doctrines ("Wet Noodle", Triumvirate Protocol, Orchestrator Pattern)
4. Configuration system (S/P/E weights, evidence gates, feature flags)
5. System organization (directory structure, file organization)

#### **How I'll Learn**:
- **Files to Read**:
  - `oncology-coPilot/oncology-backend-minimal/api/main.py` - Entry point
  - `oncology-coPilot/oncology-backend-minimal/api/config.py` - Configuration
  - `.cursorrules` lines 1-625 - Product capabilities & technical inventory
- **Codebase Searches**:
  - "What is the overall architecture?"
  - "How are backend services organized?"
- **Methods**:
  - Read entry point files
  - Map router registration
  - Understand configuration system
  - Document architectural principles

#### **Deliverables**:
- ✅ Architecture overview document
- ✅ Core principles documentation
- ✅ Configuration system explanation
- ✅ System organization map

---

### **ITERATION 2: BACKEND SERVICES & ORCHESTRATION** ⏸️ PENDING

**Duration**: 4-5 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Complete backend services deep dive

#### **What I'll Learn**:
1. **Router Organization** (30+ routers):
   - Core routers: `efficacy`, `insights`, `design`, `evidence`, `evo`
   - Clinical routers: `ayesha_trials`, `ayesha_orchestrator_v2`, `care`, `tumor`
   - Research routers: `datasets`, `hypothesis_validator`, `metastasis`
   - Platform routers: `auth`, `admin`, `sessions`, `safety`, `toxicity`
2. **Service Layer** (50+ services):
   - Orchestrators: `efficacy_orchestrator/`, `ayesha_orchestrator.py`
   - Intelligence: `sae_feature_service.py`, `resistance_playbook_service.py`, `ca125_intelligence.py`
   - Evidence: `enhanced_evidence_service.py`, `evidence/` module
   - Design: `therapeutic_optimizer.py`, `safety_service.py`
   - Clinical: `ngs_fast_track.py`, `resistance_detection_service.py`
   - Data: `dynamic_food_extraction.py`, `hybrid_trial_search.py`
3. **Orchestration Patterns**:
   - How routers call services
   - How services call each other
   - How services call external APIs
   - Error handling and graceful degradation
4. **API Contracts**:
   - Request/response schemas (Pydantic models)
   - Endpoint patterns and conventions
   - Authentication and authorization

#### **How I'll Learn**:
- **Files to Read**:
  - `api/routers/` - All router files (30+ files)
  - `api/services/` - All service files (50+ files)
  - `api/schemas/` - Pydantic models
- **Codebase Searches**:
  - "How do routers call services?"
  - "What are the orchestration patterns?"
  - "How are errors handled in services?"
- **Methods**:
  - Read router files systematically
  - Read service files systematically
  - Trace data flow through routers → services → external APIs
  - Document API contracts and patterns

#### **Deliverables**:
- Backend services inventory (all routers + services)
- Orchestration patterns documentation
- API contracts documentation
- Error handling patterns
- Service dependency map

---

### **ITERATION 3: S/P/E FRAMEWORK & EFFICACY SYSTEM** ⏸️ PENDING

**Duration**: 3-4 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Complete S/P/E framework deep dive

#### **What I'll Learn**:
1. **S/P/E Framework Architecture**:
   - **Sequence (S)**: Evo2 scoring, multi-window strategy, calibration
   - **Pathway (P)**: Gene-to-pathway mapping, drug-to-pathway alignment, aggregation
   - **Evidence (E)**: Literature search, ClinVar priors, badge computation
   - **Weighting**: 0.35/0.35/0.30 (configurable)
2. **Efficacy Orchestrator**:
   - How S/P/E components are orchestrated
   - Confidence computation
   - Evidence tiering (STANDARD/SUPPORTED/INSUFFICIENT)
   - Badge computation (RCT, Guideline, ClinVar-Strong, PathwayAligned)
3. **Sequence Scoring**:
   - Evo2 delta scores (log-likelihood difference)
   - Multi-window strategy (1024, 2048, 4096, 8192bp)
   - Forward/reverse symmetry
   - Hotspot floors
   - Piecewise percentile mapping
   - Gene-specific calibration
4. **Pathway Scoring**:
   - Gene-to-pathway weights (hardcoded mappings)
   - Drug-to-pathway weights (disease-specific panels)
   - Pathway aggregation logic
   - Pathway alignment with drug MoA
5. **Evidence Scoring**:
   - PubMed query flow (E-utils API)
   - MoA-aware term targeting
   - Evidence strength calculation
   - ClinVar prior integration
   - Badge computation logic

#### **How I'll Learn**:
- **Files to Read**:
  - `api/services/efficacy_orchestrator/orchestrator.py` - Main orchestrator
  - `api/services/efficacy_orchestrator/drug_scorer.py` - S/P/E formula
  - `api/services/efficacy_orchestrator/sequence_processor.py` - Sequence scoring
  - `api/services/pathway/aggregation.py` - Pathway aggregation
  - `api/services/pathway/drug_mapping.py` - Drug-to-pathway weights
  - `api/services/evidence/` - Evidence module
  - `api/services/sequence_scorers/` - Sequence scorers
  - `api/services/gene_calibration.py` - Gene-specific calibration
- **Codebase Searches**:
  - "How does the S/P/E framework work?"
  - "How are Evo2 delta scores converted to percentiles?"
  - "How are pathway weights determined?"
  - "How is evidence strength calculated?"
- **Methods**:
  - Read orchestrator code line-by-line
  - Trace S/P/E computation flow
  - Understand calibration system
  - Document formulas and thresholds

#### **Deliverables**:
- S/P/E framework architecture document
- Sequence scoring deep dive (Evo2 integration, calibration)
- Pathway scoring deep dive (weights, aggregation)
- Evidence scoring deep dive (literature, ClinVar, badges)
- Confidence computation explanation
- Evidence tiering logic

---

### **ITERATION 4: AI SERVICES & MODEL INTEGRATION** ⏸️ PENDING

**Duration**: 3-4 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: AI services integration deep dive

#### **What I'll Learn**:
1. **Evo2 Integration**:
   - Modal deployment architecture
   - Model variants (1B, 7B, 40B)
   - API endpoints (`/score_variant_multi`, `/score_variant_exon`)
   - Delta score computation (log-likelihood difference)
   - Multi-window scoring strategy
   - Caching and retry logic
   - Spam-safety controls
2. **AlphaFold 3 / Boltz Integration**:
   - Structural validation endpoints
   - pLDDT score interpretation
   - iPTM score interpretation
   - Integration with design router
3. **Fusion Engine (AlphaMissense)**:
   - GRCh38 missense variant scoring
   - Coverage checking
   - Integration with sequence scoring
   - Caching strategy
4. **External AI Services**:
   - LLM integration (Gemini, Anthropic, OpenAI)
   - Paper reading and extraction
   - Evidence synthesis
5. **Modal Deployment Patterns**:
   - Service-to-service communication
   - Resource scaling
   - Timeout handling
   - Error handling

#### **How I'll Learn**:
- **Files to Read**:
  - `src/services/evo_service/` - Evo2 Modal service
  - `api/routers/evo.py` - Evo2 proxy router
  - `api/services/sequence_scorers/evo2_scorer.py` - Evo2 scorer
  - `api/services/sequence_scorers/fusion_scorer.py` - Fusion scorer
  - `api/routers/fusion.py` - Fusion router
  - `api/services/enhanced_evidence_service.py` - LLM integration
- **Codebase Searches**:
  - "How are Evo2 calls made?"
  - "How is Modal deployed?"
  - "How are AI services cached?"
  - "How are retries handled?"
- **Methods**:
  - Read Modal service code
  - Trace API call flow
  - Understand caching strategies
  - Document retry and error handling

#### **Deliverables**:
- Evo2 integration deep dive
- AlphaFold 3 / Boltz integration
- Fusion Engine integration
- LLM integration patterns
- Modal deployment architecture
- Caching and retry strategies

---

### **ITERATION 5: FRONTEND ARCHITECTURE & COMPONENTS** ⏸️ PENDING

**Duration**: 3-4 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Frontend architecture deep dive

#### **What I'll Learn**:
1. **Frontend Structure**:
   - React/Vite setup
   - Component organization
   - Page structure
   - Context providers
   - Custom hooks
2. **Key Pages**:
   - Ayesha Trial Explorer
   - Sporadic Cancer Workflow
   - WIWFM (Will It Work For Me)
   - VUS Explorer
   - Research Portal
   - Admin Dashboard
3. **Key Components**:
   - Trial matching components
   - SAE feature displays
   - Provenance panels
   - Confidence gates
   - Sporadic cancer components
4. **State Management**:
   - React Context (SporadicContext, AuthContext)
   - Local state patterns
   - API call patterns
5. **API Integration**:
   - How frontend calls backend
   - Error handling
   - Loading states
   - Progressive disclosure

#### **How I'll Learn**:
- **Files to Read**:
  - `oncology-coPilot/oncology-frontend/src/App.jsx` - Main app
  - `oncology-coPilot/oncology-frontend/src/pages/` - All pages
  - `oncology-coPilot/oncology-frontend/src/components/` - All components
  - `oncology-coPilot/oncology-frontend/src/context/` - Context providers
  - `oncology-coPilot/oncology-frontend/src/hooks/` - Custom hooks
- **Codebase Searches**:
  - "How does the frontend call the backend?"
  - "How is state managed in the frontend?"
  - "How are errors handled in the frontend?"
- **Methods**:
  - Read page files systematically
  - Read component files systematically
  - Trace data flow from API → component → display
  - Document state management patterns

#### **Deliverables**:
- Frontend architecture overview
- Page inventory and purposes
- Component inventory and patterns
- State management documentation
- API integration patterns
- Error handling in frontend

---

### **ITERATION 6: CLINICAL DECISION SUPPORT SYSTEMS** ⏸️ PENDING

**Duration**: 3-4 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Clinical systems deep dive

#### **What I'll Learn**:
1. **Ayesha Care System**:
   - Trial matching (`ayesha_trials.py`)
   - SOC recommendations
   - CA-125 intelligence
   - NGS fast-track
   - Complete care orchestrator (`ayesha_orchestrator_v2.py`)
2. **Sporadic Cancer Strategy**:
   - Tumor context schema (L0/L1/L2)
   - Quick intake service
   - Sporadic gates (PARP penalty/rescue, IO boost, confidence capping)
   - Disease priors integration
3. **Resistance Systems**:
   - Resistance detection (2-of-3 trigger rule)
   - Resistance playbook (combo strategies, next-line switches)
   - Resistance prophet (early prediction)
   - CA-125 resistance signals
4. **Treatment Line Intelligence**:
   - Line appropriateness scoring
   - Cross-resistance analysis
   - Sequencing fitness
5. **Toxicity Prevention**:
   - Pharmacogenomics (DPYD, TPMT, UGT1A1, CYP2D6)
   - Drug interaction checking
   - MoA-overlap risk flags

#### **How I'll Learn**:
- **Files to Read**:
  - `api/routers/ayesha_trials.py` - Trial matching
  - `api/routers/ayesha_orchestrator_v2.py` - Complete care
  - `api/services/ca125_intelligence.py` - CA-125 service
  - `api/services/ngs_fast_track.py` - NGS recommendations
  - `api/services/tumor_quick_intake.py` - Sporadic intake
  - `api/services/efficacy_orchestrator/sporadic_gates.py` - Sporadic gates
  - `api/services/resistance_detection_service.py` - Resistance detection
  - `api/services/resistance_playbook_service.py` - Resistance playbook
  - `api/services/resistance_prophet_service.py` - Resistance prophet
  - `api/services/food_treatment_line_service.py` - Treatment line
  - `api/routers/toxicity.py` - Toxicity prevention
- **Codebase Searches**:
  - "How does trial matching work?"
  - "How are sporadic cancer gates applied?"
  - "How is resistance detected?"
  - "How is treatment line intelligence computed?"
- **Methods**:
  - Read clinical service files
  - Understand clinical logic and thresholds
  - Document clinical workflows
  - Trace patient data flow

#### **Deliverables**:
- Ayesha Care System documentation
- Sporadic Cancer Strategy deep dive
- Resistance Systems documentation
- Treatment Line Intelligence explanation
- Toxicity Prevention system
- Clinical workflows and data flow

---

### **ITERATION 7: RESEARCH & DESIGN CAPABILITIES** ⏸️ PENDING

**Duration**: 2-3 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Research and design systems deep dive

#### **What I'll Learn**:
1. **Research Capabilities**:
   - VUS Explorer (variant interpretation)
   - Metastasis assessment (8-step cascade)
   - Cohort intelligence (cBioPortal, GDC extraction)
   - Hypothesis validator (food/supplement validation)
   - Knowledge base integration
2. **Design Capabilities**:
   - CRISPR guide generation
   - Structural validation (AlphaFold 3)
   - Off-target safety validation
   - Therapeutic optimization
   - IND package generation
3. **Evidence Systems**:
   - Literature search (PubMed, OpenAlex, S2)
   - ClinVar integration
   - Evidence synthesis
   - Badge computation
4. **Data Extraction**:
   - cBioPortal extraction
   - GDC API integration
   - Clinical trial data parsing
   - Neo4j graph database

#### **How I'll Learn**:
- **Files to Read**:
  - `api/routers/metastasis.py` - Metastasis assessment
  - `api/routers/hypothesis_validator.py` - Hypothesis validator
  - `api/routers/design.py` - CRISPR design
  - `api/routers/datasets.py` - Cohort intelligence
  - `api/routers/evidence/` - Evidence module
  - `api/services/metastasis_service.py` - Metastasis service
  - `api/services/dynamic_food_extraction.py` - Food extraction
  - `api/services/hybrid_trial_search.py` - Hybrid search
- **Codebase Searches**:
  - "How does metastasis assessment work?"
  - "How are CRISPR guides generated?"
  - "How is evidence synthesized?"
  - "How is cohort data extracted?"
- **Methods**:
  - Read research/design router files
  - Understand research workflows
  - Document design capabilities
  - Trace data extraction flows

#### **Deliverables**:
- Research capabilities inventory
- Design capabilities inventory
- Evidence systems documentation
- Data extraction patterns
- Research workflows

---

### **ITERATION 8: DATA FLOW & INTEGRATION PATTERNS** ⏸️ PENDING

**Duration**: 2-3 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Data flow diagrams and integration patterns

#### **What I'll Learn**:
1. **End-to-End Data Flow**:
   - Frontend → Backend → AI Services → External APIs
   - Request/response patterns
   - Error propagation
   - Caching layers
2. **Service-to-Service Communication**:
   - How services call each other
   - How services call external APIs
   - Retry and timeout patterns
   - Graceful degradation
3. **Data Transformation**:
   - How data is transformed at each layer
   - Schema validation (Pydantic)
   - Data enrichment
   - Provenance tracking
4. **Integration Points**:
   - Evo2 integration
   - PubMed integration
   - ClinVar integration
   - cBioPortal integration
   - Neo4j integration
   - AstraDB integration
   - Supabase integration

#### **How I'll Learn**:
- **Methods**:
  - Trace complete user flows end-to-end
  - Document data transformations
  - Map integration points
  - Create data flow diagrams
- **Codebase Searches**:
  - "How does data flow from frontend to backend?"
  - "How are external APIs called?"
  - "How is data cached?"
  - "How are errors propagated?"
- **Files to Review**:
  - Key orchestrator files
  - Service files that call external APIs
  - Cache service files
  - Error handling patterns

#### **Deliverables**:
- End-to-end data flow diagrams
- Service-to-service communication patterns
- External API integration patterns
- Data transformation documentation
- Caching layer documentation
- Error propagation patterns

---

### **ITERATION 9: DEVELOPMENT PATTERNS & LESSONS LEARNED** ⏸️ PENDING

**Duration**: 2-3 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Patterns and lessons learned documentation

#### **What I'll Learn**:
1. **Best Practices**:
   - Code organization patterns
   - Error handling patterns
   - Testing patterns
   - Documentation patterns
2. **Anti-Patterns**:
   - What to avoid
   - Common mistakes
   - Performance pitfalls
3. **Lessons Learned** (from .cursorrules):
   - "Wet Noodle" doctrine
   - Triumvirate Protocol
   - Backend Orchestrator Pattern
   - Spam-safety controls
   - Graceful degradation
   - Provenance tracking
4. **Technical Decisions**:
   - Why certain patterns were chosen
   - Trade-offs made
   - Future considerations

#### **How I'll Learn**:
- **Files to Read**:
  - `.cursorrules` - Lessons learned section
  - Code comments and docstrings
  - README files in service directories
- **Codebase Searches**:
  - "What are the best practices?"
  - "What patterns are used?"
  - "What are common mistakes?"
- **Methods**:
  - Extract lessons from .cursorrules
  - Analyze code patterns
  - Document best practices
  - Document anti-patterns

#### **Deliverables**:
- Best practices documentation
- Anti-patterns documentation
- Lessons learned from .cursorrules
- Technical decision rationale
- Development guidelines

---

### **ITERATION 10: PRODUCT CAPABILITIES & POSITIONING** ⏸️ PENDING

**Duration**: 2-3 hours  
**Status**: ⏸️ **PENDING**  
**Deliverable**: Product capabilities and positioning documentation

#### **What I'll Learn**:
1. **6 Capability Groups**:
   - Clinical Decision Support
   - Research Acceleration
   - Therapeutic Design
   - Platform Intelligence
   - Conversational AI
   - Enterprise Platform
2. **Competitive Advantages**:
   - Transparent reasoning
   - Deterministic confidence
   - Action-ready outputs
   - Multi-modal validation
   - Proactive resistance detection
   - Unified care plans
3. **Customer Value**:
   - For oncologists
   - For biotechs
   - For researchers
   - For clinical trial teams
4. **Product Positioning**:
   - Who we are
   - What we do
   - Why it matters
   - How we're different

#### **How I'll Learn**:
- **Files to Read**:
  - `.cursorrules` lines 219-625 - Product capabilities
  - Product documentation files
  - Marketing materials
- **Codebase Searches**:
  - "What are the product capabilities?"
  - "What is the competitive advantage?"
  - "What is the customer value?"
- **Methods**:
  - Extract product information from .cursorrules
  - Map capabilities to code
  - Document customer value
  - Create positioning summary

#### **Deliverables**:
- 6 capability groups documentation
- Competitive advantages analysis
- Customer value proposition
- Product positioning summary
- Capability-to-code mapping

---

## 📊 EXECUTION STRATEGY

### **Iteration Execution Order**:
1. ✅ **I1**: Complete (Architecture overview)
2. **I2**: Backend Services (Foundation for everything)
3. **I3**: S/P/E Framework (Core intelligence system)
4. **I4**: AI Services (How models are integrated)
5. **I5**: Frontend (How users interact)
6. **I6**: Clinical Systems (Primary use case)
7. **I7**: Research & Design (Secondary use cases)
8. **I8**: Data Flow (How everything connects)
9. **I9**: Patterns & Lessons (Best practices)
10. **I10**: Product Capabilities (Business context)

### **Learning Methods**:
1. **Read Code**: Systematic file reading
2. **Codebase Search**: Semantic searches for understanding
3. **Trace Flows**: Follow data through system
4. **Document**: Write comprehensive documentation
5. **Question**: Ask clarifying questions when needed

### **Quality Criteria**:
- ✅ **Completeness**: All major systems covered
- ✅ **Depth**: Understanding "how" and "why"
- ✅ **Clarity**: Documentation is clear and actionable
- ✅ **Accuracy**: Information is verified against code
- ✅ **Organization**: Information is well-structured

---

## 📝 PROGRESS TRACKING

### **Completed Iterations**:
- ✅ **I1**: Overall Architecture & Core Principles (2-3h) - **COMPLETE**

### **In Progress**:
- None

### **Pending Iterations**:
- ⏸️ **I2**: Backend Services & Orchestration (4-5h)
- ⏸️ **I3**: S/P/E Framework & Efficacy System (3-4h)
- ⏸️ **I4**: AI Services & Model Integration (3-4h)
- ⏸️ **I5**: Frontend Architecture & Components (3-4h)
- ⏸️ **I6**: Clinical Decision Support Systems (3-4h)
- ⏸️ **I7**: Research & Design Capabilities (2-3h)
- ⏸️ **I8**: Data Flow & Integration Patterns (2-3h)
- ⏸️ **I9**: Development Patterns & Lessons (2-3h)
- ⏸️ **I10**: Product Capabilities & Positioning (2-3h)

### **Total Progress**: 1/10 iterations (10% complete)

---

## 🎯 FINAL DELIVERABLE

**Master Learning Document**: `.cursor/ayesha/Deliverables/NYX-v2_COMPLETE_APPLICATION_LEARNING.md`

**Contents**:
- Complete architecture understanding
- All systems documented
- All capabilities explained
- All patterns documented
- All lessons learned captured
- Ready for confident future development

---

**Status**: 📋 **PLAN COMPLETE** - Ready for iterative execution  
**Next Step**: Begin Iteration 2 (Backend Services & Orchestration)

