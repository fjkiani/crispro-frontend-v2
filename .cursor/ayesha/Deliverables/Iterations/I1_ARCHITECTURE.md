# 📊 ITERATION 1: OVERALL ARCHITECTURE & CORE PRINCIPLES

**Status**: ✅ **COMPLETE**  
**Duration**: 2-3 hours  
**Created**: January 14, 2025

---

## 📊 SEQUENCE 1: OVERALL ARCHITECTURE & CORE PRINCIPLES

### **1.1 THREE-TIER BACKEND ARCHITECTURE**

**Key Finding**: The platform has **THREE distinct backend systems** working together:

#### **Tier 1: Minimal Backend** (`oncology-coPilot/oncology-backend-minimal/`)
- **Purpose**: Vercel-deployable demo with production-ready endpoints
- **Status**: ✅ **PRODUCTION** - Real business logic, not mocks
- **Architecture**: FastAPI with modular routers
- **Key Features**:
  - 30+ operational endpoints
  - Modular router pattern (domain-specific routers)
  - Service layer separation (business logic in `services/`)
  - Feature flags for different operational profiles
  - Graceful degradation patterns
  - Complete provenance tracking

#### **Tier 2: Main Backend** (`oncology-coPilot/oncology-backend/`)
- **Purpose**: Full-featured orchestration system with AI agents
- **Status**: ⚠️ **LEGACY** - May have agent system but minimal backend is primary
- **Agents**: 15+ specialized AI agents (if still active)
- **Features**: Complete patient management, workflow orchestration

#### **Tier 3: AI Services Backend** (`src/services/`)
- **Purpose**: Production AI model inference services
- **Status**: ✅ **PRODUCTION** - Real Evo2, AlphaFold 3, Boltz-2 models
- **Deployment**: GPU-powered Modal deployments
- **Capabilities**: Actual biological AI inference

**Architectural Flow**:
```
Frontend (React/Vite)
    ↓
Minimal Backend (FastAPI) - Primary orchestration
    ↓
AI Services (Modal) - Evo2, AlphaFold 3, Boltz-2
    ↓
External APIs (PubMed, ClinVar, cBioPortal, etc.)
```

**Why This Architecture?**
- **Separation of Concerns**: Business logic separate from AI inference
- **Scalability**: AI services on Modal can scale independently
- **Cost Efficiency**: Minimal backend on Vercel (serverless), AI services on-demand
- **Development Speed**: Frontend can be built against minimal backend while AI services develop

---

### **1.2 CORE ARCHITECTURAL PRINCIPLES**

#### **Principle 1: Modular Router Pattern**
- **What**: Each router handles a specific domain (efficacy, insights, design, evidence, etc.)
- **Why**: 
  - Clean separation of concerns
  - Easy to add new capabilities
  - Clear ownership and testing boundaries
- **Example**: `api/routers/efficacy.py`, `api/routers/insights.py`, `api/routers/design.py`

#### **Principle 2: Service Layer Separation**
- **What**: Business logic in `services/`, routers are thin endpoints
- **Why**:
  - Reusable business logic
  - Easier testing (test services independently)
  - Clear API contracts
- **Example**: `api/services/efficacy_orchestrator/`, `api/services/sae_feature_service.py`

#### **Principle 3: Feature Flags**
- **What**: Environment-based toggles for different operational profiles
- **Why**:
  - Graceful degradation (disable features if services unavailable)
  - A/B testing capabilities
  - Demo vs production modes
- **Example**: `EVO_FORCE_MODEL`, `EVO_USE_DELTA_ONLY`, `DISABLE_FUSION`

#### **Principle 4: Graceful Degradation**
- **What**: Fallback chains, placeholder values, non-blocking integration
- **Why**:
  - System remains operational even if external services fail
  - Better user experience (partial results vs complete failure)
  - Resilience to network issues
- **Example**: If Evo2 unavailable → return placeholder scores with provenance

#### **Principle 5: Provenance Tracking**
- **What**: Complete audit trails (run IDs, profiles, methods, citations)
- **Why**:
  - Reproducibility (can rerun exact same analysis)
  - Transparency (users see how results were computed)
  - Compliance (audit trail for clinical use)
- **Example**: Every response includes `provenance` field with run_id, profile, methods

---

### **1.3 KEY TECHNICAL DOCTRINES**

#### **Doctrine 1: The "Wet Noodle" Problem**
- **What**: A DNA sequence that is grammatically correct in 1D (`delta_score`) can still translate into a physically useless protein that fails to fold correctly in 3D (`pLDDT` score)
- **Solution**: Multi-dimensional validation process:
  - **Phase I: The Forge** - Generate candidates
  - **Phase II: The Sieve** - Use sequence-level likelihood scores as fast filter
  - **Phase III: The Gauntlet** - Use 3D structural prediction as final arbiter
- **Why**: Prevents generating biologically invalid sequences

#### **Doctrine 2: The Triumvirate Protocol**
- **What**: Multi-layered approach for variant assessment
- **Components**:
  1. **Truncation Check** - Deterministic bioinformatic translation of CDS
  2. **Evo2 Deep Learning** - Only for non-truncating variants
  3. **Pathway Analysis** - Context-aware impact assessment
- **Why**: Evo2 has blind spot for frameshift/nonsense mutations - need deterministic check first

#### **Doctrine 3: Backend Orchestrator Pattern**
- **What**: Single powerful orchestrator endpoint manages entire multi-stage workflow
- **Why**: 
  - Simplifies frontend (one endpoint vs many)
  - Decouples frontend from backend implementation details
  - Enables progressive enhancement (mock unimplemented services)
- **Example**: `/api/ayesha/complete_care_v2` orchestrates trials + SOC + CA-125 + WIWFM + food + resistance

#### **Doctrine 4: Generative vs Inference Paradigm**
- **What**: Platform uses **generative and predictive paradigm** (`Digital Twin -> Predict -> Generate`)
- **Why**: 
  - Inference-based methods (competitors) deconstruct noisy data to *guess* at reality
  - Generative approach creates ground truth from high-fidelity sequencing
  - More accurate and verifiable
- **Impact**: Revolutionary replacement, not incremental improvement

---

### **1.4 SYSTEM ORGANIZATION**

#### **Backend Structure**:
```
oncology-coPilot/oncology-backend-minimal/
├── api/
│   ├── main.py                    # FastAPI app initialization
│   ├── config.py                  # Feature flags, weights, env vars
│   ├── routers/                   # 30+ domain-specific routers
│   ├── services/                  # Business logic (~100-150 lines each)
│   ├── schemas/                   # Pydantic models for validation
│   └── startup.py                 # Background services
```

#### **Frontend Structure**:
```
oncology-coPilot/oncology-frontend/
├── src/
│   ├── pages/                     # Route-level pages
│   ├── components/                # Reusable React components
│   ├── context/                   # React Context providers
│   └── hooks/                     # Custom React hooks
```

#### **AI Services Structure**:
```
src/services/
├── evo_service/                   # Evo2 inference (Modal)
├── boltz_service/                 # AlphaFold 3/Boltz-2 (Modal)
└── ...
```

---

### **1.5 KEY INSIGHTS FROM .cursorrules**

#### **Product Positioning**:
- **Who**: CrisPRO.ai - AI-powered precision medicine platform
- **What**: Transforms genomic data into actionable therapeutic intelligence
- **Why**: >90% of clinical trials fail - we eliminate uncertainty before expensive lab work
- **How**: Multi-modal AI validation (S/P/E framework) with transparent confidence

#### **Target Customers**:
1. **Oncologists**: Personalized treatment recommendations
2. **Biotechs**: De-risk drug development with in-silico validation
3. **Researchers**: Accelerate hypothesis validation
4. **Clinical Trial Teams**: Match patients to trials with biomarker intelligence

#### **Core Value Proposition**:
- **Transparent Reasoning**: Every recommendation shows WHY, not just WHAT
- **Deterministic Confidence**: 90-100% confidence from checkboxes, not AI magic
- **Action-Ready Outputs**: Clinician-ready dossiers with contacts, checklists, protocols
- **Multi-Modal Validation**: S/P/E framework (70-85% accuracy vs 50-60% single-metric)

---

## 📝 NOTES & QUESTIONS

### **Questions for Next Sequences**:
1. How exactly does the S/P/E framework work? (Sequence 2)
2. How are Evo2 calls made and cached? (Sequence 4)
3. How does the frontend handle sporadic cancer context? (Sequence 3)
4. What are the exact API contracts between services? (Sequence 2)

### **Key Files to Review**:
- `oncology-coPilot/oncology-backend-minimal/api/main.py` - Entry point
- `oncology-coPilot/oncology-backend-minimal/api/config.py` - Feature flags
- `.cursorrules` lines 219-625 - Product capabilities
- `.cursorrules` lines 666-1406 - Technical inventory

---

### **1.6 BACKEND ENTRY POINT ANALYSIS** (`api/main.py`)

#### **Router Registration Pattern**:
- **30+ Routers Registered**: Each router handles a specific domain
- **Conditional Registration**: Some routers only included if feature flags enabled
- **Example Routers**:
  - `health`, `myeloma`, `evo`, `evidence`, `efficacy`
  - `fusion`, `guidance`, `datasets`, `insights`, `design`
  - `sessions`, `auth`, `admin`, `toxicity`, `safety`
  - `metastasis`, `acmg`, `pharmgkb`, `clinical_trials`
  - `trials`, `trials_graph`, `trials_agent`, `resistance`, `nccn`
  - `clinical_genomics`, `kg`, `hypothesis_validator`
  - `ayesha_twin_demo`, `ayesha`, `tumor`, `care`
  - `ayesha_trials`, `ayesha_orchestrator_v2`, `ayesha_dossiers`
  - `dossiers`, `agents`

#### **Startup/Shutdown Tasks**:
- **Background Services**: Calibration preload, agent scheduler
- **Graceful Degradation**: Startup failures don't block app
- **Agent Scheduler**: Autonomous intelligence system

#### **CORS Configuration**:
- **Allow All Origins**: `allow_origins=["*"]` (development mode)
- **Full Access**: All methods and headers allowed

#### **Key Endpoints** (in main.py):
- `/api/workflow/run_seed_soil_analysis` - Metastasis pathway analysis
- `/api/twin/run` - Digital Twin orchestrator (warm model, run scoring)
- `/api/twin/submit` - Async job submission
- `/api/twin/status` - Job status check
- `/api/analytics/dashboard` - Analytics aggregation
- `/api/safety/ensembl_context` - Ensembl context for locus
- `/api/safety/clinvar_context` - ClinVar context for variant
- `/api/patients/{patient_id}` - Patient mutations data
- `/api/research/mutation-analysis` - VUS Explorer analysis

---

### **1.7 CONFIGURATION SYSTEM** (`api/config.py`)

#### **S/P/E Weight Configuration**:
- **Default Weights**: Sequence (0.35), Pathway (0.35), Evidence (0.30)
- **Why**: Balanced multi-modal approach, slightly favoring sequence/pathway
- **Configurable**: Via environment variables

#### **Evidence Gate Thresholds**:
- **Evidence Gate**: 0.7 (conservative clinical default)
- **ClinVar Strong**: 0.8 (high confidence threshold)
- **Pathway Alignment**: 0.2 (minimum pathway match)
- **Insufficient Signal**: 0.02 (below this = insufficient)

#### **Feature Flags System**:
- **Evo2 Control**: `DISABLE_EVO2` - Can disable Evo2 globally
- **Literature Control**: `DISABLE_LITERATURE` - Can disable literature search
- **Fusion Control**: `DISABLE_FUSION` - Can disable AlphaMissense integration
- **API Exposure**: `ENABLE_INSIGHTS_API`, `ENABLE_DESIGN_API`, `ENABLE_COMMAND_CENTER`
- **Operational Mode**: `OPERATIONAL_MODE` - "clinical" or "research"

#### **Evo Spam-Safety Controls**:
- **EVO_SPAM_SAFE**: Default true (prevents excessive API calls)
- **EVO_MAX_MODELS**: Default 1 (if spam-safe), else 3
- **EVO_MAX_FLANKS**: Default 1 (if spam-safe), else 5
- **EVO_DISABLE_TRANSCRIPT_SWEEP**: Default true (if spam-safe)
- **EVO_DISABLE_SYMMETRY**: Default true (if spam-safe)
- **EVO_USE_DELTA_ONLY**: Default true (prevents upstream fan-out)

**Why These Controls?**:
- **Cost Management**: Evo2 calls are expensive (GPU inference)
- **Performance**: Too many parallel calls can timeout
- **Reliability**: Prevents overwhelming external services

#### **Model URL Mapping**:
- **Dynamic Fallback**: Runtime URL resolution with fallback chain
- **Fallback Logic**: `evo2_1b` → `evo2_7b` → `evo2_40b`
- **Why**: Ensures service availability even if preferred model unavailable

#### **Calibration Configuration**:
- **TTL**: 24 hours (how long calibration data is valid)
- **Refresh Interval**: 6 hours (how often to refresh)
- **Preload Genes**: 25 common MM genes preloaded on startup

---

### **1.8 KEY ARCHITECTURAL INSIGHTS**

#### **Why Modular Routers?**
- **Separation of Concerns**: Each router owns a domain (efficacy, insights, design)
- **Easy Testing**: Can test routers independently
- **Clear Ownership**: Easy to find code for specific features
- **Progressive Enhancement**: Can add routers without breaking existing code

#### **Why Feature Flags?**
- **Graceful Degradation**: System works even if services unavailable
- **Demo Mode**: Can disable expensive features for demos
- **A/B Testing**: Can enable features for specific users
- **Cost Control**: Can disable expensive AI services in development

#### **Why Spam-Safety Controls?**
- **Cost**: Evo2 calls are expensive (GPU inference on Modal)
- **Performance**: Too many parallel calls can cause timeouts
- **Reliability**: Prevents overwhelming external services
- **User Experience**: Prevents long wait times from excessive API calls

#### **Why Dynamic Model URL Fallback?**
- **Resilience**: System works even if preferred model unavailable
- **Cost Optimization**: Can use cheaper models (1B) when 40B not needed
- **Flexibility**: Can switch models based on use case

---

**Status**: 🔄 **ITERATION 2 IN PROGRESS** - Backend Services & Orchestration  
**Next**: Continue I2, then move to I3 (S/P/E Framework)

---