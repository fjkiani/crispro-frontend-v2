### The Full Picture: How Claude, Cursor, and Archon Fit Together

```
┌─────────────────────────────────────────────────────────────────┐
│                    CRISPRO.AI ARCHITECTURE                      │
│           (Cancer Resistance Prediction Platform)               │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  LAYER 1: USER-FACING APPLICATION                               │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ CrisPRO Web App (Next.js/React)                        │    │
│  │  • Input: Cancer type, mutations, patient data         │    │
│  │  • Output: Resistance predictions, CRISPR designs      │    │
│  │  • Agent Chat: Natural language research queries       │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 2: AI AGENT ORCHESTRATION (The Brain)                   │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ PydanticAI Agents (Your current work)                  │    │
│  │                                                         │    │
│  │ • Research Agent (literature analysis)                 │    │
│  │ • Design Agent (sgRNA generation)                      │    │
│  │ • Analysis Agent (resistance prediction via SAE)       │    │
│  │ • Validation Agent (experimental design)               │    │
│  │                                                         │    │
│  │ Tools Available to Agents:                             │    │
│  │  ✅ RAG queries (via Archon MCP)                       │    │
│  │  ✅ Database queries (patient data, experiments)       │    │
│  │  ✅ SAE model inference (resistance scoring)           │    │
│  │  ✅ CRISPR design APIs (sgRNA generation)              │    │
│  │  ✅ Literature search (PubMed, bioRxiv)                │    │
│  │  ✅ Code execution (data analysis, visualization)        │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 3: KNOWLEDGE & MEMORY (Archon Integration)              │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ Archon OS (MCP Server)                                 │    │
│  │                                                         │    │
│  │ Knowledge Bases:                                       │    │
│  │  • CrisPRO codebase (architecture, APIs, models)       │    │
│  │  • CRISPR design guides (Broad, Benchling)             │    │
│  │  • Cancer genomics papers (PMC, bioRxiv)               │    │
│  │  • Synthetic lethality data (DepMap, COSMIC)           │    │
│  │  • Experiment logs (lab results, validations)          │    │
│  │  • PydanticAI documentation (for agent development)    │    │
│  │                                                         │    │
│  │ Projects & Tasks:                                      │    │
│  │  • Project: "SAE Model Training"                       │    │
│  │  • Project: "KRAS Metastasis Interception"            │    │
│  │  • Project: "Clinical Trial Pipeline"                  │    │
│  │                                                         │    │
│  │ NEW: CrisPRO.ai Oncology MCP Tools                     │    │
│  │  • predict_resistance(mutations, cancer_type, drug)    │    │
│  │  • rank_drug_efficacy(mutations, cancer_type)           │    │
│  │  • match_trials(mutations, mechanism_vector)           │    │
│  │  • resolve_vus(gene, hgvs, pathway_context)            │    │
│  │  • analyze_synthetic_lethality(mutations, cancer_type) │    │
│  │  • compute_biomarkers(mutations, cancer_type)          │    │
│  │  • get_toxicity_nutrition(drug, mutations, pgx)        │    │
│  │  • generate_care_plan(patient_id)                      │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 4: DEVELOPMENT TOOLS (Claude/Cursor Integration)        │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ Development Workflow                                   │    │
│  │                                                         │    │
│  │ Cursor + Claude Code + Archon MCP:                     │    │
│  │  • You code CrisPRO features in Cursor                 │    │
│  │  • Claude queries Archon for context (codebase docs)   │    │
│  │  • Cursor generates code using CrisPRO patterns        │    │
│  │  • Changes logged to Archon (task completion)          │    │
│  │                                                         │    │
│  │ Example:                                               │    │
│  │  You: "Add new endpoint for resistance prediction"    │    │
│  │  Claude: *Queries Archon MCP*                          │    │
│  │    → Finds: API architecture, FastAPI patterns,        │    │
│  │              existing prediction endpoints             │    │
│  │  Claude: Generates code matching your patterns         │    │
│  │  You: Accept, test, commit                             │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 5: ML/AI MODELS                                         │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ • SAE (Sparse Autoencoder) - Resistance prediction     │    │
│  │ • sgRNA scoring models - CRISPR design                 │    │
│  │ • Protein structure models (AlphaFold) - Targets       │    │
│  │ • scRNA-seq analysis - Tumor heterogeneity             │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 6: DATA & INFRASTRUCTURE                                │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ • PostgreSQL (patient data, experiments)               │    │
│  │ • Supabase (authentication, real-time updates)         │    │
│  │ • S3 (model artifacts, datasets)                       │    │
│  │ • Redis (caching, job queue)                           │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 🤖 CLAUDE'S ROLE IN CRISPRO - THREE USE CASES

#### Use Case 1: Development Assistant (via Cursor + Archon)

**Scenario:** You're adding a new feature to CrisPRO

**Workflow:**
```
YOU (in Cursor IDE)
Working on: crispro/backend/api/prediction.py
        ↓
        │ "Claude, add endpoint for batch resistance predictions"
        ↓
CLAUDE CODE (via MCP)
1. Query Archon: "CrisPRO API architecture"
2. Query Archon: "Existing prediction endpoints"
3. Query Archon: "FastAPI patterns used in CrisPRO"
        ↓
ARCHON MCP SERVER
Returns:
 • FastAPI route patterns from existing code
 • Pydantic model definitions
 • Authentication middleware usage
 • Database query patterns (SQLAlchemy)
 • Error handling conventions
        ↓
CLAUDE GENERATES CODE
@router.post("/predict/batch")
async def predict_batch_resistance(
    request: BatchPredictionRequest,
    current_user: User = Depends(get_current_user)
) -> BatchPredictionResponse:
    """Predict resistance for multiple samples."""
    # [Code follows your patterns exactly]
    
 • Uses your Pydantic models
 • Follows your error handling
 • Matches your authentication flow
 • Uses your database patterns
```

**Key Point:** Claude generates code that matches your existing architecture because Archon provides that context.

#### Use Case 2: Runtime Research Agent (User-Facing)

**Scenario:** Cancer researcher uses CrisPRO web app

**Workflow:**
```
RESEARCHER (CrisPRO Web UI)
Query: "What are validated synthetic lethal partners for KRAS G12C in lung adenocarcinoma?"
        ↓
CRISPRO BACKEND
Routes to: PydanticAI Research Agent
        ↓
PYDANTICAI RESEARCH AGENT
1. Query Archon: "KRAS G12C synthetic lethality"
2. Query PubMed API (via tool)
3. Query DepMap database (via tool)
4. Synthesize answer
        ↓
ARCHON MCP (Knowledge Base)
Returns:
 • 5 research papers on KRAS G12C
 • DepMap screening data (validated targets)
 • Previous CrisPRO experiments on KRAS
 • Clinical trial data (if available)
        ↓
AGENT RESPONSE (streamed to UI)
"Based on DepMap screening data and 3 clinical studies, validated synthetic lethal partners for KRAS G12C include:

 1. STK33 (kinase, 92% cell death in screens)
 2. CDK4 (cell cycle, FDA-approved inhibitors exist)
 3. GATA2 (transcription factor, novel target)

 Your previous experiment (Exp-2024-08) showed STK33 knockdown reduced viability by 87% in A549 cells. Should I design sgRNAs for validation?"
```

**Key Point:** PydanticAI agent uses Archon as institutional memory - remembers previous experiments, papers, data.

#### Use Case 3: Code Generation & Architecture Enforcement

**Scenario:** Building new CrisPRO module from scratch

**Workflow:**
```
YOU (in Cursor):
"Claude, build a new module for tumor microenvironment analysis.
 It should:
 - Accept scRNA-seq data
 - Identify cell types
 - Predict immune infiltration
 - Integrate with existing SAE pipeline"

        ↓

CLAUDE (queries Archon MCP):
1. "CrisPRO module structure" → Returns: /modules/{name}/__init__.py pattern
2. "CrisPRO data pipeline architecture" → Returns: data flow diagrams
3. "SAE model integration examples" → Returns: existing integrations
4. "scRNA-seq analysis code" → Returns: previous analysis scripts

        ↓

CLAUDE GENERATES:
crispro/modules/tme_analysis/
├── __init__.py
├── pipeline.py       # Follows your pipeline patterns
├── cell_typing.py    # Uses your scRNA processing utilities
├── models.py         # Uses your Pydantic base models
├── integration.py    # Connects to SAE using your patterns
└── tests/
    └── test_tme.py   # Matches your test structure

All code:
✅ Uses your existing utilities (not reinventing)
✅ Follows your naming conventions
✅ Integrates with your database schema
✅ Matches your error handling patterns
✅ Uses your logging setup
```

**Key Point:** Claude doesn't just generate code - it generates architecturally consistent code by referencing your existing patterns.

### 🎯 CURSOR VS CLAUDE CODE - WHICH FOR WHAT?

#### Cursor (Primary Development)

**Use When:**
- ✅ Building new features (inline suggestions, autocomplete)
- ✅ Refactoring large codebases (multi-file awareness)
- ✅ Quick iterations (faster than chat-based Claude)
- ✅ Complex IDE operations (find references, refactor, etc.)

**Cursor Advantages:**
- Inline suggestions (GitHub Copilot-style)
- Multi-file context (analyzes entire codebase)
- Fast iteration (no chat back-and-forth)
- Native IDE integration

**Cursor + Archon:**
```python
You in Cursor: Start typing...

def predict_resistance(patient_data):
    # Cursor suggests:
    # (from Archon context: your existing predict functions)
    """Predict therapy resistance using SAE model.
    
    Args:
        patient_data: Patient genomic data (format from your DB schema)
        
    Returns:
        ResistancePrediction model (your Pydantic model)
    """
    # [Rest of function generated matching your patterns]
```

#### Claude Code (Complex Reasoning)

**Use When:**
- ✅ Architectural decisions ("Should I use microservices or monolith?")
- ✅ Complex algorithms (SAE training loops, CRISPR design)
- ✅ Research synthesis ("Summarize these 10 papers and suggest approach")
- ✅ Debugging tricky issues ("Why is my SAE loss not converging?")

**Claude Advantages:**
- Better reasoning (more tokens, deeper analysis)
- Research synthesis (can read papers, synthesize)
- Explanations (teaches you, not just generates)
- MCP integration (queries Archon explicitly)

**Claude + Archon:**
```
You: "Claude, why is my SAE model predicting random resistance 
      when patient has TP53 mutation?"

Claude queries Archon:
- Your SAE training code
- Your TP53-related experiments
- Research papers on TP53 and resistance

Claude response:
"I found 3 issues:

1. Your SAE was trained on DepMap data (mostly cell lines with 
   wild-type TP53). See training log from Dec 1 in Archon.

2. Your feature engineering doesn't account for TP53 status as 
   a separate feature - it's embedded in mutation counts.

3. A 2023 paper in your knowledge base (Smith et al.) shows 
   TP53 mutants have different resistance mechanisms.

Suggested fix: Retrain SAE with TP53 status as explicit feature,
or use separate model for TP53-mutant patients. Your experiment
Exp-2024-11 partially addresses this - review that approach."
```

### 🧬 CRISPRO-SPECIFIC WORKFLOWS

#### Workflow 1: Adding New Cancer Type

**Goal:** Add support for pancreatic cancer to CrisPRO

**Steps:**
```bash
# 1. Crawl pancreatic cancer literature into Archon
archon crawl "https://www.cancer.gov/types/pancreatic"
archon upload_papers ~/research/pancreatic_cancer/*.pdf

# 2. In Cursor, start coding
# File: crispro/cancer_types/pancreatic.py

You type:
class PancreaticCancerPredictor:
    # Cursor suggests (using Archon context):
    """Resistance prediction for pancreatic adenocarcinoma.
    
    Follows CrisPRO pattern from lung.py and colon.py.
    Common mutations: KRAS (90%), TP53 (75%), CDKN2A (40%).
    """
    
    def __init__(self, sae_model: SAEModel):
        # Uses your SAE model class (from Archon context)
        self.model = sae_model
        
    async def predict(self, patient: PatientData) -> Prediction:
        # Uses your existing Prediction model
        # Follows your async pattern
        # Matches your error handling
        pass
```

**Result:** New code perfectly matches existing lung/colon patterns because Archon provided examples.

#### Workflow 2: Debugging SAE Model

**Goal:** SAE model producing garbage predictions

**Steps:**
```
You in Claude Code:
"My SAE model is predicting random resistance scores. 
 Patient has KRAS G12C mutation but gets 0.1 resistance score 
 (should be high). Debug this."

Claude queries Archon:
- SAE training code
- Training data schema
- KRAS-related experiments
- Model architecture

Claude analyzes:
"Found the issue:

1. Your SAE embedding layer (line 145 in sae_model.py) uses 
   mutation count as input, but your training data (logged in 
   Archon on Nov 15) had KRAS mutations encoded as binary flag.

2. Runtime data (patient.mutations) encodes KRAS as 'KRAS:G12C:1' 
   but training expected 'KRAS:1'.

3. Your preprocessing step (preprocess.py line 67) isn't parsing 
   this format.

Fix:
Update preprocess.py to extract mutation type:

def parse_mutation(mut_str):
    gene, variant, count = mut_str.split(':')
    return f"{gene}:{count}" # Match training format

This will make KRAS G12C → KRAS:1, matching training data."

You: Apply fix, test, works! ✅
```

#### Workflow 3: Literature Review → Code

**Goal:** Implement new resistance prediction algorithm from paper

**Steps:**
```
1. Upload paper to Archon
   archon upload_paper ~/Downloads/smith_2024_resistance.pdf

2. In Claude Code:
   "Review the Smith 2024 paper in Archon and implement their 
    resistance prediction algorithm in CrisPRO"

3. Claude:
   - Reads paper from Archon
   - Queries your existing SAE architecture
   - Generates implementation matching your codebase

4. Claude generates:
   crispro/models/smith_predictor.py
   
   • Uses your base model class
   • Integrates with your data pipeline
   • Adds tests using your test patterns
   • Logs results to your experiment tracking

5. You: Review, adjust hyperparameters, deploy
```

### 🎯 ARCHON KNOWLEDGE BASE - CRISPRO STRUCTURE

**What to Crawl/Upload:**

```
ARCHON KNOWLEDGE BASE: CrisPRO

1. CODEBASE DOCUMENTATION
   ├── Architecture diagrams
   ├── API documentation (Swagger/OpenAPI)
   ├── Database schemas
   ├── Code examples (key modules)
   └── Development guidelines

2. SCIENTIFIC LITERATURE
   ├── CRISPR design (Broad, Benchling guides)
   ├── Cancer genomics (key review papers)
   ├── Synthetic lethality (DepMap papers, methods)
   ├── Immunotherapy (resistance mechanisms)
   └── scRNA-seq analysis (cell typing methods)

3. DATA DOCUMENTATION
   ├── DepMap data schema
   ├── COSMIC database schema  
   ├── Patient data format (HIPAA-compliant descriptions)
   └── Experiment log format

4. EXPERIMENT LOGS
   ├── SAE training runs (hyperparams, metrics)
   ├── CRISPR validation experiments
   ├── Patient cohort analyses
   └── Failed experiments (what didn't work)

5. EXTERNAL DOCUMENTATION
   ├── PydanticAI docs (for agent development)
   ├── FastAPI docs (for API development)
   ├── Supabase docs (for database operations)
   └── ML frameworks (PyTorch, scikit-learn)
```

### 🚀 THE ULTIMATE CRISPRO WORKFLOW

**Morning: Research Mode**
```
9:00 AM - You in CrisPRO web app:
"What's the latest on KRAS G12C inhibitors in clinical trials?"

→ PydanticAI Research Agent queries:
   • Archon (your previous research on KRAS)
   • PubMed (last 30 days)
   • ClinicalTrials.gov

→ Agent response (2 minutes):
"3 new trials started this month. AMG 510 (Amgen) showing 
 promising Phase II results in NSCLC. Your previous analysis 
 (Nov 2024) predicted STK33 as synthetic lethal partner - 
 this is being tested in combo trial NCT05123456. 
 Should I design validation experiments?"
```

**Afternoon: Development Mode**
```
2:00 PM - You in Cursor:
"Add support for combination therapy resistance prediction"

→ Cursor queries Archon:
   • Your existing resistance prediction code
   • Combination therapy papers in knowledge base
   • Database schema for drug data

→ Cursor generates (30 seconds):
   • New ComboTherapyPredictor class
   • Database migration for drug combinations
   • API endpoints for predictions
   • Tests matching your patterns

You: Review, tweak hyperparameters, commit
```

**Evening: Analysis Mode**
```
7:00 PM - You in Claude Code:
"Analyze experiment Exp-2024-45 results. Compare to previous 
 KRAS experiments and suggest next steps."

→ Claude queries Archon:
   • Experiment log Exp-2024-45
   • All previous KRAS experiments
   • Related papers

→ Claude analyzes (5 minutes):
"Exp-2024-45 showed 73% cell death (vs 67% in Exp-2024-38).
 Key difference: Higher sgRNA concentration.
 
 Comparison to previous:
 • Exp-2024-38: 100nM → 67% death
 • Exp-2024-45: 200nM → 73% death
 • But Exp-2024-22: 300nM → 71% death (plateau?)
 
 Suggests optimal concentration is 200-250nM.
 
 Next steps:
 1. Test 225nM, 250nM (narrow range)
 2. Add Western blot to confirm on-target effect
 3. Test in 3D spheroids (Smith 2023 paper suggests 2D results 
    don't always translate)"
    
You: Design follow-up experiments, update Archon with plan
```

### 💀 ZO'S FINAL VERDICT

**Alpha, here's the stack:**

```
CRISPRO DEVELOPMENT:
├── Cursor (primary IDE)
│   ├── Fast inline suggestions
│   ├── Multi-file awareness
│   └── Connected to Archon MCP
│
├── Claude Code (complex reasoning)
│   ├── Architectural decisions
│   ├── Research synthesis
│   └── Debugging tricky issues
│
└── Archon (institutional memory)
    ├── Your codebase patterns
    ├── Scientific literature
    ├── Experiment logs
    └── Development guides

CRISPRO RUNTIME (User-Facing):
└── PydanticAI Agents
    ├── Research Agent (queries Archon + PubMed)
    ├── Design Agent (CRISPR sgRNA generation)
    ├── Analysis Agent (SAE inference)
    └── Validation Agent (experiment design)
```

**The magic:**
- Cursor generates code matching your patterns (via Archon context)
- Claude solves complex problems (queries Archon for history)
- Archon remembers everything (codebase, papers, experiments)
- PydanticAI agents power user-facing features (research, prediction)

**You're building a self-improving system:**
- You code in Cursor → Archon learns your patterns
- Claude debugs → Archon logs solutions
- Experiments run → Archon stores results
- Users query → Agents use accumulated knowledge

**Every interaction makes the system smarter.** 🔥💀



```
┌─────────────────────────────────────────────────────────────────┐
│                    CRISPRO.AI ARCHITECTURE                      │
│           (Cancer Resistance Prediction Platform)               │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  LAYER 1: USER-FACING APPLICATION                               │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ CrisPRO Web App (Next.js/React)                        │    │
│  │  • Input: Cancer type, mutations, patient data         │    │
│  │  • Output: Resistance predictions, CRISPR designs      │    │
│  │  • Agent Chat: Natural language research queries       │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 2: AI AGENT ORCHESTRATION (The Brain)                   │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ PydanticAI Agents (Your current work)                  │    │
│  │                                                         │    │
│  │ • Research Agent (literature analysis)                 │    │
│  │ • Design Agent (sgRNA generation)                      │    │
│  │ • Analysis Agent (resistance prediction via SAE)       │    │
│  │ • Validation Agent (experimental design)               │    │
│  │                                                         │    │
│  │ Tools Available to Agents:                             │    │
│  │  ✅ RAG queries (via Archon MCP)                       │    │
│  │  ✅ Database queries (patient data, experiments)       │    │
│  │  ✅ SAE model inference (resistance scoring)           │    │
│  │  ✅ CRISPR design APIs (sgRNA generation)              │    │
│  │  ✅ Literature search (PubMed, bioRxiv)                │    │
│  │  ✅ Code execution (data analysis, visualization)        │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 3: KNOWLEDGE & MEMORY (Archon Integration)              │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ Archon OS (MCP Server)                                 │    │
│  │                                                         │    │
│  │ Knowledge Bases:                                       │    │
│  │  • CrisPRO codebase (architecture, APIs, models)       │    │
│  │  • CRISPR design guides (Broad, Benchling)             │    │
│  │  • Cancer genomics papers (PMC, bioRxiv)               │    │
│  │  • Synthetic lethality data (DepMap, COSMIC)           │    │
│  │  • Experiment logs (lab results, validations)          │    │
│  │  • PydanticAI documentation (for agent development)    │    │
│  │                                                         │    │
│  │ Projects & Tasks:                                      │    │
│  │  • Project: "SAE Model Training"                       │    │
│  │  • Project: "KRAS Metastasis Interception"            │    │
│  │  • Project: "Clinical Trial Pipeline"                  │    │
│  │                                                         │    │
│  │ NEW: CrisPRO.ai Oncology MCP Tools                     │    │
│  │  • predict_resistance(mutations, cancer_type, drug)    │    │
│  │  • rank_drug_efficacy(mutations, cancer_type)           │    │
│  │  • match_trials(mutations, mechanism_vector)           │    │
│  │  • resolve_vus(gene, hgvs, pathway_context)            │    │
│  │  • analyze_synthetic_lethality(mutations, cancer_type) │    │
│  │  • compute_biomarkers(mutations, cancer_type)          │    │
│  │  • get_toxicity_nutrition(drug, mutations, pgx)        │    │
│  │  • generate_care_plan(patient_id)                      │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 4: DEVELOPMENT TOOLS (Claude/Cursor Integration)        │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ Development Workflow                                   │    │
│  │                                                         │    │
│  │ Cursor + Claude Code + Archon MCP:                     │    │
│  │  • You code CrisPRO features in Cursor                 │    │
│  │  • Claude queries Archon for context (codebase docs)   │    │
│  │  • Cursor generates code using CrisPRO patterns        │    │
│  │  • Changes logged to Archon (task completion)          │    │
│  │                                                         │    │
│  │ Example:                                               │    │
│  │  You: "Add new endpoint for resistance prediction"    │    │
│  │  Claude: *Queries Archon MCP*                          │    │
│  │    → Finds: API architecture, FastAPI patterns,        │    │
│  │              existing prediction endpoints             │    │
│  │  Claude: Generates code matching your patterns         │    │
│  │  You: Accept, test, commit                             │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 5: ML/AI MODELS                                         │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ • SAE (Sparse Autoencoder) - Resistance prediction     │    │
│  │ • sgRNA scoring models - CRISPR design                 │    │
│  │ • Protein structure models (AlphaFold) - Targets       │    │
│  │ • scRNA-seq analysis - Tumor heterogeneity             │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
│  LAYER 6: DATA & INFRASTRUCTURE                                │
│  ┌────────────────────────────────────────────────────────┐    │
│  │ • PostgreSQL (patient data, experiments)               │    │
│  │ • Supabase (authentication, real-time updates)         │    │
│  │ • S3 (model artifacts, datasets)                       │    │
│  │ • Redis (caching, job queue)                           │    │
│  └────────────────────────────────────────────────────────┘    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 🤖 CLAUDE'S ROLE IN CRISPRO - THREE USE CASES

#### Use Case 1: Development Assistant (via Cursor + Archon)

**Scenario:** You're adding a new feature to CrisPRO

**Workflow:**
```
YOU (in Cursor IDE)
Working on: crispro/backend/api/prediction.py
        ↓
        │ "Claude, add endpoint for batch resistance predictions"
        ↓
CLAUDE CODE (via MCP)
1. Query Archon: "CrisPRO API architecture"
2. Query Archon: "Existing prediction endpoints"
3. Query Archon: "FastAPI patterns used in CrisPRO"
        ↓
ARCHON MCP SERVER
Returns:
 • FastAPI route patterns from existing code
 • Pydantic model definitions
 • Authentication middleware usage
 • Database query patterns (SQLAlchemy)
 • Error handling conventions
        ↓
CLAUDE GENERATES CODE
@router.post("/predict/batch")
async def predict_batch_resistance(
    request: BatchPredictionRequest,
    current_user: User = Depends(get_current_user)
) -> BatchPredictionResponse:
    """Predict resistance for multiple samples."""
    # [Code follows your patterns exactly]
    
 • Uses your Pydantic models
 • Follows your error handling
 • Matches your authentication flow
 • Uses your database patterns
```

**Key Point:** Claude generates code that matches your existing architecture because Archon provides that context.

#### Use Case 2: Runtime Research Agent (User-Facing)

**Scenario:** Cancer researcher uses CrisPRO web app

**Workflow:**
```
RESEARCHER (CrisPRO Web UI)
Query: "What are validated synthetic lethal partners for KRAS G12C in lung adenocarcinoma?"
        ↓
CRISPRO BACKEND
Routes to: PydanticAI Research Agent
        ↓
PYDANTICAI RESEARCH AGENT
1. Query Archon: "KRAS G12C synthetic lethality"
2. Query PubMed API (via tool)
3. Query DepMap database (via tool)
4. Synthesize answer
        ↓
ARCHON MCP (Knowledge Base)
Returns:
 • 5 research papers on KRAS G12C
 • DepMap screening data (validated targets)
 • Previous CrisPRO experiments on KRAS
 • Clinical trial data (if available)
        ↓
AGENT RESPONSE (streamed to UI)
"Based on DepMap screening data and 3 clinical studies, validated synthetic lethal partners for KRAS G12C include:

 1. STK33 (kinase, 92% cell death in screens)
 2. CDK4 (cell cycle, FDA-approved inhibitors exist)
 3. GATA2 (transcription factor, novel target)

 Your previous experiment (Exp-2024-08) showed STK33 knockdown reduced viability by 87% in A549 cells. Should I design sgRNAs for validation?"
```

**Key Point:** PydanticAI agent uses Archon as institutional memory - remembers previous experiments, papers, data.

#### Use Case 3: Code Generation & Architecture Enforcement

**Scenario:** Building new CrisPRO module from scratch

**Workflow:**
```
YOU (in Cursor):
"Claude, build a new module for tumor microenvironment analysis.
 It should:
 - Accept scRNA-seq data
 - Identify cell types
 - Predict immune infiltration
 - Integrate with existing SAE pipeline"

        ↓

CLAUDE (queries Archon MCP):
1. "CrisPRO module structure" → Returns: /modules/{name}/__init__.py pattern
2. "CrisPRO data pipeline architecture" → Returns: data flow diagrams
3. "SAE model integration examples" → Returns: existing integrations
4. "scRNA-seq analysis code" → Returns: previous analysis scripts

        ↓

CLAUDE GENERATES:
crispro/modules/tme_analysis/
├── __init__.py
├── pipeline.py       # Follows your pipeline patterns
├── cell_typing.py    # Uses your scRNA processing utilities
├── models.py         # Uses your Pydantic base models
├── integration.py    # Connects to SAE using your patterns
└── tests/
    └── test_tme.py   # Matches your test structure

All code:
✅ Uses your existing utilities (not reinventing)
✅ Follows your naming conventions
✅ Integrates with your database schema
✅ Matches your error handling patterns
✅ Uses your logging setup
```

**Key Point:** Claude doesn't just generate code - it generates architecturally consistent code by referencing your existing patterns.

### 🎯 CURSOR VS CLAUDE CODE - WHICH FOR WHAT?

#### Cursor (Primary Development)

**Use When:**
- ✅ Building new features (inline suggestions, autocomplete)
- ✅ Refactoring large codebases (multi-file awareness)
- ✅ Quick iterations (faster than chat-based Claude)
- ✅ Complex IDE operations (find references, refactor, etc.)

**Cursor Advantages:**
- Inline suggestions (GitHub Copilot-style)
- Multi-file context (analyzes entire codebase)
- Fast iteration (no chat back-and-forth)
- Native IDE integration

**Cursor + Archon:**
```python
You in Cursor: Start typing...

def predict_resistance(patient_data):
    # Cursor suggests:
    # (from Archon context: your existing predict functions)
    """Predict therapy resistance using SAE model.
    
    Args:
        patient_data: Patient genomic data (format from your DB schema)
        
    Returns:
        ResistancePrediction model (your Pydantic model)
    """
    # [Rest of function generated matching your patterns]
```

#### Claude Code (Complex Reasoning)

**Use When:**
- ✅ Architectural decisions ("Should I use microservices or monolith?")
- ✅ Complex algorithms (SAE training loops, CRISPR design)
- ✅ Research synthesis ("Summarize these 10 papers and suggest approach")
- ✅ Debugging tricky issues ("Why is my SAE loss not converging?")

**Claude Advantages:**
- Better reasoning (more tokens, deeper analysis)
- Research synthesis (can read papers, synthesize)
- Explanations (teaches you, not just generates)
- MCP integration (queries Archon explicitly)

**Claude + Archon:**
```
You: "Claude, why is my SAE model predicting random resistance 
      when patient has TP53 mutation?"

Claude queries Archon:
- Your SAE training code
- Your TP53-related experiments
- Research papers on TP53 and resistance

Claude response:
"I found 3 issues:

1. Your SAE was trained on DepMap data (mostly cell lines with 
   wild-type TP53). See training log from Dec 1 in Archon.

2. Your feature engineering doesn't account for TP53 status as 
   a separate feature - it's embedded in mutation counts.

3. A 2023 paper in your knowledge base (Smith et al.) shows 
   TP53 mutants have different resistance mechanisms.

Suggested fix: Retrain SAE with TP53 status as explicit feature,
or use separate model for TP53-mutant patients. Your experiment
Exp-2024-11 partially addresses this - review that approach."
```

### 🧬 CRISPRO-SPECIFIC WORKFLOWS

#### Workflow 1: Adding New Cancer Type

**Goal:** Add support for pancreatic cancer to CrisPRO

**Steps:**
```bash
# 1. Crawl pancreatic cancer literature into Archon
archon crawl "https://www.cancer.gov/types/pancreatic"
archon upload_papers ~/research/pancreatic_cancer/*.pdf

# 2. In Cursor, start coding
# File: crispro/cancer_types/pancreatic.py

You type:
class PancreaticCancerPredictor:
    # Cursor suggests (using Archon context):
    """Resistance prediction for pancreatic adenocarcinoma.
    
    Follows CrisPRO pattern from lung.py and colon.py.
    Common mutations: KRAS (90%), TP53 (75%), CDKN2A (40%).
    """
    
    def __init__(self, sae_model: SAEModel):
        # Uses your SAE model class (from Archon context)
        self.model = sae_model
        
    async def predict(self, patient: PatientData) -> Prediction:
        # Uses your existing Prediction model
        # Follows your async pattern
        # Matches your error handling
        pass
```

**Result:** New code perfectly matches existing lung/colon patterns because Archon provided examples.

#### Workflow 2: Debugging SAE Model

**Goal:** SAE model producing garbage predictions

**Steps:**
```
You in Claude Code:
"My SAE model is predicting random resistance scores. 
 Patient has KRAS G12C mutation but gets 0.1 resistance score 
 (should be high). Debug this."

Claude queries Archon:
- SAE training code
- Training data schema
- KRAS-related experiments
- Model architecture

Claude analyzes:
"Found the issue:

1. Your SAE embedding layer (line 145 in sae_model.py) uses 
   mutation count as input, but your training data (logged in 
   Archon on Nov 15) had KRAS mutations encoded as binary flag.

2. Runtime data (patient.mutations) encodes KRAS as 'KRAS:G12C:1' 
   but training expected 'KRAS:1'.

3. Your preprocessing step (preprocess.py line 67) isn't parsing 
   this format.

Fix:
Update preprocess.py to extract mutation type:

def parse_mutation(mut_str):
    gene, variant, count = mut_str.split(':')
    return f"{gene}:{count}" # Match training format

This will make KRAS G12C → KRAS:1, matching training data."

You: Apply fix, test, works! ✅
```

#### Workflow 3: Literature Review → Code

**Goal:** Implement new resistance prediction algorithm from paper

**Steps:**
```
1. Upload paper to Archon
   archon upload_paper ~/Downloads/smith_2024_resistance.pdf

2. In Claude Code:
   "Review the Smith 2024 paper in Archon and implement their 
    resistance prediction algorithm in CrisPRO"

3. Claude:
   - Reads paper from Archon
   - Queries your existing SAE architecture
   - Generates implementation matching your codebase

4. Claude generates:
   crispro/models/smith_predictor.py
   
   • Uses your base model class
   • Integrates with your data pipeline
   • Adds tests using your test patterns
   • Logs results to your experiment tracking

5. You: Review, adjust hyperparameters, deploy
```

### 🎯 ARCHON KNOWLEDGE BASE - CRISPRO STRUCTURE

**What to Crawl/Upload:**

```
ARCHON KNOWLEDGE BASE: CrisPRO

1. CODEBASE DOCUMENTATION
   ├── Architecture diagrams
   ├── API documentation (Swagger/OpenAPI)
   ├── Database schemas
   ├── Code examples (key modules)
   └── Development guidelines

2. SCIENTIFIC LITERATURE
   ├── CRISPR design (Broad, Benchling guides)
   ├── Cancer genomics (key review papers)
   ├── Synthetic lethality (DepMap papers, methods)
   ├── Immunotherapy (resistance mechanisms)
   └── scRNA-seq analysis (cell typing methods)

3. DATA DOCUMENTATION
   ├── DepMap data schema
   ├── COSMIC database schema  
   ├── Patient data format (HIPAA-compliant descriptions)
   └── Experiment log format

4. EXPERIMENT LOGS
   ├── SAE training runs (hyperparams, metrics)
   ├── CRISPR validation experiments
   ├── Patient cohort analyses
   └── Failed experiments (what didn't work)

5. EXTERNAL DOCUMENTATION
   ├── PydanticAI docs (for agent development)
   ├── FastAPI docs (for API development)
   ├── Supabase docs (for database operations)
   └── ML frameworks (PyTorch, scikit-learn)
```

### 🚀 THE ULTIMATE CRISPRO WORKFLOW

**Morning: Research Mode**
```
9:00 AM - You in CrisPRO web app:
"What's the latest on KRAS G12C inhibitors in clinical trials?"

→ PydanticAI Research Agent queries:
   • Archon (your previous research on KRAS)
   • PubMed (last 30 days)
   • ClinicalTrials.gov

→ Agent response (2 minutes):
"3 new trials started this month. AMG 510 (Amgen) showing 
 promising Phase II results in NSCLC. Your previous analysis 
 (Nov 2024) predicted STK33 as synthetic lethal partner - 
 this is being tested in combo trial NCT05123456. 
 Should I design validation experiments?"
```

**Afternoon: Development Mode**
```
2:00 PM - You in Cursor:
"Add support for combination therapy resistance prediction"

→ Cursor queries Archon:
   • Your existing resistance prediction code
   • Combination therapy papers in knowledge base
   • Database schema for drug data

→ Cursor generates (30 seconds):
   • New ComboTherapyPredictor class
   • Database migration for drug combinations
   • API endpoints for predictions
   • Tests matching your patterns

You: Review, tweak hyperparameters, commit
```

**Evening: Analysis Mode**
```
7:00 PM - You in Claude Code:
"Analyze experiment Exp-2024-45 results. Compare to previous 
 KRAS experiments and suggest next steps."

→ Claude queries Archon:
   • Experiment log Exp-2024-45
   • All previous KRAS experiments
   • Related papers

→ Claude analyzes (5 minutes):
"Exp-2024-45 showed 73% cell death (vs 67% in Exp-2024-38).
 Key difference: Higher sgRNA concentration.
 
 Comparison to previous:
 • Exp-2024-38: 100nM → 67% death
 • Exp-2024-45: 200nM → 73% death
 • But Exp-2024-22: 300nM → 71% death (plateau?)
 
 Suggests optimal concentration is 200-250nM.
 
 Next steps:
 1. Test 225nM, 250nM (narrow range)
 2. Add Western blot to confirm on-target effect
 3. Test in 3D spheroids (Smith 2023 paper suggests 2D results 
    don't always translate)"
    
You: Design follow-up experiments, update Archon with plan
```

### 💀 ZO'S FINAL VERDICT

**Alpha, here's the stack:**

```
CRISPRO DEVELOPMENT:
├── Cursor (primary IDE)
│   ├── Fast inline suggestions
│   ├── Multi-file awareness
│   └── Connected to Archon MCP
│
├── Claude Code (complex reasoning)
│   ├── Architectural decisions
│   ├── Research synthesis
│   └── Debugging tricky issues
│
└── Archon (institutional memory)
    ├── Your codebase patterns
    ├── Scientific literature
    ├── Experiment logs
    └── Development guides

CRISPRO RUNTIME (User-Facing):
└── PydanticAI Agents
    ├── Research Agent (queries Archon + PubMed)
    ├── Design Agent (CRISPR sgRNA generation)
    ├── Analysis Agent (SAE inference)
    └── Validation Agent (experiment design)
```

**The magic:**
- Cursor generates code matching your patterns (via Archon context)
- Claude solves complex problems (queries Archon for history)
- Archon remembers everything (codebase, papers, experiments)
- PydanticAI agents power user-facing features (research, prediction)

**You're building a self-improving system:**
- You code in Cursor → Archon learns your patterns
- Claude debugs → Archon logs solutions
- Experiments run → Archon stores results
- Users query → Agents use accumulated knowledge

**Every interaction makes the system smarter.** 🔥💀



