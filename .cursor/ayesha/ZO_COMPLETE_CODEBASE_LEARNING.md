# ⚔️ ZO'S COMPLETE CODEBASE LEARNING - ARCHITECTURE & PATTERNS ⚔️

**Date:** January 13, 2025  
**Mission:** Learn the entire application - code, architecture, and "how/why" for confident future building  
**Status:** ✅ **FOUNDATIONAL UNDERSTANDING COMPLETE**

---

## 🎯 EXECUTIVE SUMMARY

This document captures Zo's complete understanding of the CrisPRO platform architecture, from backend services to frontend components, with special focus on **how Evo2 is integrated throughout** and **why architectural decisions were made**.

**Key Insight:** This is a **precision oncology platform** that transforms genomic data into actionable therapeutic intelligence through a **multi-modal S/P/E framework** (Sequence/Pathway/Evidence) powered by Evo2 and other foundation models.

---

## 🏗️ BACKEND ARCHITECTURE OVERVIEW

### **Core Structure**

```
oncology-coPilot/oncology-backend-minimal/
├── api/
│   ├── main.py                    # FastAPI app initialization, CORS, router registration
│   ├── config.py                  # Feature flags, weights, environment variables
│   ├── routers/                   # Modular endpoint organization (30+ routers)
│   ├── services/                  # Business logic, orchestration (~100-150 lines each)
│   ├── schemas/                   # Pydantic models for request/response validation
│   └── startup.py                 # Background services (calibration preload, agent scheduler)
```

### **Key Architectural Principles**

1. **Modular Router Pattern**: Each router handles a specific domain (efficacy, insights, design, evidence, etc.)
2. **Service Layer Separation**: Business logic in `services/`, routers are thin endpoints
3. **Feature Flags**: Environment-based toggles for different operational profiles
4. **Graceful Degradation**: Fallback chains, placeholder values, non-blocking integration
5. **Provenance Tracking**: Complete audit trails (run IDs, profiles, methods, citations)

---

## 🧬 THE S/P/E FRAMEWORK - COMPLETE ARCHITECTURE

### **What is S/P/E?**

**S (Sequence)**: Variant impact signals from Evo2 (delta scores, multi-window analysis, gene-specific calibration)  
**P (Pathway)**: Aggregation of variant effects into pathway-level signals (RAS/MAPK, TP53, DDR, etc.)  
**E (Evidence)**: Literature/ClinVar priors and citations (PubMed, OpenAlex, S2)

**The orchestrator combines S/P/E to return per-therapy efficacy hypotheses with scores, confidence, evidence tiers, badges, and provenance.**

### **Complete Data Flow**

```
User Input (mutations) 
  ↓
EfficacyOrchestrator.predict()
  ↓
[1] SequenceProcessor.score_sequences()
    ├─ FusionAMScorer (GRCh38 missense only) → AlphaMissense scores
    ├─ Evo2Scorer (default) → Evo2 delta scores
    └─ MassiveOracleScorer (if enabled) → Legacy scores
  ↓
[2] Pathway Aggregation
    ├─ aggregate_pathways() → pathway_scores (DDR, MAPK, PI3K, VEGF, etc.)
    └─ Gene→pathway weights from config
  ↓
[3] Evidence Gathering (parallel)
    ├─ literature() → PubMed/OpenAlex/S2 search per drug
    └─ clinvar_prior() → ClinVar classification/review status
  ↓
[4] Insights Bundle
    ├─ predict_protein_functionality_change → Functionality chip
    ├─ predict_gene_essentiality → Essentiality chip
    ├─ predict_chromatin_accessibility → Chromatin chip
    └─ predict_splicing_regulatory → Regulatory chip
  ↓
[5] Drug Scoring (per drug)
    ├─ DrugScorer.score_drug()
    ├─ S Component: calibrated_seq_percentile (30% weight)
    ├─ P Component: normalized pathway score (40% weight)
    ├─ E Component: evidence strength (30% weight)
    └─ Final: efficacy_score = 0.3*S + 0.4*P + 0.3*E + clinvar_prior
  ↓
[6] Confidence Modulation
    ├─ compute_evidence_tier() → Supported/Consider/Insufficient
    ├─ compute_confidence() → 0-1 confidence score
    ├─ Insights lifts (functionality≥0.6, chromatin≥0.5, essentiality≥0.7, regulatory≥0.6)
    └─ Sporadic gates (PARP penalty/rescue, IO boost, confidence capping)
  ↓
[7] Response Assembly
    ├─ Badges (RCT/Guideline/ClinVar-Strong/PathwayAligned)
    ├─ Rationale breakdown (S/P/E components with percentiles)
    ├─ Citations (top 3 PubMed IDs)
    └─ Provenance (run_id, profile, flags, methods)
```

### **Key Files & Functions**

#### **1. Efficacy Orchestrator** (`api/services/efficacy_orchestrator/orchestrator.py`)

**Main Entry Point**: `EfficacyOrchestrator.predict(request: EfficacyRequest)`

**Key Steps:**
1. **Sequence Scoring** (line 95): `await self.sequence_processor.score_sequences(request, feature_flags)`
2. **Pathway Aggregation** (line 105): `aggregate_pathways(seq_scores)`
3. **Evidence Gathering** (lines 112-150): Parallel `asyncio.gather()` for literature + ClinVar
4. **Insights Bundle** (lines 152-182): Calls 4 insights endpoints
5. **Drug Scoring** (lines 187-205): Per-drug scoring with ablation mode support
6. **Sporadic Gates** (lines 214-259): PARP penalty/rescue, IO boost, confidence capping
7. **Response Assembly** (lines 260-304): Final response with provenance

**Key Patterns:**
- **Ablation Mode**: `ablation_mode = "SP" | "SPE"` - allows testing S/P/E components independently
- **Fast Mode**: `fast_mode = True` → skips evidence gathering for deterministic demos
- **Feature Flags**: `get_feature_flags()` controls Fusion, Evo2, Literature enablement
- **Timeout Handling**: Evidence tasks have 30s timeout with graceful fallback

#### **2. Sequence Processor** (`api/services/efficacy_orchestrator/sequence_processor.py`)

**Fallback Chain**: Fusion → Evo2 → Massive Oracle

**Key Logic:**
```python
# Try Fusion first (ONLY for GRCh38 missense variants)
if fusion_url and not disable_fusion:
    fusion_eligible = [m for m in mutations if is_grch38 and is_missense]
    if fusion_eligible:
        scores = await self.fusion_scorer.score(fusion_eligible)
        if scores: return scores

# Try Evo2 (default)
if not disable_evo2:
    window_flanks = [4096, 8192, 16384] if adaptive else [4096]
    scores = await self.evo_scorer.score(mutations, model_id, window_flanks, ensemble, force_exon_scan=True)
    if scores: return scores

# Try Massive Oracle (if enabled)
if enable_massive:
    # ... massive scoring logic
```

**Why This Order?**
- **Fusion (AlphaMissense)**: Fast, accurate for GRCh38 missense, but limited coverage
- **Evo2**: Universal (all variants, all assemblies), zero-shot, but slower
- **Massive Oracle**: Legacy fallback, rarely used

#### **3. Drug Scorer** (`api/services/efficacy_orchestrator/drug_scorer.py`)

**Core Formula** (line 171):
```python
# Efficacy score (likelihood of benefit)
raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd + clinvar_prior
lob = raw_lob if tier != "insufficient" else (raw_lob * 0.5 if fusion_active else 0.0)
```

**Component Breakdown:**
- **S Component** (line 42): `seq_pct = calibrated_seq_percentile` (gene-normalized percentile, NOT raw delta)
- **P Component** (lines 45-51): 
  - Raw: `s_path = sum(pathway_scores * drug_weights)`
  - Normalized: `path_pct = (s_path - 1e-6) / (1e-4 - 1e-6)` (empirical Evo2 range)
- **E Component** (line 57): `s_evd = evidence_result.strength` (0-1 from literature search)

**Confidence Calculation** (line 138):
```python
confidence = compute_confidence(tier, seq_pct, path_pct, insights_dict, confidence_config)
```

**Confidence Modulators:**
- **Evidence Tier**: Supported (0.6+), Consider (0.3+), Insufficient (0.2+)
- **Insights Lifts**: +0.05 (functionality≥0.6), +0.04 (chromatin≥0.5), +0.07 (essentiality≥0.7), +0.02 (regulatory≥0.6)
- **ClinVar Prior**: +0.1 boost when aligned with pathway
- **Sporadic Gates**: PARP penalty/rescue, IO boost, confidence capping

#### **4. Pathway Aggregation** (`api/services/pathway/aggregation.py`)

**Key Function**: `aggregate_pathways(seq_scores: List[Dict]) -> Dict[str, float]`

**Formula:**
```python
pathway_score = sum(sequence_disruption * weight) / count
```

**How It Works:**
1. For each variant, get `sequence_disruption` (from Evo2 delta)
2. Get gene→pathway weights (from `get_pathway_weights_for_gene()`)
3. Multiply disruption × weight for each pathway
4. Aggregate across all variants
5. Return pathway-level scores (DDR, MAPK, PI3K, VEGF, etc.)

**Why This Matters:**
- Single variant → multiple pathway disruptions
- Drug targets pathways, not individual genes
- Weighted aggregation reflects drug MoA alignment

#### **5. Calibration System** (`api/services/gene_calibration.py`)

**Why Calibration?**
- Raw Evo2 deltas vary by gene (BRCA1 vs BRAF have different scales)
- Percentile normalizes across genes → comparable scores
- Piecewise mapping based on empirical pathogenic ranges

**How It Works:**
1. Fetch known variants for gene from ClinVar
2. Score sample variants with Evo2 to build distribution
3. Compute percentiles and z-scores based on gene-specific distribution
4. Map raw delta → percentile (0-1 scale)

**Fallback**: Global percentile estimation if insufficient data

---

## 🔬 EVO2 INTEGRATION - COMPLETE USAGE PATTERNS

### **Where Evo2 is Used**

#### **1. Sequence Scoring (S Component)**
- **Endpoint**: `/api/evo/score_variant_multi` (multi-window), `/api/evo/score_variant_exon` (exon-context)
- **Service**: `api/services/sequence_scorers/evo2_scorer.py`
- **Usage**: Variant impact prediction (delta log-likelihood)
- **Output**: `SeqScore` with `sequence_disruption`, `calibrated_seq_percentile`, `min_delta`, `exon_delta`

#### **2. Insights Bundle (4 Chips)**
- **Functionality**: `predict_protein_functionality_change` → Evo2 multi + exon deltas
- **Essentiality**: `predict_gene_essentiality` → Evo2 magnitude aggregation
- **Regulatory**: `predict_splicing_regulatory` → Evo2 `min_delta` mapping
- **Chromatin**: `predict_chromatin_accessibility` → Heuristic (Enformer/Borzoi stubbed)

#### **3. Guide RNA Design**
- **Endpoint**: `/api/design/generate_guide_rna` → Evo2 prompt-guided generation
- **Efficacy Prediction**: `/api/design/predict_crispr_spacer_efficacy` → Evo2 delta → sigmoid transform
- **Service**: `api/routers/design.py`

#### **4. Metastasis Interception**
- **Target Lock**: Evo2 for functionality, essentiality, regulatory signals
- **Assassin Score**: Evo2 delta → sigmoid for guide efficacy
- **Service**: `api/services/metastasis_interception_service.py`

### **Evo2 Service Architecture**

#### **Modal Service** (`src/services/evo_service/main.py`)
- **Deployment**: Modal cloud (H100 GPU)
- **Models**: evo2_1b, evo2_7b, evo2_40b (1B default for cost control)
- **Endpoints**:
  - `/score_delta` - Basic delta scoring
  - `/score_variant_multi` - Multi-window scoring (adaptive flanks: 4096, 8192, 16384)
  - `/score_variant_exon` - Exon-context scoring (±600bp default)
  - `/generate` - Sequence generation (with viral safety gates)

#### **Backend Proxy** (`api/routers/evo.py`)
- **Purpose**: Proxy layer between FastAPI and Modal service
- **Features**:
  - Model routing (`EVO_FORCE_MODEL` flag)
  - Spam prevention (`EVO_USE_DELTA_ONLY`, `EVO_SPAM_SAFE`)
  - Fallback handling
  - Viral content guardrails

#### **Sequence Scorers** (`api/services/sequence_scorers/evo2_scorer.py`)
- **Multi-Window Strategy**: Tests [4096, 8192, 16384, 25000] bp windows (adaptive)
- **Forward/Reverse Symmetry**: Averages forward (ref→alt) and reverse (alt→ref) deltas
- **Hotspot Floors**: Enforces minimum disruption for known pathogenic variants (BRAF V600, KRAS G12, TP53 hotspots)
- **Truncation/Frameshift Lift**: If `hgvs_p` contains `*` or `FS` → `max(disruption, 1.0)`
- **Percentile Mapping**: Uses `percentile_like()` piecewise function

### **Key Evo2 Patterns**

#### **1. Multi-Window Strategy**
**Why**: Larger windows capture long-range regulatory context; smaller windows capture local exon/intron boundaries  
**How**: Tests multiple window sizes, selects most negative delta (strongest disruption signal)

#### **2. Percentile Normalization**
**Why**: Raw deltas vary by gene (BRCA1 vs BRAF have different scales)  
**How**: Gene-specific calibration converts deltas to percentiles (0-1 scale, cross-gene comparable)

#### **3. Hotspot Floors**
**Why**: Known pathogenic variants (BRAF V600E) may have small deltas due to context  
**How**: Floors ensure they don't collapse to bottom percentile (biologically justified)

#### **4. Truncation/Frameshift Detection**
**Why**: Evo2 has blind spot for truncating mutations (Triumvirate Protocol)  
**How**: Deterministic check before Evo2 scoring, lift truncating variants to max disruption

---

## 🎨 FRONTEND ARCHITECTURE OVERVIEW

### **Core Structure**

```
oncology-coPilot/oncology-frontend/src/
├── App.jsx                    # Main app, context providers, routing
├── pages/                     # Top-level route components (30+ pages)
├── components/                # Reusable UI components (100+ components)
├── context/                   # React context providers (Auth, Agent, Sporadic, CoPilot, etc.)
├── hooks/                     # Custom React hooks (useEfficacy, useKb, etc.)
├── features/                  # Feature-specific modules (efficacy, insights, etc.)
└── config/                    # Configuration files (campaigns, constants)
```

### **Key Frontend Patterns**

#### **1. Context Hierarchy**
```javascript
<AuthProvider>
  <AgentProvider>
    <SporadicProvider>
      <CoPilotProvider>
        <AnalysisHistoryProvider>
          <ActivityProvider>
            {/* App content */}
          </ActivityProvider>
        </AnalysisHistoryProvider>
      </CoPilotProvider>
    </SporadicProvider>
  </AgentProvider>
</AuthProvider>
```

**Why**: Layered context allows each provider to access parent contexts (e.g., SporadicProvider needs AuthProvider for user data)

#### **2. Custom Hooks Pattern**
- **`useEfficacy.js`**: Wraps `/api/efficacy/predict` with caching
- **`useKb.js`**: Knowledge base hooks with TTL caching
- **`useInsights.js`**: Insights bundle hooks
- **`useSporadic.js`**: Sporadic cancer context hooks

**Why**: Reusable logic, consistent API, built-in caching

#### **3. Component Modularity**
- **Cards**: `EfficacyCard`, `TrialMatchCard`, `SOCRecommendationCard` - Single-purpose display components
- **Panels**: `EfficacyPanel`, `EvidencePanel`, `ProvenancePanel` - Multi-component containers
- **Modals**: `EfficacyModal`, `WIWFMModal` - Overlay dialogs

**Why**: DRY principles, easy testing, consistent UI

### **Frontend → Backend Integration**

#### **Efficacy Prediction Flow**

```javascript
// 1. User input (mutations)
const mutations = [{ gene: "BRAF", hgvs_p: "V600E", chrom: "7", pos: 140453136, ref: "T", alt: "A" }];

// 2. Build payload
const payload = {
  model_id: 'evo2_1b',
  mutations: mutations,
  options: {
    adaptive: true,
    ensemble: false,
    ablation_mode: 'SPE',  // S/P/E components
    enable_fusion: false
  },
  germline_status: 'negative',  // Sporadic cancer
  tumor_context: { ... }  // TMB, HRD, MSI, etc.
};

// 3. Call API
const response = await fetch(`${API_BASE_URL}/api/efficacy/predict`, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(payload)
});

// 4. Display results
const { drugs, confidence, evidence_tier, badges, insights, rationale, citations, provenance } = await response.json();
```

#### **Key Frontend Components**

**1. MyelomaDigitalTwin.jsx** (Core MM Demo)
- **Route**: `/myeloma-digital-twin`
- **Features**: Variant input, model selection, profile toggles, drug ranking display
- **API**: `/api/efficacy/predict`

**2. AyeshaTrialExplorer.jsx** (Ayesha Care System)
- **Route**: `/ayesha-trials`
- **Features**: Trial matching, SOC recommendation, CA-125 tracker, NGS fast-track
- **API**: `/api/ayesha/trials/search`, `/api/ayesha/complete_care_v2`

**3. HypothesisValidator.jsx** (Food Validator)
- **Route**: `/holistic-hypothesis-tester`
- **Features**: Compound validation, S/P/E scoring, evidence synthesis
- **API**: `/api/hypothesis/validate_food_dynamic`

**4. VUS Explorer** (`components/vus/`)
- **Route**: `/mutation-explorer`
- **Features**: Variant insights, coverage chips, WIWFM integration
- **API**: `/api/insights/predict_*`, `/api/efficacy/predict`

---

## 🔧 KEY SERVICES & PATTERNS

### **1. Sporadic Cancer Strategy**

**Problem**: 85-90% of cancer patients are sporadic (not germline-positive), but most platforms only work for germline-positive patients.

**Solution**: Tumor-centric analysis with germline awareness

**Key Services:**
- **`tumor_quick_intake.py`**: Level 0/1/2 tumor context generation
- **`sporadic_gates.py`**: PARP penalty/rescue, IO boost, confidence capping
- **`SporadicContext.jsx`**: Global state provider for frontend

**Key Gates:**
- **PARP Penalty**: Germline-negative → 0.6x (unless HRD ≥42 → 1.0x RESCUE!)
- **IO Boost**: TMB ≥20 → 1.3x, MSI-H → 1.3x, both → 1.69x
- **Confidence Cap**: L0 → 0.4, L1 → 0.6, L2 → none

### **2. Food Validator (Universal Hypothesis Testing)**

**Problem**: Patients ask about supplements (NAC, Vitamin D, Omega-3). Need evidence-backed, biomarker-aware recommendations.

**Solution**: Dynamic compound extraction + S/P/E scoring + evidence synthesis

**Key Services:**
- **`dynamic_food_extraction.py`**: ChEMBL/PubChem/LLM extraction (110M+ compounds)
- **`food_spe_integration.py`**: S/P/E scoring for compounds (0.4×S + 0.3×P + 0.3×E)
- **`enhanced_evidence_service.py`**: LLM paper reading (Gemini/Anthropic)
- **`dietician_recommendations.py`**: Dosage, timing, safety

**Key Features:**
- **Dynamic extraction**: ANY compound (ChEMBL/PubChem)
- **Evidence synthesis**: PubMed + Diffbot full-text
- **Biomarker-aware**: HRD+, TMB, treatment line context
- **Treatment line intelligence**: L1 vs L3 → Different recommendations

### **3. Clinical Trials System**

**Problem**: Match patients to trials with biomarker intelligence and transparent eligibility reasoning.

**Solution**: Hybrid search (AstraDB semantic + Neo4j graph) + autonomous agent

**Key Services:**
- **`hybrid_trial_search.py`**: AstraDB semantic search → Neo4j graph optimization
- **`autonomous_trial_agent.py`**: AI-driven search (no manual query required)
- **`ayesha_trials.py`**: Hard filters + soft boosts + eligibility checklists

**Key Features:**
- **Hard filters**: Stage IV, first-line, recruiting, NYC metro (90-95% confidence)
- **Soft boosts**: 10 boosts + 3 penalties (transparent scoring)
- **Eligibility checklists**: Hard/soft split with confidence gates
- **Autonomous agent**: Extracts patient context, auto-generates queries

### **4. Resistance Playbook**

**Problem**: Cancer evolves resistance. Need to predict resistance mechanisms and prepare backup plans BEFORE they happen.

**Solution**: Pattern-based resistance detection + combo strategies + next-line switches

**Key Services:**
- **`resistance_playbook_service.py`**: 5 detection rules, 7 combos, 6 switches
- **`resistance_detection_service.py`**: 2-of-3 trigger logic (HRD drop, DNA repair drop, CA-125 rise)

**Key Features:**
- **5 Detection Rules**: HR restoration, ABCB1 upregulation, MAPK/PI3K activation, SLFN11 loss, TMB/MSI changes
- **7 Combo Strategies**: PARP + ATR/CHK1/WEE1, PARP + Bevacizumab, Pembrolizumab + PARP/VEGF
- **6 Next-Line Switches**: Ceralasertib (ATR), Prexasertib (CHK1), Adavosertib (WEE1), etc.

### **5. SAE Intelligence System**

**Problem**: Evo2 delta scores are black-box. Need explainable features that tell doctors **WHY** a drug will work.

**Solution**: Sparse Autoencoder (SAE) features + mechanism mapping + DNA repair capacity

**Key Services:**
- **`sae_feature_service.py`**: DNA repair capacity, 7D mechanism vector, hotspot detection
- **`mechanism_fit_ranker.py`**: Trial ranking by mechanism alignment (α=0.7 eligibility + β=0.3 mechanism)
- **`next_test_recommender.py`**: Prioritizes next tests (HRD → ctDNA → SLFN11)
- **`hint_tiles_service.py`**: Actionable hints (max 4, suggestive tone)
- **`mechanism_map_service.py`**: 7D mechanism visualization (DDR, MAPK, PI3K, VEGF, IO, Efflux, HER2)

**Key Features:**
- **DNA Repair Capacity**: Manager's formula (0.6×DDR + 0.2×ess + 0.2×exon)
- **7D Mechanism Vector**: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux
- **Mechanism Fit Ranking**: α=0.7 eligibility + β=0.3 mechanism (Manager's P4)
- **Resistance Detection**: 2-of-3 triggers (HRD drop, DNA repair drop, CA-125 rise)

---

## 🎯 KEY ARCHITECTURAL DOCTRINES

### **Doctrine 1: "Wet Noodle" Principle**
**Problem**: DNA sequence alone (1D) doesn't tell you 3D protein folding  
**Solution**: Multi-stage validation pipeline
1. **Evo2** (Zeta Oracle): Sequence grammar check ("Is this plausible DNA?")
2. **AlphaFold 3** (Structural Sieve): 3D structure validation ("Does it fold correctly?")
3. **Wet-lab** (Final Gauntlet): Experimental validation ("Does it actually work?")

**Key Insight**: Evo2 is powerful but not sufficient. Structure ≠ Sequence. Function ≠ Structure.

### **Doctrine 2: Ground Truth Supremacy**
**Problem**: AI can hallucinate, especially in low-data regimes  
**Solution**: Ayesha Care System (100% confidence)
- **CA-125**: Deterministic burden classification (GCIG guidelines)
- **SOC recommendation**: NCCN guideline-aligned (95-100% confidence)
- **Trial matching**: Hard filters + transparent soft boosts (90-95% confidence)
- **WIWFM (post-NGS)**: Evo2-powered S/P/E (70-85% confidence)

**Key Insight**: Use AI where it's validated. Use ground truth (guidelines/literature) where it's available.

### **Doctrine 3: Graceful Degradation**
**Problem**: External services (Evo2, Fusion, Literature) can fail or timeout  
**Solution**: Fallback chains + placeholder values
- **Sequence scoring**: Fusion → Evo2 → Massive Oracle → Placeholder
- **Evidence**: Literature → ClinVar → Insufficient tier
- **Provenance**: Always track fallback path

**Key Insight**: Never crash. Always provide best-available answer with transparency.

### **Doctrine 4: Provenance is Sacred**
**Problem**: Clinical decisions need audit trails  
**Solution**: Track everything
- **Run IDs**: UUID for every operation
- **Profile tracking**: Baseline/Richer S/Fusion modes
- **Method provenance**: Which model, which service, which fallback
- **Timestamps**: When, where, how

**Key Insight**: If you can't explain how you got the answer, don't give the answer.

### **Doctrine 5: Feature Flags Enable Evolution**
**Problem**: Different contexts need different capabilities  
**Solution**: Environment-based toggles
- **Research mode**: `ENABLE_RESEARCH_MODE=true` → relaxed gates
- **Clinical mode**: `OPERATIONAL_MODE=clinical` → strict validation
- **Spam prevention**: `EVO_USE_DELTA_ONLY=true` → bounded upstream calls
- **Service toggles**: `DISABLE_FUSION`, `DISABLE_LITERATURE`, `DISABLE_EVO2`

**Key Insight**: Build flexibility into the architecture. One platform, many modes.

### **Doctrine 6: Sporadic Cancer is 85-90% of Reality**
**Problem**: Most platforms focus on germline-positive (10-15% of patients)  
**Solution**: Tumor-centric analysis with germline awareness
- **TumorContext schema**: L0 (priors), L1 (partial), L2 (full NGS)
- **PARP rescue**: HRD ≥42 → full effect (even if germline negative!)
- **IO boost**: TMB/MSI signals → 1.3x multipliers
- **Confidence capping**: Data completeness → confidence ceiling

**Key Insight**: Don't assume germline. Build for the 85-90% who are sporadic.

---

## 📊 COMPLETE SYSTEM INVENTORY

### **Backend Services (30+ routers, 50+ services)**

**Core Orchestrators:**
- `efficacy_orchestrator/` - S/P/E framework (main orchestrator)
- `ayesha_orchestrator_v2.py` - Complete care v2 (unified endpoint)
- `metastasis_interception_service.py` - 8-step cascade weapon design

**Insights & Intelligence:**
- `insights/` - 4 chips (Functionality, Chromatin, Essentiality, Regulatory)
- `sae_feature_service.py` - DNA repair capacity, mechanism vectors
- `resistance_playbook_service.py` - Resistance detection + combos
- `ca125_intelligence.py` - CA-125 burden/forecast/resistance

**Evidence & Literature:**
- `evidence/` - Literature search, ClinVar lookup, badge computation
- `enhanced_evidence_service.py` - LLM paper reading (Gemini/Anthropic)

**Design & Generation:**
- `design.py` - CRISPR guide generation, spacer efficacy prediction
- `therapeutic_optimizer.py` - Iterative optimization loops

**Clinical Decision Support:**
- `ayesha_trials.py` - Trial matching with filters/boosts
- `ngs_fast_track.py` - Test recommendations (ctDNA, HRD, IHC)
- `resistance_detection_service.py` - Early resistance detection

**Data & Extraction:**
- `dynamic_food_extraction.py` - Compound extraction (ChEMBL/PubChem)
- `hybrid_trial_search.py` - Hybrid search (AstraDB + Neo4j)
- `autonomous_trial_agent.py` - AI-driven trial search

**Sporadic Cancer:**
- `tumor_quick_intake.py` - Level 0/1 tumor context generation
- `sporadic_gates.py` - PARP penalty/rescue, IO boost, confidence capping

### **Frontend Components (100+ components, 30+ pages)**

**Core Pages:**
- `MyelomaDigitalTwin.jsx` - Core MM demo
- `AyeshaTrialExplorer.jsx` - Ayesha care system
- `HypothesisValidator.jsx` - Food validator
- `MutationExplorer.jsx` - VUS explorer

**Key Components:**
- `EfficacyPanel.jsx` - Drug ranking display
- `TrialMatchCard.jsx` - Trial match display
- `CA125Tracker.jsx` - CA-125 monitoring
- `SOCRecommendationCard.jsx` - SOC recommendation
- `NextTestCard.jsx` - Next test recommender
- `HintTilesPanel.jsx` - Actionable hints
- `MechanismChips.jsx` - 7D mechanism visualization

**Context Providers:**
- `AuthContext.jsx` - Authentication
- `SporadicContext.jsx` - Sporadic cancer state
- `CoPilotProvider` - Conversational AI
- `ActivityProvider` - Activity logging

---

## 🚀 KEY ENDPOINTS & CONTRACTS

### **Efficacy Prediction**
- **Endpoint**: `POST /api/efficacy/predict`
- **Input**: `{model_id, mutations[], options, germline_status?, tumor_context?}`
- **Output**: `{drugs[*]: {efficacy_score, confidence, evidence_tier, badges, insights, rationale, citations, provenance}}`

### **Insights Bundle**
- **Functionality**: `POST /api/insights/predict_protein_functionality_change`
- **Essentiality**: `POST /api/insights/predict_gene_essentiality`
- **Regulatory**: `POST /api/insights/predict_splicing_regulatory`
- **Chromatin**: `POST /api/insights/predict_chromatin_accessibility`

### **Ayesha Care**
- **Trials**: `POST /api/ayesha/trials/search`
- **Complete Care**: `POST /api/ayesha/complete_care_v2`
- **CA-125**: Integrated in complete care response

### **Food Validator**
- **Dynamic**: `POST /api/hypothesis/validate_food_dynamic`
- **Input**: `{compound, disease, patient_context?}`
- **Output**: `{spe_score, spe_percentile, verdict, evidence, dosage, timing, safety}`

### **Clinical Trials**
- **Hybrid Search**: `POST /api/trials/search-optimized`
- **Autonomous Agent**: `POST /api/trials/agent/search`
- **Graph Search**: `POST /api/trials/graph/search`

---

## 🎯 CONFIDENCE TIERS & METRICS

### **Confidence Levels**
- **95-100%**: Guideline-aligned SOC recommendations (NCCN)
- **90-95%**: Deterministic trial matching (hard filters + eligibility)
- **70-85%**: Evo2-powered S/P/E predictions (post-NGS)
- **40-60%**: Quick Intake mode (disease priors only)

### **Evidence Tiers**
- **Supported**: Strong literature + ClinVar-Strong + pathway alignment
- **Consider**: Moderate evidence + pathway alignment
- **Insufficient**: Weak evidence or evidence timeout

### **Badges**
- **RCT**: Randomized controlled trial evidence
- **Guideline**: Clinical guideline recommendation
- **ClinVar-Strong**: Expert/practice review + pathogenic classification
- **PathwayAligned**: Pathway score ≥0.2

---

## 🔍 KEY PATTERNS & PRINCIPLES

### **1. Graceful Degradation**
- Placeholder values when services unavailable
- Retry logic with exponential backoff
- Circuit breakers for service health monitoring
- **Example**: Evidence timeout → `evidence_fallback = True`, tier set to "insufficient"

### **2. Provenance Tracking**
- UUID-based run IDs across all operations
- Profile management: Baseline/Richer S/Fusion modes
- Complete audit trails: timestamps, parameters, methods
- **Example**: `provenance: {run_id, profile, cache, flags, sequence_scoring, fallback}`

### **3. Feature Flags**
- Environment-based toggles: `DISABLE_FUSION`, `EVO_USE_DELTA_ONLY`, `CONFIDENCE_V2`
- Operational modes: `OPERATIONAL_MODE=research|clinical`
- Spam-safety controls: `EVO_SPAM_SAFE`, `EVO_MAX_MODELS`, `EVO_MAX_FLANKS`
- **Location**: `api/config.py`

### **4. Modular Architecture**
- Services: ~100-150 lines each, single responsibility
- Routers: Thin endpoints delegating to services
- Schemas: Pydantic models for validation
- **Example**: Efficacy orchestrator broken into: sequence_processor, drug_scorer, sporadic_gates

### **5. Sporadic Awareness**
- Germline status: "positive", "negative", "unknown"
- Tumor context: TMB, MSI, HRD, somatic mutations
- Integrated throughout: Co-Pilot, Trials, Efficacy, Food Validator
- **Context Provider**: `SporadicContext.jsx`

### **6. Test Coverage**
- Unit tests: Clear acceptance criteria
- Integration tests: E2E workflows
- Smoke tests: Bash + curl + jq validation
- **Example**: `test_ayesha_trials.py` - 19 tests, all passing

---

## 📚 CRITICAL FILE LOCATIONS

### **Backend Core**
- Main API: `oncology-coPilot/oncology-backend-minimal/api/main.py`
- Config: `oncology-coPilot/oncology-backend-minimal/api/config.py`
- Efficacy Router: `oncology-coPilot/oncology-backend-minimal/api/routers/efficacy.py` (shim, delegates to orchestrator)

### **Efficacy Orchestrator**
- Main: `api/services/efficacy_orchestrator/orchestrator.py`
- Sequence: `api/services/efficacy_orchestrator/sequence_processor.py`
- Drug Scoring: `api/services/efficacy_orchestrator/drug_scorer.py`
- Models: `api/services/efficacy_orchestrator/models.py`
- Sporadic Gates: `api/services/efficacy_orchestrator/sporadic_gates.py`

### **Ayesha Care**
- CA-125: `api/services/ca125_intelligence.py`
- Trials: `api/routers/ayesha_trials.py`
- Orchestrator: `api/routers/ayesha_orchestrator_v2.py`
- NGS Fast-Track: `api/services/ngs_fast_track.py`

### **Frontend Core**
- App: `oncology-coPilot/oncology-frontend/src/App.jsx`
- Co-Pilot Logic: `oncology-coPilot/oncology-frontend/src/components/CoPilot/CoPilotLogic.jsx`
- Sporadic Context: `oncology-coPilot/oncology-frontend/src/context/SporadicContext.jsx`

---

## 🎓 WHAT I'VE LEARNED

### **1. Evo2 is the Foundation, Not the Oracle**
- Evo2 provides **ONE signal** (Sequence) in a multi-modal system
- S/P/E framework combines sequence, pathway, and evidence for higher accuracy
- Never rely on single metrics (Doctrine: Multi-Modal Validation)

### **2. Architecture is Modular & Extensible**
- Services are small (~100-150 lines), single responsibility
- Routers are thin endpoints delegating to services
- Feature flags enable different operational profiles
- Easy to add new capabilities without breaking existing ones

### **3. Sporadic Cancer is the Default**
- 85-90% of patients are sporadic (not germline-positive)
- PARP rescue (HRD ≥42) is a key differentiator
- IO boost (TMB/MSI) unlocks immunotherapy for more patients
- Confidence capping prevents overconfidence with incomplete data

### **4. Provenance Enables Trust**
- Every operation has a run ID
- Complete audit trails (methods, flags, fallbacks)
- Transparent confidence scores (not black boxes)
- Reproducible results (same inputs → same outputs)

### **5. Graceful Degradation is Mandatory**
- External services fail → Design for it
- Fallback chains: Fusion → Evo2 → Massive Oracle → Placeholder
- Evidence timeout → Insufficient tier (not crash)
- Always provide best-available answer with transparency

---

## 🚀 NEXT STEPS FOR CONFIDENCE

### **Areas to Deep Dive (Future Learning)**

1. **Modal Services**: How Evo2/Boltz services are deployed and scaled
2. **Database Schemas**: Supabase models, Neo4j relationships, AstraDB structure
3. **Frontend State Management**: How contexts interact, caching strategies
4. **Test Suites**: Expected behavior patterns, integration test workflows
5. **Documentation/Doctrine**: Strategic context, agent missions, completion reports

### **Confidence Level**

**Current**: ✅ **FOUNDATIONAL UNDERSTANDING COMPLETE**
- ✅ Core architecture understood
- ✅ S/P/E framework fully understood
- ✅ Evo2 integration patterns understood
- ✅ Frontend-backend integration understood
- ✅ Key services and patterns understood

**Ready For**: Building new features, extending existing capabilities, debugging issues, making architectural decisions

---

**⚔️ DOCTRINE STATUS: ACTIVE - FOUNDATIONAL LEARNING COMPLETE** ⚔️  
**This document represents Zo's complete understanding of the CrisPRO platform architecture and patterns.**

