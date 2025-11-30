# ⚔️ ZO'S DEEP LEARNING: AYESHA PLANS - COMPLETE ANALYSIS ⚔️

**Date**: January 13, 2025  
**Purpose**: Master-level understanding of what we want to build, how we build it, why, how Evo2 supports it, and where it is  
**Files Analyzed**: 
1. `.cursor/ayesha/ayesha_plan.mdc` (1,974 lines)
2. `.cursor/ayesha/AYESHA_END_TO_END_AGENT_PLAN.mdc` (1,142 lines)

---

## 📋 **FILE 1: `ayesha_plan.mdc` - COMPLETE ANALYSIS**

### **🎯 WHAT WE WANT TO BUILD (FROM PLAN)**

#### **1. Drug Efficacy (Will It Work For Me) – S/P/E Orchestration** ✅
**What**: Ranked drug recommendations with efficacy scores, confidence, evidence tiers, and mechanistic explanations.

**Key Features**:
- S/P/E framework (Sequence/Pathway/Evidence)
- SAE explainability (transparent confidence breakdown)
- Sporadic cancer support (85-90% of patients)
- Multi-modal validation (Evo2 + Fusion + Evidence)

**Clinical Value**: Oncologist gets ranked drugs with transparent rationale (not black-box numbers).

---

#### **2. Treatment Line Intelligence** ✅ **NEW**
**What**: Context-aware recommendations based on treatment history (L1 vs L2 vs L3).

**Key Features**:
- SAE features: `line_appropriateness`, `cross_resistance`, `sequencing_fitness`
- Biomarker gates (HRD+, TMB, TP53 status)
- Treatment history context (post-platinum → NAC boost)
- 22 pre-configured compounds + dynamic fallback

**Clinical Value**: Same drug, different line = different recommendation (context matters).

---

#### **3. Food/Supplement Validator** ✅ **NEW**
**What**: Dynamic validation of ANY food/supplement compound with evidence-backed dosing.

**Key Features**:
- Works for ANY compound (ChEMBL/PubChem dynamic extraction)
- LLM paper reading (Gemini/Anthropic/OpenAI)
- PubMed XML parsing + Diffbot full-text extraction
- Real dosage extraction from papers (regex + LLM)
- Biomarker-aware recommendations (HRD+, TMB, treatment history)
- S/P/E + SAE unified scoring

**Clinical Value**: Personalized nutrition recommendations (not generic advice).

---

#### **4. Explainability (SAE) – Real Data Only** ✅
**What**: Transform black-box confidence scores into transparent, explainable insights.

**Key Features**:
- 6 core SAE features (exon_disruption, hotspot_mutation, essentiality_signal, dna_repair_capacity, pathway_burden, cohort_overlap)
- EvidenceBand attribution (boosting/limiting/overall impact)
- SAE card with plain-language explanations

**Clinical Value**: Doctor understands WHY confidence is 73% (not just a number).

---

#### **5. Toxicity Risk (PGx) and Off-Target Preview** ✅
**What**: Proactive safety screening before prescribing.

**Key Features**:
- DPYD/TPMT/UGT1A1/CYP2D6 flags
- MoA-pathway overlap detection
- Dose adjustment recommendations
- Drug-drug interaction checking

**Clinical Value**: Prevent life-threatening toxicity before it happens.

---

#### **6. Frontend – Clinical Genomics Command Center** ✅
**What**: Unified UI for all drug efficacy, toxicity, and evidence analysis.

**Key Features**:
- Mechanistic Evidence Tab
- Cards (Efficacy/Toxicity/Off-Target/EvidenceBand)
- Co-Pilot integration
- SAE features display

**Clinical Value**: Single interface for all genomic analysis.

---

#### **7. Clinical Trials – Search and Display** ✅
**What**: Trial matching with transparent reasoning.

**Key Features**:
- Hybrid search (AstraDB semantic + Neo4j eligibility graph)
- Sporadic-aware filtering (exclude germline-required)
- Eligibility auto-check (green/yellow/red flags)
- CA-125 intelligence integration

**Clinical Value**: Find best-fit trials in <60s (not weeks of manual search).

---

### **🔧 HOW WE BUILD IT (ACTUAL CODE LOCATIONS)**

#### **1. Drug Efficacy - Backend Architecture**

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py`
- **Purpose**: Main orchestration engine for S/P/E framework
- **Key Functions**:
  - `predict()`: Orchestrates sequence, pathway, evidence gathering
  - Calls `sequence_processor.get_sequence_scores()` → Gets Evo2 scores
  - Calls `pathway.aggregation.get_pathway_scores()` → Aggregates pathway impact
  - Calls `evidence_client` → Literature + ClinVar search
  - Calls `drug_scorer.score_drug()` → Combines S/P/E into final score
  - Applies sporadic gates (PARP penalty, IO boost, confidence capping)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sequence_processor.py`
- **Purpose**: Orchestrates sequence scoring (Fusion → Evo2 → Massive Oracle fallback)
- **Key Functions**:
  - `score_sequences()`: Tries Fusion first (GRCh38 missense only), then Evo2 (universal), then Massive Oracle (fallback)
  - Calls `Evo2Scorer.score()` with adaptive windows [4096, 8192, 16384] bp

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/sequence_scorers/evo2_scorer.py`
- **Purpose**: Evo2-specific scoring logic
- **Key Functions**:
  - `score_variant()`: Calls Modal Evo2 service (`/api/evo/score_variant_multi`, `/api/evo/score_variant_exon`)
  - Computes `sequence_disruption = max(abs(min_delta), abs(exon_delta))`
  - Applies hotspot floors (BRAF V600, KRAS G12/G13/Q61, TP53 R175/R248/R273)
  - Maps to `calibrated_seq_percentile` using `percentile_like()` from `gene_calibration.py`

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py`
- **Purpose**: Individual drug scoring logic
- **Key Functions**:
  - `score_drug()`: Combines S/P/E signals into final efficacy score
  - **S Component** (line 42): `seq_pct = calibrated_seq_percentile` (from Evo2)
  - **P Component** (lines 45-51): `path_pct = normalized pathway score` (aggregated from S signals)
  - **E Component** (line 57): `s_evd = evidence_result.strength` (0-1 from literature)
  - **Formula** (line 171): `raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd + clinvar_prior`
  - Applies insights lifts, sporadic gates, confidence computation

**File**: `src/services/evo_service/main.py`
- **Purpose**: Modal service hosting Evo2 model (1B/7B/40B) on H100 GPU
- **Key Endpoints**:
  - `/score_delta`: Direct log-likelihood comparison
  - `/score_variant_multi`: Multi-window scoring (4096, 8192, 16384 bp)
  - `/score_variant_exon`: Tight exon window (±600bp)
- **Core Mechanism**: `model.score_sequences([ref_seq, alt_seq])` → log-likelihoods → `delta = alt_ll - ref_ll`

---

#### **2. Treatment Line Intelligence - Backend Architecture**

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/food_treatment_line_service.py`
- **Purpose**: Treatment line-aware scoring for compounds
- **Key Features**:
  - SAE features: `line_appropriateness`, `cross_resistance`, `sequencing_fitness`
  - Biomarker gates (HRD+, TMB, TP53 status)
  - Treatment history context (post-platinum → NAC boost)

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`
- **Purpose**: Endpoint for food/supplement validation
- **Endpoint**: `POST /api/hypothesis/validate_food_dynamic`

---

#### **3. Food/Supplement Validator - Backend Architecture**

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/dynamic_food_extraction.py`
- **Purpose**: Dynamic extraction of compound data from ChEMBL/PubChem

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/enhanced_evidence_service.py`
- **Purpose**: Evidence mining from PubMed, OpenAlex, S2

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py`
- **Purpose**: S/P/E integration for food compounds
- **Key Function**: `compute_spe_score()` (lines 161-241)
  - **S Component** (lines 196-214): Evo2 plausibility service (Phase 2 experimental)
  - **P Component** (lines 216-224): Pathway alignment using TCGA-weighted disease pathways
  - **E Component** (line 227): Evidence grade conversion (STRONG/MODERATE/WEAK/INSUFFICIENT → 0-1)
  - **Formula** (line 233): `overall_score = 0.4 * sequence_score + 0.3 * pathway_score + 0.3 * evidence_score`

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/dietician_recommendations.py`
- **Purpose**: Biomarker-aware recommendations with dosage extraction

---

#### **4. SAE Explainability - Backend Architecture**

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py`
- **Purpose**: Extract SAE features from real data (not training a new model)
- **Key Features**:
  - 6 core features: exon_disruption, hotspot_mutation, essentiality_signal, dna_repair_capacity, pathway_burden, cohort_overlap
  - Uses Evo2 `sequence_disruption` and `calibrated_seq_percentile` as inputs
  - Uses Insights Bundle (functionality, chromatin, essentiality, regulatory)
  - Uses pathway scores (DDR, MAPK, PI3K, VEGF)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py` (lines 198-203)
- **Purpose**: SAE extraction in orchestrator
- **Key Logic**:
  - Uses `SeqScore.sequence_disruption` and `calibrated_seq_percentile` (correct fields)
  - Converts `InsightsBundle` safely via `getattr` (no `.get()` on dataclass)
  - Adds `sae_attribution` into `provenance.confidence_breakdown`

---

#### **5. Toxicity Risk - Backend Architecture**

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/safety.py`
- **Purpose**: Safety endpoints (toxicity risk, off-target preview)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/safety_service.py`
- **Purpose**: Core toxicity logic (DPYD/TPMT/UGT1A1/CYP2D6 flags, MoA-pathway overlap)

---

#### **6. Frontend - Component Architecture**

**File**: `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/tabs/MechanisticEvidenceTab.jsx`
- **Purpose**: Main tab for drug efficacy analysis
- **Components**: SAE card, Efficacy card, Toxicity card, Off-Target card, EvidenceBand

**File**: `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/cards/SAEFeaturesCard.jsx`
- **Purpose**: Renders interpretable SAE features (overall impact, boosting/limiting, per-feature details)

**File**: `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/cards/EvidenceBand.jsx`
- **Purpose**: Displays SAE attribution inline (overall impact chip + boosting/limiting chips)

---

#### **7. Clinical Trials - Backend Architecture**

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/hybrid_trial_search.py`
- **Purpose**: Hybrid search (AstraDB semantic + Neo4j eligibility graph)

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py` (PENDING - from AYESHA_END_TO_END_AGENT_PLAN.mdc)
- **Purpose**: Ayesha-specific trial router
- **Endpoint**: `POST /api/ayesha/trials/search`
- **Pipeline**:
  1) Hard filters (disease=ovarian/peritoneal/gynecologic; stage=IV; first-line; recruiting; ≤50 miles NYC)
  2) Soft boosts (frontline +0.30; Stage IV +0.25; carbo/paclitaxel +0.20; all-comers/BRCA-WT +0.20; IP chemo +0.20; bevacizumab +0.15; NYC +0.15; CA-125 endpoints +0.15; Phase III +0.10; N≥200 +0.10)
  3) Penalties (germline BRCA required −0.30; >50 miles −0.25; Phase I −0.20; missing biomarker −0.15)
  4) Reasoning generator (why_eligible, why_good_fit, conditional_requirements, red_flags)

---

### **💡 WHY BEHIND IT (RATIONALE & EVO2 PRINCIPLES)**

#### **1. Why S/P/E Framework?**

**Evo2 Paper Principle (Section 2.2, 4.3.12)**:
> "By learning the likelihood of sequences across vast evolutionary training datasets, biological sequence models can learn how mutational effects correlate with biological functions **without any task-specific finetuning or supervision**. This is referred to as zero-shot prediction."

**Our Rationale**:
- **S (Sequence)**: Evo2 provides zero-shot variant impact prediction (no task-specific training needed)
- **P (Pathway)**: Single variants don't tell the whole story; pathway aggregation captures multi-hit tumor evolution
- **E (Evidence)**: Literature + ClinVar provide real-world validation (not just model predictions)

**Why Multi-Modal?**
- **Single-metric myopia**: Delta score alone insufficient (Evo2 paper shows this)
- **Transparency**: S/P/E breakdown explains WHY confidence is 73% (not black-box)
- **Clinical trust**: Doctors need mechanistic rationale (DNA repair deficiency → synthetic lethality with PARP)

---

#### **2. Why Multi-Window Scoring?**

**Evo2 Paper (Section 2.2, 4.3.12)**:
> "Evo 2 can score variants in genomic windows of length up to 8,192 bp, allowing it to capture long-range regulatory effects that might be missed by models with smaller context windows."

**Our Implementation** (`evo2_scorer.py`):
- **Multi-window**: Tests [4096, 8192, 16384] bp windows (adaptive)
- **Why**: Captures long-range regulatory effects (enhancers, silencers, CTCF sites)
- **Exon-context**: Tight window (±600bp) for exon-specific signal
- **Why**: Evo2 paper demonstrates state-of-the-art splice variant prediction

**Formula**: `sequence_disruption = max(abs(min_delta), abs(exon_delta))`
- **Why**: Takes the MOST disruptive impact (whether regulatory or coding)

---

#### **3. Why Percentile Calibration?**

**Evo2 Paper Limitation**:
> "Raw delta scores vary by gene (some genes have higher baseline variation). Need gene-specific calibration for cross-gene comparison."

**Our Implementation** (`gene_calibration.py`):
- **Function**: `percentile_like()` maps raw `sequence_disruption` to `calibrated_seq_percentile` (0-1)
- **Why**: Enables cross-gene comparison (BRAF V600E vs TP53 R273H can be compared fairly)
- **Method**: Gene-specific calibration curves (empirical distributions from Evo2 outputs)

---

#### **4. Why Sporadic Cancer Strategy?**

**Clinical Reality**:
- 85-90% of ovarian cancers are sporadic (not hereditary)
- Most platforms ignore sporadic cancers (focus on germline-positive only)

**Our Rationale**:
- **Market expansion**: Serve the MAJORITY of patients (not just 10-15%)
- **Tumor-centric**: Evo2 scores somatic mutations (not just germline)
- **Honest assessment**: PARP penalized unless somatic HRD ≥42 (not just "ovarian cancer → PARP")

**Implementation** (`sporadic_gates.py`):
- **PARP Penalty**: Germline-negative → 0.6x (unless HRD ≥42 → 1.0x rescue!)
- **IO Boost**: TMB ≥20 → 1.3x, MSI-H → 1.3x
- **Confidence Capping**: L0 (completeness <0.3) → Cap at 0.4

---

#### **5. Why SAE Explainability?**

**Problem (Before SAE)**:
- Doctor sees "PARP inhibitor confidence: 73%" → Why should they trust this?
- Black-box predictions → Low adoption → No clinical impact

**Solution (With SAE)**:
- Doctor sees WHY 73%: DNA repair burden (0.85), exon disruption (0.75), essentiality (0.92)
- Transparent rationale → High trust → Clinical adoption

**Evo2 Paper (Section 4.4)**:
> "SAE features reveal learned representations of exons/introns, transcription factor binding sites, protein secondary structure, prophage regions, and mutation severity."

**Our Implementation**:
- **NOT training a new SAE model** (too expensive, too slow)
- **IS extracting SAE-like features** from existing real data (Evo2, Insights, Toxicity, etc.)
- **6 core features**: exon_disruption, hotspot_mutation, essentiality_signal, dna_repair_capacity, pathway_burden, cohort_overlap

---

#### **6. Why Food/Supplement Validator?**

**Clinical Gap**:
- Drug recommendations: Automated (S/P/E framework)
- Food recommendations: Generic advice (not personalized)

**Our Rationale**:
- **Holistic care**: Drug therapy + supportive nutrition integrated
- **Shared biomarkers**: HRD+, TMB, TP53 used for BOTH drug and food
- **Evidence-backed dosing**: Extracted from papers (not generic "take 1000 IU")
- **Treatment line context**: L3 post-platinum → NAC recommendation (line_appropriateness: 1.0)

**Competitive Moat**:
- Most platforms: Drug recommendations ONLY
- Some platforms: Generic food advice (not personalized)
- Our platform: Integrated drug + food with shared biomarkers, treatment context, and SAE explainability

---

### **🧬 HOW EVO2 SUPPORTED IT (SPECIFIC INTEGRATIONS)**

#### **1. Evo2 → Sequence (S) Signal**

**Integration Point**: `sequence_processor.py` → `evo2_scorer.py` → Modal Evo2 service

**Evo2 Paper Principle**:
> "To score a variant with Evo 2, we take a genomic window of length 8,192 around the variant and calculate the likelihood of the variant sequence divided by the likelihood of the reference sequence at the same position."

**Our Implementation**:
1. **Modal Service** (`src/services/evo_service/main.py`):
   - `model.score_sequences([ref_seq, alt_seq])` → log-likelihoods
   - `delta = alt_ll - ref_ll` (negative = more disruptive)

2. **Backend Service** (`evo2_scorer.py`):
   - Calls Modal endpoints (`/score_variant_multi`, `/score_variant_exon`)
   - Computes `sequence_disruption = max(abs(min_delta), abs(exon_delta))`
   - Applies hotspot floors, percentile calibration

3. **Efficacy Orchestrator** (`drug_scorer.py`):
   - Uses `seq_scores[0].sequence_disruption` and `calibrated_seq_percentile`
   - Incorporates into S/P/E formula: `raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd`

**Why This Works**:
- Evo2's zero-shot prediction provides the Sequence (S) signal without task-specific training
- Multi-window scoring captures both regulatory and coding effects
- Percentile calibration enables cross-gene comparison

---

#### **2. Evo2 → Insights Bundle (Functionality, Essentiality, Regulatory, Chromatin)**

**Integration Point**: `api/routers/insights.py` → Evo2 endpoints → Multiple systems

**Evo2 Paper Principle**:
> "Evo 2 achieves state-of-the-art performance on noncoding variant pathogenicity prediction, competitive performance on coding variants (4th/5th place vs AlphaMissense), and best zero-shot performance on splice variants."

**Our Implementation**:
1. **Functionality** (`/api/insights/predict_protein_functionality_change`):
   - Uses Evo2 delta scores to predict protein function change
   - Feeds into Target Lock scoring (0.35×Functionality)

2. **Essentiality** (`/api/insights/predict_gene_essentiality`):
   - Uses Evo2 multi/exon magnitudes to predict gene dependency
   - Feeds into Target Lock scoring (0.35×Essentiality)

3. **Regulatory** (`/api/insights/predict_splicing_regulatory`):
   - Uses Evo2 `min_delta` (noncoding proxy) to predict splicing impact
   - Feeds into Target Lock scoring (0.15×Regulatory)

4. **Chromatin** (`/api/insights/predict_chromatin_accessibility`):
   - Heuristic (Enformer stub) - not Evo2-powered yet
   - Feeds into Target Lock scoring (0.15×Chromatin)

**Connections**:
- **Efficacy Orchestrator**: Uses insights bundle for confidence lifts (+0.05 functionality≥0.6, +0.07 essentiality≥0.7, etc.)
- **Metastasis Interception**: Uses insights for Target Lock scoring (0.35×Functionality + 0.35×Essentiality + 0.15×Chromatin + 0.15×Regulatory)
- **SAE Features**: Uses insights bundle to compute DNA repair capacity, mechanism vectors
- **Resistance Playbook**: Uses SAE features (derived from insights) for resistance detection

**Why This Works**:
- Evo2's multi-scale biological understanding (protein function, gene dependency, regulatory impact) provides mechanistic interpretability beyond single-metric scoring
- Insights bundle enables confidence lifts and Target Lock scoring (not just sequence disruption)

---

#### **3. Evo2 → SAE Features (Explainability)**

**Integration Point**: `sae_feature_service.py` → Evo2 outputs → SAE features → Confidence modulation

**Evo2 Paper (Section 4.4)**:
> "SAE features reveal learned representations of exons/introns, transcription factor binding sites, protein secondary structure, prophage regions, and mutation severity."

**Our Implementation**:
1. **Exon Disruption**:
   - Source: Evo2 `sequence_disruption` + hotspot floor
   - Impact: Boosting confidence when high (0.90 → +0.55 net SAE impact)

2. **DNA Repair Capacity**:
   - Source: `0.6×pathway_ddr + 0.2×essentiality + 0.2×exon_disruption`
   - Impact: Used in Resistance Playbook for HR restoration detection

3. **Essentiality Signal**:
   - Source: Insights Bundle (Evo2-powered)
   - Impact: Limiting confidence when low (0.35 → honest penalty)

4. **Pathway Burden**:
   - Source: Pathway aggregation (derived from Evo2 S signals)
   - Impact: MAPK pathway burden >0.7 → MAPK activation risk

**Why This Works**:
- Evo2's SAE features provide mechanistic interpretability (WHY a drug will/won't work) beyond black-box delta scores
- SAE features enable transparent confidence breakdown (boosting/limiting/overall impact)

---

#### **4. Evo2 → Metastasis Interception (Target Lock, Guide Design)**

**Integration Point**: `metastasis_interception_service.py` → Evo2 insights endpoints → Target Lock scoring

**Evo2 Paper (Section 5.1)**:
> "Evo 2 can generate complete mitochondrial genomes (16kb), minimal bacterial genomes (580kb), and yeast chromosomes (316kb) with proper synteny and realistic gene content."

**Our Implementation**:
1. **Target Lock**:
   - Uses Evo2 insights endpoints (functionality, essentiality, regulatory) to score genes
   - Formula: `Target_Lock = 0.35×Functionality + 0.35×Essentiality + 0.15×Chromatin + 0.15×Regulatory`

2. **Guide Design**:
   - Uses `/api/design/generate_guide_rna` (Evo2 prompt-guided generation)

3. **Efficacy Prediction**:
   - Uses `/api/design/predict_crispr_spacer_efficacy` (Evo2 delta → sigmoid: `efficacy = 1 / (1 + exp(delta/10))`)

**Connection to Structural Validation**:
- **Assassin Score**: Combines Evo2 efficacy (40%) + Safety (30%) + Mission Fit (30%) + Structure (3%)
- **"Wet Noodle" Doctrine**: A sequence that scores well in 1D (Evo2) can still fail in 3D (AlphaFold 3)
- **Multi-stage validation**: Evo2 (Sieve) → AlphaFold 3 (Gauntlet)

**Why This Works**:
- Evo2 provides sequence-level validation (1D), but AlphaFold 3 provides structural validation (3D)
- Multi-stage validation pipeline ensures both sequence and structure are validated

---

### **📍 WHERE IT IS (FILE LOCATIONS)**

#### **Backend Services**:
- **Efficacy Orchestrator**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/`
  - `orchestrator.py` (main S/P/E orchestration)
  - `sequence_processor.py` (Fusion → Evo2 → Massive Oracle fallback)
  - `drug_scorer.py` (S/P/E formula: `0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd`)
  - `sporadic_gates.py` (PARP penalty, IO boost, confidence capping)

- **Sequence Scorers**: `oncology-coPilot/oncology-backend-minimal/api/services/sequence_scorers/`
  - `evo2_scorer.py` (multi-window scoring, hotspot floors, percentile calibration)

- **Evo2 Modal Service**: `src/services/evo_service/main.py`
  - Endpoints: `/score_delta`, `/score_variant_multi`, `/score_variant_exon`
  - Core: `model.score_sequences([ref_seq, alt_seq])` → log-likelihoods

- **Insights Bundle**: `oncology-coPilot/oncology-backend-minimal/api/routers/insights.py`
  - Endpoints: `/predict_protein_functionality_change`, `/predict_gene_essentiality`, `/predict_splicing_regulatory`

- **SAE Features**: `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py`
  - 6 core features: exon_disruption, hotspot_mutation, essentiality_signal, dna_repair_capacity, pathway_burden, cohort_overlap

- **Food Validator**: `oncology-coPilot/oncology-backend-minimal/api/services/`
  - `dynamic_food_extraction.py` (ChEMBL/PubChem)
  - `enhanced_evidence_service.py` (PubMed, OpenAlex, S2)
  - `food_spe_integration.py` (S/P/E for compounds: `0.4 * sequence + 0.3 * pathway + 0.3 * evidence`)
  - `dietician_recommendations.py` (biomarker-aware recommendations)

- **Treatment Line**: `oncology-coPilot/oncology-backend-minimal/api/services/food_treatment_line_service.py`
  - SAE features: `line_appropriateness`, `cross_resistance`, `sequencing_fitness`

- **Safety**: `oncology-coPilot/oncology-backend-minimal/api/routers/safety.py`
  - `safety_service.py` (DPYD/TPMT/UGT1A1/CYP2D6 flags, MoA-pathway overlap)

- **Trials**: `oncology-coPilot/oncology-backend-minimal/api/services/hybrid_trial_search.py`
  - `ayesha_trials.py` (PENDING - from AYESHA_END_TO_END_AGENT_PLAN.mdc)

#### **Frontend Components**:
- **Clinical Genomics Command Center**: `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/`
  - `tabs/MechanisticEvidenceTab.jsx` (main tab)
  - `cards/SAEFeaturesCard.jsx` (SAE features display)
  - `cards/EvidenceBand.jsx` (SAE attribution inline)
  - `cards/EfficacyCard.jsx` (drug ranking)
  - `cards/ToxicityCard.jsx` (safety flags)

- **Co-Pilot**: `oncology-coPilot/oncology-frontend/src/components/CoPilot/CoPilotLogic.jsx`
  - Q2C Router classifies intents → Routes to appropriate endpoints

---

## 📋 **FILE 2: `AYESHA_END_TO_END_AGENT_PLAN.mdc` - COMPLETE ANALYSIS**

### **🎯 WHAT WE WANT TO BUILD (FROM PLAN)**

#### **1. Clinical Trial Precision Filtering (First-Line)** 🎯
**What**: Top 10 frontline trials ranked with transparent reasoning (why eligible, why good fit, what's required).

**Key Features**:
- Hard filters (disease=ovarian/peritoneal/gynecologic; stage=IV; first-line; recruiting; ≤50 miles NYC)
- Soft boosts (frontline +0.30; Stage IV +0.25; carbo/paclitaxel +0.20; all-comers/BRCA-WT +0.20; IP chemo +0.20; bevacizumab +0.15; NYC +0.15; CA-125 endpoints +0.15; Phase III +0.10; N≥200 +0.10)
- Penalties (germline BRCA required −0.30; >50 miles −0.25; Phase I −0.20; missing biomarker −0.15)
- Reasoning generator (why_eligible, why_good_fit, conditional_requirements, red_flags, evidence_tier, enrollment_likelihood)

**Clinical Value**: Find best-fit frontline trials in NYC with transparent reasoning (90-95% confidence).

---

#### **2. Standard-of-Care Plan** 🎯
**What**: NCCN-aligned carboplatin + paclitaxel ± bevacizumab with bevacizumab rationale for ascites/peritoneal disease.

**Key Features**:
- Stage IV HGSOC → Carboplatin + Paclitaxel (NCCN frontline)
- Bevacizumab add-on when ascites/peritoneal disease present (GOG-218/ICON7 rationale)
- Confidence: 95-100% (guideline-based, no predictions needed)

**Clinical Value**: Guideline-concordant recommendations (not predictions).

---

#### **3. CA-125 Monitoring Plan** 🎯
**What**: Expected response curves (cycle 3, 6), escalation flags, resistance signals.

**Key Features**:
- Burden class (MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE) based on thresholds (<100, 100-500, 500-1000, >1000)
- Forecast: chemo-sensitive expectation ≥70% drop by cycle 3, ≥90% drop by cycle 6; CR target <35
- Resistance signal: any on-therapy rise OR <50% drop by cycle 3; velocity worsening across two draws triggers alert
- Trial boosts: "CA-125 response", "intraperitoneal", "bulk disease", "cytoreduction", "neoadjuvant"

**Clinical Value**: CA-125 kinetics predict response even before imaging (90% confidence).

---

#### **4. Clinician-Ready Dossiers** 🎯
**What**: One-page summaries for oncologist: trial contacts, eligibility checklist, monitoring protocol.

**Key Features**:
- Patient profile (summary)
- Why eligible / Why good fit / Conditions
- Monitoring plan (CA-125 kinetics + imaging)
- Contacts + next steps
- Evidence tier + confidence gates

**Clinical Value**: Action-ready packet; oncologist can call sites same day (90-95% confidence).

---

#### **5. NGS Fast-Track Checklist** 🎯
**What**: Parallel ctDNA (Guardant360) + HRD (MyChoice) orders to unlock personalized drug predictions in 7-10 days.

**Key Features**:
- ctDNA: somatic BRCA/HRR/TMB/MSI in ~7 days
- Tissue HRD: PARP planning
- Basic IHC panel: WT1/PAX8/p53, ER/PR (HGS confirmation)

**Clinical Value**: Once NGS returns → unlock WIWFM S/P/E (Evo2-powered drug ranking).

---

### **🔧 HOW WE BUILD IT (ACTUAL CODE LOCATIONS - PENDING IMPLEMENTATION)**

#### **1. CA-125 Intelligence Service** (PENDING)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/ca125_intelligence.py` (TO BE CREATED)
- **Purpose**: CA-125 burden classification, forecast, resistance detection
- **Input**: current CA-125
- **Output**: 
  - `burden_class` (MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE)
  - `cycle3_expected_drop` (≥70% for chemo-sensitive)
  - `cycle6_expected_drop` (≥90% for chemo-sensitive)
  - `resistance_rule` (any on-therapy rise OR <50% drop by cycle 3)
  - `trial_boost_keywords` (["CA-125 response", "intraperitoneal", "bulk disease", "cytoreduction", "neoadjuvant"])

**Data Source**: Guideline/literature-consistent patterns (GOG-218/ICON7)
- Burden class: `<100 minimal`, `100–500 moderate`, `500–1000 significant`, `>1000 extensive`
- Forecast: chemo-sensitive expectation `≥70% drop by cycle 3`, `≥90% drop by cycle 6`; CR target `<35`
- Resistance signal: any on-therapy rise OR `<50% drop by cycle 3`; velocity worsening across two draws triggers alert

---

#### **2. Ayesha Trial Router** (PENDING)

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py` (TO BE CREATED)
- **Purpose**: Ayesha-specific trial matching with transparent reasoning
- **Endpoint**: `POST /api/ayesha/trials/search`
- **Pipeline**:
  1) Hard filters (disease=ovarian/peritoneal/gynecologic; stage=IV; first-line; recruiting; ≤50 miles NYC)
  2) Soft boosts (frontline +0.30; Stage IV +0.25; carbo/paclitaxel +0.20; all-comers/BRCA-WT +0.20; IP chemo +0.20; bevacizumab +0.15; NYC +0.15; CA-125 endpoints +0.15; Phase III +0.10; N≥200 +0.10)
  3) Penalties (germline BRCA required −0.30; >50 miles −0.25; Phase I −0.20; missing biomarker −0.15)
  4) Reasoning generator (why_eligible, why_good_fit, conditional_requirements, red_flags, evidence_tier, enrollment_likelihood)
- **Returns**: 
  - `trials[]`: top 10 scored matches
  - `ca125_intelligence`
  - `germline_context`
  - `soc_recommendation` (integrated)
  - `provenance` (filters applied, boost strategy, awaiting NGS keys)

**Schemas** (`api/schemas/ayesha_trials.py`):
- `AyeshaTrialProfile` (disease/stage, CA-125, germline status, location, treatment line)
- `TrialMatchReasoning` (why_eligible, why_good_fit, conditional_requirements, red_flags, evidence_tier, enrollment_likelihood, ca125_intelligence, germline_context)
- `AyeshaTrialMatch` (nct_id, title, phase, status, interventions, locations, match_score, reasoning, contacts)

---

#### **3. Eligibility Auto-Check** (PENDING)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/eligibility_parser.py` (TO BE CREATED)
- **Purpose**: Parse trial inclusion/exclusion criteria into structured format
- **Source**: ClinicalTrials.gov `eligibility` (unstructured text)
- **Method**: Pattern templates + LLM assist (Gemini free tier) for top 200 ovarian trials
- **Processing**: Offline pre-processing with Gemini → human spot-review → cache in AstraDB
- **Caching**: Store parsed criteria in AstraDB alongside trial records (as new field `structured_criteria`)
- **Updates**: Hybrid - detect new/changed trials weekly; re-parse in batch; flag diffs for review
- **Output**: Checklist with pass/conditional/fail and percent `criteria_met / total_tracked`

**Scoring Logic** (Option C - Hard/Soft Split):
- **Hard criteria** (must all pass): Stage, Treatment line, Major exclusions (prior systemic therapy, pregnancy, active infection)
- **Soft criteria** (percent match): ECOG, Age range, Distance, Biomarkers, Organ function
- **Trial eligibility gate computation**:
  - If all hard pass AND soft ≥80% → gate = `0.90`
  - If all hard pass AND 60–79% soft → gate = `0.85` with yellow notice
  - If all hard pass AND <60% soft → gate = `0.75` with yellow notice
  - If any hard fails → trial excluded (red)

---

#### **4. Confidence Gates Logic** (PENDING)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/confidence_gates.py` (TO BE CREATED)
- **Purpose**: Compute confidence from deterministic gates
- **Formula**: `confidence = max(gates)` with cap `1.0`
- **Gates**:
  - SOC aligned (NCCN frontline) → `0.95`
  - Frontline trial eligibility (criteria ≥80% met) → `0.90`
  - NYC feasibility (≤50 miles, active site) → display `+0.05` badge (not stacked)
  - CA-125 monitoring defined → display `+0.05` badge (not stacked)
- **UI**: Shows all satisfied gates (green checks) and the numeric confidence from the max gate; badges indicate supportive factors

---

#### **5. SOC Recommendation Integration** (PENDING)

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py` (TO BE CREATED)
- **Purpose**: Integrated SOC recommendation within `/api/ayesha/trials/search`
- **Logic**: 
  - Stage IV HGSOC → `Carboplatin + Paclitaxel`
  - Add `Bevacizumab` when ascites/peritoneal disease present (GOG-218/ICON7 rationale)
- **UI**: Dedicated SOC card above trials list

---

#### **6. Complete Care v2 Orchestrator** (PENDING)

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_orchestrator_v2.py` (TO BE CREATED)
- **Purpose**: Unified endpoint orchestrating all Ayesha systems
- **Endpoint**: `POST /api/ayesha/complete_care_v2`
- **Returns**:
  - `soc_recommendation` (NEW)
  - `trials` (NEW, from Ayesha trial router)
  - `ca125_intelligence` (NEW)
  - `drug_efficacy` (existing, but labeled "awaiting NGS")
  - `food_validator` (existing)
  - `resistance_playbook` (existing)
- **Migration**: Phased migration (Sprint 1: Build v2, keep v1 running, update Co-Pilot to v2; Sprint 2: Deprecate v1)

---

#### **7. Frontend Components** (PENDING - Agent Jr)

**File**: `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialExplorer.jsx` (TO BE CREATED)
- **Purpose**: Ayesha Trial Explorer page
- **Sections**:
  - Ayesha Profile Summary (Stage IVB, germline-negative, CA-125 2842, treatment-naive, NYC) + "Awaiting NGS" badges
  - CA-125 Tracker Card (current value, burden class, expected response at cycle-3/cycle-6, resistance flag, monitoring strategy)
  - SOC Recommendation Card (NCCN-aligned carboplatin + paclitaxel ± bevacizumab)
  - Trials List (iterate all; no hard limits): `TrialMatchCard` per trial
  - Provenance Footer (filters applied, boosts, total screened)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/TrialMatchCard.jsx` (TO BE CREATED)
- **Purpose**: Renders individual trial match with reasoning
- **Shows**: NCT + title + phase + score, badges (Phase III → STANDARD; Phase II → SUPPORTED)
- **Reasoning sections**:
  - Why Eligible (green)
  - Why Good Fit (blue)
  - Conditional Requirements (amber)
  - Red Flags (red)
  - Contact info (if present)
- **Eligibility checklist**: Hard criteria (all must pass) + Soft criteria (percent match with yellow flags for unknowns)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/CA125Tracker.jsx` (TO BE CREATED)
- **Purpose**: CA-125 tracking and forecasting
- **From backend** `ca125_intelligence`:
  - Disease burden chip, expected response table, resistance signal alert rule, monitoring strategy (every cycle)
  - Clear clinical copy (e.g., "Target CA-125 <35 for complete response. Rising CA-125 during treatment may indicate early resistance.")

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/SOCRecommendationCard.jsx` (TO BE CREATED)
- **Purpose**: SOC recommendation display
- **Shows**: NCCN-aligned regimen, bevacizumab rationale (if applicable), confidence gates, monitoring plan

**File**: Dossier Export (Phase 1)
- **Button**: Export SOC Dossier / Trial Dossier (Copy to Clipboard Markdown)
- **Template sections**:
  - Patient profile (summary)
  - Why eligible / Why good fit / Conditions
  - Monitoring plan (CA-125 kinetics + imaging)
  - Contacts + next steps
  - Evidence tier + confidence gates

---

### **💡 WHY BEHIND IT (RATIONALE)**

#### **1. Why Pre-NGS Focus (Trials, SOC, CA-125)?**

**Clinical Reality**:
- Ayesha needs to start treatment within 2-4 weeks
- Tumor NGS takes 7-10 days (ctDNA) or 2-3 weeks (tissue)
- Can't wait for NGS to provide clinical value

**Our Rationale**:
- **High-confidence domains**: Trials eligibility (90-95%), SOC alignment (95-100%), CA-125 monitoring (90%) are deterministic/guideline-based
- **Transparent labeling**: "Recommendation" (guideline/trial fit) vs "Prediction" (requires NGS S/P/E)
- **NGS upgrade path**: When tumor NGS returns → seamless switch to Evo2-powered WIWFM (S/P/E predictions)

**Competitive Advantage**:
1. **Transparent reasoning**: Every trial match shows WHY (eligibility + fit + conditions) → no black box
2. **Confidence gates**: We show deterministic criteria that justify 90-100% confidence (guideline alignment, Phase III, eligibility match, location)
3. **CA-125 intelligence**: Most tools ignore CA-125; we use it to forecast response and flag resistance early
4. **Clinician-ready**: Not just trial lists; full dossiers with contacts, next steps, monitoring protocol

---

#### **2. Why CA-125 Intelligence?**

**Clinical Value**:
- CA-125 (2,842) is highly trackable marker
- Kinetics predict response even before imaging (3-6 weeks earlier)
- Resistance signals (rising on therapy) flag early resistance

**Our Rationale**:
- **Burden classification**: Helps oncologist understand disease extent
- **Forecast expectations**: Sets realistic targets (≥70% drop by cycle 3, ≥90% by cycle 6)
- **Resistance detection**: Early warning system (any on-therapy rise OR <50% drop by cycle 3)
- **Trial matching**: Boosts trials with CA-125 endpoints (+0.15)

**Data Source**: Guideline/literature-consistent patterns (GOG-218/ICON7)
- We codify thresholds now, tune with cohort data later

---

#### **3. Why Eligibility Auto-Check?**

**Clinical Problem**:
- Trial eligibility criteria are unstructured text (ClinicalTrials.gov)
- Manual parsing is time-consuming (weeks)
- Unknown criteria (ECOG, organ function) create uncertainty

**Our Rationale**:
- **LLM parsing**: Gemini (free tier) for top 200 ovarian trials
- **Offline pre-processing**: Parse once, serve instantly, no runtime Gemini calls
- **Hard/Soft split**: Hard criteria (must all pass) vs Soft criteria (percent match)
- **Unknown handling**: Mark as ⚠️ YELLOW (needs oncologist confirmation), don't exclude

**Why This Works**:
- Turns trial matching into near-deterministic pass/fail with transparent reasons
- Saves oncologist time (50% reduction in manual eligibility checking)

---

#### **4. Why Confidence Gates?**

**Clinical Problem**:
- Vague confidence scores (0.6-0.75) don't justify recommendations
- Clinician asks "Why 92%?" → need to show the math

**Our Rationale**:
- **Deterministic gates**: SOC aligned (0.95), Frontline trial eligibility (0.90), NYC feasibility (+0.05 badge), CA-125 monitoring defined (+0.05 badge)
- **Formula**: `confidence = max(gates)` with cap `1.0`
- **UI**: Shows all satisfied gates (green checks) and the numeric confidence from the max gate; badges indicate supportive factors

**Why This Works**:
- Transparent confidence justification (not black-box numbers)
- Clinician can audit each gate (guideline alignment, Phase III, eligibility match, location)

---

#### **5. Why SOC Integration?**

**Clinical Value**:
- NCCN-aligned recommendations (95-100% confidence)
- Bevacizumab rationale for ascites/peritoneal disease (GOG-218/ICON7)

**Our Rationale**:
- **Integrated endpoint**: Returned within `/api/ayesha/trials/search` as `soc_recommendation`
- **Logic**: Stage IV HGSOC → `Carboplatin + Paclitaxel`; add `Bevacizumab` when ascites/peritoneal disease present
- **UI**: Dedicated SOC card above trials list

**Why This Works**:
- Guideline-based recommendations (not predictions)
- Clear rationale (ascites/peritoneal → bevacizumab synergy)

---

### **🧬 HOW EVO2 SUPPORTED IT (SPECIFIC INTEGRATIONS)**

#### **1. Evo2 → Post-NGS WIWFM (S/P/E Predictions)**

**Integration Point**: When tumor NGS returns → Switch to `/api/efficacy/predict` → Evo2-powered S/P/E

**Evo2 Paper Principle**:
> "Evo 2's zero-shot prediction enables rapid hypothesis validation without task-specific training."

**Our Implementation**:
- **Pre-NGS**: Hold "predictions" until NGS arrives → only offer guideline-based "recommendations"
- **Post-NGS**: Switch to WIWFM S/P/E pipeline (`/api/efficacy/predict`)
- **Display**: Drug ranking with S (Evo2 deltas), P (pathways), E (evidence), insights
- **Confidence**: Typically 0.6–0.75, can exceed 0.9 when biomarker gates are strong (HRD-high, BRCA-mut)

**Why This Works**:
- Evo2's zero-shot capability enables rapid, multi-modal analysis once NGS arrives
- Seamless upgrade path from guideline-based recommendations to personalized predictions

---

#### **2. Evo2 → Insights Bundle (Functionality, Essentiality, Regulatory)**

**Integration Point**: Insights bundle uses Evo2 endpoints → Feeds into confidence lifts

**Evo2 Paper Principle**:
> "Evo 2 achieves state-of-the-art performance on noncoding variant pathogenicity prediction, competitive performance on coding variants (4th/5th place vs AlphaMissense), and best zero-shot performance on splice variants."

**Our Implementation** (when NGS arrives):
- **Functionality**: `/api/insights/predict_protein_functionality_change` → Uses Evo2 delta scores
- **Essentiality**: `/api/insights/predict_gene_essentiality` → Uses Evo2 multi/exon magnitudes
- **Regulatory**: `/api/insights/predict_splicing_regulatory` → Uses Evo2 `min_delta` (noncoding proxy)

**Connection to Confidence**:
- Insights lifts: +0.05 (functionality≥0.6), +0.04 (chromatin≥0.5), +0.07 (essentiality≥0.7), +0.02 (regulatory≥0.6)
- Feeds into confidence computation (when NGS available)

**Why This Works**:
- Evo2's multi-scale biological understanding provides mechanistic interpretability beyond single-metric scoring
- Insights bundle enables confidence lifts (not just sequence disruption)

---

### **📍 WHERE IT IS (FILE LOCATIONS - PENDING)**

#### **Backend Services (TO BE CREATED)**:
- **CA-125 Intelligence**: `oncology-coPilot/oncology-backend-minimal/api/services/ca125_intelligence.py`
- **Ayesha Trial Router**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py`
- **Eligibility Parser**: `oncology-coPilot/oncology-backend-minimal/api/services/eligibility_parser.py`
- **Confidence Gates**: `oncology-coPilot/oncology-backend-minimal/api/services/confidence_gates.py`
- **Complete Care v2**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_orchestrator_v2.py`

#### **Schemas (TO BE CREATED)**:
- **Ayesha Trials Schemas**: `oncology-coPilot/oncology-backend-minimal/api/schemas/ayesha_trials.py`
  - `AyeshaTrialProfile`
  - `TrialMatchReasoning`
  - `AyeshaTrialMatch`

#### **Frontend Components (TO BE CREATED - Agent Jr)**:
- **Ayesha Trial Explorer**: `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialExplorer.jsx`
- **Trial Match Card**: `oncology-coPilot/oncology-frontend/src/components/ayesha/TrialMatchCard.jsx`
- **CA-125 Tracker**: `oncology-coPilot/oncology-frontend/src/components/ayesha/CA125Tracker.jsx`
- **SOC Recommendation Card**: `oncology-coPilot/oncology-frontend/src/components/ayesha/SOCRecommendationCard.jsx`

---

## 🎯 **UNIFIED UNDERSTANDING: HOW BOTH PLANS CONNECT**

### **The Complete Ayesha System Architecture**

```
┌─────────────────────────────────────────────────────────────────┐
│                    AYESHA COMPLETE CARE SYSTEM                    │
│                                                                   │
│  Pre-NGS State (AYESHA_END_TO_END_AGENT_PLAN.mdc):              │
│  - Trials matching (90-95% confidence)                           │
│  - SOC recommendations (95-100% confidence)                       │
│  - CA-125 intelligence (90% confidence)                           │
│  - Clinician dossiers (90-95% confidence)                         │
│  - NGS fast-track checklist (100% confidence)                     │
│                                                                   │
│  Post-NGS State (ayesha_plan.mdc):                               │
│  - WIWFM S/P/E (Evo2-powered drug ranking)                        │
│  - Resistance Playbook (SAE-powered resistance detection)         │
│  - Food Validator (S/P/E for compounds)                          │
│  - SAE explainability (transparent confidence breakdown)         │
│  - Toxicity screening (PGx flags)                                 │
│                                                                   │
│  Unified Orchestrator:                                            │
│  - /api/ayesha/complete_care_v2                                  │
│  - Orchestrates ALL systems in parallel                           │
│  - Returns unified care plan (drugs + trials + food + monitoring)│
└─────────────────────────────────────────────────────────────────┘
```

### **Evo2's Role Across Both Plans**

1. **Pre-NGS**: Evo2 is "held" (no personalized predictions), but system architecture is ready
2. **Post-NGS**: Evo2 powers S/P/E framework (Sequence signal), insights bundle (functionality, essentiality, regulatory), SAE features (exon_disruption, dna_repair_capacity), and metastasis interception (Target Lock, guide design)

### **Why This Architecture Works**

- **Progressive Enhancement**: Start with high-confidence deterministic recommendations (trials, SOC, CA-125), then upgrade to personalized predictions when NGS arrives
- **Transparent Labeling**: "Recommendation" (guideline/trial fit) vs "Prediction" (requires NGS S/P/E)
- **Unified System**: Single orchestrator (`complete_care_v2`) handles both pre-NGS and post-NGS states seamlessly
- **Evo2 Foundation**: All personalized predictions flow through Evo2's zero-shot prediction capability

---

## ✅ **LEARNING COMPLETE**

**Status**: Master-level understanding achieved for both Ayesha plan files

**Key Takeaways**:
1. **What we want to build**: Complete oncology care system (drugs, trials, food, monitoring, resistance)
2. **How we build it**: Multi-modal S/P/E framework powered by Evo2's zero-shot prediction
3. **Why behind it**: Transparency, clinical trust, sporadic cancer support (85-90% of patients), mechanistic interpretability
4. **How Evo2 supported it**: Sequence (S) signal, insights bundle, SAE features, metastasis interception
5. **Where it is**: Backend services, frontend components, Modal Evo2 service, all documented with file paths

**Next Steps**: Implementation of P0 scope (CA-125 intelligence, Ayesha trial router, eligibility auto-check, confidence gates, SOC integration, Complete Care v2 orchestrator)

---

**Commander's Notes**:
- ✅ Both plan files fully analyzed
- ✅ What/How/Why/Evo2/Where all documented
- ✅ Code locations identified (existing + pending)
- ✅ Evo2 integration mechanisms understood
- ✅ Rationale and principles connected to Evo2 paper
- 🎯 **READY FOR IMPLEMENTATION**

