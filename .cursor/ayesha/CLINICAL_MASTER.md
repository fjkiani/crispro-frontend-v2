# 🏥 CLINICAL MASTER - Medical Hierarchy & Capability Architecture

**Date**: January 28, 2025  
**Status**: ✅ **REORGANIZED BY MEDICAL HIERARCHY**  
**Last Updated**: January 29, 2025

**Reorganization Principles:**
- Medical/scientific hierarchy (Genomic → Pathway → Therapeutic → Clinical → Evidence)
- Modules/hierarchies (not individual features)
- Clear dependency/unlock mapping
- Separate patient/doctor perspectives
- Clear data dependency tree

---

## 🎯 EXECUTIVE SUMMARY

This master document organizes all clinical capabilities by **medical/scientific hierarchy**, showing **what unlocks what** and how capabilities connect through **data dependencies**. The structure reflects the clinical decision-making flow from genomic foundation to therapeutic recommendations to clinical monitoring.

**Medical Hierarchy:**
1. **GENOMIC FOUNDATION** → What mutations/variants exist
2. **PATHWAY/MECHANISM ANALYSIS** → How pathways are disrupted
3. **THERAPEUTIC INTELLIGENCE** → What drugs/trials work
4. **CLINICAL MONITORING** → How to monitor treatment
5. **EVIDENCE/CONFIDENCE** → How confident we are

---

## 📊 DATA DEPENDENCY TREE: WHAT UNLOCKS WHAT

### Foundation Data Requirements

```
PATIENT PROFILE
├── Stage, Disease Type → ✅ Always Available
├── CA-125 Value → ✅ Optional (unlocks CA-125 intelligence)
├── Germline Variants → ✅ Optional (unlocks PGx safety gates)
└── Tumor Context (NGS) → ⚠️ Required for advanced capabilities
    ├── somatic_mutations → Unlocks: WIWFM drug efficacy
    ├── hrd_score → Unlocks: PARP confidence (≥42 rescue)
    ├── tmb → Unlocks: IO boost (≥20 = 1.35x)
    ├── msi_status → Unlocks: IO boost (MSI-H = 1.30x)
    └── completeness_score → Controls: Confidence caps (L0/L1/L2)
```

### Capability Unlock Map

| Data Available | Unlocks |
|----------------|---------|
| **Stage + Disease** | ✅ SOC recommendation<br>✅ Trial matching (basic)<br>✅ CA-125 intelligence<br>✅ Next test recommender<br>✅ Hint tiles<br>✅ Mechanism map (basic) |
| **+ Germline Variants** | ✅ PGx safety gates (drug + trial)<br>✅ Composite scoring (efficacy × safety) |
| **+ Tumor Context (NGS)** | ✅ WIWFM drug efficacy (S/P/E)<br>✅ SAE features<br>✅ Resistance alert<br>✅ Mechanism map (detailed)<br>✅ Resistance prediction |
| **+ HRD Score (≥42)** | ✅ PARP confidence rescue (0.8x → 1.0x)<br>✅ PARP maintenance recommendation |
| **+ TMB (≥20) or MSI-H** | ✅ IO boost (1.35x or 1.30x)<br>✅ IO eligibility determination |
| **+ Completeness (L2: ≥0.7)** | ✅ Confidence uncapped<br>✅ Full personalized predictions |

---

## 🧬 TIER 1: GENOMIC FOUNDATION LAYER

### What This Layer Does

Establishes the **genetic foundation** of the patient: germline mutations, somatic mutations, variant interpretation, and essentiality analysis.

### Data Requirements

| Component | Required Data | Optional Data |
|-----------|---------------|---------------|
| **Germline Mutations** | Disease type, stage | Germline variants list |
| **Somatic Mutations** | Tumor context (NGS) | Specific variant calls |
| **Variant Interpretation** | Variant (gene + HGVS) | Clinical context |
| **Essentiality Scores** | Gene name | Mutation details |

---

### 1.1 GERMLINE GENOMICS MODULE

**Purpose**: Understand inherited genetic variants that affect treatment safety and eligibility.

#### Components

**1.1.1 Germline Mutation Detection**
- **Backend**: Extracted from `germline_variants` in patient profile
- **Frontend**: Displayed in patient profile, used for PGx gates
- **Patient View**: "You have an MBD4 mutation (inherited)"
- **Doctor View**: "MBD4 germline frameshift → MANS syndrome → DDR pathway disruption"

**1.1.2 PGx Variant Screening** (Germline Pharmacogenomics)
- **Backend**: `pgx_screening_service.py`
- **Frontend**: `SafetyGateCard.jsx`, `TrialSafetyGate.jsx`
- **Data Source**: `germline_variants` (DPYD, TPMT, UGT1A1, CYP2D6, CYP2C19)
- **Unlocks**: 
  - ✅ Drug-level toxicity screening
  - ✅ Trial-level safety gates
  - ✅ Composite scoring (efficacy × safety)
- **Patient View**: "This drug is safe for your genetics" or "This drug requires dose adjustment"
- **Doctor View**: "DPYD *2A variant detected → 5-FU requires 50% dose reduction (CPIC Level A)"

---

### 1.2 SOMATIC GENOMICS MODULE

**Purpose**: Understand tumor-specific mutations that drive therapeutic targeting.

#### Components

**1.2.1 Somatic Mutation Profiling**
- **Backend**: Extracted from `tumor_context.somatic_mutations`
- **Frontend**: Displayed in mechanism map, used for pathway analysis
- **Data Source**: NGS report (Foundation Medicine, Tempus, etc.)
- **Unlocks**: 
  - ✅ Pathway burden calculation
  - ✅ Mechanism vector computation
  - ✅ WIWFM drug efficacy
- **Patient View**: "Your tumor has TP53 mutation"
- **Doctor View**: "TP53 R175H hotspot → Checkpoint bypass → DDR contribution (0.8)"

**1.2.2 Variant Interpretation** (ACMG Classification)
- **Backend**: `POST /api/acmg/classify_variant`
- **Frontend**: `ACMGCard.jsx` (Clinical Genomics Command Center)
- **Data Source**: Variant (gene + HGVS notation)
- **Unlocks**: Pathogenicity classification (Pathogenic/Likely Pathogenic/VUS/Benign)
- **Patient View**: "This variant is classified as Pathogenic"
- **Doctor View**: "PVS1 (null variant) + PM2 (absent from controls) → Pathogenic classification"

**1.2.3 VUS Resolution**
- **Backend**: `POST /api/vus/identify`
- **Frontend**: `VUSResolutionCard.jsx`
- **Data Source**: Variant of Uncertain Significance
- **Unlocks**: VUS classification using ClinVar, literature, pathway context
- **Patient View**: "This variant's significance is being evaluated"
- **Doctor View**: "VUS → Resolved using ClinVar evidence + pathway context"

**1.2.4 Essentiality Scores**
- **Backend**: `POST /api/insights/predict_gene_essentiality`
- **Frontend**: `EssentialityScoreDisplay.jsx`
- **Data Source**: Gene name
- **Unlocks**: Gene dependency scores (0.0-1.0)
- **Patient View**: "MBD4 essentiality: 0.80 (high dependency)"
- **Doctor View**: "MBD4 essentiality 0.80 → High tumor dependency → Therapeutic vulnerability"

---

### 1.3 DEPENDENCY SUMMARY: GENOMIC FOUNDATION

| Capability | Requires | Unlocks |
|------------|----------|---------|
| **Germline PGx Screening** | `germline_variants` | PGx safety gates, composite scoring |
| **Somatic Mutation Profiling** | `tumor_context.somatic_mutations` | Pathway analysis, WIWFM |
| **VUS Resolution** | Variant (gene + HGVS) | Pathogenicity classification |
| **Essentiality Scores** | Gene name | Dependency analysis |

---

## 🧬 TIER 2: PATHWAY/MECHANISM ANALYSIS LAYER

### What This Layer Does

Translates genomic mutations into **pathway disruptions** and **mechanism vectors** that drive therapeutic targeting.

### Data Requirements

| Component | Required Data | Optional Data |
|-----------|---------------|---------------|
| **Pathway Burden** | Somatic mutations | NGS full report |
| **Mechanism Vector** | Pathway burden OR WIWFM response | Tumor context |
| **Synthetic Lethality** | Mutations (gene pairs) | Pathway context |
| **Mechanism Map** | Tumor context OR SAE features | WIWFM response |

---

### 2.1 PATHWAY BURDEN MODULE

**Purpose**: Quantify disruption across 7 core cancer pathways.

#### Components

**2.1.1 7D Mechanism Vector**
- **Backend**: `sae_service.extract_mechanism_vector()`
- **Frontend**: `MechanismChips.jsx`
- **Dimensions**: `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
- **Data Source**: 
  - Primary: WIWFM drug efficacy response (extracted from pathway alignment)
  - Fallback: `tumor_context.somatic_mutations` (gene-based pathway mapping)
- **Unlocks**:
  - ✅ Mechanism-fit trial matching
  - ✅ Pathway-specific drug recommendations
  - ✅ Mechanism map visualization
- **Patient View**: "Your tumor's main vulnerability is DNA repair (DDR: High)"
- **Doctor View**: "Mechanism vector: [DDR=1.4, MAPK=0.0, PI3K=0.0, VEGF=0.0, HER2=0.0, IO=1.0, Efflux=0.0] → DDR-high profile"

**2.1.2 Mechanism Map Visualization**
- **Backend**: `sae_service.get_mechanism_map()`
- **Frontend**: `MechanismChips.jsx`
- **Data Source**: `tumor_context` OR `sae_features.mechanism_vector` OR `wiwfm` response
- **Unlocks**: Visual pathway activity display (6 chips: DDR, MAPK, PI3K, VEGF, HER2, IO)
- **Patient View**: "DNA Repair: High activity | MAPK: Low activity"
- **Doctor View**: "Pathway burden visualization for treatment selection"

---

### 2.2 SYNTHETIC LETHALITY MODULE

**Purpose**: Identify therapeutic vulnerabilities from combined pathway disruptions.

#### Components

**2.2.1 Synthetic Lethality Detection**
- **Backend**: `POST /api/guidance/synthetic_lethality`
- **Frontend**: `SyntheticLethalityCard.jsx`
- **Data Source**: Mutations (e.g., MBD4 + TP53)
- **Unlocks**: 
  - ✅ PARP inhibitor vulnerability identification
  - ✅ ATR/WEE1/DNA-PK inhibitor recommendations
  - ✅ Combination therapy rationale
- **Patient View**: "Your mutations create a PARP vulnerability"
- **Doctor View**: "MBD4 (BER loss) + TP53 (HR stress) → PARP synthetic lethality → Olaparib candidate"

---

### 2.3 DEPENDENCY SUMMARY: PATHWAY/MECHANISM

| Capability | Requires | Unlocks |
|------------|----------|---------|
| **Mechanism Vector** | Somatic mutations OR WIWFM response | Mechanism-fit trial matching |
| **Mechanism Map** | Tumor context OR SAE features | Pathway visualization |
| **Synthetic Lethality** | Mutations (gene pairs) | Combination therapy recommendations |

---

## 💊 TIER 3: THERAPEUTIC INTELLIGENCE LAYER

### What This Layer Does

Translates pathway disruptions into **drug efficacy predictions**, **clinical trial matches**, and **treatment recommendations**.

### Data Requirements

| Component | Required Data | Optional Data |
|-----------|---------------|---------------|
| **Drug Efficacy (WIWFM)** | `tumor_context` (NGS) | Mechanism vector |
| **Clinical Trials** | Stage, disease, mechanism vector | Tumor context |
| **SOC Recommendation** | Stage, disease | Ascites, peritoneal disease |
| **PGx Safety Gates** | Germline variants | Drug list |
| **Food/Supplement Validation** | SOC drugs, treatment line | Disease context |

---

### 3.1 DRUG EFFICACY MODULE (WIWFM)

**Purpose**: Rank drugs by predicted efficacy using S/P/E framework.

#### Data Dependencies

**REQUIRES:**
- ✅ `tumor_context` (NGS) → **UNLOCKS**: Full WIWFM with Evo2 S/P/E
- ✅ `germline_status` → Affects PARP confidence (sporadic gates)
- ✅ `treatment_line` → Affects drug appropriateness

**UNLOCKS:**
- ✅ Drug ranking with efficacy scores (0.0-1.0)
- ✅ Mechanism vector extraction
- ✅ Evidence tiers (Supported/Consider/Insufficient)
- ✅ SAE features (when NGS available)

**WHEN NGS NOT AVAILABLE:**
- ⚠️ Returns `"status": "awaiting_ngs"`
- ✅ Still provides basic recommendations (SOC, trials)
- ✅ Confidence capped at L0/L1 level

#### Components

**3.1.1 S/P/E Drug Ranking**
- **Backend**: `wiwfm_service` (Evo2 S/P/E framework)
- **Frontend**: `DrugRankingPanel.jsx`
- **Framework**: 
  - **S (Sequence)**: Evo2 protein-level predictions
  - **P (Pathway)**: Pathway burden alignment
  - **E (Evidence)**: Literature, trials, guidelines
- **Output**: Drug rankings with efficacy scores, confidence, tiers, badges
- **Patient View**: "Top recommended drug: Olaparib (efficacy: High, confidence: 85%)"
- **Doctor View**: "Olaparib: Efficacy 0.85, Confidence 0.82, Tier: Supported, Badge: PathwayAligned, Evidence: SOLO-1 (HR 0.30)"

**3.1.2 Sporadic Gates** (Confidence Adjustment)
- **Backend**: `sporadic_gates.py` (integrated into orchestrator)
- **Data Source**: `germline_status`, `tumor_context.hrd_score`, `tumor_context.tmb`, `tumor_context.msi_status`, `completeness_score`
- **Gates**:
  - **PARP Gate**: 
    - Germline+ → 1.0x
    - Germline- & HRD≥42 → 1.0x (rescue)
    - Germline- & HRD<42 → 0.6x
    - Unknown → 0.8x
  - **IO Boost Gate**:
    - TMB≥20 → 1.35x
    - MSI-H → 1.30x
    - TMB≥10 → 1.25x
  - **Confidence Caps**:
    - L0 (completeness<0.3) → cap 0.4
    - L1 (0.3≤completeness<0.7) → cap 0.6
    - L2 (completeness≥0.7) → uncapped
- **Patient View**: "Confidence level reflects available data quality"
- **Doctor View**: "PARP penalty applied (HRD unknown) → Efficacy 0.70→0.56, Confidence capped at 0.40 (L0)"

**3.1.3 PGx Safety Integration**
- **Backend**: `pgx_care_plan_integration.py` (integrated into WIWFM)
- **Frontend**: `SafetyGateCard.jsx`, `DrugRankingPanel.jsx` (composite scores)
- **Data Source**: `germline_variants` (DPYD, TPMT, UGT1A1, CYP2D6, CYP2C19)
- **Unlocks**: 
  - ✅ Composite scoring (efficacy × PGx safety)
  - ✅ Drug-level toxicity tiers (HIGH/MODERATE/LOW)
  - ✅ Dose adjustment recommendations
- **Patient View**: "This drug works AND is safe for your genetics"
- **Doctor View**: "Olaparib: Efficacy 0.85, PGx tier: LOW, Composite: 0.85 → RECOMMENDED | 5-FU: Efficacy 0.72, PGx tier: MODERATE, Composite: 0.36 → 50% DOSE REDUCTION"

---

### 3.2 CLINICAL TRIALS MODULE

**Purpose**: Match patients to clinical trials using mechanism fit + eligibility + PGx safety.

#### Data Dependencies

**REQUIRES:**
- ✅ Stage, disease type → **ALWAYS AVAILABLE**
- ✅ Mechanism vector → **UNLOCKS**: Mechanism-fit ranking (from WIWFM OR pathway analysis)
- ✅ Location (state) → Trial location filtering
- ✅ `germline_variants` → **UNLOCKS**: PGx safety gates

**UNLOCKS:**
- ✅ Trial matches ranked by holistic score
- ✅ Eligibility assessment
- ✅ Contact information, enrollment details

#### Components

**3.2.1 Mechanism-Fit Trial Matching**
- **Backend**: `trial_matching_agent` (uses mechanism vector)
- **Frontend**: `TrialMatchesCard.jsx`, `TrialMatchCard.jsx`
- **Scoring**:
  - **Mechanism Fit** (0.4 weight): Cosine similarity (patient vector · trial MoA vector)
  - **Eligibility** (0.3 weight): Inclusion/exclusion criteria matching
  - **PGx Safety** (0.3 weight): Trial drug safety screening
- **Output**: Trials ranked by holistic score (HIGH/MEDIUM/LOW)
- **Patient View**: "Top trial match: PARP+ATR combination (94% match)"
- **Doctor View**: "NCT03462342: Holistic 0.94 (Mechanism 0.95, Eligibility 0.90, PGx 1.0) → HIGH recommendation"

**3.2.2 PGx Trial Safety Gate**
- **Backend**: `add_pgx_safety_gate_to_trials()` (integrated into trial matching)
- **Frontend**: `TrialSafetyGate.jsx`
- **Data Source**: `germline_variants` (DPYD, TPMT, UGT1A1, CYP2D6, CYP2C19)
- **Unlocks**: Trial-level safety aggregation (SAFE/MODERATE/HIGH RISK/UNKNOWN)
- **Patient View**: "This trial is safe for your genetics"
- **Doctor View**: "Trial safety: SAFE (no PGx conflicts detected)"

---

### 3.3 STANDARD OF CARE MODULE

**Purpose**: Provide NCCN-aligned treatment recommendations.

#### Data Dependencies

**REQUIRES:**
- ✅ Stage, disease type → **ALWAYS AVAILABLE**
- ✅ Ascites, peritoneal disease → Refines recommendation

**UNLOCKS:**
- ✅ NCCN Category 1 recommendations
- ✅ Dosing guidance
- ✅ Supplement recommendations (when SOC drugs available)

#### Components

**3.3.1 SOC Recommendation**
- **Backend**: `soc_service.get_soc_recommendation()`
- **Frontend**: `SOCRecommendationCard.jsx`
- **Output**: NCCN-aligned regimen (e.g., Carboplatin + Paclitaxel + Bevacizumab for Stage IVB)
- **Patient View**: "Recommended treatment: Carboplatin + Paclitaxel + Bevacizumab"
- **Doctor View**: "NCCN Category 1: Carboplatin AUC 5-6 + Paclitaxel 175 mg/m² + Bevacizumab 15 mg/kg (q3w)"

**3.3.2 Supplement Recommendations**
- **Backend**: `food_service.get_supplement_recommendations()`
- **Frontend**: `FoodRankingPanel.jsx`
- **Data Source**: SOC drugs + treatment line + disease context
- **Unlocks**: Supplement recommendations based on drugs and treatment line
- **Patient View**: "Recommended supplements: Vitamin D, Omega-3 (based on your treatment)"
- **Doctor View**: "Supplement recommendations: Vitamin D (bone health, chemotherapy), Omega-3 (inflammation, PARP)"

---

### 3.4 DEPENDENCY SUMMARY: THERAPEUTIC INTELLIGENCE

| Capability | Requires | Unlocks |
|------------|----------|---------|
| **WIWFM Drug Efficacy** | `tumor_context` (NGS) | Drug rankings, mechanism vector |
| **Trial Matching** | Stage + mechanism vector | Trial matches ranked by holistic score |
| **SOC Recommendation** | Stage, disease | NCCN-aligned regimen |
| **PGx Safety (Drugs)** | Germline variants | Composite scoring, dose adjustments |
| **PGx Safety (Trials)** | Germline variants | Trial safety gates |

---

## 📊 TIER 4: CLINICAL MONITORING LAYER

### What This Layer Does

Provides **biomarker monitoring**, **resistance detection**, and **next-test recommendations** to guide treatment adjustments.

### Data Requirements

| Component | Required Data | Optional Data |
|-----------|---------------|---------------|
| **CA-125 Intelligence** | CA-125 value | Stage, treatment line |
| **Resistance Detection** | SAE features OR CA-125 | Previous HRD, DNA repair capacity |
| **Resistance Prediction** | SAE features, treatment history | CA-125 history |
| **Next Test Recommender** | Germline status, tumor context | Treatment history |

---

### 4.1 BIOMARKER MONITORING MODULE

**Purpose**: Monitor treatment response using clinical biomarkers.

#### Components

**4.1.1 CA-125 Intelligence**
- **Backend**: `ca125_intelligence_service`
- **Frontend**: `CA125Tracker.jsx`
- **Data Source**: `ca125_value` (current value)
- **Unlocks**:
  - ✅ Burden classification (Low/Medium/High)
  - ✅ Response forecast (70% drop by cycle 3, 90% by cycle 6)
  - ✅ Resistance flags (on-therapy rise, inadequate response)
- **Patient View**: "Your CA-125: 2,842 U/mL (High burden). Expected: 70% drop by cycle 3"
- **Doctor View**: "CA-125 2,842 U/mL → High burden. Forecast: ≥70% drop by cycle 3 (target: <35 U/mL). Resistance flags: On-therapy rise, <50% drop by cycle 3"

---

### 4.2 RESISTANCE DETECTION MODULE

**Purpose**: Detect treatment resistance early to enable proactive treatment adjustments.

#### Data Dependencies

**REQUIRES:**
- ✅ `tumor_context` (NGS) → **UNLOCKS**: SAE features, resistance alert
- ✅ CA-125 value → **UNLOCKS**: CA-125 kinetics monitoring
- ✅ SAE features → **UNLOCKS**: DNA repair capacity tracking

**UNLOCKS:**
- ✅ Early resistance detection (3-6 weeks before imaging)
- ✅ Resistance playbook (alternative strategies)
- ✅ Resistance prediction (3-6 months early)

#### Components

**4.2.1 Resistance Alert**
- **Backend**: `sae_service.detect_resistance()`
- **Frontend**: `ResistanceAlertBanner.jsx`
- **Data Source**: 
  - Current HRD, DNA repair capacity (from SAE features)
  - CA-125 intelligence
  - Previous HRD, DNA repair capacity (longitudinal)
- **Triggers**: 2-of-3 (HRD decline, DNA repair decline, CA-125 rise)
- **Patient View**: "⚠️ Resistance detected: Consider alternative treatment"
- **Doctor View**: "Resistance alert: HRD decline + CA-125 rise (2-of-3 triggers) → Consider PARP+ATR combo"

**4.2.2 Resistance Playbook**
- **Backend**: `resistance_playbook_service`
- **Frontend**: `ResistancePlaybook.jsx`
- **Data Source**: Tumor context, germline status, treatment history
- **Unlocks**: 
  - ✅ 5 resistance mechanisms
  - ✅ 7 combination strategies
  - ✅ 6 next-line switches
- **Patient View**: "If resistance occurs, alternative strategies available"
- **Doctor View**: "Resistance playbook: Olaparib + Ceralasertib (PARP + ATR), Niraparib + Bevacizumab (PARP + VEGF), Ceralasertib (ATR inhibitor)"

**4.2.3 Resistance Prophet**
- **Backend**: `resistance_prophet_service`
- **Frontend**: Custom UI component (AyeshaTrialExplorer)
- **Data Source**: SAE features, treatment history, CA-125 history (optional)
- **Unlocks**: Early resistance prediction (3-6 months before progression)
- **Patient View**: "Early warning system for treatment resistance"
- **Doctor View**: "Resistance Prophet: Predicts resistance 3-6 months early using SAE features + clinical patterns"

---

### 4.3 NEXT-TEST RECOMMENDER MODULE

**Purpose**: Prioritize diagnostic tests to unlock advanced capabilities.

#### Components

**4.3.1 Next Test Recommendations**
- **Backend**: `sae_service.get_next_test_recommendations()`
- **Frontend**: `NextTestCard.jsx`
- **Data Source**: Germline status, tumor context, treatment history
- **Priorities**:
  - **Priority 1: HRD Test** → Unlocks PARP confidence (HRD≥42 → 95%)
  - **Priority 2: ctDNA Panel** → Unlocks TMB, MSI, somatic mutations
  - **Priority 3: SLFN11 IHC** → Refines PARP sensitivity
- **Patient View**: "Next recommended test: HRD test (unlocks personalized predictions)"
- **Doctor View**: "Priority 1: HRD test (MyChoice CDx, 10 days, $4-6K) → Unlocks: PARP confidence, WIWFM"

---

### 4.4 DEPENDENCY SUMMARY: CLINICAL MONITORING

| Capability | Requires | Unlocks |
|------------|----------|---------|
| **CA-125 Intelligence** | CA-125 value | Burden classification, response forecast |
| **Resistance Alert** | SAE features (NGS) OR CA-125 | Early resistance detection |
| **Resistance Playbook** | Tumor context, treatment history | Alternative strategies |
| **Resistance Prophet** | SAE features, treatment history | Early prediction (3-6 months) |
| **Next Test Recommender** | Germline status, tumor context | Test prioritization |

---

## 📈 TIER 5: EVIDENCE/CONFIDENCE LAYER

### What This Layer Does

Provides **confidence calibration**, **evidence tiers**, **SAE features**, and **provenance tracking** to enable transparent clinical decision-making.

### Data Requirements

| Component | Required Data | Optional Data |
|-----------|---------------|---------------|
| **Evidence Tiers** | Drug/trial recommendation | Evidence sources |
| **Confidence Levels** | Completeness score (L0/L1/L2) | Evidence strength |
| **SAE Features** | Tumor context (NGS) | WIWFM response |
| **Provenance** | All orchestration steps | Timestamps, methods |

---

### 5.1 CONFIDENCE CALIBRATION MODULE

**Purpose**: Calibrate confidence based on data completeness and evidence strength.

#### Components

**5.1.1 Completeness-Based Confidence Caps**
- **Backend**: `sporadic_gates.py` (confidence caps)
- **Data Source**: `completeness_score` (L0/L1/L2)
- **Caps**:
  - **L0** (completeness < 0.3): Confidence capped at 0.4
  - **L1** (0.3 ≤ completeness < 0.7): Confidence capped at 0.6
  - **L2** (completeness ≥ 0.7): Confidence uncapped
- **Patient View**: "Confidence level reflects available data quality"
- **Doctor View**: "Completeness 0.2 (L0) → Confidence capped at 0.40 (conservative estimate)"

**5.1.2 Evidence Tiers**
- **Backend**: Integrated into drug/trial recommendations
- **Tiers**: 
  - **Supported**: Strong evidence (RCTs, guidelines)
  - **Consider**: Moderate evidence (pathway alignment, literature)
  - **Insufficient**: Limited evidence
- **Patient View**: "This recommendation is Supported (strong evidence)"
- **Doctor View**: "Olaparib: Tier: Supported, Badges: PathwayAligned, RCT, Guideline"

---

### 5.2 SAE FEATURES MODULE (Explainability)

**Purpose**: Provide interpretable features from genomic data to explain predictions.

#### Data Dependencies

**REQUIRES:**
- ✅ `tumor_context` (NGS) → **UNLOCKS**: SAE features computation
- ✅ WIWFM response → Enhances SAE features
- ✅ CA-125 intelligence → Enhances resistance signals

**WHEN NGS NOT AVAILABLE:**
- ⚠️ Returns `"status": "awaiting_ngs"`

#### Components

**5.2.1 SAE Features Extraction**
- **Backend**: `sae_service.compute_sae_features()`
- **Frontend**: `AyeshaSAEFeaturesCard.jsx`
- **Features** (9 total):
  1. **exon_disruption** → Evo2 delta + hotspot floor
  2. **hotspot_mutation** → AlphaMissense / ClinVar / Hotspot DB
  3. **essentiality_signal** → Insights essentiality endpoint
  4. **DNA_repair_capacity** → Toxicity pathway overlap
  5. **seed_region_quality** → CRISPR guide quality
  6. **cohort_overlap** → Cohort validation
  7. **line_appropriateness** → Treatment line fit
  8. **cross_resistance_risk** → Resistance risk
  9. **sequencing_fitness** → Sequencing score
- **Output**: DNA repair capacity, pathway burden, mechanism vector, boosting/limiting features
- **Patient View**: "DNA Repair Capacity: High (0.85) | Pathway Burden: DDR-high"
- **Doctor View**: "SAE Features: DNA repair capacity 0.85, Mechanism vector [DDR=1.4, ...], Boosting features: 3, Limiting features: 0"

---

### 5.3 PROVENANCE MODULE

**Purpose**: Track data sources, methods, and timestamps for auditability.

#### Components

**5.3.1 Provenance Tracking**
- **Backend**: Integrated into all orchestrator responses
- **Frontend**: Displayed in component metadata
- **Tracks**:
  - Data sources (NGS report, CA-125 labs, etc.)
  - Methods (S/P/E framework, mechanism fit, etc.)
  - Timestamps (when data was processed)
  - Confidence levels (how confident we are)
- **Patient View**: "Recommendations based on your test results from [date]"
- **Doctor View**: "Provenance: NGS from Foundation Medicine (2025-01-15), WIWFM using S/P/E framework, Confidence: 0.85 (L2, uncapped)"

---

### 5.4 DEPENDENCY SUMMARY: EVIDENCE/CONFIDENCE

| Capability | Requires | Unlocks |
|------------|----------|---------|
| **Confidence Caps** | Completeness score (L0/L1/L2) | Calibrated confidence levels |
| **Evidence Tiers** | Drug/trial recommendation | Evidence strength classification |
| **SAE Features** | Tumor context (NGS) | Interpretable feature explanations |
| **Provenance** | All orchestration steps | Auditability, transparency |

---

## 🔄 CLINICAL WORKFLOW: DATA DEPENDENCY PROGRESSION

### Week 1: Initial Intake (L0 Data Completeness)

**Available Data:**
- ✅ Stage: IVB
- ✅ Disease: Ovarian cancer (HGSOC)
- ✅ Germline: MBD4 mutation (positive)
- ⚠️ Tumor Context: **NOT AVAILABLE** (awaiting NGS)
- ⚠️ CA-125: **NOT AVAILABLE** (optional)

**Unlocked Capabilities:**
1. ✅ **SOC Recommendation** → NCCN Category 1 (Carboplatin + Paclitaxel + Bevacizumab)
2. ✅ **Trial Matching** (basic) → Stage-based matching, mechanism-fit limited
3. ✅ **Next Test Recommender** → HRD test (Priority 1), ctDNA (Priority 2)
4. ✅ **Hint Tiles** → Actionable hints (max 4)
5. ✅ **Mechanism Map** (basic) → Pathway visualization (limited without NGS)
6. ✅ **PGx Safety Gates** → If germline variants provided
7. ⚠️ **WIWFM Drug Efficacy** → `"status": "awaiting_ngs"`
8. ⚠️ **SAE Features** → `"status": "awaiting_ngs"`
9. ⚠️ **Resistance Alert** → `"status": "awaiting_ngs"`

**Confidence Level**: **L0** (completeness < 0.3) → Confidence capped at **0.40**

**Patient View**: "We can provide basic recommendations. To unlock personalized predictions, order HRD test + ctDNA panel."

**Doctor View**: "L0 intake → Confidence capped at 0.40. Recommendations are guideline-based. Order HRD test (Priority 1) to unlock PARP confidence adjustment."

---

### Week 2-3: NGS Fast-Track (L1 → L2 Data Completeness)

**Available Data (After NGS Returns):**
- ✅ Stage: IVB
- ✅ Disease: Ovarian cancer (HGSOC)
- ✅ Germline: MBD4 mutation (positive)
- ✅ **Tumor Context (NGS)**: 
  - Somatic mutations: TP53 R175H
  - HRD score: 58 (≥42 → PARP rescue)
  - TMB: 25 (≥20 → IO boost)
  - MSI: MSS
  - Completeness: 0.9 (L2)
- ✅ CA-125: 2,842 U/mL

**Unlocked Capabilities (NEW):**
1. ✅ **WIWFM Drug Efficacy** → Full S/P/E rankings (Olaparib #1, efficacy 0.85)
2. ✅ **SAE Features** → DNA repair capacity, mechanism vector, boosting/limiting features
3. ✅ **Resistance Alert** → Early resistance detection (if triggers present)
4. ✅ **Mechanism Map** (detailed) → Full pathway visualization
5. ✅ **PARP Confidence Rescue** → HRD≥42 → PARP multiplier 1.0x (no penalty)
6. ✅ **IO Boost** → TMB≥20 → IO multiplier 1.35x
7. ✅ **CA-125 Intelligence** → Burden classification, response forecast

**Confidence Level**: **L2** (completeness ≥ 0.7) → Confidence **UNCAPPED** (0.85-0.90 for PARP)

**Patient View**: "Personalized predictions now available. Top recommendation: Olaparib (85% confidence, supported by strong evidence)."

**Doctor View**: "L2 intake → Confidence uncapped. PARP rescue (HRD 58 ≥42) → Efficacy 0.85, Confidence 0.85. IO boost (TMB 25 ≥20) → IO multiplier 1.35x."

---

## 👥 DUAL PERSPECTIVE: PATIENT vs DOCTOR PORTAL

### PATIENT PORTAL VIEW

**Focus**: Clear, actionable information with plain-language explanations.

**Display Priorities:**
1. **What this means for me** (biological intelligence in simple terms)
2. **What to do next** (next test recommendations, treatment options)
3. **What to watch for** (monitoring, resistance signals)
4. **How confident we are** (confidence levels explained simply)

**Language**: Plain language, avoid jargon, focus on actionable steps.

**Example (MBD4 Mutation):**
- **Patient View**: "Your MBD4 mutation creates a PARP vulnerability - this is good news for treatment. Your tumor is vulnerable to PARP inhibitors (like Olaparib), which work by blocking DNA repair. This is similar to how PARP inhibitors work for BRCA mutations."

---

### DOCTOR PORTAL VIEW

**Focus**: Clinical-ready intelligence with evidence citations and detailed reasoning.

**Display Priorities:**
1. **Biological mechanism** (pathway analysis, synthetic lethality)
2. **Evidence base** (RCTs, guidelines, literature citations)
3. **Confidence reasoning** (why we're confident, what data supports it)
4. **Action items** (order sets, monitoring protocols, trial contacts)

**Language**: Clinical terminology, evidence citations, detailed reasoning.

**Example (MBD4 Mutation):**
- **Doctor View**: "MBD4 frameshift → BER deficiency → DDR pathway disruption (1.0). Combined with TP53 R175H → Additional DDR contribution (0.8). Synthetic lethality: BER loss + HR pathway stress → PARP inhibitor vulnerability. Evidence: SOLO-1 (HR 0.30), PRIMA (HRD+ benefit). Confidence: 85% (HRD 58 ≥42 → PARP rescue). Recommendation: Olaparib maintenance (NCCN Category 1)."

---

## 🖥️ FRONTEND WIRING STATUS: MODULE-BASED ORGANIZATION

### PAGE 1: AyeshaCompleteCare (`/ayesha-complete-care`)

**Backend API**: `POST /api/ayesha/complete_care_v2`

**Module Organization** (by medical hierarchy):

#### GENOMIC FOUNDATION LAYER
- **VUSResolutionCard** → VUS resolution (separate API)
- **EssentialityScoreDisplay** → Essentiality scores (separate API)

#### PATHWAY/MECHANISM LAYER
- **SyntheticLethalityCard** → Synthetic lethality (separate API)
- **MechanismChips** → Mechanism map visualization

#### THERAPEUTIC INTELLIGENCE LAYER
- **SOCRecommendationCard** → SOC recommendation
- **DrugRankingPanel** → WIWFM drug efficacy (conditional on NGS)
- **TrialMatchesCard** → Clinical trials
- **IOSafestSelectionCard** → IO selection (RUO)
- **FoodRankingPanel** → Food/supplement validation

#### CLINICAL MONITORING LAYER
- **CA125Tracker** → CA-125 intelligence
- **NextTestCard** → Next test recommendations
- **ResistancePlaybook** → Resistance strategies

#### EVIDENCE/CONFIDENCE LAYER
- **HintTilesPanel** → Hint tiles
- **IntegratedConfidenceBar** → Integrated confidence (conditional)

---

### PAGE 2: AyeshaTrialExplorer (`/ayesha-trials`)

**Backend API**: `POST /api/ayesha/complete_care_v2`

**Tab Organization** (by clinical workflow):

#### TAB 0: OVERVIEW
- **GENOMIC FOUNDATION**: (none displayed)
- **PATHWAY/MECHANISM**: MechanismChips, AyeshaSAEFeaturesCard (conditional)
- **THERAPEUTIC INTELLIGENCE**: SOCRecommendationCard
- **CLINICAL MONITORING**: NextTestCard
- **EVIDENCE/CONFIDENCE**: HintTilesPanel

#### TAB 1: TRIALS
- **THERAPEUTIC INTELLIGENCE**: TrialMatchCard (top 10), TrialSafetyGate (PGx)

#### TAB 2: TREATMENT
- **THERAPEUTIC INTELLIGENCE**: SOCRecommendationCard, DrugRankingPanel (conditional on NGS), FoodRankingPanel

#### TAB 3: MONITORING
- **CLINICAL MONITORING**: CA125Tracker, NextTestCard

#### TAB 4: RESISTANCE
- **CLINICAL MONITORING**: ResistanceAlertBanner (conditional), ResistancePlaybook (conditional), Resistance Prophet UI (conditional)
- **EVIDENCE/CONFIDENCE**: AyeshaSAEFeaturesCard (conditional)

---

## 📋 CAPABILITY MODULE MATRIX

| Module | Backend Service | Frontend Component | Data Requirement | Unlocks |
|--------|----------------|-------------------|------------------|---------|
| **GENOMIC FOUNDATION** |
| Germline PGx | `pgx_screening_service` | `SafetyGateCard`, `TrialSafetyGate` | `germline_variants` | Composite scoring, trial safety gates |
| VUS Resolution | `vus_identification_service` | `VUSResolutionCard` | Variant (gene + HGVS) | Pathogenicity classification |
| Essentiality | `insights_service` | `EssentialityScoreDisplay` | Gene name | Dependency analysis |
| **PATHWAY/MECHANISM** |
| Mechanism Vector | `sae_service.extract_mechanism_vector()` | `MechanismChips` | Somatic mutations OR WIWFM | Mechanism-fit trial matching |
| Synthetic Lethality | `synthetic_lethality_service` | `SyntheticLethalityCard` | Mutations (gene pairs) | Combination therapy recommendations |
| **THERAPEUTIC INTELLIGENCE** |
| WIWFM Drug Efficacy | `wiwfm_service` | `DrugRankingPanel` | `tumor_context` (NGS) | Drug rankings, mechanism vector |
| Trial Matching | `trial_matching_agent` | `TrialMatchesCard`, `TrialMatchCard` | Stage + mechanism vector | Trial matches (holistic score) |
| SOC Recommendation | `soc_service` | `SOCRecommendationCard` | Stage, disease | NCCN-aligned regimen |
| Food/Supplement | `food_service` | `FoodRankingPanel` | SOC drugs, treatment line | Supplement recommendations |
| **CLINICAL MONITORING** |
| CA-125 Intelligence | `ca125_intelligence_service` | `CA125Tracker` | CA-125 value | Burden classification, response forecast |
| Resistance Alert | `sae_service.detect_resistance()` | `ResistanceAlertBanner` | SAE features (NGS) OR CA-125 | Early resistance detection |
| Resistance Playbook | `resistance_playbook_service` | `ResistancePlaybook` | Tumor context, treatment history | Alternative strategies |
| Resistance Prophet | `resistance_prophet_service` | Custom UI | SAE features, treatment history | Early prediction (3-6 months) |
| Next Test Recommender | `sae_service.get_next_test_recommendations()` | `NextTestCard` | Germline status, tumor context | Test prioritization |
| **EVIDENCE/CONFIDENCE** |
| SAE Features | `sae_service.compute_sae_features()` | `AyeshaSAEFeaturesCard` | `tumor_context` (NGS) | Interpretable features |
| Hint Tiles | `sae_service.get_hint_tiles()` | `HintTilesPanel` | Context-dependent | Actionable hints (max 4) |

---

## 🎯 SUMMARY: MEDICAL HIERARCHY ORGANIZATION

**Tier 1: Genomic Foundation** → Establishes genetic basis (germline, somatic, variants)  
**Tier 2: Pathway/Mechanism** → Translates mutations to pathway disruptions  
**Tier 3: Therapeutic Intelligence** → Converts pathway insights to treatment recommendations  
**Tier 4: Clinical Monitoring** → Tracks treatment response and resistance  
**Tier 5: Evidence/Confidence** → Calibrates confidence and provides explainability  

**Key Dependencies:**
- **NGS (`tumor_context`)** → Unlocks WIWFM, SAE Features, Resistance Prediction
- **HRD Score (≥42)** → Unlocks PARP confidence rescue
- **TMB (≥20) or MSI-H** → Unlocks IO boost
- **Completeness (L2: ≥0.7)** → Unlocks confidence uncapping
- **Germline Variants** → Unlocks PGx safety gates

**Patient vs Doctor Perspectives:**
- **Patient Portal**: Plain language, actionable steps, confidence explained simply
- **Doctor Portal**: Clinical terminology, evidence citations, detailed reasoning

---

## 🎨 FRONTEND ARCHITECTURE: CURRENT STATE & IMPROVEMENT PLAN

### Current State

**Two Monolithic Pages:**
1. **AyeshaCompleteCare.jsx** (962 lines) - Unified page for Ayesha's complete care plan
2. **UniversalCompleteCare.jsx** (695 lines) - Universal complete care plan page

**Shared Backend Endpoint:**
- Both pages use `POST /api/ayesha/complete_care_v2` (unified orchestrator)

**Problems Identified:**
- ❌ **Code Duplication**: Similar API request building, data transformation, error handling
- ❌ **Monolithic Structure**: All logic (API calls, state, rendering) in single component
- ❌ **No Shared Logic**: Request building, data transformation duplicated
- ❌ **Mixed Concerns**: API logic, state management, UI rendering intertwined
- ❌ **Hard to Test**: Large components difficult to unit test
- ❌ **Hard to Reuse**: Can't reuse logic for other patient pages

### Improvement Plan

**See**: `FRONTEND_ARCHITECTURE_IMPROVEMENT_PLAN.md` for detailed refactoring plan.

**Key Improvements:**
1. **Extract Shared Logic** → Custom hooks (`useCompleteCarePlan`, `useCompleteCareRequest`)
2. **Modular Components** → Organize by medical hierarchy (5 section components)
3. **Container/Presenter Pattern** → Separate data fetching from presentation
4. **Reusable Building Blocks** → Shared components, hooks, utilities
5. **Maintainability** → Easier to test, modify, extend

**Expected Outcomes:**
- ✅ **70% code reuse** (hooks, utilities)
- ✅ **80%+ test coverage** (hooks, utilities, components)
- ✅ **Modular architecture** (medical hierarchy organization)
- ✅ **Thin containers** (pages reduced to ~100-150 lines each)

---

**Last Updated**: January 29, 2025  
**Reorganization**: Medical hierarchy structure with dependency mapping and dual perspectives  
**Frontend Plan**: See `FRONTEND_ARCHITECTURE_IMPROVEMENT_PLAN.md`
