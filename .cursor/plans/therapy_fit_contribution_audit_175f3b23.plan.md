---
name: Therapy Fit Contribution Audit
overview: Comprehensive audit of the therapy_fit_contribution.mdc document, understanding the S/P/E framework mechanics and insights chips system for drug efficacy ranking.
todos:
  - id: audit_spe_framework
    content: "Understand S/P/E framework formula: efficacy_score = 0.3×S + 0.4×P + 0.3×E + ClinVar_prior"
    status: pending
  - id: audit_sequence_component
    content: "Understand Sequence (S) component: Evo2 adaptive multi-window scoring (30% weight)"
    status: pending
  - id: audit_pathway_component
    content: "Understand Pathway (P) component: Gene-to-pathway mapping (40% weight - highest)"
    status: pending
  - id: audit_evidence_component
    content: "Understand Evidence (E) component: Literature + ClinVar (30% weight)"
    status: pending
  - id: audit_insights_chips
    content: "Understand 4 insights chips: Functionality, Chromatin, Essentiality, Regulatory with threshold-based lifts"
    status: pending
  - id: audit_confidence_computation
    content: "Understand confidence computation: tier-based or linear S/P/E formula with insights lifts"
    status: pending
  - id: audit_evidence_tiers
    content: "Understand evidence tiers: Supported/Consider/Insufficient with criteria and confidence ranges"
    status: pending
  - id: audit_badges
    content: "Understand badges: RCT, Guideline, ClinVar-Strong, PathwayAligned"
    status: pending
---

# Therapy Fit Contribution Audit: S/P/E Framework & Insights Chips

## Executive Summary

The Therapy Fit capability ranks drugs by efficacy using a transparent **S/P/E (Sequence/Pathway/Evidence) Framework** with **30%/40%/30% weights**. The system includes **4 insights chips** (Functionality, Chromatin, Essentiality, Regulatory) that provide threshold-based confidence lifts to differentiate drugs with similar S/P/E scores.---

## 1. S/P/E Framework Architecture

### Core Formula

```javascript
efficacy_score = 0.3 × Sequence + 0.4 × Pathway + 0.3 × Evidence + ClinVar_prior
```

**Location**: `drug_scorer.py:171` (hardcoded, not using config defaults)**Weights Breakdown**:

- **Sequence (S)**: 30% - Evo2 adaptive multi-window scoring
- **Pathway (P)**: 40% - Gene-to-pathway mapping (HIGHEST WEIGHT - mechanism alignment importance)
- **Evidence (E)**: 30% - Literature + ClinVar classification
- **ClinVar_prior**: Additive boost [-0.2, +0.2] based on ClinVar classification

### Sequence (S) Component - 30% Weight

**Implementation**: `sequence_processor.py`, `evo2_scorer.py`**How It Works**:

1. **Evo2 Adaptive Multi-Window Scoring**:

- Windows: 4096, 8192, 16384, 25000 bp
- Adaptive: Select best window per variant
- Gene-specific calibration (hotspot-aware)
- Output: Variant disruption percentile (0-1)

2. **Fallback Chain**:

- Priority 1: Fusion Engine (AlphaMissense) - ONLY for GRCh38 missense variants
- Priority 2: Evo2 Adaptive - Multi-window, ensemble models
- Priority 3: Massive Oracle - Synthetic or real-context (if enabled)

**Key Files**:

- `api/services/sequence_scorers/evo2_scorer.py`
- `api/services/sequence_scorers/sequence_processor.py`

### Pathway (P) Component - 40% Weight (HIGHEST)

**Implementation**: Pathway aggregation with gene-to-pathway mapping**How It Works**:

1. **Gene-to-Pathway Mapping**:

- Map variants to pathways: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux
- Examples:
    - BRCA1 truncation → HR pathway (DDR)
    - MBD4 frameshift → BER pathway (DDR)
    - TP53 R175H → Checkpoint pathway (DDR)

2. **Drug-Pathway Alignment**:

- PARP inhibitors → DDR pathway (0.9 alignment)
- MEK inhibitors → MAPK pathway (0.9 alignment)
- Proteasome inhibitors → TP53 pathway (0.8 alignment)

3. **Output**: Pathway alignment percentile (0-1)

**Why 40%?**: Highest weight reflects mechanism alignment importance - validates that drug mechanism matches patient pathway burden.

### Evidence (E) Component - 30% Weight

**Implementation**: Literature search + ClinVar classification**How It Works**:

1. **ClinVar Classification**:

- Pathogenic/Likely Pathogenic → ClinVar_prior boost
- Benign/Likely Benign → ClinVar_prior penalty
- Review status affects confidence

2. **Literature Search**:

- PubMed/OpenAlex/S2 search per drug
- MoA-filtered queries
- Evidence strength score (0-1)

3. **Output**: Evidence strength score (0-1)

---

## 2. Insights Chips System

### Four Chips Overview

The system computes **4 biological insight scores** that provide **threshold-based confidence lifts**:| Chip | Threshold | Lift (Legacy) | Lift (V2) | Purpose ||------|-----------|---------------|-----------|---------|| **Functionality** | ≥0.6 | +0.05 | +0.04 | Protein function change prediction || **Chromatin** | ≥0.5 | +0.03 | +0.04 | Regulatory impact assessment || **Essentiality** | ≥0.7 | +0.07 | +0.02 | Gene dependency scoring || **Regulatory** | ≥0.6 | +0.02 | +0.02 | Splicing and non-coding impact |**Implementation**: `api/services/confidence/insights_lifts.py`

### How Chips Work

1. **Insight Computation**:

- Each chip is computed via separate API endpoints:
    - `/api/insights/predict_protein_functionality_change`
    - `/api/insights/predict_chromatin_accessibility`
    - `/api/insights/predict_gene_essentiality`
    - `/api/insights/predict_splicing_regulatory`

2. **Threshold-Based Lifts**:

- Lifts only applied when chip scores exceed defined thresholds
- Modest boosts (+0.02 to +0.07) help differentiate drugs with similar S/P/E scores
- **Purpose**: When biological signals align, confidence increases modestly

3. **Total Lift Cap** (V2 only):

- Cap total lifts at +0.08
- If sum exceeds 0.08, proportionally scale down all lifts

**Key Files**:

- `api/services/confidence/insights_lifts.py`
- `api/routers/insights.py`

### Real-World Example (Ayesha MBD4+TP53)

**Insight Chips**:

- Functionality = 0.60 (≥0.6 threshold → +0.05 lift)
- Chromatin = 0.60 (≥0.5 threshold → +0.04 lift)
- Essentiality = 0.35 (<0.7 threshold → no lift)
- Regulatory = 0.12 (<0.6 threshold → no lift)

**Total Lifts Applied**: +0.09 (Functionality + Chromatin)---

## 3. Confidence Computation

### Confidence Formula

**Legacy** (tier-based):

- **Supported**: 0.6 + 0.2 × max(seq_pct, path_pct)
- **Consider**: 0.3 + 0.1 × seq_pct + 0.1 × path_pct
- **Insufficient**: 0.20 + 0.35 × max_sp + 0.15 × min_sp

**V2** (linear S/P/E formula):

```javascript
confidence = clamp01(0.5 × S + 0.2 × P + 0.3 × E + lifts)
```

Where:

- **S**: Sequence percentile (0-1)
- **P**: Pathway percentile (0-1)
- **E**: Evidence tier score (Supported: +0.05, Consider: +0.02, Insufficient: +0.00)
- **lifts**: Sum of insights lifts (capped at +0.08 in V2)

**Implementation**: `api/services/confidence/confidence_computation.py`---

## 4. Evidence Tiers & Badges

### Evidence Tiers

| Tier | Criteria | Confidence Range | Example ||------|----------|------------------|---------|| **Supported** | Strong evidence (≥0.7) OR ClinVar-Strong + pathway alignment | 0.6+ | KRAS G12D → MEK inhibitor (0.85 confidence) || **Consider** | Moderate evidence (≥0.3) with some pathway alignment | 0.3-0.6 | BRAF V600E → BRAF inhibitor (0.51 confidence) || **Insufficient** | Low scores (<0.3) across S/P/E | <0.3 | Fast-path mode or weak signals |

### Badges (Evidence Quality Indicators)

| Badge | Description | Example ||-------|-------------|---------|| **RCT** | Randomized controlled trial evidence | PARP inhibitors for BRCA-mutant ovarian cancer || **Guideline** | NCCN/FDA guideline recommendation | BRAF V600E → BRAF inhibitor (melanoma) || **ClinVar-Strong** | ClinVar Pathogenic/Likely Pathogenic | TP53 R248W → Proteasome inhibitor || **PathwayAligned** | Drug mechanism matches patient pathway burden | DDR-high → PARP inhibitors |---

## 5. Complete Data Flow

```javascript
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
    └─ Final: efficacy_score = 0.3×S + 0.4×P + 0.3×E + clinvar_prior
  ↓
[6] Confidence Modulation
    ├─ compute_evidence_tier() → Supported/Consider/Insufficient
    ├─ compute_insights_lifts() → Functionality/Chromatin/Essentiality/Regulatory lifts
    └─ compute_confidence() → Final confidence (0-1)
  ↓
[7] Output: Ranked drugs with:
    ├─ efficacy_score (0-1)
    ├─ confidence (0-1)
    ├─ evidence_tier (Supported/Consider/Insufficient)
    ├─ badges[] (RCT, Guideline, ClinVar-Strong, PathwayAligned)
    ├─ insights {functionality, chromatin, essentiality, regulatory}
    ├─ rationale[] (S/P/E breakdown per drug)
    ├─ citations[] (literature PMIDs)
    └─ provenance {run_id, profile, methods, flags}
```

---

## 6. Key Implementation Files

### Core S/P/E Framework

- `api/services/drug_scorer.py` - Drug scoring with S/P/E formula (line 171)
- `api/services/efficacy_orchestrator.py` - Main orchestration
- `api/services/sequence_scorers/sequence_processor.py` - Sequence (S) component
- `api/services/sequence_scorers/evo2_scorer.py` - Evo2 scoring

### Insights Chips

- `api/services/confidence/insights_lifts.py` - Confidence lift computation
- `api/routers/insights.py` - 4 insights endpoints (Functionality, Chromatin, Essentiality, Regulatory)
- `api/services/confidence/confidence_computation.py` - Confidence formula

### Documentation

- `.cursor/lectures/drugDevelopment/therapy_fit_contribution.mdc` - This document
- `.cursor/ayesha/Deliverables/Iterations/I3_SPE_FRAMEWORK.md` - Detailed implementation notes

---

## 7. Validated Metrics (Production-Ready)

### S/P/E Framework Performance

- **S/P/E Weights**: 30% S, 40% P, 30% E
- **Pathway Alignment**: 100% (Multiple Myeloma - 5/5 MAPK variants correctly matched)
- **Overall Accuracy**: 70-85%
- **Confidence Range**: 0.45-0.85 (MM baseline)
- **Evidence Tier Promotions**: 10-20% Consider→Supported

### Real-World Validation Examples

**Multiple Myeloma (100% Pathway Alignment)**:

- KRAS G12D → MEK inhibitor (0.85 confidence, **supported** tier, PathwayAligned badge)
- NRAS Q61R → MEK inhibitor (0.83 confidence, **supported** tier, PathwayAligned badge)
- TP53 R248W → Proteasome inhibitor (0.84 confidence, **supported** tier, ClinVar-Strong badge)

**Ayesha (MBD4+TP53 HGSOC)**:

- PARP inhibitors ranked #1-3 (0.800 efficacy)
- **Supported** tier, PathwayAligned + ClinVar-Strong badges
- Insight chips: Functionality=0.60, Chromatin=0.60, Essentiality=0.35, Regulatory=0.12