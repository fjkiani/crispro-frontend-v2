---
name: Drug Efficacy Chips Audit
overview: Comprehensive audit of the Drug Efficacy Prediction (S/P/E Framework) and Insight Chips mechanism, covering how chips are computed, integrated into confidence scoring, and displayed in the UI.
todos:
  - id: audit-complete
    content: Complete comprehensive audit of drug efficacy framework and insight chips mechanism
    status: completed
  - id: understand-chips
    content: Understand how four insight chips (Functionality, Regulatory, Essentiality, Chromatin) are computed and integrated
    status: completed
  - id: trace-confidence
    content: Trace how chips modulate confidence scores through threshold-based lifts
    status: completed
  - id: document-flow
    content: Document complete data flow from chip computation → confidence integration → UI display
    status: completed
---

# Drug Efficacy Framework & Insight Chips - Complete Audit

## Executive Summary

The Drug Efficacy Prediction system uses an **S/P/E Framework** (Sequence/Pathway/Evidence) to rank drugs by mechanism alignment. Four **Insight Chips** (Functionality, Regulatory, Essentiality, Chromatin) provide transparent mechanistic signals that modulate confidence scores. This audit covers the complete pipeline from chip computation to UI display.---

## 1. S/P/E Framework Architecture

### Core Formula

```javascript
efficacy_score = 0.3×S + 0.4×P + 0.3×E + clinvar_prior
```

**Components:**

- **Sequence (S)**: Evo2 multi-window variant scoring (4096, 8192, 16384, 25000 bp windows)
- Output: `seq_pct` (calibrated percentile, 0-1)
- Location: [`drug_scorer.py:42-43`](oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py)
- **Pathway (P)**: Aggregated pathway burden mapping variants to 7 pathways (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
- Output: `path_pct` (normalized pathway score, 0-1)
- Normalization: `path_pct = min(1.0, s_path / 0.005)` where `s_path` is weighted pathway score
- Location: [`drug_scorer.py:45-57`](oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py)
- **Evidence (E)**: Literature search + ClinVar priors
- Output: `s_evd` (evidence strength, 0-1)
- Location: [`drug_scorer.py:68-74`](oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py)

**Key Finding**: Pathway (P) component is **ESSENTIAL** — 100% accuracy with P, 40% without (validated ablation study).---

## 2. Insight Chips System

### Four Chips Overview

| Chip | Description | Typical Range | API Endpoint | Threshold | Confidence Lift ||------|-------------|---------------|--------------|-----------|----------------|| **Functionality** | Evo2-based protein disruption impact | 0.60 (moderate) | `POST /api/insights/predict_protein_functionality_change` | ≥0.6 | +0.05 (legacy) / +0.04 (V2) || **Regulatory** | Splice/UTR disruption signals | 0.12 (low) | `POST /api/insights/predict_splicing_regulatory` | ≥0.6 | +0.02 || **Essentiality** | Gene-level dependency | 0.35 (moderate) | `POST /api/insights/predict_gene_essentiality` | ≥0.7 | +0.07 (legacy) / +0.02 (V2) || **Chromatin** | Regulatory accessibility (Enformer-ready) | 0.58-0.60 (moderate-high) | `POST /api/insights/predict_chromatin_accessibility` | ≥0.5 | +0.04 (legacy) / +0.02 (V2) |

### Chip Computation Pipeline

**Step 1: Insights Bundle Collection**

- Location: [`I3_SPE_FRAMEWORK.md:529-533`](.cursor/ayesha/Deliverables/Iterations/I3_SPE_FRAMEWORK.md)
- Process: Parallel API calls to four insight endpoints
- Output: `InsightsBundle` object with four scores (0-1 range)

**Step 2: Integration into Drug Scoring**

- Location: [`drug_scorer.py:148-154`](oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py)
- Process: Chips extracted into `insights_dict` and passed to `compute_confidence()`
- Code:
  ```python
          insights_dict = {
              "functionality": insights.functionality or 0.0,
              "chromatin": insights.chromatin or 0.0,
              "essentiality": insights.essentiality or 0.0,
              "regulatory": insights.regulatory or 0.0,
          }
          confidence = compute_confidence(tier, seq_pct, path_pct, insights_dict, confidence_config)
  ```


---

## 3. Confidence Calculation with Chips

### Legacy Implementation (Default)

**Location**: [`confidence_computation.py:9-108`](oncology-coPilot/oncology-backend-minimal/api/services/confidence/confidence_computation.py)**Base Confidence by Tier:**

- **Supported**: `0.6 + 0.2 * max(seq_pct, path_pct)`
- **Consider**: `0.3 + 0.1 * seq_pct + 0.1 * path_pct` (or `0.5 + 0.2 * max(seq_pct, path_pct)` if Fusion active and S/P ≥ 0.7)
- **Insufficient**: `0.20 + 0.35 * max_sp + 0.15 * min_sp` (with 0.25 floor if Fusion active)

**Chip-Based Confidence Lifts** (lines 94-98):

```python
confidence += 0.05 if func >= 0.6 else 0.0      # Functionality
confidence += 0.04 if chrom >= 0.5 else 0.0    # Chromatin
confidence += 0.07 if ess >= 0.7 else 0.0       # Essentiality
confidence += 0.02 if reg >= 0.6 else 0.0       # Regulatory
```

**Additional Modulators:**

- Alignment margin boost: +0.05 if `abs(seq_pct - path_pct) >= 0.2`
- ClinVar prior boost: +0.1 (capped) when aligned with pathway
- Gene-drug MoA tie-breaker: +0.01 for specific matches (e.g., BRAF → BRAF inhibitor)

### V2 Implementation (CONFIDENCE_V2=1)

**Location**: [`confidence_computation.py:111-160`](oncology-coPilot/oncology-backend-minimal/api/services/confidence/confidence_computation.py)**Formula**: `confidence = clamp01(0.5·S + 0.2·P + 0.3·E + lifts)`**Chip Lifts** (lines 147-150):

- Functionality ≥0.6: +0.04
- Chromatin ≥0.5: +0.02
- Essentiality ≥0.7: +0.02
- Regulatory ≥0.6: +0.02
- **Total lifts capped at +0.08** (proportionally scaled if exceeded)

**Location**: [`insights_lifts.py:61-99`](oncology-coPilot/oncology-backend-minimal/api/services/confidence/insights_lifts.py)---

## 4. Frontend Display

### InsightChips Component

**Location**: [`InsightChips.jsx`](oncology-coPilot/oncology-frontend/src/components/vus/InsightChips.jsx)**Visual Logic:**

- **High** (≥threshold): Green chip (`bg-green-600`) with ✓ icon
- **Medium** (≥mid threshold): Yellow chip (`bg-yellow-600`) with ~ icon
- **Low** (<mid threshold): Gray chip (`bg-gray-600`) with — icon

**Data Flow:**

1. Chips received via `insights` prop: `{functionality, chromatin, essentiality, regulatory}`
2. Values compared against `INSIGHT_THRESHOLDS` from constants
3. Color/icon determined by threshold comparison
4. Helper tooltips display gene/variant-specific explanations

**Key Features:**

- KB integration for gene/variant-specific helper text
- Provenance display for auditability
- Threshold legend showing High/Medium/Low ranges

---

## 5. Validated Examples

### Multiple Myeloma (100% Pathway Alignment)

- **KRAS G12D → MEK inhibitor**: 0.85 confidence, **supported** tier
- **Insight chips**: Functionality=0.60, Chromatin=0.60, Essentiality=0.35, Regulatory=0.12
- **Confidence calculation**: Base (supported tier) + Functionality lift (+0.05) + Chromatin lift (+0.04) = **0.85**

### Ayesha (MBD4+TP53 HGSOC)

- **PARP inhibitors ranked #1-3**: 0.800 efficacy
- **Mechanism vector**: [0.88, 0.12, 0.15, 0.10, 0.05, 0.0, 0.0] (DDR-high)
- **Insight chips**: Functionality=0.60, Chromatin=0.60, Essentiality=0.35, Regulatory=0.12
- **Evidence tier**: Supported (for PARP inhibitors with HRD-high)

---

## 6. Key Implementation Files

### Backend

- **Drug Scoring**: [`drug_scorer.py`](oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py) (lines 148-154, 187, 220)
- **Confidence Computation**: [`confidence_computation.py`](oncology-coPilot/oncology-backend-minimal/api/services/confidence/confidence_computation.py) (lines 94-98, 147-150)
- **Insights Lifts**: [`insights_lifts.py`](oncology-coPilot/oncology-backend-minimal/api/services/confidence/insights_lifts.py) (lines 49-56, 82-89)

### Frontend

- **Chip Display**: [`InsightChips.jsx`](oncology-coPilot/oncology-frontend/src/components/vus/InsightChips.jsx) (lines 41-51, 56-92)

### Documentation

- **Main Contribution Doc**: [`drug_efficacy_contribution.mdc`](.cursor/lectures/drugDevelopment/drug_efficacy_contribution.mdc) (lines 36-43, 219-223, 331-343)
- **Framework Flow**: [`I3_SPE_FRAMEWORK.md`](.cursor/ayesha/Deliverables/Iterations/I3_SPE_FRAMEWORK.md) (lines 529-546)

---

## 7. Current Limitations

| Limitation | Status | Impact on Chips ||-----------|--------|----------------|| Chromatin predictions | Stubs (Enformer-ready) | Chromatin chip uses heuristic v0, not Enformer || Calibration Error (ECE) | 0.479 (SPE), 0.529 (SP) | Chip lifts may not be optimally calibrated || Small sample size | n=5 MAPK variants (MM) | Chip thresholds validated on limited data |---

## 8. Summary: How Chips Work

1. **Computation**: Four parallel API calls compute chip scores (0-1 range)
2. **Integration**: Chips bundled into `insights_dict` and passed to confidence calculator
3. **Confidence Modulation**: Chips provide "small, transparent lifts" when thresholds exceeded:

- Functionality ≥0.6: +0.05/+0.04 lift
- Chromatin ≥0.5: +0.04/+0.02 lift
- Essentiality ≥0.7: +0.07/+0.02 lift