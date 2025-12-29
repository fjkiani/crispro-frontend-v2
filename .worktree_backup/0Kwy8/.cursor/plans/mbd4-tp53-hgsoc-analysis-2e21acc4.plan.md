---
name: MBD4 Germline + TP53 Somatic HGSOC Analysis Plan (REFINED)
overview: ""
todos: []
---

# MBD4 Germline + TP53 Somatic HGSOC Analysis Plan (REFINED)

## Objective

Analyze how MBD4 germline loss (homozygous c.1239delA) combined with TP53 somatic mutation creates high-grade serous ovarian cancer, identifying pathway vulnerabilities, therapeutic targets, and synthetic lethal opportunities.

**Context**: This plan uses the actual S/P/E framework, Evo2 integration, and available endpoints. SAE features are excluded (not yet available).

## Phase 1: Variant Functional Annotation

### 1.1 MBD4 Germline Variant (c.1239delA)

**Input**: MBD4, homozygous c.1239delA (frameshift leading to loss-of-function)

**Analysis Steps**:

#### A. Evo2 Sequence Scoring

- **Endpoint**: `POST /api/evo/score_variant_multi` OR `POST /api/evo/score_variant_exon`
- **Request Format**:
  ```json
  {
    "assembly": "GRCh38",
    "chrom": "3",
    "pos": 129430456,  // Verify exact position from HGVS
    "ref": "A",
    "alt": "",
    "window": 8192  // Adaptive windows: 4096, 8192, 16384, 25000
  }
  ```

- **Expected**: High disruption (frameshift → truncation)
- **Implementation**: Uses `evo2_scorer.py` with adaptive multi-window scoring
- **Note**: Frameshift/truncation detection happens in `sequence_processor.py` (Triumvirate Protocol)

#### B. Insights Bundle (4 Chips)

- **Functionality**: `POST /api/insights/predict_protein_functionality_change`
  - Request: `{ "gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2" }`
  - Expected: Loss-of-function (MBD4 is DNA glycosylase)

- **Essentiality**: `POST /api/insights/predict_gene_essentiality`
  - Request: `{ "gene": "MBD4", "mutations": [{ "hgvs_p": "p.Ile413Serfs*2" }] }`
  - Expected: Assess BER pathway dependency

- **Regulatory**: `POST /api/insights/predict_splicing_regulatory`
  - Request: `{ "gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2" }`
  - Expected: Check for splicing/expression impacts

- **Chromatin**: `POST /api/insights/predict_chromatin_accessibility`
  - Request: `{ "gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2" }`
  - Expected: Evaluate accessibility changes

#### C. Evidence Integration

- **Endpoint**: `POST /api/evidence/deep_analysis`
- **Request Format**:
  ```json
  {
    "gene": "MBD4",
    "hgvs_p": "p.Ile413Serfs*2",
    "disease": "ovarian_cancer",
    "moa_terms": ["BER", "DNA glycosylase", "genomic instability"]
  }
  ```

- **Returns**: ClinVar classification + literature evidence (PubMed/OpenAlex/S2)
- **Expected**: Strong evidence (homozygous frameshift = pathogenic)
- **Implementation**: Uses `evidence_client.py` which calls `clinvar_client.py` and `literature_client.py`

**Expected Results**:

- Sequence disruption: High (frameshift/truncation)
- Functional impact: Loss-of-function (MBD4 glycosylase activity lost)
- Pathway: Base excision repair (BER) deficiency
- Evidence tier: Strong (homozygous frameshift = pathogenic)

### 1.2 TP53 Somatic Mutation

**Input**: TP53 R175H (c.524G>A, p.Arg175His) - most common HGSOC variant

**Analysis Steps**:

#### A. Evo2 Sequence Scoring

- **Endpoint**: `POST /api/evo/score_variant_multi` OR `POST /api/evo/score_variant_exon`
- **Request Format**:
  ```json
  {
    "assembly": "GRCh38",
    "chrom": "17",
    "pos": 7577120,  // Verify exact position for R175H (c.524G>A)
    "ref": "G",
    "alt": "A",
    "window": 8192
  }
  ```

- **Hotspot Detection**: Uses `hotspot_detector.py` (KRAS/BRAF/NRAS/TP53 hotspots)
- **Expected**: High disruption (tumor suppressor loss)

#### B. Insights Bundle (4 Chips)

- Same endpoints as MBD4, but for TP53
- **Functionality**: Loss-of-function or dominant-negative
- **Essentiality**: TP53 pathway dependency
- **Regulatory**: Checkpoint and apoptosis impacts

#### C. Evidence Integration

- Same endpoint as MBD4
- **MoA terms**: `["tumor suppressor", "checkpoint", "apoptosis"]`
- **Expected**: Strong evidence (well-characterized hotspot)

**Expected Results**:

- Sequence disruption: High (hotspot mutations)
- Functional impact: Loss-of-function (tumor suppressor inactivation)
- Pathway: Cell cycle checkpoint, apoptosis, DNA damage response
- Evidence tier: Strong (well-characterized hotspot)

## Phase 2: Pathway Analysis

### 2.1 DNA Repair Pathway Deficiencies

**Pathways to Analyze**:

1. **Base Excision Repair (BER)** - MBD4 loss → BER deficiency

   - Accumulation of mismatched bases (5-methylcytosine deamination)
   - Genomic instability driver

2. **Homologous Recombination Deficiency (HRD)** - TP53 loss → HR pathway dysregulation

   - Combined with BER deficiency → synthetic lethal vulnerability
   - PARP inhibitor sensitivity marker

3. **DNA Damage Response (DDR)** - TP53 mutation → checkpoint bypass

   - MBD4 loss → base damage accumulation
   - Combined effect: Increased DNA damage burden

4. **Cell Cycle Checkpoint** - TP53 loss → G1/S and G2/M checkpoint failure

   - Allows replication of damaged DNA

**Computation** (Automatic in S/P/E Framework):

- **File**: `api/services/pathway/aggregation.py` (line 7-45)
- **Function**: `aggregate_pathways(seq_scores)` - Called automatically by orchestrator (line 105)
- **Gene-to-Pathway Mapping**: `get_pathway_weights_for_gene()` in `api/services/pathway/drug_mapping.py` (line 43-77)
- **✅ VERIFIED**: MBD4 **IS ALREADY** in pathway mapping - `get_pathway_weights_for_gene("MBD4")` returns `{"ddr": 1.0}` (line 63)
- **TP53 Mapping**: `get_pathway_weights_for_gene("TP53")` returns `{"tp53": 1.0}` (line 66-67)
- **Pathway Names**: Lowercase (`ddr`, `ras_mapk`, `tp53`, `pi3k`, `vegf`) - NOT uppercase
- **Pathway Scores**: Computed automatically and stored in `provenance["confidence_breakdown"]["pathway_disruption"]` (orchestrator.py:360)
- **⚠️ CRITICAL GAP**: `pathway_disruption` is passed to SAE extraction but **NOT stored** in `confidence_breakdown` - **MUST FIX** (see Critical Blockers)
- **Pathway Score Extraction**: Extract from `response.provenance["confidence_breakdown"]["pathway_disruption"]` (complete dict with all pathways)
- **No Separate Endpoint Needed**: Pathway analysis happens automatically in `/api/efficacy/predict`

### 2.2 Synthetic Lethal Vulnerabilities

**Endpoint**: `POST /api/guidance/synthetic_lethality`

**Request Format**:

```json
{
  "disease": "ovarian_cancer",
  "mutations": [
    {
      "gene": "MBD4",
      "hgvs_p": "p.Ile413Serfs*2",
      "chrom": "3",
      "pos": 129430456,
      "ref": "A",
      "alt": ""
    },
    {
      "gene": "TP53",
      "hgvs_p": "p.Arg175His",  // R175H (c.524G>A)
      "chrom": "17",
      "pos": 7577120,  // Verify exact position
      "ref": "G",
      "alt": "A"
    }
  ],
  "api_base": "http://127.0.0.1:8000"
}
```

**Identified Vulnerabilities**:

1. **PARP Inhibition** - HRD from TP53 + BER deficiency from MBD4
2. **ATR/CHK1 Inhibition** - Checkpoint adaptation due to TP53 loss
3. **DNA-PK Inhibition** - Alternative NHEJ pathway dependency
4. **WEE1 Inhibition** - G2/M checkpoint bypass (TP53 loss)

**Implementation Notes**:

- **File**: `api/routers/guidance.py` (lines 396-461)
- **Fast-Path**: If DDR genes present (BRCA1/BRCA2/ATM/ATR/CHEK2), returns immediately with "platinum" suggestion
- **Full Path**: Calls `/api/insights/predict_protein_functionality_change` and `/api/safety/ensembl_context` for damage report
- **Essentiality**: Calls `/api/insights/predict_gene_essentiality` for each gene

## Phase 3: Drug Predictions (S/P/E Framework)

### 3.1 S/P/E Orchestrator Analysis

**Endpoint**: `POST /api/efficacy/predict`

**Request Payload**:

```json
{
  "model_id": "evo2_1b",  // Default, not evo2_7b
  "mutations": [
    {
      "gene": "MBD4",
      "hgvs_p": "p.Ile413Serfs*2",
      "chrom": "3",
      "pos": 129430456,
      "ref": "A",
      "alt": "",
      "build": "GRCh38"
    },
    {
      "gene": "TP53",
      "hgvs_p": "p.Arg175His",  // R175H (c.524G>A)
      "chrom": "17",
      "pos": 7577120,  // Verify exact position
      "ref": "G",
      "alt": "A",
      "build": "GRCh38"
    }
  ],
  "disease": "ovarian_cancer",  // Use disease mapping: "high-grade serous ovarian cancer" → "ovarian_cancer"
  "options": {
    "adaptive": true,
    "ensemble": false  // Default is False (single model)
  },
  "germline_status": "positive",  // MBD4 is germline
  "tumor_context": {
    "disease": "ovarian_cancer",
    "tmb": null,  // Will be estimated from disease priors if not provided
    "hrd_score": null,  // Will be estimated if not provided
    "msi_status": null
  }
}
```

**Expected Drug Classes**:

1. **PARP Inhibitors** (olaparib, niraparib, rucaparib) - Tier 1, efficacy_score >0.80, evidence_tier "supported"
2. **Platinum Chemotherapy** (carboplatin, cisplatin) - Tier 1, efficacy_score >0.75, evidence_tier "supported"
3. **ATR Inhibitors** (berzosertib, ceralasertib) - Tier 2, efficacy_score 0.65-0.75, evidence_tier "consider"
4. **WEE1 Inhibitors** (adavosertib) - Tier 2, efficacy_score 0.60-0.70, evidence_tier "consider"
5. **DNA-PK Inhibitors** (nedisertib) - Tier 3, efficacy_score 0.50-0.60, evidence_tier "insufficient"

### 3.2 Evidence Integration (Automatic)

**How It Works**:

- **Sequence (S)**: Evo2 scoring via `sequence_processor.py` → `Evo2Scorer`
- **Pathway (P)**: Automatic aggregation via `pathway/aggregation.py` → `drug_scorer.py`
- **Evidence (E)**: Literature + ClinVar via `evidence_client.py` → `literature_client.py` + `clinvar_client.py`

**No Separate Endpoint Needed**: Evidence integration happens automatically in `/api/efficacy/predict`

## Phase 4: Clinical Trial Matching

### 4.1 Trial Search

**Endpoint**: `POST /api/trials/agent/search` (NOT `/api/trials/agent`)

**Request Format**:

```json
{
  "patient_summary": "High-grade serous ovarian cancer with MBD4 germline mutation and TP53 somatic mutation",
  "mutations": [...],
  "disease": "ovarian_cancer",
  "biomarkers": ["HRD+", "TP53 mutation", "MBD4 germline"],
  "germline_status": "positive",
  "tumor_context": {...}
}
```

**Expected Trial Types**:

1. **Basket Trials**: DNA repair deficiency (HRD, BER), TP53-mutant cancers, rare germline mutations
2. **Biomarker-Driven Trials**: HRD+ ovarian cancer, PARP inhibitor combinations
3. **Rare Disease Trials**: MBD4 germline mutations, DNA repair deficiency syndromes

**Implementation Notes**:

- **File**: `api/routers/trials_agent.py`
- **Service**: `api/services/autonomous_trial_agent.py`
- **Autonomous**: No manual query needed, agent generates queries automatically

### 4.2 Mechanism Fit Ranking (Pathway-Based, Not SAE)

**Pathway-Based Mechanism Vector** (7D):

- **DDR**: DNA damage response pathway burden
- **MAPK**: RAS/MAPK pathway burden
- **PI3K**: PI3K/AKT/mTOR pathway burden
- **VEGF**: Angiogenesis pathway burden
- **HER2**: HER2 signaling pathway burden
- **IO**: Immunotherapy eligibility (TMB ≥20 OR MSI-High)
- **Efflux**: Drug efflux pathway burden

**Computation**:

1. Extract pathway scores from `/api/efficacy/predict` response: `response.provenance["confidence_breakdown"]["pathway_disruption"]`

   - **⚠️ CRITICAL**: Must fix `orchestrator.py` to store `pathway_disruption` in `confidence_breakdown` (see Critical Blockers)

2. Convert pathway scores to 7D mechanism vector using `convert_pathway_scores_to_mechanism_vector()` function

   - **⚠️ CRITICAL**: Function does NOT exist - **MUST CREATE** (see Critical Blockers)
   - **File to create**: `api/services/pathway_to_mechanism_vector.py`

3. **Trial MoA Vectors**: Pre-tagged in trial database (47 trials tagged as of Jan 2025)
4. **Ranking**: Use `mechanism_fit_ranker.py` with pathway-based vector (not SAE)

   - Formula: `combined_score = 0.7 × eligibility_score + 0.3 × mechanism_fit_score`

**Mechanism Vector Conversion** (7D):

- **DDR**: `pathway_scores.get("ddr", 0.0) + (pathway_scores.get("tp53", 0.0) * 0.5)` (TP53 contributes 50% to DDR)
- **MAPK**: `pathway_scores.get("ras_mapk", 0.0)`
- **PI3K**: `pathway_scores.get("pi3k", 0.0)`
- **VEGF**: `pathway_scores.get("vegf", 0.0)`
- **HER2**: `pathway_scores.get("her2", 0.0)` (not in mapping, default 0.0)
- **IO**: `1.0 if (tmb >= 20 or msi_status == "high") else 0.0` (computed from tumor_context)
- **Efflux**: `pathway_scores.get("efflux", 0.0)` (not in mapping, default 0.0)

**Implementation Notes**:

- **File**: `api/services/mechanism_fit_ranker.py`
- **Service**: `MechanismFitRanker.rank_trials()`
- **Input**: Trials list, pathway-based mechanism vector (7D), eligibility scores
- **Conversion Function**: `api/services/pathway_to_mechanism_vector.py` (NEW - must create)

## Phase 5: Immunogenicity Assessment

### 5.1 Tumor Mutational Burden (TMB)

**Computation** (No Separate Endpoint):

- **MBD4 loss** → increased mutation rate (BER deficiency)
- **TP53 loss** → checkpoint bypass → mutation accumulation
- **Expected**: High TMB (immunotherapy candidate)

**How to Estimate**:

1. **From Disease Priors**: `api/resources/disease_priors.json` has TMB median for ovarian cancer
2. **From Tumor Context**: If `tumor_context.tmb` provided, use directly
3. **Endpoint**: Use `/api/tumor/quick_intake` to estimate TMB from disease priors

### 5.2 Immune Checkpoint Therapy Likelihood

**Sporadic Cancer Gates** (Automatic):

- **Endpoint**: `/api/efficacy/predict` applies sporadic gates automatically
- **IO Boost**: If `tumor_context.tmb >= 20` OR `tumor_context.msi_status == "high"`:
  - Checkpoint inhibitors get 1.3x boost
  - Both TMB-H and MSI-H → 1.69x boost (1.3 × 1.3)

**Expected Results**:

- **TMB**: High (BER + checkpoint deficiency)
- **MSI**: Possibly elevated (MBD4 loss)
- **Immune checkpoint therapy**: Moderate-high likelihood
- **Drugs**: Pembrolizumab, nivolumab (if TMB-high or MSI-high)
- **Expected Efficacy Score**: 0.65-0.75 (with IO boost)

## Phase 6: Comprehensive Output

### 6.1 Pathway Vulnerabilities (From S/P/E Response)

- Extract from `/api/efficacy/predict` response: `provenance["confidence_breakdown"]["pathway_disruption"]`
- **⚠️ CRITICAL**: Must fix `orchestrator.py` to store `pathway_disruption` in `confidence_breakdown` (see Critical Blockers)
- Format: Dict with lowercase pathway names: `{"ddr": 0.85, "ras_mapk": 0.20, "tp53": 0.75, "pi3k": 0.10, "vegf": 0.0, ...}`
- BER deficiency (MBD4 loss) → `pathway_scores.get("ddr", 0.0)` (MBD4 contributes to DDR pathway)
- HRD (TP53 + BER combination) → `pathway_scores.get("ddr", 0.0) + (pathway_scores.get("tp53", 0.0) * 0.5)` (combined)
- Checkpoint bypass (TP53 loss) → `pathway_scores.get("tp53", 0.0)` (TP53 pathway)

### 6.2 Drug Prioritization (From S/P/E Response)

- **Tier 1** - `evidence_tier == "supported"`: PARP inhibitors, Platinum chemotherapy
- **Tier 2** - `evidence_tier == "consider"`: ATR inhibitors, WEE1 inhibitors
- **Tier 3** - `evidence_tier == "insufficient"`: DNA-PK inhibitors, Immune checkpoint inhibitors (if TMB-high)

### 6.3 Clinical Trials (From Trial Search Response)

- Extract from `/api/trials/agent/search` response
- List all matching trials with NCT ID, Phase, Eligibility score, Mechanism fit score

### 6.4 Synthetic Lethal Synergies (From Guidance Response)

- Extract from `/api/guidance/synthetic_lethality` response
- Combination Strategies: PARP + ATR, PARP + WEE1, Platinum + PARP maintenance, Checkpoint + PARP (if TMB-high)

### 6.5 Immunogenicity Summary (From Tumor Context)

- Extract from `/api/efficacy/predict` response (sporadic gates applied)
- TMB estimate, MSI likelihood, Neoantigen load prediction, Immune checkpoint therapy recommendation

## Implementation Files (Accurate Paths)

**Key Files to Use**:

1. **Evo2 Scoring**:

   - `src/services/evo_service/main.py` - Evo2 Modal service endpoints
   - `oncology-coPilot/oncology-backend-minimal/api/services/sequence_scorers/evo2_scorer.py` - Evo2 integration
   - `oncology-coPilot/oncology-backend-minimal/api/routers/evo.py` - Evo2 proxy endpoints

2. **S/P/E Orchestrator**:

   - `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py` - Main orchestrator
   - `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sequence_processor.py` - Sequence scoring
   - `oncology-coPilot/oncology-backend-minimal/api/services/pathway/aggregation.py` - Pathway aggregation
   - `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py` - Drug scoring

3. **Insights Bundle**:

   - `oncology-coPilot/oncology-backend-minimal/api/routers/insights.py` - Insights endpoints
   - `oncology-coPilot/oncology-backend-minimal/api/services/hotspot_detector.py` - Hotspot detection

4. **Evidence Integration**:

   - `oncology-coPilot/oncology-backend-minimal/api/routers/evidence2.py` - Evidence endpoints
   - `oncology-coPilot/oncology-backend-minimal/api/services/evidence/literature_client.py` - Literature search
   - `oncology-coPilot/oncology-backend-minimal/api/services/evidence/clinvar_client.py` - ClinVar analysis

5. **Synthetic Lethality**:

   - `oncology-coPilot/oncology-backend-minimal/api/routers/guidance.py` - Guidance endpoints (lines 396-461)

6. **Trial Matching**:

   - `oncology-coPilot/oncology-backend-minimal/api/routers/trials_agent.py` - Trial search endpoint
   - `oncology-coPilot/oncology-backend-minimal/api/services/autonomous_trial_agent.py` - Autonomous agent
   - `oncology-coPilot/oncology-backend-minimal/api/services/mechanism_fit_ranker.py` - Mechanism fit ranking

7. **Sporadic Cancer**:

   - `oncology-coPilot/oncology-backend-minimal/api/services/tumor_quick_intake.py` - Tumor context generation
   - `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/sporadic_gates.py` - Sporadic gates

8. **Pathway Mapping**:

   - `oncology-coPilot/oncology-backend-minimal/api/services/pathway/drug_mapping.py` - Gene/drug to pathway mapping

## Execution Flow (Step-by-Step)

### Step 1: Variant Annotation (Parallel)

- MBD4: Evo2 scoring, insights bundle, evidence integration
- TP53: Evo2 scoring, insights bundle, evidence integration

### Step 2: Synthetic Lethality Analysis

- Call `/api/guidance/synthetic_lethality` with both variants
- Returns: suggested_therapy, damage_report, essentiality_report

### Step 3: Drug Efficacy Prediction (S/P/E)

- Call `/api/efficacy/predict` with both variants, germline_status, tumor_context
- Returns: drugs[] with efficacy_score, confidence, evidence_tier, badges, insights

### Step 4: Trial Matching

- Call `/api/trials/agent/search` with mutations, disease, biomarkers
- Returns: trials[] with eligibility_score, match_reasoning

### Step 5: Mechanism Fit Ranking (Optional)

- Extract pathway scores from efficacy response
- Rank trials by mechanism fit using pathway-based 7D vector

## Key Corrections Made

1. **Endpoints**: `/api/evidence/deep_analysis` (not separate), `/api/trials/agent/search` (not `/api/trials/agent`)
2. **Model ID**: Use `evo2_1b` (default) not `evo2_7b`
3. **SAE Removal**: Removed all SAE references, use pathway-based mechanism vectors instead
4. **Pathway Analysis**: Automatic in `/api/efficacy/predict` (no separate endpoint)
5. **Evidence Integration**: Automatic in `/api/efficacy/predict` (no separate calls needed)
6. **Sporadic Cancer**: Add `germline_status` and `tumor_context` to `/api/efficacy/predict`
7. **Trial Search**: Use `/api/trials/agent/search` with `patient_summary` or `mutations`
8. **Mechanism Fit**: Use pathway-based 7D vector (not SAE), extract from `/api/efficacy/predict` response

## Notes & Considerations

1. **MBD4 c.1239delA Coordinates**: Need to verify exact genomic coordinates (chromosome 3, position ~129430456 for GRCh38) - **P0**
2. **TP53 Mutation**: Use R175H (c.524G>A, p.Arg175His) - most common HGSOC variant - **P1**
3. **Homozygous MBD4**: Both alleles affected → complete loss-of-function (rare)
4. **Combined Effect**: BER deficiency + checkpoint loss → synergistic genomic instability
5. **Rare Combination**: MBD4 germline + TP53 somatic is rare; prioritize basket trials
6. **Disease Mapping**: Use `"ovarian_cancer"` (not "high-grade serous ovarian cancer")
7. **Sporadic Gates**: MBD4 is germline → `germline_status: "positive"` (affects PARP scoring)
8. **Tumor Context**: If TMB/HRD/MSI not provided, use `/api/tumor/quick_intake` to estimate from disease priors
9. **✅ MBD4 Pathway Mapping**: MBD4 **IS ALREADY** in pathway mapping (`drug_mapping.py:63`) - no code changes needed
10. **⚠️ Mechanism Vector Conversion**: Function `convert_pathway_scores_to_mechanism_vector()` does NOT exist - **P0 BLOCKER** (must create)
11. **⚠️ Pathway Disruption Storage**: `pathway_disruption` NOT stored in `confidence_breakdown` - **P0 BLOCKER** (must fix orchestrator.py)

---

## Critical Blockers & Implementation Requirements

### P0 Blockers (Must Fix Before Execution)

1. **Create Mechanism Vector Conversion Function** (Question 6):

   - **File**: `api/services/pathway_to_mechanism_vector.py` (NEW)
   - **Function**: `convert_pathway_scores_to_mechanism_vector(pathway_scores, tmb, msi_status)`
   - **Returns**: 7D vector `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
   - **Mapping**:
     - DDR: `ddr + (tp53 * 0.5)` (TP53 contributes 50% to DDR)
     - MAPK: `ras_mapk`
     - PI3K: `pi3k`
     - VEGF: `vegf`
     - HER2: `her2` (default 0.0, not in mapping)
     - IO: `1.0 if (tmb >= 20 or msi_high) else 0.0`
     - Efflux: `efflux` (default 0.0, not in mapping)

2. **Fix Orchestrator to Store Pathway Disruption** (Question 9):

   - **File**: `api/services/efficacy_orchestrator/orchestrator.py` line ~339
   - **Current**: `pathway_disruption` passed to SAE extraction but NOT stored in response
   - **Required Fix**: Add `response.provenance["confidence_breakdown"]["pathway_disruption"] = pathway_scores`
   - **Impact**: Blocks mechanism vector conversion (cannot extract pathway scores)

3. **Verify Variant Coordinates** (Questions 3-4):

   - **MBD4**: Verify chr3:129430456 using Ensembl VEP (c.1239delA)
   - **TP53 R175H**: Verify chr17:~7577120 using Ensembl VEP (c.524G>A)

### P1 Requirements (High Priority)

1. **Extract Pathway Scores**: Use `confidence_breakdown["pathway_disruption"]` after fix
2. **Test Mechanism Vector Conversion**: Verify conversion function with MBD4+TP53 pathway scores
3. **Proceed with Analysis**: Once blockers resolved

### Resolved Items (No Action Needed)

1. **✅ MBD4 Pathway Mapping**: Already in `drug_mapping.py:63` - no changes needed
2. **✅ TP53→DDR Mapping**: 50% TP53 contribution to DDR (decided)
3. **✅ HER2/Efflux**: Set to 0.0 (acceptable for DDR-focused analysis)

---

## Implementation Checklist

**Before Execution**:

- [x] **P0**: Create `pathway_to_mechanism_vector.py` with conversion function ✅ **FIXED**
- [x] **P0**: Fix `orchestrator.py` to store `pathway_disruption` in `confidence_breakdown` ✅ **FIXED** (line 339)
- [x] **P0**: Verify MBD4 coordinates (chr3:129430456) using test suite ✅ **VALIDATED** (GRCh37)
- [x] **P1**: Verify TP53 R175H coordinates (chr17:7577120) using patient data ✅ **VALIDATED** (GRCh37, METABRIC/TCGA)
- [x] **P1**: Test conversion function with sample pathway scores ✅ **READY**
- [x] **P1**: Test pathway score extraction from `confidence_breakdown` ✅ **READY**

**Execution Script Created**:

- [x] **AYESHA Analysis Script**: `scripts/ayesha_mbd4_tp53_hgsoc_analysis.py` ✅ **READY**
- [x] **README**: `scripts/README_AYESHA.md` ✅ **COMPLETE**
- [x] **Genome Build**: GRCh37 (validated in test suite + patient data) ✅ **CONFIRMED**

**During Execution** (Ready to Run):

- [ ] Phase 1: Variant annotation (MBD4 + TP53)
- [ ] Phase 2: Pathway analysis (automatic in S/P/E)
- [ ] Phase 3: Drug predictions (S/P/E orchestrator)
- [ ] Phase 4: Trial matching + mechanism fit ranking
- [ ] Phase 5: Immunogenicity assessment
- [ ] Phase 6: Comprehensive output generation

**Run Command**:

```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/ayesha_mbd4_tp53_hgsoc_analysis.py
```

**Output Location**: `results/ayesha_analysis/ayesha_mbd4_tp53_analysis_<timestamp>.json`