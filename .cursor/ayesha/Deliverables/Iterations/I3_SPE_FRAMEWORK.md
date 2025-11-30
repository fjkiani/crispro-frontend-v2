# 📊 ITERATION 3: S/P/E FRAMEWORK & EFFICACY SYSTEM

**Status**: ✅ **COMPLETE**  
**Duration**: 3-4 hours  
**Created**: January 14, 2025

---

### **3.1 THE CORE S/P/E FORMULA**

#### **Exact Formula** (`drug_scorer.py:171`):
```python
efficacy_score = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * evidence + clinvar_prior
```

**⚠️ IMPORTANT DISCREPANCY**: 
- **Config defaults** (`config.py:20-22`): `0.35 / 0.35 / 0.30` (Sequence/Pathway/Evidence)
- **Actual hardcoded** (`drug_scorer.py:171`): `0.3 / 0.4 / 0.3`
- **Current implementation uses hardcoded values** (0.3/0.4/0.3), NOT config values
- Config values exist but are not currently used in the formula

**Components**:
- **seq_pct**: Calibrated sequence percentile [0, 1]
- **path_pct**: Normalized pathway percentile [0, 1]
- **evidence**: Evidence strength score [0, 1]
- **clinvar_prior**: ClinVar prior boost [-0.2, +0.2]

**Insufficient Tier Penalty** (`drug_scorer.py:172`):
```python
lob = raw_lob if tier != "insufficient" else (raw_lob * 0.5 if confidence_config.fusion_active else 0.0)
```
- If tier is "insufficient": Penalize by 50% (if Fusion active) or set to 0.0 (if Fusion inactive)

---

### **3.2 SEQUENCE (S) SCORING - COMPLETE DEEP DIVE**

#### **3.2.1 Scorer Fallback Chain**
**Priority Order** (`sequence_processor.py:29-86`):
1. **Fusion Engine** (AlphaMissense) - ONLY for GRCh38 missense variants
2. **Evo2 Adaptive** - Multi-window, ensemble models
3. **Massive Oracle** - Synthetic or real-context (if enabled)

**Fusion Gate** (`sequence_processor.py:32-40`):
- Only eligible if: `build == "GRCh38"` AND `consequence == "missense"`
- If not eligible, skips Fusion and goes to Evo2

#### **3.2.2 Evo2 Scoring - Complete Flow**

**Step 1: Model Selection** (`evo2_scorer.py:95-109`):
- Priority 1: `EVO_FORCE_MODEL` env var (forces single model)
- Priority 2: `EVO_ALLOWED_MODELS` env var (filters allowed models)
- Priority 3: Default candidates: `["evo2_1b", "evo2_7b", "evo2_40b"]` if ensemble, else `[model_id]`
- **Best Model Selection**: Picks model with strongest `abs(min_delta)`

**Step 2: Multi-Window Strategy** (`evo2_scorer.py:241-267`):
- Default windows: `[4096, 8192, 16384, 25000]` bp
- Tests each window size, picks best `exon_delta`
- **Spam-safety**: `EVO_MAX_FLANKS` limits number of windows tested

**Step 3: Forward/Reverse Symmetry** (`evo2_scorer.py:276-332`):
- **Forward**: `ref → alt` (standard variant)
- **Reverse**: `alt → ref` (symmetry check)
- **Averaging**: `avg_min_delta = (forward_min + reverse_min) / 2.0`
- **Feature Flag**: `evo_disable_symmetry` (default: True, symmetry disabled)
- **Why**: Evo2 is strand-agnostic; averaging reduces noise

**Step 4: Sequence Disruption Calculation** (`evo2_scorer.py:126-136`):
```python
sequence_disruption = max(abs(min_delta), abs(exon_delta))
```
- Uses **stronger signal** between multi-window (`min_delta`) and exon-context (`exon_delta`)
- Both are averaged forward/reverse if symmetry enabled

**Step 5: Special Handling** (`evo2_scorer.py:138-157`):
- **Hotspot Floors** (raw delta level): `HOTSPOT_FLOOR = 1e-4`
  - BRAF V600 → `max(disruption, 1e-4)`
  - KRAS/NRAS/HRAS G12/G13/Q61 → `max(disruption, 1e-4)`
  - TP53 R175/R248/R273 → `max(disruption, 1e-4)`
- **Truncation/Frameshift Lift**: If `"*" in hgvs_p` or `"FS" in hgvs_p` → `max(disruption, 1.0)`
  - Stop codons and frameshifts are highly disruptive → enforce maximum

**Step 6: Percentile Calibration** (`evo2_scorer.py:159-171`):
- **Primary Method**: `percentile_like(sequence_disruption)` - Piecewise mapping
- **Hotspot Floors** (percentile level):
  - BRAF V600 → `max(pct, 0.90)`
  - KRAS/NRAS/HRAS G12/G13/Q61 → `max(pct, 0.80)`
  - TP53 R175/R248/R273 → `max(pct, 0.80)`

**Step 7: Gene-Specific Calibration** (`gene_calibration.py`):
- **Service Exists**: `GeneCalibrationService` can compute gene-specific percentiles from ClinVar
- **NOT CURRENTLY USED**: Main flow uses `percentile_like()` (simple piecewise)
- **Future Enhancement**: Could replace `percentile_like()` with gene-specific calibration
- **Method**: Fetches ClinVar variants for gene, scores with Evo2, builds distribution, computes percentiles

#### **3.2.3 Percentile Mapping Function** (`utils.py:7-23`):
```python
def percentile_like(value: float) -> float:
    if v <= 0.005: return 0.05
    if v <= 0.01:  return 0.10
    if v <= 0.02:  return 0.20
    if v <= 0.05:  return 0.50
    if v <= 0.10:  return 0.80
    return 1.0
```
- **Conservative mapping**: Small deltas → low percentiles (avoids false positives)
- **Empirical thresholds**: Based on observed pathogenic ranges from ClinVar
- **Gene-agnostic**: Same mapping for all genes (limitation)

#### **3.2.4 Delta-Only Mode** (`evo2_scorer.py:224-234`):
- **Feature Flag**: `evo_use_delta_only`
- **When enabled**: Skips exon scanning entirely, uses only `min_delta` from multi-window
- **Why**: Faster execution, lower cost (fewer API calls)

#### **3.2.5 Fusion Scorer** (`fusion_scorer.py`):
- **Gate**: ONLY for GRCh38 missense variants
- **Method**: Calls AlphaMissense Fusion Engine API
- **Fallback**: If external Fusion unavailable, tries local `/api/fusion/score_variant`
- **Caching**: Redis-based caching with TTL

#### **3.2.6 Massive Oracle Scorer** (`massive_scorer.py`):
- **Two Modes**:
  - **Synthetic**: Contrasting sequences for massive impact detection
  - **Real-Context**: Real GRCh38 context with configurable flank size
- **Feature Flag**: `enable_massive_modes` (default: false)
- **Why**: Research tool for detecting very large disruptions

---

### **3.3 PATHWAY (P) SCORING - COMPLETE DEEP DIVE**

#### **3.3.1 Pathway Aggregation** (`aggregation.py:7-45`):
```python
def aggregate_pathways(seq_scores: List[Dict[str, Any]]) -> Dict[str, float]:
    pathway_totals = {}
    pathway_counts = {}
    
    for score in seq_scores:
        pathway_weights = score.get("pathway_weights", {})
        sequence_disruption = float(score.get("sequence_disruption", 0.0))
        
        for pathway, weight in pathway_weights.items():
            pathway_totals[pathway] += sequence_disruption * weight
            pathway_counts[pathway] += 1
    
    # Compute average scores
    pathway_scores = {}
    for pathway in pathway_totals:
        pathway_scores[pathway] = pathway_totals[pathway] / pathway_counts[pathway]
    
    return pathway_scores
```

**Formula**: `pathway_score = sum(sequence_disruption * gene_weight) / count`

#### **3.3.2 Gene-to-Pathway Weights** (`drug_mapping.py:16-32`):
```python
def get_pathway_weights_for_gene(gene_symbol: str) -> Dict[str, float]:
    g = gene_symbol.strip().upper()
    if g in {"BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1"}:
        return {"ras_mapk": 1.0}
    if g in {"TP53", "MDM2", "ATM", "ATR", "CHEK2", "BRCA1", "BRCA2", "PTEN", "RAD51"}:
        return {"tp53": 1.0}  # DNA repair bucket
    return {}
```

**Current State**: Hardcoded, simplified mapping
- **MAPK drivers** → `ras_mapk` pathway
- **DNA repair genes** → `tp53` pathway (DNA repair bucket)

#### **3.3.3 Drug-to-Pathway Weights** (`panel_config.py:8-14`):
```python
DEFAULT_MM_PANEL = [
    {"name": "BRAF inhibitor", "pathway_weights": {"ras_mapk": 0.8, "tp53": 0.2}},
    {"name": "MEK inhibitor", "pathway_weights": {"ras_mapk": 0.9, "tp53": 0.1}},
    {"name": "IMiD", "pathway_weights": {"ras_mapk": 0.2, "tp53": 0.3}},
    {"name": "Proteasome inhibitor", "pathway_weights": {"ras_mapk": 0.3, "tp53": 0.4}},
    {"name": "Anti-CD38", "pathway_weights": {"ras_mapk": 0.1, "tp53": 0.1}},
]
```

**Current State**: Hardcoded, Multiple Myeloma panel only

#### **3.3.4 Pathway Score Calculation** (`drug_scorer.py:44-51`):
```python
# Step 1: Get drug's pathway weights
drug_weights = get_pathway_weights_for_drug(drug_name)

# Step 2: Weighted sum of pathway scores
s_path = sum(pathway_scores.get(pathway, 0.0) * weight 
             for pathway, weight in drug_weights.items())

# Step 3: Normalize to [0, 1] based on empirical Evo2 ranges
if s_path > 0:
    path_pct = min(1.0, max(0.0, (s_path - 1e-6) / (1e-4 - 1e-6)))
else:
    path_pct = 0.0
```

**Normalization Formula**: 
- **Empirical Range**: Pathogenic deltas ~`1e-6` to `1e-4`
- **Linear mapping**: `(s_path - 1e-6) / (1e-4 - 1e-6)`
- **Clamp**: `[0, 1]`

---

### **3.4 EVIDENCE (E) SCORING - COMPLETE DEEP DIVE**

#### **3.4.1 Literature Evidence** (`literature_client.py:51-135`):

**Query Flow**:
1. **PubMed E-utils API** (`esearch.fcgi`): Search for PMIDs
2. **PubMed E-utils API** (`esummary.fcgi` or `efetch.fcgi`): Get paper details
3. **MoA Filtering**: Prefer papers mentioning drug name or MoA in title/abstract
4. **Strength Calculation**: Based on publication types

**Strength Scoring** (`literature_client.py:18-48`):
```python
def _score_evidence_from_results(top_results: List[Dict[str, Any]]) -> float:
    score = 0.0
    for r in top_results[:3]:
        pub_types = " ".join([_safe_lower(t) for t in (r.get("publication_types") or [])])
        title = _safe_lower(r.get("title"))
        
        if "randomized" in pub_types or "randomized" in title:
            score += 0.5  # RCT
        elif "guideline" in pub_types or "practice" in title:
            score += 0.35  # Guideline
        elif "review" in pub_types or "meta" in title:
            score += 0.25  # Review
        else:
            score += 0.15  # Other
    
    return float(min(1.0, score))
```

**MoA Boost** (`literature_client.py:110-120`):
- **Base strength**: From publication types
- **MoA hits**: Count papers mentioning MoA in title/abstract
- **Final strength**: `min(1.0, base_strength + 0.10 * moa_hits)`
- **Max boost**: +0.30 (3 MoA hits)

#### **3.4.2 ClinVar Prior** (`clinvar_client.py:11-99`):

**Prior Calculation** (`clinvar_client.py:54-59`):
```python
if cls in ("pathogenic", "likely_pathogenic"):
    prior = 0.2 if strong else (0.1 if moderate else 0.05)
elif cls in ("benign", "likely_benign"):
    prior = -0.2 if strong else (-0.1 if moderate else -0.05)
```

**Review Status Tiers**:
- **Strong**: `"expert" in review` OR `"practice" in review` → ±0.2
- **Moderate**: `"criteria" in review` → ±0.1
- **Weak**: Other → ±0.05

**Research-Mode Fallback** (`clinvar_client.py:68-94`):
- **Feature Flag**: `RESEARCH_USE_CLINVAR_CANONICAL`
- **Canonical Hotspots**: BRAF V600E/V600K, KRAS G12D/G12V, NRAS Q61K
- **If ClinVar empty**: Assigns `prior = 0.2` for canonical hotspots

#### **3.4.3 Evidence Tier Computation** (`tier_computation.py:9-68`):

**Legacy Tier Logic** (default):
```python
# Evidence gate: strong evidence OR ClinVar-Strong + pathway alignment
evidence_gate = (
    s_evd >= config.evidence_gate_threshold or  # Default: 0.7
    ("ClinVar-Strong" in badges and s_path >= config.pathway_alignment_threshold)  # Default: 0.2
)

# Insufficient signal: low sequence, pathway, and evidence
insufficient = (
    s_seq < config.insufficient_signal_threshold and  # Default: 0.02
    s_path < 0.05 and 
    s_evd < 0.2
)

if evidence_gate:
    return "supported"
elif insufficient:
    return "insufficient"
else:
    return "consider"
```

**V2 Tier Logic** (`CONFIDENCE_V2=1`):
```python
# Tier I (supported): FDA on‑label OR ≥1 RCT OR (ClinVar‑Strong AND pathway_aligned)
if ("FDA-OnLabel" in badges or 
    "RCT" in badges or 
    ("ClinVar-Strong" in badges and s_path >= 0.2)):
    return "supported"

# Tier II (consider): ≥2 human studies MoA‑aligned OR 1 strong study + pathway_aligned
if (s_evd >= 0.6 or  # Strong evidence (proxy for strong study)
    (s_evd >= 0.4 and s_path >= 0.2)):  # Moderate evidence + pathway alignment
    return "consider"

# Tier III (insufficient): else
return "insufficient"
```

**Default Thresholds** (`config.py:25-28`):
- `EVIDENCE_GATE_THRESHOLD = 0.7`
- `PATHWAY_ALIGNMENT_THRESHOLD = 0.2`
- `INSUFFICIENT_SIGNAL_THRESHOLD = 0.02`

---

### **3.5 CONFIDENCE COMPUTATION - COMPLETE DEEP DIVE**

#### **3.5.1 Legacy Confidence** (`confidence_computation.py:69-108`):

**Tier-Based Approach**:
```python
if tier == "supported":
    confidence = 0.6 + 0.2 * max(seq_pct, path_pct)
elif tier == "consider":
    if config.fusion_active and max(seq_pct, path_pct) >= 0.7:
        confidence = 0.5 + 0.2 * max(seq_pct, path_pct)
    else:
        confidence = 0.3 + 0.1 * seq_pct + 0.1 * path_pct
else:  # insufficient
    max_sp = max(seq_pct, path_pct)
    min_sp = min(seq_pct, path_pct)
    base = 0.20 + 0.35 * max_sp + 0.15 * min_sp
    if config.fusion_active:
        confidence = max(0.25, base)
    else:
        confidence = base

# Insights modulation
confidence += 0.05 if func >= 0.6 else 0.0
confidence += 0.04 if chrom >= 0.5 else 0.0
confidence += 0.07 if ess >= 0.7 else 0.0
confidence += 0.02 if reg >= 0.6 else 0.0

# Alignment margin boost
margin = abs(seq_pct - path_pct)
if margin >= 0.2:
    confidence += 0.05
```

#### **3.5.2 V2 Confidence** (`CONFIDENCE_V2=1`) (`confidence_computation.py:111-160`):

**Linear S/P/E Formula**:
```python
# Convert tier to evidence score (E component)
if tier == "supported":
    e_score = 0.05
elif tier == "consider":
    e_score = 0.02
else:  # insufficient
    e_score = 0.00

# Calculate lifts
lifts = 0.0
lifts += 0.04 if func >= 0.6 else 0.0      # Functionality
lifts += 0.02 if chrom >= 0.5 else 0.0     # Chromatin
lifts += 0.02 if ess >= 0.7 else 0.0       # Essentiality
lifts += 0.02 if reg >= 0.6 else 0.0      # Regulatory
lifts = min(lifts, 0.08)  # Cap total lifts at +0.08

# Linear S/P/E formula
confidence = 0.5 * seq_pct + 0.2 * path_pct + 0.3 * e_score + lifts
confidence = clamp01(confidence)  # Clamp to [0, 1]
return round(confidence, 2)
```

**Formula**: `confidence = 0.5·S + 0.2·P + 0.3·E + lifts`
- **S weight**: 0.5 (sequence percentile)
- **P weight**: 0.2 (pathway percentile)
- **E weight**: 0.3 (evidence tier score: 0.05/0.02/0.00)
- **Lifts**: Max +0.08 from insights

#### **3.5.3 Additional Confidence Boosts** (`drug_scorer.py:140-168`):

**ClinVar Prior Boost**:
```python
if clinvar_prior > 0 and path_pct >= 0.2:
    confidence += min(0.1, clinvar_prior)
```

**Gene-Drug MoA Tie-Breaker**:
```python
if seq_scores and path_pct >= 0.2:
    primary_gene = (seq_scores[0].variant or {}).get("gene", "").upper()
    if primary_gene == "BRAF" and drug_name == "BRAF inhibitor":
        confidence += 0.01  # Target engagement bonus
    elif primary_gene in {"KRAS", "NRAS"} and drug_name == "MEK inhibitor":
        confidence += 0.01  # Downstream effector bonus
```

**Research-Mode Pathway Prior** (`RESEARCH_USE_PATHWAY_PRIOR=1`):
```python
if primary_gene in {"KRAS", "NRAS"} and drug_name == "MEK inhibitor":
    confidence += 0.02
if primary_gene == "BRAF" and drug_name == "BRAF inhibitor":
    confidence += 0.02
```

---

### **3.6 SPORADIC CANCER GATES - COMPLETE INTEGRATION**

#### **3.6.1 Gate 1: PARP Inhibitor Penalty** (`sporadic_gates.py:56-127`):

**Logic**:
- **Germline positive** → Full effect (1.0x)
- **Germline negative + HRD ≥42** → **RESCUED** (1.0x) ⚔️
- **Germline negative + HRD <42** → Reduced (0.6x)
- **Unknown germline + unknown HRD** → Conservative (0.8x)

**Formula**: `efficacy_score *= parp_penalty`

#### **3.6.2 Gate 2: Immunotherapy Boost** (`sporadic_gates.py:128-196`):

**Priority Order** (mutually exclusive, highest wins):
1. **TMB ≥20** → 1.35x boost (highest priority)
2. **MSI-High** → 1.30x boost (second priority)
3. **TMB ≥10 but <20** → 1.25x boost (lowest priority)

**Formula**: `efficacy_score *= io_boost_factor` (single factor, not multiplicative)

#### **3.6.3 Gate 3: Confidence Capping** (`sporadic_gates.py:197-232`):

**By Completeness Level**:
- **Level 0** (completeness <0.3): Cap at 0.4
- **Level 1** (0.3 ≤ completeness <0.7): Cap at 0.6
- **Level 2** (completeness ≥0.7): No cap

**Formula**: `confidence = min(confidence, cap)`

---

### **3.7 ABLATION MODES**

#### **Ablation Support** (`orchestrator.py:188-198`):
```python
ablation = (request.ablation_mode or ("SP" if fast_mode else "SPE")).upper()
use_S = "S" in ablation
use_P = "P" in ablation
use_E = "E" in ablation

# Shallow copies/masks
masked_seq_scores = seq_scores if use_S else []
masked_pathway_scores = pathway_scores if use_P else {}
masked_evidence = (evidence_results[i] if (use_E and i < len(evidence_results)) else None)
```

**Modes**:
- **"SPE"**: Full S/P/E (default in normal mode)
- **"SP"**: Sequence + Pathway only (default in fast mode)
- **"S"**: Sequence only
- **"P"**: Pathway only
- **"E"**: Evidence only

**Why**: Allows testing individual component contributions

---

### **3.8 KEY INSIGHTS & EDGE CASES**

#### **Sequence Scoring**:
1. **Hotspot floors applied TWICE**: Once at raw delta (1e-4), once at percentile (0.80-0.90)
2. **Truncation/frameshift**: Enforces `disruption = 1.0` for stop codons/frameshifts
3. **Forward/reverse symmetry**: Can be disabled via `evo_disable_symmetry` (default: disabled)
4. **Gene-specific calibration exists but not used**: `GeneCalibrationService` available but main flow uses `percentile_like()`
5. **Delta-only mode**: `evo_use_delta_only` skips exon scanning for speed

#### **Pathway Scoring**:
1. **Hardcoded weights**: Gene-to-pathway and drug-to-pathway weights are hardcoded
2. **Simplified mapping**: Only 2 pathways (ras_mapk, tp53) currently
3. **Normalization**: Based on empirical Evo2 ranges (1e-6 to 1e-4)

#### **Evidence Scoring**:
1. **MoA boost**: +0.10 per MoA hit, max +0.30
2. **Research-mode fallback**: Canonical hotspots get `prior = 0.2` if ClinVar empty
3. **Timeout handling**: Evidence timeout → `tier = "insufficient"`, `s_evd = 0.0`

#### **Confidence Computation**:
1. **Two versions**: Legacy (tier-based) and V2 (linear S/P/E)
2. **Feature flag**: `CONFIDENCE_V2=1` enables V2
3. **Lifts capped**: Total insights lifts capped at +0.08
4. **ClinVar boost**: +0.1 max when aligned (path_pct ≥ 0.2)

#### **Sporadic Gates**:
1. **PARP rescue**: HRD ≥42 rescues PARP for germline-negative patients
2. **IO boost priority**: TMB ≥20 > MSI-H > TMB ≥10 (mutually exclusive)
3. **Confidence capping**: Based on tumor context completeness (L0/L1/L2)

---

## **3.9 S/P/E FRAMEWORK ARCHITECTURE - COMPLETE OVERVIEW** (Cycle 1 Deep Dive)

### **3.9.1 Master Doctrine Summary**

**Core Mission**: Provide a unified framework for variant analysis that combines sequence-level signals, pathway-level aggregation, and evidence-based validation to deliver transparent, auditable therapeutic recommendations.

**Core Principles**:
- **Multi-Modal Integration**: Never rely on single metrics
- **Transparent Confidence**: All rationale and provenance tracked
- **Evidence-Based**: Literature + ClinVar + pathway alignment
- **Research Use Only**: All outputs clearly labeled RUO

### **3.9.2 Complete Data Flow**

```
User Input (mutations + options)
  ↓
EfficacyOrchestrator.predict()
  ↓
[1] SequenceProcessor.score_sequences()
    ├─ FusionAMScorer (GRCh38 missense only) → AlphaMissense scores
    ├─ Evo2Scorer (default) → Evo2 delta scores (multi-window, adaptive)
    └─ MassiveOracleScorer (if enabled) → Synthetic/real-context scores
  ↓
[2] Pathway Aggregation (aggregate_pathways)
    ├─ Gene→pathway weights (get_pathway_weights_for_gene)
    ├─ Weighted aggregation: pathway_score = sum(seq_disruption * weight) / count
    └─ Pathway scores: {ras_mapk: 0.3, tp53: 0.2, ...}
  ↓
[3] Evidence Gathering (parallel execution)
    ├─ literature() → PubMed E-utils API → Evidence strength [0, 1]
    └─ clinvar_prior() → /api/evidence/deep_analysis → Prior [-0.2, +0.2]
  ↓
[4] Insights Bundle (bundle_insights)
    ├─ predict_protein_functionality_change → Functionality [0, 1]
    ├─ predict_gene_essentiality → Essentiality [0, 1]
    ├─ predict_chromatin_accessibility → Chromatin [0, 1]
    └─ predict_splicing_regulatory → Regulatory [0, 1]
  ↓
[5] Drug Scoring (per drug) - DrugScorer.score_drug()
    ├─ Sequence Component: seq_pct (calibrated percentile) [0, 1]
    ├─ Pathway Component: path_pct (normalized pathway score) [0, 1]
    ├─ Evidence Component: s_evd (evidence strength) [0, 1]
    ├─ ClinVar Prior: clinvar_prior [-0.2, +0.2]
    └─ Efficacy Formula: 0.3*seq_pct + 0.4*path_pct + 0.3*s_evd + clinvar_prior
  ↓
[6] Confidence Computation
    ├─ Tier Determination: compute_evidence_tier() → "supported"/"consider"/"insufficient"
    ├─ Confidence Calculation: compute_confidence() → [0, 1]
    │   └─ V2 Formula (if CONFIDENCE_V2=1): 0.5*S + 0.2*P + 0.3*E + lifts
    └─ Badge Computation: compute_evidence_badges() → ["RCT", "ClinVar-Strong", ...]
  ↓
[7] Sporadic Cancer Gates (if tumor_context provided)
    ├─ PARP Penalty/Rescue: Based on germline_status + HRD score
    ├─ IO Boost: Based on TMB/MSI status (mutually exclusive)
    └─ Confidence Capping: Based on completeness level (L0/L1/L2)
  ↓
[8] Treatment Line Modulation (if treatment_history provided)
    ├─ Line Appropriateness: Is drug appropriate for current line?
    ├─ Cross-Resistance Risk: Overlap with failed therapies
    └─ Sequencing Fitness: Optimal ordering
  ↓
[9] SAE Features Extraction (if include_sae_features=true)
    ├─ DNA Repair Capacity: 0.6*pathway_DDR + 0.2*essentiality_HRR + 0.2*exon_disruption
    └─ Mechanism Vector: 7D vector from 6 core SAE features
  ↓
Response: EfficacyResponse
  ├─ drugs[]: Ranked by confidence
  ├─ run_signature: UUID for provenance
  ├─ evidence_tier: Overall tier
  ├─ sae_features: SAE features (if requested)
  └─ provenance: Complete audit trail
```

### **3.9.3 Orchestrator Workflow** (`orchestrator.py:48-412`)

**Step-by-Step Flow**:

1. **Initialize** (lines 58-79):
   - Generate `run_id` (UUID)
   - Initialize `EfficacyResponse` with provenance
   - Set feature flags from config

2. **Get Drug Panel** (lines 82-92):
   - Load default panel (`get_default_panel()`)
   - Optional: Limit panel size for fast mode

3. **Sequence Scoring** (lines 94-102):
   - Call `SequenceProcessor.score_sequences()`
   - Returns `List[SeqScore]` with disruption and percentile

4. **Pathway Aggregation** (lines 104-105):
   - Convert `SeqScore` to dict format
   - Call `aggregate_pathways()` → pathway scores dict

5. **Evidence Gathering** (lines 112-156):
   - **Parallel execution**: `asyncio.gather()` for literature calls
   - **Timeout handling**: 30s for literature, 10s for ClinVar
   - **Fast mode**: Skips evidence if `fast=True`

6. **Insights Bundle** (lines 157-162):
   - Calls `bundle_insights()` for functionality/chromatin/essentiality/regulatory
   - **Fast mode**: Skips insights

7. **Drug Scoring Loop** (lines 181-308):
   - For each drug in panel:
     - Apply ablation mode (SP/SPE masking)
     - Call `DrugScorer.score_drug()` → S/P/E formula
     - Apply cohort lifts (if enabled)
     - **Apply sporadic gates** (if tumor_context provided)
     - **Apply treatment line modulation** (if treatment_history provided)

8. **Sort & Finalize** (lines 310-411):
   - Sort drugs by confidence (descending)
   - Extract SAE features (if requested)
   - Build provenance with confidence breakdown

### **3.9.4 Sporadic Cancer Integration**

#### **Platform Integration** (`03_PLATFORM_INTEGRATION.mdc`):

**Sequence (S) for Sporadic**:
- **Input**: Tumor NGS panel results (VCF or report)
- **Process**: Evo2 scoring for **SOMATIC mutations** (not germline)
- **Context**: `context=tumor` flag in insights endpoints
- **Key Difference**: Somatic mutations are acquired, not inherited

**Pathway (P) for Sporadic**:
- **Multi-Hit Aggregation**: Often multiple somatic hits across pathways
- **Aggregate Across**:
  - DNA repair (HRD, MMR, BER)
  - Cell cycle (TP53, CDKN2A, RB1)
  - Growth signaling (KRAS, PIK3CA, EGFR)
  - Immune evasion (PD-L1, TMB-high)

**Evidence (E) for Sporadic**:
- **Focus**: Tumor-based literature
  - Somatic mutation clinical trials
  - Drug approvals for somatic alterations
  - NCCN guidelines for tumor genomics
  - Published response rates in similar profiles

**SAE Features for Sporadic**:
- **Somatic-Specific Features**:
  - `tumor_mutational_burden` (TMB) - immunotherapy candidate
  - `microsatellite_instability` (MSI) - checkpoint response
  - `homologous_recombination_deficiency` (HRD) - somatic PARP
  - `oncogene_amplifications` (HER2, MYC, FGFR2)
  - `tumor_suppressor_loss` (TP53, PTEN, RB1)

#### **WIWFM for Sporadic** (`09_WIWFM_SPORADIC.mdc`):

**Inputs**:
- `germline_status`: "negative" | "positive" | "unknown"
- `tumor_context` (when germline negative):
  - `somatic_mutations[]`
  - `tmb`, `msi_status`, `hrd_score`
  - `copy_number_alterations[]`
- `treatment_history`
- `options` (profile flags)

**Scoring Flow**:
1. **Sequence (S)**: Score somatic mutations with Evo2, merge with CNAs
2. **Pathway (P)**: Aggregate multi-hit disruptions, weight by disease + line
3. **Evidence (E)**: Prefer tumor-based literature, apply badges
4. **SAE Features**: TMB, MSI, HRD, amplifications, losses
5. **Treatment Lines**: `line_appropriateness`, `cross_resistance_risk`, `sequencing_fitness`
6. **Germline Gating**: 
   - If germline negative + no HRD → penalize PARP
   - If TMB ≥10-15 or MSI-high → boost checkpoint

**Output**: `drugs[*]` with efficacy_score, confidence, evidence_tier, sae_features, treatment_line, rationale, provenance

### **3.9.5 Code-Level Implementation Details**

#### **Orchestrator Integration Points**:

**Sporadic Gates Integration** (`orchestrator.py:214-259`):
```python
# Extract germline status and tumor context from request
germline_status = getattr(request, 'germline_status', 'unknown')
tumor_context_data = getattr(request, 'tumor_context', None)

# Apply sporadic gates
adjusted_efficacy, adjusted_confidence, sporadic_rationale = apply_sporadic_gates(
    drug_name=drug["name"],
    drug_class=drug.get("class", ""),
    moa=drug.get("moa", ""),
    efficacy_score=drug_result.efficacy_score,
    confidence=drug_result.confidence,
    germline_status=germline_status,
    tumor_context=tumor_context_dict
)
```

**Treatment Line Integration** (`orchestrator.py:261-299`):
```python
if TREATMENT_LINE_AVAILABLE and request.treatment_history:
    treatment_hist = TreatmentHistory(**request.treatment_history)
    treatment_line_features = compute_treatment_line_features(...)
    modulated_confidence, rationale = modulate_confidence_with_treatment_line(...)
```

#### **Evidence Client Implementation**:

**Literature Client** (`literature_client.py:51-135`):
- **API**: Calls `/api/evidence/literature` endpoint
- **Timeout**: 60 seconds
- **MoA Filtering**: Prefers papers mentioning drug name or MoA
- **Strength Calculation**: RCT (0.5) > Guideline (0.35) > Review (0.25) > Other (0.15)
- **MoA Boost**: +0.10 per MoA hit, max +0.30

**ClinVar Client** (`clinvar_client.py:11-99`):
- **API**: Calls `/api/evidence/deep_analysis` endpoint
- **Timeout**: 40 seconds
- **Prior Calculation**: Pathogenic + expert/practice → +0.2, Benign + expert/practice → -0.2
- **Research-Mode Fallback**: Canonical hotspots get `prior = 0.2` if ClinVar empty

#### **Insights Bundle Implementation** (`bundle_client.py:10-149`):

**Parallel Execution**:
- **Conditional**: Only calls endpoints when required data is available
- **Sequential Execution**: Actually executes sequentially (not parallel) to avoid None task errors
- **Timeout**: 40 seconds per endpoint
- **Error Resilience**: Individual endpoint failures don't break the bundle

**Payload Requirements**:
- **Functionality**: Requires `gene` + `hgvs_p`
- **Chromatin**: Requires `chrom` + `pos` + `radius`
- **Essentiality**: Requires `gene` + variant info with `consequence`
- **Regulatory**: Requires `chrom` + `pos` + `ref` + `alt`

### **3.9.6 Formula Reference Table**

| Formula | Location | Formula | Notes |
|---------|----------|---------|-------|
| **Efficacy Score** | `drug_scorer.py:171` | `0.3*S + 0.4*P + 0.3*E + clinvar_prior` | Hardcoded, not from config |
| **Confidence V2** | `confidence_computation.py:156` | `0.5*S + 0.2*P + 0.3*E + lifts` | Only if `CONFIDENCE_V2=1` |
| **Pathway Normalization** | `drug_scorer.py:48-51` | `(s_path - 1e-6) / (1e-4 - 1e-6)` | Empirical Evo2 ranges |
| **Pathway Aggregation** | `aggregation.py:34` | `sum(seq_disruption * weight) / count` | Weighted average |
| **Literature Strength** | `literature_client.py:18-48` | `RCT(0.5) + Guideline(0.35) + Review(0.25) + Other(0.15)` | Top 3 results |
| **ClinVar Prior** | `clinvar_client.py:56-59` | `Pathogenic+Strong: +0.2, Pathogenic+Moderate: +0.1, Pathogenic+Weak: +0.05` | Review status tiers |
| **PARP Penalty** | `sporadic_gates.py:68-127` | `1.0x (germline+), 1.0x (HRD≥42), 0.6x (HRD<42), 0.8x (unknown)` | Germline gating |
| **IO Boost** | `sporadic_gates.py:138-196` | `1.35x (TMB≥20), 1.30x (MSI-H), 1.25x (TMB≥10)` | Mutually exclusive |
| **Confidence Cap** | `sporadic_gates.py:206-230` | `L0: 0.4, L1: 0.6, L2: no cap` | By completeness |

### **3.9.7 Decision Points & Gates**

| Decision Point | Location | Logic | Impact |
|----------------|----------|-------|--------|
| **Fusion Eligibility** | `sequence_processor.py:32-40` | `GRCh38 AND missense` | Uses Fusion if eligible, else Evo2 |
| **Evidence Gate** | `tier_computation.py:51-54` | `s_evd ≥ 0.7 OR (ClinVar-Strong AND path ≥ 0.2)` | Tier = "supported" |
| **Insufficient Signal** | `tier_computation.py:57-61` | `seq < 0.02 AND path < 0.05 AND evd < 0.2` | Tier = "insufficient" |
| **PARP Rescue** | `sporadic_gates.py:85-94` | `germline negative AND HRD ≥ 42` | Efficacy = 1.0x (rescued!) |
| **IO Boost Priority** | `sporadic_gates.py:156-189` | `TMB≥20 > MSI-H > TMB≥10` | Highest priority wins |
| **Confidence Cap** | `sporadic_gates.py:206-230` | `L0: 0.4, L1: 0.6, L2: no cap` | Limits confidence by data quality |

### **3.9.8 Key Architectural Patterns**

1. **Orchestrator Pattern**: Single powerful orchestrator endpoint (`/api/efficacy/predict`) coordinates all services
2. **Fallback Chain**: Fusion → Evo2 → Massive Oracle (hierarchical scoring)
3. **Parallel Evidence**: `asyncio.gather()` for concurrent literature calls
4. **Graceful Degradation**: Evidence timeout → `tier = "insufficient"`, continues with S/P only
5. **Provenance Tracking**: Complete audit trail in `response.provenance`
6. **Feature Flags**: `CONFIDENCE_V2`, `DISABLE_FUSION`, `EVO_USE_DELTA_ONLY`, etc.
7. **Ablation Modes**: SP/SPE masking for component testing
8. **Fast Mode**: Skips evidence and insights for speed

### **3.9.9 Sporadic Cancer Workflow Integration**

**Complete Flow for Sporadic Cases**:
1. **Tumor NGS Ingestion**: Parse Foundation/Tempus reports → `TumorContext` (L0/L1/L2)
2. **Germline Status Check**: `germline_status = "negative"` (85-90% of cases)
3. **S/P/E Scoring**: Score somatic mutations (not germline)
4. **Sporadic Gates Applied**:
   - PARP: Check HRD score for rescue (HRD ≥42 → 1.0x)
   - IO: Check TMB/MSI for boost (TMB≥20 → 1.35x, MSI-H → 1.30x)
   - Confidence: Cap by completeness level (L0: 0.4, L1: 0.6)
5. **Treatment Line Intelligence**: Compute `line_appropriateness`, `cross_resistance_risk`
6. **SAE Features**: Extract DNA repair capacity, mechanism vector
7. **Response**: Ranked drugs with sporadic-adjusted scores

**Clinical Impact**:
- **85-90% of cases** are sporadic (germline negative)
- **PARP rescue** enables PARP inhibitors for HRD-high sporadic patients
- **IO boost** identifies immunotherapy candidates (TMB-high/MSI-H)
- **Confidence capping** reflects data quality (L0/L1/L2)

---

---

## **3.10 CODE IMPLEMENTATION DEEP DIVE** (Cycle 2 - Technical Details)

### **3.10.1 Evo2 Scoring Algorithm - Complete Implementation**

#### **Multi-Window Strategy** (`evo2_scorer.py:200-274`):

**Algorithm Flow**:
1. **Multi-Window Call** (`score_variant_multi`):
   - Calls `/api/evo/score_variant_multi` with GRCh38 coordinates
   - Returns `min_delta` (best delta across all windows)
   - **No exon context** - just sequence-level delta

2. **Exon-Context Scanning** (if not delta-only mode):
   - Tests multiple flank sizes: `[4096, 8192, 16384, 25000]` bp
   - For each flank, calls `/api/evo/score_variant_exon`
   - Picks **best `exon_delta`** (strongest absolute value)
   - Records `best_window_bp` (best flank size)

3. **Delta-Only Mode** (`evo_use_delta_only` flag):
   - **Skips exon scanning entirely** for speed/cost
   - Uses only `min_delta` from multi-window call
   - **When enabled**: `force_exon_scan=False` → skip exon loop

**Code** (`evo2_scorer.py:223-234`):
```python
# Spam-safe: if delta-only mode is enabled, skip exon loop entirely
if feature_flags.get("evo_use_delta_only", False) and not force_exon_scan:
    return {
        "min_delta": (j_multi or {}).get("min_delta"),
        "exon_delta": None,
        "best_window_bp": None,
        "windows_tested": [],
    }
```

#### **Forward/Reverse Symmetry** (`evo2_scorer.py:276-332`):

**Algorithm**:
1. **Forward Direction**: `ref → alt` (standard variant)
2. **Reverse Direction**: `alt → ref` (symmetry check)
3. **Averaging**: `avg_min_delta = (forward_min + reverse_min) / 2.0`
4. **Feature Flag**: `evo_disable_symmetry` (default: `True` - symmetry **disabled**)

**Why Symmetry**:
- Evo2 is strand-agnostic (doesn't care about forward/reverse)
- Averaging reduces noise and improves robustness
- **Default disabled** for performance (doubles API calls)

**Code** (`evo2_scorer.py:307-314`):
```python
# Average the delta scores for symmetry
forward_min = forward_result.get("min_delta") or 0.0
reverse_min = reverse_result.get("min_delta") or 0.0
avg_min_delta = (forward_min + reverse_min) / 2.0

forward_exon = forward_result.get("exon_delta") or 0.0
reverse_exon = reverse_result.get("exon_delta") or 0.0
avg_exon_delta = (forward_exon + reverse_exon) / 2.0
```

#### **Model Selection** (`evo2_scorer.py:95-109`):

**Priority Order**:
1. **`EVO_FORCE_MODEL`** env var → Forces single model (highest priority)
2. **`EVO_ALLOWED_MODELS`** env var → Filters allowed models
3. **Default**: `["evo2_1b", "evo2_7b", "evo2_40b"]` if ensemble, else `[model_id]`

**Best Model Selection**:
- Tests all candidate models
- Picks model with **strongest `abs(min_delta)`**
- **Why**: Stronger delta = more confident prediction

**Code** (`evo2_scorer.py:118-122`):
```python
# Keep the best result across models
if result.get("min_delta") is not None:
    if (best["min_delta"] is None or 
        abs(float(result["min_delta"])) > abs(float(best["min_delta"] or 0))):
        best.update(result)
        best["model"] = model
```

#### **Sequence Disruption Calculation** (`evo2_scorer.py:126-136`):

**Formula**: `sequence_disruption = max(abs(min_delta), abs(exon_delta))`

**Why**: Uses **stronger signal** between:
- Multi-window `min_delta` (sequence-level)
- Exon-context `exon_delta` (exon-level with adaptive flanks)

**Code**:
```python
exon_abs = abs(float(best.get("exon_delta") or 0.0))
min_abs = abs(float(best.get("min_delta") or 0.0))
sequence_disruption = max(min_abs, exon_abs)
```

#### **Hotspot Floors** (`evo2_scorer.py:138-149`):

**Raw Delta Level** (`HOTSPOT_FLOOR = 1e-4`):
- **BRAF V600**: `max(disruption, 1e-4)`
- **KRAS/NRAS/HRAS G12/G13/Q61**: `max(disruption, 1e-4)`
- **TP53 R175/R248/R273**: `max(disruption, 1e-4)`

**Why**: Known pathogenic hotspots should never score below threshold, even if Evo2 gives low delta

**Percentile Level** (`evo2_scorer.py:159-171`):
- **BRAF V600**: `max(pct, 0.90)` (90th percentile floor)
- **KRAS/NRAS/HRAS G12/G13/Q61**: `max(pct, 0.80)` (80th percentile floor)
- **TP53 R175/R248/R273**: `max(pct, 0.80)` (80th percentile floor)

**Why**: Hotspots should always be in high percentile range, regardless of raw delta

#### **Truncation/Frameshift Lift** (`evo2_scorer.py:150-157`):

**Logic**: If `"*" in hgvs_p` (stop codon) OR `"FS" in hgvs_p` (frameshift) → `disruption = 1.0`

**Why**: Stop codons and frameshifts are **highly disruptive** → enforce maximum disruption

**Code**:
```python
if "*" in hgvs_p or "FS" in hgvs_p:
    sequence_disruption = max(sequence_disruption, 1.0)
```

### **3.10.2 Fusion Scorer Implementation** (`fusion_scorer.py`)

#### **Variant Format Handling** (`fusion_scorer.py:74-85`):

**Multiple Format Attempts**:
```python
candidates = [
    f"chr{chrom}:{pos}:{ref}:{alt}",      # chr7:140453136:T:A
    f"{chrom}:{pos}:{ref}:{alt}",         # 7:140453136:T:A
    f"chr{chrom}:{pos}:{alt}:{ref}",      # chr7:140453136:A:T (reverse)
    f"{chrom}:{pos}:{alt}:{ref}",        # 7:140453136:A:T (reverse)
]
```

**Why**: Different APIs expect different formats - try all variants

#### **Fallback Chain** (`fusion_scorer.py:70-100`):

1. **External Fusion Engine** (`FUSION_AM_URL`):
   - Calls external Fusion Engine API
   - Prefers fused scores, falls back to AlphaMissense

2. **Local Router** (`/api/fusion/score_variant`):
   - If external unavailable, tries local router
   - Same API contract, different endpoint

3. **Graceful Degradation**:
   - If both fail, returns empty list
   - Sequence processor falls back to Evo2

#### **Caching** (`fusion_scorer.py:26-39, 62-65`):

**Cache Key**: `fusion_am:{hash(tuple(mutation_keys))}`

**Why**: Fusion scoring is expensive - cache results for 1 hour (TTL)

### **3.10.3 Pathway Panel Configuration** (`panel_config.py`)

#### **Default MM Panel** (`panel_config.py:8-14`):

**Hardcoded Configuration**:
```python
DEFAULT_MM_PANEL = [
    {"name": "BRAF inhibitor", "moa": "MAPK blockade", "pathway_weights": {"ras_mapk": 0.8, "tp53": 0.2}},
    {"name": "MEK inhibitor", "moa": "MAPK downstream blockade", "pathway_weights": {"ras_mapk": 0.9, "tp53": 0.1}},
    {"name": "IMiD", "moa": "immunomodulatory", "pathway_weights": {"ras_mapk": 0.2, "tp53": 0.3}},
    {"name": "Proteasome inhibitor", "moa": "proteostasis stress", "pathway_weights": {"ras_mapk": 0.3, "tp53": 0.4}},
    {"name": "Anti-CD38", "moa": "antibody", "pathway_weights": {"ras_mapk": 0.1, "tp53": 0.1}},
]
```

**Pathway Weights Rationale**:
- **BRAF/MEK inhibitors**: High `ras_mapk` weight (0.8-0.9) - directly target MAPK pathway
- **IMiD/Proteasome**: Lower `ras_mapk`, higher `tp53` - indirect effects
- **Anti-CD38**: Low weights (0.1) - antibody, pathway-agnostic

**Limitation**: Currently hardcoded for Multiple Myeloma only - needs extension for other diseases

### **3.10.4 Sporadic Integration - NGS Ingestion**

#### **Tumor NGS Ingestion** (`04_TUMOR_NGS_INGESTION.mdc`):

**Endpoint**: `POST /api/tumor/ingest_ngs`

**Input**: PDF/JSON from Foundation Medicine or Tempus

**Output**: `TumorContext` JSON with:
- `somatic_mutations[]`
- `tmb` (mutations/megabase)
- `msi_status` (MSI-high/MSS)
- `hrd_score` (somatic HRD)
- `copy_number_alterations[]`
- `purity`, `ploidy`

**Implementation Status**: **Planned** (Day 3) - Not yet implemented

**Acceptance Criteria**: Both Foundation and Tempus reports produce valid `TumorContext`

### **3.10.5 Sporadic Integration - Trial Matching**

#### **Clinical Trial Filtering** (`05_CLINICAL_TRIAL_MATCHING.mdc`):

**Exclude** (germline-required):
- ❌ "germline BRCA required" trials
- ❌ "Lynch syndrome required" trials

**Include** (somatic-focused):
- ✅ Tumor-agnostic (TMB-high, MSI-high)
- ✅ Somatic biomarker trials (somatic HRD, TP53)
- ✅ Mechanism-based (VEGF, PARP combo, IO)
- ✅ Post-progression trials (L2+, L3+)

**Ranking Priority**:
1. Biomarker match (TMB, MSI, HRD)
2. Mechanism alignment (pathway match)
3. Line-eligibility (L1/L2/L3)

### **3.10.6 Sporadic Integration - Backend Contracts**

#### **TumorContext Schema** (`08_BACKEND_CONTRACTS.mdc`):

```python
class TumorContext(BaseModel):
    somatic_mutations: List[Mutation]
    tmb: Optional[float] = None  # mutations/megabase
    msi_status: Optional[str] = None  # MSI-high/MSS
    hrd_score: Optional[float] = None  # somatic HRD
    copy_number_alterations: Optional[List[CNA]] = None
    purity: Optional[float] = None
    ploidy: Optional[float] = None
    completeness_score: float  # L0/L1/L2 completeness
    level: str  # "L0" | "L1" | "L2"
```

#### **Updated Efficacy Endpoint** (`/api/efficacy/predict`):

**New Inputs**:
- `germline_status`: "negative" | "positive" | "unknown"
- `tumor_context?`: `TumorContext` (optional)

**Behavior**:
- PARP penalty when `germline_status == "negative"` (unless HRD ≥42)
- TMB boost when `tmb >= 20` (checkpoint inhibitors)
- MSI boost when `msi_status == "MSI-H"` (checkpoint inhibitors)
- Confidence capping by `completeness_score` (L0: 0.4, L1: 0.6, L2: no cap)

### **3.10.7 Edge Cases & Error Handling**

#### **Evo2 Scoring Edge Cases**:

1. **Missing Coordinates**:
   - **Handling**: `safe_str()`, `safe_int()` utilities return defaults
   - **Result**: Skips invalid mutations, continues with valid ones

2. **API Timeout**:
   - **Timeout**: 180 seconds for Evo2 calls
   - **Handling**: Exception caught, continues to next mutation/model
   - **Result**: Partial results (some mutations scored, others skipped)

3. **Cache Failures**:
   - **Handling**: `get_cache()` returns `None` on failure
   - **Result**: Falls back to API call (no blocking)

4. **Invalid Delta Values**:
   - **Handling**: `safe_float()` returns `0.0` on invalid input
   - **Result**: Conservative default (no disruption)

#### **Evidence Gathering Edge Cases**:

1. **Literature Timeout**:
   - **Timeout**: 60 seconds
   - **Handling**: `asyncio.TimeoutError` caught
   - **Result**: `evidence_strength = 0.0`, `tier = "insufficient"`

2. **ClinVar Timeout**:
   - **Timeout**: 40 seconds
   - **Handling**: `asyncio.TimeoutError` caught
   - **Result**: `clinvar_prior = 0.0`, continues with S/P only

3. **Empty Results**:
   - **Handling**: `_score_evidence_from_results([])` returns `0.0`
   - **Result**: Conservative default (no evidence)

#### **Sporadic Gates Edge Cases**:

1. **Missing Tumor Context**:
   - **Handling**: `tumor_context is None` → uses defaults
   - **Result**: PARP penalty (0.8x), no IO boost, no confidence cap

2. **Missing HRD Score**:
   - **Handling**: `hrd_score is None` → assumes unknown
   - **Result**: PARP penalty (0.8x) for germline-negative

3. **Invalid TMB/MSI**:
   - **Handling**: `tmb is None` or `msi_status is None` → no boost
   - **Result**: No IO boost applied

### **3.10.8 Performance Optimizations**

#### **Caching Strategy**:

1. **Evo2 Caching**:
   - **Key**: `evo2:{model_id}:{hgvs_p}:{hash(window_flanks)}:{ensemble}`
   - **TTL**: 1 hour (3600 seconds)
   - **Location**: Redis (if configured) or in-memory

2. **Fusion Caching**:
   - **Key**: `fusion_am:{hash(tuple(mutation_keys))}`
   - **TTL**: 1 hour
   - **Location**: Redis (if configured)

3. **Evidence Caching**:
   - **Not currently implemented** (future optimization)

#### **Parallel Execution**:

1. **Evidence Gathering**:
   - **Method**: `asyncio.gather(*evidence_tasks)`
   - **Benefit**: Concurrent literature calls (one per drug)
   - **Timeout**: 30 seconds total

2. **Insights Bundle**:
   - **Method**: Sequential (not parallel) - see `bundle_client.py:111-121`
   - **Why**: Conditional execution (only calls when data available)
   - **Future**: Could be parallelized with `asyncio.gather()`

#### **Spam-Safety Controls**:

1. **Delta-Only Mode**:
   - **Flag**: `EVO_USE_DELTA_ONLY=1`
   - **Effect**: Skips expensive exon scanning
   - **Benefit**: Faster, lower cost

2. **Model Limiting**:
   - **Flag**: `EVO_FORCE_MODEL=evo2_1b`
   - **Effect**: Forces single model (no ensemble)
   - **Benefit**: Faster, lower cost

3. **Window Limiting**:
   - **Flag**: `EVO_MAX_FLANKS=2`
   - **Effect**: Limits number of windows tested
   - **Benefit**: Faster, lower cost

---

**Status**: ✅ **CYCLE 2 COMPLETE** - S/P/E Code Implementation Deep Dive  
**Next**: Cycle 3 (SC-I1) - Sporadic Cancer Strategic Foundation

---