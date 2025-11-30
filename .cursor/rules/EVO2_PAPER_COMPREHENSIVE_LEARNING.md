# 🧬 EVO2 PAPER - COMPREHENSIVE LEARNING & CODEBASE CONNECTION

**Date:** 2025-01-14  
**Purpose:** Master-level understanding of Evo2 paper principles and their manifestation in our codebase  
**Approach:** Section-by-section analysis with explicit codebase connections

---

## 📋 EXECUTIVE SUMMARY

**Evo2 Core Principle:** A biological foundation model trained on 9.3T DNA base pairs that learns evolutionary constraints to predict variant impact and generate sequences **without task-specific finetuning** (zero-shot).

**Key Innovation:** Single-nucleotide resolution with 1M-token context enables:
1. **Zero-shot variant effect prediction** across coding/noncoding/splice regions
2. **Mechanistic interpretability** via sparse autoencoders (SAEs) revealing biological features
3. **Genome-scale generation** with controllable epigenomic design via inference-time search

**Our Implementation:** We use Evo2 for:
- Sequence disruption scoring in S/P/E framework (Multiple Myeloma, Ayesha)
- Variant impact prediction via delta likelihoods
- Pathway aggregation based on sequence scores
- Hotspot-aware floors and percentile calibration

---

## 📖 SECTION-BY-SECTION DEEP DIVE

### **1. ABSTRACT & INTRODUCTION**

#### **Paper Key Points:**
- **Scale:** 7B/40B parameters, 9.3T tokens, 1M-token context window
- **Training Data:** OpenGenome2 spanning all domains of life (bacteria, archaea, eukarya, phage)
- **Zero-shot capability:** Predicts functional impacts without variant-specific training
- **Interpretability:** SAEs reveal exon-intron boundaries, TF motifs, protein structure, prophage regions
- **Generation:** Genome-scale sequences (mitochondrial, prokaryotic, eukaryotic) with higher naturalness
- **Controllability:** Inference-time search enables epigenomic design (first inference-time scaling in biology)

#### **Codebase Connection:**
```python
# Our Evo2 Service (src/services/evo_service/main.py)
# Uses Evo2 1B model (smaller, faster for production)
# Implements zero-shot scoring via delta likelihoods

@self.fastapi_app.post("/score_variant")
def score_variant(item: dict):
    # Extract 8,192 bp window around variant
    ref_sequence = seq
    alt_sequence = seq[:idx] + alt + seq[idx+1:]
    
    # Zero-shot scoring: delta = alt_ll - ref_ll
    ll = self.model.score_sequences([ref_sequence, alt_sequence])
    delta = alt_ll - ref_ll  # Negative = deleterious
```

**Why This Matters:**
- **No training needed:** We can score any variant immediately
- **Noncoding coverage:** Unlike AlphaMissense (protein-only), Evo2 handles noncoding variants
- **Splice awareness:** Evo2 learned exon-intron boundaries, enabling splice variant prediction

---

### **2.1 MODEL ARCHITECTURE & TRAINING**

#### **Paper Key Points:**
- **Architecture:** StripedHyena 2 (multi-hybrid: input-dependent convolutions + attention)
- **Two-phase training:**
  1. **Pretraining:** 8,192 token context, genic-weighted data (functional elements)
  2. **Midtraining:** Extend to 1M tokens, whole genomes (long-range dependencies)
- **Efficiency:** 1.3× throughput at 16k, ~3× at 1M vs Transformers
- **Data:** OpenGenome2 (8.8T nt), excludes eukaryotic viruses (safety)

#### **Codebase Connection:**
```python
# Our Adaptive Window Strategy (evo2_scorer.py)
# Tests multiple window sizes to find optimal context

window_flanks = [4096, 8192, 16384, 25000]  # Adaptive windows

# Rationale: Different variants need different context
# - Exonic variants: 4K-8K (gene-level)
# - Regulatory variants: 16K-25K (long-range enhancers)
# - Splice variants: 8K-16K (exon-intron boundaries)
```

**Why This Matters:**
- **Context matters:** 8K for functional elements, 1M for long-range regulation
- **Our adaptation:** We test multiple windows to find optimal context per variant
- **Efficiency trade-off:** We use 1B model (faster) vs 7B/40B (more accurate) for production

---

### **2.2 MUTATIONAL EFFECTS PREDICTION**

#### **Paper Key Points:**
- **Zero-shot across modalities:** DNA, RNA, protein fitness prediction
- **Canonical signals:** Start/stop codon sensitivity, triplet periodicity, RBS/Kozak patterns
- **Experimental concordance:** Correlates with DMS fitness, mRNA decay rates
- **Eukaryotic features:** Exon-intron classifier outperforms baselines
- **Organismal fitness:** Gene essentiality prediction (bacteria, phage, human lncRNA)

#### **Codebase Connection:**
```python
# Our Sequence Disruption Calculation (evo2_scorer.py:136)
sequence_disruption = max(abs(min_delta), abs(exon_delta))

# Paper Principle: Use maximum of multiple windows
# - min_delta: Minimum window delta (local impact)
# - exon_delta: Exon-context delta (functional element impact)
# Rationale: Capture both local and exon-level disruption
```

**Why This Matters:**
- **Multi-window strategy:** Paper shows 4-8kb context needed for some tasks (e.g., ciliate stop codons)
- **Our implementation:** We use adaptive windows + max aggregation to capture both local and functional impact
- **Hotspot floors:** We apply known oncogenic mutation floors (BRAF V600*, KRAS G12/G13) based on paper's finding that Evo2 captures evolutionary constraints

---

### **2.3 HUMAN CLINICAL VARIANT PREDICTION**

#### **Paper Key Points:**
- **ClinVar:** Competitive on coding SNVs (4th/5th), SOTA for noncoding and non-SNVs
- **SpliceVarDB:** Best zero-shot on exonic and intronic splice variants
- **BRCA1/2:** SOTA on noncoding SNVs, best combined coding+noncoding
- **No human variant training:** Uses multi-species constraints (only human reference in training)
- **Supervised boost:** Evo2 embeddings → BRCA1 classifier (AUROC 0.95)

#### **Codebase Connection:**
```python
# Our Variant Scoring Pipeline (sequence_processor.py)
# 1. Try Fusion first (GRCh38 missense only)
# 2. Fall back to Evo2 (all variant types)
# 3. Apply hotspot-aware floors

# Fusion = AlphaMissense (protein-only, GRCh38 missense)
# Evo2 = All variants (coding/noncoding/splice)

# Rationale: Use best tool for each variant type
```

**Why This Matters:**
- **Variant type coverage:** Evo2 handles indels, noncoding, splice (where Fusion/AlphaMissense fail)
- **Our strategy:** Fusion for missense (more accurate), Evo2 for everything else (broader coverage)
- **Clinical relevance:** Paper shows Evo2 is SOTA for noncoding variants (98% of genome)

---

### **2.4 MECHANISTIC INTERPRETABILITY (SAEs)**

#### **Paper Key Points:**
- **SAE setup:** Batch-TopK SAE on layer-26, 1B tokens (balanced euk/prok)
- **Features discovered:**
  - Mobile elements (prophage, CRISPR spacers)
  - Genome organization (ORFs, tRNAs, rRNAs, intergenic)
  - Protein structure (alpha-helices, beta-sheets)
  - Regulatory (TF motifs, exon-intron boundaries)
  - Mutational severity (frameshift/stop-sensitive features)
- **Generalization:** Features work on woolly mammoth (not in training)

#### **Codebase Connection:**
```python
# Our SAE Integration (sae_feature_service.py)
# We use SAE features for:
# - DNA repair capacity (pathway_ddr)
# - Mechanism vector (7D: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
# - Resistance signals
# - Hotspot detection

# Paper Principle: SAE features capture biological concepts
# Our Application: Use SAE features for pathway-level predictions
```

**Why This Matters:**
- **Feature interpretability:** SAE features align with known biology (exons, motifs, structure)
- **Our use case:** We extract SAE features for mechanism fit ranking (trial matching)
- **Future potential:** Could use SAE features for regulatory variant annotation

---

### **2.5 GENOME-SCALE GENERATION**

#### **Paper Key Points:**
- **Gene completion:** High sequence recovery across diverse species
- **Mitochondrial genomes:** 16 kb with correct gene counts, synteny, plausible structures
- **Prokaryotic genomes:** 600 kb M. genitalium with ~70% Pfam hits (vs 18% Evo1)
- **Eukaryotic chromosomes:** 330 kb yeast with tRNAs, promoters, introns
- **Safety:** Poor generation for human viral proteins (by design)

#### **Codebase Connection:**
```python
# Our Generation Use Case (metastasis_interception_service.py)
# We use Evo2 for guide RNA design:
# - Target Lock Score (multi-modal aggregation)
# - Assassin Score (efficacy + safety + structural confidence)
# - Structural validation (AlphaFold 3 for RNA:DNA complexes)

# Paper Principle: Evo2 can generate realistic sequences
# Our Application: Design guide RNAs for CRISPR editing
```

**Why This Matters:**
- **Design capability:** Evo2 can generate sequences that "look natural" to cellular machinery
- **Our application:** Guide RNA design for metastasis interception (CRISPR-based therapy)
- **Validation:** We use AlphaFold 3 (like paper) to validate generated structures

---

### **2.6 GENERATIVE EPIGENOMICS (INFERENCE-TIME SEARCH)**

#### **Paper Key Points:**
- **Goal:** Control chromatin accessibility patterns (position + width)
- **Method:** Beam search + Enformer/Borzoi scoring
  - Sample K chunks (128 bp each)
  - Score with Enformer/Borzoi ensemble
  - Keep top K' candidates, iterate
- **Inference-time scaling:** More compute → better designs (first scaling law in biology)
- **Demonstrations:** Morse code patterns ("LO", "ARC", "EVO2") as accessibility peaks

#### **Codebase Connection:**
```python
# Our Inference-Time Guidance (metastasis_interception_service.py)
# We use similar principles:
# - Multi-modal scoring (Evo2 + Enformer + AlphaFold 3)
# - Beam search for guide RNA candidates
# - Structural validation as "scoring function"

# Paper Principle: Guide generation with sequence-to-function models
# Our Application: Guide RNA design with structural validation
```

**Why This Matters:**
- **Controllability:** Inference-time search enables targeted design
- **Our application:** Could extend to chromatin accessibility design for gene therapy
- **Scaling law:** More compute → better designs (actionable knob for future work)

---

## 🔗 CODEBASE IMPLEMENTATION PATTERNS

### **Pattern 1: Zero-Shot Delta Likelihoods**

**Paper Principle:**
> "Variants assigned a more negative log-likelihood change from reference are considered to be more deleterious."

**Our Implementation:**
```python
# src/services/evo_service/main.py:199-202
ll = self.model.score_sequences([ref_sequence, alt_sequence])
delta = alt_ll - ref_ll  # Negative = deleterious
```

**Usage:**
- `sequence_processor.py`: Calls Evo2 for variant scoring
- `evo2_scorer.py`: Adaptive windows + ensemble support
- `drug_scorer.py`: Converts sequence scores to pathway scores

---

### **Pattern 2: Multi-Window Adaptive Strategy**

**Paper Principle:**
> "4 to 8 kb context windows around the premature stop codons were necessary to correctly identify the ciliate stop codon code."

**Our Implementation:**
```python
# evo2_scorer.py:48
window_flanks = [4096, 8192, 16384, 25000]  # Adaptive windows

# Rationale: Different variants need different context
# - Exonic: 4K-8K (gene-level)
# - Regulatory: 16K-25K (long-range)
# - Splice: 8K-16K (exon-intron)
```

**Usage:**
- Tests multiple windows, selects best score
- Forward/reverse symmetry averaging
- Exon-context scanning for functional elements

---

### **Pattern 3: Hotspot-Aware Floors**

**Paper Principle:**
> "Evo 2 learns evolutionary constraints from 9.3T tokens, capturing known biological patterns."

**Our Implementation:**
```python
# evo2_scorer.py:143-148
HOTSPOT_FLOORS = {
    "BRAF": 0.50,  # V600* mutations
    "KRAS": 0.40,  # G12/G13/Q61
    "NRAS": 0.40,
    "TP53": 0.35   # R175/R248/R273
}

# Apply floors to known oncogenic mutations
sequence_disruption = max(sequence_disruption, HOTSPOT_FLOOR)
```

**Usage:**
- Ensures known oncogenic mutations get high scores
- Prevents false negatives for clinically important variants
- Based on paper's finding that Evo2 captures evolutionary constraints

---

### **Pattern 4: Percentile Calibration**

**Paper Principle:**
> "Evo 2 likelihoods correlate with diverse definitions of fitness for both prokaryotic and eukaryotic protein and ncRNA molecules."

**Our Implementation:**
```python
# gene_calibration.py:percentile_like
# Converts raw sequence_disruption to percentile (0-1 scale)
# Uses gene-specific calibration curves

pct = percentile_like(sequence_disruption)
# Hotspot percentile clamping: max(pct, 0.90) for BRAF
```

**Usage:**
- Normalizes Evo2 scores to 0-1 scale
- Gene-specific calibration (different genes have different score distributions)
- Hotspot percentile clamping (known oncogenic mutations → high percentiles)

---

### **Pattern 5: Pathway Aggregation**

**Paper Principle:**
> "Evo 2 learns sequence features that contribute to molecular fitness across both domains."

**Our Implementation:**
```python
# pathway/aggregation.py:aggregate_pathways
# Aggregates sequence scores by pathway weights
# Gene → Pathway mapping (BRAF/KRAS → ras_mapk, TP53 → tp53)

pathway_scores = {}
for pathway, weight in pathway_weights.items():
    pathway_totals[pathway] += sequence_disruption * weight
```

**Usage:**
- Converts variant-level scores to pathway-level scores
- Drug-specific pathway weights (e.g., BRAF inhibitor → ras_mapk weight 0.8)
- Used in S/P/E framework for drug efficacy prediction

---

## 🎯 KEY INSIGHTS FOR OUR APPLICATION

### **1. Zero-Shot Capability = No Training Needed**

**Paper Finding:**
- Evo2 predicts variant impact without variant-specific training
- Uses evolutionary constraints learned from 9.3T tokens

**Our Benefit:**
- Can score any variant immediately (no model training)
- Works for rare variants (no need for training data)
- Handles noncoding variants (98% of genome)

---

### **2. Noncoding Variant Coverage = Clinical Advantage**

**Paper Finding:**
- Evo2 is SOTA for noncoding variants (ClinVar, SpliceVarDB)
- Outperforms AlphaMissense/GPN-MSA on non-SNVs

**Our Benefit:**
- Can score splice variants (intronic + exonic)
- Can score regulatory variants (enhancers, promoters)
- Can score indels (where protein models fail)

---

### **3. Multi-Window Strategy = Better Context**

**Paper Finding:**
- 4-8kb context needed for some tasks (e.g., ciliate stop codons)
- 1M context enables long-range regulatory relationships

**Our Benefit:**
- Adaptive windows find optimal context per variant
- Exon-context scanning captures functional element impact
- Forward/reverse symmetry improves robustness

---

### **4. SAE Features = Interpretable Biology**

**Paper Finding:**
- SAE features align with known biology (exons, motifs, structure)
- Features generalize to unseen genomes (woolly mammoth)

**Our Benefit:**
- Use SAE features for mechanism fit ranking (trial matching)
- DNA repair capacity from SAE pathway features
- Future: Regulatory variant annotation via SAE features

---

### **5. Inference-Time Search = Controllable Design**

**Paper Finding:**
- Beam search + scoring function enables targeted design
- More compute → better designs (scaling law)

**Our Benefit:**
- Guide RNA design with structural validation
- Future: Chromatin accessibility design for gene therapy
- Actionable knob: Increase compute for better designs

---

## 🚨 LIMITATIONS & CONSIDERATIONS

### **1. Model Size Trade-off**

**Paper:**
- 7B/40B models more accurate, but slower
- 1B model faster, but less accurate

**Our Choice:**
- Use 1B model for production (speed)
- Could upgrade to 7B for research (accuracy)

---

### **2. Context Window Limits**

**Paper:**
- 1M-token context enables long-range relationships
- But requires significant compute

**Our Implementation:**
- Use adaptive windows (4K-25K) for efficiency
- Could extend to longer windows for regulatory variants

---

### **3. Population Bias**

**Paper:**
- Evo2 has similar ancestry bias as other population-free methods
- Non-European variants scored as more pathogenic

**Our Consideration:**
- Need to validate across populations
- Consider population-specific calibration

---

### **4. Safety Exclusions**

**Paper:**
- Eukaryotic viruses excluded from training (safety)
- Evo2 performs poorly on viral sequences (by design)

**Our Consideration:**
- Cannot use Evo2 for viral genomics
- Need alternative models for viral variant prediction

---

## 📚 REFERENCES & FURTHER READING

### **Paper Sections:**
- **Abstract:** Core capabilities and scale
- **2.1:** Architecture and training details
- **2.2:** Mutational effects prediction
- **2.3:** Clinical variant prediction (ClinVar, BRCA1/2)
- **2.4:** SAE interpretability
- **2.5:** Genome-scale generation
- **2.6:** Inference-time search (epigenomics)
- **4. Methods:** Technical details

### **Codebase Files:**
- `src/services/evo_service/main.py`: Evo2 Modal service
- `oncology-coPilot/oncology-backend-minimal/api/services/sequence_scorers/evo2_scorer.py`: Adaptive scoring
- `oncology-coPilot/oncology-backend-minimal/api/services/pathway/aggregation.py`: Pathway aggregation
- `oncology-coPilot/oncology-backend-minimal/api/services/gene_calibration.py`: Percentile calibration
- `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py`: SAE features

### **Related Documents:**
- `.cursor/concept/evo2-paper-review.md`: Original section-by-section review
- `.cursor/rules/EVO2_DEEP_DIVE_ANALYSIS.md`: 10+ iteration deep dive
- `.cursor/rules/WIWFMSPE_MM_MASTER.mdc`: S/P/E framework (uses Evo2)

---

## ✅ VALIDATION CHECKLIST

- [x] Understand zero-shot prediction principle
- [x] Understand multi-window strategy rationale
- [x] Understand hotspot-aware floors
- [x] Understand percentile calibration
- [x] Understand pathway aggregation
- [x] Understand SAE feature extraction
- [x] Understand inference-time search
- [x] Connect paper principles to codebase implementation
- [x] Identify limitations and considerations
- [x] Document key insights for our application

---

**Status:** ✅ **COMPREHENSIVE UNDERSTANDING ACHIEVED**

**Next Steps:**
1. Review S/P/E framework integration (WIWFMSPE_MM_MASTER.mdc)
2. Review SAE feature extraction (sae_feature_service.py)
3. Review mechanism fit ranking (mechanism_fit_ranker.py)
4. Review metastasis interception (metastasis_interception_service.py)



