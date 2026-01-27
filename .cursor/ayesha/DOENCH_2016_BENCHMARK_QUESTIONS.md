# Questions for Manager: Doench 2016 Benchmark Implementation
**Date**: January 2025  
**Context**: Preparing Evo2 benchmark against Doench 2016 dataset for publication

---

## Critical Questions

### 1. Genomic Context Mapping (BLOCKING)

**Issue**: Doench 2016 dataset provides guide sequences (30mer) and gene names, but no genomic coordinates. Evo2 requires 8,192 bp genomic context for scoring.

**Questions**:
- Is there a pre-computed mapping of Doench guide sequences to genomic coordinates (chromosome, position)?
- If not, what's the recommended approach?
  - Option A: Use Ensembl API to map guide sequences to genomic positions (may be slow for 4,390 guides)
  - Option B: Use existing CrisPRO infrastructure if we have gene-to-coordinate mapping
  - Option C: Construct synthetic context around guide sequences (less ideal)
- How do we handle guides that map to multiple genomic locations (e.g., repetitive regions)?
- For genes with multiple transcripts, which transcript should we use for context extraction?

**Impact**: This is the primary bottleneck. Without genomic coordinates, we cannot extract proper Evo2 context.

---

### 2. Evo2 Scoring Approach

**Issue**: Need to determine the correct Evo2 scoring method for guide RNA efficacy prediction.

**Questions**:
- **Delta scoring**: Should we compute delta (mutant vs reference) where:
  - Reference = wild-type genomic sequence
  - Mutant = genomic sequence with guide RNA binding site?
- **Direct scoring**: Or should we score the guide RNA sequence directly using Evo2's sequence likelihood?
- **Context construction**: For delta scoring, how do we construct the "mutant" sequence?
  - Do we simulate guide binding (DNA double-strand break site)?
  - Or do we score the guide RNA sequence itself as a variant?
- **Spacer efficacy endpoint**: We have `/api/design/predict_crispr_spacer_efficacy` - should we use this instead of raw Evo2 endpoints?

**Impact**: Determines which Evo2 endpoint(s) to use and how to interpret scores.

---

### 3. Activity Score Interpretation

**Issue**: Doench dataset provides `score_drug_gene_rank` (0-1 normalized percentile). Need to understand relationship to cutting efficiency.

**Questions**:
- What does `score_drug_gene_rank` represent exactly?
  - Is it directly proportional to cutting efficiency (% indels)?
  - Or is it a percentile rank (top 20% = active)?
- How should we interpret the correlation?
  - Spearman ρ = 0.60 means what in practical terms?
  - What's the expected correlation for a "good" method?
- Rule Set 2 baseline: What's the exact Spearman correlation of Rule Set 2 on this dataset?
  - Paper mentions ~0.60-0.65, but need exact number for comparison

**Impact**: Determines success criteria and how to present results.

---

### 4. Dataset Selection and Splitting

**Issue**: Doench dataset has 4,390 guides total (1,841 FC + 2,549 RES). Need to decide which subset to use.

**Questions**:
- Should we benchmark on:
  - Full dataset (4,390 guides)?
  - FC dataset only (1,841 guides) - matches original paper focus?
  - Both separately (report FC and RES correlations separately)?
- For supervised training (Evo2 embeddings):
  - Train/test split: 80/20 or 70/30?
  - Should we stratify by gene to avoid gene-specific bias?
  - Cross-validation: Should we use k-fold or hold-out?
- Cell line specificity: Doench used HEK293T, A375, K562. Our metastasis guides target different contexts. Does this matter?

**Impact**: Affects benchmark methodology and results interpretation.

---

### 5. Baseline Comparisons

**Issue**: Need to compare Evo2 to established methods (GC-only, CHOPCHOP, Rule Set 2).

**Questions**:
- **GC-only baseline**: How should we compute this?
  - Simple GC content (0.4 ≤ GC ≤ 0.6 → high efficacy)?
  - Or use existing GC-based scoring formula?
- **CHOPCHOP**: 
  - Do we have CHOPCHOP scores pre-computed for Doench dataset?
  - Or do we need to run CHOPCHOP ourselves?
  - Which CHOPCHOP version/parameters?
- **Rule Set 2 (Azimuth)**:
  - The CSV already includes Rule Set 2 predictions - should we use those?
  - Or recompute for consistency?
- **Other baselines**: Should we include other methods (e.g., DeepSpCas9, DeepHF)?

**Impact**: Determines which comparisons are feasible and meaningful.

---

### 6. Guide Sequence Format

**Issue**: Doench provides 30mer (4bp upstream + 23bp guide + 3bp PAM). Need to extract correct sequence for Evo2.

**Questions**:
- For Evo2 scoring, should we use:
  - 23bp guide sequence only?
  - 30mer (with upstream context and PAM)?
  - Extended context (±150bp around target site)?
- PAM handling: Do we include PAM in Evo2 scoring, or remove it?
- Upstream context: The 4bp upstream - is this needed for Evo2, or just the 23bp guide?

**Impact**: Affects how we prepare sequences for Evo2 scoring.

---

### 7. Evo2 Endpoint Selection

**Issue**: Multiple Evo2 endpoints available. Need to determine which to use.

**Questions**:
- Should we use:
  - `/score_delta` (zero-shot variant effect)?
  - `/score_variant` (with Ensembl context fetch)?
  - `/api/design/predict_crispr_spacer_efficacy` (our custom endpoint)?
- For supervised approach:
  - Which Evo2 block/layer for embeddings? (Paper suggests block 20 for BRCA1)
  - Should we test multiple blocks and report best?
- Model size: 7B or 40B? (40B better but slower)

**Impact**: Determines implementation approach and expected performance.

---

### 8. Computational Resources

**Issue**: Evo2 scoring for 4,390 guides may be computationally intensive.

**Questions**:
- What's the expected runtime for 4,390 guides?
  - Zero-shot: ~X minutes/hours?
  - Supervised (embeddings): ~X minutes/hours?
- Do we have sufficient Evo2 service capacity?
- Should we batch process or run sequentially?
- Rate limiting: Any API rate limits we need to respect?

**Impact**: Affects timeline and feasibility.

---

### 9. Results Presentation

**Issue**: Need to determine how to present benchmark results in manuscript.

**Questions**:
- Should we report:
  - Single correlation number (best method)?
  - Both zero-shot and supervised (with train/test split)?
  - Per-gene correlations (to show consistency)?
- Statistical testing:
  - How do we test if Evo2 is significantly better than baselines?
  - Bootstrap confidence intervals?
  - Permutation testing?
- Figure/Table format:
  - Scatter plot (Evo2 predictions vs experimental)?
  - Correlation comparison bar chart?
  - Per-gene breakdown?

**Impact**: Affects manuscript presentation and reviewer evaluation.

---

### 10. Validation Against Our Guides

**Issue**: Doench guides target different genes (CD33, CD13, etc.) than our metastasis genes.

**Questions**:
- Does performance on Doench dataset predict performance on our metastasis guides?
- Should we also benchmark on a subset of our own guides (if we have experimental data)?
- How do we address potential domain shift (Doench genes vs metastasis genes)?

**Impact**: Affects generalizability claims in manuscript.

---

## Priority Ranking

### **P0 (Blocking Implementation)**:
1. Genomic context mapping (#1)
2. Evo2 scoring approach (#2)
3. Guide sequence format (#6)

### **P1 (Affects Results Quality)**:
4. Activity score interpretation (#3)
5. Dataset selection and splitting (#4)
6. Baseline comparisons (#5)

### **P2 (Optimization)**:
7. Evo2 endpoint selection (#7)
8. Computational resources (#8)
9. Results presentation (#9)
10. Validation against our guides (#10)

---

## Recommended Next Steps

1. **Immediate**: Resolve genomic context mapping (#1) - this is blocking
2. **Week 1**: Clarify Evo2 scoring approach (#2) and guide sequence format (#6)
3. **Week 2**: Implement benchmark with resolved approach
4. **Week 3**: Analyze results and prepare manuscript figures

---

**Status**: Awaiting manager input on P0 questions before proceeding with implementation.
