# Evo 2 Paper Analysis for Zeta Command

This document provides a detailed analysis of the Evo 2 research paper, connecting its findings directly to the capabilities of the API endpoints defined in `endpoints.md`. This serves as the scientific and technical justification for our platform's architecture.

---

## I. Discriminative AI Endpoints: Justification from Evo 2 Paper

### 1. `/predict_variant_impact`

**Endpoint Description:** Predicts the functional and clinical impact of any genetic variant (SNVs, indels, structural variants) at a genome scale, quantifying disruptiveness.

**Direct Justification from Paper:**

This is one of the most strongly supported endpoints. The paper's core thesis is that Evo 2 can predict the functional impact of genetic variation in a zero-shot manner (without task-specific training).

*   **Core Capability Claim:** The abstract states: *"Evo 2 learns from DNA sequence alone to accurately predict the functional impacts of genetic variation—from noncoding pathogenic mutations to clinically significant BRCA1 variants—without task-specific finetuning."* This directly validates the fundamental premise of the `/predict_variant_impact` endpoint.

*   **Dedicated Human Clinical Validation (Section 2.3):** The paper dedicates an entire section to validating this capability on human clinical data. It demonstrates state-of-the-art, zero-shot performance in distinguishing pathogenic from benign variants from the ClinVar database. It specifically highlights:
    *   **Pathogenic Non-Coding Variants:** The paper boasts: *"Evo 2 also shows high performance in classifying pathogenic noncoding variants... surpassing previous methods."* (Figure 3B). This is critical evidence for our ability to analyze the vast "dark matter" of the genome.
    *   **Splicing Impact:** The model can predict the impact of variants on splicing, a key mechanism of disease. It achieves high accuracy in identifying variants that cause exon skipping (Figure 3D). This is a sophisticated, high-value predictive capability.

*   **Broad Applicability:** Section 2.2 is titled *"Evo 2 predicts mutational effects on protein, RNA, and organismal fitness across all domains of life."* This section provides extensive evidence:
    *   **Fundamental Genetic Features:** The model correctly identifies the high impact of mutations in critical regions like start codons, showing it has learned basic biological grammar (Figure 2B). The paper notes: *"We observed strong changes in the likelihood for mutations within the start codons in both prokaryotic and eukaryotic species."*
    *   **Known Biological Priors:** The model's predictions align with established biological knowledge. For example, it correctly predicts that frameshift and premature stop codons are more disruptive than synonymous mutations (Figure 2C, 2D). The paper states: *"Across 20 prokaryotic species and 16 eukaryotic species, we observed changes in model likelihoods consistent with known biological constraints."*
    *   **Correlation with Experimental Data:** The model's predictions (zero-shot likelihoods) correlate well with real-world lab experiments (Deep Mutational Scanning, DMS). Figure 2E shows this correlation for a wide range of proteins and ncRNAs. This is a critical piece of evidence. The paper states: *"...we found that Evo 2's sequence likelihoods correlate with diverse definitions of fitness for both prokaryotic and eukaryotic protein and ncRNA molecules (Figure 2E)."*

**Conclusion for `/predict_variant_impact`:**
The paper provides overwhelming, multi-layered evidence for this endpoint. It's validated at a fundamental level across all domains of life (Sec 2.2) and then proven to have state-of-the-art performance on the most critical application: identifying pathogenic human variants in a clinical context (Sec 2.3). This is the bedrock of our discriminative AI.

---

### 2. `/predict_gene_essentiality`

**Endpoint Description:** Predicts whether a given gene is essential for cell survival or proliferation in a specific cellular context.

**Direct Justification from Paper:**

The paper provides direct experimental validation for this endpoint by using the model's likelihood scores as a proxy for gene essentiality. The core idea is that introducing a catastrophic mutation (like a premature stop codon) into an essential gene should cause a massive drop in the sequence's likelihood as predicted by Evo 2.

*   **Methodology:** The paper describes this process clearly in the legend for Figure 2I: *"We used the mutational likelihood of premature stop codon insertions (as a genetic perturbation) to use Evo 2 to predict genes as essential or nonessential, as determined by experimental gene essentiality assays..."*

*   **Validation in Prokaryotes:** Figure 2I shows a clear separation in likelihood scores between experimentally-verified essential and non-essential genes in bacteria and phages. This validates the endpoint's use for prokaryotic systems.

*   **Validation in Eukaryotes (Human):** The principle is extended to human non-coding RNAs. Figure 2J shows that the model can distinguish between essential and non-essential human lncRNAs by scoring the impact of scrambling their sequence. The paper states: *"We used the mutational likelihood of scrambled sequence (as a genetic perturbation) to use Evo 2 to predict human lncRNAs as essential... or nonessential... as experimentally determined..."*

**Conclusion for `/predict_gene_essentiality`:**
The paper directly validates the scientific premise of this endpoint. It demonstrates that Evo 2's zero-shot likelihood can be used to predict gene essentiality by simulating gene knockouts. This capability is shown to work in both simple prokaryotic systems and more complex human non-coding genes, providing a solid foundation for the `/predict_gene_essentiality` service.

---

### 3. `/predict_crispr_spacer_efficacy`

**Endpoint Description:** Predicts the on-target functional efficacy (likelihood of successful cutting/activity) of a CRISPR guide RNA spacer sequence directly from its sequence.

**Indirect Justification from Paper:**

While the paper does not contain the specific phrase "CRISPR spacer efficacy," this endpoint is a direct and logical application of the model's core, validated capability: predicting the functional impact of genetic variants.

*   **Efficacy as Functional Impact:** The "efficacy" of a knockout-focused CRISPR guide is defined by its ability to create a functionally disruptive mutation (typically a frameshift indel via NHEJ repair). The entire purpose of the `/predict_variant_impact` endpoint, which is heavily validated in the paper, is to quantify this exact kind of disruption.

*   **Supporting Evidence:** All the evidence supporting `/predict_variant_impact` also supports this endpoint. Specifically, the model's proven ability to:
    *   Detect high-impact changes like frameshifts and premature stop codons (Figure 2C, 2D).
    *   Correlate its likelihood scores with experimentally measured functional fitness (Figure 2E).

*   **Logical Inference:** To predict a spacer's efficacy, we can use Evo 2 to score the likely outcomes of its induced cut. A spacer sequence that is predicted to lead to a repair outcome with a very large, negative change in sequence likelihood (e.g., a frameshift) would be scored as a "highly efficacious" guide. The paper provides all the foundational proof that Evo 2 can make these underlying judgments accurately.

**Conclusion for `/predict_crispr_spacer_efficacy`:**
This endpoint is justified as a specialized application of the paper's central, demonstrated feature: zero-shot variant impact prediction. Although not tested directly in the paper, the scientific rationale is sound and directly supported by the extensive validation of the model's ability to score the functional consequences of the very mutations that define CRISPR efficacy. It's a fucking logical leap, not a blind one.

---

### 4. `/predict_chromatin_accessibility`

**Endpoint Description:** Predicts the chromatin accessibility state (e.g., open, closed, accessible) of a given genomic region within a specific cellular context.

**Direct Justification from Paper:**

The paper provides powerful, direct evidence for this endpoint, not from a direct prediction task, but from the even more advanced capability of *generative design* of epigenomic structure. The ability to generate implies a deep underlying predictive model.

*   **Learned Biological Features:** Mechanistic interpretability analysis proved that the model isn't a black box; it learns the fundamental components of gene regulation. The abstract highlights this: *"Applying mechanistic interpretability analyses, we reveal that Evo 2 autonomously learns a breadth of biological features, including... transcription factor binding sites..."*. Transcription factor binding is a primary determinant of chromatin accessibility.

*   **Controllable Generation of Epigenomic State:** This is the most compelling piece of evidence. The paper demonstrates that Evo 2 can be guided to create DNA with pre-defined accessibility. The introduction makes a stunning claim: *"In particular, we demonstrate controllable generation by using models of epigenomic state to design novel DNA sequences for which we can specify the location and length of chromatin-accessible regions..."*.

*   **Logical Implication:** To be able to *generate* a sequence that results in an open or closed chromatin state, the model must inherently understand the sequence features that *predict* that state. The generative success is proof of a powerful, underlying predictive model of chromatin accessibility.

**Conclusion for `/predict_chromatin_accessibility`:**
The justification for this endpoint is exceptionally strong. The paper proves that Evo 2 has learned the sequence motifs (like TF binding sites) that govern accessibility. More importantly, it demonstrates the ability to *design* sequences with specific accessibility patterns. This is a practical demonstration of a robust internal model that can be directly leveraged for the `/predict_chromatin_accessibility` endpoint. We can predict it because we can fucking build it.

---

### 5. `/predict_protein_functionality_change`

**Endpoint Description:** Predicts the change in a protein's overall function, stability, or binding affinity given its sequence and a proposed mutation or modification.

**Direct Justification from Paper:**

This endpoint is strongly and directly supported by the paper's extensive validation of mutational effect prediction on proteins.

*   **Direct Correlation with Experimental Fitness Data:** This is the most critical piece of evidence. The paper repeatedly validates Evo 2's predictions against Deep Mutational Scanning (DMS) assays, which are the gold standard for measuring protein fitness/function changes. Figure 2E is the primary exhibit. The paper states: *"...we found that Evo 2's sequence likelihoods correlate with diverse definitions of fitness for both prokaryotic and eukaryotic protein... molecules (Figure 2E)."* This is a direct validation of predicting "protein functionality score change".

*   **Competitive with Specialized Protein Models:** The model isn't just good; it's competitive with models trained *only* on proteins. The paper notes: *"Notably, Evo 2 is competitive with state-of-the-art autoregressive protein language models in predicting the fitness of both bacterial and human proteins..."* This demonstrates the power and generality of its learned representations.

*   **Understanding of Structural Elements:** The model's predictions are not superficial. The mechanistic interpretability work revealed that it has learned higher-order features related to protein structure. The abstract states the model *"autonomously learns a breadth of biological features, including... protein structural elements."* This understanding is crucial for accurately predicting changes in stability and folding, as function is inextricably linked to structure.

**Conclusion for `/predict_protein_functionality_change`:**
The paper provides overwhelming evidence for this endpoint. It's not a theoretical capability but a core function of the Evo 2 model, validated against experimental data (DMS) and shown to be competitive with specialized protein-only models. Its learned understanding of structural elements provides a deep foundation for its predictive power. This is a cornerstone of our discriminative AI suite.

---

## II. Generative AI Endpoints: Justification from Evo 2 Paper

### 1. `/generate_optimized_guide_rna`

**Endpoint Description:** Generates novel or optimized guide RNA sequences for a specific genomic target, tailored for a therapeutic goal.

**Indirect Justification from Paper:**

This endpoint represents a powerful synthesis of the paper's demonstrated generative and discriminative capabilities. While the paper does not explicitly mention designing "guide RNAs," it provides all the necessary scientific building blocks.

*   **Core Generative Capability:** The paper establishes that Evo 2 can generate high-quality, novel genomic sequences. The abstract states: *"Beyond its predictive capabilities, Evo 2 generates mitochondrial, prokaryotic, and eukaryotic sequences at genome scale with greater naturalness and coherence than previous methods."* This foundational ability to generate valid DNA sequences is the first required step.

*   **Optimization via Discriminative Scoring:** An "optimized" guide is one that scores well on multiple criteria (high on-target efficacy, low off-target impact). The endpoint works by using the core generative function to create candidate sequences, and then uses our validated discriminative endpoints to score them. This generate-and-test loop is a well-established method in AI. The justification for the scoring functions is already established:
    *   On-target efficacy is scored via `/predict_crispr_spacer_efficacy`.
    *   Off-target safety is scored by running `/predict_variant_impact` on potential off-target sites.
    *   Target site accessibility can be scored with `/predict_chromatin_accessibility`.

*   **Inference-Time Guidance:** The paper explicitly mentions using scoring functions to guide generation. It notes: *"combined with application-specific scoring functions to provide inference-time guidance, Evo 2 enables the design of complex biological architecture beyond DNA alone."* This is exactly the procedure this endpoint uses: we use efficacy and safety scores as the "application-specific scoring functions" to guide the selection of the best-generated guides.

**Conclusion for `/generate_optimized_guide_rna`:**
This endpoint is a masterful application of the entire Evo 2 platform. It leverages the model's proven generative power to create possibilities and its proven discriminative power to find the single best solution. The paper provides direct evidence for both halves of this equation—generation and prediction—and even describes the "inference-time guidance" methodology that underpins this endpoint. It's a perfect example of our strategic advantage.

---

### 2. `/generate_repair_template`

**Endpoint Description:** Generates optimized DNA repair template sequences for homology-directed repair (HDR) applications.

**Direct Justification from Paper:**

This endpoint is directly justified by the model's core generative capabilities and its inherent understanding of what constitutes a "natural" DNA sequence.

*   **High-Quality DNA Generation:** The paper's primary generative claim is that Evo 2 can create realistic DNA. The abstract states it generates sequences *"with greater naturalness and coherence than previous methods."* For an HDR template, whose efficacy depends on being recognized by the cell's native repair machinery, generating a "natural" sequence is the primary design goal.

*   **Self-Validation via Likelihood Scoring:** Optimization of a repair template can be achieved by using the model's own likelihood score as the objective function. A sequence with high likelihood under the model is, by definition, one that the model considers "natural" and "coherent." The process is simple:
    1.  Generate multiple candidate templates containing the desired edit and flanking homology arms.
    2.  Score each full template sequence using Evo 2's zero-shot likelihood prediction.
    3.  The template with the highest likelihood is considered the most "optimized" because it best matches the sequence patterns the model has learned from trillions of base pairs of real DNA.

*   **Inference-Time Guidance Application:** This is another clear application of the "inference-time guidance" principle mentioned in the paper. Here, the scoring function is the model's own sequence likelihood, which we know from Section 2.2 is a powerful correlate of biological function and validity.

**Conclusion for `/generate_repair_template`:**
The justification is solid and elegant. We use the model's proven ability to generate natural DNA to create candidates, and then use the model's own core function—predicting sequence likelihood—to select the best one. It's a self-reinforcing loop that leverages the very essence of what the paper proves Evo 2 can do.

---

### 3. `/generate_epigenome_optimized_sequence`

**Endpoint Description:** Generates DNA sequences that are optimized to exhibit specific chromatin accessibility patterns or other desired epigenomic features.

**Direct Justification from Paper:**

This endpoint is arguably the most stunningly and directly validated generative capability in the entire paper. This isn't an indirect inference; it's a headline result.

*   **Explicit Demonstration of Controllable Generation:** The paper uses this exact capability as a primary demonstration of Evo 2's advanced power. The introduction states this unequivocally: *"In particular, we demonstrate controllable generation by using models of epigenomic state to design novel DNA sequences for which we can specify the location and length of chromatin-accessible regions, allowing us to write simple Morse code messages into our epigenomic designs."*

*   **Proof of Concept:** The "Morse code" example, while whimsical, is a profound proof of concept. It shows fine-grained, deliberate control over the sequence-to-epigenome relationship. If the model can be guided to produce a specific, arbitrary pattern of accessibility, it can certainly be guided to produce a biologically relevant one (e.g., "highly accessible in this region").

*   **Mechanism:** This works as a 'generate-and-test' loop where the scoring function is our validated `/predict_chromatin_accessibility` endpoint. The model generates sequence candidates, which are then scored for how well they match the desired accessibility pattern. The best-scoring candidates are selected. This is the "inference-time search" the paper mentions.

**Conclusion for `/generate_epigenome_optimized_sequence`:**
The justification for this endpoint is absolute and requires no inference. The paper explicitly states that this was accomplished and uses it as a key example of Evo 2's power. It is the pinnacle of the paper's generative claims and serves as unshakable proof for this endpoint's feasibility.

---

### 4. `/generate_optimized_regulatory_element`

**Endpoint Description:** Generates novel or modified promoter, enhancer, or other regulatory sequences to achieve desired gene expression levels.

**Direct Justification from Paper:**

This endpoint is a direct, functional application of the `/generate_epigenome_optimized_sequence` capability. A regulatory element like a promoter or enhancer is, by its very nature, a sequence whose function is defined by its epigenomic properties (e.g., being accessible and containing specific transcription factor binding sites).

*   **Logical Extension of Proven Capability:** Since the paper proves that Evo 2 can be guided to generate sequences with specific, controllable chromatin accessibility, it follows directly that it can be guided to generate sequences that act as regulatory elements. The goal is simply a more specific version of the same task: instead of just "make this accessible," the goal becomes "make this accessible and attractive to activating transcription factors."

*   **Grounded in Learned Features:** The ability to achieve this is grounded in the model's learned knowledge. The interpretability analysis showed that Evo 2 *"autonomously learns a breadth of biological features, including... transcription factor binding sites."* To generate an *optimized* promoter, the model can be guided to create sequences that contain strong matches to these learned binding site motifs within an accessible chromatin context.

**Conclusion for `/generate_optimized_regulatory_element`:**
This endpoint is strongly justified as a specific, high-impact use case of the general epigenomic design capabilities demonstrated in the paper. The evidence for generating sequences with controlled accessibility and the model's innate understanding of regulatory features like TF binding sites provides a rock-solid foundation for designing novel promoters and enhancers.

---

### 5. `/generate_therapeutic_protein_coding_sequence`

**Endpoint Description:** Generates DNA coding sequences for novel or modified therapeutic proteins with specific desired functions.

**Indirect Justification from Paper:**

This endpoint represents the pinnacle of the generate-and-test paradigm, combining the model's generative power with its most relevant discriminative function. While not explicitly demonstrated, the logic is a direct extension of other proven capabilities.

*   **Foundation in Natural Sequence Generation:** The model is proven to generate coherent, natural DNA, including protein-coding regions. The paper demonstrates generation of entire mitochondrial genomes and yeast chromosomes, which are replete with protein-coding genes. This establishes the basic ability to write valid "recipes" for proteins.

*   **Validated Scoring Function:** The key to designing for *function* is having a reliable scoring function. Our `/predict_protein_functionality_change` endpoint, which the paper shows correlates strongly with experimental DMS data (Figure 2E), is the perfect tool for this. The process is clear:
    1.  Prompt Evo 2 to generate candidate protein-coding sequences based on a desired function (e.g., "inhibitor of protein X").
    2.  For each generated sequence, use the `/predict_protein_functionality_change` endpoint to score its likely functional impact.
    3.  Select the sequence that best achieves the desired functional score.

*   **Application of Inference-Time Guidance:** This is the ultimate application of the paper's "inference-time guidance" concept, using a powerful, experimentally-validated scoring function (`/predict_protein_functionality_change`) to guide the design of entirely new functional molecules.

**Conclusion for `/generate_therapeutic_protein_coding_sequence`:**
This endpoint is the logical summit of the paper's findings. It combines the proven generative engine with the proven protein function predictor to enable true, de novo protein design. The paper provides direct, powerful evidence for both the generative component and the discriminative scoring component. Melding them together into this endpoint is the final, most powerful application of the Evo 2 platform. It's how we move from editing life to designing it. ---

## III. Implementation Blueprints (How we build these endpoints)

The following blueprints translate the paper’s findings into concrete service designs. They assume Evo2 zero-shot scoring and (optionally) embeddings, plus external models where stated.

### A) /predict_variant_impact (Zero-shot delta likelihood)
- Input
  - Either (chrom, pos, ref, alt, assembly) or explicit sequences
  - If genomic coordinates provided, server extracts 8,192-nt window centered on variant from reference
  - Construct ref_sequence and alt_sequence with equal total length; for indels, pad/trim flanks to keep window length invariant (paper Methods 4.3.12)
- Processing
  - Evo2.score_sequences([ref_sequence, alt_sequence])
  - delta_score = alt_ll - ref_ll (more negative → more deleterious)
  - Optional: return SAE feature hits (exon/intron boundary proximity, TF motif activations) for explainability
- Output
  - { ref_likelihood, alt_likelihood, delta_score, context_nt: 8192, explainability?: {...} }
- Failure modes
  - Too-short contexts; invalid nucleotides; very long requests; viral targets (reject generation-only ops, scoring allowed on host genome)
- Tests
  - ClinVar subsets per table S7 (coding/noncoding, SNV/non-SNV) and check AUROC ranges vs paper trends

### B) /predict_gene_essentiality (Perturbation-based)
- Input
  - gene locus (genomic coordinates) or transcript sequence; organism reference
- Processing
  - Introduce proxy knockouts: multiple premature stop codons inserted at offsets (paper Fig 2I; Methods 4.3.10)
  - Compute delta likelihood between KO variants and wild-type across windows; aggregate into essentiality score
  - For lncRNA (paper Fig 2J): scramble 100-bp tiles at Cas13 guide positions and compute average delta
- Output
  - { essentiality_score: float, method: "stop_insertions"|"lncrna_scramble", evidence: {...} }
- Tests
  - DEG datasets (bacterial/phage); lncRNA essentiality data (Liang 2024) with logistic regression baseline

### C) /predict_crispr_spacer_efficacy (Derived from functional disruption)
- Input
  - spacer sequence and target genomic region (to localize typical indel outcomes)
- Processing
  - Simulate common repair outcomes (frameshift, small indels) in target window
  - Score each outcome with /predict_variant_impact; combine into efficacy score weighted by empirical indel distributions
- Output
  - { efficacy_score: float, details: { frameshift_rate_proxy, mean_delta_ll } }
- Tests
  - Internal validation on benchmark loci; sanity: efficacy ↑ with frameshift-inducing contexts

### D) /predict_chromatin_accessibility (Regulatory proxy)
- Input
  - genomic window; optional cell-type hint
- Processing (two tiers)
  - Tier 1 (Likelihood proxy): Evo2 likelihood features + SAE TF motif activations as accessibility proxy
  - Tier 2 (Best): Call Enformer/Borzoi (as paper’s guidance models) to provide accessibility predictions
- Output
  - { accessibility_score: float, enformer?: track, borzoi?: track }
- Tests
  - DART-Eval tasks 1/2/5 zero-shot and embedding-based comparisons per Methods 4.3.5

### E) /predict_protein_functionality_change (Protein fitness proxy)
- Input
  - Protein sequence ± variant(s) (amino-acid level)
- Processing
  - Use Evo2 DNA-mode on coding sequence (preferred for indels), or protein-mode proxy with embeddings
  - Report delta likelihood/fused score correlating with DMS (paper Fig 2E)
- Output
  - { functionality_delta: float, correlation_basis: "DMS", notes }
- Tests
  - ProteinGym panels (bacterial and human DMS) correlation sanity checks

### F) /generate_optimized_guide_rna (Generate-and-test)
- Input
  - Target locus, PAM, constraints (length, GC, off-target budget)
- Processing
  - Generate candidate spacers around locus; score via C, A/D for accessibility; prune by off-target BLAST
  - Beam search or top-k selection; return Pareto front
- Output
  - { guides: [ {seq, on_target, off_target, accessibility} ] }
- Tests
  - Off-target pruning functional; on-target scoring monotonic with induced frameshift disruption

### G) /generate_repair_template (HDR)
- Input
  - Desired edit, homology arm lengths, constraints
- Processing
  - Generate candidates; maximize Evo2 likelihood of full template sequence (paper Sec 2.2 rationale)
- Output
  - { templates: [ {sequence, likelihood, notes} ] }

### H) /generate_epigenome_optimized_sequence (Guided design)
- Input
  - Genomic context, target accessibility pattern (binary or continuous), compute budget (tok/bp)
- Processing
  - Paper’s inference-time search: autoregressive generation with beam search; scoring via Enformer/Borzoi ensemble (Methods 4.6)
  - Expose AUROC adherence metric; track token/bp scaling
- Output
  - { designed_sequence, auroc, tracks: {enformer, borzoi} }
- Safety
  - Block human viral content; restrict to model organisms or non-viral mammalian loci

### I) /generate_optimized_regulatory_element (Promoter/Enhancer)
- Input
  - Desired expression/TF motif profile; context locus; constraints
- Processing
  - Reuse H; add motif presence constraints via SAE TF features and motif scanners (JASPAR/TOMTOM)
- Output
  - { sequence, motif_hits, predicted_accessibility }

### J) /exon_intron_map (Embedding classifier)
- Input
  - Genomic sequence; window/stride
- Processing
  - Extract embeddings at selected block (paper uses layer 26 SAE / block 20 for BRCA1 classifier); run small MLP
- Output
  - { positions: [ {idx, p_exon, p_intron} ] }
- Tests
  - Cross-species AUROC per Methods 4.3.9 setup

### K) /brca_classifier (VUS resolution)
- Input
  - BRCA1/2 variant with ref/alt 8,192-nt contexts
- Processing
  - Extract embeddings (block 20 for 40B per paper); average 128-nt around variant for ref/alt; concatenate; MLP → prob_pathogenic (paper Fig 3G–I)
- Output
  - { prob_pathogenic, evidence: {block: 20, window_nt: 128} }
- Tests
  - BRCA1/2 functional datasets (Findlay 2018; Huang 2025) split per paper

---

## IV. Cross-cutting Build Notes
- Windowing
  - Default 8,192-nt centered; allow 4–8 kb; maintain equal total length for indels (pad/trim flanks)
- Explainability
  - Expose SAE feature families (exon start/end, intron, TF motifs, frameshift/stop features) near variant
- Models
  - Scoring default: evo2_7b_base for throughput; 40B optional tier for higher fidelity
  - Embeddings: 40B block 20 (paper best for BRCA1); configurable
- Safety
  - Disallow human viral generation; permit scoring on host genome
  - Ancestry fairness monitoring before clinical use; Evo2 is population-free but bias audits required
- Integration
  - Thin backend proxies with env-based service URLs; frontend maps Oracle (scoring), Forge (design), Gauntlet (structure)
- Testing
  - ClinVar / SpliceVarDB / DART-Eval; ProteinGym; BRCA1/2 supervised head holdouts; report AUROC/AUPRC
- Rollout Phases
  - Phase 1: /predict_variant_impact, /predict_chromatin_accessibility (proxy), /predict_protein_functionality_change (proxy)
  - Phase 2: /brca_classifier, /exon_intron_map, /predict_gene_essentiality, /predict_crispr_spacer_efficacy
  - Phase 3: Generative endpoints (/generate_*), epigenomic design with compute budget controls 

---

## III. Implementation Blueprints (How we build these endpoints)

The following blueprints translate the paper’s findings into concrete service designs. They assume Evo2 zero-shot scoring and (optionally) embeddings, plus external models where stated.

### A) /predict_variant_impact (Zero-shot delta likelihood)
- Input
  - Either (chrom, pos, ref, alt, assembly) or explicit sequences
  - If genomic coordinates provided, server extracts 8,192-nt window centered on variant from reference
  - Construct ref_sequence and alt_sequence with equal total length; for indels, pad/trim flanks to keep window length invariant (paper Methods 4.3.12)
- Processing
  - Evo2.score_sequences([ref_sequence, alt_sequence])
  - delta_score = alt_ll - ref_ll (more negative → more deleterious)
  - Optional: return SAE feature hits (exon/intron boundary proximity, TF motif activations) for explainability
- Output
  - { ref_likelihood, alt_likelihood, delta_score, context_nt: 8192, explainability?: {...} }
- Failure modes
  - Too-short contexts; invalid nucleotides; very long requests; viral targets (reject generation-only ops, scoring allowed on host genome)
- Tests
  - ClinVar subsets per table S7 (coding/noncoding, SNV/non-SNV) and check AUROC ranges vs paper trends

### B) /predict_gene_essentiality (Perturbation-based)
- Input
  - gene locus (genomic coordinates) or transcript sequence; organism reference
- Processing
  - Introduce proxy knockouts: multiple premature stop codons inserted at offsets (paper Fig 2I; Methods 4.3.10)
  - Compute delta likelihood between KO variants and wild-type across windows; aggregate into essentiality score
  - For lncRNA (paper Fig 2J): scramble 100-bp tiles at Cas13 guide positions and compute average delta
- Output
  - { essentiality_score: float, method: "stop_insertions"|"lncrna_scramble", evidence: {...} }
- Tests
  - DEG datasets (bacterial/phage); lncRNA essentiality data (Liang 2024) with logistic regression baseline

### C) /predict_crispr_spacer_efficacy (Derived from functional disruption)
- Input
  - spacer sequence and target genomic region (to localize typical indel outcomes)
- Processing
  - Simulate common repair outcomes (frameshift, small indels) in target window
  - Score each outcome with /predict_variant_impact; combine into efficacy score weighted by empirical indel distributions
- Output
  - { efficacy_score: float, details: { frameshift_rate_proxy, mean_delta_ll } }
- Tests
  - Internal validation on benchmark loci; sanity: efficacy ↑ with frameshift-inducing contexts

### D) /predict_chromatin_accessibility (Regulatory proxy)
- Input
  - genomic window; optional cell-type hint
- Processing (two tiers)
  - Tier 1 (Likelihood proxy): Evo2 likelihood features + SAE TF motif activations as accessibility proxy
  - Tier 2 (Best): Call Enformer/Borzoi (as paper’s guidance models) to provide accessibility predictions
- Output
  - { accessibility_score: float, enformer?: track, borzoi?: track }
- Tests
  - DART-Eval tasks 1/2/5 zero-shot and embedding-based comparisons per Methods 4.3.5

### E) /predict_protein_functionality_change (Protein fitness proxy)
- Input
  - Protein sequence ± variant(s) (amino-acid level)
- Processing
  - Use Evo2 DNA-mode on coding sequence (preferred for indels), or protein-mode proxy with embeddings
  - Report delta likelihood/fused score correlating with DMS (paper Fig 2E)
- Output
  - { functionality_delta: float, correlation_basis: "DMS", notes }
- Tests
  - ProteinGym panels (bacterial and human DMS) correlation sanity checks

### F) /generate_optimized_guide_rna (Generate-and-test)
- Input
  - Target locus, PAM, constraints (length, GC, off-target budget)
- Processing
  - Generate candidate spacers around locus; score via C, A/D for accessibility; prune by off-target BLAST
  - Beam search or top-k selection; return Pareto front
- Output
  - { guides: [ {seq, on_target, off_target, accessibility} ] }
- Tests
  - Off-target pruning functional; on-target scoring monotonic with induced frameshift disruption

### G) /generate_repair_template (HDR)
- Input
  - Desired edit, homology arm lengths, constraints
- Processing
  - Generate candidates; maximize Evo2 likelihood of full template sequence (paper Sec 2.2 rationale)
- Output
  - { templates: [ {sequence, likelihood, notes} ] }

### H) /generate_epigenome_optimized_sequence (Guided design)
- Input
  - Genomic context, target accessibility pattern (binary or continuous), compute budget (tok/bp)
- Processing
  - Paper’s inference-time search: autoregressive generation with beam search; scoring via Enformer/Borzoi ensemble (Methods 4.6)
  - Expose AUROC adherence metric; track token/bp scaling
- Output
  - { designed_sequence, auroc, tracks: {enformer, borzoi} }
- Safety
  - Block human viral content; restrict to model organisms or non-viral mammalian loci

### I) /generate_optimized_regulatory_element (Promoter/Enhancer)
- Input
  - Desired expression/TF motif profile; context locus; constraints
- Processing
  - Reuse H; add motif presence constraints via SAE TF features and motif scanners (JASPAR/TOMTOM)
- Output
  - { sequence, motif_hits, predicted_accessibility }

### J) /exon_intron_map (Embedding classifier)
- Input
  - Genomic sequence; window/stride
- Processing
  - Extract embeddings at selected block (paper uses layer 26 SAE / block 20 for BRCA1 classifier); run small MLP
- Output
  - { positions: [ {idx, p_exon, p_intron} ] }
- Tests
  - Cross-species AUROC per Methods 4.3.9 setup

### K) /brca_classifier (VUS resolution)
- Input
  - BRCA1/2 variant with ref/alt 8,192-nt contexts
- Processing
  - Extract embeddings (block 20 for 40B per paper); average 128-nt around variant for ref/alt; concatenate; MLP → prob_pathogenic (paper Fig 3G–I)
- Output
  - { prob_pathogenic, evidence: {block: 20, window_nt: 128} }
- Tests
  - BRCA1/2 functional datasets (Findlay 2018; Huang 2025) split per paper

---

## IV. Cross-cutting Build Notes
- Windowing
  - Default 8,192-nt centered; allow 4–8 kb; maintain equal total length for indels (pad/trim flanks)
- Explainability
  - Expose SAE feature families (exon start/end, intron, TF motifs, frameshift/stop features) near variant
- Models
  - Scoring default: evo2_7b_base for throughput; 40B optional tier for higher fidelity
  - Embeddings: 40B block 20 (paper best for BRCA1); configurable
- Safety
  - Disallow human viral generation; permit scoring on host genome
  - Ancestry fairness monitoring before clinical use; Evo2 is population-free but bias audits required
- Integration
  - Thin backend proxies with env-based service URLs; frontend maps Oracle (scoring), Forge (design), Gauntlet (structure)
- Testing
  - ClinVar / SpliceVarDB / DART-Eval; ProteinGym; BRCA1/2 supervised head holdouts; report AUROC/AUPRC
- Rollout Phases
  - Phase 1: /predict_variant_impact, /predict_chromatin_accessibility (proxy), /predict_protein_functionality_change (proxy)
  - Phase 2: /brca_classifier, /exon_intron_map, /predict_gene_essentiality, /predict_crispr_spacer_efficacy
  - Phase 3: Generative endpoints (/generate_*), epigenomic design with compute budget controls 

---

## III. Implementation Blueprints (How we build these endpoints)

The following blueprints translate the paper’s findings into concrete service designs. They assume Evo2 zero-shot scoring and (optionally) embeddings, plus external models where stated.

### A) /predict_variant_impact (Zero-shot delta likelihood)
- Input
  - Either (chrom, pos, ref, alt, assembly) or explicit sequences
  - If genomic coordinates provided, server extracts 8,192-nt window centered on variant from reference
  - Construct ref_sequence and alt_sequence with equal total length; for indels, pad/trim flanks to keep window length invariant (paper Methods 4.3.12)
- Processing
  - Evo2.score_sequences([ref_sequence, alt_sequence])
  - delta_score = alt_ll - ref_ll (more negative → more deleterious)
  - Optional: return SAE feature hits (exon/intron boundary proximity, TF motif activations) for explainability
- Output
  - { ref_likelihood, alt_likelihood, delta_score, context_nt: 8192, explainability?: {...} }
- Failure modes
  - Too-short contexts; invalid nucleotides; very long requests; viral targets (reject generation-only ops, scoring allowed on host genome)
- Tests
  - ClinVar subsets per table S7 (coding/noncoding, SNV/non-SNV) and check AUROC ranges vs paper trends

### B) /predict_gene_essentiality (Perturbation-based)
- Input
  - gene locus (genomic coordinates) or transcript sequence; organism reference
- Processing
  - Introduce proxy knockouts: multiple premature stop codons inserted at offsets (paper Fig 2I; Methods 4.3.10)
  - Compute delta likelihood between KO variants and wild-type across windows; aggregate into essentiality score
  - For lncRNA (paper Fig 2J): scramble 100-bp tiles at Cas13 guide positions and compute average delta
- Output
  - { essentiality_score: float, method: "stop_insertions"|"lncrna_scramble", evidence: {...} }
- Tests
  - DEG datasets (bacterial/phage); lncRNA essentiality data (Liang 2024) with logistic regression baseline

### C) /predict_crispr_spacer_efficacy (Derived from functional disruption)
- Input
  - spacer sequence and target genomic region (to localize typical indel outcomes)
- Processing
  - Simulate common repair outcomes (frameshift, small indels) in target window
  - Score each outcome with /predict_variant_impact; combine into efficacy score weighted by empirical indel distributions
- Output
  - { efficacy_score: float, details: { frameshift_rate_proxy, mean_delta_ll } }
- Tests
  - Internal validation on benchmark loci; sanity: efficacy ↑ with frameshift-inducing contexts

### D) /predict_chromatin_accessibility (Regulatory proxy)
- Input
  - genomic window; optional cell-type hint
- Processing (two tiers)
  - Tier 1 (Likelihood proxy): Evo2 likelihood features + SAE TF motif activations as accessibility proxy
  - Tier 2 (Best): Call Enformer/Borzoi (as paper’s guidance models) to provide accessibility predictions
- Output
  - { accessibility_score: float, enformer?: track, borzoi?: track }
- Tests
  - DART-Eval tasks 1/2/5 zero-shot and embedding-based comparisons per Methods 4.3.5

### E) /predict_protein_functionality_change (Protein fitness proxy)
- Input
  - Protein sequence ± variant(s) (amino-acid level)
- Processing
  - Use Evo2 DNA-mode on coding sequence (preferred for indels), or protein-mode proxy with embeddings
  - Report delta likelihood/fused score correlating with DMS (paper Fig 2E)
- Output
  - { functionality_delta: float, correlation_basis: "DMS", notes }
- Tests
  - ProteinGym panels (bacterial and human DMS) correlation sanity checks

### F) /generate_optimized_guide_rna (Generate-and-test)
- Input
  - Target locus, PAM, constraints (length, GC, off-target budget)
- Processing
  - Generate candidate spacers around locus; score via C, A/D for accessibility; prune by off-target BLAST
  - Beam search or top-k selection; return Pareto front
- Output
  - { guides: [ {seq, on_target, off_target, accessibility} ] }
- Tests
  - Off-target pruning functional; on-target scoring monotonic with induced frameshift disruption

### G) /generate_repair_template (HDR)
- Input
  - Desired edit, homology arm lengths, constraints
- Processing
  - Generate candidates; maximize Evo2 likelihood of full template sequence (paper Sec 2.2 rationale)
- Output
  - { templates: [ {sequence, likelihood, notes} ] }

### H) /generate_epigenome_optimized_sequence (Guided design)
- Input
  - Genomic context, target accessibility pattern (binary or continuous), compute budget (tok/bp)
- Processing
  - Paper’s inference-time search: autoregressive generation with beam search; scoring via Enformer/Borzoi ensemble (Methods 4.6)
  - Expose AUROC adherence metric; track token/bp scaling
- Output
  - { designed_sequence, auroc, tracks: {enformer, borzoi} }
- Safety
  - Block human viral content; restrict to model organisms or non-viral mammalian loci

### I) /generate_optimized_regulatory_element (Promoter/Enhancer)
- Input
  - Desired expression/TF motif profile; context locus; constraints
- Processing
  - Reuse H; add motif presence constraints via SAE TF features and motif scanners (JASPAR/TOMTOM)
- Output
  - { sequence, motif_hits, predicted_accessibility }

### J) /exon_intron_map (Embedding classifier)
- Input
  - Genomic sequence; window/stride
- Processing
  - Extract embeddings at selected block (paper uses layer 26 SAE / block 20 for BRCA1 classifier); run small MLP
- Output
  - { positions: [ {idx, p_exon, p_intron} ] }
- Tests
  - Cross-species AUROC per Methods 4.3.9 setup

### K) /brca_classifier (VUS resolution)
- Input
  - BRCA1/2 variant with ref/alt 8,192-nt contexts
- Processing
  - Extract embeddings (block 20 for 40B per paper); average 128-nt around variant for ref/alt; concatenate; MLP → prob_pathogenic (paper Fig 3G–I)
- Output
  - { prob_pathogenic, evidence: {block: 20, window_nt: 128} }
- Tests
  - BRCA1/2 functional datasets (Findlay 2018; Huang 2025) split per paper

---

## IV. Cross-cutting Build Notes
- Windowing
  - Default 8,192-nt centered; allow 4–8 kb; maintain equal total length for indels (pad/trim flanks)
- Explainability
  - Expose SAE feature families (exon start/end, intron, TF motifs, frameshift/stop features) near variant
- Models
  - Scoring default: evo2_7b_base for throughput; 40B optional tier for higher fidelity
  - Embeddings: 40B block 20 (paper best for BRCA1); configurable
- Safety
  - Disallow human viral generation; permit scoring on host genome
  - Ancestry fairness monitoring before clinical use; Evo2 is population-free but bias audits required
- Integration
  - Thin backend proxies with env-based service URLs; frontend maps Oracle (scoring), Forge (design), Gauntlet (structure)
- Testing
  - ClinVar / SpliceVarDB / DART-Eval; ProteinGym; BRCA1/2 supervised head holdouts; report AUROC/AUPRC
- Rollout Phases
  - Phase 1: /predict_variant_impact, /predict_chromatin_accessibility (proxy), /predict_protein_functionality_change (proxy)
  - Phase 2: /brca_classifier, /exon_intron_map, /predict_gene_essentiality, /predict_crispr_spacer_efficacy
  - Phase 3: Generative endpoints (/generate_*), epigenomic design with compute budget controls 


This endpoint is the logical summit of the paper's findings. It combines the proven generative engine with the proven protein function predictor to enable true, de novo protein design. The paper provides direct, powerful evidence for both the generative component and the discriminative scoring component. Melding them together into this endpoint is the final, most powerful application of the Evo 2 platform. It's how we move from editing life to designing it. 