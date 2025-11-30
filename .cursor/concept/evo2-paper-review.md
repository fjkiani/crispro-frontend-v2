# Evo 2 Paper — Section-by-Section Learnings

## Abstract
- **Scale and scope**: Evo 2 is trained on 9.3T base pairs with 7B/40B parameters and a 1M-token context at single-nucleotide resolution.
- **Zero-shot predictions**: Accurately predicts functional impacts across coding and noncoding variants (including BRCA1) without task-specific finetuning.
- **Mechanistic signals**: Learns exon–intron boundaries, TF motifs, protein structural elements, and prophage regions; validated via sparse autoencoders.
- **Generation**: Produces genome-scale sequences (mitochondrial, prokaryotic, eukaryotic) with higher naturalness and coherence than prior models.
- **Controllability**: Inference-time search enables controllable epigenomic structure; first inference-time scaling results in biology are demonstrated.
- **Open release**: Models, training/inference code, and the OpenGenome2 dataset are fully open to accelerate research.

## 1. Introduction
- **Motivation**: Building a generalist model for biology requires a representation that captures complexity across all domains of life, beyond human intuition, leveraging AI scaling laws and diverse data.
- **Diversity of genomes**: Prokaryotes vs eukaryotes differ dramatically (noncoding content, splicing, epigenomic layers), necessitating long-context learning.
- **Evo 2 proposal**: Train 7B/40B parameter models on 9.3T tokens with 1M-token context to learn cross-scale patterns and enable retrieval across long genomic distances.
- **Capabilities**: Zero-shot variant impact prediction across DNA/RNA/protein and domains of life; state-of-the-art on noncoding variant pathogenicity and strong BRCA1 classification with supervised embeddings.
- **Interpretability**: SAEs reveal features corresponding to exon/intron boundaries, TF motifs, protein structure, and mobile elements, enabling discovery tasks.
- **Generation & control**: Genome/chromosome-scale generation (mitochondrial, minimal bacteria, yeast chromosome) and inference-time guided designs for epigenomic structure (e.g., Morse code accessibility patterns).
- **Open science**: Release of model weights, code, and OpenGenome2 dataset to accelerate community work.

## 2.1. Evo 2 model architecture, training procedure, and data
- **Scaling & context**: 7B/40B models trained with two-phase strategy: short-context pretraining (8,192 tokens, genic-weighted) then multi-stage midtraining to 1M tokens for long-range dependencies.
- **Architecture**: StripedHyena 2 multi-hybrid (input-dependent convolutions + attention) delivers 1.3× throughput at 16k and ~3× at 1M vs optimized Transformers; better loss scaling than SH1/Transformers.
- **Context extension**: Rotary embeddings enable effective 1M-token context; successful recall in a 1M-length needle-in-a-haystack task.
- **Data**: OpenGenome2 curated dataset spanning bacteria, archaea, eukarya, phage (8.8T nt); excludes eukaryotic viruses (validated by high perplexity on excluded domain).
- **Openness**: Release of weights, training/inference code, and data positions Evo 2 among the largest fully open AI models.

## 2.2. Evo 2 predicts mutational effects on protein, RNA, and organismal fitness
- **Zero-shot across modalities**: Likelihood deltas track mutational impact over DNA/RNA/protein and across domains of life, without finetuning.
- **Canonical signals captured**: Start/stop codon sensitivity, triplet periodicity, RBS/Kozak-like upstream patterns; stronger disruption for nonsyn, frameshift, stop; greater sensitivity in constrained elements (tRNA/rRNA).
- **Genetic code awareness**: Distinguishes organism-specific stop codon codes; requires kb-scale context to resolve ciliate code.
- **Experimental concordance**: Correlates with DMS fitness for proteins and ncRNAs; negative association with mRNA decay; no signal on excluded viral proteins (by design).
- **Eukaryotic features**: Embedding-based exon classifier outperforms baselines across species; accurate exon-intron boundary detection on human locus example.
- **Organismal fitness**: Matches Evo 1 on prokaryotic gene essentiality; outperforms baselines on human lncRNA essentiality predictions via perturbation scoring.
- **Takeaway**: Evo 2’s zero-shot likelihoods and embeddings generalize broadly and align with measured biology from molecule to organism.

## 2.3. Human clinical variant effect prediction
- **ClinVar benchmarks**: Competitive on coding SNVs (ranked 4th/5th), state-of-the-art for coding non-SNVs and for noncoding SNVs/non-SNVs.
- **Splicing (SpliceVarDB)**: Best zero-shot performance on both exonic and intronic splice variants.
- **BRCA1/2**: Zero-shot sets SOTA on BRCA1 noncoding SNVs; combined coding+noncoding beats baselines; extends to BRCA2.
- **No human variant training**: Achieves performance using multi-species constraints; only human genome in training is the reference.
- **Supervised boost**: Evo 2 embeddings power a BRCA1 classifier that surpasses zero-shot and specialized baselines (AUROC ~0.95), showing layer-wise embedding utility.
- **Takeaway**: Evo 2 is a robust zero-shot VEP for diverse variant classes and an excellent foundation for supervised clinical predictors.

## 2.4. Feature interpretation from molecular to genome scale
- **SAE setup**: Batch-TopK SAE trained on layer-26 representations over 1B tokens (balanced euk/prok) to extract sparse, interpretable features via contrastive feature search.
- **Mobile element feature**: Feature f/19746 activates on prophages and CRISPR spacers; associates spacers with phage-derived sequences rather than memorizing phage genomes; flags unannotated prophage-like regions.
- **Genome organization features**: Features align with ORFs, intergenic regions, tRNAs, rRNAs; protein-level secondary structure features (alpha-helices, beta-sheets) emerge, indicating multimodal signals.
- **Human regulatory and severity features**: Frameshift/stop-sensitive feature (f/24278) and TF-motif-like activations; exon/intron architecture features generalize to the woolly mammoth genome.
- **Tooling**: Public web viewer released to explore 104-genome SAE features.
- **Takeaway**: Evo 2 latent space encodes rich, biologically grounded concepts enabling annotation, discovery, and potential steering.

## 2.5. Genome-scale generation across the domains of life
- **Prompted gene completion**: High sequence recovery across six diverse species; scales with model size and remains stable through context extension.
- **Safety by exclusion**: Poor generation for human viral proteins consistent with deliberate data exclusions.
- **Organelle generation**: Generates 16 kb mitochondrial genomes with correct gene counts and synteny; AF3 supports plausible multimer folds and interactions.
- **Prokaryotic genome generation**: 600 kb M. genitalium-scale sequences with ~70% Pfam hits (vs 18% for Evo 1); protein properties match natural distributions; AF3 structures plausible.
- **Eukaryotic chromosome generation**: 330 kb yeast chromosome-like sequences with tRNAs, promoters, intron-bearing genes; structural/sequence similarity to natural yeast proteins.
- **Limitations & paths forward**: Feature density lower than native in yeast; unconstrained decoding suggests room for improved inference strategies.
- **Takeaway**: Evo 2 can author long, coherent genomic sequences with realistic coding and regulatory content across domains of life.

## 2.6. Generative epigenomics via inference-time search
- **Goal**: Control the position and width of chromatin-accessible regions in multi-kb sequences.
- **Method**: Guide Evo 2 generation with Enformer+Borzoi scoring inside a beam search; rescore after each 128-bp chunk and continue from top candidates.
- **Inference-time scaling**: Wider beam (more tokens sampled) predictably improves AUROC toward ~0.9; clear relationship between compute and design quality.
- **Demonstrations**: Encode Morse code patterns (“LO”, “ARC”, “EVO2”) as accessibility peaks; designs retain natural dinucleotide frequencies and ensemble consensus.
- **Ablation**: Replacing Evo 2 with uniform proposal degrades scaling, naturalness, and ensemble agreement, hinting at adversarial sequences.
- **General paradigm**: Pair a capable generator with sequence-to-function/structure models to steer complex biological properties beyond accessibility.
- **Takeaway**: Establishes the first inference-time scaling law for biological LMs and a practical recipe for controllable epigenomic design.

## 3. Discussion
- **Unified capability**: A single model achieves generalist prediction (mutational effects, variant pathogenicity) and design (genomes, epigenomic patterns) across domains of life.
- **Interpretability**: SAE features span coding, regulatory, structural, and mobile-element signals, enabling annotation and potential activation steering.
- **Engineering & data**: Achieved via datacenter-scale training, new multi-hybrid architecture, long-context extension, and high-information-content data curation.
- **Safety & bias**: Eukaryotic virus exclusion reduces misuse risk and degrades viral capability as intended; population analyses suggest comparable performance across ancestries.
- **Future directions**: Add population-scale variation, leverage activation steering, supervised/RL with wet-lab feedback, and extend inference-time design to other properties.

## 4. Methods (selected)
- **4.1 Training**: Next-token prediction on byte-tokenized DNA; two-phase schedule (pretraining at 1k→8k context, then midtraining to 1M) with sequence packing. SH2 with RoPE; mixed precision with 3D parallelism; FP8 where applicable. Reweighted CE loss downweights repetitive DNA (0.1). Needle-in-haystack long-context retrieval evaluation defined with categorical Jacobian metric.
- **4.2 Data (OpenGenome2)**: Expanded to 8.84T nt. Rigorous curation for eukaryotes (clustering via Mash distances, filtering ambiguous bases/short contigs), metagenomes (ORF-filtered, redundancy-reduced), organelles, mRNA/ncRNA, promoters. Augmentations: stitched genic windows, exon overhangs, strand flips; special tokens ‘@’/‘#’; phylogenetic tags during midtraining.
- **4.3 Evaluations**: Standardized zero-shot scoring: delta log-likelihood with window normalization; comprehensive benchmarks: start/stop contexts, prok/euk mutation sweeps, genetic code tests, DART-Eval regulatory tasks, DMS protein/ncRNA fitness, mRNA decay, exon/intron classification, gene and lncRNA essentiality, ClinVar, SpliceVarDB, BRCA1/2.
- **4.4 SAEs**: Batch-TopK SAE trained on layer-26 activations; 1B-activation corpus; metrics include activation density and UMAP embeddings; contrastive feature search pipelines for prophage, genome organization, protein secondary structure, frameshift/stop sensitivity, TF motifs, and exon/intron boundaries; public viewer built on igv.js.

## 4.5. Unconstrained generation — methods details
- **Gene completion setup**: Prompts = 1 kb upstream + 0.5–1 kb gene start; decoding with temperature 0.7, top-k 4; evaluate amino-acid recovery across 10 samples per prompt.
- **Mitochondria generation**: 3 kb prompts; 250×16 kb sequences (temp 1.0/0.7, top-k 4); annotate with MitoZ; synteny via LoVis4u; BLASTp diversity; AF3 multimer plausibility.
- **M. genitalium generation**: Prompt 10.5 kb; 35×580 kb sequences (temp 1.0, top-k 4); ORFs via Prodigal; Pfam via HHpred; protein quality via ESMFold pLDDT, DSSP; structures visualized in ChimeraX.
- **Yeast chromosome generation**: 20×330 kb sequences (temp 1.0, top-k 4); genes via GeneMark-S/Prot-hint; promoters via Promoter 2.0; tRNAs via YGAP; fold with ESMFold.

## 4.6.1. Beam search algorithm (formal)
- **Procedure**: Sample K chunks per step from the LM; score partial sequences with an external function f; append best (or keep K′ beams) and iterate until length L.
- **Chunking**: Fixed chunk length C; concatenative growth with rescoring at each step; supports multi-beam K′ ≤ K.

## 4.6.2. Beam search for generative epigenomics
- **Parameters**: C=128, L=19,968, K=42, K′=2; Evo 2 7B, temp 1.0, top-k 4; mouse chrX context upstream/downstream to meet model input lengths.
- **Scoring**: Enformer (896 bins @128 bp) and Borzoi (6144 bins @32 bp) normalized; loss = average L1 to desired binary peak pattern; Borzoi uses ensemble with variance penalty.

## 4.6.3. Design patterns and tasks
- **Peak patterns**: Long/medium/short square waves; Morse code “LO” (dot=768 bp), “ARC” (384 bp), “EVO2” (384 bp); dots/dashes=open chromatin, spaces=closed.

## 4.6.4. Inference-time scaling experiments
- **tok/bp sweep**: Configs from 1 to 60 tok/bp by varying K and K′; AUROC averaged across Enformer and Borzoi predictions vs desired patterns.
- **Result**: More inference-time compute → better adherence to target patterns; uniform proposal baseline underperforms.

## 5. Data availability
- **OpenGenome2**: https://huggingface.co/datasets/arcinstitute/opengenome2

## 6. Code and model availability
- **Code**: Top-level repo, training (Savanna), inference (Vortex), Evo Designer UI, SAE viewer, NVIDIA BioNeMo links.
- **Weights**: Evo 2 40B/7B and base variants, plus 1B base on Hugging Face.

## 7. Acknowledgments
- **Note**: Extensive academic/industry support and funding acknowledgments.

## 8. Author contributions
- **Note**: Project conception, architecture/training/inference, data/SAE curation, evaluations, UI tools, safety analyses, and manuscript writing roles delineated.

## 9. Competing interests
- **Note**: Disclosures provided (various affiliations, advisory, and startup roles); others declare none.

## What this means for disease (plain language)
- **Think of Evo 2 as a “DNA grammar checker + autocomplete”**: It learned what healthy DNA looks like across many species. When it sees a change (a variant), it can tell how “surprising” that change is. Big surprises often mean the change could break function.
- **Faster variant interpretation**: Many patients carry “variants of uncertain significance” (VUS). Evo 2 helps rank which variants are likely harmful, including in the 98% of the genome that doesn’t code for proteins (noncoding) and in splice regions that control how genes are cut and pasted.
- **Better splicing and noncoding insight**: The model’s ability to recognize exon/intron boundaries and transcription factor motifs means it can flag variants that disrupt gene on/off switches or RNA processing—common culprits in rare disease and cancer.
- **From prediction to design**: Beyond reading DNA, Evo 2 can write realistic DNA. Guided by a “coach” (models like Enformer/Borzoi), it can design sequences predicted to open or close chromatin at the right places—useful for gene therapies and cell engineering.
- **Safety by design**: The team intentionally removed human-infecting viral genomes from training. Evo 2 performs poorly on those viruses—which is good for biosecurity.

## How each capability maps to medicine
- **Zero-shot likelihoods (surprise score)**: A quick, alignment-free way to triage variants in clinical reports, including indels and noncoding changes.
- **Exon/intron and motif detection**: Helps annotate poorly understood genes/regions, spot splicing defects, and prioritize follow-up assays (e.g., minigene tests, RNA-seq).
- **BRCA1/2 results**: Shows promise for hereditary cancer genetics; similar pipelines could extend to other clinically important genes.
- **Mechanistic interpretability (SAEs)**: Reveals what the model “notices” (e.g., motifs, exon boundaries). That transparency builds trust and enables feature-level tools (e.g., scanning for risky promoter motifs).
- **Genome-scale generation**: Aids rapid prototyping of constructs, libraries, and test sequences that still “look natural” to cellular machinery.
- **Controllable epigenomic design**: Early path to tuning expression without changing the underlying gene—potentially safer, more targeted control knobs for future therapies.

## Figure guide (what to look for, in plain terms)
- **Figure 1 (overview)**: How the model is built and trained to read long DNA (up to a million letters) and still remember important bits far apart—important for diseases where distant DNA affects genes.
- **Figure 2 (mutation effects across life)**: The model rediscovers core genetics (start/stop patterns, codon periodicity) and shows predictions align with lab measurements. This builds confidence that disease variant scores are meaningful.
- **Figure 3 (human variants)**: Benchmark results showing strong performance on noncoding and splice variants and competitive coding predictions. BRCA tests demonstrate clinical relevance.
- **Figure 4 (interpretability)**: Shows the model’s “concepts” lines up with biology (motifs, exon-intron, protein structure) and even helps annotate extinct genomes—evidence of generalization.
- **Figure 5 (generation)**: The model can write entire mitochondria, minimal bacteria, and yeast-like chromosomes that look realistic to analysis tools—useful for design.
- **Figure 6 (epigenomic design)**: Using a steering strategy (beam search + scoring), more compute yields better designs that match desired accessibility patterns—an actionable knob for future therapeutic design.

## Practical disease workflows this can enable
- **Clinical genetics**: Prioritize VUS in coding and noncoding regions; add splicing-aware scoring to reports; route top candidates to functional assays.
- **Rare disease research**: Rapidly annotate genes and regulatory elements in understudied organisms or patient-derived assemblies; generate hypotheses for regulatory variants.
- **Oncology**: Improve triage of tumor variants including noncoding drivers; design regulatory elements for cell models and screens.
- **Therapeutic design**: Prototype regulatory DNA to achieve cell-type-specific expression; explore safer control via epigenomic patterning instead of editing the gene body.

## Limitations and care points
- **Lab validation still needed**: Model scores are hypothesis generators, not diagnoses. Use as decision support alongside assays (e.g., RNA-seq, MPRA, CRISPR tests).
- **Eukaryotic viruses excluded**: By design, viral predictions/generation are weak—good for safety, but not a tool for viral genomics.
- **Bias and coverage**: Early population analyses look comparable to other methods, but equitable performance needs continued evaluation with broader data.

## If you’re new to this area, here’s the gist
- **Evo 2 learned the “rules” of DNA** from many species, so it can spot changes that likely break those rules (possible disease-causing variants).
- **It can also write DNA** that follows those rules and be steered to create sequences with desired activity patterns, a step toward programmable biology.
- **It’s open** (data, code, models), so clinicians, researchers, and tool builders can test, validate, and improve it for real medical impact. 