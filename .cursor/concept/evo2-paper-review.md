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