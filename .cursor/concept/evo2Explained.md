Evo2 paper — pass 1 (lines 1–300): key takeaways for oncology, R&D de‑risking, and VUS
Capabilities
7B/40B model trained on 9.3T tokens with 1M-token context; single-nucleotide resolution.
Zero-shot prediction of variant impact across DNA/RNA/protein; strong on noncoding; competitive on coding.
Mechanistic interpretability via SAEs: captures exon/intron, TF motifs, protein structural elements, prophage regions.
Genome-scale generation; controllable design via inference-time search (e.g., shaping chromatin accessibility).
Architecture: StripedHyena 2 multi-hybrid; efficient long-context recall to 1M tokens.
Open release: models, training/inference code, and OpenGenome2 dataset.
Security: excludes eukaryotic viruses; poor performance on viral proteins by design.
Oncology application angles
Zero-shot pathogenicity for both coding and noncoding human variants → quick target plausibility checks.
Splice-altering variant prediction (later sections) and noncoding regulatory sensitivity → regulatory VUS relevance.
Embeddings power small supervised models; used to classify BRCA1 VUS at SOTA in their results.
R&D de-risking levers
“Delta likelihood” as a universal, cheap first filter for variant impact (coding and noncoding).
Long-context sensitivity (up to 1M bp) enables context-aware assessment near start/stop codons, UTRs, TF sites.
Mechanistic interpretability features support auditability and rationale for decisions in regulatory settings.
Sequence generation + inference-time control can prototype CRISPR targets or regulatory edits with desired accessibility patterns.
VUS resolution pathway (implied by paper)
Step 1: Zero-shot pathogenicity/splicing likelihoods to triage.
Step 2: Use embeddings to train task-specific classifiers (e.g., BRCA1) on curated assay labels.
Step 3: Combine with orthogonal evidence (gene essentiality, chromatin context) for higher-confidence calls.


Evo2 VEP performance
Coding SNVs: ranks 4th/5th (behind AlphaMissense, ESM-1b, GPN-MSA).
Coding non-SNVs (indels): outperforms others zero-shot.
Noncoding SNVs and non-SNVs: state-of-the-art zero-shot.
Splicing (SpliceVarDB): best zero-shot for exonic and intronic splice variants.
BRCA1/BRCA2: strong zero-shot; best when coding+noncoding combined. Using embeddings, a supervised BRCA1 classifier hits AUROC 0.94–0.95.
Implications for oncology/VUS
Zero-shot triage across all variant classes (coding, noncoding, indels, splice) reduces wet-lab validation burden.
Embeddings enable high-accuracy supervised VUS classifiers for specific genes (e.g., BRCA1).
Long-context modeling captures start/stop codon, UTR, and regulatory context for more reliable calls.
Particularly valuable for noncoding and splice VUS, where standard tools are weaker.
Interpretability
SAE features align with exon/intron boundaries, TF motifs, protein structural signatures, prophage/mobile elements.
Mutation-severity features (frameshift/stop) identified; supports explainability and regulatory traceability.
Generative capabilities
Genome-scale generation: realistic mitochondria, minimal bacterial genomes, yeast chromosome-scale sequences; structural plausibility validated with AlphaFold 3.
Generative epigenomics: inference-time search guided by Enformer/Borzoi designs specified chromatin-accessibility patterns; more compute → better adherence (inference-time scaling).
Security boundary: fails on human virus proteins (by design).
R&D de-risking use
Workflow: zero-shot triage → supervised refinement on embeddings for high-stakes genes → interpretable evidence → design regulatory/guide sequences with targeted accessibility profiles.
Addresses noncoding/splice VUS at scale and informs CRISPR targetability via designed accessibility.


---

## Evo2 paper — pass 2 (lines 301–600): clinical VEP, splicing, BRCA1/2, generation primer
- Coding SNVs (ClinVar): Evo2 40B/7B rank ~4th/5th; competitive but not top vs AlphaMissense/ESM/GPN-MSA.
- Non-SNV (indels) and Noncoding SNVs: Evo2 sets the pace zero-shot; this fills a major tooling gap.
- SpliceVarDB: Best zero-shot for exonic and intronic splice-altering variants → strong splicing VUS coverage.
- BRCA1/BRCA2 (functional data): zero-shot strong; when combining coding+noncoding variants Evo2 leads. Using embeddings, a supervised BRCA1 head achieves AUROC ~0.94–0.95.
- Oncology/R&D impact:
  - Triages VUS across coding/noncoding/indel/splice, lowering wet-lab burden.
  - For high-stakes genes (e.g., BRCA1), train small supervised heads on Evo2 embeddings to resolve borderline VUS.

## Evo2 paper — pass 3 (lines 601–1000): discussion, safety, open resources, future work
- Open release: models, training + inference code, OpenGenome2 dataset; UI tools for design and SAE feature exploration.
- Safety & ethics: Eukaryotic viral genomes excluded; Evo2 intentionally weak on human viral proteins (generation + fitness) → minimizes dual-use risk.
- Population bias: “Population-free” design generalizes comparably across ancestries in tested benchmarks; still monitor bias.
- Future directions: integrate population-scale variation; steer generations via SAE features; RL/finetuning with experimental feedback; inference-time compute for complex properties (e.g., cell-type specificity).
- Oncology/R&D impact:
  - Safe default: block viral requests in generation endpoints.
  - Plan for ancestry-aware audits when using outputs clinically.

## Evo2 paper — pass 4 (lines 1001–1400): data curation, tokenization, windowing, zero-shot protocol
- Curated training data across organelles, transcripts, ncRNA, promoters; clustering to reduce redundancy.
- Augmentations: stitched gene windows, promoter/exon/splice overhangs; special tokens ‘@’, ‘#’; reverse-complementing.
- Midtraining: stitched long sequences + phylogenetic tags → long-context learning to 1M tokens.
- Zero-shot scoring: use 8,192-nt windows centered on variant; for indels keep total window length constant by trimming/padding flanks.
- Some tasks require ≥4 kb context (e.g., genetic code recognition in ciliates).
- Oncology/R&D impact:
  - Our endpoints should fetch 8,192-nt context; support 4–8 kb adjustable windows.
  - Always normalize indel windows to equal length to keep scores comparable.

## Evo2 paper — pass 5 (lines 1401–1800): mechanistic interpretability (SAE)
- Batch-TopK SAE trained on Evo2 activations reveals features for: prophages/mobile elements, CRISPR spacers, ORFs, tRNA/rRNA, protein secondary structures, TF motifs, exon–intron boundaries, and mutation-severity (frameshift/stop) signatures.
- Public viewer for 104 genomes exposes feature activations.
- Oncology/R&D impact:
  - Add explainability artifacts: report which SAE feature families activate near a variant (e.g., TF motif proximity, exon boundary).
  - Use these features as orthogonal evidence in VUS adjudication and dossier narratives.

## Evo2 paper — pass 6 (lines 1801–2200): metadata, references
- Acknowledgments and references; confirms integration points (Enformer/Borzoi, AlphaFold 3, ProteinGym, DART‑Eval).

## Evo2 paper — pass 7 (lines 2201–2600): supplements, ablations, data composition
- Repeat down-weighting (0.1×) improves ClinVar AUROC in ablation; adopted in final training.
- Data composition: enriching genic windows/augmented transcripts improves zero-shot across ClinVar/SpliceVarDB/BRCA1 (esp. combined coding+noncoding); whole-genome-heavy ablation underperforms.
- Oncology/R&D impact:
  - Reinforces our focus on function-rich context windows when scoring; discourage scoring on low-information, repeat-heavy windows.

## Evo2 paper — pass 8 (lines 2601–4273): detailed tables/metrics
- Large results tables corroborate headline claims:
  - Noncoding + non-SNV SOTA zero-shot; splice SOTA; BRCA1 supervised head on embeddings outperforms baselines.
  - Confirms window sizes, variant class coverage, and evaluation splits.

---

## Consolidated implementation guidance (for later coding)
- Scoring endpoints
  - Delta likelihood: `/evo2/score_delta`, `/evo2/score_batch` with 8,192-nt windows; constant-length handling for indels.
  - Splice-specific: `/evo2/score_splice` (exonic/intronic flag) using same windowing.
  - Regulatory: `/evo2/score_regulatory` for cCRE/TF tasks; likelihood or embedding-based scores.
- VUS resolution
  - Embedding head for BRCA1/2: extract block-20 (40B) or best-performing block; average 128-nt windows around variant; small MLP head returns `prob_pathogenic`.
  - Include SAE-derived feature hits and TF motif proximity as explainability.
- Design (guarded)
  - Prompted completion (non-viral); inference-time guided design for chromatin accessibility via Enformer/Borzoi; expose compute budget controls; return AUROC of adherence.
- Safety & fairness
  - Reject viral protein generation; log ancestry-sensitive contexts; expose disclaimers that Evo2 is population-free and perform bias audits before clinical use.