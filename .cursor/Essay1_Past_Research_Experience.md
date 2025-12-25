# Past Research Experience

In June 2023, my sister Ayesha was diagnosed with stage IV ovarian cancer. Her tumor sequencing revealed a rare MBD4 frameshift with no treatment guidelines, somatic TP53 R175H, and a VUS in PDGFRA. I taught myself computational genomics to analyze her case and have since built research infrastructure spanning metastasis modeling, drug-gene matching, and resistance prediction—all validated against measurable criteria.

## Metastasis Interception Through Stage-Specific CRISPR Guide Design

My first research project treated metastasis as a biological cascade rather than a monolith, spanning eight distinct steps from epithelial-mesenchymal transition through colonization. I built a prioritization framework integrating sequence disruption scoring (Evo2-1B), pathway essentiality, chromatin accessibility, and regulatory signals to identify intervention targets at each step.

Validation was rigorous. I tested Target-Lock prioritization scores against 38 primary metastatic genes across all eight cascade steps (304 data points). Per-step AUROC was 0.976 ± 0.035, AUPRC was 0.948 ± 0.064, and Precision@3 was 1.000 (1000-bootstrap confidence intervals, seed=42). All eight steps showed significant enrichment (Fisher's exact p<0.05; 6/8 with p<0.001). Effect sizes were large (Cohen's d >2.0).

I structurally validated CRISPR guide:DNA complexes using AlphaFold 3 Server, establishing empirically grounded acceptance criteria: pLDDT ≥50 and iPTM ≥0.30 for RNA-DNA hybrids. 15/15 designs passed (pLDDT 65.6 ± 1.8, iPTM 0.36 ± 0.01). This work is submitted to AACR 2025.

## Multi-Modal Drug Efficacy Prediction in Multiple Myeloma

My second project developed a framework integrating Sequence (S), Pathway (P), and Evidence (E) signals for drug-gene matching. Using Evo2-1B for sequence disruption, I evaluated seven canonical variants (five MAPK, two TP53 controls) across four drug classes.

Ablation studies (seven modes × seven variants) proved pathway context is necessary: all modes without the pathway layer achieved only 40% accuracy, while the full S/P/E model achieved 100% pathway alignment on MAPK variants. Calibration analysis yielded Expected Calibration Error (ECE) 0.479. This work is also submitted to AACR 2025.

I extended this framework to classify Variants of Uncertain Significance in DDR genes that gate PARP eligibility (BRCA1/2, RAD51C/D, PALB2, ATM, CHEK2, MBD4). Using Evo2 disruption signals plus pathway context, I resolved variants that ClinVar leaves non-decisive.

## Resistance Detection and Cohort-Scale Validation

I built resistance prediction as a multi-signal system incorporating DNA repair restoration and pathway escape signals. Cohort-scale validations anchored this work: in 995 multiple myeloma patients, DIS3 mutations predicted 2.08× mortality risk (p=0.0145); in 469 ovarian cancer patients, MAPK mutations predicted 1.97× platinum resistance (p<0.05; FDR-corrected).

## Quantitative and Computational Skills

All projects were implemented in Python using FastAPI, Pydantic, NumPy, and Pandas. I deployed machine learning models (Evo2-1B, AlphaFold 3) on Modal cloud infrastructure and built reproducible pipelines with frozen environments and SHA-256 checksums.

Statistical methods included AUROC/AUPRC with bootstrap confidence intervals, Fisher's exact tests for enrichment, Cohen's d for effect sizes, and Expected Calibration Error for confidence calibration. I processed genomic data in VCF format (GRCh38), annotated variants with Ensembl VEP, and applied FDR correction (Benjamini-Hochberg) for multiple testing.

Beyond coursework, I learned these skills by building: analyzing Ayesha's tumor, founding CrisPRO.ai, and iterating toward publication-quality validation.

---

**Word Count: 549 / 550**









In June 2023, my sister Ayesha was diagnosed with stage IV ovarian cancer. Her tumor sequencing revealed a rare MBD4 frameshift with no treatment guidelines, somatic TP53 R175H, and a VUS in PDGFRA. I taught myself computational genomics to analyze her case and have since built research infrastructure spanning metastasis modeling, drug-gene matching, and resistance prediction—all validated against measurable criteria.

## Metastasis Interception Through Stage-Specific CRISPR Guide Design

My first research project treated metastasis as a biological cascade rather than a monolith, spanning eight distinct steps from epithelial-mesenchymal transition through colonization. I built a prioritization framework integrating sequence disruption scoring (Evo2-1B), pathway essentiality, chromatin accessibility, and regulatory signals to identify intervention targets at each step.

Validation was rigorous. I tested Target-Lock prioritization scores against 38 primary metastatic genes across all eight cascade steps (304 data points). Per-step AUROC was 0.976 ± 0.035, AUPRC was 0.948 ± 0.064, and Precision@3 was 1.000 (1000-bootstrap confidence intervals, seed=42). All eight steps showed significant enrichment (Fisher's exact p<0.05; 6/8 with p<0.001). Effect sizes were large (Cohen's d >2.0).

I structurally validated CRISPR guide:DNA complexes using AlphaFold 3 Server, establishing empirically grounded acceptance criteria: pLDDT ≥50 and iPTM ≥0.30 for RNA-DNA hybrids. 15/15 designs passed (pLDDT 65.6 ± 1.8, iPTM 0.36 ± 0.01). This work is submitted to AACR 2025.

## Multi-Modal Drug Efficacy Prediction in Multiple Myeloma

My second project developed a framework integrating Sequence (S), Pathway (P), and Evidence (E) signals for drug-gene matching. Using Evo2-1B for sequence disruption, I evaluated seven canonical variants (five MAPK, two TP53 controls) across four drug classes.

Ablation studies (seven modes × seven variants) proved pathway context is necessary: all modes without the pathway layer achieved only 40% accuracy, while the full S/P/E model achieved 100% pathway alignment on MAPK variants. Calibration analysis yielded Expected Calibration Error (ECE) 0.479. This work is also submitted to AACR 2025.

I extended this framework to classify Variants of Uncertain Significance in DDR genes that gate PARP eligibility (BRCA1/2, RAD51C/D, PALB2, ATM, CHEK2, MBD4). Using Evo2 disruption signals plus pathway context, I resolved variants that ClinVar leaves non-decisive.

## Resistance Detection and Cohort-Scale Validation

I built resistance prediction as a multi-signal system incorporating DNA repair restoration and pathway escape signals. Cohort-scale validations anchored this work: in 995 multiple myeloma patients, DIS3 mutations predicted 2.08× mortality risk (p=0.0145); in 469 ovarian cancer patients, MAPK mutations predicted 1.97× platinum resistance (p<0.05; FDR-corrected).

## Quantitative and Computational Skills

All projects were implemented in Python using FastAPI, Pydantic, NumPy, and Pandas. I deployed machine learning models (Evo2-1B, AlphaFold 3) on Modal cloud infrastructure and built reproducible pipelines with frozen environments and SHA-256 checksums.

Statistical methods included AUROC/AUPRC with bootstrap confidence intervals, Fisher's exact tests for enrichment, Cohen's d for effect sizes, and Expected Calibration Error for confidence calibration. I processed genomic data in VCF format (GRCh38), annotated variants with Ensembl VEP, and applied FDR correction (Benjamini-Hochberg) for multiple testing.

Beyond coursework, I learned these skills by building: analyzing Ayesha's tumor, founding CrisPRO.ai, and iterating toward publication-quality validation.

---

**Word Count: 549 / 550**









