# Methods

## Overview
We present a stage‑aware CRISPR guide design framework (Interception) that maps genomic signals to the metastatic cascade and ranks intervention candidates by on‑target efficacy, genome‑wide safety, and mission fit. The system is implemented as a modular FastAPI backend with Pydantic schemas and a React frontend. All responses include full provenance (run id, methods, model ids, feature flags) and are packaged with scripts for complete reproducibility.

## Data Sources and Inputs
- Variants: GRCh38 (chrom, pos, ref, alt), optionally gene and HGVS protein.
- Model endpoints: Evo2 proxy for sequence‑level signals; Enformer/Borzoi (local stubs) for chromatin accessibility.
- Gene sets and weights: Config files in `api/config/` with versioning.
- Reference genome: GRCh38 for off‑target alignment.

## Backend Architecture
The backend exposes routers for insights, design, fusion (optional), and metastasis interception:
- Insights (`api/routers/insights.py`): functionality change, essentiality, regulatory/splicing proxy, chromatin accessibility.
- Design (`api/routers/design.py`): `POST /api/design/predict_crispr_spacer_efficacy` (Evo2 delta→sigmoid).
- Interception (`api/routers/metastasis_intercept.py`): orchestrates target lock → candidate design → safety preview → assassin ranking.

Services are modularized under `api/services/` and configuration is versioned; provenance includes ruleset version and profile.

## Target Lock (Mission Fit)
For a selected metastatic step, a ruleset maps the mission to one or more curated gene sets (e.g., MAPK for primary growth, EMT_TF for local invasion, MMP for invasion). For each candidate gene, we query four biological signals and compute a weighted score:

Target‑Lock formula (default weights):
```
target_lock = 0.35 * functionality
             + 0.35 * essentiality
             + 0.15 * chromatin
             + 0.15 * regulatory
```
Each signal is normalized to [0,1] and threshold‑gated (default ≥0.6); contributions are summed and clipped to [0,1]. The top gene by `target_lock` is the validated target for design.

## Functionality (Evo2 multi + exon context)
We estimate variant‑level functionality change using Evo2 sequence likelihood shifts with two complementary calls:
1) `score_variant_multi` (broad, multi‑window context) – returns `min_delta`.
2) `score_variant_exon` (deep exon context with large flanks; default flank=8192 bp) – returns `exon_delta`.

We combine absolute magnitudes with a bias toward exon context for coding variants. Hotspot domain hints (e.g., BRAF V600, RAS G12/G13/Q61) add a small lift. Contract accepts either top‑level fields or `mutations[0]` for flexibility.

## Essentiality (Evo2 gene‑level proxy)
We compute essentiality as an Evo2 magnitude aggregation across relevant exons/windows for the gene. Scores are normalized to [0,1] with conservative calibration; used as a relative prioritization signal rather than binary lethality.

## Regulatory / Splicing Proxy
We use Evo2 sequence delta in non‑coding/regulatory contexts as a proxy for splice/regulatory impact and normalize to [0,1]. This signal is conservatively weighted in target lock.

## Chromatin Accessibility (Enformer/Borzoi)
When available, we query a local Enformer stub (`ENFORMER_URL`) to estimate accessibility in a window centered at the variant (default ±1–5 kb depending on request). If Enformer is unavailable, we attempt a Borzoi stub; otherwise we fall back to a simple heuristic. Responses include `provenance.method ∈ {enformer, borzoi, heuristic}` and a confidence hint (0.55–0.60 for models, 0.40 for heuristic).

## Design (Guide Candidate Generation)
For the validated target locus, PAM‑aware candidates (NGG) are enumerated in 20 bp spacers within the target window. Each candidate is annotated with GC%, homopolymer flags, and a preliminary safety heuristic if alignment services are unavailable. The number of candidates per target is configurable (default n=3 for demos; n=10 for validation datasets).

## Efficacy (Evo2 delta → sigmoid)
On‑target efficacy is estimated by scoring the local context with Evo2 and transforming the delta via a sigmoid:
```
efficacy = 1 / (1 + exp(delta / scale))
```
Default `scale=10.0`. Context can be provided directly, fetched via Ensembl using (chrom, pos, ref, alt) with configurable ±window size (default 150 bp flanks; total 300 bp), or fall back to guide‑only with reduced confidence.

## Safety (Genome‑wide Off‑Target Search)
We perform hierarchical alignment against GRCh38:
1) Primary: minimap2 short‑read preset (`-x sr`, `-N 100`).
2) Fallback: BLAST (blastn-short) with tabular output.

Hits are parsed to extract mismatch counts (0–3 mismatches) and summarized as total off‑targets and per‑mismatch distribution. Safety is mapped via exponential decay:
```
safety = exp(-0.5 * total_hits)
```
This penalizes promiscuous guides aggressively; 0 hits → 1.0, 1 hit → ~0.606, 5 hits → ~0.082.

## Assassin Score (Composite Ranking)
Candidates are ranked by a composite score balancing efficacy, safety, and mission fit:
```
assassin = 0.40 * efficacy + 0.30 * safety + 0.30 * target_lock
```
Weights reflect a pragmatic priority: guides must cut efficiently, avoid off‑targets, and align with the chosen mission biology. Scores are clipped to [0,1] and surfaced with per‑candidate provenance (methods, model ids, inputs, context length).

## Publication Figures and Tables
We provide scripts to generate all figures and tables from live endpoints:
- `scripts/generate_target_lock_data_v2.py`: heatmap (synthetic fast path).
- `scripts/generate_publication_data_real.py`: real ClinVar‑like variants (heatmap + guides).
- `scripts/generate_guide_validation_data.py`: F3–F5 and Table 2 from ranked candidates.

Generated artifacts are copied to `publication/figures/`, `publication/tables/`, and `publication/data/` and referenced in `PUBLICATION_OUTPUT_SUMMARY.mdc`.

## Reproducibility and Environment
We provide `publication/environment.yml` and a stubbed chromatin setup (`tools/chromatin/enformer_server.py`, `.../borzoi_server.py`). Required environment variables: `ENFORMER_URL`, `BORZOI_URL`. The `SUBMISSION_PACKAGE.md` includes exact commands, file manifest, and SHA‑256 checksums.

## Statistics and Plotting
We compute summary statistics (mean ± SD, min, max, median) for efficacy, safety, and assassin scores across mission steps and render violin/box plots using matplotlib/seaborn. Decimal precision in figure labels matches Table 2.

## Limitations and RUO Compliance
Chromatin models are stubbed for deterministic development; production deployments should point to validated services. Functionality and essentiality are conservative, pending wet‑lab calibration. All results are Research Use Only; we do not make clinical claims. Provenance is provided for auditability and reproducibility.

## Availability
All code paths, configurations, and figure/data generation scripts are included in the repository under version control. A submission bundle (`publication_submission_v1.zip`) contains the full publication package.
