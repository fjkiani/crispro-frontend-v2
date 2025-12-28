# Additional validation cohorts (shortlist)

## What we need
- Ovarian cancer cohort with **expression** and **platinum response/sensitivity** labels (binary or ordinal).
- Prefer n≥80 with ≥20 resistant.

## Confirmed usable (already integrated)
- **GSE63885 (GPL570)**: platinum sensitivity labels; n=101 (34 resistant / 67 sensitive). Primary external validation.

## High-priority candidates to screen next

### 1) GSE32062 (n≈270)
- Metadata indicates a `platinum` field exists in sample characteristics.
- Next step: determine whether it encodes sensitivity/response (not just “received platinum”).

### 2) curatedOvarianData (Bioconductor)
- A curated collection of clinically annotated ovarian expression datasets.
- Next step: identify which included datasets have **platinum response** annotations and extract a harmonized response label.

### 3) Additional GEO cohorts with explicit response fields
We observed GEO cohorts where sample characteristics include response fields such as:
- “clinical status post 1st line chemotherapy” (CR/PR/SD/P)
- “platinum sensitivity” with time-to-relapse based definitions

Next step: programmatically scan a short list of ovarian GEO series for these fields and sample counts.

## Acceptance criteria (for CCR)
- At least **one additional independent validation cohort** reproducing MFAP4 directionality with AUROC ≥0.65.
- Ideally demonstrate robustness across platform (microarray vs RNA-seq) and across label definitions.
