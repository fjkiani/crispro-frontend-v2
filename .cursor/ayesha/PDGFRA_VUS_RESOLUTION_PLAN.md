# PDGFRA p.S755P (c.2263T>C) VUS Resolution Plan

## Objective
Eliminate the VUS classification for PDGFRA p.S755P by running it through the complete 8-step VUS identification workflow.

## Variant Information
- **Gene**: PDGFRA
- **Protein Change**: p.S755P (Serine to Proline at position 755)
- **cDNA Change**: c.2263T>C
- **Classification**: Currently VUS (Variant of Unknown Significance)
- **Zygosity**: Heterozygous

## 8-Step VUS Identification Workflow

### Step 1: Inputs and Normalization
- **Input**: HGVS protein notation `p.S755P` and cDNA `c.2263T>C`
- **Action**: Resolve to GRCh38 coordinates (chromosome, position, ref, alt)
- **Output**: Normalized variant coordinates

### Step 2: Priors (Decide "Already Not VUS?")
- **ClinVar Lookup**: Check for existing classification
  - Endpoint: `/api/evidence/clinvar` or `/api/evidence/deep_analysis`
- **AlphaMissense Coverage**: Check Fusion Engine eligibility
  - Endpoint: `/api/fusion/coverage?chrom=...&pos=...&ref=...&alt=...`
- **Output**: Prior classification status (if available)

### Step 3: Triumvirate Protocol Gate
- **Check**: Is this a truncating/frameshift variant?
- **Logic**: p.S755P is a missense variant (Ser→Pro), so NOT truncating
- **Output**: PASS (continue to scoring)

### Step 4: S/P/E Core Signals
- **Sequence (S)**: Evo2 delta scores
  - Endpoint: `/api/evo/score_variant_multi` or via insights
- **Pathway (P)**: Variant-to-pathway mapping
  - PDGFRA is in RTK/MAPK pathway
- **Evidence (E)**: Literature strength, ClinVar priors
  - Endpoint: `/api/evidence/deep_analysis`
- **Output**: S/P/E scores with rationale

### Step 5: Insights Bundle (Real Signals)
- **Functionality**: Protein function change prediction
  - Endpoint: `/api/insights/predict_protein_functionality_change`
- **Chromatin**: Accessibility and regulatory impact
  - Endpoint: `/api/insights/predict_chromatin_accessibility`
- **Essentiality**: Gene dependency scoring
  - Endpoint: `/api/insights/predict_gene_essentiality`
- **Regulatory**: Splicing and non-coding impact
  - Endpoint: `/api/insights/predict_splicing_regulatory`
- **Output**: 4 insight scores with thresholds

### Step 6: Fusion Gating (AlphaMissense)
- **Coverage Check**: GRCh38 missense variants only
  - p.S755P is missense → eligible
- **Fusion Integration**: AlphaMissense scores when available
  - Endpoint: `/api/fusion/score_variant`
- **Output**: Fused S scores with provenance

### Step 7: SAE Interpretable Features (Real Data Only)
- **Feature Extraction**: Sparse autoencoder features from Evo2
- **Interpretation**: Feature attribution and activation steering
- **Output**: SAE features with interpretation

### Step 8: VUS Triage Scoring and Decision Rubric
- **Confidence Calculation**: S/P/E + insights + Fusion + SAE
- **Tier Classification**: Supported/Consider/Insufficient
- **Decision Logic**: Clear pass/fail/conditional verdicts
- **Output**: Final VUS triage verdict with confidence

## Expected Resolution Paths

### Path 1: Benign/Likely Benign
- **If**: Low functionality score (<0.3), ClinVar benign, low essentiality
- **Action**: Reclassify as "Benign" or "Likely Benign"
- **Outcome**: VUS eliminated

### Path 2: Pathogenic/Likely Pathogenic
- **If**: High functionality score (>0.7), ClinVar pathogenic, high essentiality
- **Action**: Reclassify as "Pathogenic" or "Likely Pathogenic"
- **Outcome**: VUS eliminated

### Path 3: Remains VUS
- **If**: Mixed signals, insufficient evidence
- **Action**: Document evidence gaps, recommend additional testing
- **Outcome**: VUS remains but with clear evidence summary

## Implementation Script
See `resolve_pdgra_vus.py` for automated execution.

## Success Criteria
1. ✅ Complete 8-step workflow executed
2. ✅ All insights endpoints called successfully
3. ✅ ClinVar and Fusion coverage checked
4. ✅ Final classification determined (Benign/Pathogenic/VUS)
5. ✅ Evidence summary generated
6. ✅ VUS eliminated OR clear path forward documented




## Objective
Eliminate the VUS classification for PDGFRA p.S755P by running it through the complete 8-step VUS identification workflow.

## Variant Information
- **Gene**: PDGFRA
- **Protein Change**: p.S755P (Serine to Proline at position 755)
- **cDNA Change**: c.2263T>C
- **Classification**: Currently VUS (Variant of Unknown Significance)
- **Zygosity**: Heterozygous

## 8-Step VUS Identification Workflow

### Step 1: Inputs and Normalization
- **Input**: HGVS protein notation `p.S755P` and cDNA `c.2263T>C`
- **Action**: Resolve to GRCh38 coordinates (chromosome, position, ref, alt)
- **Output**: Normalized variant coordinates

### Step 2: Priors (Decide "Already Not VUS?")
- **ClinVar Lookup**: Check for existing classification
  - Endpoint: `/api/evidence/clinvar` or `/api/evidence/deep_analysis`
- **AlphaMissense Coverage**: Check Fusion Engine eligibility
  - Endpoint: `/api/fusion/coverage?chrom=...&pos=...&ref=...&alt=...`
- **Output**: Prior classification status (if available)

### Step 3: Triumvirate Protocol Gate
- **Check**: Is this a truncating/frameshift variant?
- **Logic**: p.S755P is a missense variant (Ser→Pro), so NOT truncating
- **Output**: PASS (continue to scoring)

### Step 4: S/P/E Core Signals
- **Sequence (S)**: Evo2 delta scores
  - Endpoint: `/api/evo/score_variant_multi` or via insights
- **Pathway (P)**: Variant-to-pathway mapping
  - PDGFRA is in RTK/MAPK pathway
- **Evidence (E)**: Literature strength, ClinVar priors
  - Endpoint: `/api/evidence/deep_analysis`
- **Output**: S/P/E scores with rationale

### Step 5: Insights Bundle (Real Signals)
- **Functionality**: Protein function change prediction
  - Endpoint: `/api/insights/predict_protein_functionality_change`
- **Chromatin**: Accessibility and regulatory impact
  - Endpoint: `/api/insights/predict_chromatin_accessibility`
- **Essentiality**: Gene dependency scoring
  - Endpoint: `/api/insights/predict_gene_essentiality`
- **Regulatory**: Splicing and non-coding impact
  - Endpoint: `/api/insights/predict_splicing_regulatory`
- **Output**: 4 insight scores with thresholds

### Step 6: Fusion Gating (AlphaMissense)
- **Coverage Check**: GRCh38 missense variants only
  - p.S755P is missense → eligible
- **Fusion Integration**: AlphaMissense scores when available
  - Endpoint: `/api/fusion/score_variant`
- **Output**: Fused S scores with provenance

### Step 7: SAE Interpretable Features (Real Data Only)
- **Feature Extraction**: Sparse autoencoder features from Evo2
- **Interpretation**: Feature attribution and activation steering
- **Output**: SAE features with interpretation

### Step 8: VUS Triage Scoring and Decision Rubric
- **Confidence Calculation**: S/P/E + insights + Fusion + SAE
- **Tier Classification**: Supported/Consider/Insufficient
- **Decision Logic**: Clear pass/fail/conditional verdicts
- **Output**: Final VUS triage verdict with confidence

## Expected Resolution Paths

### Path 1: Benign/Likely Benign
- **If**: Low functionality score (<0.3), ClinVar benign, low essentiality
- **Action**: Reclassify as "Benign" or "Likely Benign"
- **Outcome**: VUS eliminated

### Path 2: Pathogenic/Likely Pathogenic
- **If**: High functionality score (>0.7), ClinVar pathogenic, high essentiality
- **Action**: Reclassify as "Pathogenic" or "Likely Pathogenic"
- **Outcome**: VUS eliminated

### Path 3: Remains VUS
- **If**: Mixed signals, insufficient evidence
- **Action**: Document evidence gaps, recommend additional testing
- **Outcome**: VUS remains but with clear evidence summary

## Implementation Script
See `resolve_pdgra_vus.py` for automated execution.

## Success Criteria
1. ✅ Complete 8-step workflow executed
2. ✅ All insights endpoints called successfully
3. ✅ ClinVar and Fusion coverage checked
4. ✅ Final classification determined (Benign/Pathogenic/VUS)
5. ✅ Evidence summary generated
6. ✅ VUS eliminated OR clear path forward documented









