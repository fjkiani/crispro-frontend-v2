---
title: "Tumor Board Packet — Bundle Only (AK, L1, 2026-02-12)"
source_bundle: "artifacts/payload.md (embedded JSON bundle)"
derived_clean_json: "artifacts/payload.extracted.json (machine-extracted copy)"
ruo: true
contract:
  - "This document uses ONLY what the captured bundle contains."
  - "Every factual claim must cite a JSON field path present in the bundle."
  - "If the bundle is silent, we say Unknown and name the missing field."
  - "No re-runs. No backend code. No equations. No inference beyond what fields encode."
---

## 0) What this bundle is (and why it exists)

This packet is a **single captured API response** intended to be a quote-safe “tumor board packet.”

- **Endpoint captured**: `POST /api/ayesha/therapy-fit/bundle?level=l1&include_synthetic_lethality=true` (bundle header statement: see `## 0` in the legacy scribe; the bundle object itself is stored at `artifacts/payload.extracted.json`)
- **What it contains (at L1)**:
  - inputs used (`levels.L1.inputs_used`)
  - WIWFM efficacy ranking (`levels.L1.efficacy`)
  - Synthetic Lethality panel (`levels.L1.synthetic_lethality`)
  - completeness + missing fields (`levels.L1.completeness`, `tests_needed`)

## 1) Bundle contents — minimal factual synopsis (no interpretation)

### 1.1 Identity / timing

- **Contract**: `contract_version = "v2.0"` (`contract_version`)
- **Patient ID**: `patient_id = "AK"` (`patient_id`)
- **Generated at**: `generated_at = "2026-02-12T16:51:53.836997"` (`generated_at`)
- **Levels**: `requested_levels = ["L1"]` (`requested_levels`) and `levels` contains `L1` (`levels.L1`)

### 1.2 Inputs actually used at L1

Mutations (`levels.L1.inputs_used.mutations`):
- **MBD4**: HGVS only; no chrom/pos/ref/alt present (`levels.L1.inputs_used.mutations[0]`)
- **PDGFRA**: HGVS only; no chrom/pos/ref/alt present (`levels.L1.inputs_used.mutations[1]`)
- **TP53**: allele-complete with `chrom="17" pos=7577120 ref="G" alt="A"` (`levels.L1.inputs_used.mutations[2]`)

Tumor context (`levels.L1.inputs_used.tumor_context`):
- `completeness_score = 0.55` (`levels.L1.inputs_used.tumor_context.completeness_score`)
- `msi_status = "MSS"` (`levels.L1.inputs_used.tumor_context.msi_status`)
- `pd_l1_status = "POSITIVE"`, `pd_l1_cps = 10` (`levels.L1.inputs_used.tumor_context.pd_l1_status`, `levels.L1.inputs_used.tumor_context.pd_l1_cps`)
- `er_status = "WEAKLY_POSITIVE"`, `er_percent = 50` (`levels.L1.inputs_used.tumor_context.er_status`, `levels.L1.inputs_used.tumor_context.er_percent`)

### 1.3 L1 completeness / missing data

- `levels.L1.completeness.completeness_score = 0.55` (`levels.L1.completeness.completeness_score`)
- Missing fields listed explicitly: `["HRD score","TMB score","RNA expression data","CA-125 lab values"]` (`levels.L1.completeness.missing`)
- Confidence cap stated: `levels.L1.completeness.confidence_cap = 0.6` (`levels.L1.completeness.confidence_cap`)

### 1.4 WIWFM efficacy panel summary (what it actually scored)

- Drug rows are present under `levels.L1.efficacy.drugs` (10 rows). (`levels.L1.efficacy.drugs`)
- Evidence fields are empty/zeroed across those rows:
  - `evidence_strength = 0.0` (`levels.L1.efficacy.drugs[*].evidence_strength`)
  - `citations = []` (`levels.L1.efficacy.drugs[*].citations`)
  - `citations_count = 0` (`levels.L1.efficacy.drugs[*].citations_count`)
  - `evidence_manifest.pubmed_query = null` (`levels.L1.efficacy.drugs[*].evidence_manifest.pubmed_query`)
  - `ruo_reasons` includes `"no-citations"` (`levels.L1.efficacy.drugs[*].ruo_reasons`)

Sequence scoring provenance for efficacy (`levels.L1.efficacy.provenance.sequence_scoring`):
- `engine_used = "evo2"` (`levels.L1.efficacy.provenance.sequence_scoring.engine_used`)
- `count = 1` (`levels.L1.efficacy.provenance.sequence_scoring.count`)
- `variants_sent_to_engine`: only TP53 (`levels.L1.efficacy.provenance.sequence_scoring.variants_sent_to_engine`)
- `variants_excluded`: MBD4 + PDGFRA with `reason="missing_chrom_or_pos"` (`levels.L1.efficacy.provenance.sequence_scoring.variants_excluded`)

### 1.5 Synthetic Lethality (SL) panel summary

SL status + headline:
- `synthetic_lethality_detected = true` (`levels.L1.synthetic_lethality.synthetic_lethality_detected`)
- `double_hit_description = "Base Excision Repair pathway loss"` (`levels.L1.synthetic_lethality.double_hit_description`)
- `suggested_therapy = "Ceralasertib"` (`levels.L1.synthetic_lethality.suggested_therapy`)

Broken pathways (`levels.L1.synthetic_lethality.broken_pathways`):
- `BER`: `status="non_functional"`, `genes_affected=["MBD4"]`, `disruption_score=0.654` (`levels.L1.synthetic_lethality.broken_pathways[?pathway_id=="BER"]`)
- `CHECKPOINT`: `status="compromised"`, `genes_affected=["TP53"]`, `disruption_score=0.55` (`levels.L1.synthetic_lethality.broken_pathways[?pathway_id=="CHECKPOINT"]`)
- `UNKNOWN`: `status="functional"`, `genes_affected=["PDGFRA"]`, `disruption_score=0.35` (`levels.L1.synthetic_lethality.broken_pathways[?pathway_id=="UNKNOWN"]`)

Essential backup pathways (`levels.L1.synthetic_lethality.essential_pathways`):
- `ATR` and `WEE1` appear with `disruption_score=0.15` and descriptions that explicitly mention “DepMap lineage grounding … +0.15 confidence boost.” (`levels.L1.synthetic_lethality.essential_pathways[?pathway_id=="ATR"]`, `levels.L1.synthetic_lethality.essential_pathways[?pathway_id=="WEE1"]`)
- `HR` and `PARP` appear with `disruption_score=0.0` (`levels.L1.synthetic_lethality.essential_pathways[?pathway_id=="HR"]`, `levels.L1.synthetic_lethality.essential_pathways[?pathway_id=="PARP"]`)

SL recommended drugs (`levels.L1.synthetic_lethality.recommended_drugs`):
- `Ceralasertib` (ATR) `confidence=0.85` (`levels.L1.synthetic_lethality.recommended_drugs[0]`)
- `Adavosertib` (WEE1) `confidence=0.85` (`levels.L1.synthetic_lethality.recommended_drugs[1]`)
- `Olaparib` (HR) `confidence=0.7` (`levels.L1.synthetic_lethality.recommended_drugs[2]`)
- `Niraparib` (HR) `confidence=0.7` (`levels.L1.synthetic_lethality.recommended_drugs[3]`)
- `Rucaparib` (PARP) `confidence=0.7` (`levels.L1.synthetic_lethality.recommended_drugs[4]`)

SL provenance contains auditable receipts:
- `provenance.status = "ok"` (`levels.L1.synthetic_lethality.provenance.status`)
- `true_scoring_required = true` (`levels.L1.synthetic_lethality.provenance.true_scoring_required`)
- `hgvs_resolution[]` receipts include a left-pad deletion normalization note for MBD4 (`levels.L1.synthetic_lethality.provenance.hgvs_resolution[*]`)
- `sequence_scoring.variants_sent_to_engine` includes Evo2 receipts for MBD4/PDGFRA/TP53 (`levels.L1.synthetic_lethality.provenance.sequence_scoring.variants_sent_to_engine[*].evo2_receipt`)

---

## 2) Tumor Board Scribe — 18 precise questions an ovarian oncologist might ask (with “why this question arises”)

These are questions you should expect **based on what this specific bundle shows**.

### A) Why are specific drugs ranked as they are?

1) **Why are paclitaxel and adavosertib at 0.32 while most others are 0.12?**
   - Why it arises: the table contains two efficacy bands (0.32 vs 0.12) without an explicit explanatory flag.
   - Supporting paths: `levels.L1.efficacy.drugs[*].name`, `levels.L1.efficacy.drugs[*].efficacy_score`

2) **Why is olaparib’s efficacy 0.12 even though the PARP germline gate says FULL_EFFECT?**
   - Why it arises: the gate indicates “no penalty,” yet the headline efficacy is still low.
   - Supporting paths: `levels.L1.efficacy.drugs[?name=="olaparib"].efficacy_score`, `levels.L1.efficacy.drugs[?name=="olaparib"].sporadic_gates_provenance`

3) **Which variants were actually scored by the WIWFM sequence engine in this run?**
   - Why it arises: ranking can only reflect variants that were scored.
   - Supporting paths: `levels.L1.efficacy.provenance.sequence_scoring.variants_sent_to_engine`, `...engine_used`, `...count`

4) **Why were MBD4 and PDGFRA excluded from WIWFM sequence scoring?**
   - Why it arises: they are present in `inputs_used.mutations` but absent from the engine inputs.
   - Supporting paths: `levels.L1.inputs_used.mutations`, `levels.L1.efficacy.provenance.sequence_scoring.variants_excluded[*]`

5) **What pathway(s) contributed to the WIWFM panel, and are we essentially seeing a TP53-only pathway model?**
   - Why it arises: pathway_scores appear limited to TP53 (`{"tp53":0.4}`).
   - Supporting paths: `levels.L1.efficacy.pathway_scores`, `levels.L1.efficacy.provenance.sequence_scoring.variants_sent_to_engine`

6) **Did any drug meet an evidence gate (or are all “insufficient” because evidence is 0)?**
   - Why it arises: `evidence_strength=0.0` and `citations_count=0` appear across drugs.
   - Supporting paths: `levels.L1.efficacy.drugs[*].meets_evidence_gate`, `...evidence_strength`, `...citations_count`, `...evidence_tier`

7) **Is “evidence enabled” but returning empty, or is evidence not being queried?**
   - Why it arises: the run flags `evidence_enabled=true`, but `pubmed_query=null` and citations are empty.
   - Supporting paths: `levels.L1.efficacy.provenance.flags.evidence_enabled`, `levels.L1.efficacy.drugs[*].evidence_manifest.pubmed_query`

8) **Were “insights” included, or skipped due to fast mode?**
   - Why it arises: insights can modulate confidence and interpretation; this bundle says it was skipped.
   - Supporting paths: `levels.L1.efficacy.provenance.insights`

### B) How did SL arrive at dependencies and recommendations?

9) **What exact criteria triggered `synthetic_lethality_detected = true`?**
   - Why it arises: tumor board will ask “what rule are you using?”
   - Supporting paths: `levels.L1.synthetic_lethality.synthetic_lethality_detected`, `levels.L1.synthetic_lethality.provenance.decision_rule.sl_detected_criteria`

10) **Which pathways are considered broken, and what are their disruption scores and gene drivers?**
   - Why it arises: this is the core mechanistic claim.
   - Supporting paths: `levels.L1.synthetic_lethality.broken_pathways[*].pathway_id`, `...status`, `...disruption_score`, `...genes_affected`

11) **Why are ATR and WEE1 listed as essential backups, and what is the evidence for that in this bundle?**
   - Why it arises: ATR/WEE1 are the top recommended axis in SL.
   - Supporting paths: `levels.L1.synthetic_lethality.essential_pathways[*].pathway_id`, `...description`, `...disruption_score`

12) **Why are Ceralasertib and Adavosertib recommended at confidence 0.85, and what is the basis of that confidence?**
   - Why it arises: tumor board will challenge the confidence number.
   - Supporting paths: `levels.L1.synthetic_lethality.recommended_drugs[*].drug_name`, `...confidence`, `...rationale`

13) **Does SL use rules, a model, or both — and does it state its detection method?**
   - Why it arises: “black box” vs “rule-based” matters.
   - Supporting paths: `levels.L1.synthetic_lethality.provenance.detection_method`, `...signals_used`

14) **Did SL truly score the variants with Evo2, and can we audit the receipts?**
   - Why it arises: this is the “proof chain” requirement.
   - Supporting paths: `levels.L1.synthetic_lethality.provenance.sequence_scoring.variants_sent_to_engine[*].evo2_receipt`

15) **How were HGVS-only variants normalized into genomic alleles for Evo2 scoring?**
   - Why it arises: MBD4 is an indel; normalization can change alleles.
   - Supporting paths: `levels.L1.synthetic_lethality.provenance.hgvs_resolution[*].receipt`, `...note`

### C) What data is missing, and what should be gathered next?

16) **What does the bundle itself claim is missing at L1, and what confidence cap does it apply?**
   - Why it arises: the system is explicit about missing fields and caps.
   - Supporting paths: `levels.L1.completeness.missing`, `levels.L1.completeness.confidence_cap`, `levels.L1.completeness.completeness_score`

17) **Which tests are explicitly recommended, and what downstream analyses do they unlock?**
   - Why it arises: this maps to next clinical actions.
   - Supporting paths: `tests_needed[*].test`, `tests_needed[*].unlocks`, `tests_needed[*].why`

18) **Which biomarkers are present today (PD-L1, ER, MSI) and which are absent (HRD/TMB/CA-125/expression)?**
   - Why it arises: clarifies what conclusions can/can’t be drawn today.
   - Supporting paths: `levels.L1.inputs_used.tumor_context.*`, `levels.L1.completeness.has_hrd`, `...has_tmb`, `...has_ca125`, `...has_expression`

---

## 3) Clarifier answers — full explanations using only what the payload can prove

This section answers likely follow-ups **as you would speak in a tumor board**, but without importing external facts. When the payload can’t answer, we say what’s missing.

### 3.1 Why are paclitaxel/adavosertib 0.32 while most others are 0.12?

From the bundle, we can say only this:
- The efficacy table contains at least two efficacy-score clusters: `0.32` (e.g., paclitaxel/adavosertib) and `0.12` (many others).
  - Evidence: `levels.L1.efficacy.drugs[*].efficacy_score`, `levels.L1.efficacy.drugs[*].name`
- The bundle does **not** provide an explicit per-drug formula or a “why this drug got 0.32” trace beyond the per-drug `rationale` components and gate provenance.
  - Evidence: `levels.L1.efficacy.drugs[*].rationale`, `levels.L1.efficacy.drugs[*].sporadic_gates_provenance`

**Unknown from this payload**: the precise rule that creates the 0.32 vs 0.12 split.  
**Would need**: an explicit score breakdown (e.g., `levels.L1.efficacy.drugs[*].score_breakdown`) or the scoring equation with intermediate values exposed (not present).

### 3.2 Why is olaparib 0.12 even though PARP_GERMLINE is FULL_EFFECT?

The payload supports the following:
- The PARP gate logic applied to olaparib explicitly says `PARP_GERMLINE` with `verdict="FULL_EFFECT"` and `penalty=1.0`, i.e., it did **not** penalize PARP.
  - Evidence: `levels.L1.efficacy.drugs[?name=="olaparib"].sporadic_gates_provenance`
- The drug’s returned S/P/E rationale includes component values showing:
  - sequence percentile is `0.8`
  - pathway percentile is `0.0`
  - evidence strength is `0.0`
  - Evidence: `levels.L1.efficacy.drugs[?name=="olaparib"].rationale`

So, the payload itself explains “low efficacy” in terms of returned component values (no evidence + zero pathway percentile) despite “no PARP penalty.”

**Unknown from this payload**: why pathway percentile is 0.0 for olaparib in this run.  
**Would need**: an explicit pathway-to-drug mapping trace for olaparib in this run (not present) or explicit per-drug intermediate scoring fields.

### 3.3 Which variants actually drove WIWFM scoring here?

From the bundle:
- `engine_used="evo2"`, `count=1`, and `variants_sent_to_engine` lists TP53 only.
  - Evidence: `levels.L1.efficacy.provenance.sequence_scoring.engine_used`, `...count`, `...variants_sent_to_engine`
- MBD4 and PDGFRA were excluded due to `missing_chrom_or_pos`.
  - Evidence: `levels.L1.efficacy.provenance.sequence_scoring.variants_excluded[*].reason`
- Pathway scores show only `{"tp53":0.4}`.
  - Evidence: `levels.L1.efficacy.pathway_scores`

### 3.4 Why are there zero citations and `evidence_strength=0.0` for all drugs?

From the bundle:
- Every drug row shows `evidence_strength=0.0`, `citations=[]`, `citations_count=0`.
  - Evidence: `levels.L1.efficacy.drugs[*].evidence_strength`, `...citations`, `...citations_count`
- Evidence manifest has `pubmed_query=null` and `citations=[]`.
  - Evidence: `levels.L1.efficacy.drugs[*].evidence_manifest.pubmed_query`, `...evidence_manifest.citations`
- Run flags show `evidence_enabled=true`.
  - Evidence: `levels.L1.efficacy.provenance.flags.evidence_enabled`

**Unknown from this payload**: whether evidence retrieval failed, was skipped, or is stubbed, because the payload has no explicit evidence-provider status/error object.  
**Would need**: an evidence provenance object such as `levels.L1.efficacy.provenance.evidence_provider_status` / `...evidence_errors` / a structured `evidence_run_log` (not present).

### 3.5 How does SL arrive at ATR/WEE1 dependencies and why are ATR/WEE1 top recommendations?

This bundle provides a self-contained chain (as encoded by fields):
- SL detection is on: `synthetic_lethality_detected=true`.
  - Evidence: `levels.L1.synthetic_lethality.synthetic_lethality_detected`
- SL labels BER as non-functional driven by MBD4, and CHECKPOINT as compromised driven by TP53.
  - Evidence: `levels.L1.synthetic_lethality.broken_pathways[*]`
- SL lists ATR and WEE1 under `essential_pathways` with descriptions stating they are dependencies due to checkpoint loss, including “DepMap lineage grounding … +0.15 confidence boost.”
  - Evidence: `levels.L1.synthetic_lethality.essential_pathways[*].pathway_id`, `...description`, `...disruption_score`
- The top SL recommended drugs match those pathways:
  - `Ceralasertib` targets `ATR` with `confidence=0.85`
  - `Adavosertib` targets `WEE1` with `confidence=0.85`
  - Evidence: `levels.L1.synthetic_lethality.recommended_drugs[*].drug_name`, `...target_pathway`, `...confidence`
- Suggested therapy is `Ceralasertib`.
  - Evidence: `levels.L1.synthetic_lethality.suggested_therapy`

SL is also more “input-complete” than WIWFM in this bundle:
- SL sets `true_scoring_required=true` and includes HGVS→GRCh38 receipts and Evo2 receipts for all three genes in its provenance.
  - Evidence: `levels.L1.synthetic_lethality.provenance.true_scoring_required`, `...hgvs_resolution[*]`, `...sequence_scoring.variants_sent_to_engine[*].evo2_receipt`

---

## 4) Tumor-board note (paste-ready; RUO; bundle-only)

**Patient/Bundle**: Ayesha (`patient_id="AK"`), Therapy-Fit L1 bundle (`contract_version="v2.0"`) generated at `generated_at="2026-02-12T16:51:53.836997"`.  
Evidence: `patient_id`, `contract_version`, `generated_at`, `levels.L1`

### Tumor biology represented in this payload

- L1 completeness is `0.55` with an explicit confidence cap of `0.6`, and missing HRD score, TMB score, RNA expression data, and CA-125 lab values.  
  - Evidence: `levels.L1.completeness.completeness_score`, `levels.L1.completeness.confidence_cap`, `levels.L1.completeness.missing`
- Mutations listed include MBD4 (HGVS-only), PDGFRA (HGVS-only), and TP53 R175H with allele-complete coordinates.  
  - Evidence: `levels.L1.inputs_used.mutations[*]`
- Tumor context includes MSS, PD-L1 positive with CPS 10, and weakly positive ER at 50%.  
  - Evidence: `levels.L1.inputs_used.tumor_context.*`

### Synthetic lethality signal and mechanistic rationale (bundle-defined)

- SL is detected with BER-pathway loss described in `double_hit_description`.  
  - Evidence: `levels.L1.synthetic_lethality.synthetic_lethality_detected`, `levels.L1.synthetic_lethality.double_hit_description`
- Broken pathways include BER (MBD4) and CHECKPOINT (TP53).  
  - Evidence: `levels.L1.synthetic_lethality.broken_pathways[*]`
- Essential backup pathways include ATR and WEE1 (`disruption_score=0.15`) with descriptions explicitly mentioning “DepMap lineage grounding … +0.15 confidence boost.”  
  - Evidence: `levels.L1.synthetic_lethality.essential_pathways[*]`
- SL recommended therapies prioritize ATR/WEE1 axis:
  - Ceralasertib (ATR) `confidence=0.85`
  - Adavosertib (WEE1) `confidence=0.85`
  - Suggested therapy: Ceralasertib  
  - Evidence: `levels.L1.synthetic_lethality.recommended_drugs[*]`, `levels.L1.synthetic_lethality.suggested_therapy`

### WIWFM efficacy panel limitations in this payload (why ranking is constrained)

- The efficacy panel’s sequence scoring is bounded by variants sent to the engine; use `sequence_scoring.count`, `variants_sent_to_engine[]`, and `variants_excluded[]` to see what was scored vs. skipped.
  - Evidence: `levels.L1.efficacy.provenance.sequence_scoring.*`
- Evidence/literature is absent in this efficacy panel (`evidence_strength=0.0`, `citations_count=0`, `pubmed_query=null`) even though the flags say evidence is enabled.
  - Evidence: `levels.L1.efficacy.drugs[*].evidence_strength`, `...citations_count`, `...evidence_manifest.pubmed_query`, `levels.L1.efficacy.provenance.flags.evidence_enabled`

### Follow-up tests explicitly recommended by the bundle

- HRD assay.  
  - Evidence: `tests_needed[0].*`
- Comprehensive genomic profiling for TMB.  
  - Evidence: `tests_needed[1].*`
- RNA sequencing/transcriptome.  
  - Evidence: `tests_needed[2].*`
- Serum CA-125 baseline.  
  - Evidence: `tests_needed[3].*`

---

## 5) Practical “how to answer in the room” crib notes (bundle-only)

### If asked “why do we believe the SL panel more than the WIWFM panel here?”

You can say (bundle-only):
- SL includes HGVS→GRCh38 resolution receipts and scored MBD4/PDGFRA/TP53 (`true_scoring_required=true` + `hgvs_resolution[]` + `variants_sent_to_engine[]`).
  - Evidence: `levels.L1.synthetic_lethality.provenance.true_scoring_required`, `...hgvs_resolution[*]`, `...sequence_scoring.variants_sent_to_engine[*]`
- WIWFM may exclude variants from sequence scoring depending on normalization/availability; quote the actual `variants_excluded[*].reason` and `count` from the run you’re presenting.
  - Evidence: `levels.L1.efficacy.provenance.sequence_scoring.variants_excluded[*].reason`, `levels.L1.efficacy.provenance.sequence_scoring.count`

### If asked “why are there no citations?”

You can say (bundle-only):
- This run returns `citations_count=0` and `pubmed_query=null` for all efficacy drug rows, and RUO reasons include “no-citations.”
  - Evidence: `levels.L1.efficacy.drugs[*].citations_count`, `levels.L1.efficacy.drugs[*].evidence_manifest.pubmed_query`, `levels.L1.efficacy.drugs[*].ruo_reasons`
- The payload does not include an evidence-provider status/error object; we can’t say whether evidence retrieval failed or was skipped.
  - Missing field example: `levels.L1.efficacy.provenance.evidence_provider_status` (not present)

