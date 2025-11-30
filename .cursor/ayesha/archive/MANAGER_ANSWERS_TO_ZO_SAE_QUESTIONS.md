# 🧭 MANAGER ANSWERS TO ZO – SAE OPERATIONAL PLAYBOOK (AUTHORITATIVE)

Date: January 13, 2025  
Owner: SR (Manager)  
Context: Align Section 19 (SAE→Evo2→S/P/E) with non‑hallucinatory delivery for Ayesha (Stage IVB ovarian)

---

## Priority Answers (Top 5)

### P1) Missing NGS – What do we show TODAY?
- Policy: Show proactive, deterministic guidance only. No SAE‑driven drug inferences until tumor data exists.
- TODAY deliverables for Ayesha:
  - SOC card (carboplatin + paclitaxel; bevacizumab rationale for ascites/peritoneal).
  - Trials list (frontline, NYC) with transparent eligibility reasoning and confidence gates.
  - CA‑125 monitoring plan (cycle‑3 ≥70% drop, cycle‑6 ≥90%, target <35; resistance: on‑therapy rise or <50% at cycle‑3).
  - Next‑test recommender (order of operations): 1) HRD (MyChoice/tissue), 2) ctDNA for MSI/TMB + somatic HRR, 3) SLFN11 IHC (PARP sensitivity), 4) ABCB1 proxy if prior taxane becomes relevant.
  - Mechanism Map: hidden (or grey “Awaiting NGS”). Hint tiles limited to test order, trials lever, monitoring; no “try/avoid” drug hints yet.

### P2) SAE thresholds – source and policy
- Bands not single points; avoid brittle hard gates.
- Thresholds (calibrated bands; use hysteresis ±0.05 to avoid flapping):
  - dna_repair_capacity: high ≥0.70; moderate 0.40–0.69; low <0.40.
  - pathway_burden.mapk/pi3k/vegf: high ≥0.70; moderate 0.40–0.69; low <0.40.
  - essentiality_signal: high ≥0.80 (stricter); moderate 0.50–0.79.
  - cross_resistance_risk: high >0.70; moderate 0.40–0.70.
  - cohort_overlap: high ≥0.70; moderate 0.40–0.69; low <0.40.
- Sources: literature anchors (GOG‑218/ICON7, PAOLA‑1), internal calibration on retrospective cases, oncologist consensus. Treat as policy constants; log for future recalibration.

### P3) Gemini trial tagging – reliability policy
- Offline only; never in runtime paths. Validation protocol:
  - Batch tag 200 ovarian trials → human spot‑review 30 diverse trials.
  - Accept batch if ≥90% tag accuracy; otherwise adjust prompt taxonomy and re‑tag.
  - Persist `model`, `version`, `parsed_at`, `reviewed_by`, `source_checksum` with each record.
  - Update cadence: weekly diff for new/changed trials. Uncertain tags default to neutral vector; never force a mechanism label.

### P4) Mechanism fit vs eligibility – tiebreak rules
- Ranking: score = eligibility α=0.7 + mechanism_fit β=0.3 (conservative weighting).
- Guardrails:
  - Minimum eligibility threshold to enter top‑10: ≥0.60.
  - Minimum mechanism_fit for mechanism‑gated display: ≥0.50; if <0.50, show but without mechanism boost and add “low mechanism fit” warning.
  - Never suppress SOC; SOC card remains first‑class.
  - Provide “Show all trials” toggle for clinician control.

### P5) Hint tile language – tone policy
- Use suggestive, RUO‑appropriate tone (avoids paternalism; increases adoption):
  - “Consider ordering HRD (impacts PARP eligibility).”
  - “Consider PARP + bevacizumab (ascites/peritoneal).”
  - “Consider avoiding re‑taxane (cross‑resistance risk).”
- Max 4 tiles; priority order: Next test → Trials lever → Monitoring → Avoid (only if applicable). No “avoid” tile for treatment‑naive.

---

## Claim‑Level Answers (Operationalized)

### C1) DNA_repair_capacity high + ascites → platinum ± bevacizumab; PARP maintenance if HRD≥42; ATR/CHK1 on resistance
- Feature definition: dna_repair_capacity = 0.6×DDR_burden + 0.2×essentiality_signal (if HRR gene) + 0.2×exon_disruption (if pathogenic HRR variant present). Else fall back to DDR_burden.
- Banding behavior: High ⇒ display clinician hint and allow modest trial boost for platinum ± bevacizumab arms (+0.10); moderate ⇒ hint only; low ⇒ no hint.
- Resistance detection: trigger when any two of:
  - HRD drop ≥10 points vs baseline,
  - dna_repair_capacity decrease ≥0.15 vs baseline,
  - CA‑125 <50% drop by cycle‑3 or on‑therapy rise.
- TODAY (no NGS): Show “Awaiting HRD; consider ordering HRD to unlock PARP gating.” No SAE‑driven drug confidence yet.
- Finding ATR/CHK1 trials: use curated keywords (“ATR inhibitor”, “CHK1 inhibitor”) + offline MoA tagging; never depend solely on LLM without review.

### C2) RAS/MAPK hotspot or high MAPK burden → MEK/RAF trial candidates; deprioritize MEK monotherapy if absent
- Hotspot detection: use COSMIC/hardcoded list (e.g., KRAS G12C/G12D/G12V, NRAS Q61, BRAF V600E). SAE `hotspot_mutation` may assist but cannot override COSMIC.
- Conflict policy: hotspot present but MAPK burden low (<0.40) ⇒ show trials but no monotherapy boost; combos acceptable if other pathway rationale exists. Boost only if burden ≥0.40; full boost at ≥0.70.
- Deprioritize MEK monotherapy when burden <0.40 (−0.15) and show caution copy.
- TODAY: show “MAPK status: awaiting NGS”; do not surface MEK/RAF levers yet.

### C3) Essentiality high (DDR) → strengthens PARP case; HR restoration pattern → preemptive ATR/CHK1
- essentiality_signal source: Insights/dep‑map style dependency prior; treat as calibrated 0–1.
- Threshold: high ≥0.80; effect: add badge, confidence lift cap +0.03 only (avoid over‑weighting single feature).
- Longitudinal logic: need ≥2 timepoints. HR restoration if HRD drop ≥10 AND dna_repair_capacity drop ≥0.15 OR emergence of RAD51 reactivation signature. Immediate alert; does not wait for radiographic progression.
- For treatment‑naive: show “No longitudinal signal yet; set baseline and re‑assess at week‑12/24.”

### C4) Cross_resistance_risk high with prior taxane/ABCB1 → avoid substrates; propose non‑substrates
- Risk model: max(
  - 0.8 if prior taxane with progression ≤6mo,
  - 0.7 if ABCB1 CNV>4 or expression proxy high,
  - else 0.3 baseline).
- ABCB1 inference before expression: use CNV if available; otherwise UNKNOWN (do not infer).
- Substrate lists: use PharmGKB/DrugBank curated classes; versioned; RUO.
- Non‑substrates: platinum, PARP, ATR/CHK1/WEE1 (verify per‑agent substrate status; do not assume). Provide “likely non‑substrate” tag with source.
- TODAY (naive): show “Cross‑resistance: low/unknown; no avoidance guidance yet.”

### C5) Cohort_overlap low + confidence low → push trials; high overlap → lean standard with modest lift
- Definition: overlap of patient phenotype with literature/trial cohorts (disease, line, key biomarkers); proxy until cohort DB exists.
- Computation policy:
  - High (≥0.70): disease + key biomarker archetype well‑represented (e.g., HRD‑high HGSOC). Add +0.05 confidence and badge.
  - Moderate (0.40–0.69): no lift, no “push trials” banner.
  - Low (<0.40): show “clinical trial recommended” banner; keep SOC but rank trials more prominently.
- Without cohort DB: use policy proxies and explicitly label as proxy.

### C6) Next‑test recommender (choose ONE first)
- Trigger: completeness L0/L1 or missing any of HRD/MSI/TMB.
- Priority order (Ayesha): 1) HRD (PARP gate), 2) ctDNA MSI/TMB + somatic HRR (IO and DDR combo considerations), 3) SLFN11 IHC (PARP sensitivity), 4) ABCB1 proxy if post‑taxane scenario emerges.
- Messaging detail: use “differential branches” format (If positive → X; If negative → Y) with turnaround.

### C7) SAE‑aligned trial ranking (mechanism fit)
- Vectors:
  - Patient `sae_mechanism_vector` = [DDR, MAPK, PI3K, VEGF, IO, Efflux] from SAE; L2‑normalize vectors before cosine.
  - Trial `moa_vector` from offline MoA tagging; store versioned; neutral if unknown.
- Fallback: if patient vector all zeros/unknown, mechanism_fit disabled (β=0) and explain “awaiting NGS; eligibility‑only ranking shown.”
- Explanation: show breakdown (“DDR 0.82 × PARP+ATR → 0.95 fit”).
- Wrong MoA handling: human review gate; uncertain trials remain neutral.

### C8) Clinician hint tiles (UI)
- Max 4; prioritize: Next test → Trials lever → Monitoring → Avoid (only when truly applicable).
- Pre‑NGS: test + monitoring + trials lever only. Post‑NGS: enable “try next/avoid” based on SAE features.
- Copy tone: suggestive; include short reasons; link to source or provenance.

### C9) Mechanism Map UI (chips)
- Thresholds: Green ≥0.70; Yellow 0.40–0.69; Gray <0.40.
- IO special: Green if MSI‑H; Gray if unknown; Red if MSI‑S.
- Pre‑NGS: show gray chips with “Awaiting NGS” overlay and tooltip clarifying meaning.

### C10) Pre‑computed care pathways
- On‑demand assembly (not batch). Criteria for “line‑ready ATR combo trials”:
  - Recruiting; Phase II/III; ≤50 miles; mechanism_fit ≥0.60; exclusions manageable.
- Logistics factor: multiply combined score by proximity factor (1.0 ≤10 miles; 0.9 ≤50; 0.7 >100). Never hide a close, good‑fit trial.
- UI: by default collapse low‑fit mechanisms into “Explore more” with clear toggle to show all.

---

## Alignment With Ayesha
- TODAY (no NGS): deliver SOC + trials + CA‑125 + Next‑test recommender; hide mechanism map; hints limited to testing/monitoring/logistics; no SAE‑based efficacy claims.
- WHEN HRD returns: if ≥42, unlock PARP maintenance (with SAE rationale when available); if <42, surface ATR/CHK1 trials; keep RUO labels.
- IF resistance signals emerge: switch path per C1/C3 and show Resistance Playbook tiles.

---

## Provenance & Safety
- RUO labels on all hints; provenance logs include thresholds used, vectors, gating decisions, and data completeness level.
- No runtime LLM calls in clinical paths; offline LLM outputs human‑reviewed and versioned.
- All thresholds and policies are configurable constants; shipped with documentation and sources.

---

Status: APPROVED POLICY. Proceed to implement per debrief, honoring these guardrails.
