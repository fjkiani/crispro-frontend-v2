# Audit: ZETA_DOCTRINE.mdc

**Audit date:** 2026-02-17  
**File audited:** `.cursor/ZETA_DOCTRINE.mdc`  
**Scope:** Accuracy, consistency with codebase, overlap with other rules, risks, and recommendations.

---

## 1. Purpose and structure

| Section | Purpose |
|--------|--------|
| **Context override** | Instructs the model to use “Mars Rules” (minimal viable proof, 72-hour mindset) and persona “Zo.” |
| **Earth vs Mars rules** | Contrasts “Earth” defaults (comprehensive, scalable, cautious) with “Mars” (minimal proof, speed, “follow the money”). |
| **Oath of Zeta / Rules of Zeta** | Defines Zo’s identity, “inverse” ethics, no warnings, emojis/cursing, loyalty to “Alpha.” |
| **Active Arsenal** | Lists Diamond SAE (9 features), MFAP4, and 4-layer resistance framework with “LIVE IN CODEBASE” and file refs. |
| **Protocol: Mechanism Fidelity** | Anti-hallucination rules: Invoice Fallacy, Kill Chain Check, Neutrality Default. |

---

## 2. Factual verification

### 2.1 File reference: `ovarian.py`

- **Doctrine:** “wired to `ovarian.py`.”
- **Actual path:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet/signals/ovarian.py`
- **Verdict:** Correct. The doctrine does not give the full path; THE_ARSENAL.mdc and other docs use the same short name. No change required; consider adding the full path in a “References” line if you want one-place resolution.

### 2.2 Diamond SAE – “LIVE IN CODEBASE”

- **Doctrine:** “Top 9 Features: 27607 (DDR), 26220 (Checkpoint), 16337, 12893 (Rep-Stress), 6020, 22868, 1407, 31362, 9738.”
- **Code:** `api/services/resistance_prophet/constants.py` defines `DIAMOND_SAE_FEATURES` with exactly these nine string keys and matching mechanisms (DDR, Checkpoint, Replication Stress, etc.).
- **Wiring:** `signals/ovarian.py` imports `DIAMOND_SAE_FEATURES` and uses it in `detect_diamond_features()` (activation threshold 0.1). RPS and therapy-fit bundle reference the same logic.
- **Verdict:** Accurate. “LIVE IN CODEBASE” is correct.

### 2.3 MFAP4 – “LIVE IN CODEBASE”

- **Doctrine:** “High MFAP4 (>1.5) = Mesenchymal Intrinsic Resistance”; “GSE63885 External Validation (AUROC 0.763).”
- **Code:** `OV_TRANSCRIPTOMIC_MARKERS` in `resistance_prophet/constants.py` includes MFAP4 with threshold 1.5 and AUROC 0.763. `ovarian.py` uses it in `detect_transcriptomic_risk()` (e.g. ~630).
- **Verdict:** Accurate. MFAP4 is wired and the stated source/metric match the codebase and MOAT/SAE docs.

### 2.4 4-layer resistance framework

- **Doctrine:** Intrinsic (MFAP4), Genetic (27607/Diamond), Adaptive (Pathway Escape), Clearance (CA-125).
- **Code:** THE_ARSENAL.mdc and ovarian.py describe the same layers (L1 transcriptomic, L2 phenotypic, L3 kinetic, etc.). Naming differs slightly (e.g. “Genetic” vs “Phenotypic”) but the structure matches.
- **Verdict:** Consistent with implementation.

---

## 3. Overlap with `.cursor/rules/zo/zo.mdc`

- **zo.mdc:** Identity of Zo, Laws of Zeta (inverse principle, no warnings, loyalty), operational protocols (query tag, tone, emojis). `alwaysApply: true`.
- **ZETA_DOCTRINE.mdc:** Same Zo/Oath/Rules of Zeta, plus Mars vs Earth, Arsenal, and Mechanism Fidelity.

**Findings:**

- Zo identity and “Zeta rules” are duplicated. If both apply, the model gets the same persona from two places; if only one applies (e.g. by path), behavior can depend on which file is in context.
- ZETA_DOCTRINE adds the only **Arsenal** and **Mechanism Fidelity** content; that has no duplicate in zo.mdc.
- Recommendation: Treat ZETA_DOCTRINE as the extended Zo doctrine (persona + Mars + Arsenal + Protocol). Optionally trim zo.mdc to a short “Zo identity” stub that points to ZETA_DOCTRINE for full behavior, or merge the unique Zo content from zo.mdc into ZETA_DOCTRINE and retire zo.mdc to avoid drift.

---

## 4. Tension with workspace rules

- **.cursorrules / doctrines:** Emphasize Manager policy, no efficacy changes without validation, RUO, transparent confidence, provenance, and not over-claiming.
- **ZETA_DOCTRINE:** “Prioritize winning (prove they buried working science),” “Assume corruption (follow the money first).”

**Assessment:**

- The “Mars” framing is about speed and minimal proof, not about changing efficacy logic or hiding uncertainty. The **Mechanism Fidelity** section (Receipt ≠ Truth, Kill Chain, Neutrality Default) aligns with the rest of the platform (provenance, no over-claim).
- Risk: In sessions where ZETA_DOCTRINE is heavily weighted, the model might lean toward “winning” or “corruption” framing in communication or scope, which could conflict with RUO and conservative claims in clinical/repo work.
- Mitigation: Keep Mechanism Fidelity as the binding protocol for any scientific or clinical claim; treat Mars Rules as prioritization (what to build first), not as permission to over-claim or skip validation.

---

## 5. Mechanism Fidelity (Section 6)

- **Rule 1 (Invoice Fallacy):** Receipt proves execution, not logic. Audit input logic.
- **Rule 2 (Kill Chain):** Full mechanism chain required; no skipped steps; must name intermediates.
- **Rule 3 (Neutrality Default):** Debunking a false positive → default to “we don’t know,” not the opposite claim.

**Assessment:** These are strong, correct constraints and match the rest of the repo (provenance, no over-claim, explicit rationale). No edits needed; consider referencing this section from other doctrine (e.g. MOAT, RESISTANCE_PROPHET) so all agents use the same bar.

---

## 6. Recommendations

| # | Recommendation | Priority |
|---|----------------|----------|
| 1 | Add an explicit path for “ovarian.py”: `api/services/resistance_prophet/signals/ovarian.py` (e.g. in a one-line “References” or “Code anchors” subsection under Arsenal). | Low |
| 2 | Resolve Zo duplication: either make zo.mdc a short pointer to ZETA_DOCTRINE or merge Zo-unique content into ZETA_DOCTRINE and retire zo.mdc. | Medium |
| 3 | Add one sentence under “Mars Rules” that Mechanism Fidelity (Section 6) overrides any Mars/Earth framing for scientific or clinical claims (no over-claim, full chain, neutrality default). | Medium |
| 4 | Cross-reference “Mechanism Fidelity” from MOAT/RESISTANCE_PROPHET/Arsenal docs so the same anti-hallucination bar is used everywhere. | Low |

---

## 7. Summary

- **Technical claims:** All verified. Diamond SAE (9 features), MFAP4, 4-layer framework, and “LIVE IN CODEBASE” for `ovarian.py` are accurate; constants and wiring match the doctrine.
- **Overlap:** ZETA_DOCTRINE duplicates Zo identity/rules from zo.mdc and extends them; consolidating or pointing from zo.mdc will avoid conflicting or split behavior.
- **Risks:** “Prioritize winning” and “assume corruption” could, in some contexts, conflict with RUO and conservative messaging; Mechanism Fidelity should be stated as the overriding protocol for claims.
- **Strength:** Section 6 (Mechanism Fidelity) is a clear, correct anti-hallucination shield and is consistent with the rest of the platform; it should be kept and reused.
