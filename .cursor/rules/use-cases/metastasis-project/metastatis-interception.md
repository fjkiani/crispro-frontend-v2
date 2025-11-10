Doctrine: The Metastasis Interception Arsenal
1.0 Mission Briefing
The old guard designs weapons against a primary tumor. We will design weapons against the invasion itself. This doctrine outlines the operational protocol for our CommandCenter to translate a high-level strategic objective (e.g., "Prevent Angiogenesis") into a validated, AI-forged therapeutic blueprint.

This is not about simply finding a gene and designing a guide. This is a multi-stage in silico campaign that fuses our entire arsenal to create stage-specific anti-metastatic weapons.

**CLARIFICATION QUESTIONS (Added during review):**
Q1: How does the user select a "mission" (e.g., "Disrupt Angiogenesis")? Is this a dropdown in UI mapped to the 8 cascade steps, or free-form text that we parse?
Q2: What is "Digital_Twin_ID_007"? Is this a patient identifier we store, or is this conceptual? Do we need a patient/session context system?

2.0 The Kill Chain: From Strategic Intent to Validated Weapon
The CommandCenter will execute the following workflow when commanded to design a metastasis interception therapeutic.

Phase I: Target Validation & Prioritization (The "Why")
Before we forge a weapon, we must be certain we are aiming at the right fucking target.

Receive Strategic Command: The workflow is initiated with a high-level command from the user via the CrisPRO Studio UI. For example:

Mission: Disrupt Angiogenesis

Patient: Digital_Twin_ID_007

Pathway Intelligence: The CommandCenter consults its internal knowledge base to identify the primary biological pathway associated with the mission (e.g., the VEGF signaling pathway).

**Q3: Where is the "internal knowledge base" that maps mission→pathway→gene_sets? Is this the `metastasis_interception_rules.json` config, or a separate KB we need to build?**

Contextual Target Analysis (Unleash the Oracle): This is the fucking kill shot. The CommandCenter does not blindly target VEGF. It first uses the Zeta Oracle to analyze the patient's specific Digital Twin and answer a critical question: "Which gene in the VEGF pathway is the primary driver of angiogenesis in this specific tumor?"

It uses /predict_gene_essentiality to score the importance of each gene in the pathway (VEGFA, VEGFR1, VEGFR2, etc.).

**Q4: For `/predict_gene_essentiality`, what are the exact inputs we need? Do we need gene coordinates, or just gene symbols? What "context" do we provide (cell line, disease, patient mutations)?**

It analyzes the patient's specific mutations in these genes with /predict_variant_impact.

**Q5: Do we currently have `/predict_variant_impact` endpoint? I see we have `/api/insights/predict_protein_functionality_change` - is that the same, or do we need to create a new wrapper?**

Target Lock: The CommandCenter synthesizes this intelligence and achieves a "target lock." It identifies the single most vulnerable and critical gene for that specific patient's angiogenic strategy (e.g., VEGFA).

**Q6: How do we "synthesize" multiple signals (essentiality, variant impact, functionality, regulatory) into a single "target lock" score? What is the weighting formula? Is this in the config as `weights.target_lock`?**
**Q7: What happens if multiple genes tie for top rank? Do we return all of them, or use a tiebreaker (e.g., prefer genes present in patient mutations)?**

Phase II: Weapon Forging & Validation (The "How")
With a validated, high-value target, the CommandCenter now has a green light to forge the weapon.

Initiate Design & Validation Workflow: The CommandCenter now makes an internal call to our established /workflow/design_and_validate_guides endpoint.

**Q8: Does `/workflow/design_and_validate_guides` currently exist? Or do we need to build this as part of interception? If it doesn't exist, should we just call `/api/design/generate_guide_rna` directly?**

Input: It passes the precise genomic locus of the validated target gene (e.g., the coordinates for VEGFA).

**Q9: How do we obtain "precise genomic locus" for a gene? Do we need a gene→coordinates lookup service, or do we assume the user provides coords in `mutations`? What if coords are missing?**
**Q10: What happens if the target gene has multiple isoforms/transcripts? Do we target all exons, or just one canonical transcript?**

Execute the Kill Chain: The /design_and_validate_guides workflow executes its own automated, multi-step campaign:

Forge Candidates: The Zeta Forge is commanded to generate a suite of optimal gRNA candidates.

**Q11: "Zeta Forge" = `/api/design/generate_guide_rna`? What are the exact parameters (target_sequence length, PAM, num candidates, model_id)?**
**Q12: How many candidates do we request? Is this configurable in the interception_rules.json (e.g., `num_candidates_per_target: 5`)?**

Validate Efficacy: Each candidate is passed to the Zeta Oracle to calculate a "Zeta Score," quantifying its predicted functional impact.

**Q13: What is "Zeta Score" in this context? Is this the `efficacy_proxy` from design endpoint, or a separate scoring step using Evo2 delta likelihood?**
**Q14: Do we need to predict the on-target efficacy for each guide? Is there a `/predict_crispr_spacer_efficacy` endpoint, or do we use heuristics (GC content, homopolymer runs)?**

Validate Safety: Each candidate is passed to our BLAST Service for a comprehensive, genome-wide off-target analysis.

**Q15: "BLAST Service" = `api/services/safety_service.py`? What function do we call (e.g., `preview_off_targets`, `assess_safety`)?**
**Q16: How long does off-target analysis take per guide? Do we need timeouts/retries? What if BLAST is unavailable - do we fail gracefully with placeholder safety scores?**

Synthesize the Blueprint: The CommandCenter receives the final, rank-ordered list of fully validated gRNAs. It packages this intelligence into a complete therapeutic blueprint.

**Q17: What is the ranking formula for the final "Assassin Score"? Is it `w_eff*efficacy + w_safe*safety + w_fit*mission_fit`, where mission_fit is based on target_lock score?**
**Q18: What data structure is the "therapeutic blueprint"? Is this just the response schema from `/api/metastasis/intercept`, or do we need to generate a PDF/report artifact?**

3.0 The Deliverable: A Mission-Specific Arsenal
The final output is not a generic list of guides. It is a "Validated Anti-Angiogenic Payload," a dossier containing:

Mission Objective: Disrupt Angiogenesis for Patient Digital_Twin_ID_007.

**Q19: How is "Mission Objective" displayed in the UI? Is this a formatted string, or do we need to generate human-readable text from `mission_step` and patient context?**

Validated Target: VEGFA, identified as the key vulnerability.

**Q20: Do we return JUST the top-ranked target, or do we also show runner-ups (e.g., top 3 targets with their scores) for transparency?**

The Arsenal: A rank-ordered list of gRNA candidates, each with its own:

Sequence

Predicted Efficacy Score (Zeta Score)

Off-Target Safety Score

Composite "Assassin Score"

**Q21: What is the minimum number of candidates we must return? What if design only generates 1 candidate - do we still proceed, or error out?**
**Q22: Do we need to include provenance per-candidate (which model generated it, which safety check ran), or is overall provenance sufficient?**
**Q23: Should candidates include PAM site, GC%, position in target sequence, or just the 20bp guide sequence?**

This doctrine transforms our platform from a simple design tool into a true strategic weapon system. We don't just build weapons; we ensure they are the right fucking weapons for the right fucking target, for every stage of the war. 🚀

---

## SUMMARY OF CLARIFICATION QUESTIONS (18 total)

**Mission & Patient Context (Q1-Q2):**
- How does user select mission? Dropdown or free-form?
- What is Digital_Twin_ID? Do we need patient/session management?

**Target Lock Phase (Q3-Q7):**
- Where is mission→pathway→gene_set mapping stored?
- What inputs does `/predict_gene_essentiality` need?
- Do we have `/predict_variant_impact` or use existing insights endpoints?
- What is the weighting formula for target_lock score?
- How do we handle ties in target ranking?

**Weapon Forging Phase (Q8-Q14):**
- Does `/workflow/design_and_validate_guides` exist, or build new?
- How to obtain genomic locus for target gene?
- How to handle multiple transcripts/isoforms?
- Confirm "Zeta Forge" = `/api/design/generate_guide_rna`?
- How many candidates to request (configurable)?
- What is "Zeta Score" - efficacy_proxy or separate Evo2 scoring?
- Do we have `/predict_crispr_spacer_efficacy` or use heuristics?

**Safety Validation (Q15-Q16):**
- Confirm "BLAST Service" = `safety_service.py`? Which function?
- What are timeouts/retries for off-target analysis?

**Final Blueprint (Q17-Q23):**
- What is the Assassin Score formula?
- Is "therapeutic blueprint" just JSON response or PDF artifact?
- How to format "Mission Objective" text?
- Return top target only, or show runner-ups?
- Minimum number of candidates required?
- Per-candidate provenance vs overall provenance?
- What metadata to include for each candidate (PAM, GC%, position)?

---

## ANSWERS TO CLARIFICATION QUESTIONS

A1: Mission is a fixed dropdown mapped to the 8 cascade steps; free‑form is not supported in v1.

A2: `Digital_Twin_ID_*` is a patient/session identifier string. We will accept an optional `patient_id` and echo it into provenance; storage is handled by the Sessions API/FE `SessionContext` (local persistence acceptable for v1).

A3: Use `metastasis_interception_rules.json` as the internal KB for mission→pathway→gene_sets mapping in v1. No separate KB required.

A4: For `/api/insights/predict_gene_essentiality`, provide `{ gene, variants: [mutation?], model_id }`. Coordinates are optional; if present in `mutations`, include them. Default `model_id = evo2_1b`.

A5: We do not have `/predict_variant_impact`. Use existing insights: `predict_protein_functionality_change` (coding proxy) and `predict_splicing_regulatory` (noncoding proxy). Combine as needed.

A6: Target‑lock score = weighted sum over available signals using `weights.target_lock` in config; enforce per‑signal thresholds before aggregation. Normalize to [0,1].

A7: On ties, pick the gene present in `mutations` first; then higher essentiality; then stable alphabetical. Return the single `validated_target` plus a `considered_targets` array for transparency.

A8: `/workflow/design_and_validate_guides` does not exist. Call `/api/design/generate_guide_rna` directly in v1; keep the workflow name reserved for v2.

A9: Prefer locus from provided `mutations` for the target gene. If absent, use canonical transcript (MANE/longest CDS) from local gene DB if available; else require a user‑supplied 30–100bp `target_sequence` (return a clear error if none available).

A10: Use one canonical transcript by default (MANE/longest CDS). Allow override via request `options.transcript_id` in v2.

A11: Yes—“Zeta Forge” = `POST /api/design/generate_guide_rna` with `{ target_sequence: 30–100bp, pam: "NGG", num: K, model_id: "evo2_1b" }`.

A12: Configure `K` in `metastasis_interception_rules.json` (e.g., `num_candidates_per_target: 3–5`). Default to 3.

A13: In v1, use the design endpoint’s `efficacy_proxy` as the efficacy component. Evo2 delta‑based re‑scoring can be added in v2 as an optional enhancement.

A14: No first‑class `/predict_crispr_spacer_efficacy` yet. Use heuristics already returned (GC%, homopolymer penalties) as the on‑target proxy in v1.

A15: Use `safety_service.preview_off_targets(...)` if available; otherwise call the existing safety preview helper used by design flows. Return `{ off_target_hits, red_flags[] }`.

A16: Set per‑guide safety timeout to 30–60s with 1 retry (exponential backoff). On timeout/failure, return placeholder safety `{ off_target_hits: null, red_flags: ["safety_unavailable"] }` and record in provenance.

A17: `assassin_score = w_eff * efficacy_proxy + w_safe * safety_score + w_fit * mission_fit`, with weights from config. Map safety to [0,1] where fewer predicted off‑targets ⇒ higher `safety_score`.

A18: v1 returns a JSON blueprint only (the `/api/metastasis/intercept` response). PDF/report artifact generation is deferred to v2.

A19: FE derives a human‑readable “Mission Objective” from `mission_step` and optional `patient_id`. BE also returns a formatted string for convenience.

A20: Return one `validated_target` and a `considered_targets` list (top 3) with scores and brief rationale.

A21: Target ≥2 candidates. If only 1 (or 0) are generated, return what’s available with a `provenance.status_warning` and do not error; FE displays a warning.

A22: Provide both overall provenance and per‑candidate provenance: `{ design_method, safety_method/status }`. Keep methods consistent with v1 naming.

A23: Include `{ sequence (20bp spacer), pam, gc, position_in_window?, efficacy_proxy }`. `position_in_window` is reported when the locus window is known.

---

## FOLLOW-UP CLARIFICATIONS NEEDED (7 questions)

**Gene Coordinate Lookup (Critical for A9):**
F1: No comprehensive gene DB exists. For v1, should we:
  - (Option A) **Require coords in mutations** - fail with clear error if target gene lacks coords
  - (Option B) **Hardcode canonical targets** - add VEGFA, MMP2, MMP9, SNAI1, etc. with known exon sequences (~10-15 genes)
  - (Option C) **Use Ensembl API** - add `/api/safety/ensembl_context` call to fetch transcript on-demand (slower, external dependency)
  - **Recommended: Option A** for v1 simplicity; Option B for demo polish; Option C for v2

**Safety Service Integration (A15 verification):**
F2: `safety_service.preview_off_targets()` returns heuristics ONLY (GC, homopolymer, heuristic_score), NOT actual off-target hit counts. Should we:
  - (Option A) Use `heuristic_score` as the safety component (higher = better)
  - (Option B) Return placeholder `{ off_target_hits: null, red_flags: [] }` with note "Genome alignment not available in v1"
  - **Recommended: Option A** - map `heuristic_score` to safety_score directly

**Design Endpoint Return (A13 clarification):**
F3: `/api/design/generate_guide_rna` already returns `spacer_efficacy_heuristic` per candidate. Confirm we use this directly as `efficacy_proxy` in v1?
  - **Verified**: Yes - the design endpoint returns `{ sequence, pam, gc, spacer_efficacy_heuristic }` which we'll use as-is

**Mission Step Gene Sets (A3 expansion):**
F4: The intervention config has gene_sets (MAPK, HRR, EMT_TF, MMP, HOMING, IMMUNE, DORMANCY). Do we need additional gene_sets for interception-specific targeting? E.g.:
  - `ANGIO`: [VEGFA, VEGFR1, VEGFR2, FGF2, etc.] for angiogenesis mission
  - `ADHESION`: [CDH1, ICAM1, VCAM1, etc.] for EMT mission
  - Or reuse existing gene_sets from intervention?
  - **Recommended: Reuse + add ANGIO** for v1; full expansion in v2

**Config Weights (A6/A17 details):**
F5: Confirm default weights in `metastasis_interception_rules.json`:
  - `weights.target_lock`: { functionality: 0.35, essentiality: 0.35, regulatory: 0.15, chromatin: 0.15 }
  - `weights.assassin`: { efficacy: 0.4, safety: 0.3, mission_fit: 0.3 }
  - `thresholds`: { functionality: 0.6, essentiality: 0.6, regulatory: 0.6, chromatin: 0.6 }
  - Are these reasonable? Adjust?

**Mission Fit Calculation (A17 expansion):**
F6: `mission_fit` in assassin_score formula: should this be:
  - (Option A) The target_lock score of the validated target (0.0-1.0)
  - (Option B) Binary: 1.0 if gene in mission gene_set, 0.5 otherwise
  - (Option C) A combination: target_lock * gene_set_match_bonus
  - **Recommended: Option A** - use target_lock score directly as mission_fit

**Candidate Metadata (A23 final detail):**
F7: The design endpoint doesn't return `position_in_window`. Should we:
  - (Option A) Omit this field in v1 response
  - (Option B) Track PAM positions from our window extraction and include it
  - **Recommended: Option A** - defer to v2 when we have precise locus tracking



   CLARIFYING QUESTIONS BY TASK
Task 1: Expand design window to ±300bp
Q1.1: What is the current window size in /api/design/generate_guide_rna? Is it just the 20bp PAM windows?
Q1.2: Should the config knob be in metastasis_interception_rules.json as design.window_size: 300, or a new env var DESIGN_WINDOW_SIZE?
Q1.3: For tests, do you want unit tests for window extraction logic, or integration tests with real gene sequences?
Task 2: Enable design API via feature flag
Q2.1: Is the design API currently disabled? I don't recall seeing any feature flags around it.
Q2.2: Proposed env var: ENABLE_GUIDE_RNA_DESIGN=true? Should it default to false (opt-in) or true (opt-out)?
Q2.3: Where to document - in main README, or create a separate DESIGN_API.md guide?
Task 3: Harden Ensembl fetch
Q3.1: Do we currently have Ensembl integration? I don't see it in the existing codebase.
Q3.2: Is this for fetching gene coordinates/transcripts to enable auto-design when coords are missing? (relates to F1)
Q3.3: Cache TTL and storage - Redis with 24hr TTL, or in-memory LRU with 1hr TTL?
Q3.4: Graceful fallback - return clear error message, or fall back to hardcoded coords from Task 4?
Task 4: Hardcode ANGIO gene exon windows
Q4.1: Which specific genes? Proposed: [VEGFA, VEGFR1, VEGFR2, FGF2, HIF1A, ANGPT1, ANGPT2] - confirm or adjust?
Q4.2: Storage location - add gene_loci section to metastasis_interception_rules.json, or create separate gene_sequences.json?
Q4.3: Data format per gene:
  {
    "VEGFA": {
      "chrom": "6",
      "exon_windows": [
        {"exon": 1, "start": 43737945, "end": 43738245, "sequence": "ACGT...300bp..."}
      ]
    }
  }
  Confirm this structure, or provide alternative?
Task 5: Real off-target search (BLAST/minimap2)
Q5.1: Do we have BLAST/minimap2 installed? Is this a Modal service, local binary, or external API?
Q5.2: Genome reference - use GRCh38.fa, or a pre-indexed reference?
Q5.3: Timeout per guide - 60s as per F2, with 1 retry? Or different?
Q5.4: Safety score mapping formula - proposed: safety_score = exp(-0.1 * off_target_hits) where 0 hits = 1.0, 10 hits ≈ 0.37? Adjust?
Task 6: Add /api/design/predict_crispr_spacer_efficacy
Q6.1: Which algorithm? Options:
(A) Evo2 delta scoring (use existing Evo proxy, most consistent with S/P/E)
(B) Rule-based heuristics (GC%, position, homopolymers - already done)
(C) DeepSpCas9/Azimour (external model, more accurate but new dependency)
Recommended: Option A for v1
Q6.2: Should this replace spacer_efficacy_heuristic in response, or add as separate field evo2_efficacy_score?
Q6.3: Inputs - just guide_sequence (20bp), or guide_sequence + target_context (±100bp)?
Task 7: Polish FE interception UX
Q7.1: "Disable without coords" - should we:
(A) Hide the entire interception panel?
(B) Show panel but disable "Design Guides" button with tooltip: "Requires variant coordinates"?
Recommended: Option B for transparency
Q7.2: Per-candidate provenance detail level:
(A) Just method names: { design_method: "evo2_1b", safety_method: "heuristic" }
(B) Full provenance: { design_run_id, safety_run_id, timestamps }
Recommended: Option A for v1
Q7.3: Toast library - do we already have one installed? If not, prefer react-hot-toast or sonner?
Task 8: E2E test (variant → target + candidates)
Q8.1: E2E framework - Playwright (recommended), Cypress, or manual bash script?
Q8.2: Test against local backend (port 8000) or deployed staging?
Q8.3: Expected duration - is <2 min acceptable, or need <30s?
Task 9: Containerize & deploy to staging
Q9.1: Where is "staging" - Modal (already using it), AWS ECS, GCP Cloud Run, or local Docker Compose?
Q9.2: Do we have Dockerfiles? If not, should I create:
oncology-backend-minimal/Dockerfile
oncology-frontend/Dockerfile
Root docker-compose.yml for orchestration?
Q9.3: Required env vars for staging:
Backend: EVO_URL_1B, ENABLE_GUIDE_RNA_DESIGN, ENSEMBL_API_URL?
Frontend: VITE_API_ROOT?
Others?
Task 10: Update docs & demo script
Q10.1: Which docs to update:
(A) Update existing metastatis-interception.md with implementation notes?
(B) Create new METASTASIS_INTERCEPTION_USER_GUIDE.md?
(C) Both?
Recommended: Option C
Q10.2: Demo script format:
(A) Bash script with curl commands (smoke_test_interception.sh)
(B) Jupyter notebook with narrative walkthrough
(C) Video script with UI screenshots
Recommended: Option A + Option C
Q10.3: Figures needed:
F1: Architecture diagram (mission → target lock → design → safety → assassin score)
F2: UI screenshot (MetastasisInterceptionPanel with results)
F3: Sample output JSON (validated_target + candidates with scores)
Confirm or add/remove?
Q10.4: v2 roadmap location - add section to metastatis-interception.md, or create ROADMAP_V2.md?
