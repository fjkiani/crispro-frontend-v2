## Evo2 Endpoint Refactor Plan

### Why this document
- We need a clean, reliable way to expose Evo2-powered capabilities without bloating the Command Center or blowing deployment quotas (OOM on Vercel).
- This plan standardizes service names, defines endpoints grounded in the Evo2 paper, and describes a thin-proxy backend pattern for Vercel.

### Current pain points
- Heavy monoliths (full backend + models) cause OOM during Vercel builds.
- Command Center has accumulated orchestration + business logic + compute that belong in specialized services.
- Modal URLs are inconsistent and hard to maintain.

### Target architecture (clean separation of concerns)
- Command Center (FastAPI on Modal): Orchestrator only
  - Accepts high-level workflow requests (e.g., “Target Dossier”)
  - Calls specialized services (Oracle/Evo2, Forge, Gauntlet, Boltz, Fusion)
  - Fuses results and returns a single response
- Specialized compute microservices (Modal):
  - Evo Service (Evo2-based scoring/generation)
  - Zeta Oracle (thin wrapper combining rules with Evo2 scoring)
  - Zeta Forge (CRISPR/inhibitor design, orchestrates Gauntlet/Boltz as needed)
  - Gauntlet (validation and simulated trials)
  - Boltz Service (structure/affinity simulation proxy)
  - Fusion Engine (pipelines/transformers for multi-modal fusion)
- Public Backend (Vercel): Thin proxy
  - Small, stateless endpoints that forward to Modal webhooks
  - No heavy ML deps; no long-running tasks

### Service naming conventions (Modal)
- zeta-oracle → app name: `zeta-oracle`
- zeta-forge → app name: `zeta-forge`
- gauntlet → app name: `gauntlet`
- boltz-service → app name: `boltz-service`
- fusion-engine → app name: `fusion-engine`
- evo-service → app name: `evo-service`
- hunter-analyst → app name: `hunter-analyst`

The resulting webhook URLs follow Modal’s pattern after deploy; example:
- `https://crispro--zeta-oracle-zetaoracle-api.modal.run/invoke`
- `https://crispro--evo-service-evoservice-api.modal.run/score_delta`

### Endpoint suite (grounded in Evo2 paper)

#### Phase 1 (Must-have, low-risk)
- Evo Service (Evo2 scoring)
  - POST `/score_delta`
    - Input: `{ref_sequence, alt_sequence}` or `{assembly, chrom, pos, ref, alt}`
    - Output: `{ref_likelihood, alt_likelihood, delta_likelihood_score}`
    - Notes: 8,192-nt window centered on variant; keep total length constant for indels
  - POST `/score_batch`
    - Input: list of records; batched version of above
- Vercel Backend (proxy only)
  - POST `/api/evo/score_delta` → Modal `/score_delta`
  - POST `/api/evo/score_batch` → Modal `/score_batch`

#### Phase 2 (Oncology R&D core)
- Evo Service (extend)
  - POST `/score_splice` (splicing VUS; Evo2 zero-shot on splice windows)
  - POST `/score_regulatory` (TF motif/regulatory assessments; Evo2 zero-shot on noncoding windows)
  - POST `/embeddings` (for supervised heads; e.g., BRCA1 classifier)
- Oracle (rules + Evo2 fusion)
  - POST `/predict_variant_impact` (wraps Evo2 delta + rules)
  - POST `/predict_gene_essentiality` (KO/scramble proxies + Evo2 signatures)
- Forge
  - POST `/generate_crispr_guides` (design + off-target proxy)
  - POST `/generate_inhibitors` (design + rank; orchestrate Boltz)
- Gauntlet
  - POST `/run_trials` (validate CRISPR/inhibitors; return efficacy metrics)
- Vercel Backend (proxy only)
  - Mirror the above with `/api/*` routes that forward to Modal

#### Phase 3 (Advanced)
- Evo Service
  - POST `/design_accessibility` (inference-time control for chromatin accessibility)
  - POST `/generate_sequence` (long-form generation with safety guardrails)
- Boltz Service
  - POST `/v1/predict_structure`
  - POST `/v1/binding_affinity`

### Request/Response schemas (examples)

#### Evo2 delta scoring
Request
```json
{
  "ref_sequence": "ACGT...",
  "alt_sequence": "ACGt..."
}
```
Response
```json
{
  "ref_likelihood": -1234.56,
  "alt_likelihood": -1248.90,
  "delta_likelihood_score": -14.34
}
```
Notes
- Negative delta → disruptive; magnitudes context-dependent
- Validate `[ACGTN]` only; cap sequences (e.g., <= 16,384)

#### Splicing/regulatory
Request
```json
{
  "assembly": "GRCh38",
  "chrom": "17",
  "pos": 43045736,
  "ref": "A",
  "alt": "G",
  "context": "splice" | "regulatory",
  "window": 8192
}
```
Response
```json
{
  "delta_likelihood_score": -3.21,
  "feature_scores": {
    "exon_boundary": 0.87,
    "tf_motif_disruption": 0.42
  },
  "explainability": {
    "saliency": [ ... ],
    "top_features": ["exon_start", "splice_acceptor"]
  }
}
```

### Implementation notes (Evo2 specifics)
- Windowing: default 8,192-nt centered on variant; pad/trim for indels to keep constant length
- Indels: represent minimal ref/alt; reflow flanks to preserve total length
- Batching: use small batches by token length; 7B for speed, 40B optional for higher fidelity
- Explainability: expose SAE-derived features where available (exon/intron boundaries, TF motifs, mutation severity)

### Deployment

#### Modal (services)
```bash
# One-time auth
venv/bin/modal token new

# Verify CLI
venv/bin/modal --version

# Deploy services (examples)
venv/bin/modal deploy src/services/evo_service/main.py
venv/bin/modal deploy src/services/oracle/main.py
venv/bin/modal deploy src/services/forge/main.py
venv/bin/modal deploy src/services/gauntlet/main.py
venv/bin/modal deploy src/services/boltz_service/main.py
venv/bin/modal deploy src/services/fusion_engine/main.py
```

#### Vercel (thin proxy backend)
- Layout: `oncology-coPilot/oncology-backend-minimal/api/index.py`
- Config: `vercel.json` routes all to `api/index.py`
- Requirements: `requirements-minimal.txt` → copy to `requirements.txt`
```bash
cd oncology-coPilot/oncology-backend-minimal
cp requirements-minimal.txt requirements.txt
npx -y vercel --yes --prod
```

### Environment variables
- In Vercel
  - `EVO_SERVICE_URL` → Modal evo-service base URL
  - Add others as services go live (`ORACLE_URL`, `FORGE_URL`, `GAUNTLET_URL`, etc.)
- In services
  - Use `modal.Secret.from_dotenv()` and `.env` for internal API keys

### Testing (smoke tests)
```bash
# Evo2 delta (Modal)
curl -X POST "$EVO_SERVICE_URL/score_delta" \
  -H 'content-type: application/json' \
  -d '{"ref_sequence":"ACGT", "alt_sequence":"AGGT"}'

# Evo2 delta (Vercel proxy)
curl -X POST "$VERCEL_BASE/api/evo/score_delta" \
  -H 'content-type: application/json' \
  -d '{"ref_sequence":"ACGT", "alt_sequence":"AGGT"}'
```

### Rollout plan
- Phase 1
  - Ship `/score_delta`, `/score_batch` in `evo-service` (Modal)
  - Add `/api/evo/score_delta`, `/api/evo/score_batch` to Vercel proxy
  - Wire Command Center to call Modal directly (keep Vercel for public/demo)
- Phase 2
  - Add splicing/regulatory scoring + embeddings
  - Implement Oracle/Forge/Gauntlet proxy routes
- Phase 3
  - Add generation endpoints with strong safety guardrails
  - Integrate Boltz structure/affinity simulation

### Command Center cleanup (scope reduction)
- Only orchestrate workflows
  - Transform inputs → call services → fuse outputs → return response
- Remove embedded compute/business logic that belongs in Evo2/Forge/Gauntlet
- Centralize downstream URLs as constants; use environment-driven config

### Success criteria
- Vercel backend stays <50MB and deploys in <60s
- Modal services expose stable, standardized URLs
- Command Center has <500 LOC of orchestration code, zero heavy deps
- End-to-end demo runs fully against Modal services with thin Vercel proxy

### Future considerations
- Add caching for repeated variant windows
- Add dataset-backed validations (ClinVar/SpliceVarDB) for CI
- Log schema versions in responses for robust clients 
 

### Why this document
- We need a clean, reliable way to expose Evo2-powered capabilities without bloating the Command Center or blowing deployment quotas (OOM on Vercel).
- This plan standardizes service names, defines endpoints grounded in the Evo2 paper, and describes a thin-proxy backend pattern for Vercel.

### Current pain points
- Heavy monoliths (full backend + models) cause OOM during Vercel builds.
- Command Center has accumulated orchestration + business logic + compute that belong in specialized services.
- Modal URLs are inconsistent and hard to maintain.

### Target architecture (clean separation of concerns)
- Command Center (FastAPI on Modal): Orchestrator only
  - Accepts high-level workflow requests (e.g., “Target Dossier”)
  - Calls specialized services (Oracle/Evo2, Forge, Gauntlet, Boltz, Fusion)
  - Fuses results and returns a single response
- Specialized compute microservices (Modal):
  - Evo Service (Evo2-based scoring/generation)
  - Zeta Oracle (thin wrapper combining rules with Evo2 scoring)
  - Zeta Forge (CRISPR/inhibitor design, orchestrates Gauntlet/Boltz as needed)
  - Gauntlet (validation and simulated trials)
  - Boltz Service (structure/affinity simulation proxy)
  - Fusion Engine (pipelines/transformers for multi-modal fusion)
- Public Backend (Vercel): Thin proxy
  - Small, stateless endpoints that forward to Modal webhooks
  - No heavy ML deps; no long-running tasks

### Service naming conventions (Modal)
- zeta-oracle → app name: `zeta-oracle`
- zeta-forge → app name: `zeta-forge`
- gauntlet → app name: `gauntlet`
- boltz-service → app name: `boltz-service`
- fusion-engine → app name: `fusion-engine`
- evo-service → app name: `evo-service`
- hunter-analyst → app name: `hunter-analyst`

The resulting webhook URLs follow Modal’s pattern after deploy; example:
- `https://crispro--zeta-oracle-zetaoracle-api.modal.run/invoke`
- `https://crispro--evo-service-evoservice-api.modal.run/score_delta`

### Endpoint suite (grounded in Evo2 paper)

#### Phase 1 (Must-have, low-risk)
- Evo Service (Evo2 scoring)
  - POST `/score_delta`
    - Input: `{ref_sequence, alt_sequence}` or `{assembly, chrom, pos, ref, alt}`
    - Output: `{ref_likelihood, alt_likelihood, delta_likelihood_score}`
    - Notes: 8,192-nt window centered on variant; keep total length constant for indels
  - POST `/score_batch`
    - Input: list of records; batched version of above
- Vercel Backend (proxy only)
  - POST `/api/evo/score_delta` → Modal `/score_delta`
  - POST `/api/evo/score_batch` → Modal `/score_batch`

#### Phase 2 (Oncology R&D core)
- Evo Service (extend)
  - POST `/score_splice` (splicing VUS; Evo2 zero-shot on splice windows)
  - POST `/score_regulatory` (TF motif/regulatory assessments; Evo2 zero-shot on noncoding windows)
  - POST `/embeddings` (for supervised heads; e.g., BRCA1 classifier)
- Oracle (rules + Evo2 fusion)
  - POST `/predict_variant_impact` (wraps Evo2 delta + rules)
  - POST `/predict_gene_essentiality` (KO/scramble proxies + Evo2 signatures)
- Forge
  - POST `/generate_crispr_guides` (design + off-target proxy)
  - POST `/generate_inhibitors` (design + rank; orchestrate Boltz)
- Gauntlet
  - POST `/run_trials` (validate CRISPR/inhibitors; return efficacy metrics)
- Vercel Backend (proxy only)
  - Mirror the above with `/api/*` routes that forward to Modal

#### Phase 3 (Advanced)
- Evo Service
  - POST `/design_accessibility` (inference-time control for chromatin accessibility)
  - POST `/generate_sequence` (long-form generation with safety guardrails)
- Boltz Service
  - POST `/v1/predict_structure`
  - POST `/v1/binding_affinity`

### Request/Response schemas (examples)

#### Evo2 delta scoring
Request
```json
{
  "ref_sequence": "ACGT...",
  "alt_sequence": "ACGt..."
}
```
Response
```json
{
  "ref_likelihood": -1234.56,
  "alt_likelihood": -1248.90,
  "delta_likelihood_score": -14.34
}
```
Notes
- Negative delta → disruptive; magnitudes context-dependent
- Validate `[ACGTN]` only; cap sequences (e.g., <= 16,384)

#### Splicing/regulatory
Request
```json
{
  "assembly": "GRCh38",
  "chrom": "17",
  "pos": 43045736,
  "ref": "A",
  "alt": "G",
  "context": "splice" | "regulatory",
  "window": 8192
}
```
Response
```json
{
  "delta_likelihood_score": -3.21,
  "feature_scores": {
    "exon_boundary": 0.87,
    "tf_motif_disruption": 0.42
  },
  "explainability": {
    "saliency": [ ... ],
    "top_features": ["exon_start", "splice_acceptor"]
  }
}
```

### Implementation notes (Evo2 specifics)
- Windowing: default 8,192-nt centered on variant; pad/trim for indels to keep constant length
- Indels: represent minimal ref/alt; reflow flanks to preserve total length
- Batching: use small batches by token length; 7B for speed, 40B optional for higher fidelity
- Explainability: expose SAE-derived features where available (exon/intron boundaries, TF motifs, mutation severity)

### Deployment

#### Modal (services)
```bash
# One-time auth
venv/bin/modal token new

# Verify CLI
venv/bin/modal --version

# Deploy services (examples)
venv/bin/modal deploy src/services/evo_service/main.py
venv/bin/modal deploy src/services/oracle/main.py
venv/bin/modal deploy src/services/forge/main.py
venv/bin/modal deploy src/services/gauntlet/main.py
venv/bin/modal deploy src/services/boltz_service/main.py
venv/bin/modal deploy src/services/fusion_engine/main.py
```

#### Vercel (thin proxy backend)
- Layout: `oncology-coPilot/oncology-backend-minimal/api/index.py`
- Config: `vercel.json` routes all to `api/index.py`
- Requirements: `requirements-minimal.txt` → copy to `requirements.txt`
```bash
cd oncology-coPilot/oncology-backend-minimal
cp requirements-minimal.txt requirements.txt
npx -y vercel --yes --prod
```

### Environment variables
- In Vercel
  - `EVO_SERVICE_URL` → Modal evo-service base URL
  - Add others as services go live (`ORACLE_URL`, `FORGE_URL`, `GAUNTLET_URL`, etc.)
- In services
  - Use `modal.Secret.from_dotenv()` and `.env` for internal API keys

### Testing (smoke tests)
```bash
# Evo2 delta (Modal)
curl -X POST "$EVO_SERVICE_URL/score_delta" \
  -H 'content-type: application/json' \
  -d '{"ref_sequence":"ACGT", "alt_sequence":"AGGT"}'

# Evo2 delta (Vercel proxy)
curl -X POST "$VERCEL_BASE/api/evo/score_delta" \
  -H 'content-type: application/json' \
  -d '{"ref_sequence":"ACGT", "alt_sequence":"AGGT"}'
```

### Rollout plan
- Phase 1
  - Ship `/score_delta`, `/score_batch` in `evo-service` (Modal)
  - Add `/api/evo/score_delta`, `/api/evo/score_batch` to Vercel proxy
  - Wire Command Center to call Modal directly (keep Vercel for public/demo)
- Phase 2
  - Add splicing/regulatory scoring + embeddings
  - Implement Oracle/Forge/Gauntlet proxy routes
- Phase 3
  - Add generation endpoints with strong safety guardrails
  - Integrate Boltz structure/affinity simulation

### Command Center cleanup (scope reduction)
- Only orchestrate workflows
  - Transform inputs → call services → fuse outputs → return response
- Remove embedded compute/business logic that belongs in Evo2/Forge/Gauntlet
- Centralize downstream URLs as constants; use environment-driven config

### Success criteria
- Vercel backend stays <50MB and deploys in <60s
- Modal services expose stable, standardized URLs
- Command Center has <500 LOC of orchestration code, zero heavy deps
- End-to-end demo runs fully against Modal services with thin Vercel proxy

### Future considerations
- Add caching for repeated variant windows
- Add dataset-backed validations (ClinVar/SpliceVarDB) for CI
- Log schema versions in responses for robust clients 
 
 
 