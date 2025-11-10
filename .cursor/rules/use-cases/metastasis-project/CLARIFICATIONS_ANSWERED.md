# Clarifications Answered - Metastasis Interception Dominance Strategy

**Date:** October 13, 2025  
**Purpose:** Answer all technical questions from peer review to enable Day 1 execution

---

## Q1: Where is the canonical `metastasis_rules.json` for step labels?

**Answer:**
- **Path:** `oncology-coPilot/oncology-backend-minimal/api/rules/metastasis_rules.json`
- **Version:** v1.0.0 (will pin commit hash on Day 1)
- **Schema:**
```json
{
  "steps": {
    "local_invasion": {
      "genes": ["MMP2", "MMP9", "SNAI1", "TWIST1", ...],
      "pmids": ["12345678", "23456789", ...]
    },
    "intravasation": { ... },
    ...
  }
}
```
- **Expansion Plan:** Add 10 genes from ongoing trials (VEGFR2, AXL, TGFβR1, etc.) → 24 total
- **Provenance:** Commit hash + version in all figure legends

---

## Q2: What sequence context length for Enformer inference?

**Answer:**
- **Context:** ±32kb (64kb total) centered on variant position
- **Rationale:** Enformer expects 131kb input; we use 64kb centered window for speed/memory
- **Strand handling:** Both strands averaged (forward + reverse complement)
- **Tracks aggregated:** DNase-seq, CAGE, ATAC-seq → scalar accessibility [0,1]
- **Method:** Mean signal across tracks, normalized to [0,1]

---

## Q3: Which AF3 runner parameters (recycles, ensembles, templates)?

**Answer:**
- **Recycles:** 3 (balance quality vs speed)
- **Model:** model_1_multimer_v3 (single model, no ensembling)
- **Templates:** OFF (de novo prediction, no PDB templates)
- **Seed:** 42 (reproducibility)
- **Complex:** Multimer (guide RNA 20nt + target DNA 23nt + PAM 3nt)
- **Rationale:** Standard AlphaFold3 params for nucleic acid complexes

---

## Q4: Artifact storage target (S3/GCS bucket, path, retention)?

**Answer:**
- **Storage:** AWS S3 `s3://crispro-structures/`
- **Path structure:** `s3://crispro-structures/{job_id}/`
  - `{job_id}/structure.pdb` - PDB file
  - `{job_id}/metrics.json` - pLDDT, PAE, validation results
  - `{job_id}/provenance.json` - model, params, timestamps
- **Retention:** 1 year in S3, permanent mirror on Zenodo
- **Access:** Public read after publication, IAM write-only for service
- **Zenodo export:** Batch upload all 40 structures as supplementary dataset

---

## Q5: Fusion default - confirm "off" globally, auto-on for hotspots?

**Answer:**
- **Default:** OFF globally (determinism, no AlphaMissense latency)
- **Auto-enable:** Canonical hotspots only (BRAF V600E, KRAS G12D/G13C/Q61H, TP53 R175H/R248W/R273H)
- **Banner:** `✓ Fusion active (AlphaMissense)` when enabled
- **Coverage check:** Display badge when AM has GRCh38 missense coverage
- **Provenance:** Always log `fusion_enabled: true/false` in response metadata
- **User override:** Allow manual toggle in UI (profile selection)

---

## Q6: FE changes - which components consume structural metrics?

**Answer:**

### Components to Update:
1. **`StructurePredictionResults.jsx`**
   - **Props:** `{job_id, status, plddt_mean, pae_interface, pdb_url, pass}`
   - **Display:** Status badge, pLDDT bar chart, PAE heatmap, 3D viewer link
   
2. **`StructureViewer3D.jsx`**
   - **Props:** `{pdb_url, plddt_per_residue, highlights}`
   - **Display:** Mol* 3D viewer, pLDDT coloring, interaction highlights

3. **`MetastasisInterceptionPanel.jsx`**
   - **New section:** "Structural Validation" card
   - **Display:** Pass/fail badge, pLDDT mean, link to viewer

### API Schema:
```typescript
interface StructureJob {
  id: string;
  status: 'queued' | 'running' | 'succeeded' | 'failed' | 'retry';
  plddt_mean: number;        // Overall confidence
  pae_interface: number;     // Interface confidence
  clash_count: number;       // Stereochemistry
  molprobity_score: number;  // Overall quality
  pass: boolean;             // pLDDT≥70 AND PAE≤10
  pdb_url: string;           // S3 link
  created_at: string;
  completed_at: string;
  provenance: {
    model: string;
    recycles: number;
    seed: number;
    commit: string;
  }
}
```

### Polling Pattern:
```javascript
const pollJob = async (jobId) => {
  const response = await fetch(`/api/structure/${jobId}`);
  const job = await response.json();
  if (job.status === 'succeeded' || job.status === 'failed') {
    return job;
  }
  await new Promise(r => setTimeout(r, 5000)); // 5s
  return pollJob(jobId);
};
```

---

## Implementation Priority

### Immediate (Day 1):
1. ✅ Lock `metastasis_rules.json` path and version
2. ✅ Expand gene set (14 → 24)
3. ✅ Start per-step validation script

### Day 3:
4. ✅ Deploy Enformer with ±32kb context
5. ✅ Integrate backend with fallback

### Week 2:
6. ✅ Deploy AF3 with specified params
7. ✅ Wire FE components
8. ✅ S3 storage + Zenodo export

---

**Status:** ✅ ALL CLARIFICATIONS ANSWERED - READY TO EXECUTE  
**Next Step:** Start Day 1 per-step label derivation script  
**Command:** `./scripts/metastasis/derive_step_labels.py --expand-to 24 --output publication/data/step_labels_v1.0.0.json`

