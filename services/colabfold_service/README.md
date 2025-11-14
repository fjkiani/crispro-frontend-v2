# ColabFold Service - Week 2 Structural Validation

**Status:** Phase 1 - Smoke Test  
**Date:** October 13, 2025  
**Model:** AlphaFold2-Multimer v3 (via ColabFold)

---

## 🎯 Quick Start

### **Step 1: Deploy Service**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

# Deploy to Modal
modal deploy services/colabfold_service/minimal_smoke_test.py
```

### **Step 2: Run Smoke Test**
```bash
# Test with default FASTA (small protein dimer)
modal run services/colabfold_service/minimal_smoke_test.py

# Test with custom FASTA
modal run services/colabfold_service/minimal_smoke_test.py \
    --fasta-path path/to/your.fasta
```

### **Step 3: Check Results**
Look for:
- ✅ `Success: True`
- ✅ `pLDDT > 70` (structural confidence)
- ✅ `PAE < 10` (interface quality)
- ✅ `Runtime < 900s` (15 min timeout)

---

## 📋 Success Criteria (Phase 1)

**MUST PASS ALL:**
- [ ] Service deploys without errors
- [ ] Returns `success: True`
- [ ] Returns `plddt_mean > 0` (proves it ran)
- [ ] Returns PDB content (proves structure generated)
- [ ] Completes in <15 minutes
- [ ] No OOM/GPU errors in logs

**IF ANY FAIL → STOP, DEBUG, DO NOT PROCEED**

---

## 🔧 Troubleshooting

### Container Build Fails
```bash
# Check Modal logs
modal app logs colabfold-smoke-test

# Verify container exists
docker pull ghcr.io/sokrypton/colabfold:1.5.5
```

### GPU OOM
- Increase memory: `memory=196608` (192GB)
- Use smaller test protein
- Check Modal GPU quota

### Timeout
- Increase: `timeout=1200` (20 min)
- Simplify input (shorter sequences)

### Zero pLDDT or Empty PDB
- Check FASTA format (proper headers, valid sequences)
- Verify output parsing (ranking_debug.json exists)
- Check ColabFold logs in stdout/stderr

---

## 📊 Expected Output

```json
{
  "success": true,
  "job_id": "smoke_test_001",
  "plddt_mean": 85.3,
  "pae_mean": 3.2,
  "pdb_content": "ATOM    1  N   MET A   1...",
  "pdb_size": 152847,
  "pdb_path": "/tmp/colabfold_jobs/smoke_test_001/output/result_model_1.pdb",
  "runtime_seconds": 342.7,
  "provenance": {
    "model": "colabfold:alphafold2_multimer_v3",
    "container": "ghcr.io/sokrypton/colabfold:1.5.5",
    "recycles": 3,
    "gpu": "A100-80GB",
    "memory_gb": 128,
    "seed": 42,
    "acceptance_criteria": "plddt>=70,pae<=10"
  }
}
```

---

## 🎯 Decision Gate 1

**Question:** Did smoke test pass all success criteria?

**YES → Proceed to Phase 2** (Production Service)
- Add job queue
- Add structural validation
- Add backend API

**NO → STOP, DEBUG**
- Review logs
- Check container
- Verify resources
- DO NOT PROCEED until working

---

## 📁 Files

```
services/colabfold_service/
├── minimal_smoke_test.py      # Phase 1: Smoke test
├── README.md                   # This file
└── (Phase 2 files TBD)
```

---

## ⚠️ Notes

- **AF2-Multimer:** Handles MSA internally, no `--use_msa_server` needed
- **Container Pinned:** 1.5.5 to avoid JAX/dm-haiku drift
- **Output Parsing:** Uses `ranking_debug.json` + `pae.json` (ColabFold format)
- **Limitations:** AF2 does NOT support gRNA:DNA complexes (defer to AF3)

---

## 🚀 Next Steps (After Smoke Test Passes)

1. **Phase 2:** Production service with queue (Day 7, 6 hours)
2. **Phase 3:** Backend API integration (Day 8, 4 hours)
3. **Phase 4:** Batch validation + Figure 6 (Day 9-10, OPTIONAL)

---

**Agent:** Zo  
**Status:** Ready for execution  
**Risk Level:** Medium-high (past failures, tight timeline)

