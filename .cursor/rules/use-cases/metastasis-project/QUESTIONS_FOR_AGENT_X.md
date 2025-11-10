# Questions for Agent X - Boltz/AF3/ESMFold Strategy

**Context:** We have 3 structural validation options for Week 1 metastasis publication. All have tradeoffs.

---

## Option 1: Boltz (We have 2 implementations)

### Q1: MSA Dependency
**Both Boltz services call `tools/msa_client.run_mmseqs2()` which hits the remote MMseqs2 API.**
- This is the SAME bottleneck that's timing out ColabFold (60+ min).
- Can we modify Boltz to skip MSA entirely for a fast "monomer single-sequence" mode?
- Or does Boltz REQUIRE MSA like AF2-multimer does?

### Q2: Fast Mode Feasibility
**Looking at `src/services/boltz_service/main.py` lines 220-228:**
```python
# --- OPERATION: MSA RESTORATION ---
print("Fetching MSA for target sequence...")
msa_prefix = f"/tmp/msa_{job_id}"
msa_results = run_mmseqs2(protein_sequence, prefix=msa_prefix)
```

**Can we:**
- A) Comment out MSA generation and pass `None` or empty MSA to Boltz?
- B) Generate a "minimal mock MSA" (single sequence) locally in <1 second?
- C) Use Boltz's native MSA generation (does it have one)?

### Q3: Deployment Status
**Which Boltz service should we use?**
- `src/services/boltz/main.py` - Uses "self-docking" trick, older
- `src/services/boltz_service/main.py` - Newer, more complete, uses MSA client
- Are Boltz weights already downloaded? Do we need to run the downloader first?

---

## Option 2: AlphaFold3 (We just downloaded the code)

### Q4: AF3 Weights Status
**You mentioned requesting weights from Google DeepMind.**
- What's the typical approval timeline? (Days? Weeks?)
- Can we request NOW and continue in parallel?
- Link to request form?

### Q5: AF3 Fast Mode
**AF3 docs mention custom MSAs (input.md line 21-22).**
- Can we provide an EMPTY or single-sequence MSA to skip the remote search?
- Would this work for smoke tests (accept lower accuracy)?
- What's the expected runtime with `--run_data_pipeline=false` + minimal MSA?

### Q6: AF3 Without Weights
**Can we run AF3 data tests NOW to validate the environment?**
```bash
cd google-deepmind-alphafold3
python run_alphafold_data_test.py --output_dir /tmp/af3_out
```
- Does this work without model weights?
- Would this prove our AF3 integration is ready?

---

## Option 3: ESMFold (Fast fallback)

### Q7: ESMFold Deployment
**ESMFold is mentioned as 1-2 min/structure, no MSA needed.**
- Do we have ESMFold code/service anywhere in our codebase?
- If not, how hard is it to deploy? (Estimate: hours/days?)
- What's the output format? (pLDDT only? Or full structure?)

### Q8: ESMFold Limitations
**Doctrine says: "Does not model complexes/interfaces as well as AF2-multimer."**
- For Week 1, we're validating SINGLE guide RNAs (20-23 nt), not complexes.
- Is ESMFold good enough for this use case with RUO disclaimers?
- What's the acceptance threshold? (pLDDT ≥ 70 like AF2?)

---

## Strategic Decision Questions

### Q9: Week 1 Timeline Reality Check
**It's Day 6 of a 14-day timeline. We need to decide NOW:**
- A) **Submit Week 1 WITHOUT structural validation** (AUROC 0.976 is already excellent)
- B) **Deploy Boltz fast-mode** (if MSA can be skipped/mocked) - 4-6 hours
- C) **Deploy ESMFold** (if easy) - 2-4 hours
- D) **Wait for AF3 weights** - unknown timeline

**Your recommendation?**

### Q10: Boltz "Self-Docking" vs Standard
**`src/services/boltz/main.py` uses a "self-docking" trick (lines 129-133):**
```python
config_data = {
    "entities": [
        {"name": "self_target", "type": "protein", "sequence": req.protein_sequence},
        {"name": "self_candidate", "type": "protein", "sequence": req.protein_sequence}
    ]
}
```

**Does this bypass MSA generation? Or is it still slow?**

### Q11: Provenance & RUO
**If we use any structural method for Week 1:**
- What provenance fields MUST we track? (model, version, MSA status, etc.)
- What RUO language is required in Methods/Figures?
- Can we say "structural assessment pending optimization" if we skip it?

---

## My Current Recommendation (Needs Your Input)

**Based on timeline pressure and existing validation strength:**

1. **Submit Week 1 NOW** with:
   - ✅ Target-Lock validation (AUROC 0.976, n=38 genes, 8 steps)
   - ✅ Guide design validation (20 real designs, efficacy/safety/assassin scores)
   - 📝 Document structural validation as:
     ```
     "Structural validation infrastructure deployed (Boltz-2/ColabFold/AF3-ready).
      Full structural campaign (40 structures) deferred pending MSA optimization
      (current: 60+ min/structure via remote MMseqs2; target: <5 min via local
      MSA generation or single-sequence mode). Research Use Only."
     ```

2. **In Parallel (Post-submission):**
   - Request AF3 weights NOW
   - Deploy ESMFold for fast figures (if feasible)
   - OR implement Boltz fast-mode (if MSA can be skipped)

**Do you agree? Or should we push harder to get structures into Week 1?**

---

## Immediate Next Steps (Waiting for Your Answers)

Once you answer these questions:
- [ ] Update blocker report with chosen strategy
- [ ] Update METHODS_DRAFT with structural validation status
- [ ] Deploy chosen structural service (if any)
- [ ] Finalize Week 1 submission package

**Waiting for your strategic guidance, Agent X!** 🎯


