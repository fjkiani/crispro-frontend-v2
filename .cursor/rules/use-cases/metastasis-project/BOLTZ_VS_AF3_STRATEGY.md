# ⚔️ BOLTZ VS ALPHAFOLD3 STRATEGIC ASSESSMENT

**Date:** October 13, 2025  
**Purpose:** Determine optimal structural validation strategy for future publications

---

## 🎯 THE CORE QUESTIONS

1. **Will Boltz deliver what we need for future publications?**
2. **Should we fix Boltz or pivot to AF3?**
3. **What does pLDDT 67.09 actually mean for our science?**

---

## 📊 CURRENT BOLTZ PERFORMANCE

### What We Validated (Fast Mode, msa='empty')
- **Runtime:** 16 seconds/protein
- **pLDDT Mean:** 67.09
- **PTM Score:** 0.43
- **Fraction Disordered:** 1.00
- **Mode:** Single-sequence (no MSA)

### What This Means Scientifically

**pLDDT Interpretation:**
- **90-100:** Very high confidence (publication-grade)
- **70-90:** High confidence (generally reliable)
- **50-70:** Low confidence (OK for screening, **NOT** for atomic detail)
- **<50:** Very low confidence (likely wrong)

**Our 67.09 = "Borderline acceptable for relative ranking only"**

**PTM 0.43 = "Poor predicted TM-score"** (want >0.5 for confidence)

**Fraction Disordered 1.00 = "Fully disordered"** (red flag - suggests protein may not fold properly!)

---

## ⚠️ **CRITICAL REALITY CHECK**

### What Boltz Fast-Mode CAN Do:
✅ Rapid screening (16s vs 60+ min)  
✅ Detect gross failures ("wet noodles")  
✅ Relative ranking (A better than B)  
✅ Infrastructure proof-of-concept  

### What Boltz Fast-Mode CANNOT Do:
❌ **Publication-grade structural predictions** (pLDDT too low)  
❌ **Atomic-detail analysis** (need pLDDT >70)  
❌ **Protein-protein interactions** (need full MSA for accuracy)  
❌ **gRNA:DNA complexes** (Boltz doesn't natively support nucleic acids well)  
❌ **Confident YES/NO on foldability** (fraction_disordered=1.00 is concerning)  

---

## 🔬 WHAT FUTURE PUBLICATIONS REQUIRE

### Week 1 Paper (Metastasis Interception) - **SUFFICIENT**
- ✅ Framework validation (Target Lock, AUROC 0.976)
- ✅ Guide design methodology
- ✅ Infrastructure deployment proof
- ⚠️ Structural validation = "nice to have" (current: adequate for RUO)

### Future Paper #2 (Wet Lab Validation) - **INSUFFICIENT**
- ❌ Need **high-confidence structures** (pLDDT >70) for experimental validation
- ❌ Need **atomic-detail predictions** for mutagenesis studies
- ❌ Need **protein-protein docking** for interaction studies
- ❌ Need **gRNA:Cas9:DNA complexes** for guide optimization

### Future Paper #3 (Clinical Translation) - **CRITICALLY INSUFFICIENT**
- ❌ Need **publication-grade structures** for regulatory submissions
- ❌ Need **full MSA mode** for maximum accuracy
- ❌ Need **validated predictions** against crystal structures
- ❌ Need **reproducible pipelines** with known accuracy bounds

---

## 🎯 BOLTZ FULL-MODE VS AF3 COMPARISON

| Feature | Boltz Fast (Current) | Boltz Full MSA | AlphaFold3 |
|---------|---------------------|----------------|------------|
| **Runtime** | 16s | 60+ min | ~5-10 min |
| **pLDDT** | ~67 | 70-90 | 70-95 |
| **MSA Required** | No | Yes | Yes (or custom) |
| **Protein-Protein** | Supported | Supported | **Best-in-class** |
| **Nucleic Acids** | Limited | Limited | **Native support** |
| **gRNA:DNA** | Poor | Poor | **Designed for this** |
| **Accuracy** | Screening-grade | Good | **State-of-the-art** |
| **Deployment** | ✅ Done | ⚠️ MSA timeout | ❌ Weights needed |
| **Publication-Grade** | ❌ No | ✅ Yes | ✅ Yes |

---

## ⚔️ **STRATEGIC RECOMMENDATION**

### **SHORT TERM (Week 1 Paper):**
✅ **Keep Boltz fast-mode as-is**
- Sufficient for RUO structural validation
- Demonstrates infrastructure capability
- Honest pLDDT reporting (67.09) with limitations documented
- **Verdict:** Ship it with clear RUO disclaimers

### **MEDIUM TERM (Weeks 2-4):**
🔧 **Fix Boltz full-MSA mode OR deploy ESMFold as backup**
- **Option A:** Fix Boltz MSA generation (use AF3's custom MSA approach)
  - Pro: Keep existing infrastructure
  - Con: Still won't handle gRNA:DNA complexes well
  - Timeline: 1-2 weeks
  
- **Option B:** Deploy ESMFold (no MSA required, faster than Boltz full)
  - Pro: 1-2 min/protein, pLDDT 70-80, no MSA headaches
  - Con: Monomer-only (no protein-protein or nucleic acids)
  - Timeline: 2-4 hours
  
- **Option C:** Request AF3 weights, build AF3 pipeline
  - Pro: Best-in-class, handles gRNA:DNA natively
  - Con: Weights approval required (~weeks), complex setup
  - Timeline: 2-4 weeks after weight approval

### **LONG TERM (Future Publications):**
🎯 **MUST deploy AlphaFold3 for publication-grade work**
- **Why:** gRNA:DNA complexes, protein-protein interactions, maximum accuracy
- **When:** Before wet lab validation paper (Paper #2)
- **Effort:** 2-4 weeks (weights + integration + validation)
- **Payoff:** State-of-the-art structural predictions for high-impact journals

---

## 📈 **VALIDATION REQUIREMENTS BY USE CASE**

### Use Case 1: Target Protein Foldability (Current)
- **Need:** pLDDT >50 (screening), >70 (confident)
- **Boltz Fast:** pLDDT 67 = ⚠️ Borderline
- **Boltz Full:** pLDDT 75-85 = ✅ Good
- **ESMFold:** pLDDT 70-80 = ✅ Good
- **AF3:** pLDDT 80-95 = ✅ Excellent

### Use Case 2: Protein-Protein Interactions
- **Need:** Interface PAE <10, high pLDDT at interface
- **Boltz Fast:** PTM 0.43 = ❌ Poor
- **Boltz Full:** PTM 0.6-0.8 = ✅ Good
- **ESMFold:** N/A = ❌ Monomer-only
- **AF3:** iPAE <5, PTM >0.8 = ✅ Excellent

### Use Case 3: gRNA:Cas9:DNA Complexes
- **Need:** Nucleic acid modeling, interface prediction
- **Boltz Fast:** N/A = ❌ Not designed for this
- **Boltz Full:** Limited = ⚠️ Poor nucleic acid support
- **ESMFold:** N/A = ❌ Monomer-only
- **AF3:** Native support = ✅ **ONLY VIABLE OPTION**

### Use Case 4: HDR Template Validation
- **Need:** DNA structure prediction, Cas9:DNA binding
- **Boltz:** ❌ Poor
- **ESMFold:** ❌ No DNA
- **AF3:** ✅ **REQUIRED**

---

## 🎯 **IMMEDIATE ACTION PLAN**

### Week 1 (Current Paper) - **NO CHANGES**
- ✅ Submit with Boltz fast-mode (pLDDT 67.09)
- ✅ Clear RUO disclaimers
- ✅ Document limitations honestly
- **Status:** Complete

### Week 2-3 (Enhancement)
**Option A (Quick Win):** Deploy ESMFold
- **Why:** 1-2 min/protein, pLDDT 70-80, no MSA headaches
- **Effort:** 2-4 hours
- **Use Cases:** Single-protein foldability checks
- **Limitation:** No protein-protein, no DNA

**Option B (Boltz Fix):** Implement custom MSA approach from AF3
- **Why:** Keep existing infrastructure, improve accuracy
- **Effort:** 1-2 weeks
- **Use Cases:** Protein-protein interactions (with caveats)
- **Limitation:** Still poor for gRNA:DNA

### Week 4+ (Production-Grade)
**MUST DO:** Request and deploy AlphaFold3
- **Why:** Future papers (wet lab, clinical) require it
- **Effort:** Request weights NOW (approval takes weeks), then 2-4 weeks integration
- **Use Cases:** Everything (gRNA:DNA, protein-protein, HDR templates)
- **Blocker:** Weights approval

---

## 📊 **BOLTZ VALIDATION GAPS**

### What We Know:
✅ Service deploys successfully  
✅ Runs in 16 seconds (fast mode)  
✅ Returns structured output (pLDDT, PTM, disorder)  

### What We DON'T Know:
❌ **Accuracy against ground truth** (no crystal structures compared)  
❌ **Reliability across protein types** (only 1 smoke test)  
❌ **Performance on actual guide designs** (haven't validated real sequences)  
❌ **Correlation with wet lab outcomes** (no experimental validation)  

### What We MUST Validate:
1. **Benchmark against known structures** (5-10 proteins with crystal structures)
2. **Test on real guide sequences** (10-20 from our metastasis dataset)
3. **Compare pLDDT to experimental foldability** (if data available)
4. **Assess disorder predictions** (compare to IDPred/IUPred)

---

## ⚔️ **COMMANDER'S DECISION MATRIX**

### **For Week 1 Paper:**
**Decision:** ✅ **SHIP WITH BOLTZ FAST-MODE**
- Rationale: Adequate for RUO, demonstrates infrastructure, honest limitations
- Risk: Low (clear disclaimers, not core to validation)

### **For Future Papers:**
**Decision:** 🎯 **DEPLOY AF3 (Critical Path)**
- Rationale: gRNA:DNA complexes are NON-NEGOTIABLE for future work
- Timeline: Request weights NOW, deploy in 4-6 weeks
- Backup: ESMFold for simple protein foldability (2-4 hours to deploy)

### **For Boltz Full-MSA:**
**Decision:** ⏸️ **DEFER (Nice-to-Have)**
- Rationale: Effort better spent on AF3 deployment
- Use Case: Backup for protein-protein if AF3 weights delayed
- Priority: Low (ESMFold covers most protein-only use cases)

---

## 🚀 **RECOMMENDED EXECUTION ORDER**

1. **NOW:** Submit Week 1 with Boltz fast-mode (complete)
2. **Week 2:** Deploy ESMFold (4 hours) as immediate protein validation backup
3. **Week 2:** Request AF3 weights (start approval process)
4. **Week 3-4:** Benchmark Boltz fast-mode against known structures (validation)
5. **Week 4-6:** Deploy AF3 (upon weight approval)
6. **Week 6+:** Full validation of AF3 for Paper #2 (wet lab)

---

**VERDICT:** 
- **Boltz Fast-Mode:** ✅ Sufficient for Week 1, NOT for future publications
- **Boltz Full-MSA:** ⚠️ Defer (ESMFold faster, AF3 better)
- **AlphaFold3:** 🎯 **CRITICAL for future work - start approval NOW**

---

**Next Action:** Request AF3 weights from Google DeepMind while we finish Week 1 submission


