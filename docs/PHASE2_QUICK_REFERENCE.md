# Phase 2 Data Hunt: Quick Reference

**For:** Parallel search agent  
**Mission:** Find serial SAE datasets  
**Time:** 1.5 hours

---

## 🎯 WHAT I NEED

**Find datasets with:**
- ✅ Paired samples (baseline + progression)
- ✅ Genomic data (WES mutations OR RNA-seq expression)
- ✅ Outcome data (TTP, PFS, response)
- ✅ n ≥ 20 patients with pairs
- ✅ Accessible data (public OR requestable)

---

## 🔍 WHERE TO SEARCH

### 1. PubMed (30 min)
**Queries:**
- `("ovarian cancer" OR "ovarian carcinoma") AND ("serial sequencing" OR "longitudinal sequencing" OR "paired samples")`
- `("ovarian cancer") AND ("ctDNA") AND ("tumor sequencing")`
- `("ovarian cancer") AND ("organoid") AND ("drug testing" OR "resistance")`

**Specific papers:**
- Goranova et al. (ovarian ctDNA)
- Parkinson et al. (ovarian serial plasma)
- Patch et al. (recurrent ovarian WGS)

### 2. GEO (20 min)
- Search: "ovarian cancer serial" OR "ovarian cancer longitudinal"
- Filter: Multiple samples per patient

### 3. cBioPortal (20 min)
- Filter: Ovarian cancer studies
- Check: Studies with "recurrent" samples

### 4. dbGaP (10 min)
- Search: Ovarian cancer studies
- Check: Serial samples

---

## 📋 FOR EACH DATASET, REPORT:

1. **Citation:** Author, Journal, Year, PMID
2. **Cohort:** Sample size, timepoints, treatment
3. **Data:** Types (WES/RNA-seq), availability (public/restricted)
4. **SAE Compatible?** YES/NO + reasoning
5. **Outcomes?** TTP, PFS, response available?
6. **Feasibility:** HIGH/MEDIUM/LOW

---

## 📝 UPDATE THIS FILE

**File:** `docs/SERIAL_SAE_DATA_SOURCES.md`

**Add sections for each dataset found using this template:**

```markdown
### Dataset X: [Study Name]

**Citation:** Author et al., Journal Year (PMID: xxxxx)

**Cohort:**
- Cancer type: Ovarian cancer
- Sample size: n=50
- Serial timepoints: Baseline, progression (n=30 pairs)
- Treatment: PARP inhibitor

**Data Availability:**
- Status: Public / Restricted
- Repository: GEO / dbGaP / SRA
- Accession: GSE12345
- Data types: WES, RNA-seq

**SAE Compatibility:**
- Can compute pathway scores: YES / NO
- Reasoning: [Why yes/no?]
- Feasibility: HIGH / MEDIUM / LOW

**Outcome Data:**
- Available: YES / NO
- Outcomes: TTP, PFS, response

**Notes:** [Any additional info]
```

---

## ✅ SUCCESS = FIND 1-2 DATASETS WITH:
- n ≥ 20 paired samples
- WES OR RNA-seq data
- Outcome data (TTP/PFS/response)
- Public OR requestable access

---

**Full specification:** See `PHASE2_DATA_HUNT_SPECIFICATION.md`  
**Agent instructions:** See `PHASE2_AGENT_INSTRUCTIONS.md`
