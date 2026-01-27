# Phase 2 Data Hunt: Agent Instructions

**For:** Parallel search agent  
**Mission:** Find serial SAE datasets (baseline + progression pairs)  
**Timeline:** 1.5 hours  
**Deliverable:** Updated `docs/SERIAL_SAE_DATA_SOURCES.md`

---

## 🎯 YOUR MISSION

Find datasets with **paired samples** (baseline + progression) that allow serial SAE validation.

**Minimum viable:** n ≥ 20 patients with paired samples

---

## 📋 EXACT SEARCH STRATEGY

### Step 1: PubMed Search (30 minutes)

**Query 1: Ovarian Serial Sequencing**
```
("ovarian cancer" OR "ovarian carcinoma") 
AND 
("serial sequencing" OR "longitudinal sequencing" OR "paired samples" OR "primary recurrent")
AND 
("whole exome" OR "WES" OR "RNA-seq")
```

**Query 2: Ovarian ctDNA + Tumor**
```
("ovarian cancer" OR "ovarian carcinoma")
AND
("ctDNA" OR "circulating tumor DNA")
AND
("tumor sequencing" OR "tissue sequencing")
```

**Query 3: Ovarian Organoid Serial**
```
("ovarian cancer" OR "ovarian carcinoma")
AND
("organoid" OR "patient-derived organoid")
AND
("drug testing" OR "resistance" OR "serial")
```

**For each paper found, extract:**
- Authors, journal, year, PMID
- Sample size, timepoints, treatment
- Data types (WES, RNA-seq, etc.)
- Data availability (public/restricted)
- SAE compatibility assessment

### Step 2: Specific Studies Check (20 minutes)

**Check these specific papers:**

1. **Goranova et al.**
   - Search: `"Goranova" AND "ovarian" AND "ctDNA"`
   - Question: Does it have paired tumor samples (not just plasma)?

2. **Parkinson et al.**
   - Search: `"Parkinson" AND "ovarian" AND "serial" AND "plasma"`
   - Question: Does it have tumor sequencing data?

3. **Patch et al.**
   - Search: `"Patch" AND "recurrent ovarian" AND "whole genome"`
   - Question: Does it have primary + recurrent pairs?

### Step 3: Repository Search (30 minutes)

**GEO (Gene Expression Omnibus):**
- Search: "ovarian cancer serial" OR "ovarian cancer longitudinal"
- Filter: Studies with multiple samples per patient
- Check: Expression data + clinical outcomes

**cBioPortal:**
- Filter: Ovarian cancer studies
- Check: Studies with "recurrent" or "progression" samples
- Question: Which have paired primary+recurrent?

**dbGaP:**
- Search: Ovarian cancer studies
- Check: Which have serial samples?
- Note: May require application

### Step 4: Document Findings (10 minutes)

**For each dataset found, create entry:**

```markdown
### Dataset X: [Study Name]

**Citation:** Author et al., Journal Year (PMID: xxxxx)

**Cohort:**
- Cancer type: Ovarian cancer
- Sample size: n=50
- Serial timepoints: Baseline, progression (n=30 pairs)
- Treatment: PARP inhibitor

**Data Availability:**
- Status: Public / Restricted / Not available
- Repository: GEO / dbGaP / SRA / Other
- Accession: GSE12345 or phs000123
- Data types: WES, RNA-seq, WGS, targeted_panel

**SAE Compatibility:**
- Can compute pathway scores: YES / NO
- Reasoning: Why yes/no?
- Required data: What data needed?
- Missing data: What's missing?
- Feasibility: HIGH / MEDIUM / LOW

**Outcome Data:**
- Available: YES / NO
- Outcomes: OS, PFS, TTP, response, resistance
- Time to progression: Available?

**Notes:** Any additional relevant information
```

---

## 🔍 WHAT TO LOOK FOR

### ✅ GOOD DATASETS (Prioritize These)

- **Has paired samples:** Baseline + progression (or mid-treatment)
- **Has genomic data:** WES mutations OR RNA-seq expression
- **Has outcomes:** TTP, PFS, response status, resistance labels
- **Is accessible:** Public download OR requestable
- **Sample size:** n ≥ 20 patients with pairs

### ❌ SKIP THESE (Don't Waste Time)

- Only plasma ctDNA (no tumor sequencing)
- Only single timepoint
- <10 paired samples (too small)
- No outcome data
- Restricted access with no contact info
- Only CNV data (no mutations/expression)

---

## 📊 REPORTING FORMAT

**Update:** `docs/SERIAL_SAE_DATA_SOURCES.md`

**Add sections:**
1. One section per dataset found (use template above)
2. Summary table at end
3. Priority ranking (Tier 1 = best, Tier 4 = future work)

**Summary Table Format:**
| Study | n (Paired) | Data Types | SAE Compatible | Access | Feasibility |
|-------|------------|------------|----------------|--------|-------------|
| Goranova 2017 | ? | ? | ? | ? | ? |
| ... | ... | ... | ... | ... | ... |

---

## ✅ SUCCESS CRITERIA

**Minimum Success:**
- Find 1-2 datasets with n ≥ 20 paired samples
- Data accessible
- SAE computable

**Ideal Success:**
- Find 3-5 datasets with n ≥ 20 paired samples
- Mix of ovarian + other cancers
- All publicly accessible

---

## 🚨 CRITICAL INFORMATION NEEDED

For each dataset, I need to know:

1. **Can we compute SAE?** (mutations OR expression available?)
2. **How many paired samples?** (n ≥ 20?)
3. **Is data accessible?** (public OR requestable?)
4. **What outcomes?** (TTP, PFS, response?)
5. **What's the feasibility?** (HIGH/MEDIUM/LOW)

**If you find a dataset but can't answer these questions, note what's missing.**

---

**Status:** ⏳ **READY FOR SEARCH**  
**Report back:** Within 1.5 hours with findings in specified format
