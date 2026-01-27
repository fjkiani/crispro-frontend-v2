# Phase 2 Data Hunt: Detailed Specification

**Date:** January 13, 2026  
**Mission:** Find serial SAE datasets (baseline + progression pairs)  
**For:** Parallel agent search  
**Timeline:** 1.5 hours

---

## 🎯 OBJECTIVE

Find datasets with **paired samples** (baseline + progression) that allow:
1. Compute SAE pathway scores at multiple timepoints
2. Correlate pathway changes (ΔSAE) with clinical outcomes
3. Validate serial monitoring hypothesis

**Minimum viable:** n ≥ 20 patients with paired samples

---

## 📋 DATA SOURCE REQUIREMENTS

### Required Information Per Dataset

For each dataset found, extract:

```json
{
  "study": {
    "authors": "First Author et al.",
    "journal": "Journal Name",
    "year": 2020,
    "pmid": "12345678",
    "doi": "10.xxxx/xxxx",
    "title": "Full paper title"
  },
  "cohort": {
    "cancer_type": "ovarian cancer / breast cancer / etc.",
    "sample_size": 50,
    "serial_timepoints": ["baseline", "progression", "mid-treatment?"],
    "timepoint_dates": "Are actual dates available?",
    "treatment": "What treatment was given? (e.g., PARP inhibitor, platinum)"
  },
  "data_availability": {
    "status": "public / restricted / not available",
    "repository": "GEO / dbGaP / SRA / other",
    "accession_id": "GSE12345 or phs000123",
    "access_required": "Open / dbGaP application / Contact authors",
    "data_types": ["WES", "RNA-seq", "WGS", "targeted_panel"]
  },
  "sae_compatibility": {
    "can_compute_pathway_scores": "yes/no",
    "reasoning": "Why yes/no?",
    "required_data": "What data needed for SAE?",
    "missing_data": "What's missing?",
    "feasibility": "high / medium / low"
  },
  "outcome_data": {
    "available": "yes/no",
    "outcomes": ["OS", "PFS", "TTP", "response", "resistance"],
    "time_to_progression": "available?",
    "resistance_labels": "available?"
  },
  "notes": "Any additional relevant information"
}
```

---

## 🔍 SEARCH STRATEGY

### Search 1: Ovarian Cancer Serial Sequencing

**PubMed Queries:**
1. `("ovarian cancer" OR "ovarian carcinoma") AND ("serial sequencing" OR "longitudinal sequencing" OR "paired samples" OR "primary recurrent")`
2. `("ovarian cancer" OR "ovarian carcinoma") AND ("serial ctDNA" OR "longitudinal ctDNA" OR "ctDNA monitoring")`
3. `("ovarian cancer" OR "ovarian carcinoma") AND ("recurrent tumor" OR "progression") AND ("whole exome" OR "WES" OR "RNA-seq")`

**Specific Studies to Check:**
1. **Goranova et al.**
   - Search: `"Goranova" AND "ovarian" AND "ctDNA"`
   - Check: Does it have paired tumor samples (not just plasma)?

2. **Parkinson et al.**
   - Search: `"Parkinson" AND "ovarian" AND "serial" AND "plasma"`
   - Check: Does it have tumor sequencing data?

3. **Patch et al.**
   - Search: `"Patch" AND "recurrent ovarian" AND "whole genome"`
   - Check: Does it have primary + recurrent pairs?

4. **TCGA Recurrent Analysis Papers**
   - Search: `"TCGA" AND "ovarian" AND "recurrent" AND "sequencing"`
   - Check: Any papers analyzing recurrent samples separately?

### Search 2: Breast Cancer Serial Sequencing

**Why:** Cross-cancer validation, may have more serial data

**PubMed Queries:**
1. `("breast cancer" OR "breast carcinoma") AND ("serial sequencing" OR "longitudinal sequencing" OR "paired samples")`
2. `("breast cancer" OR "breast carcinoma") AND ("primary recurrent" OR "baseline progression") AND ("WES" OR "RNA-seq")`

**Specific Studies:**
1. **METABRIC Serial Analysis**
   - Check: Does METABRIC have serial samples?

2. **TRACERx (Lung, but methodology)**
   - Check: Similar methodology for ovarian?

### Search 3: Organoid Serial Testing

**PubMed Queries:**
1. `("ovarian cancer" OR "ovarian carcinoma") AND "organoid" AND ("drug testing" OR "resistance" OR "serial")`
2. `("patient-derived organoid" OR "PDO") AND "ovarian" AND ("longitudinal" OR "serial" OR "time course")`

**What to Look For:**
- Organoid datasets with drug exposure over time
- Pathway changes during drug treatment
- Genomic data (mutations, expression) at multiple timepoints

### Search 4: Clinical Trial Datasets

**PubMed Queries:**
1. `("clinical trial" OR "phase 2" OR "phase 3") AND "ovarian" AND ("biopsy" OR "tissue") AND ("baseline" AND "progression")`
2. `("PARP inhibitor" OR "niraparib" OR "olaparib") AND "ovarian" AND ("serial" OR "paired" OR "biopsy")`

**What to Look For:**
- Trials that collected biopsies at multiple timepoints
- Trials with published genomic data
- Trials with resistance outcomes

### Search 5: Data Repositories

**GEO (Gene Expression Omnibus):**
- Search: `"ovarian cancer" AND ("serial" OR "longitudinal" OR "time course")`
- Filter: Studies with multiple samples per patient
- Check: Expression data + clinical outcomes

**dbGaP (Database of Genotypes and Phenotypes):**
- Search: Ovarian cancer studies
- Check: Which have serial samples?
- Note: May require dbGaP application

**SRA (Sequence Read Archive):**
- Search: Ovarian cancer + serial/longitudinal
- Check: WES/RNA-seq data availability

**cBioPortal:**
- Check: Studies with "recurrent" or "progression" samples
- Filter: Ovarian cancer studies
- Check: Which have paired primary+recurrent?

---

## 🔬 SAE COMPATIBILITY ASSESSMENT

### What SAE Needs

**Required Data:**
1. **Somatic Mutations** (for pathway scoring)
   - VCF files OR MAF files OR mutation calls
   - Genes: DDR (BRCA1/2, ATM, etc.), MAPK (KRAS, BRAF, etc.), PI3K (PIK3CA, PTEN, etc.)

2. **OR Expression Data** (alternative pathway scoring)
   - RNA-seq counts OR microarray expression
   - Can infer pathway activity from expression

3. **Clinical Outcomes** (for validation)
   - Time to progression (TTP)
   - Progression-free survival (PFS)
   - Response status (sensitive/resistant)
   - Overall survival (OS)

### Compatibility Checklist

For each dataset, assess:

- [ ] **Has somatic mutations?** (WES, targeted panel, or mutation calls)
- [ ] **Has expression data?** (RNA-seq or microarray)
- [ ] **Has paired samples?** (baseline + progression)
- [ ] **Has outcome data?** (TTP, PFS, response)
- [ ] **Is data accessible?** (public download or requestable)
- [ ] **Sample size adequate?** (n ≥ 20 for validation)

**Feasibility Scoring:**
- **High:** Has mutations + outcomes + paired samples + accessible
- **Medium:** Has expression + outcomes + paired samples + accessible
- **Low:** Missing critical data or not accessible

---

## 📊 REPORTING FORMAT

### For Each Dataset Found

Create entry in `docs/SERIAL_SAE_DATA_SOURCES.md`:

```markdown
### Dataset X: [Study Name]

**Citation:** Author et al., Journal Year (PMID: xxxxx)

**Cohort:**
- Cancer type: Ovarian cancer
- Sample size: n=50
- Serial timepoints: Baseline, progression (n=30 pairs)
- Treatment: PARP inhibitor (niraparib)

**Data Availability:**
- Status: Public
- Repository: GEO
- Accession: GSE12345
- Data types: WES, RNA-seq

**SAE Compatibility:**
- Can compute pathway scores: ✅ YES
- Reasoning: Has WES mutations + RNA-seq expression
- Required data: ✅ All available
- Feasibility: HIGH

**Outcome Data:**
- Available: ✅ YES
- Outcomes: TTP, PFS, response status
- Resistance labels: ✅ YES

**Notes:** 
- 30 patients with paired baseline+progression samples
- All have WES and RNA-seq data
- Publicly available via GEO
- Perfect for serial SAE validation
```

### Summary Table

Create summary table at end of document:

| Study | n (Paired) | Data Types | SAE Compatible | Access | Feasibility |
|-------|------------|------------|----------------|--------|-------------|
| Goranova 2017 | ? | ? | ? | ? | ? |
| Parkinson 2018 | ? | ? | ? | ? | ? |
| ... | ... | ... | ... | ... | ... |

---

## 🎯 PRIORITY RANKING

### Tier 1: Highest Priority
- Ovarian cancer + paired samples + WES/RNA-seq + outcomes + public access
- **Target:** n ≥ 20 patients

### Tier 2: Medium Priority
- Ovarian cancer + paired samples + expression only + outcomes + public access
- **Target:** n ≥ 20 patients

### Tier 3: Lower Priority
- Other cancers + paired samples + WES/RNA-seq + outcomes + public access
- **Use:** Cross-cancer validation

### Tier 4: Future Work
- Restricted access datasets (dbGaP, contact authors)
- **Use:** If Tier 1-3 insufficient

---

## 🔍 SPECIFIC SEARCH TERMS

### PubMed Advanced Search

**Query 1: Ovarian Serial Sequencing**
```
("ovarian cancer"[Title/Abstract] OR "ovarian carcinoma"[Title/Abstract]) 
AND 
("serial sequencing"[Title/Abstract] OR "longitudinal sequencing"[Title/Abstract] 
OR "paired samples"[Title/Abstract] OR "primary recurrent"[Title/Abstract])
AND 
("whole exome"[Title/Abstract] OR "WES"[Title/Abstract] OR "RNA-seq"[Title/Abstract])
```

**Query 2: Ovarian ctDNA + Tumor**
```
("ovarian cancer"[Title/Abstract] OR "ovarian carcinoma"[Title/Abstract])
AND
("ctDNA"[Title/Abstract] OR "circulating tumor DNA"[Title/Abstract])
AND
("tumor sequencing"[Title/Abstract] OR "tissue sequencing"[Title/Abstract])
```

**Query 3: Ovarian Organoid Serial**
```
("ovarian cancer"[Title/Abstract] OR "ovarian carcinoma"[Title/Abstract])
AND
("organoid"[Title/Abstract] OR "patient-derived organoid"[Title/Abstract])
AND
("drug testing"[Title/Abstract] OR "resistance"[Title/Abstract] OR "serial"[Title/Abstract])
```

### GEO Search

**Search Terms:**
- "ovarian cancer serial"
- "ovarian cancer longitudinal"
- "ovarian cancer time course"
- "ovarian cancer paired"

**Filters:**
- Expression profiling by array OR Expression profiling by high throughput sequencing
- Multiple samples per patient

### cBioPortal Search

**Study Filters:**
- Cancer Type: Ovarian Cancer
- Sample Type: Contains "recurrent" OR "progression"
- Data Types: Mutations, Expression

---

## 📝 DATA EXTRACTION TEMPLATE

For each paper found, extract:

1. **Paper Metadata**
   - Authors, journal, year, PMID, DOI, title

2. **Cohort Details**
   - Cancer type, sample size, timepoints, treatment

3. **Data Information**
   - What data types? (WES, RNA-seq, WGS, panel)
   - How many patients have paired samples?
   - What timepoints? (baseline, 3mo, 6mo, progression?)

4. **Access Information**
   - Public? Repository? Accession ID?
   - Restricted? How to access?

5. **SAE Assessment**
   - Can we compute pathway scores? Why/why not?
   - What's missing?
   - Feasibility score

6. **Outcome Information**
   - What outcomes available?
   - TTP? PFS? Response? Resistance?

---

## ✅ SUCCESS CRITERIA

**Minimum Success:**
- Find 1-2 datasets with n ≥ 20 paired samples
- Data accessible (public or requestable)
- SAE computable (mutations or expression)
- Outcomes available

**Ideal Success:**
- Find 3-5 datasets with n ≥ 20 paired samples
- Mix of ovarian + other cancers
- Mix of WES + RNA-seq
- All publicly accessible

**Bonus:**
- Find organoid datasets with serial drug testing
- Find clinical trial datasets with serial biopsies
- Find datasets with >50 paired samples

---

## 🚨 RED FLAGS (Skip These)

**Don't Waste Time On:**
- ❌ Studies with only plasma ctDNA (no tumor sequencing)
- ❌ Studies with only single timepoint
- ❌ Studies with <10 paired samples (too small)
- ❌ Studies with no outcome data
- ❌ Studies with restricted access AND no contact info
- ❌ Studies with only CNV data (no mutations/expression)

**Focus On:**
- ✅ Paired samples (baseline + progression)
- ✅ Genomic data (mutations OR expression)
- ✅ Outcome data (TTP, PFS, response)
- ✅ Accessible data (public OR requestable)

---

## 📋 DELIVERABLE

**File:** `docs/SERIAL_SAE_DATA_SOURCES.md` (update existing)

**Format:** 
- One section per dataset found
- Summary table at end
- Priority ranking
- Next steps for data access

**Timeline:** Report back within 1.5 hours with findings

---

**Status:** ⏳ **READY FOR PARALLEL SEARCH**  
**Agent Instructions:** Follow this specification exactly, report all findings in specified format
