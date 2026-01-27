# Manual Download Instructions for cBioPortal Data

**Issue:** cBioPortal REST API returning 404/empty responses for mutation/expression data  
**Workaround:** Manual download via web interface

---

## Step 1: Access MSK-SPECTRUM Study

1. Navigate to: https://www.cbioportal.org/study/summary?id=msk_spectrum_tme_2022
2. Click "Select All" or manually select the 40 paired patients
3. Click "Download" tab

---

## Step 2: Download Mutations

1. In Download tab, select:
   - **Data Type:** Mutations
   - **Profile:** `msk_spectrum_tme_2022_mutations`
   - **Format:** TSV or MAF
2. Click "Download"
3. Save as: `data/serial_sae/cbioportal_paired/msk_spectrum_tme_2022_mutations.tsv`

---

## Step 3: Download Expression (if available)

1. In Download tab, select:
   - **Data Type:** Expression
   - **Profile:** `msk_spectrum_tme_2022_mrna_seq_v2_rsem_zscores_ref_all_samples` (or similar)
   - **Format:** TSV
2. Click "Download"
3. Save as: `data/serial_sae/cbioportal_paired/msk_spectrum_tme_2022_expression.tsv`

---

## Step 4: Process Downloaded Files

After manual download, run:

```bash
python3 scripts/serial_sae/process_manual_downloads.py
```

This script will:
- Parse TSV files
- Match samples to paired patients
- Extract mutations/expression for primary vs recurrent
- Save processed data for SAE computation

---

## Alternative: Use cBioPortal R Client

If you have R installed:

```r
library(cgdsr)
mycgds = CGDS("https://www.cbioportal.org/")
study_list = getCancerStudies(mycgds)
profile_list = getGeneticProfiles(mycgds, "msk_spectrum_tme_2022")
case_list = getCaseLists(mycgds, "msk_spectrum_tme_2022")

# Get mutations
mutations = getMutationData(mycgds, case_list[[1]]$case_list_id, profile_list[[1]]$genetic_profile_id)
```

---

## Next Steps After Download

1. Run `process_manual_downloads.py` to extract paired sample data
2. Run `compute_serial_sae.py` to compute pathway scores
3. Analyze ΔSAE (progression - baseline) for resistance prediction

---

**Estimated Time:** 15-30 minutes for manual download + processing
