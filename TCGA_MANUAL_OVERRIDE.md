# 🛑 Zeta Protocol: Manual Override Required

Alpha, the automatic data extraction systems (Python, Curl) were blocked by the target's defense systems (403 Forbidden). 
However, the **Weapon Logic** is built and ready to fire. We just need the ammunition.

## 🧱 The Blockade
- **Target**: TCGA-OV Expression Data (Transcriptional State).
- **Status**: Access Denied (Network/SSL restrictions).
- **Impact**: We cannot run the transcriptional validation automatically.

## 🔓 The Override (Your Action)
I need you to manually download **ONE file** and place it in the directory.

1. **Download**: [Link to cBioPortal TCGA-OV Dataset](https://cbioportal-datahub.s3.amazonaws.com/ov_tcga_pan_can_atlas_2018.tar.gz)
   - *Alternative*: If that link fails in browser, go to [cBioPortal Data Sets](https://www.cbioportal.org/datasets), search "Ovarian Serous Cystadenocarcinoma (TCGA, PanCancer Atlas)", and click the download icon.

2. **Place File**:
   - Extract the `data_mrna_seq_v2_rsem.txt` file.
   - Move it to: `oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/xena_mirror/expression.tsv` (or just tell me where you put it).

## 🚀 Activation Sequence
Once the file is there, run this single command to generate the proof:

```bash
# 1. Score the Patients (DDR & E2F)
python3 scripts/validation/score_tcga_hallmarks.py --input <PATH_TO_FILE> --output scores.csv

# 2. Correlate with Survival (The Truth)
python3 scripts/validation/correlate_tcga_outcomes.py --scores scores.csv --survival oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/xena_mirror/survival.tsv.gz
```

*Nyx status: Weapon primed. Awaiting ammo.*
