# Serial SAE Analysis Scripts

**Mission:** Compute serial SAE pathway scores and validate resistance prediction

---

## Scripts

### 1. `download_cbioportal_paired.py` ✅
**Purpose:** Extract paired primary+recurrent samples from cBioPortal

**Usage:**
```bash
python3 scripts/serial_sae/download_cbioportal_paired.py
```

**Output:**
- `data/serial_sae/cbioportal_paired/paired_patients.json`
- `data/serial_sae/cbioportal_paired/summary.json`

**Status:** ✅ Complete - Found 40 paired patients from MSK-SPECTRUM

---

### 2. `download_molecular_data_mcp.py` ⚠️
**Purpose:** Download mutations and expression using cBioPortal MCP API client

**Usage:**
```bash
python3 scripts/serial_sae/download_molecular_data_mcp.py
```

**Output:**
- `data/serial_sae/cbioportal_paired/msk_spectrum_tme_2022_mutations.csv`
- `data/serial_sae/cbioportal_paired/msk_spectrum_tme_2022_expression.csv`

**Status:** ⚠️ Working but returns 0 mutations - may need sample ID verification or manual download

**Note:** Uses cBioPortal MCP API client from `tools/cbioportal-mcp/`

---

### 3. `download_molecular_data.py` (Legacy)
**Purpose:** Original script using direct REST API (had 404 issues)

**Status:** ⚠️ Replaced by `download_molecular_data_mcp.py`

---

### 4. `process_manual_downloads.py` ✅
**Purpose:** Process manually downloaded TSV files from cBioPortal web interface

**Usage:**
```bash
# After manual download from cBioPortal web interface
python3 scripts/serial_sae/process_manual_downloads.py
```

**Input:** TSV files saved to `data/serial_sae/cbioportal_paired/`
- `msk_spectrum_tme_2022_mutations.tsv`
- `msk_spectrum_tme_2022_expression.tsv`

**Output:**
- `data/serial_sae/cbioportal_paired/processed_mutations.csv`
- `data/serial_sae/cbioportal_paired/processed_expression.csv`

**Status:** ✅ Ready - Use if MCP script doesn't work

---

### 5. `compute_serial_sae.py` ⏳
**Purpose:** Compute SAE pathway scores at baseline and progression

**Usage:**
```bash
python3 scripts/serial_sae/compute_serial_sae.py
```

**Output:**
- `data/serial_sae/results/serial_sae_computation.json`

**Status:** ⏳ Waiting for mutation/expression data

---

### 6. `download_gse165897.py` ✅
**Purpose:** Download GSE165897 scRNA-seq dataset (11 paired patients)

**Usage:**
```bash
python3 scripts/serial_sae/download_gse165897.py
```

**Status:** ✅ Created - Can proceed in parallel with cBioPortal data

---

### 7. `download_tcga_ov_gdc.py` ✅
**Purpose:** Query GDC API and create manifest for TCGA-OV files

**Usage:**
```bash
python3 scripts/serial_sae/download_tcga_ov_gdc.py
```

**Output:** Manifest file with file IDs

**Status:** ✅ Complete - Manifest created (434 RNA-seq files)

---

### 8. `download_tcga_ov_api.py` ✅ (RECOMMENDED)
**Purpose:** Download TCGA-OV files directly via GDC API (no gdc-client needed)

**Usage:**
```bash
# Download all files
python3 scripts/serial_sae/download_tcga_ov_api.py

# Or limit to first N files for testing
python3 scripts/serial_sae/download_tcga_ov_api.py 10
```

**Output:** Downloaded files in `data/serial_sae/tcga_ov/downloads/`

**Status:** ✅ Ready - Direct API download, no external tools required

---

### 9. `process_tcga_ov_rnaseq.py` ✅
**Purpose:** Process downloaded RNA-seq files into expression matrix

**Usage:**
```bash
python3 scripts/serial_sae/process_tcga_ov_rnaseq.py
```

**Output:** 
- `data/serial_sae/tcga_ov/processed/tcga_ov_expression_matrix.csv`
- `data/serial_sae/tcga_ov/processed/tcga_ov_sample_metadata.csv`

**Status:** ✅ Ready - Processes TSV files into combined expression matrix

---

**Note:** TCGA-OV has mostly primary tumors, not paired samples. Use for baseline pathway score validation.

---

## Data Sources

### Immediate (Open Access)
1. **cBioPortal MSK-SPECTRUM** - 40 paired patients ✅ (identified)
2. **GSE165897** - 11 paired patients (scRNA-seq) ✅ (script ready)

### Pending
3. **BriTROC-1** - 276 paired patients (EGA access)
4. **Williams et al.** - 18 patients (data embargoed)

---

## Next Steps

1. ✅ Extract paired patient IDs (DONE - 40 patients)
2. ⚠️ Download mutations for paired samples (MCP script working but 0 results - may need manual download)
3. ⏳ Download expression data for paired samples
4. ⏳ Compute SAE pathway scores
5. ⏳ Calculate pathway kinetics (ΔSAE)
6. ⏳ Correlate with outcomes

---

## Troubleshooting

### MCP Script Returns 0 Mutations
- **Check:** Verify sample IDs match cBioPortal database
- **Try:** Manual download via web interface (see `MANUAL_DOWNLOAD_INSTRUCTIONS.md`)
- **Alternative:** Use `process_manual_downloads.py` after manual download

### API Errors
- **Solution:** Use cBioPortal MCP API client (already integrated)
- **Fallback:** Manual download via web interface

---

**Status:** ⏳ **IN PROGRESS** - Scripts created, data download in progress
