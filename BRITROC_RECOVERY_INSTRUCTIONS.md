# BriTROC-1 Data Recovery: The "R-Wall" Breaker

## 🚨 The Situation
The previous agent correctly identified that **BriTROC-1 Copy Number (CN) Signature** data is the "Gold Mine" for validating validation. 
However, they failed to mention that this data is locked inside **R Binary (`.rds`)** files and **does not exist** as CSV/TSV in the repository. 

Our current environment **does not have R installed**, and Python tools (`pyreadr`) failed to build due to missing system libraries (`lzma`).

## 🛠️ The Solution
You must extraction the data using an environment that has R installed.

### Step 1: Install R (if not present)
`brew install r` (on Mac)

### Step 2: Run this Extraction Script
Save the following as `extract_rds.R` in the root folder:

```r
# extract_rds.R
# Path to the locked file
rds_path <- "oncology-coPilot/oncology-backend-minimal/data/britoc/source_repos/britroc-1-cn-analysis/copy_number_signatures/britroc_30kb_signature_data.rds"

if (!file.exists(rds_path)) {
  stop("RDS file not found at: ", rds_path)
}

# Load the data
data <- readRDS(rds_path)

# The object is likely a matrix or list. Let's inspect and save components.
# If it's a list, save relevant dataframes.
if (is.list(data)) {
  print(names(data))
  # Assuming 'signatures' or 'contributions' is what we want
  # You may need to adjust based on structure output
  capture.output(str(data), file="rds_structure.txt")
  
  # Try to save the main element (adjust index/name if needed after inspection)
  # For now, we try to save the whole thing if it's a simple list-like structure convertible to dataframe
  # or iterate names
} 

# Attempt to save as CSV (simplest assumption: it's a dataframe/matrix)
tryCatch({
  write.csv(as.data.frame(data), "britroc_signatures_extracted.csv")
  print("Success! Data saved to britroc_signatures_extracted.csv")
}, error = function(e) {
  print("Could not save directly as CSV. See rds_structure.txt for details.")
})
```

### Step 3: Re-Run the Master CSV Builder
Once `britroc_signatures_extracted.csv` is generated, I can merge it into `BRITROC_MASTER.csv`.

## ⚠️ Alternative Plan (TCGA Only)
If you cannot run R, we must abandon BriTROC-1 CN validation for now and rely solely on **TCGA-OV RNA-seq** (which is open access and Python-friendly) to validate the "DDR State" hypothesis.
