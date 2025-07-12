import gzip
import os
import pandas as pd

# --- Configuration ---
TARGET_GENES = [
    "TP53", "KRAS", "BRAF", "EGFR", "PIK3CA", 
    "PTEN", "BRCA1", "BRCA2", "ATM"
]
# Using the new, more comprehensive data source
INPUT_FILE = "data/adjudicator_training/variant_summary.txt.gz"
OUTPUT_FILE = "data/adjudicator_training/clinvar_missense_labels.csv"

def main():
    """
    Main orchestration function to parse the ClinVar variant_summary.txt.gz 
    file and extract a high-quality, labeled dataset for training.
    """
    print("🚀 Starting Operation Adjudicator: Phase 1.1 (Final Strategy)")
    print(f"  -> Input file: {INPUT_FILE}")

    if not os.path.exists(INPUT_FILE):
        print(f"  -> FATAL ERROR: variant_summary.txt.gz file not found at '{INPUT_FILE}'. Please download it first.")
        return

    print("STEP 1: Loading data into pandas... (This may take a few minutes)")
    df = pd.read_csv(INPUT_FILE, sep='\t', compression='gzip', low_memory=False, on_bad_lines='skip')
    print(f"  -> Successfully loaded {len(df)} total records from source.")

    print("STEP 2: Filtering records based on ground truth schema...")

    # Filter 1: Keep only records for our target genes
    df_filtered = df[df['GeneSymbol'].isin(TARGET_GENES)].copy()
    print(f"  -> {len(df_filtered)} records remaining after filtering for target genes.")

    # Filter 2: CORRECTED - Keep records where 'Type' is 'single nucleotide variant'
    # and the 'Name' column indicates a protein change (p. notation). This is the
    # ground truth method for identifying missense variants in this file.
    df_filtered = df_filtered[
        (df_filtered['Type'] == 'single nucleotide variant') & 
        (df_filtered['Name'].str.contains(r'\(p\.', na=False))
    ].copy()
    print(f"  -> {len(df_filtered)} records remaining after identifying SNVs with protein changes.")

    # Filter 3: CORRECTED - Relax the review status filter to be more inclusive
    # of high-quality entries that may not have the single "expert panel" status.
    # We will accept any record where assertion criteria have been provided.
    df_filtered = df_filtered[df_filtered['ReviewStatus'].str.contains('criteria provided', na=False)].copy()
    print(f"  -> {len(df_filtered)} records remaining after filtering for 'criteria provided' review status.")

    # Filter 4: Keep only records with clear Pathogenic or Benign labels
    df_filtered['significance_mapped'] = df_filtered['ClinicalSignificance'].apply(
        lambda x: 1 if isinstance(x, str) and 'pathogenic' in x.lower() else (0 if isinstance(x, str) and 'benign' in x.lower() else -1)
    )
    # Also include 'Conflicting classifications of pathogenicity' as it contains valuable data we can sift through
    df_filtered = df_filtered[
        df_filtered['significance_mapped'].isin([0, 1]) | 
        (df_filtered['ClinicalSignificance'].str.contains('Conflicting', na=False))
    ].copy()

    # If a record has conflicting data, but contains the word 'pathogenic', we will consider it pathogenic for now.
    # This is an assumption we can refine later.
    df_filtered.loc[df_filtered['ClinicalSignificance'].str.contains('Conflicting', na=False) & df_filtered['ClinicalSignificance'].str.contains('Pathogenic', na=False), 'significance_mapped'] = 1
    df_filtered.loc[df_filtered['ClinicalSignificance'].str.contains('Conflicting', na=False) & df_filtered['ClinicalSignificance'].str.contains('Benign', na=False), 'significance_mapped'] = 0

    # Now, drop any remaining records that are not 0 or 1
    df_filtered = df_filtered[df_filtered['significance_mapped'].isin([0, 1])].copy()
    print(f"  -> {len(df_filtered)} records remaining after filtering for pathogenic/benign significance.")

    # --- Data Standardization and Saving ---
    print("\nSTEP 3: Standardizing and saving final dataset...")

    # CORRECTED: More robust regex to handle different HGVS formats
    df_filtered['hgvsp'] = df_filtered['Name'].str.extract(r'\((p\..*?)\)')[0]
    df_filtered['hgvsc'] = df_filtered['Name'].str.extract(r'\:(c\..*?)\s*\(')[0]

    # Select and rename columns for the final output
    df_final = df_filtered[[
        'GeneSymbol', 'hgvsp', 'hgvsc', 'VariationID', 
        'ClinicalSignificance', 'significance_mapped', 'ReviewStatus'
    ]].rename(columns={
        'GeneSymbol': 'gene',
        'VariationID': 'clinvar_id',
        'ClinicalSignificance': 'significance_raw',
        'ReviewStatus': 'review_status'
    })

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(OUTPUT_FILE)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_final.to_csv(OUTPUT_FILE, index=False)
    print(f"  -> Successfully saved {len(df_final)} standardized records to '{OUTPUT_FILE}'")

    # --- Data Validation ---
    print("\n--- Data Summary ---")
    print(f"Total Records: {len(df_final)}")
    if not df_final.empty:
        print("Class Distribution:")
        print(df_final['significance_mapped'].value_counts())
    print("--------------------")

    print("\n🔥 Operation Adjudicator: Phase 1.1 Complete.")

if __name__ == "__main__":
    main() 