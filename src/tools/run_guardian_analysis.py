import pandas as pd
import os
import requests
import time
from tqdm import tqdm

# Define file paths
INPUT_PATH = "results/guardian_analysis_brca1.tsv"
OUTPUT_PATH = "results/guardian_analysis_brca1_pathogenicity.tsv"

# CommandCenter API Endpoint - CORRECTED TO V2
COMMAND_CENTER_URL = "https://crispro--command-center-v2-commandcenter-api.modal.run/workflow/assess_threat"

def assess_mutation_pathogenicity(gene, protein_change):
    """
    Calls the CommandCenter to assess the pathogenicity of a single mutation.
    """
    if pd.isna(protein_change):
        return None, "No protein change", None, None

    payload = {
        "gene_symbol": gene,
        "protein_change": protein_change
    }
    try:
        print(f"  -> Assessing {gene} {protein_change} using endpoint: {COMMAND_CENTER_URL}")
        response = requests.post(COMMAND_CENTER_URL, json=payload, timeout=120) # Increased timeout
        if response.status_code == 200:
            data = response.json()
            return data.get("is_pathogenic"), data.get("assessment_source"), data.get("pathogenicity_score"), data.get("error")
        else:
            print(f"  -> API Error for {gene} {protein_change}: {response.status_code} - {response.text}")
            return "API Error", response.status_code, None, response.text
    except requests.exceptions.RequestException as e:
        print(f"  -> Request Exception for {gene} {protein_change}: {e}")
        return "Request Failed", str(e), None, str(e)

def main():
    """
    Main function to run the BRCA1 patient analysis.
    """
    print("--- Starting Operation Guardian: BRCA1 Pathogenicity Assessment (Target: v2) ---")

    try:
        df = pd.read_csv(INPUT_PATH, sep='\t', low_memory=False)
    except FileNotFoundError:
        print(f"ERROR: Input file not found at {INPUT_PATH}")
        return

    # Add columns for the assessment results
    df['is_pathogenic'] = None
    df['assessment_source'] = None
    df['pathogenicity_score'] = None
    df['error_message'] = None

    print(f"Found {len(df)} mutations to assess.")
    
    # Use tqdm for a progress bar
    tqdm.pandas(desc="Assessing Mutations")
    
    results = df.progress_apply(
        lambda row: assess_mutation_pathogenicity(row['Hugo_Symbol'], row['HGVSp_Short']), 
        axis=1
    )

    # Unpack results tuple into separate columns
    df[['is_pathogenic', 'assessment_source', 'pathogenicity_score', 'error_message']] = pd.DataFrame(results.tolist(), index=df.index)

    # Save the results to a new TSV file
    df.to_csv(OUTPUT_PATH, sep='\t', index=False)
    print(f"--- Analysis complete. Results saved to {OUTPUT_PATH} ---")

if __name__ == "__main__":
    main()