import pandas as pd
import os
import requests
import time
from tqdm import tqdm
import lifelines
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
MUTATIONS_PATH = "data/ov_data_mutations.txt"
CLINICAL_PATH = "data/ov_tcga_clinical_data.tsv"
OUTPUT_PLOT_PATH = "results/guardian_survival_analysis_brca1.png"
COMMAND_CENTER_URL = "https://crispro--command-center-v2-commandcenter-api.modal.run/workflow/assess_threat"

# --- PART 1: PATHOGENICITY ANALYSIS ---

def assess_mutation_pathogenicity(gene, protein_change):
    """
    Calls the CommandCenter to assess the pathogenicity of a single mutation.
    Returns a tuple: (is_pathogenic, assessment_source, pathogenicity_score, error_message)
    """
    if pd.isna(protein_change):
        return False, "No protein change", None, "Missing protein change data"

    payload = {"gene_symbol": gene, "protein_change": protein_change}
    try:
        response = requests.post(COMMAND_CENTER_URL, json=payload, timeout=120)
        if response.status_code == 200:
            data = response.json()
            
            # --- DEFINITIVE FIX: Parse the nested assessment dictionary ---
            assessment_data = data.get("assessment")
            if assessment_data:
                is_pathogenic = assessment_data.get("is_pathogenic", False)
                source = assessment_data.get("assessment_source", "Unknown")
                score = assessment_data.get("pathogenicity_score")
            else:
                is_pathogenic = False
                source = "No Assessment"
                score = None

            return (is_pathogenic, source, score, data.get("error"))
        else:
            return False, "API Error", None, f"{response.status_code}: {response.text}"
    except requests.exceptions.RequestException as e:
        return False, "Request Failed", None, str(e)


def run_pathogenicity_analysis(df):
    """
    Runs the pathogenicity analysis on the dataframe and adds results columns.
    """
    print("--- Starting Part 1: Pathogenicity Assessment ---")
    tqdm.pandas(desc="Assessing Mutations")
    
    results = df.progress_apply(
        lambda row: assess_mutation_pathogenicity(row['Hugo_Symbol'], row['HGVSp_Short']),
        axis=1
    )
    df[['is_pathogenic', 'assessment_source', 'pathogenicity_score', 'error_message']] = pd.DataFrame(results.tolist(), index=df.index)
    
    # --- DEBUGGING STEP ---
    print("\n--- Raw Pathogenicity Results from API ---")
    print(df[['HGVSp_Short', 'is_pathogenic', 'assessment_source', 'error_message']])
    # --- END DEBUGGING STEP ---

    print("--- Pathogenicity Assessment Complete ---")
    return df

# --- PART 2: SURVIVAL ANALYSIS ---

def run_survival_analysis(df):
    """
    Performs and plots a survival analysis on the dataframe.
    """
    print("\n--- Starting Part 2: Survival Analysis ---")
    
    df['OS_MONTHS'] = pd.to_numeric(df['Overall Survival (Months)'], errors='coerce')
    df['OS_STATUS'] = df['Overall Survival Status'].apply(lambda x: 1 if str(x).strip() == '1:DECEASED' else 0)
    
    df_survival = df.dropna(subset=['OS_MONTHS', 'OS_STATUS']).copy()
    
    print(f"Found {len(df_survival)} patients with complete survival data.")

    pathogenic_cohort = df_survival[df_survival['is_pathogenic'] == True]
    non_pathogenic_cohort = df_survival[df_survival['is_pathogenic'] == False]

    print(f"Pathogenic cohort size: {len(pathogenic_cohort)}")
    print(f"Non-Pathogenic/Other cohort size: {len(non_pathogenic_cohort)}")

    if len(pathogenic_cohort) < 2 or len(non_pathogenic_cohort) < 2:
        print("Warning: One or both cohorts have fewer than 2 patients. Cannot generate a comparative plot or run log-rank test.")
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(10, 6))

    kmf.fit(pathogenic_cohort['OS_MONTHS'], event_observed=pathogenic_cohort['OS_STATUS'], label=f"BRCA1 Pathogenic (n={len(pathogenic_cohort)})")
    kmf.plot_survival_function()

    kmf.fit(non_pathogenic_cohort['OS_MONTHS'], event_observed=non_pathogenic_cohort['OS_STATUS'], label=f"BRCA1 Non-Pathogenic / Other (n={len(non_pathogenic_cohort)})")
    kmf.plot_survival_function()

    plt.title('BRCA1 Pathogenicity and Overall Survival in Ovarian Cancer')
    plt.xlabel('Time (Months)')
    plt.ylabel('Overall Survival Probability')
    plt.grid(True)
    
    plt.savefig(OUTPUT_PLOT_PATH)
    print(f"--- Survival plot saved to {OUTPUT_PLOT_PATH} ---")
    
    results = lifelines.statistics.logrank_test(
        pathogenic_cohort['OS_MONTHS'],
        non_pathogenic_cohort['OS_MONTHS'],
        event_observed_A=pathogenic_cohort['OS_STATUS'],
        event_observed_B=non_pathogenic_cohort['OS_STATUS']
    )
    print("\n--- Log-Rank Test Results ---")
    results.print_summary()

# --- MAIN EXECUTION ---

def main():
    """
    Main function to run the full analysis pipeline.
    """
    print("====== Starting Operation Guardian: Full Unified Analysis ======")
    try:
        df_mutations = pd.read_csv(MUTATIONS_PATH, sep='\t', comment='#', low_memory=False)
        df_clinical = pd.read_csv(CLINICAL_PATH, sep='\t', low_memory=False)
    except FileNotFoundError as e:
        print(f"FATAL ERROR: Input file not found. {e}")
        return

    df_brca1 = df_mutations[df_mutations['Hugo_Symbol'] == 'BRCA1'].copy()
    
    # FIX: Create a common 'Patient ID' column to enable the merge.
    df_brca1['Patient ID'] = df_brca1['Tumor_Sample_Barcode'].str.slice(0, 12)
    
    df = pd.merge(df_brca1, df_clinical, on='Patient ID', how='inner')

    print(f"Found {len(df)} BRCA1 mutations in patients with clinical data.")

    df_with_results = run_pathogenicity_analysis(df)
    
    run_survival_analysis(df_with_results)
    
    print("\n====== Full Unified Analysis Complete ======")

if __name__ == "__main__":
    main() 