import pandas as pd
import requests
import os
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
import lifelines

# --- Constants ---
# Services
COMMAND_CENTER_URL = "https://crispro--command-center-v2-commandcenter-api.modal.run/workflow/assess_threat"

# Data Files
MUTATION_DATA_PATH = "data/data_mutations.txt"
CLINICAL_DATA_PATH = "data/tcga.tsv"

# Output Files
MASTER_ANALYSIS_FILE = "results/rad_onc_analysis_data.tsv"
SURVIVAL_PLOT_FILE = "results/rad_onc_survival_analysis.png"

# Analysis Parameters
GENE_OF_INTEREST = "TP53"

# --- Data Ingestion ---

def load_mutation_data(path: str, gene: str) -> pd.DataFrame:
    """Loads the mutation data and filters for the gene of interest."""
    print(f"🧬 Loading mutation data from {path}...")
    try:
        df = pd.read_csv(path, sep='\t', comment='#', low_memory=False)
        print(f"✅ Loaded {len(df)} mutations.")
        
        print(f"🔎 Filtering for {gene} mutations...")
        gene_mutations = df[df['Hugo_Symbol'] == gene].copy()
        print(f"Found {len(gene_mutations)} mutations for {gene}.")
        return gene_mutations
    except FileNotFoundError:
        print(f"❌ ERROR: Mutation data file not found at {path}")
        return None

def load_clinical_data(path: str) -> pd.DataFrame:
    """Loads the clinical data."""
    print(f"🏥 Loading clinical data from {path}...")
    try:
        df = pd.read_csv(path, sep='\t', comment='#', low_memory=False)
        print(f"✅ Loaded clinical data for {len(df)} patients.")
        return df
    except FileNotFoundError:
        print(f"❌ ERROR: Clinical data file not found at {path}")
        return None

# --- Triumvirate Assessment ---

def get_pathogenicity(gene_symbol: str, protein_change: str) -> bool | None:
    """
    Calls the CommandCenter to get a pathogenicity assessment for a single variant.
    Returns True if pathogenic, False if benign, None if assessment fails.
    """
    if not isinstance(protein_change, str) or not protein_change.startswith('p.'):
        return None # Skip invalid protein change formats

    # The API expects the 'p.' prefix to be removed from the protein change string.
    protein_change_for_api = protein_change.replace('p.', '')

    payload = {
        "gene_symbol": gene_symbol,
        "protein_change": protein_change_for_api
    }
    
    try:
        response = requests.post(COMMAND_CENTER_URL, json=payload, timeout=120)
        response.raise_for_status()
        
        data = response.json()
        assessment = data.get("assessment", {})
        
        if assessment.get("error"):
            print(f"⚠️ API returned an error for {protein_change}: {assessment['error']}")
            return None
            
        return assessment.get("is_pathogenic")

    except requests.exceptions.RequestException as e:
        print(f"💥 API call failed for {protein_change}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred during API call for {protein_change}: {e}")
        return None

def assess_all_mutations(mutations_df: pd.DataFrame, gene: str) -> dict:
    """
    Iterates through all unique mutations and gets their pathogenicity.
    """
    unique_protein_changes = mutations_df['HGVSp_Short'].dropna().unique()
    print(f"🔬 Assessing {len(unique_protein_changes)} unique {gene} mutations...")
    
    pathogenicity_map = {}
    
    for pc in tqdm(unique_protein_changes, desc=f"Assessing {gene} Mutations"):
        is_pathogenic = get_pathogenicity(gene, pc)
        pathogenicity_map[pc] = is_pathogenic
        time.sleep(0.1) # Small delay to be kind to the API
        
    print("✅ Assessment complete.")
    return pathogenicity_map

# --- Analysis & Output ---

def create_master_table(clinical_df: pd.DataFrame, mutations_df: pd.DataFrame, pathogenicity_map: dict) -> pd.DataFrame:
    """
    Merges clinical data with mutation and pathogenicity data to create a master table.
    """
    print("🤝 Merging dataframes to create master analysis table...")
    
    # Map the pathogenicity results to the mutations dataframe
    mutations_df['is_pathogenic'] = mutations_df['HGVSp_Short'].map(pathogenicity_map)

    # We only need patient ID, mutation, and pathogenicity status from the mutation table
    patient_mutation_info = mutations_df[['Tumor_Sample_Barcode', 'HGVSp_Short', 'is_pathogenic']].copy()
    
    # The clinical patient ID is a substring of the sample barcode
    patient_mutation_info['Patient ID'] = patient_mutation_info['Tumor_Sample_Barcode'].str[:12]
    
    # Merge the clinical data with our mutation info
    master_df = pd.merge(clinical_df, patient_mutation_info, on='Patient ID', how='left')

    # For patients without a TP53 mutation, is_pathogenic will be NaN.
    # We'll classify them as not pathogenic (False). Wild-Type is not pathogenic.
    master_df['is_pathogenic'].fillna(False, inplace=True)
    master_df['HGVSp_Short'].fillna('Wild-Type', inplace=True)
    
    print(f"✅ Master table created with {len(master_df)} records.")
    
    # Save the master table
    print(f"💾 Saving master table to {MASTER_ANALYSIS_FILE}...")
    if not os.path.exists("results"):
        os.makedirs("results")
    master_df.to_csv(MASTER_ANALYSIS_FILE, sep='\t', index=False)
    print("✅ Save complete.")
    
    return master_df

def perform_survival_analysis(master_df: pd.DataFrame):
    """
    Performs Kaplan-Meier survival analysis on the master dataframe
    and generates a plot.
    """
    print("\\n--- 🔬 Starting Phase 2: Survival Analysis ---")
    
    # Prepare data for lifelines
    master_df['OS_EVENT'] = master_df['Overall Survival Status'].apply(lambda x: 1 if x == '1:DECEASED' else 0)
    master_df['OS_MONTHS'] = pd.to_numeric(master_df['Overall Survival (Months)'], errors='coerce')
    master_df.dropna(subset=['OS_MONTHS', 'OS_EVENT'], inplace=True)
    
    # Stratify patients into the four cohorts
    master_df['cohort'] = 'TP53 Wild-Type, No Radiation'
    master_df.loc[(master_df['is_pathogenic'] == True) & (master_df['Radiation Therapy'] == 'No'), 'cohort'] = 'TP53 Pathogenic, No Radiation'
    master_df.loc[(master_df['is_pathogenic'] == False) & (master_df['Radiation Therapy'] == 'Yes'), 'cohort'] = 'TP53 Wild-Type, Received Radiation'
    master_df.loc[(master_df['is_pathogenic'] == True) & (master_df['Radiation Therapy'] == 'Yes'), 'cohort'] = 'TP53 Pathogenic, Received Radiation'
    
    print("Cohorts created:")
    print(master_df['cohort'].value_counts())

    # Perform Kaplan-Meier analysis
    kmf = lifelines.KaplanMeierFitter()
    
    plt.figure(figsize=(12, 8))
    
    for name, grouped_df in master_df.groupby('cohort'):
        kmf.fit(grouped_df['OS_MONTHS'], event_observed=grouped_df['OS_EVENT'], label=f'{name} (n={len(grouped_df)})')
        kmf.plot_survival_function()

    plt.title('Survival Analysis of Lung Cancer Patients by TP53 Status and Radiation Therapy')
    plt.xlabel('Time (Months)')
    plt.ylabel('Overall Survival Probability')
    plt.grid(True)
    plt.tight_layout()
    
    print(f"💾 Saving survival plot to {SURVIVAL_PLOT_FILE}...")
    plt.savefig(SURVIVAL_PLOT_FILE)
    print("✅ Plot saved.")
    

def main():
    """Main orchestration function for the analysis."""
    print("--- 🔥 Operation Phoenix: Radiation Oncology Analysis 🔥 ---")
    
    # Phase 1: Data Processing & Enrichment
    mutations = load_mutation_data(MUTATION_DATA_PATH, GENE_OF_INTEREST)
    if mutations is None:
        return

    clinical = load_clinical_data(CLINICAL_DATA_PATH)
    if clinical is None:
        return
        
    pathogenicity_results = assess_all_mutations(mutations, GENE_OF_INTEREST)
    
    master_table = create_master_table(clinical, mutations, pathogenicity_results)

    print("\\n--- ✅ Operation Phoenix: Phase 1 Complete ---")
    
    perform_survival_analysis(master_table)
    
    print("\\n--- ✅ Operation Phoenix: Analysis Complete ---")


if __name__ == "__main__":
    main() 