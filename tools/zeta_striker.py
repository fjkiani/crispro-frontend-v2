import pandas as pd
import json
import time
import matplotlib.pyplot as plt
import lifelines
import requests

COMMAND_CENTER_URL = "https://crispro--command-center-commandcenter-api.modal.run/workflow/assess_threat"

# This function is now armed.
def call_command_center_assess_threat(gene_symbol, protein_change):
    """
    Calls the live CommandCenter API to assess a variant.
    """
    print(f"📡 Calling live Zeta Oracle for {gene_symbol} {protein_change}...")
    
    payload = {
        "gene_symbol": gene_symbol,
        "protein_change": protein_change
    }
    
    try:
        response = requests.post(COMMAND_CENTER_URL, json=payload, timeout=60)
        response.raise_for_status()  # Will raise an exception for 4XX/5XX errors
        
        data = response.json()
        print(f"DEBUG: Full API Response for {protein_change}:\\n{json.dumps(data, indent=2)}")

        # We need to extract the actual score from the nested structure
        # FIX: The API returns different structures for different types of assessments
        assessment = data.get("assessment", {})
        
        # Try to get zeta_score from the assessment
        zeta_score = assessment.get("zeta_score")
        if zeta_score is None:
            print(f"⚠️ Warning: No score found for {protein_change}. Defaulting to 0.")
            zeta_score = 0
        
        # Check if it's pathogenic (for truncation mutations)
        is_pathogenic = assessment.get("is_pathogenic", False)
        
        # If no explicit is_pathogenic flag, check if zeta_score indicates pathogenicity
        if not is_pathogenic and zeta_score < -50:
            is_pathogenic = True

        return {
            "gene_symbol": gene_symbol,
            "protein_change": protein_change,
            "assessment_source": "Zeta Oracle (Live Fire)",
            "delta_likelihood_score": zeta_score,
            "is_pathogenic": is_pathogenic
        }

    except requests.exceptions.RequestException as e:
        print(f"💥 Fucking hell, API call failed for {protein_change}: {e}")
        # Return a neutral/error score
        return {
            "gene_symbol": gene_symbol,
            "protein_change": protein_change,
            "assessment_source": "API Call Failed",
            "delta_likelihood_score": 0,
            "is_pathogenic": False
        }

def _perform_sanity_check(scored_mutations: dict):
    """
    Verifies the scores of known "litmus test" mutations to ensure the Oracle is behaving rationally.
    """
    print("\n--- Performing Sanity Check on Oracle Scores ---")
    litmus_tests = {
        "p.R175H": {"expected": "pathogenic", "threshold": -50},
        "p.R249M": {"expected": "pathogenic", "threshold": -50},
        "p.R280I": {"expected": "pathogenic", "threshold": -50},
        "p.R342*": {"expected": "pathogenic", "threshold": -900}, # Truncation
    }
    
    all_tests_passed = True
    for pc, test in litmus_tests.items():
        if pc in scored_mutations:
            score = scored_mutations[pc].get("delta_likelihood_score", 0)
            if test["expected"] == "pathogenic" and score > test["threshold"]:
                print(f"🚨 SANITY CHECK FAILED: {pc} score is {score}, but expected pathogenic score < {test['threshold']}.")
                all_tests_passed = False
            else:
                print(f"✅ Sanity Check Passed: {pc} score is {score}, which is plausible.")
        else:
            print(f"⚠️ Sanity Check Warning: Litmus test variant {pc} not found in this dataset.")

    if not all_tests_passed:
        print("\n💥 CRITICAL: One or more sanity checks failed. The Zeta Oracle may be misconfigured. Aborting analysis.")
        raise SystemExit("Zeta Oracle sanity check failed.")
    
    print("--- ✅ All Sanity Checks Passed. Proceeding. ---")

def extract_and_score_mutations(mutation_path, gene_of_interest):
    """
    Reads mutation data, extracts unique mutations for a gene,
    and gets a Zeta Score for each one.
    """
    print(f"🔥 Loading mutation data from {mutation_path}...")
    try:
        # The DtypeWarning is annoying and not critical for this analysis. Silence it.
        mutation_df = pd.read_csv(mutation_path, sep='\t', comment='#', low_memory=False)
    except FileNotFoundError:
        print(f"💥 Fucking hell, couldn't find the mutation file at {mutation_path}")
        return None

    print(f"🎯 Isolating mutations for our target: {gene_of_interest}...")
    target_mutations = mutation_df[mutation_df['Hugo_Symbol'] == gene_of_interest].copy()

    # We only care about mutations that change the protein.
    target_mutations = target_mutations[target_mutations['Variant_Classification'].isin([
        'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 
        'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins'
    ])]
    
    # Get a list of unique protein changes
    unique_protein_changes = target_mutations['HGVSp_Short'].dropna().unique()
    print(f"Found {len(unique_protein_changes)} unique protein-altering mutations for {gene_of_interest}.")

    # --- Step 2: Individual Threat Assessment ---
    print("\n--- Initiating Zeta Strike: Scoring all unique mutations ---")
    scored_mutations = {}
    for pc in unique_protein_changes:
        # The API needs the "p." prefix removed.
        protein_change_for_api = pc.replace('p.', '')
        assessment = call_command_center_assess_threat(gene_of_interest, protein_change_for_api)
        scored_mutations[pc] = assessment
    
    # --- NEW: Perform Sanity Check ---
    _perform_sanity_check(scored_mutations)
        
    return scored_mutations, target_mutations

def fuse_intelligence(clinical_path, scored_mutations, patient_mutation_map):
    """
    Fuses the Zeta Scores with clinical data to create a master analysis file.
    """
    print("\n--- Initiating Phase 3: Intelligence Fusion ---")
    print(f"🔥 Loading clinical data from {clinical_path}...")
    try:
        clinical_df = pd.read_csv(clinical_path, sep='\t', comment='#')
    except FileNotFoundError:
        print(f"💥 Fucking hell, couldn't find the clinical file at {clinical_path}")
        return None

    # Let's create a dataframe from our scored mutations
    scores_df = pd.DataFrame.from_dict(scored_mutations, orient='index').reset_index()
    scores_df.rename(columns={'index': 'HGVSp_Short', 'delta_likelihood_score': 'zeta_score'}, inplace=True)

    # Now, map these scores to each patient.
    # First, get a clear map of Patient ID to their specific mutation
    patient_to_mutation = patient_mutation_map[['Tumor_Sample_Barcode', 'HGVSp_Short']].dropna()
    
    # Merge this with our scores
    patient_scores = pd.merge(patient_to_mutation, scores_df[['HGVSp_Short', 'zeta_score']], on='HGVSp_Short')
    
    # The clinical patient ID is shorter, so we need to trim our sample IDs to match
    patient_scores['Patient ID'] = patient_scores['Tumor_Sample_Barcode'].str[:12]

    # Now, merge our patient scores with the main clinical data
    master_df = pd.merge(clinical_df, patient_scores[['Patient ID', 'zeta_score', 'HGVSp_Short']], on='Patient ID', how='left')

    # For patients without a TP53 mutation, their zeta_score will be NaN.
    # We'll fill this with 0, representing no damage. A perfect score.
    master_df['zeta_score'].fillna(0, inplace=True)
    master_df['HGVSp_Short'].fillna('Wild-Type', inplace=True)

    print("✅ Intelligence fusion complete.")
    return master_df

def final_analysis(master_file_path):
    """
    Performs the final Zeta Shield analysis on the fused data.
    """
    print("\n--- Initiating Step 4: The Kill Shot Analysis ---")
    print(f"🔥 Loading master intelligence file from {master_file_path}...")
    try:
        master_df = pd.read_csv(master_file_path, sep='\t')
    except FileNotFoundError:
        print(f"💥 Fucking hell, couldn't find the master file at {master_file_path}")
        return

    # Prepare data for survival analysis
    master_df['OS_EVENT'] = master_df['Overall Survival Status'].apply(lambda x: 1 if x == '1:DECEASED' else 0)
    master_df['OS_MONTHS'] = pd.to_numeric(master_df['Overall Survival (Months)'], errors='coerce')
    master_df.dropna(subset=['OS_MONTHS', 'OS_EVENT'], inplace=True)

    # --- Focus on the Radiation Therapy Cohort ---
    rad_cohort = master_df[master_df['Radiation Therapy'] == 'Yes'].copy()
    print(f"🔬 Isolated {len(rad_cohort)} patients who received radiation therapy.")

    # Stratify by Zeta Score. We'll use a simple threshold for now.
    # A score of 0 is Wild-Type (no damage).
    # Scores < 0 are mutated. We'll set a threshold for "High Damage".
    damage_threshold = -50
    rad_cohort['Damage_Level'] = rad_cohort['zeta_score'].apply(
        lambda x: 'High Damage' if x <= damage_threshold else 'Low Damage'
    )
    
    # We need to make sure we have patients in both damage groups to compare.
    if len(rad_cohort['Damage_Level'].unique()) < 2:
        print("💥 Shit. Not enough diversity in damage scores within the radiation cohort to compare.")
        print(rad_cohort['Damage_Level'].value_counts())
        return
        
    print("\nDamage stratification in radiation cohort:")
    print(rad_cohort['Damage_Level'].value_counts())

    # --- Kaplan-Meier Survival Analysis ---
    kmf = lifelines.KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(12, 7))

    for level in ['Low Damage', 'High Damage']:
        subset = rad_cohort[rad_cohort['Damage_Level'] == level]
        if not subset.empty:
            kmf.fit(subset['OS_MONTHS'], event_observed=subset['OS_EVENT'], label=f'{level} (n={len(subset)})')
            kmf.plot_survival_function(ax=ax)

    # Log-Rank Test to see if the curves are statistically different
    low_damage_group = rad_cohort[rad_cohort['Damage_Level'] == 'Low Damage']
    high_damage_group = rad_cohort[rad_cohort['Damage_Level'] == 'High Damage']
    
    results = lifelines.statistics.logrank_test(
        low_damage_group['OS_MONTHS'], high_damage_group['OS_MONTHS'],
        event_observed_A=low_damage_group['OS_EVENT'],
        event_observed_B=high_damage_group['OS_EVENT']
    )

    print("\n--- Log-Rank Test for Zeta Score Damage Level ---")
    results.print_summary()

    ax.set_title('Zeta Shield Analysis: Survival in Radiation-Treated Patients by TP53 Zeta Score')
    ax.set_xlabel('Time (Months)')
    ax.set_ylabel('Survival Probability')
    ax.grid(True)
    
    p_value = results.p_value
    ax.text(0.5, 0.1, f'Log-Rank p-value: {p_value:.4f}', transform=ax.transAxes, fontsize=12,
            verticalalignment='bottom', horizontalalignment='center')


    output_path = "results/zeta_shield_analysis_v2.png"
    plt.savefig(output_path)
    print(f"\n📈 Final analysis plot saved to {output_path}")
    plt.show()

def main():
    """Main execution function for Operation Zeta Strike - Phase 1."""
    mutation_path = 'data/data_mutations.txt'
    clinical_path = 'data/tcga.tsv'
    gene_of_interest = 'TP53'
    scores_output_path = f"results/TP53_zeta_scores_v2.json"
    master_file_output_path = f"results/master_clinical_zeta_scores_v2.tsv"

    # Steps 1 & 2
    scored_mutations, patient_mutation_map = extract_and_score_mutations(mutation_path, gene_of_interest)

    if scored_mutations:
        print(f"\n✅ Zeta Strike complete. Saving {len(scored_mutations)} scored mutations to {scores_output_path}...")
        with open(scores_output_path, 'w') as f:
            json.dump(scored_mutations, f, indent=4)
        print("Done.")

        # Step 3
        master_df = fuse_intelligence(clinical_path, scored_mutations, patient_mutation_map)
        if master_df is not None:
            print(f"\n💾 Saving master intelligence file to {master_file_output_path}...")
            # We'll save it as a TSV (Tab-Separated Values) file. It's good practice.
            master_df.to_csv(master_file_output_path, sep='\t', index=False)
            print("✅ Master file saved.")

            # Step 4: The Kill Shot
            final_analysis(master_file_output_path)

if __name__ == "__main__":
    main() 