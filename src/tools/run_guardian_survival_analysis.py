import pandas as pd
import lifelines
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

# Define file paths
INPUT_PATH = "results/guardian_analysis_brca1_pathogenicity.tsv"
OUTPUT_PLOT_PATH = "results/guardian_survival_analysis_brca1.png"

def run_survival_analysis():
    """
    Performs and plots a survival analysis based on BRCA1 mutation pathogenicity.
    """
    print("--- Starting Operation Guardian: Survival Analysis ---")
    
    try:
        df = pd.read_csv(INPUT_PATH, sep='\t', low_memory=False)
    except FileNotFoundError:
        print(f"ERROR: Analysis input file not found at {INPUT_PATH}")
        return

    # Prepare data for lifelines
    # Ensure 'Overall Survival (Months)' and 'Overall Survival Status' are numeric
    df['OS_MONTHS'] = pd.to_numeric(df['Overall Survival (Months)'], errors='coerce')
    
    # Convert 'Overall Survival Status' to binary (1 for deceased, 0 for living)
    # 1:DECEASED, 0:LIVING
    df['OS_STATUS'] = df['Overall Survival Status'].apply(lambda x: 1 if str(x).strip() == '1:DECEASED' else 0)
    
    # *** FIX: Convert string 'True'/'False' to actual booleans ***
    # The value read from the TSV is a string, not a boolean. This must be converted.
    df['is_pathogenic'] = df['is_pathogenic'].apply(lambda x: str(x).lower() == 'true')

    # Drop rows where survival data is missing
    df.dropna(subset=['OS_MONTHS', 'OS_STATUS'], inplace=True)
    
    print(f"Found {len(df)} patients with complete survival data.")

    # Create two cohorts by directly comparing the string values
    pathogenic_cohort = df[df['is_pathogenic'].astype(str).str.lower() == 'true']
    non_pathogenic_cohort = df[df['is_pathogenic'].astype(str).str.lower() != 'true']

    print(f"Pathogenic cohort size: {len(pathogenic_cohort)}")
    print(f"Non-Pathogenic/Other cohort size: {len(non_pathogenic_cohort)}")

    # Perform Kaplan-Meier survival analysis
    kmf = KaplanMeierFitter()

    plt.figure(figsize=(10, 6))

    # Fit and plot for the pathogenic cohort
    if not pathogenic_cohort.empty:
        kmf.fit(pathogenic_cohort['OS_MONTHS'], event_observed=pathogenic_cohort['OS_STATUS'], label='BRCA1 Pathogenic')
        kmf.plot_survival_function()

    # Fit and plot for the non-pathogenic cohort
    if not non_pathogenic_cohort.empty:
        kmf.fit(non_pathogenic_cohort['OS_MONTHS'], event_observed=non_pathogenic_cohort['OS_STATUS'], label='BRCA1 Non-Pathogenic / Other')
        kmf.plot_survival_function()

    # Add plot details
    plt.title('BRCA1 Pathogenicity and Overall Survival in Ovarian Cancer')
    plt.xlabel('Time (Months)')
    plt.ylabel('Overall Survival Probability')
    plt.grid(True)
    
    # Save the plot
    plt.savefig(OUTPUT_PLOT_PATH)
    print(f"--- Survival plot saved to {OUTPUT_PLOT_PATH} ---")
    
    # Perform log-rank test to see if the difference is statistically significant
    if not pathogenic_cohort.empty and not non_pathogenic_cohort.empty:
        results = lifelines.statistics.logrank_test(
            pathogenic_cohort['OS_MONTHS'],
            non_pathogenic_cohort['OS_MONTHS'],
            event_observed_A=pathogenic_cohort['OS_STATUS'],
            event_observed_B=non_pathogenic_cohort['OS_STATUS']
        )
        print("\n--- Log-Rank Test Results ---")
        results.print_summary()

def main():
    run_survival_analysis()

if __name__ == "__main__":
    main() 