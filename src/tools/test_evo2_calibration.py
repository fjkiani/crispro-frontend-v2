import os
import requests
import re
import time
import json
from datetime import datetime
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import gzip

# --- Configuration ---
TEST_MODE_LIMIT = 20 # Set to 0 to run on all variants. Positive number limits the run.
# API endpoint for the real Evo2 model
MODAL_API_URL = "https://fjkiani--variant-analysis-evo2-evo2model-web.modal.run"
# Input VCF from Phase 1
FILTERED_VCF_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "clinvar_filtered_pathogenic_benign.vcf.gz")
# Output directory
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIRECTORY = os.path.join(OUTPUT_DIR, "figures")

def call_evo2_api_for_vcf_variant(chrom, pos, ref, alt):
    """
    Calls the deployed Evo2 model on Modal using VCF-style variant data.
    """
    # Ensure chromosome has 'chr' prefix, which the API expects.
    if not str(chrom).lower().startswith('chr'):
        chrom = f"chr{chrom}"

    try:
        payload = {
            "genome": "hg38",
            "chromosome": chrom,
            "variant_position": int(pos),
            "alternative": alt,
        }
        response = requests.post(MODAL_API_URL, json=payload, timeout=120)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"API call failed for variant {chrom}:{pos} {ref}>{alt}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred for variant {chrom}:{pos} {ref}>{alt}: {e}")
        return None

def run_clinvar_calibration():
    """
    Loads the filtered ClinVar VCF, gets real Evo2 delta_scores for each variant,
    and saves the results for analysis.
    """
    print("--- Starting ClinVar Calibration Experiment ---")
    
    results = []
    variants_to_process = []

    # First, read all variants from the VCF file
    print(f"Reading variants from {FILTERED_VCF_PATH}...")
    # Check if file exists before trying to open
    if not os.path.exists(FILTERED_VCF_PATH):
        print(f"ERROR: Input file not found at {FILTERED_VCF_PATH}")
        return None

    with gzip.open(FILTERED_VCF_PATH, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt_alleles, _, _, info = fields[:8]
            
            # Extract clinical significance array. It may be comma-separated.
            clinsig_match = re.search(r'CLNSIG=([^;]+)', info)
            if not clinsig_match:
                continue
            clinsig_values = clinsig_match.group(1).split(',')

            # Handle multi-allelic sites by splitting them.
            # We assume the number of alts matches the number of clinsig values.
            alts = alt_alleles.split(',')
            
            for i, alt in enumerate(alts):
                # Assign clinical significance. If lists don't match, use the first.
                clinsig_raw = clinsig_values[i] if i < len(clinsig_values) else clinsig_values[0]
                clinsig_raw = clinsig_raw.lower()

                if 'pathogenic' in clinsig_raw:
                    significance = 'Pathogenic'
                elif 'benign' in clinsig_raw:
                    significance = 'Benign'
                else:
                    continue # Skip 'uncertain', 'conflicting', etc.
            
                variants_to_process.append({
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "significance": significance
                })
    
    print(f"Found {len(variants_to_process)} total variants after splitting multi-allelic sites.")

    # Now, call the API for each variant
    variants_for_this_run = variants_to_process
    if TEST_MODE_LIMIT > 0:
        print(f"\n--- RUNNING IN TEST MODE: Analyzing up to {TEST_MODE_LIMIT} variants. ---\n")
        variants_for_this_run = variants_to_process[:TEST_MODE_LIMIT]

    for i, variant in enumerate(variants_for_this_run):
        print(f"Processing variant {i+1}/{len(variants_for_this_run)}: {variant['chrom']}:{variant['pos']}... ", end="")
        evo2_result = call_evo2_api_for_vcf_variant(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
        
        if evo2_result and 'delta_score' in evo2_result:
            results.append({
                "variant": f"{variant['chrom']}:{variant['pos']}",
                "delta_score": evo2_result['delta_score'],
                "significance": variant['significance']
            })
            print(f"Success. Delta Score: {evo2_result['delta_score']:.6f}")
        else:
            print("Failed.")
            
        time.sleep(1) # Be respectful to the API

    # Convert results to a DataFrame and save
    df = pd.DataFrame(results)
    if df.empty:
        print("No results obtained from API. Exiting.")
        return None

    output_file = os.path.join(OUTPUT_DIR, "data", "clinvar_calibration_results.tsv")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved raw results to {output_file}")
    
    return df

def generate_calibration_plot(df):
    """
    Generates a boxplot comparing the delta_score distributions for
    Pathogenic vs. Benign variants.
    """
    if df is None or df.empty:
        print("Cannot generate plot from empty dataframe.")
        return
        
    print("Generating calibration plot...")
    plt.figure(figsize=(8, 7))
    sns.set_theme(style="whitegrid")
    
    # Define order and palette for consistency
    order = ['Benign', 'Pathogenic']
    palette = {"Benign": "#7EC3E6", "Pathogenic": "#E67E7E"}

    sns.boxplot(data=df, x='significance', y='delta_score', order=order, palette=palette)
    
    plt.title('Evo2 Delta Score Calibration against ClinVar', fontsize=16, pad=20)
    plt.xlabel('ClinVar Clinical Significance', fontsize=12)
    plt.ylabel('Evo2 Model Delta Score', fontsize=12)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    
    # Perform Mann-Whitney U test
    benign_scores = df[df['significance'] == 'Benign']['delta_score'].dropna()
    pathogenic_scores = df[df['significance'] == 'Pathogenic']['delta_score'].dropna()
    
    if len(benign_scores) > 1 and len(pathogenic_scores) > 1:
        # Using 'less' because we hypothesize pathogenic scores are less (more negative) than benign scores
        stat, p_val = stats.mannwhitneyu(pathogenic_scores, benign_scores, alternative='less')
        
        # Display the p-value on the plot
        y_max = df['delta_score'].max()
        y_min = df['delta_score'].min()
        plt.text(0.5, 0.95, f'Mann-Whitney U test\np-value = {p_val:.2e}', 
                 ha='center', va='top', transform=plt.gca().transAxes, fontsize=12, 
                 bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

    plot_filename = os.path.join(FIGURES_DIRECTORY, "D_clinvar_calibration_boxplot.png")
    os.makedirs(os.path.dirname(plot_filename), exist_ok=True)
    plt.savefig(plot_filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Saved calibration plot to {plot_filename}")


def main():
    """Main function to run the calibration experiment."""
    if not os.path.exists(FIGURES_DIRECTORY):
        os.makedirs(FIGURES_DIRECTORY)
        print(f"Created directory: {FIGURES_DIRECTORY}")
        
    results_df = run_clinvar_calibration()
    if results_df is not None and not results_df.empty:
        generate_calibration_plot(results_df)
    print("\n--- ClinVar Calibration Experiment Finished ---")


if __name__ == "__main__":
    main() 