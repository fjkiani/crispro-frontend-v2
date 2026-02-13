import pandas as pd
import numpy as np
import argparse
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib.pyplot as plt

def analyze_prognostic_value(scores_file, survival_file, output_prefix):
    print(f"Loading scores from {scores_file}...")
    scores = pd.read_csv(scores_file, index_col=0)
    
    print(f"Loading survival data from {survival_file}...")
    # Xena survival file structure: sample, OS, OS.time, etc.
    surv = pd.read_csv(survival_file, sep='\t', index_col='sample')
    
    # Merge (Inner join to find patients with both)
    # Ensure indices match format (TCGA-XX-XXXX)
    # Xena often uses full barcode, while some expression uses short. 
    # Zeta Protocol: Try direct merge first.
    merged = scores.join(surv, how='inner')
    print(f"Matched {len(merged)} patients for analysis.")
    
    if len(merged) < 10:
        print("❌ Not enough patients matched! Check ID formats (Barcode vs UUID).")
        return

    # Clean data
    # OS: 1=Dead, 0=Alive. OS.time: Days
    merged = merged.dropna(subset=['OS', 'OS.time'])
    
    results = []
    
    # Iterate pathways
    for pathway in ['DDR_Score', 'E2F_Score']:
        print(f"\nAnalyzing {pathway}...")
        
        # 1. Cox Regression (Continuous)
        cph = CoxPHFitter()
        # Subset data
        data = merged[[pathway, 'OS', 'OS.time']]
        cph.fit(data, duration_col='OS.time', event_col='OS')
        
        hr = cph.hazard_ratios_[pathway]
        p_val = cph.summary.loc[pathway, 'p']
        print(f"  Hazard Ratio: {hr:.4f} (p={p_val:.4f})")
        
        results.append({
            'Pathway': pathway,
            'HR': hr,
            'P_Value': p_val
        })
        
        # 2. Kaplan-Meier (Median Split)
        median = merged[pathway].median()
        high = merged[merged[pathway] > median]
        low = merged[merged[pathway] <= median]
        
        kmf = KaplanMeierFitter()
        plt.figure(figsize=(10, 6))
        
        kmf.fit(high['OS.time'], high['OS'], label=f"High {pathway} (n={len(high)})")
        kmf.plot_survival_function()
        
        kmf.fit(low['OS.time'], low['OS'], label=f"Low {pathway} (n={len(low)})")
        kmf.plot_survival_function()
        
        plt.title(f"TCGA-OV Survival by {pathway} (Median Split)")
        plt.ylabel("Survival Probability")
        plt.xlabel("Days")
        plt.savefig(f"{output_prefix}_{pathway}_KM.png")
        print(f"  KM Plot saved to {output_prefix}_{pathway}_KM.png")
        plt.close()

    # Save summary
    res_df = pd.DataFrame(results)
    res_df.to_csv(f"{output_prefix}_summary.csv", index=False)
    print(f"✅ Analysis complete. Summary saved to {output_prefix}_summary.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, help="Pathway scores CSV")
    parser.add_argument("--survival", required=True, help="Survival data (Xena format)")
    parser.add_argument("--output", default="tcga_ov_validation", help="Output prefix")
    args = parser.parse_args()
    
    analyze_prognostic_value(args.scores, args.survival, args.output)
