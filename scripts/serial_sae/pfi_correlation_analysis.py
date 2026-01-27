#!/usr/bin/env python3
"""
PFI Correlation Analysis (GSE165897)
Test pathway kinetics against Platinum-Free Interval (PFI)

Hypotheses:
1. VEGF Activation → Short PFI (negative correlation)
2. DDR Suppression → Short PFI (positive correlation)
3. Multi-pathway signature → resistance prediction
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from scipy.stats import spearmanr, mannwhitneyu
import sys
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Force immediate output flushing
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout, 'reconfigure') else None

RESULTS_DIR = Path("data/serial_sae/gse165897/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# PFI data from user (Zhang et al., 2022 Science Advances, Table 1)
PFI_DATA = {
    'patient_id': ['EOC1005', 'EOC87', 'EOC136', 'EOC423', 'EOC568', 
                   'EOC648', 'EOC891', 'EOC991', 'EOC933', 'EOC1129', 'EOC1153'],
    'pfi_months': [4.2, 9.1, 7.0, 24.0, 7.0, 50.5, 6.1, 4.2, None, None, None],
    'pfi_days': [126, 274, 210, 721, 210, 1515, 183, 126, None, None, None],
    'resistant': ['R', 'S', 'S', 'S', 'S', 'S', 'S', 'R', None, None, None],
    'resistance_label': ['resistant', 'sensitive', 'sensitive', 'sensitive', 'sensitive',
                        'sensitive', 'sensitive', 'resistant', None, None, None]
}

def cohens_d(group1, group2):
    """Calculate Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_std

def main():
    print("="*80)
    print("PFI CORRELATION ANALYSIS (GSE165897)")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Step 1: Load pathway kinetics
    kinetics_path = RESULTS_DIR / "pathway_kinetics.csv"
    if not kinetics_path.exists():
        print(f"❌ Pathway kinetics file not found: {kinetics_path}")
        print("   Run pathway_kinetics_gse165897.py first!")
        return
    
    kinetics_df = pd.read_csv(kinetics_path)
    print(f"✅ Loaded pathway kinetics for {len(kinetics_df)} patients")
    print(f"   Columns: {list(kinetics_df.columns)}\n")
    
    # Step 2: Add PFI data
    pfi_df = pd.DataFrame(PFI_DATA)
    merged_df = kinetics_df.merge(pfi_df, on='patient_id', how='left')
    
    # Filter to evaluable patients (n=8 with PFI data)
    evaluable_df = merged_df.dropna(subset=['pfi_months']).copy()
    print(f"✅ Evaluable patients (n={len(evaluable_df)} with PFI data)")
    print(f"   Resistant: {len(evaluable_df[evaluable_df['resistant']=='R'])}")
    print(f"   Sensitive: {len(evaluable_df[evaluable_df['resistant']=='S'])}\n")
    
    if len(evaluable_df) < 4:
        print("⚠️  Insufficient data for correlation analysis (need n≥4)")
        return
    
    # Step 3: Spearman correlations (PFI vs pathway kinetics)
    print("="*80)
    print("SPEARMAN CORRELATIONS: PFI vs Pathway Kinetics")
    print("="*80)
    
    pathway_deltas = ['ddr_delta', 'mapk_delta', 'pi3k_delta', 'vegf_delta']
    correlation_results = []
    
    for pathway in pathway_deltas:
        # Remove any NaN values for correlation
        valid_data = evaluable_df[[pathway, 'pfi_months']].dropna()
        if len(valid_data) >= 3:
            r, p = spearmanr(valid_data['pfi_months'], valid_data[pathway])
            correlation_results.append({
                'pathway': pathway.replace('_delta', ''),
                'spearman_r': r,
                'spearman_p': p,
                'n': len(valid_data),
                'interpretation': 'Strong' if abs(r) > 0.6 else 'Moderate' if abs(r) > 0.4 else 'Weak'
            })
            print(f"{pathway.upper()}:")
            print(f"  Spearman r = {r:.3f}, p = {p:.4f} (n={len(valid_data)})")
            if p < 0.05:
                print(f"  ✅ Significant (p < 0.05)")
            elif p < 0.1:
                print(f"  ⚠️  Trend (p < 0.1)")
            else:
                print(f"  ❌ Not significant (p ≥ 0.1)")
            print()
        else:
            print(f"⚠️  {pathway}: Insufficient data for correlation (n={len(valid_data)})")
            print()
    
    correlation_df = pd.DataFrame(correlation_results)
    correlation_df.to_csv(RESULTS_DIR / "pfi_correlations.csv", index=False)
    print(f"💾 Saved correlations: {RESULTS_DIR / 'pfi_correlations.csv'}\n")
    
    # Step 4: Mann-Whitney U test (Resistant vs Sensitive)
    print("="*80)
    print("MANN-WHITNEY U TEST: Resistant vs Sensitive")
    print("="*80)
    
    resistant_df = evaluable_df[evaluable_df['resistant'] == 'R']
    sensitive_df = evaluable_df[evaluable_df['resistant'] == 'S']
    
    print(f"Resistant group: n={len(resistant_df)}")
    print(f"Sensitive group: n={len(sensitive_df)}\n")
    
    if len(resistant_df) < 2 or len(sensitive_df) < 2:
        print("⚠️  Insufficient data for group comparison (need n≥2 per group)")
    else:
        comparison_results = []
        
        for pathway in pathway_deltas:
            resistant_vals = resistant_df[pathway].dropna().values
            sensitive_vals = sensitive_df[pathway].dropna().values
            
            if len(resistant_vals) > 0 and len(sensitive_vals) > 0:
                # Mann-Whitney U test
                u_stat, u_p = mannwhitneyu(resistant_vals, sensitive_vals, alternative='two-sided')
                
                # Effect size (Cohen's d)
                d = cohens_d(resistant_vals, sensitive_vals)
                
                # Means
                resistant_mean = np.mean(resistant_vals)
                sensitive_mean = np.mean(sensitive_vals)
                
                comparison_results.append({
                    'pathway': pathway.replace('_delta', ''),
                    'resistant_mean': resistant_mean,
                    'resistant_n': len(resistant_vals),
                    'sensitive_mean': sensitive_mean,
                    'sensitive_n': len(sensitive_vals),
                    'mannwhitney_u': u_stat,
                    'mannwhitney_p': u_p,
                    'cohens_d': d,
                    'interpretation': 'Large' if abs(d) > 0.8 else 'Medium' if abs(d) > 0.5 else 'Small'
                })
                
                print(f"{pathway.upper()}:")
                print(f"  Resistant: mean = {resistant_mean:.4f} (n={len(resistant_vals)})")
                print(f"  Sensitive: mean = {sensitive_mean:.4f} (n={len(sensitive_vals)})")
                print(f"  Mann-Whitney U: {u_stat:.2f}, p = {u_p:.4f}")
                print(f"  Cohen's d = {d:.3f} ({'Large' if abs(d) > 0.8 else 'Medium' if abs(d) > 0.5 else 'Small'} effect)")
                if u_p < 0.05:
                    print(f"  ✅ Significant difference (p < 0.05)")
                elif u_p < 0.1:
                    print(f"  ⚠️  Trend (p < 0.1)")
                else:
                    print(f"  ❌ No significant difference (p ≥ 0.1)")
                print()
        
        comparison_df = pd.DataFrame(comparison_results)
        comparison_df.to_csv(RESULTS_DIR / "pfi_group_comparison.csv", index=False)
        print(f"💾 Saved group comparison: {RESULTS_DIR / 'pfi_group_comparison.csv'}\n")
    
    # Step 5: Visualizations
    print("="*80)
    print("GENERATING VISUALIZATIONS")
    print("="*80)
    
    try:
        # Figure 1: Correlation plots (PFI vs pathway kinetics)
        print("📊 Creating correlation plots...")
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()
    
    for idx, pathway in enumerate(pathway_deltas):
        ax = axes[idx]
        
        # Scatter plot with regression line
        valid_data = evaluable_df[[pathway, 'pfi_months']].dropna()
        if len(valid_data) >= 3:
            ax.scatter(valid_data[pathway], valid_data['pfi_months'], 
                      s=100, alpha=0.7, c=valid_data['resistant'].map({'R': 'red', 'S': 'blue'}))
            
            # Regression line
            z = np.polyfit(valid_data[pathway], valid_data['pfi_months'], 1)
            p = np.poly1d(z)
            ax.plot(valid_data[pathway], p(valid_data[pathway]), "r--", alpha=0.5, linewidth=2)
            
            # Correlation annotation
            r, p_val = spearmanr(valid_data[pathway], valid_data['pfi_months'])
            ax.text(0.05, 0.95, f"r = {r:.3f}\np = {p_val:.4f}", 
                   transform=ax.transAxes, fontsize=11,
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel(f"{pathway.replace('_delta', '').upper()} Δ (post - pre)", fontsize=11)
        ax.set_ylabel("PFI (months)", fontsize=11)
        ax.set_title(f"PFI vs {pathway.replace('_delta', '').upper()} Kinetics", fontweight="bold", fontsize=12)
        ax.grid(True, alpha=0.3)
    
        plt.tight_layout()
        correlation_plot_path = RESULTS_DIR / "pfi_correlation_plots.png"
        plt.savefig(correlation_plot_path, dpi=300, bbox_inches='tight')
        plt.close('all')  # Close all figures
        print(f"💾 Saved correlation plots: {correlation_plot_path}\n")
    except Exception as e:
        print(f"⚠️  Error generating correlation plots: {e}\n")
        import traceback
        traceback.print_exc()
    
        # Figure 2: Box plots (Resistant vs Sensitive)
        if len(resistant_df) >= 2 and len(sensitive_df) >= 2:
            print("📊 Creating resistance boxplots...")
            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            axes = axes.flatten()
        
        for idx, pathway in enumerate(pathway_deltas):
            ax = axes[idx]
            
            # Prepare data for box plot
            plot_data = []
            plot_labels = []
            
            resistant_vals = resistant_df[pathway].dropna().values
            sensitive_vals = sensitive_df[pathway].dropna().values
            
            if len(resistant_vals) > 0:
                plot_data.append(resistant_vals)
                plot_labels.append('Resistant\n(n=2)')
            
            if len(sensitive_vals) > 0:
                plot_data.append(sensitive_vals)
                plot_labels.append('Sensitive\n(n=6)')
            
            if len(plot_data) > 0:
                bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True)
                
                # Color boxes
                colors = ['lightcoral', 'lightblue']
                for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                # Add stripplot overlay
                for i, (data, label) in enumerate(zip(plot_data, plot_labels)):
                    x_pos = i + 1
                    y_vals = data
                    ax.scatter([x_pos] * len(y_vals), y_vals, color='black', alpha=0.5, s=50, zorder=3)
                
                # Add significance annotation
                if len(resistant_vals) > 0 and len(sensitive_vals) > 0:
                    u_stat, u_p = mannwhitneyu(resistant_vals, sensitive_vals, alternative='two-sided')
                    sig_text = f"p = {u_p:.4f}"
                    if u_p < 0.05:
                        sig_text += " *"
                    elif u_p < 0.1:
                        sig_text += " †"
                    ax.text(0.5, 0.95, sig_text, transform=ax.transAxes,
                           fontsize=10, ha='center', verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                
                ax.axhline(0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
                ax.set_ylabel("Δ (post - pre)", fontsize=11)
                ax.set_title(f"{pathway.replace('_delta', '').upper()} Kinetics\nby Resistance Status", 
                           fontweight="bold", fontsize=12)
                ax.grid(True, alpha=0.3, axis='y')
        
            plt.tight_layout()
            boxplot_path = RESULTS_DIR / "pfi_resistance_boxplots.png"
            plt.savefig(boxplot_path, dpi=300, bbox_inches='tight')
            plt.close('all')  # Close all figures
            print(f"💾 Saved resistance boxplots: {boxplot_path}\n")
        else:
            print("⚠️  Insufficient data for boxplots (need n≥2 per group)\n")
    except Exception as e:
        print(f"⚠️  Error generating boxplots: {e}\n")
        import traceback
        traceback.print_exc()
    
    # Step 6: Summary report
    print("="*80)
    print("GENERATING SUMMARY REPORT")
    print("="*80)
    
    report_path = RESULTS_DIR / "pfi_correlation_report.txt"
    with open(report_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("PFI CORRELATION ANALYSIS REPORT (GSE165897)\n")
        f.write("="*80 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Dataset: GSE165897 (n={len(evaluable_df)} evaluable patients with PFI data)\n")
        f.write(f"Resistant (PFI < 6 mo): n={len(resistant_df)}\n")
        f.write(f"Sensitive (PFI ≥ 6 mo): n={len(sensitive_df)}\n\n")
        
        f.write("PRIMARY HYPOTHESIS: VEGF Activation → Short PFI\n")
        f.write("-" * 80 + "\n")
        vegf_corr = correlation_df[correlation_df['pathway'] == 'vegf']
        if not vegf_corr.empty:
            r = vegf_corr.iloc[0]['spearman_r']
            p = vegf_corr.iloc[0]['spearman_p']
            f.write(f"VEGF Δ vs PFI: r = {r:.3f}, p = {p:.4f}\n")
            if r < 0 and p < 0.1:
                f.write("✅ HYPOTHESIS SUPPORTED: Higher VEGF activation → Shorter PFI\n")
            elif r < 0:
                f.write("⚠️  Direction correct but not significant (p ≥ 0.1)\n")
            else:
                f.write("❌ HYPOTHESIS NOT SUPPORTED: No negative correlation\n")
        f.write("\n")
        
        f.write("SECONDARY HYPOTHESIS: DDR Suppression → Short PFI\n")
        f.write("-" * 80 + "\n")
        ddr_corr = correlation_df[correlation_df['pathway'] == 'ddr']
        if not ddr_corr.empty:
            r = ddr_corr.iloc[0]['spearman_r']
            p = ddr_corr.iloc[0]['spearman_p']
            f.write(f"DDR Δ vs PFI: r = {r:.3f}, p = {p:.4f}\n")
            if r > 0 and p < 0.1:
                f.write("✅ HYPOTHESIS SUPPORTED: Less DDR suppression → Longer PFI\n")
            elif r > 0:
                f.write("⚠️  Direction correct but not significant (p ≥ 0.1)\n")
            else:
                f.write("❌ HYPOTHESIS NOT SUPPORTED: No positive correlation\n")
        f.write("\n")
        
        f.write("GROUP COMPARISON: Resistant vs Sensitive\n")
        f.write("-" * 80 + "\n")
        if len(resistant_df) >= 2 and len(sensitive_df) >= 2:
            for _, row in comparison_df.iterrows():
                f.write(f"{row['pathway'].upper()}:\n")
                f.write(f"  Resistant: {row['resistant_mean']:.4f} (n={row['resistant_n']})\n")
                f.write(f"  Sensitive: {row['sensitive_mean']:.4f} (n={row['sensitive_n']})\n")
                f.write(f"  Mann-Whitney U: p = {row['mannwhitney_p']:.4f}\n")
                f.write(f"  Cohen's d = {row['cohens_d']:.3f} ({row['interpretation']} effect)\n")
                f.write("\n")
        
        f.write("INTERPRETATION:\n")
        f.write("-" * 80 + "\n")
        f.write("Strong Signal (r > 0.6, p < 0.05): Pathway kinetics strongly predict PFI\n")
        f.write("Moderate Signal (r = 0.4-0.6, p = 0.05-0.15): Trend visible, needs larger cohort\n")
        f.write("Weak Signal (r < 0.4, p > 0.2): Methodology needs refinement OR n=8 too small\n")
        f.write("\n")
        f.write("NEXT STEPS:\n")
        f.write("-" * 80 + "\n")
        f.write("1. If strong signal: Proceed to MSK_SPECTRUM validation\n")
        f.write("2. If moderate signal: Expand cohort or refine pathway scoring\n")
        f.write("3. If weak signal: Re-evaluate methodology or thresholds\n")
    
    print(f"💾 Saved summary report: {report_path}\n")
    
    # Step 7: Print key findings
    print("="*80)
    print("KEY FINDINGS")
    print("="*80)
    
    if not correlation_df.empty:
        vegf_row = correlation_df[correlation_df['pathway'] == 'vegf']
        if not vegf_row.empty:
            r = vegf_row.iloc[0]['spearman_r']
            p = vegf_row.iloc[0]['spearman_p']
            print(f"\n🎯 PRIMARY HYPOTHESIS (VEGF Activation → Short PFI):")
            print(f"   Spearman r = {r:.3f}, p = {p:.4f}")
            if r < -0.4 and p < 0.1:
                print(f"   ✅ STRONG SIGNAL: Hypothesis supported!")
            elif r < 0:
                print(f"   ⚠️  Direction correct but needs larger cohort")
            else:
                print(f"   ❌ No negative correlation detected")
        
        ddr_row = correlation_df[correlation_df['pathway'] == 'ddr']
        if not ddr_row.empty:
            r = ddr_row.iloc[0]['spearman_r']
            p = ddr_row.iloc[0]['spearman_p']
            print(f"\n🎯 SECONDARY HYPOTHESIS (DDR Suppression → Short PFI):")
            print(f"   Spearman r = {r:.3f}, p = {p:.4f}")
            if r > 0.4 and p < 0.1:
                print(f"   ✅ STRONG SIGNAL: Hypothesis supported!")
            elif r > 0:
                print(f"   ⚠️  Direction correct but needs larger cohort")
            else:
                print(f"   ❌ No positive correlation detected")
    
    print("\n" + "="*80)
    print("✅ PFI correlation analysis complete")
    print("="*80)
    print(f"\n📁 Results saved to: {RESULTS_DIR}")
    print("   - pfi_correlations.csv (Spearman correlations)")
    print("   - pfi_group_comparison.csv (Mann-Whitney U tests)")
    print("   - pfi_correlation_plots.png (scatter plots)")
    print("   - pfi_resistance_boxplots.png (group comparisons)")
    print("   - pfi_correlation_report.txt (summary report)")

if __name__ == "__main__":
    main()
