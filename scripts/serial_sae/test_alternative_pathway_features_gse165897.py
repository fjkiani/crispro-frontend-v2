"""
Test alternative pathway features on GSE165897 to find PFI correlation.

Tests 5 feature sets sequentially:
1. Post-treatment scores only
2. Pathway ratios
3. Composite multi-pathway score
4. Pre-treatment scores only
5. Delta magnitude (not direction)
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from pathlib import Path
import json
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
DATA_DIR = Path("data/serial_sae/gse165897")
RESULTS_DIR = DATA_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Success criteria
STRONG_SIGNAL = {"r": 0.5, "p": 0.1}  # Immediate report
MODERATE_SIGNAL = {"r": 0.4, "p": 0.15}  # Continue testing

def load_data():
    """Load pathway scores and PFI data."""
    pathway_scores_path = DATA_DIR / "results/pathway_scores.csv"
    pfi_data_path = DATA_DIR / "pfi_data_complete.json"
    
    # Load pathway scores (index is patient_timepoint, columns are pathways)
    pathway_df = pd.read_csv(pathway_scores_path, index_col=0)
    
    # Load PFI data
    if pfi_data_path.exists():
        with open(pfi_data_path) as f:
            patient_pfi = json.load(f)
    else:
        # Fallback: try resistance_labels.json
        resistance_labels_path = DATA_DIR / "resistance_labels.json"
        if resistance_labels_path.exists():
            with open(resistance_labels_path) as f:
                pfi_data = json.load(f)
            # Extract from _known_values or direct entries
            patient_pfi = {}
            if "_known_values" in pfi_data:
                for pid, data in pfi_data["_known_values"].items():
                    if isinstance(data, dict) and "pfi_days" in data:
                        patient_pfi[pid] = data["pfi_days"]
            # Also check direct entries
            for pid, value in pfi_data.items():
                if pid.startswith("EOC") and isinstance(value, dict) and "pfi_days" in value:
                    patient_pfi[pid] = value["pfi_days"]
    
    print(f"Loaded PFI data: {len(patient_pfi)} patients")
    
    # Pathway names (from columns, excluding metadata)
    pathway_names = [col for col in pathway_df.columns if col not in ["patient_id", "treatment_phase"]]
    print(f"Pathway names: {pathway_names}")
    
    # Create pre/post pathway scores per patient
    # pathway_df index format: "EOC1005_post-NACT", "EOC1005_treatment-naive"
    patients_data = {}
    
    for patient_id, pfi_days in patient_pfi.items():
        # Find pre and post samples for this patient
        pre_rows = pathway_df[
            (pathway_df.index.str.contains(patient_id)) & 
            (pathway_df.get("treatment_phase", pd.Series()) == "treatment-naive")
        ]
        post_rows = pathway_df[
            (pathway_df.index.str.contains(patient_id)) & 
            (pathway_df.get("treatment_phase", pd.Series()) == "post-NACT")
        ]
        
        if len(pre_rows) > 0 and len(post_rows) > 0:
            pre_scores = pre_rows.iloc[0][pathway_names]
            post_scores = post_rows.iloc[0][pathway_names]
            
            patients_data[patient_id] = {
                "pfi_days": pfi_days,
                "pre_scores": pre_scores.to_dict(),
                "post_scores": post_scores.to_dict(),
                "delta_scores": (post_scores - pre_scores).to_dict()
            }
    
    print(f"Patients with paired samples + PFI: {len(patients_data)}")
    return patients_data, pathway_df

def test_feature_set(features_dict, pfi_values, feature_name):
    """
    Test a set of features against PFI.
    
    Returns:
        results_df: DataFrame with r, p, feature_name
        top_features: List of top 3 features by |r|
    """
    results = []
    
    for feature_name_key, feature_values in features_dict.items():
        # Filter to patients with both feature and PFI
        valid_indices = [
            i for i, (f, p) in enumerate(zip(feature_values, pfi_values))
            if not (np.isnan(f) or np.isnan(p))
        ]
        
        if len(valid_indices) < 3:
            continue
        
        valid_features = [feature_values[i] for i in valid_indices]
        valid_pfi = [pfi_values[i] for i in valid_indices]
        
        # Pearson correlation
        r, p = pearsonr(valid_features, valid_pfi)
        
        # Spearman (rank-based, more robust)
        rho, p_rho = spearmanr(valid_features, valid_pfi)
        
        results.append({
            "feature": feature_name_key,
            "n": len(valid_indices),
            "r_pearson": r,
            "p_pearson": p,
            "r_spearman": rho,
            "p_spearman": p_rho,
            "abs_r": abs(r),
            "abs_rho": abs(rho)
        })
    
    results_df = pd.DataFrame(results)
    if len(results_df) == 0:
        return None, []
    
    # Sort by absolute correlation
    results_df = results_df.sort_values("abs_r", ascending=False)
    
    # Get top 3
    top_features = results_df.head(3)["feature"].tolist()
    
    return results_df, top_features

def test_1_post_treatment_scores(patients_data):
    """Test 1: Post-treatment scores only."""
    print("\n" + "="*80)
    print("TEST 1: POST-TREATMENT SCORES ONLY")
    print("="*80)
    
    # Get pathway names from first patient (lowercase: ddr, mapk, pi3k, vegf)
    if len(patients_data) == 0:
        return False, None, [], {}, []
    
    first_patient = list(patients_data.values())[0]
    pathway_names = [k for k in first_patient["post_scores"].keys() if k not in ["patient_id", "treatment_phase"]]
    print(f"Pathway names found: {pathway_names}")
    
    features_dict = {}
    pfi_values = []
    patient_ids = []
    
    for patient_id, data in patients_data.items():
        pfi_values.append(data["pfi_days"])
        patient_ids.append(patient_id)
        
        for pathway in pathway_names:
            feature_key = f"post_{pathway}"
            if feature_key not in features_dict:
                features_dict[feature_key] = []
            
            post_score = data["post_scores"].get(pathway, np.nan)
            features_dict[feature_key].append(post_score)
    
    results_df, top_features = test_feature_set(features_dict, pfi_values, "Post-Treatment")
    
    if results_df is not None:
        print("\nTop correlations:")
        print(results_df.head(10).to_string(index=False))
        
        # Check for strong signal
        strong = results_df[
            (results_df["abs_r"] >= STRONG_SIGNAL["r"]) & 
            (results_df["p_pearson"] <= STRONG_SIGNAL["p"])
        ]
        
        if len(strong) > 0:
            print(f"\n🚨 STRONG SIGNAL FOUND: {len(strong)} features")
            print(strong[["feature", "r_pearson", "p_pearson"]].to_string(index=False))
            return True, results_df, top_features, features_dict, pfi_values
    
    return False, results_df, top_features, features_dict, pfi_values

def test_2_pathway_ratios(patients_data):
    """Test 2: Pathway ratios."""
    print("\n" + "="*80)
    print("TEST 2: PATHWAY RATIOS")
    print("="*80)
    
    if len(patients_data) == 0:
        return False, None, [], {}, []
    
    # Get pathway names (lowercase)
    first_patient = list(patients_data.values())[0]
    pathway_map = {k.lower(): k for k in first_patient["post_scores"].keys() if k not in ["patient_id", "treatment_phase"]}
    
    features_dict = {}
    pfi_values = []
    
    for patient_id, data in patients_data.items():
        pfi_values.append(data["pfi_days"])
        
        post = data["post_scores"]
        
        # VEGF/DDR (use lowercase keys)
        vegf_key = pathway_map.get("vegf", "vegf")
        ddr_key = pathway_map.get("ddr", "ddr")
        if post.get(vegf_key, 0) > 0 and post.get(ddr_key, 0) > 0:
            features_dict.setdefault("VEGF_DDR_ratio", []).append(
                post[vegf_key] / post[ddr_key]
            )
        else:
            features_dict.setdefault("VEGF_DDR_ratio", []).append(np.nan)
        
        # PI3K/MAPK
        pi3k_key = pathway_map.get("pi3k", "pi3k")
        mapk_key = pathway_map.get("mapk", "mapk")
        if post.get(pi3k_key, 0) > 0 and post.get(mapk_key, 0) > 0:
            features_dict.setdefault("PI3K_MAPK_ratio", []).append(
                post[pi3k_key] / post[mapk_key]
            )
        else:
            features_dict.setdefault("PI3K_MAPK_ratio", []).append(np.nan)
        
        # (VEGF + Efflux) / (DDR + PI3K) - Efflux might not exist, use 0
        efflux_key = pathway_map.get("efflux", None)
        numerator = post.get(vegf_key, 0) + (post.get(efflux_key, 0) if efflux_key else 0)
        denominator = post.get(ddr_key, 0) + post.get(pi3k_key, 0)
        if denominator > 0:
            features_dict.setdefault("resistance_sensitivity_ratio", []).append(
                numerator / denominator
            )
        else:
            features_dict.setdefault("resistance_sensitivity_ratio", []).append(np.nan)
    
    results_df, top_features = test_feature_set(features_dict, pfi_values, "Ratios")
    
    if results_df is not None:
        print("\nTop correlations:")
        print(results_df.to_string(index=False))
        
        strong = results_df[
            (results_df["abs_r"] >= STRONG_SIGNAL["r"]) & 
            (results_df["p_pearson"] <= STRONG_SIGNAL["p"])
        ]
        
        if len(strong) > 0:
            print(f"\n🚨 STRONG SIGNAL FOUND: {len(strong)} features")
            return True, results_df, top_features, features_dict, pfi_values
    
    return False, results_df, top_features, features_dict, pfi_values

def test_3_composite_score(patients_data):
    """Test 3: Composite multi-pathway score."""
    print("\n" + "="*80)
    print("TEST 3: COMPOSITE MULTI-PATHWAY SCORE")
    print("="*80)
    
    if len(patients_data) == 0:
        return False, None, [], {}, []
    
    # Get pathway names (lowercase)
    first_patient = list(patients_data.values())[0]
    pathway_map = {k.lower(): k for k in first_patient["post_scores"].keys() if k not in ["patient_id", "treatment_phase"]}
    
    # Test different weight combinations (use lowercase keys)
    weight_sets = {
        "equal": {"vegf": 1, "efflux": 0, "her2": 0, "ddr": -1, "pi3k": -1, "mapk": 0},
        "tcga_informed": {"vegf": 1.5, "efflux": 0, "her2": 0, "ddr": -1.5, "pi3k": -1.2, "mapk": 0},
        "resistance_focused": {"vegf": 2.0, "efflux": 0, "her2": 0, "ddr": -2.0, "pi3k": -1.5, "mapk": -0.5},
    }
    
    features_dict = {}
    pfi_values = []
    
    for patient_id, data in patients_data.items():
        pfi_values.append(data["pfi_days"])
        post = data["post_scores"]
        
        for weight_name, weights in weight_sets.items():
            composite = 0.0
            for pathway_lower, weight in weights.items():
                pathway_key = pathway_map.get(pathway_lower)
                if pathway_key and pathway_key in post:
                    composite += weight * post[pathway_key]
            
            features_dict.setdefault(f"composite_{weight_name}", []).append(composite)
    
    results_df, top_features = test_feature_set(features_dict, pfi_values, "Composite")
    
    if results_df is not None:
        print("\nTop correlations:")
        print(results_df.to_string(index=False))
        
        strong = results_df[
            (results_df["abs_r"] >= STRONG_SIGNAL["r"]) & 
            (results_df["p_pearson"] <= STRONG_SIGNAL["p"])
        ]
        
        if len(strong) > 0:
            print(f"\n🚨 STRONG SIGNAL FOUND: {len(strong)} features")
            return True, results_df, top_features, features_dict, pfi_values
    
    return False, results_df, top_features, features_dict, pfi_values

def test_4_pre_treatment_scores(patients_data):
    """Test 4: Pre-treatment scores only."""
    print("\n" + "="*80)
    print("TEST 4: PRE-TREATMENT SCORES ONLY")
    print("="*80)
    
    if len(patients_data) == 0:
        return False, None, [], {}, []
    
    # Get pathway names from first patient (lowercase)
    first_patient = list(patients_data.values())[0]
    pathway_names = [k for k in first_patient["pre_scores"].keys() if k not in ["patient_id", "treatment_phase"]]
    
    features_dict = {}
    pfi_values = []
    
    for patient_id, data in patients_data.items():
        pfi_values.append(data["pfi_days"])
        
        for pathway in pathway_names:
            feature_key = f"pre_{pathway}"
            if feature_key not in features_dict:
                features_dict[feature_key] = []
            
            pre_score = data["pre_scores"].get(pathway, np.nan)
            features_dict[feature_key].append(pre_score)
    
    results_df, top_features = test_feature_set(features_dict, pfi_values, "Pre-Treatment")
    
    if results_df is not None:
        print("\nTop correlations:")
        print(results_df.head(10).to_string(index=False))
        
        strong = results_df[
            (results_df["abs_r"] >= STRONG_SIGNAL["r"]) & 
            (results_df["p_pearson"] <= STRONG_SIGNAL["p"])
        ]
        
        if len(strong) > 0:
            print(f"\n🚨 STRONG SIGNAL FOUND: {len(strong)} features")
            return True, results_df, top_features, features_dict, pfi_values
    
    return False, results_df, top_features, features_dict, pfi_values

def test_5_delta_magnitude(patients_data):
    """Test 5: Delta magnitude (not direction)."""
    print("\n" + "="*80)
    print("TEST 5: DELTA MAGNITUDE (NOT DIRECTION)")
    print("="*80)
    
    if len(patients_data) == 0:
        return False, None, [], {}, []
    
    # Get pathway names from first patient (lowercase)
    first_patient = list(patients_data.values())[0]
    pathway_names = [k for k in first_patient["delta_scores"].keys() if k not in ["patient_id", "treatment_phase"]]
    
    features_dict = {}
    pfi_values = []
    
    for patient_id, data in patients_data.items():
        pfi_values.append(data["pfi_days"])
        
        for pathway in pathway_names:
            feature_key = f"abs_delta_{pathway}"
            if feature_key not in features_dict:
                features_dict[feature_key] = []
            
            delta = data["delta_scores"].get(pathway, np.nan)
            features_dict[feature_key].append(abs(delta) if not np.isnan(delta) else np.nan)
    
    results_df, top_features = test_feature_set(features_dict, pfi_values, "Delta Magnitude")
    
    if results_df is not None:
        print("\nTop correlations:")
        print(results_df.head(10).to_string(index=False))
        
        strong = results_df[
            (results_df["abs_r"] >= STRONG_SIGNAL["r"]) & 
            (results_df["p_pearson"] <= STRONG_SIGNAL["p"])
        ]
        
        if len(strong) > 0:
            print(f"\n🚨 STRONG SIGNAL FOUND: {len(strong)} features")
            return True, results_df, top_features, features_dict, pfi_values
    
    return False, results_df, top_features, features_dict, pfi_values

def plot_top_features(features_dict, pfi_values, top_features, test_name, results_dir):
    """Create scatter plots for top 3 features."""
    if len(top_features) == 0:
        return
    
    fig, axes = plt.subplots(1, min(3, len(top_features)), figsize=(15, 5))
    if len(top_features) == 1:
        axes = [axes]
    
    for idx, feature_name in enumerate(top_features[:3]):
        ax = axes[idx] if len(top_features) > 1 else axes[0]
        
        # Get feature values
        feature_values = features_dict[feature_name]
        
        # Filter to valid pairs
        valid_pairs = [
            (f, p) for f, p in zip(feature_values, pfi_values)
            if not (np.isnan(f) or np.isnan(p))
        ]
        
        if len(valid_pairs) == 0:
            continue
        
        x_vals, y_vals = zip(*valid_pairs)
        
        # Scatter plot
        ax.scatter(x_vals, y_vals, alpha=0.7, s=100)
        
        # Fit line
        z = np.polyfit(x_vals, y_vals, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(x_vals), max(x_vals), 100)
        ax.plot(x_line, p(x_line), "r--", alpha=0.8, label=f"r={pearsonr(x_vals, y_vals)[0]:.3f}")
        
        ax.set_xlabel(feature_name)
        ax.set_ylabel("PFI (days)")
        ax.set_title(f"{feature_name} vs PFI")
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_path = results_dir / f"{test_name}_top_features.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Plots saved: {plot_path}")

def generate_report(test_name, results_df, top_features, features_dict, pfi_values, patients_data):
    """Generate full report with correlation table, plots, and GO/NO-GO decision."""
    print("\n" + "="*80)
    print("GENERATING FULL REPORT")
    print("="*80)
    
    # 1. Correlation table
    print("\n📊 CORRELATION TABLE (All Features vs PFI)")
    print("="*80)
    print(results_df.to_string(index=False))
    
    # Save correlation table
    table_path = RESULTS_DIR / f"{test_name}_correlation_table.csv"
    results_df.to_csv(table_path, index=False)
    print(f"\n✅ Correlation table saved: {table_path}")
    
    # 2. Scatter plots for top 3 features
    print("\n📈 GENERATING SCATTER PLOTS (Top 3 Features)")
    print("="*80)
    plot_top_features(features_dict, pfi_values, top_features, test_name, RESULTS_DIR)
    
    # 3. Statistical summary
    print("\n📋 STATISTICAL SUMMARY")
    print("="*80)
    top_3 = results_df.head(3)
    for _, row in top_3.iterrows():
        print(f"\n{row['feature']}:")
        print(f"  Pearson r: {row['r_pearson']:.4f} (p = {row['p_pearson']:.4f})")
        print(f"  Spearman ρ: {row['r_spearman']:.4f} (p = {row['p_spearman']:.4f})")
        print(f"  Sample size: n = {int(row['n'])}")
        print(f"  Effect size: |r| = {row['abs_r']:.4f}")
    
    # 4. GO/NO-GO Decision
    print("\n" + "="*80)
    print("GO/NO-GO DECISION FOR MSK_SPECTRUM")
    print("="*80)
    
    strong_signals = results_df[
        (results_df["abs_r"] >= STRONG_SIGNAL["r"]) & 
        (results_df["p_pearson"] <= STRONG_SIGNAL["p"])
    ]
    
    if len(strong_signals) > 0:
        print(f"\n✅ GO: STRONG SIGNAL FOUND ({len(strong_signals)} features)")
        print("\nTop signals:")
        for _, row in strong_signals.iterrows():
            print(f"  - {row['feature']}: r = {row['r_pearson']:.4f}, p = {row['p_pearson']:.4f}")
        print("\n🎯 RECOMMENDATION: PROCEED with MSK_SPECTRUM dataset")
        print("   Rationale: Strong correlation (r > 0.5, p < 0.1) indicates")
        print("   pathway features can predict platinum resistance.")
        decision = "GO"
    else:
        moderate_signals = results_df[
            (results_df["abs_r"] >= MODERATE_SIGNAL["r"]) & 
            (results_df["p_pearson"] <= MODERATE_SIGNAL["p"])
        ]
        if len(moderate_signals) > 0:
            print(f"\n⚠️  CONDITIONAL GO: MODERATE SIGNAL ({len(moderate_signals)} features)")
            print("\nTop signals:")
            for _, row in moderate_signals.iterrows():
                print(f"  - {row['feature']}: r = {row['r_pearson']:.4f}, p = {row['p_pearson']:.4f}")
            print("\n🎯 RECOMMENDATION: PROCEED with MSK_SPECTRUM (with caution)")
            print("   Rationale: Moderate correlation (r > 0.4, p < 0.15) suggests")
            print("   pathway features may predict resistance, but validation needed.")
            decision = "CONDITIONAL_GO"
        else:
            print("\n❌ NO-GO: NO SIGNIFICANT SIGNAL")
            print("\n🎯 RECOMMENDATION: RE-EVALUATE approach or test additional datasets")
            print("   Rationale: No features meet correlation thresholds (r > 0.4, p < 0.15)")
            decision = "NO_GO"
    
    # Save decision
    decision_path = RESULTS_DIR / f"{test_name}_decision.txt"
    with open(decision_path, "w") as f:
        f.write(f"DECISION: {decision}\n")
        f.write(f"Test: {test_name}\n")
        f.write(f"Top Feature: {top_features[0] if top_features else 'N/A'}\n")
        f.write(f"Correlation: {results_df.iloc[0]['r_pearson']:.4f} (p = {results_df.iloc[0]['p_pearson']:.4f})\n")
    print(f"\n✅ Decision saved: {decision_path}")
    
    return decision

def main():
    """Run all tests sequentially."""
    print("="*80)
    print("ALTERNATIVE PATHWAY FEATURES TEST - GSE165897")
    print("="*80)
    
    # Load data
    patients_data, pathway_df = load_data()
    
    print(f"\n✅ Loaded {len(patients_data)} patients with paired samples + PFI")
    
    all_results = []
    found_strong_signal = False
    
    # Test 1: Post-treatment scores
    stop, results_df, top_features, features_dict, pfi_values = test_1_post_treatment_scores(patients_data)
    if results_df is not None:
        all_results.append(("Post-Treatment", results_df))
        if stop:
            found_strong_signal = True
            print("\n🚨 STRONG SIGNAL FOUND in Test 1 - Generating full report...")
            decision = generate_report("Post-Treatment", results_df, top_features, features_dict, pfi_values, patients_data)
            return
    
    # Test 2: Pathway ratios
    stop, results_df, top_features = test_2_pathway_ratios(patients_data)
    if results_df is not None:
        all_results.append(("Ratios", results_df))
        if stop:
            found_strong_signal = True
            print("\n🚨 STOPPING: Strong signal found in Test 2")
            return
    
    # Test 3: Composite scores
    stop, results_df, top_features = test_3_composite_score(patients_data)
    if results_df is not None:
        all_results.append(("Composite", results_df))
        if stop:
            found_strong_signal = True
            print("\n🚨 STOPPING: Strong signal found in Test 3")
            return
    
    # Test 4: Pre-treatment scores
    stop, results_df, top_features = test_4_pre_treatment_scores(patients_data)
    if results_df is not None:
        all_results.append(("Pre-Treatment", results_df))
        if stop:
            found_strong_signal = True
            print("\n🚨 STOPPING: Strong signal found in Test 4")
            return
    
    # Test 5: Delta magnitude
    stop, results_df, top_features = test_5_delta_magnitude(patients_data)
    if results_df is not None:
        all_results.append(("Delta Magnitude", results_df))
        if stop:
            found_strong_signal = True
            print("\n🚨 STOPPING: Strong signal found in Test 5")
            return
    
    # Compile all results
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    
    combined_results = []
    for test_name, df in all_results:
        df_copy = df.copy()
        df_copy["test"] = test_name
        combined_results.append(df_copy)
    
    if len(combined_results) > 0:
        all_results_df = pd.concat(combined_results, ignore_index=True)
        all_results_df = all_results_df.sort_values("abs_r", ascending=False)
        
        print("\nTop 10 features across all tests:")
        print(all_results_df.head(10)[["test", "feature", "r_pearson", "p_pearson", "n"]].to_string(index=False))
        
        # Save results
        results_path = RESULTS_DIR / "alternative_features_correlation.csv"
        all_results_df.to_csv(results_path, index=False)
        print(f"\n✅ All results saved: {results_path}")
        
        # Check for moderate signals
        moderate = all_results_df[
            (all_results_df["abs_r"] >= MODERATE_SIGNAL["r"]) & 
            (all_results_df["p_pearson"] <= MODERATE_SIGNAL["p"])
        ]
        
        if len(moderate) > 0:
            print(f"\n✅ MODERATE SIGNAL: {len(moderate)} features")
            print(moderate[["test", "feature", "r_pearson", "p_pearson"]].to_string(index=False))
            print("\n✅ RECOMMENDATION: PROCEED with MSK_SPECTRUM (moderate signal found)")
        else:
            print("\n⚠️  NO SIGNIFICANT SIGNAL FOUND")
            print("   Recommendation: Re-evaluate approach or test additional datasets")
    else:
        print("\n❌ NO RESULTS: Check data loading")

if __name__ == "__main__":
    main()
