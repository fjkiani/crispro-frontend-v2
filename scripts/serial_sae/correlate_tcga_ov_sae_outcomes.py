#!/usr/bin/env python3
"""
Correlate TCGA-OV SAE Pathway Scores with Clinical Outcomes

Matches computed SAE scores with clinical outcomes (OS, PFS, platinum response)
and performs statistical correlation analysis.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import sys
from scipy import stats
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

# Add backend to path for clinical extraction
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

DATA_DIR = Path("data/serial_sae/tcga_ov")
RESULTS_DIR = Path("data/serial_sae/tcga_ov/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# GDC API for metadata
GDC_API_BASE = "https://api.gdc.cancer.gov"

def extract_tcga_patient_id_from_gdc_uuid(uuid: str) -> Optional[str]:
    """
    Extract TCGA patient ID from GDC file UUID.
    Uses GDC API query endpoint to get case submitter_id.
    """
    import httpx
    
    try:
        # Query GDC cases API for cases associated with this file
        query = {
            "filters": {
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "files.file_id",
                            "value": [uuid]
                        }
                    }
                ]
            },
            "fields": "submitter_id",
            "format": "JSON",
            "size": 1
        }
        
        url = f"{GDC_API_BASE}/cases"
        response = httpx.post(url, json=query, timeout=10.0)
        
        if response.status_code == 200:
            data = response.json()["data"]
            hits = data.get("hits", [])
            
            if hits:
                submitter_id = hits[0].get("submitter_id")
                
                # Extract TCGA patient ID (e.g., TCGA-XX-XXXX from TCGA-XX-XXXX-XX)
                if submitter_id:
                    # TCGA patient ID is first 3 segments (TCGA-XX-XXXX)
                    parts = submitter_id.split("-")
                    if len(parts) >= 3:
                        return f"{parts[0]}-{parts[1]}-{parts[2]}"
        
        return None
    except Exception as e:
        print(f"  ⚠️  Error extracting patient ID for {uuid}: {e}")
        return None

def load_existing_clinical_data() -> Optional[pd.DataFrame]:
    """
    Load existing TCGA-OV clinical data if available.
    Checks multiple potential locations.
    """
    # Check existing clinical report
    clinical_report_path = Path("publications/synthetic_lethality/results/clinical/tcga_ov_clinical_report.json")
    if clinical_report_path.exists():
        print("📥 Loading existing clinical report...")
        with open(clinical_report_path) as f:
            report = json.load(f)
        
        # Extract patient-level data if available
        # Note: This report has aggregate stats, not patient-level
        # We'll need to extract fresh from cBioPortal
    
    # Try to extract fresh from cBioPortal
    try:
        from biomarker_enriched_cohorts.scripts.extract_tcga_outcomes import extract_tcga_ov_outcomes
        
        print("📡 Extracting TCGA-OV clinical outcomes from cBioPortal...")
        cohort_artifact, receipt = extract_tcga_ov_outcomes()
        
        # Convert to DataFrame
        patients = cohort_artifact.get("cohort", {}).get("patients", [])
        if not patients:
            return None
        
        clinical_records = []
        for patient in patients:
            patient_id = patient.get("patient_id", "")
            outcomes = patient.get("outcomes", {})
            
            clinical_records.append({
                "patient_id": patient_id,
                "os_days": outcomes.get("os_days"),
                "os_event": outcomes.get("os_event"),
                "pfs_days": outcomes.get("pfs_days"),
                "pfs_event": outcomes.get("pfs_event")
            })
        
        clinical_df = pd.DataFrame(clinical_records)
        print(f"✅ Loaded {len(clinical_df)} patients with clinical data")
        return clinical_df
        
    except Exception as e:
        print(f"⚠️  Error extracting clinical data: {e}")
        return None

def batch_extract_patient_ids(uuids: List[str], batch_size: int = 100) -> Dict[str, Optional[str]]:
    """
    Batch extract TCGA patient IDs from GDC file UUIDs using GDC API.
    Queries files directly and expands cases to get submitter_id.
    """
    import httpx
    
    results = {}
    
    for i in range(0, len(uuids), batch_size):
        batch = uuids[i:i+batch_size]
        print(f"  Batch {i//batch_size + 1}/{(len(uuids)-1)//batch_size + 1}: Processing {len(batch)} files...")
        
        query = {
            "filters": {
                "op": "in",
                "content": {
                    "field": "file_id",
                    "value": batch
                }
            },
            "fields": "file_id,cases.submitter_id",
            "format": "JSON",
            "size": len(batch),
            "expand": "cases"
        }
        
        try:
            url = f"{GDC_API_BASE}/files"
            response = httpx.post(url, json=query, timeout=60.0)
            
            if response.status_code == 200:
                data = response.json()["data"]
                hits = data.get("hits", [])
                
                # Build mapping from file_id to submitter_id
                for hit in hits:
                    file_id = hit.get("file_id")
                    cases = hit.get("cases", [])
                    
                    # Get submitter_id from first case (files typically have one case)
                    if cases and len(cases) > 0:
                        submitter_id = cases[0].get("submitter_id")
                        
                        # Extract TCGA patient ID (first 3 segments)
                        if submitter_id:
                            parts = submitter_id.split("-")
                            if len(parts) >= 3:
                                patient_id = f"{parts[0]}-{parts[1]}-{parts[2]}"
                            else:
                                patient_id = None
                        else:
                            patient_id = None
                    else:
                        patient_id = None
                    
                    if file_id:
                        results[file_id] = patient_id
                
                # Mark unmatched UUIDs as None
                for uuid in batch:
                    if uuid not in results:
                        results[uuid] = None
            else:
                print(f"  ⚠️  Batch {i//batch_size + 1} API error: {response.status_code}")
                # Mark entire batch as None on error
                for uuid in batch:
                    results[uuid] = None
                    
        except Exception as e:
            print(f"  ⚠️  Error in batch {i//batch_size + 1}: {e}")
            # Mark entire batch as None on error
            for uuid in batch:
                results[uuid] = None
    
    return results

def match_samples_to_patients(sample_metadata_df: pd.DataFrame) -> pd.DataFrame:
    """
    Match GDC file IDs to TCGA patient IDs using batch GDC API queries.
    First maps filename UUIDs to GDC file_ids using the manifest.
    """
    print("🔗 Matching samples to TCGA patient IDs (batch processing)...")
    
    # Load manifest to map filename UUIDs to GDC file_ids
    manifest_path = DATA_DIR / "tcga_ov_rnaseq_manifest.csv"
    if not manifest_path.exists():
        print(f"  ⚠️  Manifest not found: {manifest_path}")
        print("  Using sample_id directly (may not work if it's not a GDC file_id)")
        manifest_df = None
    else:
        manifest_df = pd.read_csv(manifest_path)
        # Extract filename UUID (first part before first dot)
        manifest_df["filename_uuid"] = manifest_df["file_name"].apply(
            lambda x: x.split('.')[0] if '.' in str(x) else str(x)
        )
        # Create mapping: filename_uuid -> file_id
        filename_to_fileid = dict(zip(manifest_df["filename_uuid"], manifest_df["file_id"]))
        print(f"  ✅ Loaded manifest: {len(filename_to_fileid)} file mappings")
    
    # Map sample_ids (filename UUIDs) to GDC file_ids
    total_samples = len(sample_metadata_df)
    gdc_file_ids = []
    
    if manifest_df is not None:
        for sample_id in sample_metadata_df["sample_id"]:
            gdc_file_id = filename_to_fileid.get(sample_id, sample_id)  # Fallback to sample_id if not found
            gdc_file_ids.append(gdc_file_id)
    else:
        gdc_file_ids = sample_metadata_df["sample_id"].tolist()
    
    print(f"  Processing {total_samples} samples in batches of 100...")
    patient_id_map = batch_extract_patient_ids(gdc_file_ids, batch_size=100)
    
    # Build matched DataFrame
    matched_data = []
    for idx, row in sample_metadata_df.iterrows():
        sample_id = row["sample_id"]
        # Get GDC file_id (either from mapping or use sample_id directly)
        if manifest_df is not None:
            gdc_file_id = filename_to_fileid.get(sample_id, sample_id)
        else:
            gdc_file_id = sample_id
        
        # Get patient_id from map (keyed by GDC file_id)
        patient_id = patient_id_map.get(gdc_file_id)
        
        matched_data.append({
            "file_id": gdc_file_id,  # Use GDC file_id for merging with pathway scores
            "patient_id": patient_id,
            "sample_id": sample_id  # Keep original sample_id for reference
        })
    
    matched_df = pd.DataFrame(matched_data)
    matched_count = matched_df["patient_id"].notna().sum()
    print(f"✅ Matched {matched_count}/{total_samples} samples to patient IDs")
    
    return matched_df

def correlate_sae_with_outcomes(
    sae_scores: Dict[str, List[float]],
    pathway_scores_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    sample_to_patient_df: pd.DataFrame
) -> Dict:
    """
    Correlate SAE pathway scores with clinical outcomes.
    """
    print("\n📊 Correlating SAE scores with clinical outcomes...")
    
    # Convert patient_id to string in both DataFrames for consistent merging
    sample_to_patient_df["patient_id"] = sample_to_patient_df["patient_id"].astype(str)
    clinical_df["patient_id"] = clinical_df["patient_id"].astype(str)
    
    # Filter out NaN/None values before merging
    sample_to_patient_df = sample_to_patient_df[sample_to_patient_df["patient_id"] != "nan"]
    sample_to_patient_df = sample_to_patient_df[sample_to_patient_df["patient_id"] != "None"]
    
    # Merge sample IDs → patient IDs → clinical data
    merged_df = sample_to_patient_df.merge(
        clinical_df,
        on="patient_id",
        how="inner"
    )
    
    # Add pathway scores (pathway_scores_df is indexed by sample_id/filename UUID)
    pathway_cols = ["ddr", "ras_mapk", "pi3k", "vegf", "her2", "io", "efflux"]
    for col in pathway_cols:
        if col in pathway_scores_df.columns:
            # Map sample_id to pathway score (pathway scores are indexed by filename UUID)
            merged_df = merged_df.merge(
                pathway_scores_df[[col]].reset_index().rename(columns={"index": "sample_id"}),
                on="sample_id",
                how="left"
            )
    
    # Add mechanism vector magnitude (L2 norm)
    mechanism_vector_cols = [col for col in pathway_cols if col in merged_df.columns]
    if mechanism_vector_cols:
        merged_df["mechanism_vector_magnitude"] = np.sqrt(
            merged_df[mechanism_vector_cols].pow(2).sum(axis=1)
        )
    
    # Filter to patients with both SAE and clinical data
    analysis_df = merged_df[
        merged_df["os_days"].notna() & 
        merged_df["os_event"].notna() &
        merged_df["ddr"].notna()
    ].copy()
    
    print(f"✅ {len(analysis_df)} patients with complete SAE + clinical data")
    
    # Statistical analyses
    results = {
        "n_patients": len(analysis_df),
        "pathway_correlations": {},
        "survival_analyses": {}
    }
    
    # 1. Pathway correlations with OS/PFS
    for pathway in pathway_cols:
        if pathway not in analysis_df.columns:
            continue
        
        # OS correlation
        os_corr = analysis_df[[pathway, "os_days"]].corr().iloc[0, 1]
        os_p = stats.pearsonr(analysis_df[pathway], analysis_df["os_days"])[1]
        
        # PFS correlation (if available)
        pfs_corr = None
        pfs_p = None
        if "pfs_days" in analysis_df.columns:
            pfs_valid = analysis_df["pfs_days"].notna()
            if pfs_valid.sum() > 10:
                pfs_corr = analysis_df.loc[pfs_valid, [pathway, "pfs_days"]].corr().iloc[0, 1]
                pfs_p = stats.pearsonr(
                    analysis_df.loc[pfs_valid, pathway],
                    analysis_df.loc[pfs_valid, "pfs_days"]
                )[1]
        
        results["pathway_correlations"][pathway] = {
            "os_correlation": float(os_corr) if not np.isnan(os_corr) else None,
            "os_p_value": float(os_p) if not np.isnan(os_p) else None,
            "pfs_correlation": float(pfs_corr) if pfs_corr is not None and not np.isnan(pfs_corr) else None,
            "pfs_p_value": float(pfs_p) if pfs_p is not None and not np.isnan(pfs_p) else None
        }
    
    # 2. Survival analysis (high vs low DDR pathway)
    if "ddr" in analysis_df.columns:
        ddr_median = analysis_df["ddr"].median()
        analysis_df["ddr_high"] = analysis_df["ddr"] > ddr_median
        
        # OS by DDR
        os_high = analysis_df[analysis_df["ddr_high"]]["os_days"].median()
        os_low = analysis_df[~analysis_df["ddr_high"]]["os_days"].median()
        
        # Survival analysis using Mann-Whitney U test as proxy for log-rank
        # (Full log-rank test would require lifelines package)
        try:
            high_times = analysis_df[analysis_df["ddr_high"]]["os_days"].values
            high_events = analysis_df[analysis_df["ddr_high"]]["os_event"].astype(int).values
            low_times = analysis_df[~analysis_df["ddr_high"]]["os_days"].values
            low_events = analysis_df[~analysis_df["ddr_high"]]["os_event"].astype(int).values
            
            # Note: scipy.stats doesn't have logrank_test, would need lifelines
            # For now, use Mann-Whitney U test as proxy
            mw_stat, mw_p = stats.mannwhitneyu(high_times, low_times, alternative='two-sided')
            
            results["survival_analyses"]["ddr_high_vs_low"] = {
                "median_os_high": float(os_high),
                "median_os_low": float(os_low),
                "mann_whitney_p": float(mw_p),
                "n_high": int(analysis_df["ddr_high"].sum()),
                "n_low": int((~analysis_df["ddr_high"]).sum())
            }
        except Exception as e:
            print(f"  ⚠️  Error in survival analysis: {e}")
    
    # 3. Mechanism vector magnitude correlation
    if "mechanism_vector_magnitude" in analysis_df.columns:
        mag_os_corr = analysis_df[["mechanism_vector_magnitude", "os_days"]].corr().iloc[0, 1]
        mag_os_p = stats.pearsonr(
            analysis_df["mechanism_vector_magnitude"],
            analysis_df["os_days"]
        )[1]
        
        results["mechanism_vector_magnitude"] = {
            "os_correlation": float(mag_os_corr) if not np.isnan(mag_os_corr) else None,
            "os_p_value": float(mag_os_p) if not np.isnan(mag_os_p) else None
        }
    
    return results, analysis_df

def main():
    print("="*80)
    print("Correlate TCGA-OV SAE Scores with Clinical Outcomes")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Load SAE scores
    sae_scores_path = RESULTS_DIR / "tcga_ov_sae_scores.json"
    pathway_scores_path = RESULTS_DIR / "tcga_ov_pathway_scores.csv"
    sample_metadata_path = DATA_DIR / "processed" / "tcga_ov_sample_metadata.csv"
    
    if not sae_scores_path.exists():
        print(f"❌ Error: SAE scores not found at {sae_scores_path}")
        print("Please run scripts/serial_sae/compute_tcga_ov_sae.py first.")
        return
    
    print("📥 Loading SAE scores...")
    with open(sae_scores_path) as f:
        sae_scores = json.load(f)
    
    pathway_scores_df = pd.read_csv(pathway_scores_path, index_col=0)
    sample_metadata_df = pd.read_csv(sample_metadata_path)
    
    print(f"✅ Loaded: {len(sae_scores)} samples with SAE scores")
    print(f"✅ Pathway scores: {pathway_scores_df.shape[0]} samples × {pathway_scores_df.shape[1]} pathways")
    
    # Load clinical data
    clinical_df = load_existing_clinical_data()
    if clinical_df is None or clinical_df.empty:
        print("❌ Error: Could not load clinical outcomes")
        return
    
    # Match samples to patients
    # Check if we already have a mapping
    mapping_path = RESULTS_DIR / "tcga_ov_sample_to_patient_mapping.csv"
    if mapping_path.exists():
        print("📥 Loading existing sample-to-patient mapping...")
        sample_to_patient_df = pd.read_csv(mapping_path)
    else:
        print("🔗 Creating sample-to-patient mapping (this may take a while)...")
        sample_to_patient_df = match_samples_to_patients(sample_metadata_df)
        sample_to_patient_df.to_csv(mapping_path, index=False)
        print(f"💾 Saved mapping: {mapping_path}")
    
    # Correlate SAE with outcomes
    correlation_results, analysis_df = correlate_sae_with_outcomes(
        sae_scores,
        pathway_scores_df,
        clinical_df,
        sample_to_patient_df
    )
    
    # Save results
    results_path = RESULTS_DIR / "tcga_ov_sae_clinical_correlation.json"
    with open(results_path, "w") as f:
        json.dump(correlation_results, f, indent=2)
    print(f"💾 Saved correlation results: {results_path}")
    
    # Save analysis dataframe
    analysis_df_path = RESULTS_DIR / "tcga_ov_sae_clinical_merged.csv"
    analysis_df.to_csv(analysis_df_path, index=False)
    print(f"💾 Saved merged data: {analysis_df_path}")
    
    # Print summary
    print("\n" + "="*80)
    print("Correlation Summary")
    print("="*80)
    print(f"Patients analyzed: {correlation_results['n_patients']}")
    print("\nPathway Correlations with OS:")
    for pathway, corr_data in correlation_results.get("pathway_correlations", {}).items():
        os_corr = corr_data.get("os_correlation")
        os_p = corr_data.get("os_p_value")
        if os_corr is not None:
            sig = "***" if os_p < 0.001 else "**" if os_p < 0.01 else "*" if os_p < 0.05 else ""
            print(f"  {pathway:12s}: r={os_corr:6.3f}, p={os_p:.4f} {sig}")
    
    if "survival_analyses" in correlation_results:
        ddr_analysis = correlation_results["survival_analyses"].get("ddr_high_vs_low", {})
        if ddr_analysis:
            print(f"\nDDR High vs Low OS:")
            print(f"  High DDR median OS: {ddr_analysis.get('median_os_high', 0):.1f} days")
            print(f"  Low DDR median OS: {ddr_analysis.get('median_os_low', 0):.1f} days")
            print(f"  Mann-Whitney p: {ddr_analysis.get('mann_whitney_p', 0):.4f}")
    
    print("\n" + "="*80)
    print("✅ Correlation analysis complete")
    print("="*80)
    print("\nNext steps:")
    print("1. Generate figures (survival curves, correlation plots)")
    print("2. Validate pathway signatures against platinum response")
    print("3. Compare with published TCGA-OV pathway analyses")

if __name__ == "__main__":
    main()
