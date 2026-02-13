
import pandas as pd
import numpy as np
import os
import sys

# Define PARP Drugs
PARP_DRUGS = {"Olaparib", "Niraparib", "Talazoparib", "Rucaparib", "Veliparib"}

def main():
    print("=== MBD4 Logic Verification ===")
    
    # Paths (relative to repo root)
    omics_path = "publications/synthetic_lethality/data/OmicsSomaticMutations.csv"
    gdsc_path = "publications/synthetic_lethality/data/GDSC2_fitted_dose_response_27Oct23.xlsx"
    model_path = "data/depmap/Model.csv"
    
    # 1. Load MBD4 Mutations
    print(f"Loading MBD4 mutations from {omics_path}...")
    # Read only necessary columns to speed up
    # Guessing columns based on standard DepMap files: ModelID, HugoSymbol, ProteinChange, VariantInfo
    # We'll read headers first to be sure, but let's assume standard DepMap
    try:
        # Load chunks to find MBD4 rows efficiently
        mbd4_mutations = []
        chunk_size = 100000
        for chunk in pd.read_csv(omics_path, chunksize=chunk_size, low_memory=False):
            # Check column names in first chunk
            if "HugoSymbol" in chunk.columns:
                gene_col = "HugoSymbol"
            elif "Gene" in chunk.columns:
                gene_col = "Gene"
            else:
                print(f"Error: Could not find gene column in {chunk.columns}")
                return

            mbd4_chunk = chunk[chunk[gene_col] == "MBD4"]
            if not mbd4_chunk.empty:
                mbd4_mutations.append(mbd4_chunk)
        
        if not mbd4_mutations:
            print("No MBD4 mutations found in Omics file.")
            return

        mbd4_df = pd.concat(mbd4_mutations)
        print(f"Found {len(mbd4_df)} MBD4 mutations.")
    except Exception as e:
        print(f"Error loading Omics: {e}")
        return

    # 2. Load Model Metadata (to map ModelID to COSMIC_ID for GDSC)
    print(f"Loading Model metadata from {model_path}...")
    model_df = pd.read_csv(model_path, usecols=["ModelID", "COSMICID", "OncotreeLineage"])
    model_df = model_df.dropna(subset=["COSMICID"])
    # Normalize COSMIC ID
    def normalize_cosmic(x):
        try:
            return str(int(float(x)))
        except:
            return None
    model_df["COSMIC_ID_str"] = model_df["COSMICID"].apply(normalize_cosmic)
    
    # 3. Load GDSC Response
    print(f"Loading GDSC2 response from {gdsc_path}...")
    gdsc_df = pd.read_excel(gdsc_path, engine="openpyxl")
    gdsc_df["COSMIC_ID_str"] = gdsc_df["COSMIC_ID"].apply(normalize_cosmic)
    
    # Filter for PARP drugs
    parp_df = gdsc_df[gdsc_df["DRUG_NAME"].isin(PARP_DRUGS)]
    
    # 4. Integrate and Analyze
    print("\nAnayzing MBD4-mutated Cell Lines...")
    
    # Get MBD4 ModelIDs
    mbd4_models = mbd4_df["ModelID"].unique()
    print(f"Unique Cell Lines with MBD4 mutations: {len(mbd4_models)}")
    
    results = []
    
    for model_id in mbd4_models:
        # Get mutation details
        muts = mbd4_df[mbd4_df["ModelID"] == model_id]
        mut_desc = "; ".join([f"{r.get('ProteinChange', '')} ({r.get('VariantInfo', '')})" for _, r in muts.iterrows()])
        
        # Get COSMIC ID
        model_meta = model_df[model_df["ModelID"] == model_id]
        if model_meta.empty:
            continue
        cosmic_id = model_meta.iloc[0]["COSMIC_ID_str"]
        lineage = model_meta.iloc[0]["OncotreeLineage"]
        
        # Get PARP Response
        resp = parp_df[parp_df["COSMIC_ID_str"] == cosmic_id]
        if resp.empty:
            avg_z = np.nan
            status = "No Data"
        else:
            avg_z = resp["Z_SCORE"].mean()
            # Threshold: < -0.8 is Sensitive, > -0.8 is Resistant
            status = "SENSITIVE" if avg_z < -0.8 else "RESISTANT"
            
        prediction = "PARP (DDR=1.0)" # Based on the hardcoded logic we saw
        
        results.append({
            "ModelID": model_id,
            "Lineage": lineage,
            "Mutations": mut_desc,
            "Avg_PARP_Z": avg_z,
            "Actual_Status": status,
            "Predicted": prediction,
            "Verdict": "FALSE POSITIVE" if status == "RESISTANT" else "TRUE POSITIVE (Maybe)"
        })
        
    # 5. Output Results
    df_res = pd.DataFrame(results).dropna(subset=["Avg_PARP_Z"]) # Only show where we have data
    df_res = df_res.sort_values("Avg_PARP_Z")
    
    print("\n=== VALIDATION RESULTS ===")
    print(df_res[["ModelID", "Lineage", "Mutations", "Actual_Status", "Verdict"]].to_string(index=False))
    
    print(f"\nSummary:")
    print(f"Total Evaluated: {len(df_res)}")
    print(f"False Positives (Resistant but Predicted PARP): {len(df_res[df_res['Actual_Status'] == 'RESISTANT'])}")
    print(f"True Positives (Sensitive and Predicted PARP): {len(df_res[df_res['Actual_Status'] == 'SENSITIVE'])}")

if __name__ == "__main__":
    main()
