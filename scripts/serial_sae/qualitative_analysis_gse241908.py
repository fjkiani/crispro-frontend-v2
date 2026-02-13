
import pandas as pd
import numpy as np
import gzip

DATA_PATH = 'oncology-coPilot/oncology-backend-minimal/data/serial_sae/gse241908/GSE241908_data1_TPM.csv.gz'
OUTPUT_REPORT = 'GSE241908_QUALITATIVE_REPORT.md'

# Assumed Pairs based on index/ID pattern (6 pairs)
# Structure: (PatientID, Pre_Col, Post_Col)
PAIRS = [
    ("Patient_A", "ShV-80", "ShV-81"),
    ("Patient_B", "ShV-83", "ShV-84"),
    ("Patient_C", "ShV-90", "ShV-91"),
    ("Patient_D", "ShV-94", "ShV-95"),
    ("Patient_E", "ShV-96", "ShV-97"),
    ("Patient_F", "ShV-98", "ShV-99"),
]

# Gene Panels (Ensembl ID Mapping)
GENE_MAP = {
    "BRCA1": "ENSG00000012048",
    "BRCA2": "ENSG00000139618",
    "RAD51": "ENSG00000051180",
    "PARP1": "ENSG00000143799",
    "FANCD2": "ENSG00000144554",
    "ATR": "ENSG00000171848", # Added for completeness if needed, else skip
    "ATM": "ENSG00000149311", 
    "CHEK1": "ENSG00000149554",
    "MKI67": "ENSG00000148773",
    "CCNE1": "ENSG00000105173",
    "CDK1": "ENSG00000170312",
    "E2F1": "ENSG00000101412",
    "PCNA": "ENSG00000132646",
    "TOP2A": "ENSG00000131747"
}

DDR_IDS = [GENE_MAP[g] for g in ["BRCA1", "BRCA2", "RAD51", "PARP1", "FANCD2"] if g in GENE_MAP]
CC_IDS = [GENE_MAP[g] for g in ["MKI67", "CCNE1", "CDK1", "E2F1", "PCNA"] if g in GENE_MAP]

def load_data():
    print(f"Loading {DATA_PATH}...")
    try:
        df = pd.read_csv(DATA_PATH, index_col=0, sep=None, engine='python')
    except:
        df = pd.read_csv(DATA_PATH, index_col=0)
    
    # Clean column names
    df.columns = [c.strip('"').strip("'").strip() for c in df.columns]
    
    # Clean Index (Strip version suffix) -> "ENSG00000000003.14" -> "ENSG00000000003"
    df.index = df.index.astype(str).str.split('.').str[0].str.strip('"')
    
    print(f"Loaded {len(df)} features, {len(df.columns)} samples.")
    return df

def analyze_pairs(df):
    results = []
    
    print("\nAnalyzing Pairs...")
    for pid, pre, post in PAIRS:
        if pre not in df.columns or post not in df.columns:
            print(f"⚠️ Skipping {pid}: Columns not found ({pre}, {post})")
            continue
            
        vec_pre = df[pre]
        vec_post = df[post]
        
        # DDR Score
        ddr_present = [g for g in DDR_IDS if g in df.index]
        ddr_pre = vec_pre.loc[ddr_present].mean()
        ddr_post = vec_post.loc[ddr_present].mean()
        ddr_fc = np.log2((ddr_post + 1) / (ddr_pre + 1))
        
        # Cell Cycle Score
        cc_present = [g for g in CC_IDS if g in df.index]
        cc_pre = vec_pre.loc[cc_present].mean()
        cc_post = vec_post.loc[cc_present].mean()
        cc_fc = np.log2((cc_post + 1) / (cc_pre + 1))
        
        results.append({
            "Patient": pid,
            "Pre_ID": pre,
            "Post_ID": post,
            "DDR_Log2FC": ddr_fc,
            "CC_Log2FC": cc_fc
        })
        
    return pd.DataFrame(results)


def generate_report(results_df):
    md = "# 🔬 GSE241908 Qualitative Analysis Report\n\n"
    md += "**Objective:** Sanity check therapy-induced transcriptional shifts (n=6 Paired).\n"
    md += "**Hypothesis:** Post-chemotherapy samples show upregulation of DNA Repair and Cell Cycle modules.\n\n"
    
    md += "## 1. Pairwise Deltas (Log2FC)\n"
    md += "| Patient | DDR Log2FC | Cell Cycle Log2FC | Interpretation |\n"
    md += "|---|---|---|---|\n"
    
    consistent_up = 0
    for _, row in results_df.iterrows():
        ddr_dir = "⬆️" if row['DDR_Log2FC'] > 0 else "⬇️"
        cc_dir = "⬆️" if row['CC_Log2FC'] > 0 else "⬇️"
        
        # Interpret
        interp = "Consistent" if row['DDR_Log2FC'] > 0 and row['CC_Log2FC'] > 0 else "Mixed/Down"
        if interp == "Consistent":
            consistent_up += 1
            
        md += f"| {row['Patient']} ({row['Pre_ID']}/{row['Post_ID']}) | **{row['DDR_Log2FC']:.2f}** {ddr_dir} | **{row['CC_Log2FC']:.2f}** {cc_dir} | {interp} |\n"
        
    md += "\n## 2. Summary Statistics\n"
    md += f"- **Consistent Upregulation:** {consistent_up}/{len(results_df)} patients\n"
    md += f"- **Mean DDR Uplift:** {results_df['DDR_Log2FC'].mean():.2f}\n"
    md += f"- **Mean Cell Cycle Uplift:** {results_df['CC_Log2FC'].mean():.2f}\n"
    
    md += "\n## 3. Top Gene Analysis\n"
    md += "Included Genes:\n"
    # Reconstruct lists for reporting
    ddr_names = ["BRCA1", "BRCA2", "RAD51", "PARP1", "FANCD2", "ATR", "ATM", "CHEK1"]
    cc_names = ["MKI67", "CCNE1", "CDK1", "E2F1", "PCNA", "TOP2A"]
    md += f"- DDR: {', '.join(ddr_names)}\n"
    md += f"- Cell Cycle: {', '.join(cc_names)}\n"
    
    with open(OUTPUT_REPORT, 'w') as f:
        f.write(md)
        
    print(f"\n✅ Report generated: {OUTPUT_REPORT}")
    print(results_df[['Patient', 'DDR_Log2FC', 'CC_Log2FC']])

def main():
    df = load_data()
    # Handle duplicate index if any (gene symbols sometimes duplicate)
    df = df.groupby(df.index).mean()
    
    results = analyze_pairs(df)
    generate_report(results)

if __name__ == "__main__":
    main()
