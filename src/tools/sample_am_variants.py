import modal
import pandas as pd
import json

# Use the same volume that holds AlphaMissense parquet
alphamissense_volume = modal.Volume.from_name("alphamissense-data")
image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas", "pyarrow")
app = modal.App("alphamissense-sampler", image=image)

@app.function(volumes={"/data": alphamissense_volume})
def sample_variants(n_pos: int = 10, n_neg: int = 10, parquet_path: str = "/data/AlphaMissense_hg38.parquet") -> str:
    df = pd.read_parquet(parquet_path, columns=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity"])  # type: ignore
    # High-confidence positives and negatives
    pos_df = df[df["am_pathogenicity"] >= 0.9].sample(min(n_pos, len(df[df["am_pathogenicity"] >= 0.9])), random_state=42)
    neg_df = df[df["am_pathogenicity"] <= 0.1].sample(min(n_neg, len(df[df["am_pathogenicity"] <= 0.1])), random_state=42)
    out = []
    for _, r in pd.concat([pos_df, neg_df], axis=0).iterrows():
        out.append({
            "chrom": str(r["#CHROM"]),
            "pos": int(r["POS"]),
            "ref": str(r["REF"]),
            "alt": str(r["ALT"]),
            "am_pathogenicity": float(r["am_pathogenicity"]),
            "label": 1 if float(r["am_pathogenicity"]) >= 0.9 else 0,
        })
    return json.dumps({"variants": out})

@app.local_entrypoint()
def main():
    import os
    n_pos = int(os.getenv("N_POS", "10"))
    n_neg = int(os.getenv("N_NEG", "10"))
    print(sample_variants.remote(n_pos=n_pos, n_neg=n_neg))



import json

# Use the same volume that holds AlphaMissense parquet
alphamissense_volume = modal.Volume.from_name("alphamissense-data")
image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas", "pyarrow")
app = modal.App("alphamissense-sampler", image=image)

@app.function(volumes={"/data": alphamissense_volume})
def sample_variants(n_pos: int = 10, n_neg: int = 10, parquet_path: str = "/data/AlphaMissense_hg38.parquet") -> str:
    df = pd.read_parquet(parquet_path, columns=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity"])  # type: ignore
    # High-confidence positives and negatives
    pos_df = df[df["am_pathogenicity"] >= 0.9].sample(min(n_pos, len(df[df["am_pathogenicity"] >= 0.9])), random_state=42)
    neg_df = df[df["am_pathogenicity"] <= 0.1].sample(min(n_neg, len(df[df["am_pathogenicity"] <= 0.1])), random_state=42)
    out = []
    for _, r in pd.concat([pos_df, neg_df], axis=0).iterrows():
        out.append({
            "chrom": str(r["#CHROM"]),
            "pos": int(r["POS"]),
            "ref": str(r["REF"]),
            "alt": str(r["ALT"]),
            "am_pathogenicity": float(r["am_pathogenicity"]),
            "label": 1 if float(r["am_pathogenicity"]) >= 0.9 else 0,
        })
    return json.dumps({"variants": out})

@app.local_entrypoint()
def main():
    import os
    n_pos = int(os.getenv("N_POS", "10"))
    n_neg = int(os.getenv("N_NEG", "10"))
    print(sample_variants.remote(n_pos=n_pos, n_neg=n_neg))









import json

# Use the same volume that holds AlphaMissense parquet
alphamissense_volume = modal.Volume.from_name("alphamissense-data")
image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas", "pyarrow")
app = modal.App("alphamissense-sampler", image=image)

@app.function(volumes={"/data": alphamissense_volume})
def sample_variants(n_pos: int = 10, n_neg: int = 10, parquet_path: str = "/data/AlphaMissense_hg38.parquet") -> str:
    df = pd.read_parquet(parquet_path, columns=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity"])  # type: ignore
    # High-confidence positives and negatives
    pos_df = df[df["am_pathogenicity"] >= 0.9].sample(min(n_pos, len(df[df["am_pathogenicity"] >= 0.9])), random_state=42)
    neg_df = df[df["am_pathogenicity"] <= 0.1].sample(min(n_neg, len(df[df["am_pathogenicity"] <= 0.1])), random_state=42)
    out = []
    for _, r in pd.concat([pos_df, neg_df], axis=0).iterrows():
        out.append({
            "chrom": str(r["#CHROM"]),
            "pos": int(r["POS"]),
            "ref": str(r["REF"]),
            "alt": str(r["ALT"]),
            "am_pathogenicity": float(r["am_pathogenicity"]),
            "label": 1 if float(r["am_pathogenicity"]) >= 0.9 else 0,
        })
    return json.dumps({"variants": out})

@app.local_entrypoint()
def main():
    import os
    n_pos = int(os.getenv("N_POS", "10"))
    n_neg = int(os.getenv("N_NEG", "10"))
    print(sample_variants.remote(n_pos=n_pos, n_neg=n_neg))



import json

# Use the same volume that holds AlphaMissense parquet
alphamissense_volume = modal.Volume.from_name("alphamissense-data")
image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas", "pyarrow")
app = modal.App("alphamissense-sampler", image=image)

@app.function(volumes={"/data": alphamissense_volume})
def sample_variants(n_pos: int = 10, n_neg: int = 10, parquet_path: str = "/data/AlphaMissense_hg38.parquet") -> str:
    df = pd.read_parquet(parquet_path, columns=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity"])  # type: ignore
    # High-confidence positives and negatives
    pos_df = df[df["am_pathogenicity"] >= 0.9].sample(min(n_pos, len(df[df["am_pathogenicity"] >= 0.9])), random_state=42)
    neg_df = df[df["am_pathogenicity"] <= 0.1].sample(min(n_neg, len(df[df["am_pathogenicity"] <= 0.1])), random_state=42)
    out = []
    for _, r in pd.concat([pos_df, neg_df], axis=0).iterrows():
        out.append({
            "chrom": str(r["#CHROM"]),
            "pos": int(r["POS"]),
            "ref": str(r["REF"]),
            "alt": str(r["ALT"]),
            "am_pathogenicity": float(r["am_pathogenicity"]),
            "label": 1 if float(r["am_pathogenicity"]) >= 0.9 else 0,
        })
    return json.dumps({"variants": out})

@app.local_entrypoint()
def main():
    import os
    n_pos = int(os.getenv("N_POS", "10"))
    n_neg = int(os.getenv("N_NEG", "10"))
    print(sample_variants.remote(n_pos=n_pos, n_neg=n_neg))








