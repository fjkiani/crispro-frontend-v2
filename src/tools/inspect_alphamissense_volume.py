import modal
import pandas as pd

# Define the same volume used by the main service
alphamissense_volume = modal.Volume.from_name("alphamissense-data")
image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas", "pyarrow")
app = modal.App("alphamissense-inspector", image=image)

@app.function(volumes={"/data": alphamissense_volume})
def inspect_volume():
    """
    A one-off function to inspect the AlphaMissense data in the Modal Volume.
    """
    data_path = "/data/AlphaMissense_hg38.parquet"
    print(f"--- 🔍 Inspecting AlphaMissense data in Modal Volume at {data_path} ---")
    
    try:
        df = pd.read_parquet(data_path, columns=['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity'])
        
        # Coordinates for BRAF gene on hg38
        chrom = 'chr7'
        start_pos = 140719327
        end_pos = 140924929

        print(f"Filtering for BRAF variants ({chrom}:{start_pos}-{end_pos})...")
        df_braf = df[(df['#CHROM'] == chrom) & (df['POS'] >= start_pos) & (df['POS'] <= end_pos)]
        
        if df_braf.empty:
            print("❌ No BRAF variants found in the specified coordinate range.")
            return

        # Find a high-confidence pathogenic variant
        pathogenic_variants = df_braf[df_braf['am_pathogenicity'] > 0.9] # High threshold for confidence
        
        if pathogenic_variants.empty:
            print("❌ No high-confidence pathogenic BRAF variants found (score > 0.9).")
        else:
            print("✅ Found high-confidence pathogenic BRAF variants. Here is one:")
            variant = pathogenic_variants.iloc[0]
            variant_str = f"{variant['#CHROM']}:{variant['POS']}:{variant['REF']}:{variant['ALT']}"
            print(f"   Variant: {variant_str}")
            print(f"   Pathogenicity Score: {variant['am_pathogenicity']}")
            print("\nUse this variant string to test the fusion-engine.")

    except Exception as e:
        print(f"💥 An error occurred during inspection: {e}")

@app.local_entrypoint()
def main():
    inspect_volume.remote() 