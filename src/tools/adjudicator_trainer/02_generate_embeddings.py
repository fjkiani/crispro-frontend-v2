import pandas as pd
import asyncio
import aiohttp
from tqdm.asyncio import tqdm_asyncio
import os
import json

# --- Configuration ---
INPUT_FILE = "data/adjudicator_training/clinvar_missense_labels.csv"
OUTPUT_FILE = "data/adjudicator_training/labeled_embeddings.parquet"
ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run/invoke"
CONCURRENT_REQUESTS = 100
SAVE_INTERVAL = 500

# --- Main Logic ---

async def get_embedding(session, variant_id, gene, hgvsp):
    """
    Asynchronously fetches an embedding for a single variant from the Zeta Oracle.
    Returns the embedding vector or None if an error occurs.
    """
    payload = {
        "action": "embed",
        "params": {
            "gene_symbol": gene,
            "protein_change": hgvsp
        }
    }
    try:
        async with session.post(ZETA_ORACLE_URL, json=payload, timeout=60) as response:
            if response.status == 200:
                data = await response.json()
                # The oracle's 'embed' action returns a dictionary with the embedding
                return data.get("embedding")
            else:
                print(f"Error for {variant_id} ({gene} {hgvsp}): Status {response.status} | Response: {await response.text()}")
                return None
    except Exception as e:
        print(f"Exception for {variant_id} ({gene} {hgvsp}): {e}")
        return None

async def main():
    """
    Main orchestration function to generate and save embeddings for the ClinVar dataset.
    """
    print("🚀 Starting Operation Adjudicator: Phase 1.3 - Generate Embeddings")
    
    # Load the source data
    if not os.path.exists(INPUT_FILE):
        print(f"FATAL: Input file not found at {INPUT_FILE}")
        return
    source_df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(source_df)} variants from {INPUT_FILE}")

    # Setup for async requests
    all_results = []
    if os.path.exists(OUTPUT_FILE):
        all_results = [pd.read_parquet(OUTPUT_FILE)]
        processed_ids = set(all_results[0]['clinvar_id'])
        source_df = source_df[~source_df['clinvar_id'].isin(processed_ids)]
        print(f"Resuming with {len(source_df)} remaining variants.")
    
    async with aiohttp.ClientSession() as session:
        tasks = []
        for index, row in tqdm_asyncio(source_df.iterrows(), total=len(source_df), desc="Dispatching jobs"):
            tasks.append(get_embedding(session, row['clinvar_id'], row['gene'], row['hgvsp']))
            
            if len(tasks) >= CONCURRENT_REQUESTS:
                # Process a chunk of tasks
                processed_embeddings = await asyncio.gather(*tasks)
                tasks = [] # Reset tasks
                
                # Save chunk to results
                chunk_df = source_df.loc[[idx for idx, emb in zip(source_df.index[:len(processed_embeddings)], processed_embeddings) if emb is not None]].copy()
                chunk_df['embedding'] = [json.dumps(emb) for emb in processed_embeddings if emb is not None]
                
                if not chunk_df.empty:
                    all_results.append(chunk_df)
                    # Incremental save
                    if len(all_results) > 1:
                         pd.concat(all_results, ignore_index=True).to_parquet(OUTPUT_FILE, index=False)
                         print(f"✅ Incremental save: {sum(len(df) for df in all_results)} records saved.")

        # Process any remaining tasks
        if tasks:
            processed_embeddings = await asyncio.gather(*tasks)
            chunk_df = source_df.loc[[idx for idx, emb in zip(source_df.index[-len(processed_embeddings):], processed_embeddings) if emb is not None]].copy()
            chunk_df['embedding'] = [json.dumps(emb) for emb in processed_embeddings if emb is not None]
            if not chunk_df.empty:
                all_results.append(chunk_df)

    # Final save
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        final_df.to_parquet(OUTPUT_FILE, index=False)
        print(f"\nSaved {len(final_df)} total labeled embeddings to {OUTPUT_FILE}")

    print("\n🔥 Operation Adjudicator: Phase 1.3 Complete.")


if __name__ == "__main__":
    asyncio.run(main()) 