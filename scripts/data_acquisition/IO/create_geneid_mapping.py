#!/usr/bin/env python3
"""
Quick GeneID → Symbol mapping for GSE91061
Uses mygene API for bulk mapping (fastest approach per ZETA DOCTRINE)
"""

import pandas as pd
import requests
import json
import time
from tqdm import tqdm

print("="*80)
print("GENEID → SYMBOL MAPPING (GSE91061)")
print("="*80)

# Load expression matrix to get all GeneIDs
print("\n[1] Loading expression matrix...")
expr = pd.read_csv('GSE91061_norm_counts_TPM_GRCh38.p13_NCBI.tsv', 
                   sep='\t', index_col=0, nrows=0)  # Just get index

# Get all GeneIDs
all_gene_ids = pd.read_csv('GSE91061_norm_counts_TPM_GRCh38.p13_NCBI.tsv', 
                           sep='\t', index_col=0).index.tolist()

print(f"Total GeneIDs: {len(all_gene_ids)}")

# Use mygene API (bulk query - 1000 at a time)
print("\n[2] Mapping GeneIDs to symbols via mygene API...")
url = 'https://mygene.info/v3/query'

mapping = {}
batch_size = 1000

for i in tqdm(range(0, len(all_gene_ids), batch_size)):
    batch = all_gene_ids[i:i+batch_size]
    
    params = {
        'q': ','.join(map(str, batch)),
        'scopes': 'entrezgene',
        'fields': 'symbol,entrezgene',
        'species': 'human'
    }
    
    try:
        response = requests.post(url, data=params, timeout=30)
        if response.ok:
            data = response.json()
            for item in data:
                if 'entrezgene' in item and 'symbol' in item:
                    gene_id = str(item['entrezgene'])
                    symbol = item['symbol']
                    mapping[gene_id] = symbol
        time.sleep(0.5)  # Rate limiting
    except Exception as e:
        print(f"⚠️  Batch {i//batch_size + 1} failed: {e}")
        continue

print(f"\n✓ Mapped {len(mapping)}/{len(all_gene_ids)} genes ({len(mapping)/len(all_gene_ids)*100:.1f}%)")

# Save mapping
output_file = 'gse91061_geneid_to_symbol.json'
with open(output_file, 'w') as f:
    json.dump(mapping, f, indent=2)

print(f"✓ Saved mapping to {output_file}")

# Create reverse mapping (symbol → GeneID) for pathway lookup
reverse_mapping = {v: k for k, v in mapping.items()}
with open('gse91061_symbol_to_geneid.json', 'w') as f:
    json.dump(reverse_mapping, f, indent=2)

print(f"✓ Saved reverse mapping to gse91061_symbol_to_geneid.json")

print("\n" + "="*80)
print("MAPPING COMPLETE")
print("="*80)
