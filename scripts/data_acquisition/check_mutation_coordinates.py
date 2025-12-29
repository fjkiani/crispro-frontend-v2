#!/usr/bin/env python3
"""Check if mutations have coordinates for TRUE-SAE extraction."""

import json
from pathlib import Path

mut_file = Path("data/external/ov_platinum_non_tcga/raw/ovarian_msk_2025_mutations.json")

with open(mut_file) as f:
    mutations = json.load(f)

print(f"Total mutations: {len(mutations)}")
print()

if mutations:
    m = mutations[0]
    print("Sample mutation keys:")
    for key in sorted(m.keys()):
        print(f"  - {key}: {type(m[key]).__name__}")
    
    print()
    print("Coordinate fields check:")
    print(f"  chromosome: {'chromosome' in m and m['chromosome']}")
    print(f"  startPosition: {'startPosition' in m and m['startPosition']}")
    print(f"  endPosition: {'endPosition' in m and m['endPosition']}")
    print(f"  referenceAllele: {'referenceAllele' in m and m['referenceAllele']}")
    print(f"  variantAllele: {'variantAllele' in m and m['variantAllele']}")
    
    print()
    print("Example mutation (coordinate fields):")
    coord_fields = ["chromosome", "startPosition", "endPosition", "referenceAllele", "variantAllele", "gene", "proteinChange"]
    example = {k: m.get(k) for k in coord_fields if k in m}
    print(json.dumps(example, indent=2))
    
    # Count mutations with coordinates
    with_coords = sum(1 for m in mutations if m.get("chromosome") and m.get("startPosition"))
    print()
    print(f"Mutations with coordinates: {with_coords} / {len(mutations)} ({100*with_coords/len(mutations):.1f}%)")

