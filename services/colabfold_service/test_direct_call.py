#!/usr/bin/env python3
"""
Direct call to deployed ColabFold service - bypass local timeout
"""
import modal

# Connect to deployed app
app = modal.App.lookup("colabfold-smoke-test", create_if_missing=False)

# Get function handle
predict_structure = modal.Function.lookup("colabfold-smoke-test", "predict_structure")

# Tiny test protein (47 residues)
fasta_content = """>test_tiny
MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV
"""

print("=" * 80)
print("🔬 DIRECT MODAL CALL - NO LOCAL TIMEOUT")
print("=" * 80)
print(f"FASTA: {len(fasta_content)} bytes")
print("Calling deployed function (this may take 20-30 min)...")
print("Check Modal dashboard: https://modal.com/apps/crispro/main")
print("=" * 80)

# Call deployed function (server-side timeout = 900s from decorator)
result = predict_structure.remote(fasta_content, job_id="direct_test")

print("\n" + "=" * 80)
print("✅ RESULT:")
print("=" * 80)
print(f"Success: {result.get('success', False)}")
print(f"pLDDT Mean: {result.get('plddt_mean', 'N/A')}")
print(f"PAE Mean: {result.get('pae_mean', 'N/A')}")
print(f"PDB Size: {result.get('pdb_size', 'N/A')} bytes")

if result.get('error'):
    print(f"\n❌ ERROR: {result['error']}")
    
if result.get('pdb_content'):
    # Save PDB
    with open("/tmp/colabfold_output.pdb", "w") as f:
        f.write(result['pdb_content'])
    print(f"\n💾 Saved PDB to: /tmp/colabfold_output.pdb")

print("=" * 80)

