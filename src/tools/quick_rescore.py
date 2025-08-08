#!/usr/bin/env python3
"""
Quick script to regenerate TP53 scores with corrected parsing logic.
"""

import json
import pandas as pd
from zeta_striker import call_command_center_assess_threat, fuse_intelligence, final_analysis

def main():
    # Load the existing mutation data to get the list of mutations
    mutation_df = pd.read_csv('data/data_mutations.txt', sep='\t', comment='#', low_memory=False)
    target_mutations = mutation_df[mutation_df['Hugo_Symbol'] == 'TP53'].copy()
    
    # Filter for protein-altering mutations
    target_mutations = target_mutations[target_mutations['Variant_Classification'].isin([
        'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 
        'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins'
    ])]
    
    # Get unique protein changes
    unique_protein_changes = target_mutations['HGVSp_Short'].dropna().unique()
    print(f"Found {len(unique_protein_changes)} unique mutations to rescore.")
    
    # Focus on a subset that we know includes both truncation and missense mutations
    priority_mutations = []
    for pc in unique_protein_changes:
        if '*' in pc or 'fs' in pc:  # Truncation/frameshift mutations
            priority_mutations.append(pc)
        elif len(priority_mutations) < 20:  # Add some missense mutations
            priority_mutations.append(pc)
        if len(priority_mutations) >= 30:  # Limit to 30 for speed
            break
    
    print(f"Rescoring {len(priority_mutations)} priority mutations...")
    
    # Rescore the mutations
    scored_mutations = {}
    for pc in priority_mutations:
        protein_change_for_api = pc.replace('p.', '')
        assessment = call_command_center_assess_threat('TP53', protein_change_for_api)
        scored_mutations[pc] = assessment
        print(f"  {pc}: score={assessment['delta_likelihood_score']}, pathogenic={assessment['is_pathogenic']}")
    
    # Save the new scores
    with open('results/TP53_zeta_scores_corrected.json', 'w') as f:
        json.dump(scored_mutations, f, indent=4)
    
    # Regenerate the master file
    clinical_df = pd.read_csv('data/tcga.tsv', sep='\t', comment='#')
    master_df = fuse_intelligence('data/tcga.tsv', scored_mutations, target_mutations)
    
    if master_df is not None:
        master_df.to_csv('results/master_clinical_zeta_scores_corrected.tsv', sep='\t', index=False)
        print("✅ Corrected master file saved.")
        
        # Run the final analysis
        final_analysis('results/master_clinical_zeta_scores_corrected.tsv')

if __name__ == "__main__":
    main() 