import pandas as pd
import os

def find_brca1_variant_pair():
    """
    Parses the Findlay et al. (2018) BRCA1 dataset to find a pair of variants
    at the same position, one pathogenic (LOF) and one benign (FUNC/INT).
    """
    print("🚀 Initiating BRCA1 Target Acquisition...")
    
    # Load the dataset as identified in the notebook
    brca1_df = pd.read_excel(
        os.path.join('scripts', 'evo2', 'evo2', 'notebooks', 'brca1', '41586_2018_461_MOESM3_ESM.xlsx'),
        header=2,
    )
    
    # Basic cleaning as per the notebook
    brca1_df.rename(columns={
        'position (hg19)': 'pos',
        'reference': 'ref',
        'alt': 'alt',
        'func.class': 'class',
    }, inplace=True)
    brca1_df['class'] = brca1_df['class'].replace(['FUNC', 'INT'], 'FUNC/INT')
    
    # Find positions that have both LOF and FUNC/INT variants
    pos_counts = brca1_df.groupby('pos')['class'].nunique()
    target_positions = pos_counts[pos_counts > 1].index
    
    if len(target_positions) == 0:
        print("❌ No positions found with both LOF and FUNC/INT variants.")
        return None, None

    # Select the first such position for our test
    target_pos = target_positions[0]
    variants_at_pos = brca1_df[brca1_df['pos'] == target_pos]
    
    lof_variant = variants_at_pos[variants_at_pos['class'] == 'LOF'].iloc[0]
    func_variant = variants_at_pos[variants_at_pos['class'] == 'FUNC/INT'].iloc[0]
    
    print("\n--- TARGETS ACQUIRED ---")
    print(f"Genomic Position (hg19): {target_pos}")
    print("\n🔥 Pathogenic (LOF) Variant:")
    print(lof_variant[['ref', 'alt', 'class']])
    print("\n✅ Benign (FUNC/INT) Variant:")
    print(func_variant[['ref', 'alt', 'class']])
    print("------------------------\n")
    
    return lof_variant.to_dict(), func_variant.to_dict()

if __name__ == "__main__":
    find_brca1_variant_pair() 
import os

def find_brca1_variant_pair():
    """
    Parses the Findlay et al. (2018) BRCA1 dataset to find a pair of variants
    at the same position, one pathogenic (LOF) and one benign (FUNC/INT).
    """
    print("🚀 Initiating BRCA1 Target Acquisition...")
    
    # Load the dataset as identified in the notebook
    brca1_df = pd.read_excel(
        os.path.join('scripts', 'evo2', 'evo2', 'notebooks', 'brca1', '41586_2018_461_MOESM3_ESM.xlsx'),
        header=2,
    )
    
    # Basic cleaning as per the notebook
    brca1_df.rename(columns={
        'position (hg19)': 'pos',
        'reference': 'ref',
        'alt': 'alt',
        'func.class': 'class',
    }, inplace=True)
    brca1_df['class'] = brca1_df['class'].replace(['FUNC', 'INT'], 'FUNC/INT')
    
    # Find positions that have both LOF and FUNC/INT variants
    pos_counts = brca1_df.groupby('pos')['class'].nunique()
    target_positions = pos_counts[pos_counts > 1].index
    
    if len(target_positions) == 0:
        print("❌ No positions found with both LOF and FUNC/INT variants.")
        return None, None

    # Select the first such position for our test
    target_pos = target_positions[0]
    variants_at_pos = brca1_df[brca1_df['pos'] == target_pos]
    
    lof_variant = variants_at_pos[variants_at_pos['class'] == 'LOF'].iloc[0]
    func_variant = variants_at_pos[variants_at_pos['class'] == 'FUNC/INT'].iloc[0]
    
    print("\n--- TARGETS ACQUIRED ---")
    print(f"Genomic Position (hg19): {target_pos}")
    print("\n🔥 Pathogenic (LOF) Variant:")
    print(lof_variant[['ref', 'alt', 'class']])
    print("\n✅ Benign (FUNC/INT) Variant:")
    print(func_variant[['ref', 'alt', 'class']])
    print("------------------------\n")
    
    return lof_variant.to_dict(), func_variant.to_dict()

if __name__ == "__main__":
    find_brca1_variant_pair() 
 
 
 
 
 