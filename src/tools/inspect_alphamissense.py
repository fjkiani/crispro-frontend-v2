import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def inspect_alphamissense_data(data_path, gene_symbol, chrom):
    """
    Inspects the AlphaMissense Parquet file to find variants for a specific gene.
    """
    try:
        logger.info(f"Loading AlphaMissense data from {data_path}...")
        df = pd.read_parquet(data_path, columns=['#CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'am_pathogenicity'])
        
        logger.info(f"Filtering for chromosome '{chrom}'...")
        df_chrom = df[df['#CHROM'] == chrom]

        if df_chrom.empty:
            logger.warning(f"No data found for chromosome '{chrom}'.")
            return

        # Unfortunately, the main file doesn't have gene symbols. 
        # We can only search by coordinate.
        # This function is therefore less useful than anticipated.
        # We will try a different approach.
        # For now, let's just dump the first few VWF variants we can find.
        # We know from the clinvar data that VWF is on chromosome 12.
        # A quick web search reveals VWF is at: chr12:5,948,013-6,126,523 on hg38.
        
        start_pos = 5948013
        end_pos = 6126523

        logger.info(f"Filtering for VWF gene coordinates: chr12:{start_pos}-{end_pos}")
        df_vwf = df_chrom[(df_chrom['POS'] >= start_pos) & (df_chrom['POS'] <= end_pos)]
        
        if df_vwf.empty:
            logger.warning("No variants found for VWF in the specified coordinate range.")
            return

        print("Found variants for VWF. Displaying first 20 pathogenic variants:")
        pathogenic_vwf = df_vwf[df_vwf['am_pathogenicity'] > 0.564] # Threshold for likely pathogenic
        print(pathogenic_vwf.head(20).to_string())

    except FileNotFoundError:
        logger.error(f"FATAL: AlphaMissense data file not found at '{data_path}'.")
        logger.error("Please ensure the file exists and the path is correct.")
    except Exception as e:
        logger.error(f"An error occurred: {e}")

if __name__ == '__main__':
    data_path = "data/alphamissense/AlphaMissense_hg38.parquet"
    inspect_alphamissense_data(data_path, "VWF", "chr12") 