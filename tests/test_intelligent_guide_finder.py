import os
import sys
from dotenv import load_dotenv

# --- Setup ---
# Load environment variables from .env file
load_dotenv()
# Add the project root to the Python path to allow importing from 'tools'
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# --- Temp Override for Testing ---
# This hardcoded endpoint is causing issues. It points to a specific,
# potentially outdated deployment. Removing it will allow the test to use the
# endpoint defined in the .env file, which is the correct behavior.
# os.environ["EVO2_GENERATIVE_ENDPOINT"] = 'https://fjkiani--evo2-sequence-generator-1b-evo2generator1b-web.modal.run'

from tools.intelligent_guide_finder import find_intelligent_guides

def test_end_to_end_guide_finder():
    """
    Tests the full end-to-end workflow of the Intelligent Guide Finder.
    
    This test will:
    1. Define a target locus (e.g., for the KRAS gene).
    2. Call the main `find_intelligent_guides` controller function.
    3. Print the results in a readable format.
    """
    print("\\n--- Testing Intelligent Guide Finder End-to-End Workflow ---")
    
    # Define a specific, small locus in a well-known gene for a focused test.
    # This targets a region in the KRAS gene.
    target_locus = "chr12:25,245,200-25,245,400" 
    target_genome = "hg38"
    
    print(f"\\nTargeting Locus: {target_locus} in Genome: {target_genome}")
    
    # Call the main controller function
    # This function orchestrates sequence retrieval, candidate generation,
    # AI-powered efficacy scoring, and off-target safety analysis.
    ranked_guides = find_intelligent_guides(locus=target_locus, genome=target_genome)
    
    print("\\n--- Results ---")
    if not ranked_guides:
        print("⚠️ Test finished, but no guides were found or an error occurred.")
        return

    print(f"Successfully found and ranked {len(ranked_guides)} guides.")
    print("Top 5 ranked guides:")
    print("-" * 50)
    # Print the top 5 guides for review
    for i, guide in enumerate(ranked_guides[:5]):
        print(f"Rank {i+1}:")
        print(f"  Sequence: {guide.get('guide_sequence')}")
        print(f"  On-Target Score (AI): {guide.get('on_target_score', 0.0):.2f}")
        print(f"  Off-Target Hits: {guide.get('off_target_hits', 'N/A')}")
        print("-" * 20)
        
    print("\\n✅ Test finished successfully. Review the ranked guides above.")

if __name__ == "__main__":
    test_end_to_end_guide_finder() 