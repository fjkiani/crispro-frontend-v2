import requests
import logging
import re
from typing import Optional
import time
import json
import os

from Bio import Align

from tools.llm_api import query_llm

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

UCSC_API_URL = "https://api.genome.ucsc.edu/getData/sequence"

def get_sequence_for_locus(locus: str, genome: str, window_size: int = 1000) -> Optional[str]:
    """
    Fetches a DNA sequence for a given genomic locus using the UCSC API.

    This function handles two types of locus formats:
    1. Gene symbol (e.g., "BRAF") - NOT YET IMPLEMENTED.
    2. Genomic coordinates (e.g., "chr7:140,750,000-140,760,000" or "chr7:140753336").

    Args:
        locus: The genomic locus (gene symbol or coordinates).
        genome: The genome build (e.g., 'hg38').
        window_size: The size of the window to fetch around a single point coordinate.

    Returns:
        The DNA sequence as a string, or None if an error occurs.
    """
    # Regular expression to match coordinate formats like chr:start-end or chr:pos
    coord_match = re.match(r'^(chr[A-Za-z0-9]+):([\d,]+)-?([\d,]*)', locus.strip())
    
    if not coord_match:
        logger.error(f"Locus format not recognized or not yet supported: {locus}")
        # Placeholder for future gene symbol to coordinate conversion
        return None

    chrom = coord_match.group(1)
    start_str = coord_match.group(2).replace(',', '')
    end_str = coord_match.group(3).replace(',', '') if coord_match.group(3) else None

    try:
        start = int(start_str)
        if end_str:
            end = int(end_str)
        else:
            # If no end is provided, create a window around the start position
            start = start - (window_size // 2)
            end = start + window_size
        
        if start < 0:
            start = 0

    except ValueError:
        logger.error(f"Invalid coordinate numbers in locus: {locus}")
        return None

    params = {
        "genome": genome,
        "chrom": chrom,
        "start": start,
        "end": end
    }

    try:
        logger.info(f"Fetching sequence from UCSC API with params: {params}")
        response = requests.get(UCSC_API_URL, params=params, verify=False) # verify=False for local dev envs
        response.raise_for_status()
        data = response.json()
        
        sequence = data.get("dna")
        if not sequence:
            logger.error(f"No DNA sequence returned from UCSC API for locus {locus}.")
            return None
            
        logger.info(f"Successfully fetched {len(sequence)}bp sequence for {locus}.")
        return sequence.upper()

    except requests.exceptions.RequestException as e:
        logger.error(f"API call to UCSC failed for locus {locus}: {e}")
        return None
    except Exception as e:
        logger.error(f"An unexpected error occurred while fetching sequence for {locus}: {e}")
        return None

def find_guide_candidates(sequence: str, pam: str = "NGG", guide_length: int = 20) -> list[str]:
    """
    Finds all potential guide RNA candidates in a DNA sequence.

    This function scans a sequence for a given PAM (e.g., "NGG") and returns
    a list of all valid guide sequences of a specified length that precede it.

    Args:
        sequence: The DNA sequence to search within.
        pam: The Protospacer Adjacent Motif to search for. 'N' is treated as any base.
        guide_length: The length of the guide RNA sequence.

    Returns:
        A list of unique 20-bp guide RNA sequences (strings).
    """
    pam_regex = pam.replace('N', '[ATCG]')
    candidates = set()
    
    # Find all occurrences of the PAM sequence
    for match in re.finditer(pam_regex, sequence, re.IGNORECASE):
        pam_start_index = match.start()
        
        # The guide is the sequence immediately preceding the PAM
        guide_start_index = pam_start_index - guide_length
        
        if guide_start_index >= 0:
            guide_sequence = sequence[guide_start_index:pam_start_index]
            candidates.add(guide_sequence)
            
    logger.info(f"Found {len(candidates)} unique guide candidates.")
    return list(candidates)

def score_sequence_importance(sequence: str) -> float:
    """
    Scores the functional importance of a DNA sequence using the "Generate & Compare" method.

    This function uses our deployed generative Evo2 model to predict an "ideal"
    sequence based on the provided guide sequence. It then aligns the original
    and generated sequences and returns the alignment score. A higher score
    implies the original sequence is highly conserved and functionally important.

    Args:
        sequence: The DNA sequence (e.g., a guide candidate) to score.

    Returns:
        The alignment score as a float. Returns 0.0 if an error occurs.
    """
    if not sequence:
        logger.warning("score_sequence_importance received an empty sequence.")
        return 0.0

    try:
        # --- 1. GENERATE: Call the Evo2 generative model directly ---
        endpoint = os.environ.get("EVO2_GENERATIVE_ENDPOINT")
        if not endpoint:
            logger.error("EVO2_GENERATIVE_ENDPOINT environment variable not set.")
            return 0.0
        
        prompt = f"Complete the following DNA sequence: {sequence}"
        payload = {"prompt": prompt}
        headers = {"Content-Type": "application/json"}
        
        logger.info(f"Calling generative model directly at {endpoint} for sequence: {sequence}")
        
        response = requests.post(endpoint, headers=headers, json=payload, timeout=90, verify=False)
        response.raise_for_status()
        
        response_json = response.json()
        generated_sequence = response_json.get("completion")

        if not generated_sequence or not isinstance(generated_sequence, str) or "Error" in generated_sequence:
            logger.error(f"Failed to get a valid generated sequence for '{sequence}': {generated_sequence}")
            return 0.0
        
        # Evo2 can return sequences with non-DNA characters or in lowercase.
        cleaned_generated_sequence = re.sub(r'[^ATCG]', '', generated_sequence.upper())

        if not cleaned_generated_sequence:
            logger.error(f"Generated sequence was empty after cleaning: '{generated_sequence}'")
            return 0.0

        # --- 2. COMPARE: Align original and generated sequences ---
        logger.info(f"Aligning original '{sequence}' with generated '{cleaned_generated_sequence}'")
        
        # Use the modern Bio.Align.PairwiseAligner
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        # Set scoring to reward matches and heavily penalize mismatches/gaps
        aligner.match_score = 1
        aligner.mismatch_score = -2
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        
        # Perform the alignment
        alignments = aligner.align(sequence.upper(), cleaned_generated_sequence)
        
        # We expect one global alignment
        if not alignments:
            logger.error(f"Alignment failed for '{sequence}'.")
            return 0.0
            
        # The score from the aligner is a raw score.
        # To make it comparable, we can normalize it by the length of the original sequence.
        # This gives a sense of "percent match quality"
        raw_score = alignments[0].score
        normalized_score = raw_score / len(sequence)

        logger.info(f"Alignment for '{sequence}': Raw Score={raw_score}, Normalized Score={normalized_score:.4f}")
        
        return float(normalized_score)

    except ImportError:
        logger.error("BioPython is not installed. Please install it with 'pip install biopython'")
        return 0.0
    except requests.exceptions.HTTPError as e:
        logger.error(f"HTTP Error during sequence scoring for '{sequence}': {e}. Response: {e.response.text}")
        return 0.0
    except requests.exceptions.RequestException as e:
        logger.error(f"Request failed during sequence scoring for '{sequence}': {e}")
        return 0.0
    except json.JSONDecodeError:
        logger.error(f"Failed to decode JSON from model response for '{sequence}'.")
        return 0.0
    except Exception as e:
        logger.error(f"An unexpected error occurred during sequence scoring for '{sequence}': {e}")
        return 0.0

def find_potential_off_targets(guide_sequence: str, genome: str = "hg38") -> dict:
    """
    Uses the NCBI BLAST API to find potential off-target sites for a guide RNA.
    This function submits a search, polls for the result, and retrieves the hits.
    """
    blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    payload = {
        'CMD': 'Put', 'PROGRAM': 'blastn', 'DATABASE': 'nt',
        'QUERY': guide_sequence, 'MEGABLAST': 'on', 'EXPECT': '1000',
        'WORD_SIZE': 7, 'NUCL_REWARD': 1, 'NUCL_PENALTY': -3,
        'FILTER': 'F', 'FORMAT_TYPE': 'JSON2'
    }

    try:
        logger.info(f"Submitting BLAST search for guide: {guide_sequence}")
        submit_response = requests.post(blast_url, data=payload, verify=False)
        submit_response.raise_for_status()
        
        rid_match = re.search(r'RID = (\w+)', submit_response.text)
        if not rid_match:
            logger.error(f"Could not find RID in BLAST submission response: {submit_response.text[:500]}")
            return {"error": "Failed to submit BLAST search."}
        rid = rid_match.group(1)
        logger.info(f"BLAST search submitted. RID: {rid}")

        # Step 2: Poll for results with a timeout
        status_payload = {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'JSON2'}
        max_retries = 24 # Poll for up to 2 minutes (24 * 5s = 120s)
        for i in range(max_retries):
            time.sleep(5)
            logger.info(f"Checking BLAST status for RID: {rid} (Attempt {i+1}/{max_retries})")
            
            status_response = requests.post(blast_url, data=status_payload, verify=False)
            status_response.raise_for_status()

            # First, try to parse JSON. If it works, results are ready.
            try:
                json_response = status_response.json()
                # Check for the specific structure of a results object.
                if json_response.get("BlastOutput2"):
                    logger.info(f"BLAST search for {rid} is ready (JSON response). Fetching results.")
                    return json_response # Results are directly in this response
                # It might be a JSON response that still says waiting.
                elif "Status=WAITING" in status_response.text:
                     logger.info(f"BLAST search for {rid} is still running (JSON status).")
                     continue
            except json.JSONDecodeError:
                # If JSON parsing fails, check the text content for status.
                if 'Status=READY' in status_response.text:
                    logger.info(f"BLAST search for {rid} is ready (Text status). Fetching results.")
                    # The results are ready, break the loop and fetch them below.
                    break 
                elif 'Status=WAITING' in status_response.text:
                    logger.info(f"BLAST search for {rid} is still running (Text status).")
                    continue # Continue polling
                else:
                    logger.warning(f"BLAST status for {rid} is unknown. Breaking loop. Response: {status_response.text[:200]}")
                    break # Unknown state, break loop
        else:
            # This 'else' belongs to the 'for' loop and executes if the loop finishes without a 'break'.
            logger.error(f"BLAST search for RID {rid} timed out after {max_retries * 5} seconds.")
            return {"error": "BLAST search timed out."}

        # Step 3: Retrieve the final results if we broke out of the loop (for READY text status)
        logger.info(f"Retrieving final results for RID: {rid}")
        final_result = requests.post(blast_url, data=status_payload, verify=False).json()
        return final_result

    except requests.exceptions.RequestException as e:
        logger.error(f"An error occurred with the BLAST API request: {e}")
        return {"error": str(e)}
    except Exception as e:
        logger.error(f"An unexpected error occurred during off-target analysis: {e}")
        return {"error": str(e)}

def find_intelligent_guides(locus: str, genome: str, pam: str = "NGG") -> list[dict]:
    """
    Main controller function for the Intelligent Guide Designer.

    This function orchestrates the entire process:
    1. Fetches the DNA sequence for a given locus.
    2. Finds all potential guide RNA candidates.
    3. Scores each candidate for on-target efficacy (functional importance).
    4. Analyzes each candidate for potential off-target effects.
    5. Returns a ranked list of guides with their scores.

    Args:
        locus: The genomic locus to target.
        genome: The genome build (e.g., 'hg38').
        pam: The PAM sequence to use for finding guides.

    Returns:
        A list of dictionaries, where each dictionary represents a guide
        and contains its sequence, on-target score, and off-target summary.
        Returns an empty list if the process fails at any critical step.
    """
    # Step 1: Get Sequence
    logger.info(f"Starting intelligent guide search for locus: {locus}")
    sequence = get_sequence_for_locus(locus, genome)
    if not sequence:
        logger.error("Guide search failed: could not retrieve sequence.")
        return []

    # Step 2: Find Candidates
    candidates = find_guide_candidates(sequence, pam)
    if not candidates:
        logger.warning("No guide candidates found for the given sequence.")
        return []
    
    logger.info(f"Found {len(candidates)} candidates. Now scoring...")
    
    results = []
    # Process a limited number of candidates for demonstration purposes
    for i, guide in enumerate(candidates[:5]): 
        logger.info(f"Processing candidate {i+1}/{len(candidates[:5])}: {guide}")
        
        # Step 3: Score On-Target Efficacy
        on_target_score = score_sequence_importance(guide)
        
        # Step 4: Analyze Off-Target Effects
        off_target_data = find_potential_off_targets(guide)
        
        num_off_targets = 0
        if "error" not in off_target_data:
            try:
                # Count the number of significant hits as a simple off-target metric
                num_off_targets = len(off_target_data["BlastOutput2"][0]["report"]["results"]["search"]["hits"])
            except (KeyError, IndexError):
                logger.warning(f"Could not parse BLAST results for guide {guide}")
        
        results.append({
            "guide_sequence": guide,
            "on_target_score": on_target_score,
            "off_target_hits": num_off_targets
        })

    # Step 5: Rank results (e.g., high on-target score, low off-target hits)
    # Simple ranking: prioritize higher on-target score, then lower off-target hits.
    ranked_results = sorted(results, key=lambda x: (x['on_target_score'] is None, -x.get('on_target_score', 0), x['off_target_hits']))
    
    logger.info("Intelligent guide search complete.")
    return ranked_results

if __name__ == '__main__':
    # Example usage for testing
    # Test with a region
    braf_locus_region = "chr7:140,753,000-140,754,000"
    print(f"Testing with region: {braf_locus_region}")
    sequence = get_sequence_for_locus(braf_locus_region, "hg38")
    if sequence:
        print(f"  - Fetched sequence length: {len(sequence)}")
    else:
        print("  - Failed to fetch sequence.")

    print("-" * 20)

    # Test with a single point from our workflow
    braf_locus_point = "chr7:140753336"
    print(f"Testing with point: {braf_locus_point}")
    sequence = get_sequence_for_locus(braf_locus_point, "hg38", window_size=500)
    if sequence:
        print(f"  - Fetched sequence length: {len(sequence)}")
        print("  - This demonstrates fetching a window around a variant for context.")
        
        # Test finding guides in the fetched sequence
        guides = find_guide_candidates(sequence)
        print(f"  - Found {len(guides)} potential guides.")
        if guides:
            print(f"  - Example guide: {guides[0]}")

            # Test scoring the first guide
            # NOTE: This requires a running LLM service accessible via query_llm
            # You may need to configure the 'provider' argument.
            # print(f"  - Scoring example guide...")
            # score = score_sequence_importance(guides[0])
            # if score is not None:
            #     print(f"  - Importance Score: {score}")
            # else:
            #     print("  - Failed to score guide.")

            # Test Off-Target Analysis
            # NOTE: This is a network-intensive call and may take a minute
            print(f"  - Performing off-target analysis for example guide...")
            off_target_results = find_potential_off_targets(guides[0])
            if "error" in off_target_results:
                print(f"  - Off-target analysis failed: {off_target_results['error']}")
            else:
                # Process and display summary of results
                try:
                    hits = off_target_results["BlastOutput2"][0]["report"]["results"]["search"]["hits"]
                    print(f"  - Found {len(hits)} potential off-target clusters.")
                    if hits:
                        top_hit = hits[0]
                        print(f"    - Top hit description: {top_hit['description'][0]['title']}")
                        print(f"    - Alignment E-value: {top_hit['hsps'][0]['evalue']}")
                except (KeyError, IndexError):
                    print("  - Could not parse off-target results format.")

    else:
        print("  - Failed to fetch sequence.")

    print("\\n" + "="*30 + "\\n")

    # --- Test the main controller ---
    print("Testing the main `find_intelligent_guides` controller function...")
    final_guides = find_intelligent_guides(braf_locus_point, "hg38")
    if final_guides:
        print(f"Successfully generated {len(final_guides)} ranked guides.")
        for guide in final_guides:
            print(f"  - Guide: {guide['guide_sequence']}, On-Target Score: {guide['on_target_score']}, Off-Target Hits: {guide['off_target_hits']}")
    else:
        print("The main controller function did not return any guides.") 