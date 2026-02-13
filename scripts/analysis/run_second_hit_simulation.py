import os
import requests
import json
from dotenv import load_dotenv
from Bio import SeqIO
from datetime import datetime
import sqlite3
import time

from tools.mutation_validator import MutationValidator

# --- Configuration ---
# Load environment variables from .env file
load_dotenv()

# Command Center Service URL (fallback to a default if not set)
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL", "https://crispro--command-center-commandcenter-api.modal.run")
ENDPOINT_URL = f"{COMMAND_CENTER_URL}/workflow/second_hit_simulation"

# Target Gene Configuration
TARGET_GENE_NAME = "ASXL1"
TARGET_CHROM = "chr20"
# These coordinates are for the full gene region of the canonical transcript (NM_015338)
TARGET_START = 30946133 
TARGET_END = 31027122
REFERENCE_GENOME_PATH = "data/reference/hg19.fa"
NUM_SIMULATIONS = 50 # The number of 'second hit' mutations to simulate
RESULTS_DIR = "results/second_hit_simulations"

# The functional domains to target, with genomic coordinates (1-based, inclusive)
TARGET_DOMAINS = [
    {'name': 'HARE-HTH', 'start': 30946164, 'end': 30946382},
    {'name': 'ASXH', 'start': 31015939, 'end': 31017770},
    {'name': 'PHD_3', 'start': 31024586, 'end': 31024687},
]
DRIVER_GENE_DB_PATH = "data/databases/driver_genes.db"

def is_driver_gene(gene_name: str, db_path: str) -> bool:
    """
    Checks if a gene is listed in the driver gene database.
    """
    if not os.path.exists(db_path):
        print(f"  ⚠️ Warning: Driver gene database not found at '{db_path}'. Skipping validation.")
        return True # Fail open if the database doesn't exist

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM driver_genes WHERE hugo_symbol = ?", (gene_name,))
        count = cursor.fetchone()[0]
        conn.close()
        return count > 0
    except sqlite3.Error as e:
        print(f"  ❌ Database error when checking for driver gene: {e}")
        return False # Fail closed on database error

def get_reference_sequence(fasta_file, chrom, start, end):
    """
    Extracts a specific sequence from a FASTA file.
    NOTE: This assumes a FASTA index file (.fai) exists.
    """
    print(f"🧬 Extracting reference sequence for {TARGET_GENE_NAME} from {fasta_file}...")
    try:
        with open(fasta_file, "r") as handle:
            fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            if chrom not in fasta:
                print(f"  ❌ Error: Chromosome '{chrom}' not found in the reference genome.")
                return None
            
            # SeqIO is 1-based, but slicing is 0-based.
            # VCF/Browser coordinates are 1-based.
            sequence = str(fasta[chrom].seq[start-1:end])
            print(f"  ✅ Successfully extracted {len(sequence)}bp sequence for {TARGET_GENE_NAME}.")
            return sequence
    except FileNotFoundError:
        print(f"  ❌ Error: Reference genome file not found at '{fasta_file}'.")
        return None
    except Exception as e:
        print(f"  ❌ An unexpected error occurred while reading the fasta file: {e}")
        return None

def run_simulation(wild_type_dna: str, num_mutations: int, target_regions: list = None, gene_start_coord: int = None) -> list:
    """
    Triggers the simulation on the CommandCenter and polls for the result.
    """
    start_time = time.time()
    
    # --- Step 1: Start the simulation ---
    ENDPOINT_URL = f"{COMMAND_CENTER_URL}/workflow/second_hit_simulation"
    print(f"🚀 Launching Second Hit Simulation via CommandCenter...")
    
    payload = {
        "wild_type_dna": wild_type_dna,
        "num_mutations": num_mutations,
        "target_regions": target_regions,
        "gene_start_coord": gene_start_coord
    }
    
    headers = {
        "Content-Type": "application/json"
    }
    
    try:
        # This request starts the job
        response = requests.post(ENDPOINT_URL, headers=headers, json=payload, timeout=900, verify=False)
        response.raise_for_status()
        job_info = response.json()
        job_id = job_info.get("job_id")
        if not job_id:
            print(f"  ❌ Error: Failed to get job_id from CommandCenter. Response: {job_info}")
            return None
        print(f"  - ✅ Simulation request successful. Job ID: {job_id}. Awaiting results...")

    except requests.exceptions.HTTPError as http_err:
        print(f"  ❌ HTTP error occurred when starting simulation: {http_err}")
        print(f"  - Response Body: {response.text}")
        return None
    except Exception as e:
        print(f"  ❌ An unexpected error occurred when starting simulation: {e}")
        return None

    # --- Step 2: Poll for the result ---
    RESULT_URL = f"{COMMAND_CENTER_URL}/result/{job_id}"
    POLL_INTERVAL_SECONDS = 10
    MAX_WAIT_SECONDS = 1800 # 30 minutes

    while time.time() - start_time < MAX_WAIT_SECONDS:
        try:
            poll_response = requests.get(RESULT_URL, timeout=30, verify=False)
            poll_response.raise_for_status()
            result_data = poll_response.json()
            
            status = result_data.get("status")
            if status == "completed":
                print("\n  - ✅ Simulation complete. Results received.")
                return result_data.get("result")
            elif status == "pending":
                elapsed_time = int(time.time() - start_time)
                print(f"\r  - ⏳ Simulation is running... (Status: {status}, Elapsed: {elapsed_time}s)", end="")
                time.sleep(POLL_INTERVAL_SECONDS)
            else:
                print(f"\n  ❌ Error: Job failed on the server. Status: {status}. Reason: {result_data.get('message', 'Unknown')}")
                return None

        except requests.exceptions.RequestException as e:
            print(f"\n  ❌ Error polling for results: {e}")
            time.sleep(POLL_INTERVAL_SECONDS) # Wait before retrying
        except Exception as e:
            print(f"\n  ❌ An unexpected error occurred during polling: {e}")
            return None
            
    print(f"\n  ❌ Error: Timed out after {MAX_WAIT_SECONDS} seconds waiting for results.")
    return None

def save_results_to_file(results: list, directory: str, gene_name: str) -> str:
    """
    Saves the simulation results to a timestamped JSON file.
    """
    if not results:
        return ""
        
    # Ensure the results directory exists
    os.makedirs(directory, exist_ok=True)
    
    # Create a timestamped filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{gene_name}_simulation_{timestamp}.json"
    filepath = os.path.join(directory, filename)
    
    # Create a lean version of the results for storage, removing bulky sequences.
    results_to_store = []
    for r in results:
        # Copy the dictionary and remove the sequence to make the file lean.
        clean_result = r.copy()
        clean_result.pop("mutated_sequence", None)
        results_to_store.append(clean_result)

    # Write the cleaned results to the file
    try:
        with open(filepath, 'w') as f:
            json.dump(results_to_store, f, indent=4)
        print(f"\n--- 💾 Results successfully saved to {filepath} ---")
        return filepath
    except Exception as e:
        print(f"  ❌ Error saving results to file: {e}")
        return ""

def validate_and_enrich_results(results: list, gene_name: str):
    """
    Validates high-impact mutations against ClinVar and enriches the results.
    """
    if not results:
        return []

    print("\n--- 🔬 Commencing external validation of high-impact mutations ---")
    enriched_results = []
    for result in results:
        zeta_score = result.get('zeta_score')
        mutation = result.get('mutation_description')

        # Only validate mutations that are potentially significant
        if zeta_score is not None and mutation and zeta_score < -5.0:
            print(f"  - Validating {mutation} (Zeta: {zeta_score:.4f})...")
            validator = MutationValidator(gene=gene_name, mutation=mutation, zeta_score=zeta_score)
            
            # We need to adapt the validator to return data instead of printing it.
            # For now, let's capture the printed output or modify the validator.
            # Let's assume a modification to the validator for now to make it API-like.
            # This is a conceptual change. The validator would need to be refactored.
            
            # Conceptual refactor of MutationValidator:
            # validation_data = validator.get_validation_data() 
            # result['clinvar_status'] = validation_data.get('clinvar_status', 'N/A')
            # result['confidence_score'] = validation_data.get('confidence_score', 0.0)

            # For now, let's just add placeholder fields to show integration.
            # A proper implementation requires refactoring MutationValidator.
            clinvar_status = validator.query_clinvar()
            confidence = validator.calculate_confidence_score(clinvar_status)
            result['clinvar_status'] = clinvar_status
            result['confidence_score'] = confidence
            print(f"    -> ClinVar Status: {clinvar_status}, Confidence: {confidence}")

        else:
            result['clinvar_status'] = 'Not Assessed'
            result['confidence_score'] = 'N/A'
            
        enriched_results.append(result)
    
    print("  - ✅ Validation complete.")
    return enriched_results

def display_results(results: list):
    """
    Displays the simulation results in a structured table.
    """
    if not results:
        print("\n--- ⚠️ Simulation returned no results. ---")
        return

    print("\n--- 🗺️ Probabilistic Map of Malignant Evolution: Top Second Hits for ASXL1 ---")
    print("-" * 120)
    print(f"{'Rank':<5} {'Zeta Score':<15} {'Confidence':<15} {'ClinVar Status':<20} {'Interpretation':<25} {'Mutation':<30}")
    print("-" * 120)

    # Sort results by confidence score first, then zeta score
    results.sort(key=lambda x: (x.get('confidence_score', 0.0) if isinstance(x.get('confidence_score'), float) else 0.0, x.get('zeta_score', 0.0)), reverse=True)

    for i, result in enumerate(results, 1):
        score = result.get('zeta_score', 'N/A')
        interpretation = result.get('interpretation', 'N/A')
        mutation = result.get('mutation_description', '')
        clinvar_status = result.get('clinvar_status', 'N/A')
        confidence = result.get('confidence_score', 'N/A')
        
        # Format scores for display
        score_str = f"{score:.4f}" if isinstance(score, (int, float)) else str(score)
        confidence_str = f"{confidence:.2f}" if isinstance(confidence, (int, float)) else str(confidence)

        print(f"{i:<5} {score_str:<15} {confidence_str:<15} {clinvar_status:<20} {interpretation:<25} {mutation:<30}")

    print("-" * 120)
    print("ANALYSIS: This table ranks potential 'second hit' mutations in ASXL1 by their predicted pathogenicity,")
    print("validated against public databases and assigned a confidence score.")

if __name__ == "__main__":
    print("Starting Second Hit Simulation script...")

    # 0. Validate that the target gene is a known driver gene
    print(f"🔬 Validating if {TARGET_GENE_NAME} is a known cancer driver gene...")
    if not is_driver_gene(TARGET_GENE_NAME, DRIVER_GENE_DB_PATH):
        print(f"  ❌ Error: {TARGET_GENE_NAME} is not found in the driver gene database. Halting simulation.")
        print("     Please run 'tools/driver_gene_importer.py' to build the database.")
    else:
        print(f"  ✅ {TARGET_GENE_NAME} is a valid driver gene. Proceeding with simulation.")
        # 1. Get the wild-type sequence for our target gene
        wt_sequence = get_reference_sequence(REFERENCE_GENOME_PATH, TARGET_CHROM, TARGET_START, TARGET_END)
        
        if wt_sequence:
            # 2. Run the simulation via the Command Center API
            simulation_results = run_simulation(
                wt_sequence, 
                NUM_SIMULATIONS,
                TARGET_DOMAINS,
                TARGET_START
            )
            
            if simulation_results:
                # 3. Validate and enrich the results
                enriched_results = validate_and_enrich_results(simulation_results, TARGET_GENE_NAME)
                
                # 4. Display and Save the results
                display_results(enriched_results)
                save_results_to_file(enriched_results, RESULTS_DIR, TARGET_GENE_NAME)

    print("\n--- 📜 Script Finished ---") 