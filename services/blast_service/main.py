import modal
import os
import subprocess
from pathlib import Path
import shutil

# --- Configuration ---
BLAST_DB_PATH = Path("/blast_db")
TMP_PATH = Path("/tmp")
# Switched from E.coli to Human Genome (Chromosome 21 for speed, full genome for production)
# Full Genome URL: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
GENOME_FASTA_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
GENOME_FILE_NAME = Path(GENOME_FASTA_URL).name
DB_NAME = "grch38"

# --- Modal Image Setup ---
# This image includes the NCBI BLAST+ toolkit and downloads the reference genome to a temporary path.
# The database itself is built on the first run and stored in the persistent volume.
blast_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("ncbi-blast+", "wget")
    .run_commands(
        f"echo 'Downloading Human reference genome to {TMP_PATH}...'",
        f"wget -q -O {TMP_PATH / GENOME_FILE_NAME} {GENOME_FASTA_URL}",
        "echo 'Download complete.'",
    )
)

app = modal.App("blast-service-human", image=blast_image)

# --- Persistent Volume for BLAST DB ---
# This volume stores the formatted BLAST database to avoid rebuilding it on every run.
blast_db_volume = modal.Volume.from_name("blast-db-volume-human-v2", create_if_missing=True)

@app.cls(
    volumes={str(BLAST_DB_PATH): blast_db_volume},
    cpu=4, 
    memory=8192, # Increased memory for human genome
    timeout=3600 # Increased timeout for db build and search
)
class BlastService:
    @modal.enter()
    def build_db_if_needed(self):
        """
        Builds the BLAST database on the first run and saves it to the persistent volume.
        Subsequent runs will be fast as they will find the pre-built database.
        """
        # Check if the database is already built by looking for a key file.
        db_marker_file = BLAST_DB_PATH / f"{DB_NAME}.nsq"
        if db_marker_file.exists():
            print("✅ Human BLAST database already exists in volume. Startup is fast.")
            return

        print("💥 Human BLAST database not found in volume. Building now... (This will take a while)")
        
        # Define paths
        source_gz_path = TMP_PATH / GENOME_FILE_NAME
        dest_gz_path = BLAST_DB_PATH / GENOME_FILE_NAME
        dest_fasta_path = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem

        try:
            # 1. Copy from the image's /tmp to the volume's /blast_db
            print(f"  - Copying genome from {source_gz_path} to {dest_gz_path}...")
            shutil.copy(str(source_gz_path), str(dest_gz_path))

            # 2. Decompress the file within the volume
            print(f"  - Decompressing {dest_gz_path}...")
            subprocess.run(["gunzip", str(dest_gz_path)], check=True, capture_output=True)
            
            # 3. Build the database from the file in the volume
            print(f"  - Building BLAST database: {DB_NAME}...")
            cmd = [
                "makeblastdb",
                "-in", str(dest_fasta_path),
                "-dbtype", "nucl",
                "-out", str(BLAST_DB_PATH / DB_NAME)
            ]
            subprocess.run(cmd, check=True, capture_output=True)

            # 4. Commit the changes to the persistent volume
            print("  - Committing database to persistent volume...")
            blast_db_volume.commit()
            
            print("✅ Human BLAST database built and saved successfully.")

        except FileNotFoundError:
            print(f"❌ CRITICAL: Genome file not found at {source_gz_path}. Image build may have failed.")
        except subprocess.CalledProcessError as e:
            print(f"❌ CRITICAL: Failed to build BLAST database. Subprocess error: {e.stderr.decode()}")
        except Exception as e:
            print(f"❌ CRITICAL: An unexpected error occurred during database setup: {e}")

    @modal.method()
    def search(self, query_sequence: str, num_mismatches: int = 3) -> dict:
        """
        Performs a BLAST search for a given sequence against the local database.
        """
        print(f"🔬 Running BLAST search for: {query_sequence}")
        
        # We need to write the query to a temporary file for blastn to use.
        query_file = BLAST_DB_PATH / "query.fa"
        with open(query_file, "w") as f:
            f.write(f">query\n{query_sequence}\n")

        output_file = BLAST_DB_PATH / "results.xml"

        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "5",  # XML format
            "-task", "blastn-short", # Optimized for short sequences
            "-word_size", "7",
            "-gapopen", "5",
            "-gapextend", "2",
            "-reward", "1",
            "-penalty", "-3",
        ]
        
        print(f"Blast command: {' '.join(cmd)}")
        
        try:
            # Execute the command
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
            
            # Read the results from the output file
            with open(output_file, "r") as f:
                results_xml = f.read()

            print("✅ BLAST search completed successfully.")
            return {"raw_blast_xml": results_xml}
        except subprocess.TimeoutExpired:
            print("❌ BLAST search timed out.")
            return {"error": "BLAST search process timed out."}
        except subprocess.CalledProcessError as e:
            print(f"❌ BLAST search failed with exit code {e.returncode}.")
            print(f"  - stderr: {e.stderr}")
            return {"error": "BLAST process failed.", "details": e.stderr}

def main():
    """A simple function to test the service locally."""
    # This requires `modal shell services/blast_service/main.py`
    # and running `main()` in the interactive shell.
    service = BlastService()
    service.build_db_if_needed()
    
    # Example gRNA sequence
    test_guide = "GATTACAGATTACAGATTAC"
    
    print(f"--- Testing BLAST search with guide: {test_guide} ---")
    results = service.search(test_guide)
    
    if "raw_blast_xml" in results and results["raw_blast_xml"]:
        print("✅ Test search successful. Received XML output.")
        # Pretty print first 200 chars of XML
        print(results["raw_blast_xml"][:200] + "...")
    else:
        print("❌ Test search failed.")
        print(results)

# --- API ---
# We can add a FastAPI wrapper here if we want to call it via HTTP.
# For inter-service calls within Modal, calling the class method directly is more efficient. 