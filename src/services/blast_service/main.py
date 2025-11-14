import modal
import os
import subprocess
from pathlib import Path
import shutil
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

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
    .pip_install("fastapi", "pydantic")
    .apt_install("ncbi-blast+", "wget", "minimap2", "samtools")
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

# --- Pydantic Models for API Input ---
class BlastRequest(BaseModel):
    query_sequence: str

class OffTargetRequest(BaseModel):
    guide_sequence: str
    mismatch_tolerance: int = 3
    max_hits: int = 100

# --- FastAPI App ---
fastapi_app = FastAPI(title="BLAST Service")

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
            process_result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
            
            # Read the results from the output file
            with open(output_file, "r") as f:
                results_xml = f.read()

            print("✅ BLAST search completed successfully.")
            return {"raw_blast_xml": results_xml}
        except subprocess.TimeoutExpired as e:
            print(f"❌ BLAST search timed out. Stderr: {e.stderr}")
            return {"error": "BLAST search process timed out.", "details": e.stderr}
        except subprocess.CalledProcessError as e:
            print(f"❌ BLAST search failed with exit code {e.returncode}.")
            print(f"  - stdout: {e.stdout}")
            print(f"  - stderr: {e.stderr}")
            return {"error": "BLAST process failed.", "details": e.stderr}
        except Exception as e:
            print(f"❌ An unexpected error occurred during BLAST search: {e}")
            return {"error": "An unexpected error occurred.", "details": str(e)}

    def search_offtargets(self, guide_sequence: str, mismatch_tolerance: int = 3, max_hits: int = 100) -> dict:
        """
        Search for off-target sites using minimap2 (fast) with BLAST fallback.
        
        Returns:
            {
                "hits": [{"chrom": str, "pos": int, "strand": str, "mismatches": int, "score": float}],
                "total_hits": int,
                "method": "minimap2" | "blast" | "error",
                "provenance": {...}
            }
        """
        print(f"🎯 Searching off-targets for guide: {guide_sequence} (tolerance: {mismatch_tolerance})")
        
        # Try minimap2 first (much faster than BLAST)
        try:
            return self._search_with_minimap2(guide_sequence, mismatch_tolerance, max_hits)
        except Exception as e:
            print(f"⚠️ minimap2 failed: {e}. Falling back to BLAST...")
            # Fallback to BLAST if minimap2 fails
            try:
                return self._search_with_blast_offtarget(guide_sequence, mismatch_tolerance, max_hits)
            except Exception as e2:
                print(f"❌ Both minimap2 and BLAST failed: {e2}")
                return {
                    "hits": [],
                    "total_hits": 0,
                    "method": "error",
                    "provenance": {"error": str(e2)}
                }
    
    def _search_with_minimap2(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Use minimap2 for fast off-target search"""
        # Create query file
        query_file = BLAST_DB_PATH / "offtarget_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        # Use the reference FASTA (not BLAST DB)
        ref_fasta = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem
        
        if not ref_fasta.exists():
            raise FileNotFoundError(f"Reference genome not found at {ref_fasta}")
        
        # Run minimap2 with parameters tuned for short sequences
        cmd = [
            "minimap2",
            "-a",  # Output SAM format
            "-x", "sr",  # Short read preset
            "-N", str(max_hits),  # Max alignments
            str(ref_fasta),
            str(query_file)
        ]
        
        print(f"minimap2 command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)
        
        # Parse SAM output
        hits = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):  # Skip header
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            
            # Extract alignment info
            chrom = fields[2]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            
            # Estimate mismatches from CIGAR and NM tag
            mismatches = 0
            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    mismatches = int(field.split(":")[-1])
                    break
            
            if mismatches <= mismatch_tolerance:
                hits.append({
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "mismatches": mismatches,
                    "score": 1.0 - (mismatches / len(guide_sequence))
                })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "minimap2",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }
    
    def _search_with_blast_offtarget(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Fallback to BLAST for off-target search"""
        # Similar to existing search() but with specific parameters for off-target detection
        query_file = BLAST_DB_PATH / "offtarget_blast_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        output_file = BLAST_DB_PATH / "offtarget_results.txt"
        
        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "6 sseqid sstart send sstrand mismatch",  # Tabular format
            "-task", "blastn-short",
            "-word_size", "7",
            "-max_target_seqs", str(max_hits),
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        
        # Parse BLAST output
        hits = []
        with open(output_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                
                chrom, start, end, strand, mismatches = fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])
                
                if mismatches <= mismatch_tolerance:
                    hits.append({
                        "chrom": chrom,
                        "pos": start,
                        "strand": strand,
                        "mismatches": mismatches,
                        "score": 1.0 - (mismatches / len(guide_sequence))
                    })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "blast",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }

    # This makes the class callable as a web endpoint
    @modal.asgi_app()
    def web_app(self):
        return fastapi_app

@fastapi_app.post("/search")
async def run_blast_search(request: BlastRequest):
    """Web endpoint to run a BLAST search."""
    blast_service = BlastService() # Create an instance
    # Note: This is a blocking call within an async endpoint.
    # For a real high-throughput service, we'd use .spawn() or .call_later()
    # But for this use case, it's acceptable.
    results = blast_service.search(request.query_sequence)
    if "error" in results:
        raise HTTPException(status_code=500, detail=results)
    return results

@fastapi_app.post("/offtarget_search")
async def run_offtarget_search(request: OffTargetRequest):
    """
    Web endpoint to search for off-target sites using minimap2 (fast) with BLAST fallback.
    
    Returns list of genomic locations where the guide RNA may bind with up to N mismatches.
    """
    blast_service = BlastService()
    results = blast_service.search_offtargets(
        request.guide_sequence,
        request.mismatch_tolerance,
        request.max_hits
    )
    
    if results.get("method") == "error":
        raise HTTPException(status_code=500, detail=results.get("provenance", {}))
    
    return results

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
            f.write(f">guide\n{guide_sequence}\n")
        
        # Use the reference FASTA (not BLAST DB)
        ref_fasta = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem
        
        if not ref_fasta.exists():
            raise FileNotFoundError(f"Reference genome not found at {ref_fasta}")
        
        # Run minimap2 with parameters tuned for short sequences
        cmd = [
            "minimap2",
            "-a",  # Output SAM format
            "-x", "sr",  # Short read preset
            "-N", str(max_hits),  # Max alignments
            str(ref_fasta),
            str(query_file)
        ]
        
        print(f"minimap2 command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)
        
        # Parse SAM output
        hits = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):  # Skip header
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            
            # Extract alignment info
            chrom = fields[2]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            
            # Estimate mismatches from CIGAR and NM tag
            mismatches = 0
            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    mismatches = int(field.split(":")[-1])
                    break
            
            if mismatches <= mismatch_tolerance:
                hits.append({
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "mismatches": mismatches,
                    "score": 1.0 - (mismatches / len(guide_sequence))
                })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "minimap2",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }
    
    def _search_with_blast_offtarget(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Fallback to BLAST for off-target search"""
        # Similar to existing search() but with specific parameters for off-target detection
        query_file = BLAST_DB_PATH / "offtarget_blast_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        output_file = BLAST_DB_PATH / "offtarget_results.txt"
        
        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "6 sseqid sstart send sstrand mismatch",  # Tabular format
            "-task", "blastn-short",
            "-word_size", "7",
            "-max_target_seqs", str(max_hits),
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        
        # Parse BLAST output
        hits = []
        with open(output_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                
                chrom, start, end, strand, mismatches = fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])
                
                if mismatches <= mismatch_tolerance:
                    hits.append({
                        "chrom": chrom,
                        "pos": start,
                        "strand": strand,
                        "mismatches": mismatches,
                        "score": 1.0 - (mismatches / len(guide_sequence))
                    })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "blast",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }

    # This makes the class callable as a web endpoint
    @modal.asgi_app()
    def web_app(self):
        return fastapi_app

@fastapi_app.post("/search")
async def run_blast_search(request: BlastRequest):
    """Web endpoint to run a BLAST search."""
    blast_service = BlastService() # Create an instance
    # Note: This is a blocking call within an async endpoint.
    # For a real high-throughput service, we'd use .spawn() or .call_later()
    # But for this use case, it's acceptable.
    results = blast_service.search(request.query_sequence)
    if "error" in results:
        raise HTTPException(status_code=500, detail=results)
    return results

@fastapi_app.post("/offtarget_search")
async def run_offtarget_search(request: OffTargetRequest):
    """
    Web endpoint to search for off-target sites using minimap2 (fast) with BLAST fallback.
    
    Returns list of genomic locations where the guide RNA may bind with up to N mismatches.
    """
    blast_service = BlastService()
    results = blast_service.search_offtargets(
        request.guide_sequence,
        request.mismatch_tolerance,
        request.max_hits
    )
    
    if results.get("method") == "error":
        raise HTTPException(status_code=500, detail=results.get("provenance", {}))
    
    return results

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
            f.write(f">guide\n{guide_sequence}\n")
        
        # Use the reference FASTA (not BLAST DB)
        ref_fasta = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem
        
        if not ref_fasta.exists():
            raise FileNotFoundError(f"Reference genome not found at {ref_fasta}")
        
        # Run minimap2 with parameters tuned for short sequences
        cmd = [
            "minimap2",
            "-a",  # Output SAM format
            "-x", "sr",  # Short read preset
            "-N", str(max_hits),  # Max alignments
            str(ref_fasta),
            str(query_file)
        ]
        
        print(f"minimap2 command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)
        
        # Parse SAM output
        hits = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):  # Skip header
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            
            # Extract alignment info
            chrom = fields[2]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            
            # Estimate mismatches from CIGAR and NM tag
            mismatches = 0
            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    mismatches = int(field.split(":")[-1])
                    break
            
            if mismatches <= mismatch_tolerance:
                hits.append({
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "mismatches": mismatches,
                    "score": 1.0 - (mismatches / len(guide_sequence))
                })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "minimap2",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }
    
    def _search_with_blast_offtarget(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Fallback to BLAST for off-target search"""
        # Similar to existing search() but with specific parameters for off-target detection
        query_file = BLAST_DB_PATH / "offtarget_blast_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        output_file = BLAST_DB_PATH / "offtarget_results.txt"
        
        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "6 sseqid sstart send sstrand mismatch",  # Tabular format
            "-task", "blastn-short",
            "-word_size", "7",
            "-max_target_seqs", str(max_hits),
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        
        # Parse BLAST output
        hits = []
        with open(output_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                
                chrom, start, end, strand, mismatches = fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])
                
                if mismatches <= mismatch_tolerance:
                    hits.append({
                        "chrom": chrom,
                        "pos": start,
                        "strand": strand,
                        "mismatches": mismatches,
                        "score": 1.0 - (mismatches / len(guide_sequence))
                    })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "blast",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }

    # This makes the class callable as a web endpoint
    @modal.asgi_app()
    def web_app(self):
        return fastapi_app

@fastapi_app.post("/search")
async def run_blast_search(request: BlastRequest):
    """Web endpoint to run a BLAST search."""
    blast_service = BlastService() # Create an instance
    # Note: This is a blocking call within an async endpoint.
    # For a real high-throughput service, we'd use .spawn() or .call_later()
    # But for this use case, it's acceptable.
    results = blast_service.search(request.query_sequence)
    if "error" in results:
        raise HTTPException(status_code=500, detail=results)
    return results

@fastapi_app.post("/offtarget_search")
async def run_offtarget_search(request: OffTargetRequest):
    """
    Web endpoint to search for off-target sites using minimap2 (fast) with BLAST fallback.
    
    Returns list of genomic locations where the guide RNA may bind with up to N mismatches.
    """
    blast_service = BlastService()
    results = blast_service.search_offtargets(
        request.guide_sequence,
        request.mismatch_tolerance,
        request.max_hits
    )
    
    if results.get("method") == "error":
        raise HTTPException(status_code=500, detail=results.get("provenance", {}))
    
    return results

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
            f.write(f">guide\n{guide_sequence}\n")
        
        # Use the reference FASTA (not BLAST DB)
        ref_fasta = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem
        
        if not ref_fasta.exists():
            raise FileNotFoundError(f"Reference genome not found at {ref_fasta}")
        
        # Run minimap2 with parameters tuned for short sequences
        cmd = [
            "minimap2",
            "-a",  # Output SAM format
            "-x", "sr",  # Short read preset
            "-N", str(max_hits),  # Max alignments
            str(ref_fasta),
            str(query_file)
        ]
        
        print(f"minimap2 command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)
        
        # Parse SAM output
        hits = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):  # Skip header
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            
            # Extract alignment info
            chrom = fields[2]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            
            # Estimate mismatches from CIGAR and NM tag
            mismatches = 0
            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    mismatches = int(field.split(":")[-1])
                    break
            
            if mismatches <= mismatch_tolerance:
                hits.append({
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "mismatches": mismatches,
                    "score": 1.0 - (mismatches / len(guide_sequence))
                })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "minimap2",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }
    
    def _search_with_blast_offtarget(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Fallback to BLAST for off-target search"""
        # Similar to existing search() but with specific parameters for off-target detection
        query_file = BLAST_DB_PATH / "offtarget_blast_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        output_file = BLAST_DB_PATH / "offtarget_results.txt"
        
        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "6 sseqid sstart send sstrand mismatch",  # Tabular format
            "-task", "blastn-short",
            "-word_size", "7",
            "-max_target_seqs", str(max_hits),
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        
        # Parse BLAST output
        hits = []
        with open(output_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                
                chrom, start, end, strand, mismatches = fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])
                
                if mismatches <= mismatch_tolerance:
                    hits.append({
                        "chrom": chrom,
                        "pos": start,
                        "strand": strand,
                        "mismatches": mismatches,
                        "score": 1.0 - (mismatches / len(guide_sequence))
                    })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "blast",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }

    # This makes the class callable as a web endpoint
    @modal.asgi_app()
    def web_app(self):
        return fastapi_app

@fastapi_app.post("/search")
async def run_blast_search(request: BlastRequest):
    """Web endpoint to run a BLAST search."""
    blast_service = BlastService() # Create an instance
    # Note: This is a blocking call within an async endpoint.
    # For a real high-throughput service, we'd use .spawn() or .call_later()
    # But for this use case, it's acceptable.
    results = blast_service.search(request.query_sequence)
    if "error" in results:
        raise HTTPException(status_code=500, detail=results)
    return results

@fastapi_app.post("/offtarget_search")
async def run_offtarget_search(request: OffTargetRequest):
    """
    Web endpoint to search for off-target sites using minimap2 (fast) with BLAST fallback.
    
    Returns list of genomic locations where the guide RNA may bind with up to N mismatches.
    """
    blast_service = BlastService()
    results = blast_service.search_offtargets(
        request.guide_sequence,
        request.mismatch_tolerance,
        request.max_hits
    )
    
    if results.get("method") == "error":
        raise HTTPException(status_code=500, detail=results.get("provenance", {}))
    
    return results

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
            f.write(f">guide\n{guide_sequence}\n")
        
        # Use the reference FASTA (not BLAST DB)
        ref_fasta = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem
        
        if not ref_fasta.exists():
            raise FileNotFoundError(f"Reference genome not found at {ref_fasta}")
        
        # Run minimap2 with parameters tuned for short sequences
        cmd = [
            "minimap2",
            "-a",  # Output SAM format
            "-x", "sr",  # Short read preset
            "-N", str(max_hits),  # Max alignments
            str(ref_fasta),
            str(query_file)
        ]
        
        print(f"minimap2 command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)
        
        # Parse SAM output
        hits = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):  # Skip header
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            
            # Extract alignment info
            chrom = fields[2]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            
            # Estimate mismatches from CIGAR and NM tag
            mismatches = 0
            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    mismatches = int(field.split(":")[-1])
                    break
            
            if mismatches <= mismatch_tolerance:
                hits.append({
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "mismatches": mismatches,
                    "score": 1.0 - (mismatches / len(guide_sequence))
                })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "minimap2",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }
    
    def _search_with_blast_offtarget(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Fallback to BLAST for off-target search"""
        # Similar to existing search() but with specific parameters for off-target detection
        query_file = BLAST_DB_PATH / "offtarget_blast_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        output_file = BLAST_DB_PATH / "offtarget_results.txt"
        
        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "6 sseqid sstart send sstrand mismatch",  # Tabular format
            "-task", "blastn-short",
            "-word_size", "7",
            "-max_target_seqs", str(max_hits),
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        
        # Parse BLAST output
        hits = []
        with open(output_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                
                chrom, start, end, strand, mismatches = fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])
                
                if mismatches <= mismatch_tolerance:
                    hits.append({
                        "chrom": chrom,
                        "pos": start,
                        "strand": strand,
                        "mismatches": mismatches,
                        "score": 1.0 - (mismatches / len(guide_sequence))
                    })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "blast",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }

    # This makes the class callable as a web endpoint
    @modal.asgi_app()
    def web_app(self):
        return fastapi_app

@fastapi_app.post("/search")
async def run_blast_search(request: BlastRequest):
    """Web endpoint to run a BLAST search."""
    blast_service = BlastService() # Create an instance
    # Note: This is a blocking call within an async endpoint.
    # For a real high-throughput service, we'd use .spawn() or .call_later()
    # But for this use case, it's acceptable.
    results = blast_service.search(request.query_sequence)
    if "error" in results:
        raise HTTPException(status_code=500, detail=results)
    return results

@fastapi_app.post("/offtarget_search")
async def run_offtarget_search(request: OffTargetRequest):
    """
    Web endpoint to search for off-target sites using minimap2 (fast) with BLAST fallback.
    
    Returns list of genomic locations where the guide RNA may bind with up to N mismatches.
    """
    blast_service = BlastService()
    results = blast_service.search_offtargets(
        request.guide_sequence,
        request.mismatch_tolerance,
        request.max_hits
    )
    
    if results.get("method") == "error":
        raise HTTPException(status_code=500, detail=results.get("provenance", {}))
    
    return results

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
            f.write(f">guide\n{guide_sequence}\n")
        
        # Use the reference FASTA (not BLAST DB)
        ref_fasta = BLAST_DB_PATH / Path(GENOME_FILE_NAME).stem
        
        if not ref_fasta.exists():
            raise FileNotFoundError(f"Reference genome not found at {ref_fasta}")
        
        # Run minimap2 with parameters tuned for short sequences
        cmd = [
            "minimap2",
            "-a",  # Output SAM format
            "-x", "sr",  # Short read preset
            "-N", str(max_hits),  # Max alignments
            str(ref_fasta),
            str(query_file)
        ]
        
        print(f"minimap2 command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)
        
        # Parse SAM output
        hits = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith("@"):  # Skip header
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            
            # Extract alignment info
            chrom = fields[2]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            
            # Estimate mismatches from CIGAR and NM tag
            mismatches = 0
            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    mismatches = int(field.split(":")[-1])
                    break
            
            if mismatches <= mismatch_tolerance:
                hits.append({
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "mismatches": mismatches,
                    "score": 1.0 - (mismatches / len(guide_sequence))
                })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "minimap2",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }
    
    def _search_with_blast_offtarget(self, guide_sequence: str, mismatch_tolerance: int, max_hits: int) -> dict:
        """Fallback to BLAST for off-target search"""
        # Similar to existing search() but with specific parameters for off-target detection
        query_file = BLAST_DB_PATH / "offtarget_blast_query.fa"
        with open(query_file, "w") as f:
            f.write(f">guide\n{guide_sequence}\n")
        
        output_file = BLAST_DB_PATH / "offtarget_results.txt"
        
        cmd = [
            "blastn",
            "-db", str(BLAST_DB_PATH / DB_NAME),
            "-query", str(query_file),
            "-out", str(output_file),
            "-outfmt", "6 sseqid sstart send sstrand mismatch",  # Tabular format
            "-task", "blastn-short",
            "-word_size", "7",
            "-max_target_seqs", str(max_hits),
        ]
        
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        
        # Parse BLAST output
        hits = []
        with open(output_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                
                chrom, start, end, strand, mismatches = fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])
                
                if mismatches <= mismatch_tolerance:
                    hits.append({
                        "chrom": chrom,
                        "pos": start,
                        "strand": strand,
                        "mismatches": mismatches,
                        "score": 1.0 - (mismatches / len(guide_sequence))
                    })
        
        return {
            "hits": hits[:max_hits],
            "total_hits": len(hits),
            "method": "blast",
            "provenance": {"mismatch_tolerance": mismatch_tolerance, "guide_length": len(guide_sequence)}
        }

    # This makes the class callable as a web endpoint
    @modal.asgi_app()
    def web_app(self):
        return fastapi_app

@fastapi_app.post("/search")
async def run_blast_search(request: BlastRequest):
    """Web endpoint to run a BLAST search."""
    blast_service = BlastService() # Create an instance
    # Note: This is a blocking call within an async endpoint.
    # For a real high-throughput service, we'd use .spawn() or .call_later()
    # But for this use case, it's acceptable.
    results = blast_service.search(request.query_sequence)
    if "error" in results:
        raise HTTPException(status_code=500, detail=results)
    return results

@fastapi_app.post("/offtarget_search")
async def run_offtarget_search(request: OffTargetRequest):
    """
    Web endpoint to search for off-target sites using minimap2 (fast) with BLAST fallback.
    
    Returns list of genomic locations where the guide RNA may bind with up to N mismatches.
    """
    blast_service = BlastService()
    results = blast_service.search_offtargets(
        request.guide_sequence,
        request.mismatch_tolerance,
        request.max_hits
    )
    
    if results.get("method") == "error":
        raise HTTPException(status_code=500, detail=results.get("provenance", {}))
    
    return results

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