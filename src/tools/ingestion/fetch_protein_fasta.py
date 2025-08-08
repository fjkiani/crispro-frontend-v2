import requests
import argparse
import os
import time
import xml.etree.ElementTree as ET

# --- NCBI API Configuration ---
EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API_KEY = os.getenv("NCBI_API_KEY") # It's good practice to use an API key

def fetch_canonical_protein_sequence(gene_symbol: str) -> str:
    """
    Fetches the canonical protein sequence for a gene from NCBI.
    """
    print(f"🔎 Searching for canonical protein for '{gene_symbol}'...")
    
    # 1. Search for the gene to get its Gene ID
    search_params = {
        "db": "gene",
        "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism] AND alive[prop]",
        "retmode": "json",
    }
    if API_KEY:
        search_params["api_key"] = API_KEY

    try:
        response = requests.get(f"{EUTILS_BASE_URL}esearch.fcgi", params=search_params)
        response.raise_for_status()
        search_data = response.json()
        gene_ids = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not gene_ids:
            print(f"❌ Could not find Gene ID for '{gene_symbol}'. Skipping.")
            return None
        
        gene_id = gene_ids[0]
        print(f"  -> Found Gene ID: {gene_id}")
        time.sleep(1) # Respect NCBI rate limits

        # 2. Fetch the gene summary
        summary_params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "xml",
        }
        if API_KEY:
            summary_params["api_key"] = API_KEY

        response = requests.get(f"{EUTILS_BASE_URL}efetch.fcgi", params=summary_params)
        response.raise_for_status()
        
        root = ET.fromstring(response.content)
        
        protein_accession = None
        
        # Robust manual iteration to find the protein accession
        for gene_commentary in root.iter("Gene-commentary"):
            commentary_type_node = gene_commentary.find("Gene-commentary_type")
            if commentary_type_node is not None and commentary_type_node.get("value") == "genomic_info":
                products = gene_commentary.find("Gene-commentary_products")
                if products is not None:
                    for product in products.iter("Gene-commentary"):
                        accession_node = product.find("Gene-commentary_accession")
                        if accession_node is not None and accession_node.text and accession_node.text.startswith("NP_"):
                            protein_accession = accession_node.text
                            break
                if protein_accession:
                    break
        
        if not protein_accession:
            for elem in root.iter():
                if elem.text and elem.text.strip().startswith("NP_"):
                    protein_accession = elem.text.strip()
                    break

        if not protein_accession:
            print(f"❌ Could not find ANY protein accession (NP_...) for '{gene_symbol}'.")
            return None
            
        print(f"  -> Found Protein Accession: {protein_accession}")
        time.sleep(1)

        # 3. Fetch the FASTA sequence
        fasta_params = {
            "db": "protein",
            "id": protein_accession,
            "rettype": "fasta",
            "retmode": "text",
        }
        if API_KEY:
            fasta_params["api_key"] = API_KEY

        response = requests.get(f"{EUTILS_BASE_URL}efetch.fcgi", params=fasta_params)
        response.raise_for_status()
        
        print(f"✅ Successfully fetched FASTA for '{gene_symbol}' ({protein_accession})")
        return response.text

    except requests.exceptions.RequestException as e:
        print(f"HTTP Error fetching data for '{gene_symbol}': {e}")
        return None
    except ET.ParseError as e:
        print(f"XML Parse Error for '{gene_symbol}': {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred for '{gene_symbol}': {e}")
        return None

def save_sequence_to_fasta(gene_symbol: str, sequence_data: str, output_dir: str):
    """Saves the sequence data to a FASTA file."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📂 Created output directory: {output_dir}")
        
    file_path = os.path.join(output_dir, f"{gene_symbol}.fasta")
    with open(file_path, 'w') as f:
        f.write(sequence_data)
    print(f"💾 Saved sequence to {file_path}")

def main():
    parser = argparse.ArgumentParser(description="Fetch canonical protein sequences from NCBI.")
    parser.add_argument("genes", nargs='+', help="A list of gene symbols to fetch.")
    parser.add_argument("--output-dir", default="data/protein_references", help="Directory to save the FASTA files.")
    
    args = parser.parse_args()
    
    print(f"--- Starting Protein Ingestion for {len(args.genes)} genes ---")
    
    for gene in args.genes:
        output_file = os.path.join(args.output_dir, f"{gene}.fasta")
        if os.path.exists(output_file):
            print(f"✅ File for '{gene}' already exists at {output_file}. Skipping.")
            continue
            
        sequence = fetch_canonical_protein_sequence(gene)
        if sequence:
            save_sequence_to_fasta(gene, sequence, args.output_dir)
        print("-" * 20)
        time.sleep(1) # Be extra cautious with rate limits

    print("--- Protein Ingestion Complete ---")

if __name__ == "__main__":
    main() 