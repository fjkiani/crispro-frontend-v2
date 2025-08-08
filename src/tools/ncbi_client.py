import requests
from loguru import logger
import re
from xml.etree import ElementTree
import os
import json
from Bio import Entrez, SeqIO
import time
from typing import Optional, Dict, List
from tenacity import retry, stop_after_attempt, wait_exponential
import xml.etree.ElementTree as ET

# --- DOCTRINE: PROVIDE A CLEAR USER AGENT ---
# This identifies our script to NCBI, which is a best practice.
Entrez.email = "crispr-assistant@example.com"
Entrez.api_key = os.environ.get("NCBI_API_KEY") # Use API key if available

class NCBIClient:
    """
    A client to interact with NCBI's E-utils for bioinformatics reconnaissance.
    It is initialized with a dynamic threat matrix to guide its prioritization.
    This client implements the Hydra Protocol for DNA sequence retrieval.
    """
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    LOCAL_GENE_DB_PATH = "data/gene_database"

    def __init__(self, threat_matrix: dict):
        """
        Initializes the client with a dynamic threat matrix.
        :param threat_matrix: A dictionary containing 'keywords' and 'cdd' lists.
        """
        self.high_value_keywords = threat_matrix.get("keywords", [])
        self.high_value_cdd = threat_matrix.get("cdd", {})
        self.threat_matrix = threat_matrix
        logger.info("Hydra Protocol Client armed. Mode: LOCAL_FIRST, NETWORK_FALLBACK.")

    def find_protein_domains(self, gene_symbol: str) -> list[dict]:
        """
        Finds functional domains for a given gene symbol by fetching its protein record.
        """
        logger.info(f"Initiating Subdomain Hunt for gene: {gene_symbol}")
        try:
            protein_id = self._get_protein_id(gene_symbol)
            if not protein_id:
                logger.warning(f"No protein ID found for {gene_symbol}. Hunt aborted.")
                return []

            return self._get_domains_from_protein_record(protein_id)
        except Exception as e:
            logger.error(f"Subdomain Hunt for {gene_symbol} failed: {e}")
            return []

    def _get_protein_id(self, gene_symbol: str) -> str | None:
        """
        Uses esearch to find the NCBI protein UID for a given gene symbol.
        """
        search_term = f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]"
        search_url = f"{self.BASE_URL}esearch.fcgi"
        params = {"db": "protein", "term": search_term, "retmode": "json", "sort": "relevance"}
        
        time.sleep(1) # Respect NCBI rate limits
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        result = response.json()
        
        id_list = result.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None
        
        logger.info(f"Found protein ID for {gene_symbol}: {id_list[0]}")
        return id_list[0]

    def _get_domains_from_protein_record(self, protein_id: str) -> list[dict]:
        """
        Uses efetch to get a protein record and parses it for domain information.
        This has been re-engineered for robustness against format changes.
        """
        fetch_url = f"{self.BASE_URL}efetch.fcgi"
        params = {"db": "protein", "id": protein_id, "rettype": "gb", "retmode": "xml"}
        
        logger.info(f"Fetching protein record for ID: {protein_id}")
        time.sleep(1) # Respect NCBI rate limits
        response = requests.get(fetch_url, params=params)
        response.raise_for_status()
        
        xml_root = ElementTree.fromstring(response.content)
        domains = []
        
        # --- DOCTRINE: RESILIENT PARSING ---
        # The NCBI XML structure can be inconsistent. We search for all possible
        # feature tags and then filter for the ones we need.
        for feature in xml_root.findall(".//GBFeature"):
            # Find the type of feature
            key_element = feature.find("GBFeature_key")
            if key_element is None or key_element.text not in ["Region", "Site"]:
                continue

            # Find the human-readable name of the region
            region_name = None
            value_text = None
            for qual in feature.findall(".//GBQualifier"):
                name = qual.find("GBQualifier_name")
                value = qual.find("GBQualifier_value")
                if value is not None:
                    value_text = value.text # Store the last seen value
                if name is not None and name.text == "region_name" and value is not None:
                    region_name = value.text
                # Also check for 'site_type' for Site features
                elif name is not None and name.text == "site_type" and value is not None:
                     region_name = value.text

            if not region_name:
                continue
            
            # Find the location string (e.g., "27..191")
            location = feature.find("GBFeature_location")
            if location is None or location.text is None:
                continue

            # Extract start and end coordinates from the location string
            match = re.search(r"(\d+)\.\.(\d+)", location.text)
            if match:
                start, end = match.groups()
                domain_name = value_text if value_text else region_name
                
                # --- DOCTRINE: THREAT PRIORITIZATION & ACTIONABLE INTELLIGENCE ---
                # Analyze the domain name for keywords and provide justification.
                priority = "Standard"
                justification = "Standard structural domain."

                # Extract the clean name and check for keywords or CDD codes
                clean_name = domain_name.split(';')[0].split('.')[0]
                
                # Check for any CDD first, as it's a strong indicator
                if "cdd:" in domain_name.lower():
                    priority = "High"
                    justification = f"High-value target: Conserved Domain ({clean_name})."

                matched_keyword = next((kw for kw in self.high_value_keywords if kw in clean_name.lower()), None)
                if matched_keyword:
                    priority = "High"
                    justification = f"High-value target: Contains keyword '{matched_keyword}'."
                
                # Allow specific, known CDDs to overwrite the generic justification
                matched_cdd_specific = next((cdd for cdd in self.high_value_cdd if cdd in domain_name.lower()), None)
                if matched_cdd_specific:
                    priority = "High"
                    justification = f"High-value target: {self.high_value_cdd[matched_cdd_specific]}."

                domain_info = {
                    "name": clean_name,
                    "start": int(start),
                    "end": int(end),
                    "priority": priority,
                    "justification": justification
                }
                logger.info(f"Found subdomain: {domain_info}")
                domains.append(domain_info)
        
        # Sort domains to show high-priority targets first
        domains.sort(key=lambda x: x['priority'] == 'High', reverse=True)
        
        if not domains:
            logger.warning(f"No domains parsed from protein ID {protein_id}. The data may lack 'Region' or 'Site' annotations.")
        
        return domains

    # --- HYDRA PROTOCOL IMPLEMENTATION ---

    def get_dna_sequence_from_gene_symbol(self, gene_symbol: str) -> Optional[str]:
        """
        Executes the Hydra Protocol to retrieve a DNA sequence.
        1. Checks the local armory.
        2. If not found, hunts on the NCBI network.
        3. Caches any network find to the local armory.
        """
        logger.debug(f"Executing Hydra Protocol for gene: {gene_symbol}")
        
        # Head 1: The Local Armory
        local_sequence = self._get_local_sequence(gene_symbol)
        if local_sequence:
            return local_sequence
        
        # Head 2: The Network Hunter
        logger.warning(f"'{gene_symbol}' not found in local armory. Engaging NCBI network hunter.")
        network_sequence = self._get_network_sequence(gene_symbol)
        
        if network_sequence:
            self._cache_sequence(gene_symbol, network_sequence)
            return network_sequence
        
        logger.error(f"Hydra Protocol failed for '{gene_symbol}'. Target not found in local armory or on network.")
        return None

    def _get_local_sequence(self, gene_symbol: str) -> Optional[str]:
        """
        (Hydra Head 1) Attempts to retrieve a DNA sequence from our local gene database.
        """
        possible_files = [
            f"{gene_symbol.upper()}.fasta", f"{gene_symbol.upper()}.fa",
            f"{gene_symbol.lower()}.fasta", f"{gene_symbol.lower()}.fa",
            f"{gene_symbol}.fasta", f"{gene_symbol}.fa"
        ]
        
        search_paths = [self.LOCAL_GENE_DB_PATH, os.path.join(self.LOCAL_GENE_DB_PATH, "reference")]
        
        for path in search_paths:
            if not os.path.exists(path): continue
            for filename in possible_files:
                filepath = os.path.join(path, filename)
                if os.path.exists(filepath):
                    logger.info(f"Local target acquired: {filepath}")
                    try:
                        with open(filepath, 'r') as f:
                            sequence_record = SeqIO.read(f, "fasta")
                            sequence = str(sequence_record.seq)
                            logger.success(f"Successfully loaded {gene_symbol} from local armory. Length: {len(sequence)}")
                            return sequence
                    except Exception as e:
                        logger.error(f"Error reading local file {filepath}: {e}")
        return None

    @retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=2, max=10))
    def _get_network_sequence(self, gene_symbol: str) -> Optional[str]:
        """
        (Hydra Head 2) Hardened network client to fetch sequence from NCBI.
        --- NEW DOCTRINE: DIRECT PATH ---
        This function bypasses the complex gene record and forges a direct
        path from gene symbol -> nuccore record -> FASTA sequence.
        """
        try:
            logger.debug(f"DIRECT PATH HUNT: Searching for Gene ID for '{gene_symbol}'.")
            
            # 1. ESearch: Find the Gene ID
            search_handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if not search_results.get("IdList"):
                logger.warning(f"No Gene ID found for '{gene_symbol}'.")
                return None
            
            gene_id = search_results["IdList"][0]
            logger.info(f"DIRECT PATH: Found Gene ID: {gene_id}")

            # 2. ELink: Find the linked Nucleotide (nuccore) ID
            # --- DOCTRINE: ADAPTIVE LINKING ---
            # First, try the most specific link. If it fails, broaden the search.
            link_names = ["gene_nuccore_refseq", "gene_nuccore"]
            linked_ids = []

            for link_name in link_names:
                logger.debug(f"Attempting to link with '{link_name}'...")
                link_handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname=link_name)
                link_results = Entrez.read(link_handle)
                link_handle.close()
                
                if not link_results or 'LinkSetDb' not in link_results[0] or not link_results[0]['LinkSetDb']:
                    continue # Try the next link name

                links = link_results[0]['LinkSetDb'][0].get('Link', [])
                if not links:
                    continue # Try the next link name
                
                linked_ids = [link['Id'] for link in links]
                if linked_ids:
                    logger.info(f"ADAPTIVE LINKING SUCCESS: Found link using '{link_name}'.")
                    break # Success, exit the loop
            
            if not linked_ids:
                logger.warning(f"No linked nuccore record ID found for Gene ID {gene_id} after all attempts.")
                return None
            
            # We will use the first linked ID, which is typically the canonical sequence.
            nuccore_id = linked_ids[0]
            logger.info(f"DIRECT PATH: Found linked Nucleotide ID: {nuccore_id}")

            # 3. EFetch: Retrieve the FASTA sequence directly
            fetch_handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="fasta", retmode="text")
            fasta_data = fetch_handle.read()
            fetch_handle.close()

            # Parse the FASTA data to get only the sequence
            sequence = "".join(fasta_data.splitlines()[1:])
            
            if not sequence:
                logger.error(f"Fetched FASTA for {nuccore_id} but it was empty.")
                return None
            
            logger.success(f"DIRECT PATH HUNT SUCCESS: Acquired sequence for {gene_symbol} (Length: {len(sequence)})")
            return sequence

        except Exception as e:
            logger.error(f"DIRECT PATH HUNT FAILED for '{gene_symbol}': {e}", exc_info=True)
            return None

    def _cache_sequence(self, gene_symbol: str, sequence: str):
        """Saves a successfully retrieved sequence to the local cache."""
        try:
            filepath = os.path.join(self.LOCAL_GENE_DB_PATH, f"{gene_symbol.upper()}.fasta")
            os.makedirs(self.LOCAL_GENE_DB_PATH, exist_ok=True)
            with open(filepath, "w") as f:
                f.write(f">{gene_symbol.upper()} sequence retrieved from NCBI\n")
                f.write(sequence)
            logger.info(f"Target cached: Saved '{gene_symbol}' sequence to {filepath}")
        except Exception as e:
            logger.error(f"Failed to cache sequence for '{gene_symbol}': {e}")
            
    def get_protein_sequence_from_gene_symbol(self, gene_symbol: str) -> Optional[str]:
        """
        --- HYDRA PROTOCOL V2 ---
        Fetches the canonical protein sequence for a given gene symbol, now with local-first caching.
        """
        logger.info(f"Predator Protocol: Hunting for protein sequence of '{gene_symbol}'")
        
        # Head 1: The Local Armory (Protein Division)
        local_sequence = self._get_local_protein_sequence(gene_symbol)
        if local_sequence:
            return local_sequence
            
        # Head 2: The Network Hunter (Protein Division)
        logger.warning(f"Protein for '{gene_symbol}' not in local armory. Engaging NCBI network hunter.")
        network_sequence = self._get_network_protein_sequence(gene_symbol)
        
        if network_sequence:
            self._cache_protein_sequence(gene_symbol, network_sequence)
            return network_sequence
            
        logger.error(f"Hydra Protocol failed for PROTEIN '{gene_symbol}'. Target not found.")
        return None

    def _get_local_protein_sequence(self, gene_symbol: str) -> Optional[str]:
        """(Hydra Head 1) Attempts to retrieve a PROTEIN sequence from the local database."""
        filepath = os.path.join(self.LOCAL_GENE_DB_PATH, f"{gene_symbol.upper()}_PROTEIN.fasta")
        if os.path.exists(filepath):
            logger.info(f"Local protein target acquired: {filepath}")
            try:
                with open(filepath, 'r') as f:
                    sequence_record = SeqIO.read(f, "fasta")
                    return str(sequence_record.seq)
            except Exception as e:
                logger.error(f"Error reading local protein file {filepath}: {e}")
        return None

    def _cache_protein_sequence(self, gene_symbol: str, sequence: str):
        """Saves a successfully retrieved protein sequence to the local cache."""
        try:
            filepath = os.path.join(self.LOCAL_GENE_DB_PATH, f"{gene_symbol.upper()}_PROTEIN.fasta")
            os.makedirs(self.LOCAL_GENE_DB_PATH, exist_ok=True)
            with open(filepath, "w") as f:
                f.write(f">{gene_symbol.upper()} protein sequence retrieved from NCBI\n")
                f.write(sequence)
            logger.info(f"Target cached: Saved '{gene_symbol}' PROTEIN sequence to {filepath}")
        except Exception as e:
            logger.error(f"Failed to cache PROTEIN sequence for '{gene_symbol}': {e}")
            
    @retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=2, max=10))
    def _get_network_protein_sequence(self, gene_symbol: str) -> Optional[str]:
        """(Hydra Head 2) Hardened network client to fetch a protein sequence."""
        try:
            # Find the Gene ID first
            search_handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            if not search_results["IdList"]:
                logger.warning(f"No Gene ID found for {gene_symbol}")
                return None
            gene_id = search_results["IdList"][0]

            # Link from Gene to Protein
            link_handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_id, linkname="gene_protein_refseq")
            link_results = Entrez.read(link_handle)
            link_handle.close()
            
            protein_ids = [link["Id"] for link in link_results[0]["LinkSetDb"][0]["Link"]]
            if not protein_ids:
                logger.warning(f"No RefSeq protein ID found for {gene_symbol}")
                return None
            
            # Fetch the protein FASTA record
            fetch_handle = Entrez.efetch(db="protein", id=protein_ids[0], rettype="fasta", retmode="text")
            fasta_record = fetch_handle.read()
            fetch_handle.close()

            sequence = "".join(fasta_record.split("\n")[1:])
            if not sequence: return None
            
            logger.success(f"Predator Protocol: Acquired protein sequence for {gene_symbol} (Length: {len(sequence)})")
            return sequence
        except Exception as e:
            logger.error(f"Failed to fetch protein sequence for {gene_symbol}: {e}", exc_info=True)
            return None

    # This method is preserved for legacy compatibility but may be deprecated.
    def get_protein_sequence(self, accession: str) -> Optional[str]:
        """
        Fetches the canonical nucleotide sequence (CDS) for a given gene symbol.
        """
        logger.info(f"Fetching nucleotide sequence for gene: {accession}")
        try:
            search_term = f"({accession}[Gene Name]) AND Homo sapiens[Organism] AND RefSeq[Filter]"
            search_url = f"{self.BASE_URL}esearch.fcgi"
            params = {"db": "nuccore", "term": search_term, "retmode": "json", "sort": "relevance"}
            
            time.sleep(1)
            response = requests.get(search_url, params=params)
            response.raise_for_status()
            result = response.json()
            
            id_list = result.get("esearchresult", {}).get("idlist", [])
            if not id_list:
                logger.warning(f"No RefSeq ID found for {accession}.")
                return None
            
            nuccore_id = id_list[0]
            fetch_url = f"{self.BASE_URL}efetch.fcgi"
            params = {"db": "nuccore", "id": nuccore_id, "rettype": "fasta", "retmode": "text"}
            
            time.sleep(1)
            response = requests.get(fetch_url, params=params)
            response.raise_for_status()
            
            fasta_content = response.text
            sequence_lines = fasta_content.splitlines()[1:]
            sequence = "".join(sequence_lines)

            logger.success(f"Successfully fetched sequence for {accession} (length: {len(sequence)}).")
            return sequence

        except json.JSONDecodeError as e:
            logger.error(f"Failed to fetch gene sequence for {accession}: {e}")
            logger.debug(f"RAW NCBI Response Text that failed parsing:\n---\n{response.text}\n---")
            return None
        except Exception as e:
            logger.error(f"Failed to fetch gene sequence for {accession}: {e}")
            return None

    @retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=2, max=10))
    def get_dna_sequence_from_accession(self, accession_id: str) -> Optional[str]:
        """
        Fetches a specific DNA sequence directly using its NCBI accession number (e.g., NM_001754.5).
        This is the most precise method for acquiring a known-good sequence.
        """
        try:
            logger.info(f"GROUND TRUTH PROTOCOL: Acquiring sequence for accession ID: {accession_id}")
            fetch_handle = Entrez.efetch(db="nuccore", id=accession_id, rettype="fasta", retmode="text")
            fasta_data = fetch_handle.read()
            fetch_handle.close()

            sequence = "".join(fasta_data.splitlines()[1:])
            
            if not sequence:
                logger.error(f"Fetched FASTA for {accession_id} but it was empty.")
                return None
            
            logger.success(f"GROUND TRUTH acquired: Sequence for {accession_id} (Length: {len(sequence)})")
            return sequence

        except Exception as e:
            logger.error(f"GROUND TRUTH protocol failed for '{accession_id}': {e}", exc_info=True)
            return None

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def get_cds_from_accession(self, accession_id: str) -> Optional[str]:
        """
        --- OPERATION: GENBANK PRECISION ---
        Fetches the full GenBank record for an accession ID, finds the annotated
        Coding Sequence (CDS), and returns just that sequence. This is the most
        robust method for getting the correct, protein-coding portion of an mRNA.
        """
        try:
            logger.info(f"GENBANK PRECISION: Acquiring GenBank record for {accession_id} to find CDS.")
            handle = Entrez.efetch(db="nuccore", id=accession_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            cds_feature = None
            for feature in record.features:
                if feature.type == "CDS":
                    cds_feature = feature
                    break
            
            if not cds_feature:
                logger.error(f"No CDS feature found in GenBank record for {accession_id}.")
                return None
            
            # Extract the sequence using the feature's location
            cds_sequence = str(cds_feature.extract(record.seq))
            
            logger.success(f"GENBANK PRECISION: Extracted CDS for {accession_id} (Length: {len(cds_sequence)})")
            return cds_sequence

        except Exception as e:
            logger.error(f"GENBANK PRECISION protocol failed for '{accession_id}': {e}", exc_info=True)
            return None
