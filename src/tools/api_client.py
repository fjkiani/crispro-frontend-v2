import os
import requests
import logging

logging.basicConfig(level=logging.INFO)

class APIClient:
    """
    A dedicated client to interact with the Evo2 backend services.
    """
    def __init__(self):
        self.discriminative_endpoint = os.environ.get("EVO2_DISCRIMINATIVE_ENDPOINT")
        self.generative_endpoint = os.environ.get("EVO2_GENERATIVE_ENDPOINT")
        self.adjudicator_endpoint = os.environ.get("ADJUDICATOR_ENDPOINT")
        self.headers = {"Content-Type": "application/json"}

        if not self.discriminative_endpoint or not self.generative_endpoint:
            logging.warning("Evo2 endpoint URLs are not set in environment variables.")
            # We can set default URLs for development if needed
            # self.discriminative_endpoint = "http://localhost:8000"
            # self.generative_endpoint = "http://localhost:8001"

    def _make_request(self, endpoint, payload):
        """Helper function to make POST requests to the API."""
        if not endpoint:
            logging.error("Endpoint URL is not configured.")
            return None
        
        try:
            response = requests.post(endpoint, json=payload, headers=self.headers)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logging.error(f"API request failed for endpoint {endpoint}: {e}")
            return None

    # --- Phase 1: Discriminative Endpoints ---

    def predict_variant_impact(self, variants):
        """
        Analyzes variants to predict their impact.
        
        Args:
            variants (list): A list of variants to analyze.
        
        Returns:
            dict: The API response.
        """
        payload = {"variants": variants}
        return self._make_request(f"{self.discriminative_endpoint}/predict_variant_impact", payload)

    def predict_gene_essentiality(self, genes, context):
        """
        Predicts the essentiality of genes in a given context.
        
        Args:
            genes (list): A list of gene identifiers.
            context (str): The biological context (e.g., cell type, tissue).
        
        Returns:
            dict: The API response.
        """
        payload = {"genes": genes, "context": context}
        return self._make_request(f"{self.discriminative_endpoint}/predict_gene_essentiality", payload)

    def predict_protein_functionality_change(self, variants):
        """
        Predicts changes in protein functionality due to variants.
        
        Args:
            variants (list): A list of variants affecting protein-coding regions.
        
        Returns:
            dict: The API response.
        """
        payload = {"variants": variants}
        return self._make_request(f"{self.discriminative_endpoint}/predict_protein_functionality_change", payload)

    def predict_chromatin_accessibility(self, genomic_regions):
        """
        Predicts chromatin accessibility for given genomic regions.
        
        Args:
            genomic_regions (list): A list of genomic regions (e.g., ["chr1:123-456"]).
        
        Returns:
            dict: The API response.
        """
        payload = {"genomic_regions": genomic_regions}
        return self._make_request(f"{self.discriminative_endpoint}/predict_chromatin_accessibility", payload)

    def classify_embedding(self, embedding: list[float]):
        """
        Sends an embedding to the Adjudicator service for classification.
        
        Args:
            embedding (list): The 8192-dimensional embedding vector.
        
        Returns:
            dict: The API response from the Adjudicator.
        """
        # The Adjudicator service expects the embedding directly as the body.
        return self._make_request(f"{self.adjudicator_endpoint}/classify", embedding)

    # --- Placeholder for Phase 2: Generative Endpoints ---
    
    def generate_optimized_guide_rna(self, target_sequence):
        """
        Generates optimized guide RNAs for a target sequence.
        (Placeholder for Phase 2)
        """
        logging.info("Placeholder for generate_optimized_guide_rna")
        return None

# Example usage:
if __name__ == "__main__":
    client = APIClient()

    # Example for predict_variant_impact
    mock_variants = [{"id": "rs123", "chr": "1", "pos": 1000, "ref": "A", "alt": "G"}]
    variant_impact = client.predict_variant_impact(mock_variants)
    if variant_impact:
        logging.info(f"Variant Impact Prediction: {variant_impact}")

    # Example for predict_gene_essentiality
    mock_genes = ["TP53", "KRAS"]
    gene_essentiality = client.predict_gene_essentiality(mock_genes, context="Lung Adenocarcinoma")
    if gene_essentiality:
        logging.info(f"Gene Essentiality Prediction: {gene_essentiality}") 