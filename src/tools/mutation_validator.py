import argparse
import logging
import requests
import json

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MutationValidator:
    """
    Validates the predicted impact of a mutation by cross-referencing it
    with external databases like ClinVar and COSMIC.
    """
    CLINVAR_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    COSMIC_API_URL = "https://cancer.sanger.ac.uk/cosmic/search" # Placeholder, requires specific API

    def __init__(self, gene, mutation, zeta_score):
        """
        Initializes the validator.

        Args:
            gene (str): The HUGO symbol of the gene (e.g., 'ASXL1').
            mutation (str): The specific mutation to validate (e.g., 'p.Gly646fs').
            zeta_score (float): The internally calculated Zeta Score for the mutation.
        """
        self.gene = gene
        self.mutation = mutation
        self.zeta_score = zeta_score
        self.session = requests.Session()

    def query_clinvar(self):
        """
        Queries the ClinVar API to get pathogenicity information.
        """
        logging.info(f"Querying ClinVar for {self.gene} {self.mutation}...")
        
        search_params = {
            'db': 'clinvar',
            'term': f'"{self.gene}"[gene] AND "{self.mutation}"[title]',
            'retmode': 'json'
        }
        try:
            # Step 1: Find the record UID using esearch
            search_response = self.session.get(self.CLINVAR_API_URL + "esearch.fcgi", params=search_params)
            search_response.raise_for_status()
            search_data = search_response.json()
            
            id_list = search_data.get('esearchresult', {}).get('idlist', [])
            if not id_list:
                logging.info(f"No exact match found in ClinVar for {self.gene} {self.mutation}.")
                return "Not Found"

            uid = id_list[0]
            logging.info(f"Found ClinVar record UID: {uid}. Fetching summary...")

            # Step 2: Fetch the record summary using the UID
            summary_params = {
                'db': 'clinvar',
                'id': uid,
                'retmode': 'json'
            }
            summary_response = self.session.get(self.CLINVAR_API_URL + "esummary.fcgi", params=summary_params)
            summary_response.raise_for_status()
            summary_data = summary_response.json()

            # Step 3: Parse the clinical significance
            try:
                # Prioritize Oncogenicity classification for our use case
                classification = summary_data['result'][uid].get('oncogenicity_classification')
                if not classification:
                    # Fallback to germline classification if oncogenicity is not available
                    classification = summary_data['result'][uid].get('germline_classification')

                if classification and 'description' in classification:
                    clinical_significance = classification['description']
                    logging.info(f"Found clinical significance: {clinical_significance}")
                    return clinical_significance
                else:
                    raise KeyError("Relevant classification not found")
                    
            except (KeyError, IndexError):
                logging.warning("Could not parse clinical significance from ClinVar record. Dumping response:")
                logging.warning(json.dumps(summary_data, indent=2))
                return "Record Found, Significance Unknown"

        except requests.exceptions.RequestException as e:
            logging.error(f"Failed to query ClinVar: {e}")
            return "API Error"
        except (KeyError, IndexError) as e:
            logging.error(f"Error parsing ClinVar response: {e}")
            return "Parsing Error"

    def query_cosmic(self):
        """
        Placeholder for querying the COSMIC database.
        Note: COSMIC has strict access policies and may not have a simple public query API.
        """
        logging.warning("COSMIC query is a placeholder. A dedicated API client or manual lookup is likely required.")
        return "Not Implemented"

    def calculate_confidence_score(self, clinvar_status):
        """
        Calculates a confidence score based on agreement between the internal
        Zeta score and the external ClinVar validation.

        Args:
            clinvar_status (str): The clinical significance from ClinVar.

        Returns:
            float: A confidence score between 0.0 and 1.0.
        """
        logging.info("Calculating confidence score based on external validation...")
        
        # Normalize ClinVar status for easier matching
        clinvar_status_lower = clinvar_status.lower()
        is_pathogenic = "pathogenic" in clinvar_status_lower or "oncogenic" in clinvar_status_lower
        is_uncertain = "uncertain" in clinvar_status_lower
        is_benign = "benign" in clinvar_status_lower
        is_not_found = "not found" in clinvar_status_lower or "error" in clinvar_status_lower

        # High confidence: Strong internal score agrees with external validation
        if self.zeta_score < -10 and is_pathogenic:
            return 0.9

        # Medium confidence: Moderate agreement or conflicting signals
        if (self.zeta_score <= -5 and self.zeta_score >= -10) and is_pathogenic:
            return 0.6
        if self.zeta_score < -10 and (is_uncertain or is_benign):
            return 0.5 # Conflicting signal
        if is_not_found:
             return 0.5 # Cannot validate externally

        # Low confidence: Strong disagreement
        if self.zeta_score > -5 and is_pathogenic:
            return 0.2
            
        # Default to medium-low confidence for other cases
        return 0.4

    def run_validation(self):
        """
        Runs the full validation process.
        """
        logging.info(f"--- Running Validation for {self.gene} {self.mutation} (Zeta: {self.zeta_score}) ---")
        clinvar_status = self.query_clinvar()
        cosmic_status = self.query_cosmic()
        
        confidence = self.calculate_confidence_score(clinvar_status)
        
        print("\n--- Validation Report ---")
        print(f"  Gene: {self.gene}")
        print(f"  Mutation: {self.mutation}")
        print(f"  Zeta Score: {self.zeta_score}")
        print(f"  ClinVar Status: {clinvar_status}")
        print(f"  COSMIC Status: {cosmic_status}")
        print(f"  Confidence Score: {confidence}")
        print("------------------------")


def main():
    parser = argparse.ArgumentParser(description="Validate a mutation's predicted impact against public databases.")
    parser.add_argument("gene", help="The HUGO symbol of the gene.")
    parser.add_argument("mutation", help="The description of the mutation (e.g., 'p.Gly646fs').")
    parser.add_argument("zeta_score", type=float, help="The internal Zeta Score of the mutation.")
    args = parser.parse_args()

    validator = MutationValidator(args.gene, args.mutation, args.zeta_score)
    validator.run_validation()

if __name__ == "__main__":
    main() 