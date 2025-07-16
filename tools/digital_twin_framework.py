import pandas as pd
from typing import List, Dict, Any, Tuple
import logging
import httpx
from abc import ABC, abstractmethod
import re
import os

# --- Service URLs ---
# Corrected to use the Zeta Oracle service, not the Evo service
# The URL must include the class name as per this project's convention.
ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run"

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DigitalTwinBase(ABC):
    """
    An abstract base class for creating disease-specific Digital Twins.
    This framework provides the core logic for analyzing patient mutations by calling the Zeta Oracle.
    """
    # --- Sequence Utilities (adapted from Zeta Oracle) ---
    AA_CODES = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
        'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
        'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
        'Tyr': 'Y', 'Val': 'V'
    }

    def _get_canonical_sequence(self, gene_symbol: str) -> str:
        """Retrieves the canonical protein sequence from a local FASTA file."""
        # This assumes a specific directory structure for reference genomes.
        fasta_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'reference', f"{gene_symbol}.fasta")
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file for gene '{gene_symbol}' not found at {fasta_path}")
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
        return "".join([line.strip() for line in lines if not line.startswith('>')])

    def _apply_hgvs_mutation(self, sequence: str, hgvsp: str) -> str:
        """Applies a missense mutation described in HGVS protein notation."""
        match = re.match(r"p\.\(?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|[A-Z])\)?", hgvsp)
        if not match:
            raise ValueError(f"Invalid or unsupported HGVS notation: {hgvsp}")
        
        groups = match.groupdict()
        ref_aa, pos_str, alt_aa = groups['ref'], groups['pos'], groups['alt']
        position = int(pos_str) - 1
        
        ref_aa_one = self.AA_CODES.get(ref_aa) if len(ref_aa) == 3 else ref_aa
        alt_aa_one = self.AA_CODES.get(alt_aa) if len(alt_aa) == 3 else alt_aa
        
        if not ref_aa_one or not alt_aa_one:
             raise ValueError(f"Invalid amino acid code in HGVS string: {hgvsp}")
        
        if not position < len(sequence):
            raise ValueError(f"Position {position+1} is out of bounds for sequence of length {len(sequence)}")
        
        if sequence[position] != ref_aa_one:
            raise ValueError(f"Reference AA at position {position+1} is '{sequence[position]}', but HGVS string specifies '{ref_aa_one}'.")
        
        mutated_sequence = list(sequence)
        mutated_sequence[position] = alt_aa_one
        return "".join(mutated_sequence)

    @abstractmethod
    def get_pathways(self) -> Dict[str, List[str]]:
        raise NotImplementedError

    @abstractmethod
    def _apply_prediction_logic(self, pathway_impact_scores: Dict[str, float]) -> str:
        raise NotImplementedError

    def _get_variant_impact_from_oracle(self, mutation: Dict[str, Any]) -> Tuple[float, Dict]:
        """
        Calls the Zeta Oracle service to get a pathogenicity score for a single variant.
        """
        try:
            gene = mutation.get("gene")
            hgvsp = mutation.get("hgvs_p")

            # 1. Get the reference sequence
            ref_seq = self._get_canonical_sequence(gene)
            
            # 2. Create the alternate sequence
            alt_seq = self._apply_hgvs_mutation(ref_seq, hgvsp)

            # 3. Call the Zeta Oracle
            payload = {
                "action": "score",
                "params": {
                    "reference_sequence": ref_seq,
                    "alternate_sequence": alt_seq,
                }
            }
            url = f"{ZETA_ORACLE_URL}/invoke"
            logger.info(f"Querying Zeta Oracle for variant: {hgvsp}")

            # Increased timeout to 300s to handle Modal cold starts for the H100 GPU.
            with httpx.Client(timeout=300.0) as client:
                response = client.post(url, json=payload)
                response.raise_for_status()
                oracle_result = response.json()
            
            logger.info(f"Zeta Oracle response: {oracle_result}")
            
            # 4. Interpret the result
            interpretation = oracle_result.get("interpretation", "Unknown").lower()
            impact_level = 1.0 # Default to low impact
            if "disruptive" in interpretation:
                impact_level = 3.0
            elif "tolerated" in interpretation:
                impact_level = 0.0

            return impact_level, oracle_result

        except Exception as e:
            logger.error(f"Error scoring variant {mutation.get('hgvs_p')}: {e}")
            return -1.0, {"error": "Failed to score variant", "details": str(e)}

    def run_analysis(self, mutations: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Runs the full analysis pipeline for a list of patient mutations.
        """
        logger.info(f"Starting analysis for patient with {len(mutations)} mutations.")
        
        pathways = self.get_pathways()
        gene_to_pathway = {gene: path_name for path_name, genes in pathways.items() for gene in genes}
        
        pathway_impact_scores = {path_name: 0.0 for path_name in pathways}
        variant_analysis = []

        for mutation in mutations:
            gene = mutation.get("gene")
            if gene not in gene_to_pathway:
                continue

            impact_level, details = self._get_variant_impact_from_oracle(mutation)
            
            if impact_level < 0:
                variant_analysis.append({
                    "variant": f"{gene} {mutation.get('hgvs_p', 'N/A')}",
                    "calculated_impact_level": "Error",
                    "evo2_result": details # Keep key name for UI compatibility
                })
                continue

            pathway_name = gene_to_pathway[gene]
            pathway_impact_scores[pathway_name] += impact_level
            
            variant_analysis.append({
                "variant": f"{gene} {mutation.get('hgvs_p')}",
                "calculated_impact_level": impact_level,
                "evo2_result": details # Keep key name for UI compatibility
            })
        
        final_prediction = self._apply_prediction_logic(pathway_impact_scores)
        
        logger.info(f"Analysis complete. Final prediction: {final_prediction}")
        
        return {
            "prediction": final_prediction,
            "pathway_scores": pathway_impact_scores,
            "variant_analysis": variant_analysis,
        } 