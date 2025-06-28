import logging
from typing import List, Dict, Any
from tools.digital_twin_framework import DigitalTwinBase

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Gene Sets for Myeloma Drug Response Analysis ---
RAS_MAPK_PATHWAY_GENES = ["KRAS", "NRAS", "BRAF"]
TP53_GENE = ["TP53"]

class MyelomaDigitalTwin(DigitalTwinBase):
    """
    A disease-specific Digital Twin for Multiple Myeloma.

    This class defines the key genetic pathways and prediction logic relevant
    to predicting drug response in Multiple Myeloma.
    """
    def get_pathways(self) -> Dict[str, List[str]]:
        """
        Defines the gene pathways relevant to Multiple Myeloma.
        The keys are used for reporting and must match UI expectations.
        """
        return {
            "summed_impact_ras_pathway": RAS_MAPK_PATHWAY_GENES,
            "summed_impact_tp53": TP53_GENE,
        }

    def _apply_prediction_logic(self, pathway_impact_scores: Dict[str, float]) -> str:
        """
        Implements the specific logic to determine a drug response prediction
        based on the calculated pathway impact scores for Myeloma.
        """
        summed_impact_ras_pathway = pathway_impact_scores.get("summed_impact_ras_pathway", 0)
        summed_impact_tp53 = pathway_impact_scores.get("summed_impact_tp53", 0)

        # --- Prediction Logic ---
        if summed_impact_ras_pathway >= 3:
            return "Likely Resistant to MAPK Pathway-Involved Therapy"
        elif summed_impact_ras_pathway >= 1:
            return "Intermediate/Potential Resistance to MAPK Pathway-Involved Therapy"
        elif summed_impact_tp53 >= 2:
            return "Potential Resistance (due to TP53 status)"
        else:
            return "Likely Sensitive to MAPK Pathway-Involved Therapy"


def predict_myeloma_drug_response(patient_mutations: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Predicts a myeloma patient's likely response to MAPK-pathway-related therapies.

    This function instantiates the MyelomaDigitalTwin, runs the analysis, and
    formats the output to be compatible with the Streamlit UI.
    """
    logger.info("Initializing Myeloma Digital Twin analysis...")
    myeloma_twin = MyelomaDigitalTwin()
    results = myeloma_twin.analyze_patient_mutations(patient_mutations)
    
    # Adapt the core analysis output to the format expected by the Streamlit page
    # This keeps the UI decoupled from the core framework's data structure
    final_output = {
        "prediction": results["prediction"],
        "detailed_analysis": results["detailed_analysis"],
    }
    # Add the pathway scores to the top level of the dictionary
    final_output.update(results["pathway_impact_scores"])

    logger.info(f"Final prediction generated: {final_output['prediction']}")
    
    return final_output 