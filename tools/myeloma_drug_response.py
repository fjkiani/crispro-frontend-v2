import requests
import json
import time
from datetime import datetime
import re
import os

# --- Configuration ---
# Set to 'local' to use the local test function, or 'remote' to call the Modal API
# API_MODE = "local"
API_MODE = "remote"
MODAL_API_URL = "https://fjkiani--variant-analysis-evo2-evo2model-web.modal.run"

# --- Tiered Impact Scoring Configuration ---
# These thresholds define how we translate the model's output into a simple impact level.
IMPACT_CONF_HIGH = 0.7  # Confidence threshold for high impact
IMPACT_CONF_MED = 0.35 # Confidence threshold for medium impact
IMPACT_DELTA_SEVERE = -0.001 # Delta score threshold for severe impact
IMPACT_DELTA_MODERATE = -0.0001 # Delta score threshold for moderate impact

# --- Gene Sets for Myeloma Drug Response Analysis ---
# Genes in the RAS/MAPK pathway are critical drivers of proliferation and often implicated
# in resistance to targeted therapies and chemotherapy in Multiple Myeloma.
RAS_MAPK_PATHWAY_GENES = [
    "KRAS",  # Frequently mutated in Myeloma, activating mutations drive resistance
    "NRAS",  # Another key RAS isoform, similar role to KRAS
    "BRAF",  # Downstream of RAS, V600E is a known resistance driver
]

# TP53 is a critical tumor suppressor. Its inactivation is associated with poor prognosis
# and broad therapy resistance. We analyze it separately.
TP53_GENE = ["TP53"]


# --- Sample Patient Data for Myeloma Drug Response ---
# We simulate three patients to test our drug response prediction logic.
myeloma_patient_data = {
    # Patient MM-201 is designed to be 'Resistant' due to a canonical activating BRAF mutation.
    "Patient_MM-201": {
        "description": "Patient with relapsed myeloma, has a known BRAF V600E mutation.",
        "somatic_mutations": [
            {"gene": "BRAF", "hgvs_c": "c.1799T>A", "hgvs_p": "p.Val600Glu", "variant_info": "chr7:140753336 A>T", "build": "hg38"},
            {"gene": "TET2", "hgvs_c": "c.4714C>T", "hgvs_p": "p.Arg1572Ter", "variant_info": "chr4:105247953 G>A", "build": "hg38"},
        ]
    },
    # Patient MM-202 is designed to be 'Sensitive' as their mutations are not in the primary resistance pathway.
    "Patient_MM-202": {
        "description": "Newly diagnosed patient, no known mutations in key drug resistance pathways.",
        "somatic_mutations": [
            {"gene": "FAM46C", "hgvs_c": "c.319C>T", "hgvs_p": "p.Arg107Cys", "variant_info": "chr1:1193988 C>T", "build": "hg38"},
            {"gene": "IGLL5", "hgvs_c": "c.509G>A", "hgvs_p": "p.Arg170His", "variant_info": "chr22:23249079 C>T", "build": "hg38"},
        ]
    },
    # Patient MM-203 is 'Intermediate' or 'Unclear'. The TP53 mutation suggests poor prognosis and potential
    # resistance, but it's not a direct MAPK pathway activator.
    "Patient_MM-203": {
        "description": "Patient with a TP53 mutation, which is a general marker of poor outcomes.",
        "somatic_mutations": [
            {"gene": "TP53", "hgvs_c": "c.524G>A", "hgvs_p": "p.Arg175His", "variant_info": "chr17:7675088 G>A", "build": "hg38"}
        ]
    }
}


def call_evo2_variant_analysis(variant_data):
    """
    Calls the deployed Evo2 model on Modal to analyze a single genetic variant.
    This function now parses the 'variant_info' string to build the correct
    JSON payload for the Evo2Model.analyze_single_variant method.
    """
    if API_MODE == "local":
        # This is a local mock for testing purposes.
        print(f"--- MOCK API CALL for {variant_data.get('gene')} ---")
        if variant_data.get("gene") == "BRAF" and "V600E" in variant_data.get("hgvs_p", ""):
            return {"classification": "Likely pathogenic", "delta_score": -0.006, "confidence": 0.95}
        if variant_data.get("gene") == "TP53":
            return {"classification": "Likely pathogenic", "delta_score": -0.005, "confidence": 0.92}
        return {"classification": "Likely benign", "delta_score": 0.0001, "confidence": 0.99}

    elif API_MODE == "remote":
        try:
            # --- Parse the variant_info string ---
            # E.g., "chr7:140753336 A>T"
            variant_string = variant_data.get('variant_info', '')
            # The regex needs to handle formats like "chr7:140753336 A>T"
            match = re.match(r'^(chr.+?):(\d+)\s+(.+?)>(.+?)$', variant_string.strip())
            if not match:
                raise ValueError(f"Invalid variant_info format: '{variant_string}'")

            chrom, pos_str, _, alt = match.groups()
            
            payload = {
                "genome": variant_data.get("build", "hg38"),
                "chromosome": chrom,
                "variant_position": int(pos_str),
                "alternative": alt,
            }
            
            # Send as a JSON body in the POST request.
            response = requests.post(MODAL_API_URL, json=payload, timeout=90)
            response.raise_for_status()
            
            # The backend script returns 'prediction' and 'classification_confidence'.
            # We need to map these to the keys our script expects.
            api_response = response.json()
            return {
                "classification": api_response.get("prediction"),
                "confidence": api_response.get("classification_confidence"),
                "delta_score": api_response.get("delta_score"),
            }

        except requests.exceptions.RequestException as e:
            print(f"API call failed for variant {variant_data.get('variant_info')}: {e}")
            return {"error": str(e), "classification": "Error", "delta_score": 0, "confidence": 0}
        except (ValueError, IndexError) as e:
            print(f"Failed to parse variant_info '{variant_data.get('variant_info')}': {e}")
            return {"error": str(e), "classification": "Error", "delta_score": 0, "confidence": 0}


def get_variant_impact_level(evo2_result):
    """
    Translates the raw Evo2 model output into a tiered impact level (0-3).
    Level 3: High confidence, severe predicted impact.
    Level 2: High confidence, moderate impact OR Medium confidence, severe impact.
    Level 1: Medium confidence, moderate impact.
    Level 0: All others (likely benign, low confidence, etc.).
    """
    classification = evo2_result.get("classification", "")
    confidence = evo2_result.get("confidence", 0)
    delta_score = evo2_result.get("delta_score", 0)

    is_pathogenic = "pathogenic" in classification.lower()

    if not is_pathogenic:
        return 0

    conf_is_high = confidence >= IMPACT_CONF_HIGH
    conf_is_med = IMPACT_CONF_MED <= confidence < IMPACT_CONF_HIGH
    delta_is_severe = delta_score <= IMPACT_DELTA_SEVERE
    delta_is_moderate = IMPACT_DELTA_MODERATE <= delta_score < IMPACT_DELTA_SEVERE

    if conf_is_high and delta_is_severe:
        return 3
    elif (conf_is_high and delta_is_moderate) or (conf_is_med and delta_is_severe):
        return 2
    elif conf_is_med and delta_is_moderate:
        return 1
    else:
        return 0


def predict_myeloma_drug_response(patient_mutations):
    """
    Predicts a myeloma patient's likely response to MAPK-pathway-related therapies
    by analyzing their somatic mutations with the Evo2 model.
    """
    analysis_results = []
    summed_impact_ras_pathway = 0
    summed_impact_tp53 = 0
    
    print("Analyzing patient mutations...")
    for variant in patient_mutations:
        print(f"  - Analyzing variant: {variant['gene']} {variant['hgvs_p']} ({variant['variant_info']})")
        evo2_result = call_evo2_variant_analysis(variant)
        time.sleep(1) # Be respectful to the API endpoint

        impact_level = get_variant_impact_level(evo2_result)
        
        analysis_details = {
            "variant": f"{variant['gene']}_{variant['hgvs_p']}",
            "evo2_result": evo2_result,
            "calculated_impact_level": impact_level
        }
        analysis_results.append(analysis_details)
        
        if variant['gene'] in RAS_MAPK_PATHWAY_GENES:
            summed_impact_ras_pathway += impact_level
        elif variant['gene'] in TP53_GENE:
            summed_impact_tp53 += impact_level

    # --- Prediction Logic ---
    # High impact in the RAS/MAPK pathway is a strong indicator of resistance.
    # A TP53 hit also suggests resistance, but we'll prioritize the direct pathway.
    final_prediction = ""
    if summed_impact_ras_pathway >= 3:
        final_prediction = "Likely Resistant to MAPK Pathway-Involved Therapy"
    elif summed_impact_ras_pathway >= 1:
        final_prediction = "Intermediate/Potential Resistance to MAPK Pathway-Involved Therapy"
    elif summed_impact_tp53 >= 2:
        final_prediction = "Potential Resistance (due to TP53 status)"
    else:
        final_prediction = "Likely Sensitive to MAPK Pathway-Involved Therapy"
        
    return {
        "prediction": final_prediction,
        "summed_impact_ras_pathway": summed_impact_ras_pathway,
        "summed_impact_tp53": summed_impact_tp53,
        "detailed_analysis": analysis_results
    }


def main():
    """
    Main function to run the myeloma drug response prediction for all sample patients.
    """
    print("="*60)
    print("  Myeloma Digital Twin: Drug Response Prediction Test")
    print("="*60)
    
    all_results = {}

    for patient_id, data in myeloma_patient_data.items():
        print(f"\n--- Running Analysis for {patient_id} ---")
        print(f"Description: {data['description']}")
        
        prediction_output = predict_myeloma_drug_response(data['somatic_mutations'])
        all_results[patient_id] = prediction_output
        
        print(f"\n>>> Final Prediction for {patient_id}: {prediction_output['prediction']}")
        print(f"  - Summed RAS/MAPK Pathway Impact: {prediction_output['summed_impact_ras_pathway']}")
        print(f"  - Summed TP53 Impact: {prediction_output['summed_impact_tp53']}")
        print("  - Detailed Variant Analysis:")
        for detail in prediction_output['detailed_analysis']:
            print(f"    - Variant: {detail['variant']}")
            print(f"      Evo2 says: '{detail['evo2_result'].get('classification')}' (Confidence: {detail['evo2_result'].get('confidence', 0):.2%}, Delta: {detail['evo2_result'].get('delta_score', 0):.6f})")
            print(f"      Calculated Impact Level: {detail['calculated_impact_level']}")

    # --- Save results to a file ---
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"myeloma_drug_response_results_{timestamp}.json"
    with open(output_filename, 'w') as f:
        json.dump(all_results, f, indent=4)
    
    print(f"\n\nResults saved to {output_filename}")
    print("="*60)


if __name__ == "__main__":
    main() 