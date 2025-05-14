import os
import json
import subprocess
import sys
import re
import time # Import time for sleep
import random # Import random for mock scores
import importlib.util # Import for dynamic imports

# Assuming get_llm_chat_response is available from tools.llm_api
from tools.llm_api import get_llm_chat_response # Already imported in original code
# Assuming CRISPRessoResultParser is available from tools.result_parser
# from tools.result_parser import CRISPRessoResultParser # Already imported in original code

# Local imports
sys.path.append(os.path.join(os.path.dirname(__file__), "tools"))
try:
    from tools.conda_manager import CondaManager
except ImportError:
    # Create a stub if conda_manager is not available
    class CondaManager:
        def __init__(self):
            self.is_conda_installed = False
            self.crispresso_env = None

        def get_installation_guide(self):
            return "Error: conda_manager.py not found. Cannot provide installation instructions."

# --- Agent State (relevant parts) ---
chopchop_config_path = "tools/chopchop/config_local.json"
# conversation_history = [] # Global history, or pass as param

# --- Global State (Simple Example) ---
# In a more complex app, this might be a class or a context manager
last_chopchop_result = None
last_crispresso_output_dir = None

# --- Placeholder for LLM Interaction ---
# Assuming this function is already defined and works as shown in the user's code
# def ask_llm_and_get_user_response(instruction_for_llm_to_ask: str, conversation_history: list):
#    ... (existing implementation) ...

# --- CHOPCHOP Configuration Functions ---
# Assuming these functions are already defined
# def generate_config_local_json(...):
#    ... (existing implementation) ...
# def handle_chopchop_config_interaction(conversation_history):
#    ... (existing implementation) ...

# --- Placeholder for other handlers ---
# Assuming these functions are already defined
# def handle_chopchop_execution_interaction(conversation_history):
#    ... (existing implementation) ...
# def handle_crispresso_execution_interaction(conversation_history):
#    ... (existing implementation) ...

# --- NEW MOCK COMPONENTS FOR THERAPY DESIGN (MULTI-GENE SYSTEM) ---

# Mock AlphaFold 3 Result Structure for a System (Updated to include CHOPCHOP metrics)
class MockSystemAlphaFoldResult:
    """Simulates parsed output from AlphaFold 3 prediction and CHOPCHOP validation for a multi-component system."""
    def __init__(self, system_sequences: dict, design_goal: str, display_callback=print):
        self.system_sequences = system_sequences
        self.design_goal = design_goal
        self.component_metrics = {} # Dictionary to hold metrics per component (including AF and CHOPCHOP)
        self.display_callback = display_callback

        self.display_callback(f"AI Research Assistant (AlphaFoldSim): Simulating structural prediction for system related to '{design_goal}'...")

        # Simulate metrics for each component in the system
        for component_name, sequence in system_sequences.items():
            self.display_callback(f"AI Research Assistant (AlphaFoldSim):   - Processing component: {component_name}")
            component_type = "protein" # Default type
            # Infer component type based on name (simple heuristic)
            if "rna" in component_name.lower() or "guide" in component_name.lower():
                component_type = "rna"
            elif "dna" in component_name.lower() or "template" in component_name.lower():
                component_type = "dna"
            self.display_callback(f"AI Research Assistant (AlphaFoldSim):     - Inferred component type: {component_type}")

            metrics = {}
            # Simulate AlphaFold Metrics
            if component_type == "protein":
                metrics["mock_plddt_score"] = round(random.uniform(0.4, 0.95), 2)
                metrics["mock_ranking_score"] = round(random.uniform(0.5, 0.98), 2)
                if design_goal == 'correct_mutation':
                    metrics["mock_structural_similarity_score"] = round(min(1.0, len(sequence) / 100.0 + random.uniform(-0.2, 0.2)), 2)
                elif design_goal == 'disrupt_oncogene':
                    metrics["mock_structural_similarity_score"] = round(random.uniform(0.1, 0.5), 2)
            elif component_type == "rna":
                metrics["mock_mfe"] = round(random.uniform(-100, -10), 2)
                metrics["mock_ranking_score"] = round(random.uniform(0.6, 0.99), 2)
            elif component_type == "dna":
                metrics["mock_ranking_score"] = round(random.uniform(0.5, 0.95), 2)
            
            self.display_callback(f"AI Research Assistant (AlphaFoldSim):     - Generated mock AF metrics for {component_name}")
            self.component_metrics[component_name] = metrics

        self.overall_mock_system_score = round(random.uniform(0.3, 0.9), 2)
        self.display_callback(f"AI Research Assistant (AlphaFoldSim):   - Generated overall mock system score: {self.overall_mock_system_score}")
        time.sleep(0.5) # Reduced sleep


# Mock AlphaFold 3 Prediction Component (Updated for Systems)
class MockAlphaFoldPredictor:
    """Simulates calling AlphaFold 3 and returning mock results for a system."""
    def predict_system(self, system_sequences: dict, design_goal: str, display_callback=print) -> MockSystemAlphaFoldResult:
        display_callback("AI Research Assistant (AlphaFoldSim): AlphaFold 3 prediction process initiated.")
        return MockSystemAlphaFoldResult(system_sequences, design_goal, display_callback=display_callback)

# Mock Evo 2 Generation Component (Updated for Systems)
class MockEvo2Generator:
    """Simulates Evo 2 sequence generation for multi-gene systems."""
    def generate_system_candidates(self, target_mutation: dict, design_goal: str, num_systems: int = 3, display_callback=print) -> list[dict]:
        display_callback(f"AI Research Assistant (Evo2Sim): Simulating Evo 2 generation for multi-gene systems targeting {target_mutation.get('hugo_gene_symbol', 'N/A')} {target_mutation.get('protein_change','N/A')} ({design_goal})...")
        mock_systems = []
        component_types = []
        if design_goal == 'correct_mutation':
            component_types = ['repair_template_dna', 'modified_cas_protein_sequence']
        elif design_goal == 'disrupt_oncogene':
            component_types = ['guide_rna_sequence', 'cas_protein_sequence']
        elif design_goal == 'design_delivery_vector':
            component_types = ['viral_packaging_signal_dna', 'therapeutic_gene_sequence', 'promoter_sequence']
        else:
            component_types = ['component_A_sequence', 'component_B_sequence']

        display_callback(f"AI Research Assistant (Evo2Sim):   - Generating {num_systems} candidates, each with components: {', '.join(component_types)}")

        for i in range(num_systems):
            system_candidate = {}
            system_candidate["system_id"] = f"MOCK_SYSTEM_{design_goal}_Candidate_{i+1}"
            system_candidate["description"] = f"Mock system candidate {i+1} for {target_mutation.get('hugo_gene_symbol', 'N/A')} {target_mutation.get('protein_change', 'N/A')}"
            display_callback(f"AI Research Assistant (Evo2Sim):     - Generating sequences for candidate {i+1} (ID: {system_candidate['system_id']})...")

            for comp_type in component_types:
                base_seq = f"MOCK_{comp_type.upper()}_{target_mutation.get('hugo_gene_symbol', 'N/A')}_{target_mutation.get('protein_change','N/A')}_v{i+1}"
                random_suffix = ''.join(random.choices('ATCG', k=random.randint(10, 50)))
                system_candidate[comp_type] = f"{base_seq}_{random_suffix}"
                display_callback(f"AI Research Assistant (Evo2Sim):       - Generated mock sequence for {comp_type}: {system_candidate[comp_type][:20]}...")
            
            display_callback(f"AI Research Assistant (Evo2Sim):       - System candidate {i+1} details: ID {system_candidate['system_id']}")
            mock_systems.append(system_candidate)
            time.sleep(0.3) # Reduced sleep

        time.sleep(0.5) # Reduced sleep
        display_callback(f"AI Research Assistant (Evo2Sim): Finished generating {len(mock_systems)} mock system candidates.")
        return mock_systems

# Mock CHOPCHOP Validation Component
def simulate_chopchop_validation(guide_rna_sequence: str, display_callback=print) -> dict:
    display_callback(f"AI Research Assistant (ChopChopSim): Simulating CHOPCHOP validation for guide: {guide_rna_sequence[:20]}...")
    time.sleep(0.2) # Reduced sleep
    mock_ontarget_score = random.uniform(0.6, 0.95)
    mock_offtarget_count = random.randint(0, 10)
    mock_mfe = random.uniform(-30, -5)
    display_callback(f"AI Research Assistant (ChopChopSim):   - Generated mock CHOPCHOP metrics: On-target {mock_ontarget_score:.2f}, Off-targets {mock_offtarget_count}, MFE {mock_mfe:.2f}")
    return {
        "chopchop_on_target_score": round(mock_ontarget_score, 2),
        "chopchop_off_target_count": mock_offtarget_count,
        "chopchop_mfe": round(mock_mfe, 2)
    }


# Mock Cancer-Specific Scoring Functions (Updated for Systems & CHOPCHOP)
def calculate_mock_system_design_score(simulated_vep_detail: dict, mock_evaluation_results: dict, design_goal: str, display_callback=print) -> dict:
    """
    Calculates a mock confidence score for a therapeutic system design.
    Now includes detailed rationale, structured simulated_findings (strengths/weaknesses),
    and structured potential_agent_suggestions for LLM-driven conversation.
    """
    display_callback("AI Research Assistant (Scoring):   Calculating mock system design score...")

    mock_evo2_confidence = simulated_vep_detail.get("evo2_confidence", 0.0)
    mock_delta_likelihood = simulated_vep_detail.get("delta_likelihood_score", 0.0)

    alphafold_results = mock_evaluation_results.get('alphafold', MockSystemAlphaFoldResult({}, ""))
    overall_af_score = alphafold_results.overall_mock_system_score if alphafold_results else 0.0
    component_metrics_af = alphafold_results.component_metrics if alphafold_results else {}

    chopchop_metrics_data = mock_evaluation_results.get('chopchop', {}).get('guide_rna_sequence', {})
    chopchop_ontarget = chopchop_metrics_data.get('chopchop_on_target_score', 0.0)
    chopchop_offtarget_count = chopchop_metrics_data.get('chopchop_off_target_count', 5)
    chopchop_mfe = chopchop_metrics_data.get('chopchop_mfe', -10.0)

    display_callback(f"AI Research Assistant (Scoring):     Inputs for scoring:")
    display_callback(f"AI Research Assistant (Scoring):       - Mock Evo2 Confidence: {mock_evo2_confidence:.2f}")
    display_callback(f"AI Research Assistant (Scoring):         L_Explanation: Simulated confidence from Evo2 for the original variant's impact. Higher suggests more confidence in the variant being pathogenic/relevant.")
    display_callback(f"AI Research Assistant (Scoring):       - Mock Evo2 Delta Likelihood: {mock_delta_likelihood:.7f}")
    display_callback(f"AI Research Assistant (Scoring):         L_Explanation: Simulated change in sequence likelihood due to original mutation. Large negative values can indicate impactful variants.")
    display_callback(f"AI Research Assistant (Scoring):       - Mock Overall AlphaFold System Score: {overall_af_score:.2f}")
    display_callback(f"AI Research Assistant (Scoring):         L_Explanation: Simulated overall structural integrity and interaction score for the entire multi-component system from AlphaFold3. Higher is better.")
    display_callback(f"AI Research Assistant (Scoring):       - Mock CHOPCHOP On-target Score: {chopchop_ontarget:.2f}")
    display_callback(f"AI Research Assistant (Scoring):         L_Explanation: Simulated CHOPCHOP on-target efficiency score for the guide RNA. Higher (towards 1.0) is better.")
    display_callback(f"AI Research Assistant (Scoring):       - Mock CHOPCHOP Off-target Count: {chopchop_offtarget_count}")
    display_callback(f"AI Research Assistant (Scoring):         L_Explanation: Simulated CHOPCHOP count of potential off-target sites. Lower is better.")
    display_callback(f"AI Research Assistant (Scoring):       - Mock CHOPCHOP MFE: {chopchop_mfe:.2f}")
    display_callback(f"AI Research Assistant (Scoring):         L_Explanation: Simulated CHOPCHOP Minimum Free Energy for guide RNA structure. More negative indicates more stable structure.")

    confidence_score = 0.0
    supporting_evidence = [] # Remains a list of strings for raw metric logging
    detailed_rationale_parts = [f"The overall simulated confidence score for this '{design_goal}' system is derived from several mock data points:"]
    
    # New structured fields
    simulated_findings = [] # List of dictionaries
    potential_agent_suggestions = [] # List of dictionaries

    # Define thresholds (examples, can be tuned)
    thresholds = {
        "evo2_confidence_strong": 0.7, "evo2_confidence_weak": 0.3,
        "overall_af_strong": 0.75, "overall_af_weak": 0.5,
        "component_af_plddt_strong": 0.8, "component_af_plddt_weak": 0.6,
        "chopchop_ontarget_strong": 0.8, "chopchop_ontarget_weak": 0.6,
        "chopchop_offtarget_low": 1, "chopchop_offtarget_high": 5, 
        "chopchop_mfe_good": -15 
    }

    # Weights remain the same
    weight_evo2_confidence = 0.2
    weight_overall_af = 0.3
    weight_af_components = 0.2
    weight_chopchop_ontarget = 0.15
    weight_chopchop_offtarget = 0.1 
    weight_chopchop_mfe = 0.05

    # --- Populate simulated_findings and detailed_rationale_parts --- 
    display_callback(f"AI Research Assistant (Scoring):     - Calculating score contribution from Evo2 original variant...")
    score_from_evo2 = mock_evo2_confidence * weight_evo2_confidence
    supporting_evidence.append(f"Mock Original Variant Evo2 Confidence: {mock_evo2_confidence:.2f}")
    supporting_evidence.append(f"Mock Original Variant Delta Likelihood: {mock_delta_likelihood:.7f}")

    if mock_evo2_confidence >= thresholds["evo2_confidence_strong"]:
        simulated_findings.append({
            'type': "strength", 'metric': "mock_evo2_confidence", 'component': "target_variant", 'value': mock_evo2_confidence,
            'description': f"High mock Evo2 confidence ({mock_evo2_confidence:.2f}) for the target variant suggests strong rationale for intervention."
        })
        detailed_rationale_parts.append(f"A high simulated Evo2 confidence ({mock_evo2_confidence:.2f}) for the target variant supports the therapeutic relevance.")
    elif mock_evo2_confidence < thresholds["evo2_confidence_weak"]:
        simulated_findings.append({
            'type': "weakness", 'metric': "mock_evo2_confidence", 'component': "target_variant", 'value': mock_evo2_confidence,
            'description': f"Low mock Evo2 confidence ({mock_evo2_confidence:.2f}) for target variant; biological validation of target's role is critical."
        })
        detailed_rationale_parts.append(f"The lower simulated Evo2 confidence ({mock_evo2_confidence:.2f}) underscores the need for thorough biological validation of the target.")

    display_callback(f"AI Research Assistant (Scoring):     - Calculating score contribution from overall AlphaFold system assessment...")
    score_from_overall_af = overall_af_score * weight_overall_af
    supporting_evidence.append(f"Mock Overall System AF Score: {overall_af_score:.2f}")
    if overall_af_score >= thresholds["overall_af_strong"]:
        simulated_findings.append({
            'type': "strength", 'metric': "mock_overall_af_score", 'component': "overall_system", 'value': overall_af_score,
            'description': f"Strong mock overall AlphaFold system score ({overall_af_score:.2f}) indicates good simulated structural integrity."
        })
        detailed_rationale_parts.append(f"The system's overall structure received a strong mock AlphaFold score ({overall_af_score:.2f}), positive for assembly/function.")
    elif overall_af_score < thresholds["overall_af_weak"]:
        finding_desc = f"Moderate/Low mock overall AlphaFold system score ({overall_af_score:.2f}) suggests potential simulated assembly/stability issues."
        simulated_findings.append({
            'type': "weakness", 'metric': "mock_overall_af_score", 'component': "overall_system", 'value': overall_af_score,
            'description': finding_desc
        })
        detailed_rationale_parts.append(f"A moderate/low mock AlphaFold score for the overall system ({overall_af_score:.2f}) suggests areas for structural optimization.")
        potential_agent_suggestions.append({
            'type': "refinement_suggestion", 'text': "Review overall system structure", 
            'details': "The overall system AlphaFold score is moderate/low. We could explore if re-evaluating component interactions or linker designs (if applicable) using more detailed structural modeling might improve the predicted stability.",
            'linked_finding': "mock_overall_af_score"
        })

    display_callback(f"AI Research Assistant (Scoring):     - Calculating score contribution from aggregated AlphaFold component metrics...")
    total_component_af_aggregate = 0
    component_af_evidence = [] # for supporting_evidence
    num_af_components_scored = 0
    for comp_name, metrics in component_metrics_af.items():
        comp_af_aggregate_score_basis = 0; num_metrics_for_comp = 0
        if "mock_plddt_score" in metrics:
            plddt = metrics["mock_plddt_score"]
            comp_af_aggregate_score_basis += plddt; num_metrics_for_comp +=1
            component_af_evidence.append(f"Comp '{comp_name}' pLDDT: {plddt:.2f}")
            if plddt >= thresholds["component_af_plddt_strong"]:
                simulated_findings.append({
                    'type': "strength", 'metric': "mock_plddt_score", 'component': comp_name, 'value': plddt,
                    'description': f"Component '{comp_name}' shows good simulated structural confidence (pLDDT: {plddt:.2f})."
                })
            elif plddt < thresholds["component_af_plddt_weak"]:
                finding_desc = f"Component '{comp_name}' has lower simulated structural confidence (pLDDT: {plddt:.2f})."
                simulated_findings.append({
                    'type': "weakness", 'metric': "mock_plddt_score", 'component': comp_name, 'value': plddt,
                    'description': finding_desc
                })
                potential_agent_suggestions.append({
                    'type': "refinement_suggestion", 'text': f"Investigate {comp_name} stability",
                    'details': f"The pLDDT for '{comp_name}' ({plddt:.2f}) is on the lower side. If its structure is critical, we could investigate reasons (e.g., intrinsic disorder, interactions not modeled) or explore sequence variants for improved simulated stability.",
                    'linked_finding': "mock_plddt_score"
                })
        if "mock_ranking_score" in metrics: component_af_evidence.append(f"Comp '{comp_name}' Rank: {metrics['mock_ranking_score']:.2f}")
        if "mock_structural_similarity_score" in metrics: component_af_evidence.append(f"Comp '{comp_name}' Similarity: {metrics['mock_structural_similarity_score']:.2f}")
        if num_metrics_for_comp > 0: 
            total_component_af_aggregate += (comp_af_aggregate_score_basis / num_metrics_for_comp); num_af_components_scored +=1
    avg_component_af_score = 0
    if num_af_components_scored > 0: avg_component_af_score = total_component_af_aggregate / num_af_components_scored
    score_from_af_components = avg_component_af_score * weight_af_components
    supporting_evidence.extend(component_af_evidence)
    if avg_component_af_score > 0: detailed_rationale_parts.append(f"The average simulated structural quality (pLDDT) of components was {avg_component_af_score:.2f}.")

    display_callback(f"AI Research Assistant (Scoring):     - Calculating score contribution from CHOPCHOP metrics...")
    score_from_chopchop_ontarget = chopchop_ontarget * weight_chopchop_ontarget
    score_from_chopchop_offtarget = (1.0 / (chopchop_offtarget_count + 1)) * weight_chopchop_offtarget * 10 
    mfe_contribution_factor = min(abs(chopchop_mfe) / 25.0, 1.0) 
    score_from_chopchop_mfe = mfe_contribution_factor * weight_chopchop_mfe * 10
    supporting_evidence.extend([
        f"Mock CHOPCHOP On-target Score: {chopchop_ontarget:.2f}",
        f"Mock CHOPCHOP Off-target Count: {chopchop_offtarget_count}",
        f"Mock CHOPCHOP MFE: {chopchop_mfe:.2f}"
    ])

    if chopchop_ontarget >= thresholds["chopchop_ontarget_strong"]:
        simulated_findings.append({
            'type': "strength", 'metric': "chopchop_on_target_score", 'component': "guide_rna_sequence", 'value': chopchop_ontarget,
            'description': f"High mock CHOPCHOP on-target score ({chopchop_ontarget:.2f}), suggesting good potential guide efficacy."
        })
        detailed_rationale_parts.append(f"The guide RNA's high simulated on-target score ({chopchop_ontarget:.2f}) is a strong positive.")
    elif chopchop_ontarget < thresholds["chopchop_ontarget_weak"]:
        finding_desc = f"Lower mock CHOPCHOP on-target score ({chopchop_ontarget:.2f}) for the guide RNA."
        simulated_findings.append({
            'type': "weakness", 'metric': "chopchop_on_target_score", 'component': "guide_rna_sequence", 'value': chopchop_ontarget,
            'description': finding_desc
        })
        detailed_rationale_parts.append(f"A lower simulated guide on-target score ({chopchop_ontarget:.2f}) is noted.")
        potential_agent_suggestions.append({
            'type': "refinement_suggestion", 'text': "Optimize guide RNA on-target score",
            'details': "The guide RNA on-target score is somewhat low. We could explore alternative guide sequences for the same target or use prediction tools considering a wider range of Cas enzymes.",
            'linked_finding': "chopchop_on_target_score"
        })

    if chopchop_offtarget_count <= thresholds["chopchop_offtarget_low"]:
        simulated_findings.append({
            'type': "strength", 'metric': "chopchop_off_target_count", 'component': "guide_rna_sequence", 'value': chopchop_offtarget_count,
            'description': f"Low mock CHOPCHOP off-target count ({chopchop_offtarget_count}), indicating good predicted specificity."
        })
        detailed_rationale_parts.append(f"The guide RNA shows good predicted specificity (off-targets: {chopchop_offtarget_count}).")
    elif chopchop_offtarget_count > thresholds["chopchop_offtarget_high"]:
        finding_desc = f"High mock CHOPCHOP off-target count ({chopchop_offtarget_count})."
        simulated_findings.append({
            'type': "weakness", 'metric': "chopchop_off_target_count", 'component': "guide_rna_sequence", 'value': chopchop_offtarget_count,
            'description': finding_desc
        })
        detailed_rationale_parts.append(f"The higher mock off-target count ({chopchop_offtarget_count}) warrants caution.")
        potential_agent_suggestions.append({
            'type': "refinement_suggestion", 'text': "Address high off-target count",
            'details': "The predicted off-target count is high. We could use advanced *in silico* tools to refine predictions or explore high-fidelity Cas variants to mitigate this risk.",
            'linked_finding': "chopchop_off_target_count"
        })
        potential_agent_suggestions.append({
            'type': "experimental_suggestion", 'text': "Prioritize off-target validation",
            'details': "Given the high predicted off-target count, comprehensive experimental validation (e.g., GUIDE-seq, CIRCLE-seq) would be crucial if this guide is pursued.",
            'linked_finding': "chopchop_off_target_count"
        })

    if chopchop_mfe <= thresholds["chopchop_mfe_good"]:
        simulated_findings.append({
            'type': "strength", 'metric': "chopchop_mfe", 'component': "guide_rna_sequence", 'value': chopchop_mfe,
            'description': f"Favorable mock CHOPCHOP MFE ({chopchop_mfe:.2f}), suggesting stable guide structure."
        })
        detailed_rationale_parts.append(f"A favorable MFE ({chopchop_mfe:.2f}) suggests good guide stability.")
    else:
        simulated_findings.append({
            'type': "weakness", 'metric': "chopchop_mfe", 'component': "guide_rna_sequence", 'value': chopchop_mfe,
            'description': f"Less favorable mock CHOPCHOP MFE ({chopchop_mfe:.2f}), could impact activity/stability."
        })
        detailed_rationale_parts.append(f"The guide MFE ({chopchop_mfe:.2f}) is less optimal.")

    confidence_score = ( score_from_evo2 + score_from_overall_af + score_from_af_components + score_from_chopchop_ontarget + score_from_chopchop_offtarget + score_from_chopchop_mfe )
    display_callback(f"AI Research Assistant (Scoring):   - Final raw confidence score before capping: {confidence_score:.2f}")
    raw_confidence_score_before_capping = confidence_score 
    confidence_score = max(0.0, min(1.0, raw_confidence_score_before_capping))
    if confidence_score != raw_confidence_score_before_capping:
        if confidence_score == 1.0 and raw_confidence_score_before_capping > 1.0:
             display_callback(f"AI Research Assistant (Scoring):   - Confidence score was capped at 1.0 (was {raw_confidence_score_before_capping:.2f})")
        elif confidence_score == 0.0 and raw_confidence_score_before_capping < 0.0:
             display_callback(f"AI Research Assistant (Scoring):   - Confidence score was floored at 0.0 (was {raw_confidence_score_before_capping:.2f})")

    # Finalize detailed_rationale (remains same logic)
    if not detailed_rationale_parts: detailed_rationale_parts.append("The score is a composite of various simulated metrics.")
    if confidence_score >= 0.8: detailed_rationale_parts.append(f"Overall, strong promise (score: {confidence_score:.2f}).")
    elif confidence_score >= 0.5: detailed_rationale_parts.append(f"Overall, moderate promise (score: {confidence_score:.2f}), areas for attention noted.")
    else: detailed_rationale_parts.append(f"Overall, lower confidence (score: {confidence_score:.2f}), likely requiring significant refinement.")
    final_detailed_rationale = " ".join(detailed_rationale_parts)

    # --- Populate potential_agent_suggestions (General/Standard Suggestions) ---
    # In Silico (General)
    potential_agent_suggestions.append({
        'type': "refinement_suggestion", 'text': "General in silico refinement", 
        'details': "We can always continue to refine component designs and explore alternative strategies using available *in silico* modeling and prediction tools. For example, iterative design-build-test cycles, even simulated, can be powerful.",
        'linked_finding': None
    })
    # Experimental Validation (Standard)
    standard_exp_validation_suggestions = [
        {"text": "Plan in vitro characterization", "details": "Key components like guide RNAs (cleavage activity) and proteins (expression, purification) should undergo *in vitro* characterization.", 'type': "experimental_suggestion"},
        {"text": "Design cell-based assays", "details": "Relevant cell lines should be used to assess on-target editing efficiency and functional outcomes (e.g., target gene knockout, protein expression changes via Western Blot/FACS).", 'type': "experimental_suggestion"},
        {"text": "Outline off-target analysis strategy", "details": "A rigorous experimental off-target analysis for the selected guide RNA is crucial (e.g., unbiased methods like GUIDE-seq or targeted deep sequencing).", 'type': "experimental_suggestion"},
        {"text": "Assess HDR efficiency (if applicable)", "details": "If Homology Directed Repair (HDR) is the intended outcome, its efficiency and the correctness of integration patterns must be assessed experimentally.", 'type': "experimental_suggestion"}
    ]
    for suggestion in standard_exp_validation_suggestions:
        suggestion['linked_finding'] = None # Ensure linked_finding is None for general ones
        # Avoid duplicates if a more specific one was already added
        if not any(s['text'] == suggestion['text'] for s in potential_agent_suggestions if s['type'] == suggestion['type']):
            potential_agent_suggestions.append(suggestion)
    
    # Broader Therapeutic Considerations (Standard)
    standard_broader_considerations = [
        {"text": "Discuss delivery strategy", "details": "A suitable delivery strategy (e.g., viral vectors like AAV, non-viral methods like LNPs, RNP delivery) needs to be defined and developed based on the therapeutic application and target tissue.", 'type': "broader_consideration"},
        {"text": "Plan pre-clinical studies", "details": "Pre-clinical efficacy and safety studies in appropriate animal models of the disease are essential next steps before any clinical translation.", 'type': "broader_consideration"},
        {"text": "Consider immunogenicity", "details": "Potential immunogenicity of system components (e.g., Cas enzyme) should be evaluated, and mitigation strategies developed if necessary.", 'type': "broader_consideration"},
        {"text": "Evaluate long-term effects", "details": "The long-term durability of the therapeutic effect and the potential for off-target effects over time need to be considered and assessed in longer-term studies.", 'type': "broader_consideration"}
    ]
    for suggestion in standard_broader_considerations:
        suggestion['linked_finding'] = None
        if not any(s['text'] == suggestion['text'] for s in potential_agent_suggestions if s['type'] == suggestion['type']):
            potential_agent_suggestions.append(suggestion)

    return {
        "confidence_score": round(confidence_score, 2),
        "rationale": final_detailed_rationale, # Using detailed_rationale as the main rationale now
        "detailed_rationale": final_detailed_rationale, # Keep for clarity if UI specifically looks for it
        "simulated_findings": simulated_findings, # New structured list
        "potential_agent_suggestions": potential_agent_suggestions, # New structured list
        "supporting_evidence": supporting_evidence 
    }


# --- NEW HANDLER FOR DESIGN THERAPY INTERACTION (UPDATED FOR SYSTEMS & CHOPCHOP) ---

def handle_design_therapy_interaction(target_mutation_details: dict, conversation_history: list, display_callback=print):
    display_callback("AI Research Assistant: Initiating mock multi-tool integrated therapy system design workflow...")
    initial_message = "Okay, I can help you explore potential multi-gene CRISPR therapy system designs for a cancer mutation using simulated Evo 2, CHOPCHOP, and AlphaFold 3 data."
    display_callback(f"AI Research Assistant: {initial_message}")
    conversation_history.append({"role": "assistant", "content": initial_message})
    mock_target_mutation = target_mutation_details
    try:
        spec = importlib.util.spec_from_file_location("mock_evo2_api", "tools/mock_evo2_api.py")
        mock_evo2_api = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mock_evo2_api)
        display_callback(f"AI Research Assistant: Retrieving mock VEP data for target mutation: {mock_target_mutation.get('hugo_gene_symbol', 'N/A')} {mock_target_mutation.get('protein_change', 'N/A')}")
        mock_original_vep_detail = mock_evo2_api.get_variant_effect_mock(
            gene_symbol=mock_target_mutation.get("hugo_gene_symbol", "N/A"),
            variant_query=mock_target_mutation.get("protein_change", "N/A"),
            variant_type=mock_target_mutation.get("variant_type", "N/A")
        )
        display_callback(f"AI Research Assistant: Mock VEP data retrieved (Evo2 confidence: {mock_original_vep_detail.get('evo2_confidence', 'N/A'):.2f}).")
    except ImportError:
        display_callback("AI Research Assistant: Warning: mock_evo2_api.py not found. Cannot generate mock VEP data for scoring.")
        mock_original_vep_detail = {"simulated_classification": "UNKNOWN", "classification_reasoning": "Mock VEP API not available.", "evo2_confidence": 0.5, "delta_likelihood_score": 0.0}
    except Exception as e:
        display_callback(f"AI Research Assistant: Error retrieving mock VEP data: {e}")
        mock_original_vep_detail = {"simulated_classification": "ERROR", "classification_reasoning": str(e), "evo2_confidence": 0.0, "delta_likelihood_score": 0.0}

    mock_design_goal = 'disrupt_oncogene'
    if mock_target_mutation.get('hugo_gene_symbol') == 'TP53': mock_design_goal = 'correct_mutation'
    elif mock_target_mutation.get('hugo_gene_symbol') == 'KRAS': mock_design_goal = 'design_delivery_vector'
    target_message = f"Focusing on the mock mutation: {mock_target_mutation.get('hugo_gene_symbol', 'N/A')} {mock_target_mutation.get('protein_change', 'N/A')} (Goal: {mock_design_goal})."
    display_callback(f"AI Research Assistant: {target_message}")
    conversation_history.append({"role": "assistant", "content": target_message})

    evo2_generator = MockEvo2Generator()
    mock_system_candidates = evo2_generator.generate_system_candidates(
        target_mutation=mock_target_mutation,
        design_goal=mock_design_goal,
        num_systems=3, # Reduced for quicker demo
        display_callback=display_callback
    )
    display_callback(f"AI Research Assistant: Mock Evo 2 generation complete. Generated {len(mock_system_candidates)} system candidates.")

    alphafold_predictor = MockAlphaFoldPredictor()
    mock_recommendations = []
    display_callback("AI Research Assistant: Simulating evaluation (CHOPCHOP & AlphaFold 3) and scoring for system candidates...")
    for i, system_candidate in enumerate(mock_system_candidates):
        display_callback(f"AI Research Assistant: Processing system candidate {i+1}/{len(mock_system_candidates)} (ID: {system_candidate.get('system_id', 'N/A')})...")
        mock_chopchop_results = {}
        if 'guide_rna_sequence' in system_candidate:
            guide_seq = system_candidate['guide_rna_sequence']
            mock_chopchop_metrics_for_guide = simulate_chopchop_validation(guide_seq, display_callback=display_callback)
            mock_chopchop_results['guide_rna_sequence'] = mock_chopchop_metrics_for_guide
            display_callback(f"AI Research Assistant: Mock CHOPCHOP validation complete for guide in candidate {i+1}.")
        else:
            display_callback(f"AI Research Assistant: No guide RNA component in candidate {i+1} for CHOPCHOP validation.")

        mock_af_result = alphafold_predictor.predict_system(system_candidate, mock_design_goal, display_callback=display_callback)
        display_callback(f"AI Research Assistant: Mock AlphaFold 3 prediction complete for candidate {i+1}.")

        mock_evaluation_results = {'alphafold': mock_af_result, 'chopchop': mock_chopchop_results}
        mock_score_details = calculate_mock_system_design_score(
            simulated_vep_detail=mock_original_vep_detail,
            mock_evaluation_results=mock_evaluation_results,
            design_goal=mock_design_goal,
            display_callback=display_callback
        )
        display_callback(f"AI Research Assistant: Mock scoring complete for candidate {i+1} (Confidence: {mock_score_details.get('confidence_score', 'N/A'):.2f}).")
        mock_recommendation = {
            "editing_type": f"Mock System Design ({mock_design_goal})",
            "recommended_approach": system_candidate.get("description", f"System Candidate {i+1}"),
            "designed_system_components": system_candidate,
            "potential_tools": ["Mock System Delivery Tool A", "Mock System Analysis Tool B"],
            "confidence_score": mock_score_details.get("confidence_score", 0.0),
            "supporting_evidence": mock_score_details.get("supporting_evidence", []),
            "mock_alphafold_system_metrics": mock_af_result.component_metrics if mock_af_result else {},
            "mock_overall_af_score": mock_af_result.overall_mock_system_score if mock_af_result else 0.0,
            "mock_chopchop_metrics": mock_chopchop_results, # This contains the raw CHOPCHOP output for the guide
            
            # Updated fields with new structure
            "detailed_rationale": mock_score_details.get("detailed_rationale", "Rationale not available."),
            "simulated_findings": mock_score_details.get("simulated_findings", []),
            "potential_agent_suggestions": mock_score_details.get("potential_agent_suggestions", [])
        }
        mock_recommendations.append(mock_recommendation)
        time.sleep(0.1) # Reduced sleep

    mock_recommendations.sort(key=lambda x: x.get('confidence_score', 0.0), reverse=True)

    if mock_recommendations:
        summary_for_llm = f"I have generated and evaluated {len(mock_recommendations)} potential multi-gene CRISPR therapy system designs for {mock_target_mutation.get('hugo_gene_symbol', 'N/A')} {mock_target_mutation.get('protein_change', 'N/A')} (Goal: {mock_design_goal}). Here are the top recommendations based on simulated Evo 2, CHOPCHOP, and AlphaFold 3 analysis:\\n\\n"
        for i, rec in enumerate(mock_recommendations[:2]): # Show top 2 for brevity in LLM prompt
            summary_for_llm += f"Recommendation {i+1} (Confidence: {rec.get('confidence_score', 0.0):.2f}):\\n"
            summary_for_llm += f"  Approach: {rec.get('recommended_approach', 'N/A')}\\n"
            system_components = rec.get('designed_system_components', {})
            summary_for_llm += f"  Key Components (Mock): {', '.join(system_components.keys()) if system_components else 'N/A'}\\n"
            summary_for_llm += f"  Rationale (Detailed): {rec.get('detailed_rationale', 'N/A')}\\\\n"
            
            # Update LLM summary to hint at structured findings/suggestions if desired, or keep high-level
            # For now, the main change is that the full structured data is available in mock_recommendation for frontend processing
            findings_summary_parts = []
            for f_item in rec.get('simulated_findings', [])[:2]: # Max 2 findings for LLM summary brevity
                # Ensure description is truncated safely and any internal quotes are handled if necessary (though unlikely for this content)
                description_snippet = f_item['description'][:30].replace("\n", " ") + "..."
                findings_summary_parts.append(f"{f_item['type'].capitalize()}: {f_item['metric']} for {f_item['component']} (value: {f_item['value']}) - '{description_snippet}'")
            if findings_summary_parts:
                summary_for_llm += f"  Key Simulated Findings (examples): { '; '.join(findings_summary_parts) }\\\\n"

            suggestions_summary_parts = []
            for s_item in rec.get('potential_agent_suggestions', [])[:2]: # Max 2 suggestions for LLM summary
                suggestions_summary_parts.append(f"{s_item['text']} ({s_item['type'].replace('_', ' ').title()})")
            if suggestions_summary_parts:
                summary_for_llm += f"  Potential Areas for Discussion: { '; '.join(suggestions_summary_parts) }\\\\n"
            
            # The original supporting_evidence_str was removed in the previous diff, let's ensure it or an equivalent is present if needed by the LLM prompt.
            # For now, focusing on the new structured summaries.
            supporting_evidence_str = '; '.join(rec.get('supporting_evidence', [])[:3]) # Show first 3 raw evidences for brevity
            summary_for_llm += f"  Supporting Evidence Snippet (Mock): {supporting_evidence_str}\\\\n\\\\n" # Add back double newline

        presentation_instruction_for_llm = (
            f"Present the following mock multi-tool integrated CRISPR therapy system recommendations to the user. "
            f"Explain that these are based on simulated Evo 2, CHOPCHOP, and AlphaFold 3 data for {mock_target_mutation.get('hugo_gene_symbol', 'N/A')} {mock_target_mutation.get('protein_change', 'N/A')} with the goal of {mock_design_goal}. "
            f"Clearly state that the confidence scores, rationale, and supporting evidence are **simulated** for demonstration purposes and represent the evaluation of an entire multi-gene system. "
            f"Explain what the confidence score conceptually represents for a system (likelihood of the overall system achieving the therapeutic goal based on the simulated data). "
            f"Summarize the top recommendations provided, mentioning the key components of each system:\\n\\n{summary_for_llm}"
            f"Finally, ask the user if they would like to explore any of these mock system recommendations further or try a different design goal."
        )
        formatted_conversation_history = []
        for item in conversation_history:
            if isinstance(item, dict) and 'role' in item and 'content' in item:
                formatted_conversation_history.append(item)
            elif isinstance(item, str):
                formatted_conversation_history.append({'role': 'user', 'content': item})
        llm_presentation_input = formatted_conversation_history + [{'role': 'system', 'content': presentation_instruction_for_llm}]
        try:
            llm_presentation = get_llm_chat_response(llm_presentation_input, provider="gemini")
            display_callback(f"AI Research Assistant (LLM): {llm_presentation}")
            conversation_history.append({"role": "assistant", "content": llm_presentation})
        except Exception as e:
            display_callback(f"AI Research Assistant (LLM Error): Failed to get LLM presentation - {e}")
            conversation_history.append({"role": "assistant", "content": "Error: Could not retrieve summary from LLM."})
    else:
        no_recommendations_message = "AI Research Assistant: No mock therapy system designs were generated for this mutation."
        display_callback(no_recommendations_message)
        conversation_history.append({"role": "assistant", "content": no_recommendations_message})
    return mock_recommendations


# --- Example of how this new handler could be called in a main loop ---
# This part is illustrative and depends on your main agent loop structure
# Assuming you have a loop that gets user input and calls handlers based on commands

# Example main loop structure (replace with your actual loop)
# def main_agent_loop():
#    conversation_history = [] # Or load from session state
#    print("AI Research Assistant: Hello! How can I help you with CRISPR today?")
#    conversation_history.append({"role": "assistant", "content": "Hello! How can I help you with CRISPR today?"})
#
#    while True:
#        user_input = input("You: ").strip()
#        conversation_history.append({"role": "user", "content": user_input})
#
#        if user_input.lower() == "exit":
#            print("AI Research Assistant: Goodbye!")
#            break
#        elif "configure chopchop" in user_input.lower():
#            handle_chopchop_config_interaction(conversation_history)
#        elif "run chopchop" in user_input.lower() or "design guides" in user_input.lower():
#             # In a real app, you'd need to determine the target from context/user input
#             # For this mock scenario, let's trigger the new design flow
#             handle_design_therapy_interaction(
#                 target_mutation_details={"hugo_gene_symbol": "BRAF", "protein_change": "V600E", "variant_type": "Missense_Mutation"}, 
#                 conversation_history=conversation_history
#             )
#        elif "run crispresso" in user_input.lower() or "analyze editing" in user_input.lower():
#            handle_crispresso_execution_interaction(conversation_history)
#        elif "design therapy" in user_input.lower() or "design system" in user_input.lower(): # New commands to trigger the mock therapy design
#             handle_design_therapy_interaction(
#                 target_mutation_details={"hugo_gene_symbol": "KRAS", "protein_change": "G12C", "variant_type": "Missense_Mutation"},
#                 conversation_history=conversation_history
#             )
#        else:
#            response = ask_llm_and_get_user_response(f"The user said: '{user_input}'. Respond helpfully based on the conversation history and available tools (CHOPCHOP, CRISPResso2, therapy design mock). If they are asking to design a therapy or guides, suggest the relevant commands. If they are asking about analysis, suggest the relevant commands. If they are asking about configuration, suggest that. Otherwise, provide a general helpful response.", conversation_history)
#            # The ask_llm_and_user_response function already appends to history
#            # print(f"AI Research Assistant: {response}") # Printed inside the function
#
# # To run the example mock scenario directly:
# if __name__ == '__main__':
#     # This block will run the mock therapy design directly for demonstration
#     print("--- Running Mock Multi-Gene Therapy System Design Scenario ---")
#     # Example call with target_mutation_details
#     sample_target_mutation = {
#         "hugo_gene_symbol": "BRAF",
#         "protein_change": "V600E",
#         "variant_type": "Missense_Mutation",
#         "genomic_coordinate_hg38": "chr7:140753336A>T",
#         "allele_frequency": 0.45,
#         "mutation_id": "MU100001"
#     }
#     handle_design_therapy_interaction(target_mutation_details=sample_target_mutation, conversation_history=[]) # Start with empty history
#     print("--- Mock Scenario Complete ---")
#     # You could add other handler calls here to test them sequentially
#     # handle_chopchop_config_interaction([])
#     # handle_chopchop_execution_interaction([])
#     # handle_crispresso_execution_interaction([])
