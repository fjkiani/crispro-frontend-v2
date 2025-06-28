import modal
import os
import requests
import re
import xml.etree.ElementTree as ET
from fastapi import FastAPI, BackgroundTasks
from pydantic import BaseModel
import time
import random

# --- Service Configuration ---
# This is a standard, non-GPU service that acts as our API Gateway and Orchestrator.
# It's built to be lean, fast, and highly available.
app_image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "requests", 
    "fastapi", 
    "uvicorn", 
    "pydantic", 
    "python-dotenv", 
    "biopython", 
    "pandas"
).apt_install("ncbi-blast+")

app = modal.App("command-center", image=app_image)

# --- Pydantic Models for API Validation ---
class DigitalTwinRequest(BaseModel):
    wild_type_dna: str
    mutated_dna: str

class GeneEssentialityResult(BaseModel):
    target_gene: str
    essentiality_score: float

class ThreatReportRequest(BaseModel):
    normal_essentiality_results: list[GeneEssentialityResult]
    stressed_essentiality_results: list[GeneEssentialityResult]
    environmental_factor: str
    vulnerability_threshold: float = 0.5

class GermlineTarget(BaseModel):
    gene: str
    variant: str
    reference_sequence: str
    mutated_sequence: str

class SomaticTarget(BaseModel):
    id: str
    ref_seq: str
    mut_seq: str

class FullPatientAssessmentRequest(BaseModel):
    patient_id: str
    genes_of_interest: list[str]
    germline_target: GermlineTarget
    somatic_targets: list[SomaticTarget]

class SecondHitRequest(BaseModel):
    wild_type_dna: str
    num_mutations: int = 50

# --- Command Center Logic ---
# All orchestration logic and the API itself live within this single, unified class.
@app.cls(cpu=1, memory=1024, timeout=1200)
class CommandCenter:
    @modal.enter()
    def setup(self):
        """Load environment variables and prepare connections."""
        from dotenv import load_dotenv
        load_dotenv()
        
        # We will fetch service URLs from environment variables.
        # FORCE OVERRIDE: Bypassing os.getenv to ensure the correct URL is used,
        # avoiding stale environment variables in the Modal deployment.
        self.zeta_oracle_endpoint = "https://crispro--zeta-oracle-service-zetaoracleservice-invoke.modal.run"
        self.alphafold_endpoint = os.getenv("ALPHAFOLD3_API_ENDPOINT", "https://alphafold.api.placeholder.com/predict")
        
        # Establish connection to our internal BlastService's search method with retries
        # to handle cases where the BlastService is performing its one-time DB setup.
        MAX_RETRIES = 12
        RETRY_DELAY_SECONDS = 15
        for i in range(MAX_RETRIES):
            try:
                print(f"Attempt {i+1}/{MAX_RETRIES} to connect to BlastService...")
                self.blast_search = modal.Function.lookup("blast-service", "BlastService.search")
                print("✅ Successfully connected to BlastService.")
                break
            except Exception as e:
                print(f"⚠️ Failed to connect to BlastService: {e}. Retrying in {RETRY_DELAY_SECONDS}s...")
                if i == MAX_RETRIES - 1:
                    print("❌ CRITICAL: Could not connect to BlastService after all retries. The service will fail.")
                    raise
                time.sleep(RETRY_DELAY_SECONDS)

        print("✅ Command Center is operational.")
        print(f"  - Zeta Oracle Target: {self.zeta_oracle_endpoint}")
        print(f"  - AlphaFold3 Target: {self.alphafold_endpoint}")
        print(f"  - Blast Service: Connected")

    def _query_zeta_oracle(self, action: str, params: dict) -> dict:
        """Sends a request to the Unified Zeta Oracle service."""
        print(f"🛰️  Dispatching '{action}' request to Unified Zeta Oracle...")
        payload = {"action": action, "params": params}
        response = requests.post(self.zeta_oracle_endpoint, json=payload, timeout=600) # Increased timeout for generation
        response.raise_for_status()
        print("✅ Zeta Oracle response received.")
        return response.json()

    def _query_alphafold3(self, protein_sequence: str) -> dict:
        """(Placeholder) Sends a request to the AlphaFold 3 service."""
        print(f"🧬 Dispatching request to AlphaFold3 for sequence: {protein_sequence[:30]}...")
        print("⚠️  (Placeholder) AlphaFold3 response generated.")
        return {"pdb_structure": f"PLACEHOLDER_PDB_FOR_{protein_sequence[:10]}"}

    @modal.method()
    def run_digital_twin_workflow(self, wild_type_dna: str, mutated_dna: str) -> dict:
        """The core orchestration logic for the Digital Twin workflow."""
        print("--- 💥 Executing Digital Twin Workflow 💥 ---")
        
        # Step 1: Get Functional Damage Score from Zeta Oracle
        oracle_result = self._query_zeta_oracle("score", {"reference_sequence": wild_type_dna, "alternate_sequence": mutated_dna})
        
        # Step 2: Get structural predictions from AlphaFold 3
        healthy_structure = self._query_alphafold3(wild_type_dna)
        mutated_structure = self._query_alphafold3(mutated_dna)
        
        # Step 3: Synthesize and return the final report data
        report = {
            "functional_damage_assessment": {
                "zeta_score": oracle_result.get('zeta_score', 'N/A'),
                "interpretation": oracle_result.get('interpretation', 'N/A'),
            },
            "structural_analysis": {
                "healthy_protein_structure": healthy_structure.get('pdb_structure'),
                "mutated_protein_structure": mutated_structure.get('pdb_structure'),
            }
        }
        print("--- ✅ Workflow Complete. Report generated. ---")
        return report

    @modal.method()
    def run_second_hit_simulation(self, wild_type_dna: str, num_mutations: int = 50) -> list:
        """The core orchestration logic for the Second Hit Simulation workflow."""
        print("--- 🔥 Executing Second Hit Simulation 🔥 ---")
        
        WINDOW_SIZE = 2048 # Define a reasonable window size to avoid OOM errors.

        generated_sequences = []
        for i in range(num_mutations):
            print(f"  - Generating mutation {i+1}/{num_mutations}...")
            
            # Step 1a: Select a random window from the full sequence
            if len(wild_type_dna) <= WINDOW_SIZE:
                start_pos = 0
                dna_window = wild_type_dna
            else:
                start_pos = random.randint(0, len(wild_type_dna) - WINDOW_SIZE)
                dna_window = wild_type_dna[start_pos : start_pos + WINDOW_SIZE]

            # Step 1b: Call the Zeta Oracle with just the raw DNA window.
            # The model is a foundational DNA model, not an instruction-following
            # chat model. The temperature setting will induce the necessary mutations.
            generation_params = {
                "prompt": dna_window,
                "gen_params": {
                    "n_tokens": WINDOW_SIZE, # We want a sequence of this length.
                    "temperature": 0.8
                }
            }
            oracle_result = self._query_zeta_oracle("generate", generation_params)
            
            # Step 1c: Process the result and splice it back
            completion = oracle_result.get('completion', '')
            # Find the longest contiguous block of DNA characters.
            dna_matches = re.findall(r'[ATCGN]+', completion.upper())
            if not dna_matches:
                print(f"    ⚠️ No DNA sequence found in Oracle response. Skipping.")
                generated_sequences.append(None)
                continue
            
            mutated_window = max(dna_matches, key=len)

            if len(mutated_window) == WINDOW_SIZE:
                # Reconstruct the full gene with the mutated window
                full_mutated_sequence = wild_type_dna[:start_pos] + mutated_window + wild_type_dna[start_pos + WINDOW_SIZE:]
                generated_sequences.append(full_mutated_sequence)
            else:
                print(f"    ⚠️ Oracle returned sequence of incorrect length ({len(mutated_window)} vs {WINDOW_SIZE}). Skipping.")
                generated_sequences.append(None)

        print("✅ Generation complete.")
        
        # --- The "Never Again" Protocol: Uniqueness Validation ---
        unique_sequences = set(filter(None, generated_sequences))
        if not unique_sequences or len(unique_sequences) < (num_mutations * 0.2): # Fails if less than 20% are unique
             error_msg = f"CRITICAL: Zeta Oracle failed to generate diverse mutations. Only {len(unique_sequences)} unique sequences were returned for {num_mutations} requests. Aborting."
             print(f"  - ❌ {error_msg}")
             raise ValueError(error_msg)
        print(f"  - ✅ Uniqueness validation passed. Found {len(unique_sequences)} unique mutations.")
        
        # Step 2: Score each generated mutation
        print("📊 Scoring generated mutations...")
        scored_mutations = []
        for mutant_sequence in unique_sequences:
            print(f"  - Scoring: {mutant_sequence[:30]}...")
            score_params = {"reference_sequence": wild_type_dna, "alternate_sequence": mutant_sequence}
            score_result = self._query_zeta_oracle("score", score_params)
            
            scored_mutations.append({
                "mutated_sequence": mutant_sequence,
                "zeta_score": score_result.get('zeta_score'),
                "interpretation": score_result.get('interpretation')
            })
            
        print("✅ Scoring complete.")
        
        # Step 3: Sort by impact (most disruptive first)
        scored_mutations.sort(key=lambda x: x.get('zeta_score', 0))
        
        print("--- ✅ Second Hit Simulation Complete. Report generated. ---")
        return scored_mutations

    @modal.method()
    def run_gene_essentiality_prediction(self, genomic_context: str, target_gene_sequence: str) -> dict:
        """
        [DEPRECATED - This workflow is conceptually flawed for most use cases.]
        The core orchestration logic for the Gene Essentiality Prediction workflow.
        
        This method calculates the impact of a subsequence's presence by comparing the stability
        of a larger sequence with and without the subsequence. It's only valid when 
        `target_gene_sequence` is a perfect substring of `genomic_context`. It should NOT
        be used to compare mutated vs. reference sequences.
        """
        print("--- 🔬 Executing Gene Essentiality Prediction 🔬 ---")
        
        # Step 1: Define the wild-type (with gene) and knockout (without gene) sequences
        wild_type_sequence = genomic_context
        knockout_sequence = genomic_context.replace(target_gene_sequence, "")
        
        if wild_type_sequence == knockout_sequence:
            raise ValueError("Target gene sequence not found in the provided genomic context.")
            
        print(f"  - Wild-Type Length: {len(wild_type_sequence)}")
        print(f"  - Knockout Length:  {len(knockout_sequence)}")
        
        # Step 2: Score the two sequences using the Zeta Oracle
        score_params = {"reference_sequence": wild_type_sequence, "alternate_sequence": knockout_sequence}
        score_result = self._query_zeta_oracle("score", score_params)
        
        # Step 3: Format and return the report
        essentiality_score = score_result.get('zeta_score')
        report = {
            "target_gene": target_gene_sequence[:30] + "...",
            "essentiality_score": essentiality_score,
            "interpretation": f"A score of {essentiality_score:.4f} indicates the gene's contribution to sequence stability."
        }
        
        print("--- ✅ Gene Essentiality Prediction Complete. ---")
        return report

    @modal.method()
    def run_threat_report_generation(self, normal_results: list[dict], stressed_results: list[dict], factor: str, threshold: float) -> dict:
        """
        Synthesizes gene essentiality results from two conditions (e.g., normal vs. stressed)
        to identify and report on synthetic lethalities.
        """
        print("--- 🛡️  Generating Validated Environmental Threat Report 🛡️  ---")
        
        normal_scores = {result['target_gene']: result['essentiality_score'] for result in normal_results}
        stressed_scores = {result['target_gene']: result['essentiality_score'] for result in stressed_results}

        synthetic_lethalities = []
        all_genes = set(normal_scores.keys()) | set(stressed_scores.keys())

        for gene in all_genes:
            # Strip off any '...' suffixes from the gene names for consistent matching
            clean_gene = gene.replace("...", "")
            normal_score = normal_scores.get(gene, 0)
            stressed_score = stressed_scores.get(gene, 0)
            
            # Vulnerability increase is measured by how much more essential (more negative) a gene becomes.
            # A large positive value indicates a significant increase in vulnerability.
            vulnerability_increase = normal_score - stressed_score 

            if vulnerability_increase > threshold:
                synthetic_lethalities.append({
                    "gene": clean_gene,
                    "essentiality_normal": f"{normal_score:.4f}",
                    "essentiality_stressed": f"{stressed_score:.4f}",
                    "vulnerability_increase": f"{vulnerability_increase:.4f}",
                    "comment": f"This gene's essentiality critically increases under {factor} stress, making it a prime therapeutic target."
                })
        
        # Sort by the most vulnerable first
        synthetic_lethalities.sort(key=lambda x: float(x['vulnerability_increase']), reverse=True)

        report = {
            "summary": f"Validated Environmental Threat Report for RUNX1-FPD under {factor} stress.",
            "environmental_factor": factor,
            "synthetic_lethalities_identified": len(synthetic_lethalities),
            "synthetic_lethalities": synthetic_lethalities,
        }
        
        print("--- ✅ Threat Report Generation Complete. ---")
        return report

    @modal.method()
    def run_full_patient_assessment(self, payload: FullPatientAssessmentRequest):
        """The main orchestration workflow for a full patient assessment."""
        print(f"--- 🩺 Kicking off Full Patient Assessment for Patient: {payload.patient_id} 🩺 ---")
        
        # --- 1. Germline Correction Blueprint ---
        print("\n--- Generating Germline Correction Blueprint for RUNX1 ---")
        germline_blueprint = self.run_germline_correction_blueprint.remote(
            # Pass the mutated sequence to be corrected
            corrected_sequence=payload.germline_target.mutated_sequence,
            # Pass the original reference sequence for context
            target_site_context=payload.germline_target.reference_sequence
        )

        # --- 2. Digital Twin Simulation ---
        print("\n--- Running Digital Twin Simulation on Germline Variant ---")
        digital_twin_report = self.run_digital_twin_workflow.remote(
            wild_type_dna=payload.germline_target.reference_sequence,
            mutated_dna=payload.germline_target.mutated_sequence
        )

        # --- 3. Somatic Clone Assassination ---
        assassination_plans = []
        for target in payload.somatic_targets:
            print(f"\n--- Initiating Assassination for Somatic Target: {target.id} ---")
            assassination_plan = self.run_clone_assassination_workflow.remote(
                genomic_context=target.ref_seq, # The wild-type context
                target_gene_sequence=target.mut_seq, # The mutated sequence to target
                patient_id=payload.patient_id
            )
            assassination_plans.append(assassination_plan)
            
        # --- 4. Synthesize Final Report ---
        # Note: .remote() calls within the same class are blocking and return the result directly.
        final_report = {
            "patient_id": payload.patient_id,
            "germline_correction_plan": germline_blueprint,
            "digital_twin_analysis": digital_twin_report,
            "somatic_assassination_plans": assassination_plans
        }
        
        print("\n--- ✅ Full Patient Assessment Workflow Complete. ---")
        return final_report

    @modal.method()
    def run_germline_correction_blueprint(self, corrected_sequence: str, target_site_context: str, arm_length: int = 500) -> dict:
        """
        Designs a full gene therapy blueprint for correcting a germline mutation.
        This includes designing guide RNAs and a homology-directed repair (HDR) template.
        """
        print("--- 🧬 Designing Germline Correction Blueprint 🧬 ---")
        
        # Step 1: Design guide RNAs using a helper method
        guides = self._design_guide_rnas(corrected_sequence)
        
        # Step 2: Design the HDR template
        # Find the center of the corrected sequence to design the homology arms
        center_pos = len(target_site_context) // 2
        
        # Define the start and end of the homology arms from the original reference sequence
        left_arm_start = max(0, center_pos - arm_length)
        left_arm_end = center_pos
        right_arm_start = center_pos
        right_arm_end = min(len(target_site_context), center_pos + arm_length)

        left_homology_arm = target_site_context[left_arm_start:left_arm_end]
        right_homology_arm = target_site_context[right_arm_start:right_arm_end]

        # The therapeutic payload is the corrected version of the gene
        therapeutic_payload = corrected_sequence
        
        hdr_template = f"{left_homology_arm}{therapeutic_payload}{right_homology_arm}"

        # Step 3: Off-target analysis for the best guide RNA
        best_guide_sequence = guides['best_guide']['sequence'] if guides['guides_found'] > 0 else None
        off_target_report = {}
        if best_guide_sequence:
            print(f"🔬 Performing off-target analysis for guide: {best_guide_sequence}")
            off_target_report = self.run_off_target_analysis.remote(guide_sequence=best_guide_sequence)
        
        blueprint = {
            "strategy": "Homology-Directed Repair (HDR) for Germline Correction",
            "target_gene": "RUNX1",
            "guide_rna_design": guides,
            "hdr_template_design": {
                "left_homology_arm": left_homology_arm[:50] + "...",
                "therapeutic_payload": therapeutic_payload[:60] + "...",
                "right_homology_arm": right_homology_arm[:50] + "...",
                "total_length": len(hdr_template),
            },
            "off_target_analysis": off_target_report if best_guide_sequence else "No guide RNA to analyze."
        }
        
        print("--- ✅ Germline Correction Blueprint Generated. ---")
        return blueprint

    @modal.method()
    def run_clone_assassination_workflow(self, genomic_context: str, target_gene_sequence: str, patient_id: str) -> dict:
        """
        Analyzes a somatic mutation and determines the best strategy to eliminate the clone.
        """
        print("--- ⚔️  Executing Clone Assassination Workflow ⚔️  ---")

        # Step 1: Score the functional damage of the somatic mutation.
        # This is the correct way to assess the mutation's impact. We compare the reference
        # sequence (genomic_context) with the sequence containing the mutation.
        print(f"🔬 Assessing functional damage of the somatic variant...")
        damage_assessment = self._query_zeta_oracle("score", {
            "reference_sequence": genomic_context,
            "alternate_sequence": target_gene_sequence
        })
        zeta_score = damage_assessment.get('zeta_score', 0)
        
        # Step 2: Based on the damage score, decide on a strategy.
        # This is a simplified decision tree. A real system would be more complex.
        strategy_report = {}
        if zeta_score < -5.0: # High damage score -> likely loss-of-function
             strategy_report = {
                 "chosen_strategy": "Rescue",
                 "details": "The mutation causes significant functional damage. A rescue strategy (e.g., re-introducing a functional copy via AAV) may be viable.",
                 "next_step": "Design AAV vector with wild-type cDNA.",
                 "zeta_score": f"{zeta_score:.4f}"
             }
        elif zeta_score > 5.0: # High positive score -> potential gain-of-function
            print("🔬 Mutation appears to be gain-of-function. Designing inhibitor.")
            # This is a placeholder for a more complex workflow.
            inhibitor_design = self.design_nanobody_inhibitor.remote(target_protein_name="Somatic Mutant", target_domain="Active Site")
            strategy_report = {
                "chosen_strategy": "Inhibition",
                "details": "The mutation may confer a gain-of-function. Designing a targeted inhibitor.",
                "next_step": "Synthesize and validate nanobody inhibitor.",
                "design_details": inhibitor_design,
                "zeta_score": f"{zeta_score:.4f}"
            }
        else: # Neutral or minor effect
             strategy_report = {
                 "chosen_strategy": "Monitor",
                 "details": "The mutation does not appear to cause a major functional change. Active monitoring is recommended.",
                 "next_step": "Continue periodic sequencing.",
                 "zeta_score": f"{zeta_score:.4f}"
             }

        print("--- ✅ Clone Assassination Workflow Complete. ---")
        return strategy_report

    @modal.method()
    def design_nanobody_inhibitor(self, target_protein_name: str, target_domain: str) -> dict:
        """
        Designs a novel nanobody sequence to inhibit a specific protein domain.
        """
        print("--- 🔬 Forging Novel Biologic Inhibitor 🔬 ---")
        print(f"  - Target Protein: {target_protein_name}")
        print(f"  - Target Domain: {target_domain}")

        # Step 1: Create a structured prompt for the generative model
        prompt = f"[TARGET_PROTEIN:{target_protein_name}] [TARGET_DOMAIN:{target_domain}] [THERAPEUTIC_ACTION:INHIBIT] [MOLECULE_TYPE:NANOBODY] Please generate a complete protein coding sequence for a nanobody that achieves this objective."
        
        # Step 2: Query the Zeta Oracle's generation endpoint
        generation_params = {"prompt": prompt, "gen_params": {"n_tokens": 700}} # Nanobodies are ~110-150 AA, so 330-450 nucleotides. 700 is a safe upper bound.
        oracle_result = self._query_zeta_oracle("generate", generation_params)
        
        nanobody_sequence = oracle_result.get("completion", "")

        if not nanobody_sequence:
            raise RuntimeError("Failed to generate nanobody sequence.")

        # Step 3: Assemble and return the report
        report = {
            "task": "Novel Biologic Inhibitor Design",
            "design_parameters": {
                "target_protein": target_protein_name,
                "target_domain": target_domain,
                "action": "INHIBIT",
                "molecule_type": "Nanobody"
            },
            "generated_sequence": {
                "dna_sequence": nanobody_sequence,
                "comment": "Sequence requires wet-lab validation for expression, binding affinity, and efficacy."
            }
        }
        print("--- ✅ Novel Biologic Inhibitor Forged. ---")
        return report

    @modal.method()
    def run_drug_repurposing_screen(self, mutant_protein_structure: str, drug_library: str = "DrugBank_approved") -> dict:
        """(Placeholder) Simulates a high-throughput virtual screening."""
        print("--- 💊 Executing In Silico Drug Repurposing Screen 💊 ---")
        print(f"  - Screening against library: {drug_library}")
        print(f"  - Placeholder: Returning top 3 dummy hits for structure: {mutant_protein_structure[:30]}...")
        
        report = {
            "task": "Drug Repurposing Screen",
            "library": drug_library,
            "top_hits": [
                {"drug_id": "DB00619", "drug_name": "Imatinib", "binding_affinity": -9.8},
                {"drug_id": "DB01254", "drug_name": "Dasatinib", "binding_affinity": -9.5},
                {"drug_id": "DB08901", "drug_name": "Ponatinib", "binding_affinity": -9.2},
            ]
        }
        print("--- ✅ Drug Repurposing Screen Complete. ---")
        return report

    @modal.method()
    def run_off_target_analysis(self, guide_sequence: str, num_mismatches_allowed: int = 3) -> dict:
        """
        Performs a comprehensive off-target analysis by querying the dedicated BlastService.
        """
        print("--- 🎯 Executing Off-Target Analysis 🎯 ---")
        print(f"  - Querying BlastService for gRNA: {guide_sequence}")

        # Step 1: Call the blast service
        # NOTE: For reasons specific to this Modal environment, this remote call appears
        # to be blocking and returns a dict directly, not a future.
        blast_result = self.blast_search.remote(
            query_sequence=guide_sequence, 
            num_mismatches=num_mismatches_allowed
        )

        if not blast_result or "raw_blast_xml" not in blast_result:
            print("❌ Error: Received invalid response from BlastService.")
            return {"error": "Invalid response from BlastService"}

        raw_xml = blast_result["raw_blast_xml"]

        # Step 2: Parse the results
        parsed_hits = self._parse_blast_results(raw_xml, guide_sequence)
        
        # Step 3: Assemble the report
        report = {
            "guide_sequence": guide_sequence,
            "total_off_targets_found": len(parsed_hits),
            "off_target_hits": parsed_hits
        }
        
        print(f"--- ✅ Off-Target Analysis Complete. Found {len(parsed_hits)} potential off-targets. ---")
        return report

    def _parse_blast_results(self, blast_xml: str, query_seq: str) -> list:
        """Parses BLAST XML output to identify and score off-target hits."""
        if not blast_xml:
            return []
            
        root = ET.fromstring(blast_xml)
        hits = []
        
        # Namespace is often present in NCBI BLAST XML
        ns = {'blast': 'http://www.ncbi.nlm.nih.gov'}
        
        for hit in root.findall('.//blast:Hit', ns):
            hit_def = hit.find('blast:Hit_def', ns).text
            for hsp in hit.findall('.//blast:Hsp', ns):
                hsp_qseq = hsp.find('blast:Hsp_qseq', ns).text
                hsp_hseq = hsp.find('blast:Hsp_hseq', ns).text
                
                mismatches = sum(1 for q, h in zip(hsp_qseq, hsp_hseq) if q != h and q != '-' and h != '-')
                
                # We only care about hits that are not the query sequence itself
                if hsp_hseq.upper() != query_seq.upper():
                    hits.append({
                        "hit_definition": hit_def,
                        "match_sequence": hsp_hseq,
                        "mismatches": mismatches,
                        "score": float(hsp.find('blast:Hsp_score', ns).text),
                        "e_value": float(hsp.find('blast:Hsp_evalue', ns).text),
                    })
        return hits

    def _design_guide_rnas(self, target_site_context: str) -> dict:
        """
        Internal helper to design and generate guide RNAs for a given DNA sequence.
        """
        print(f"🎯 Generating optimized guide RNAs for target site: {target_site_context[:30]}...")

        MIN_CONTEXT_LEN = 23  # 20 for gRNA + 3 for PAM
        if len(target_site_context) < MIN_CONTEXT_LEN:
            print(f"⚠️  Target site context is too short ({len(target_site_context)}bp) to find 20bp gRNAs. Minimum required: {MIN_CONTEXT_LEN}bp. Skipping gRNA generation.")
            candidate_gRNAs = []
            gRNA_comment = f"Target site context was too short to generate guide RNAs. A minimum of {MIN_CONTEXT_LEN}bp is recommended."
        else:
            gRNA_prompt = f"Scan the following DNA sequence and identify up to 5 optimal 20-nucleotide protospacer sequences for CRISPR targeting. A valid protospacer must be immediately followed by a 5'-NGG-3' PAM sequence. Return ONLY the 20-nucleotide DNA sequences, separated by spaces. Sequence: {target_site_context}"
            gRNA_result = self._query_zeta_oracle("generate", {"prompt": gRNA_prompt, "gen_params": {"n_tokens": 150}})  # (5 guides * 20 chars) + 4 spaces = 104. 150 is safe.

            # Use regex to find all 20-character DNA sequences to make parsing more robust
            completion_text = gRNA_result.get("completion", "")
            candidate_gRNAs = re.findall(r'[ATGC]{20}', completion_text.upper())

            print(f"✅ Found {len(candidate_gRNAs)} candidate guide RNAs.")
            gRNA_comment = "Candidates require off-target validation."

        report = {
            "guides_found": len(candidate_gRNAs),
            "best_guide": None,
            "candidates": candidate_gRNAs,
            "comment": gRNA_comment
        }

        if candidate_gRNAs:
            # For simplicity, we'll select the first candidate as the "best".
            # A real system would score them based on on-target efficacy and off-target risk.
            report["best_guide"] = {
                "sequence": candidate_gRNAs[0],
                "score": "N/A" 
            }
            report["comment"] += " The 'best_guide' is the first candidate found."

        return report

    @modal.asgi_app()
    def api(self):
        """This method serves the FastAPI application, consolidating all endpoints into one class."""
        fastapi_app = FastAPI(
            title="CrisPRO Command Center API",
            description="The central nervous system for the RUNX1-FPD Conquest.",
            version="1.0.0"
        )

        @fastapi_app.post("/workflow/digital_twin", summary="Create Digital Twin")
        async def create_digital_twin(request: DigitalTwinRequest):
            """Triggers the Digital Twin workflow."""
            # Note: Using .remote() on a method within the same class is blocking.
            report = self.run_digital_twin_workflow.remote(
                wild_type_dna=request.wild_type_dna, 
                mutated_dna=request.mutated_dna
            )
            return report

        @fastapi_app.post("/workflow/second_hit_simulation", summary="Simulate Evolutionary Trajectories")
        async def second_hit_simulation(request: SecondHitRequest):
            """Triggers the Second Hit Simulation workflow."""
            report = self.run_second_hit_simulation.remote(
                wild_type_dna=request.wild_type_dna,
                num_mutations=request.num_mutations
            )
            return report

        @fastapi_app.post("/workflow/predict_gene_essentiality")
        async def predict_gene_essentiality(request: dict):
            # This endpoint is intentionally left as-is due to the deprecation of its target method.
            raise NotImplementedError("This workflow has been deprecated due to conceptual flaws.")

        @fastapi_app.post("/workflow/clone_assassination")
        async def clone_assassination(request: dict):
            """
            Initiates the Clone Assassination workflow.
            """
            genomic_context = request.get("genomic_context")
            target_gene = request.get("target_gene")

            if not genomic_context or not target_gene:
                from fastapi.responses import JSONResponse
                return JSONResponse(content={"error": "genomic_context and target_gene are required"}, status_code=400)
            
            report = self.run_clone_assassination_workflow.remote(genomic_context, target_gene)
            return report

        @fastapi_app.post("/workflow/design_nanobody")
        async def design_nanobody(request: dict):
            """
            Initiates the Novel Biologic Inhibitor (Nanobody) design workflow.
            """
            target_protein_name = request.get("target_protein_name")
            target_domain = request.get("target_domain")

            if not target_protein_name or not target_domain:
                from fastapi.responses import JSONResponse
                return JSONResponse(content={"error": "target_protein_name and target_domain are required"}, status_code=400)

            report = self.design_nanobody_inhibitor.remote(target_protein_name, target_domain)
            return report

        @fastapi_app.post("/workflow/run_drug_repurposing")
        async def run_drug_repurposing(request: dict):
            """
            Initiates the In Silico Drug Repurposing screen.
            """
            mutant_protein_structure = request.get("mutant_protein_structure")

            if not mutant_protein_structure:
                from fastapi.responses import JSONResponse
                return JSONResponse(content={"error": "mutant_protein_structure is required"}, status_code=400)

            report = self.run_drug_repurposing_screen.remote(mutant_protein_structure)
            return report

        @fastapi_app.post("/workflow/off_target_analysis")
        async def off_target_analysis(request: dict):
            """
            Takes a guide RNA sequence and performs off-target analysis using the Blast service.
            """
            report = self.run_off_target_analysis.remote(
                guide_sequence=request.get("guide_sequence"),
                num_mismatches_allowed=request.get("num_mismatches_allowed", 3)
            )
            return report

        @fastapi_app.post("/workflow/generate_threat_report", summary="Generate Threat Report from Essentiality Data")
        async def generate_threat_report(request: ThreatReportRequest):
            """
            Identifies synthetic lethalities by comparing gene essentiality scores
            between a baseline (normal) and a stressed condition.
            """
            report = self.run_threat_report_generation.remote(
                normal_results=[r.dict() for r in request.normal_essentiality_results],
                stressed_results=[r.dict() for r in request.stressed_essentiality_results],
                factor=request.environmental_factor,
                threshold=request.vulnerability_threshold,
            )
            return report

        @fastapi_app.post("/workflow/full_patient_assessment", summary="Run Full Patient Assessment")
        async def full_patient_assessment(request: FullPatientAssessmentRequest):
            """
            Orchestrates a full therapeutic design workflow based on patient genomic data.
            """
            report = self.run_full_patient_assessment.remote(payload=request)
            return report

        @fastapi_app.post("/workflow/design_germline_correction")
        async def design_germline_correction(request: dict):
            """
            Designs an HDR template and guide RNAs for correcting a germline mutation.
            """
            corrected_sequence = request.get("corrected_sequence")
            target_site_context = request.get("target_site_context")
            arm_length = request.get("arm_length", 500)
            
            if not corrected_sequence or not target_site_context:
                from fastapi.responses import JSONResponse
                return JSONResponse(content={"error": "corrected_sequence and target_site_context are required"}, status_code=400)
            
            report = self.run_germline_correction_blueprint.remote(corrected_sequence, target_site_context, arm_length)
            return report

        @fastapi_app.post("/health", summary="Health Check")
        def health_check():
            """Confirm that the Command Center is alive and operational."""
            return {"status": "healthy", "service": "command-center"}
        
        return fastapi_app 