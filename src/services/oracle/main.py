# FORCE REDEPLOY: v1.6 (Adding debug check for BRAF)
# FORCE REDEPLOY: v1.8 (Internalizing sequence utils)
import modal
import os
import re
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Any
import asyncio
import logging
import random # Import random for confidence simulation
import sys

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
# This is CRITICAL for the local `modal deploy` command to find the `tools` module.
# It adjusts the Python path *locally* before deployment.
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


# --- DEFINITIVE FIX: Imports must reflect the 'src' subdirectory ---
# The local `sys.path` and the container's `PYTHONPATH` (set to /root)
# both expect 'src' as the top-level package for modules under src/.
from src.tools.alphamissense_client import AlphaMissenseClient
from src.tools.esm_client import ESMClient


# --- Image Definition ---
#   **Base Image**: `nvidia/cuda:12.4.0-devel-ubuntu22.04`
# *   **Python**: `3.11`
# *   **Numpy**: `1.22.0`
# *   **Transformer Engine**: `1.13`
# This is the battle-tested image definition from the proven Zeta Oracle and Evo Service.
# It correctly handles all dependencies and build complexities for Evo2.
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11" # DOCTRINE ALIGNMENT: Use Python 3.11
    )
    .apt_install(
        ["build-essential", "cmake", "ninja-build",
            "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++", "gfortran", "llvm"]
    )
    # Create a symlink to satisfy the hardcoded path in the numpy build script.
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/llvm-ar /tools/llvm/bin/llvm-ar",
        "ln -s /usr/bin/llvm-ranlib /tools/llvm/bin/llvm-ranlib"
    )
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
    })
    .run_commands(
        "echo 'Cache invalidation comment'",
        "git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .",
        # Force install numpy 1.22.0 for compatibility, then transformer_engine
        "pip install --force-reinstall numpy==1.22.0", # DOCTRINE ALIGNMENT: Pin numpy version
        "pip uninstall -y transformer-engine transformer_engine",
        "pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation"
    )
    .pip_install(
        "fastapi",
        "uvicorn[standard]",
        "loguru",
        "pydantic",
        "scikit-learn>=1.5.0",
        "func-timeout",
        "python-dotenv",
        "pandas",
        "fair-esm",
        "pyarrow",
        "flash-attn"
    )
    .env({"PYTHONPATH": "/root"})
    .add_local_dir("src", remote_path="/root/src")
    .add_local_dir("data", remote_path="/root/data")
)

# --- App Definition ---
app = modal.App("zeta-oracle", image=evo2_image)

fastapi_app = FastAPI(
    title="Zeta Oracle V2",
    description="Provides access to the Zeta Oracle for variant scoring and analysis.",
)

# --- DOCTRINE: ALLOW CROSS-ORIGIN COMMUNICATION ---
# We must add the CORS middleware to allow our frontend application (running on localhost)
# to communicate with the Zeta Oracle service. We are allowing all origins, headers,
# and methods for maximum flexibility during development. This is a critical step
# for enabling the user interface to invoke the backend services.

origins = [
    "http://localhost:5173", # Local development frontend
    "http://localhost:5174", # Local development frontend
    "https://crispro.io",      # Production domain
]

fastapi_app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# --- Pydantic Models for API Validation ---
class VariantInput(BaseModel):
    # Example: "chr7:g.140753336A>T" for AlphaMissense
    alphamissense_variant_str: str 
    # Example: "V600E" for ESM
    mutation_hgvs: str
    # The full protein sequence for ESM
    protein_sequence: str
    # Our internal identifier
    variant_id: str

class ScoringRequest(BaseModel):
    variants: List[VariantInput]

class ScoredVariant(BaseModel):
    variant_id: str
    alphamissense_score: float
    esm_score: float
    # Placeholder for our own Evo2 zero-shot score
    evo2_score: float
    # The final, fused score
    zeta_score: float

class ScoringResponse(BaseModel):
    scored_variants: List[ScoredVariant]

class OracleRequest(BaseModel):
    protein_sequence: str
    variants: List[str] # e.g., ["M1A", "D2N"] for protein changes, or "1:949696:C:T" for genomic
    variant_type: str = "protein" # or "genomic"

class OracleResponse(BaseModel):
    variant: str
    zeta_score: float
    am_score: float
    esm_score: float
    details: str
    error: str = None


# --- Models for Inhibitor Validation (Ambush Protocol Phase III) ---
class ValidatedCandidate(BaseModel):
    inhibitor_sequence: str
    stability_score: float
    commentary: str

class InhibitorValidationRequest(BaseModel):
    target_dna_sequence: str
    candidate_inhibitors: List[str]

class InhibitorValidationResponse(BaseModel):
    ranked_candidates: List[ValidatedCandidate]

# --- Models for Interaction Judging (Triumvirate Protocol) ---
class JudgingRequest(BaseModel):
    """
    The request to the Zeta Oracle's new /judge_interaction endpoint.
    It sends the PDB string of the complex and the raw AlphaFold confidence data.
    """
    job_id: str
    pdb_string: str
    alphafold_results_json: Dict[str, Any] # This will contain pAE, pLDDT, etc.

class JudgingResponse(BaseModel):
    """
    The response from the /judge_interaction endpoint.
    Returns the final, quantitative Interface Confidence Score (ICS).
    """
    job_id: str
    interface_confidence_score: float
    error: str = None


# --- App Definition ---
# --- Volume for Caching ---
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="H100:1", volumes={mount_path: volume}, scaledown_window=600, timeout=1200)
class ZetaOracle:
    # --- Internalized Sequence Utilities ---
    AA_CODES = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }

    def _get_canonical_sequence(self, gene_symbol: str) -> str:
        """Retrieves the canonical protein sequence from a local FASTA file."""
        fasta_path = f"/root/data/reference/{gene_symbol}.fasta"
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
        
        if not all([ref_aa_one, alt_aa_one, ref_aa_one in self.AA_CODES.values(), alt_aa_one in self.AA_CODES.values()]):
            raise ValueError(f"Invalid amino acid code in HGVS string: {hgvsp}")
        
        if not position < len(sequence):
            raise ValueError(f"Position {position+1} is out of bounds for sequence of length {len(sequence)}")
        
        if sequence[position] != ref_aa_one:
            raise ValueError(f"Reference AA mismatch for {hgvsp} at position {position+1}. Expected {ref_aa_one}, found {sequence[position]}.")
        
        mutated_sequence = list(sequence)
        mutated_sequence[position] = alt_aa_one
        return "".join(mutated_sequence)

    @modal.enter()
    def load_model(self):
        """
        Load the Evo2 model and other clients into memory when the container starts.
        The tier is controlled by the ORACLE_TIER environment variable.
        """
        from dotenv import load_dotenv
        load_dotenv()
        
        # Initialize clients for the fusion engine
        self.am_client = AlphaMissenseClient()
        self.esm_client = ESMClient()
        
        self.tier = os.getenv("ORACLE_TIER", "production")
        self.model = None

        print(f"--- CrisPRO Zeta Oracle ---")
        print(f"  MODE: TIER-{self.tier.upper()}")
        print(f"---------------------------")
        
        if self.tier == "production":
            from evo2 import Evo2
            print("🚀 Loading Evo2 1B model for PRODUCTION tier...")
            # Revert to the documented, correct initialization method.
            self.model = Evo2('evo2_1b_base')
            self.tokenizer = self.model.tokenizer
            self.device = 'cuda:0'
            # Note: Evo2 models handle device placement internally, no .to() needed
            print("🎉 Zeta Oracle model loaded successfully! WE ARE OPERATIONAL. 💥")
        else:
            print(f"✅ Zeta Oracle loaded in TESTING tier. Using fast, low-cost heuristics.")
            print(f"   To enable the full model, set ORACLE_TIER=production in your environment.")


    def score_variant(self, ref_sequence: str, alt_sequence: str) -> dict:
        """
        Internal method to score the likelihood difference between reference and alternate sequences.
        This version implements the DEFINITIVE SOFTWARE FIX, combining two doctrines:
        1. The Ironclad Protocol: Sequential, mixed-precision scoring to solve OOM errors.
        2. The int32 Fix: Manual tokenization with explicit dtype to solve canUse32BitIndexMath errors.
        """
        import torch
        from loguru import logger

        logger.info("--- GHOST IN THE MACHINE - FINAL CANARY v24 ---") 
        
        if self.tier == "production":
            logger.info(f"🧬 Scoring sequences (Definitive Fix): REF={len(ref_sequence)}bp ALT={len(alt_sequence)}bp")
            
            # Ironclad Protocol: Score sequentially to manage memory
            torch.cuda.empty_cache()
            logger.info("Scoring reference sequence...")
            # This is the price of dealing with 260kbp sequences.
        
            # --- DOCTRINE: DIRECT MANUAL OVERRIDE ---
            # The high-level score_sequences method has proven untrustworthy, likely due to
            # caching or an internal override. We are bypassing it and implementing the
            # log-likelihood calculation manually to guarantee a 'sum' reduction.

            def _get_total_log_likelihood(sequence: str) -> float:
                # 1. Tokenize
                input_ids = self.tokenizer.tokenize(sequence)
                if not isinstance(input_ids, torch.Tensor):
                    input_ids = torch.tensor(input_ids, dtype=torch.long) # Ensure it's a tensor
                input_ids = input_ids.to(self.device).unsqueeze(0) # Add batch dimension

                # 2. Get Raw Logits
                with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
                    model_output = self.model(input_ids)

                # Handle potential nested tuple output from the model
                logits = model_output
                while isinstance(logits, tuple):
                    logits = logits[0]

                # 3. Convert to Log-Probabilities
                # We want the log-probability of the *next* token, so we compare
                # the logits for positions 0..n-1 with the input tokens at 1..n.
                shifted_logits = logits[..., :-1, :].contiguous()
                shifted_labels = input_ids[..., 1:].contiguous()
                
                log_probs = torch.nn.functional.log_softmax(shifted_logits, dim=-1)

                # 4. Gather the log-probs of the actual tokens in the sequence
                # This gets the log-probability of the ground-truth token at each position.
                true_log_probs = torch.gather(log_probs, 2, shifted_labels.unsqueeze(-1)).squeeze(-1)
                
                # 5. Summation - this is our guaranteed 'sum' reduction.
                total_log_likelihood = true_log_probs.sum().item()
                return total_log_likelihood

            logger.info(f"Manual Override: Calculating REF likelihood for {len(ref_sequence)}bp...")
            ref_ll = _get_total_log_likelihood(ref_sequence)
            logger.info(f"REF Likelihood: {ref_ll}")
            
            logger.info(f"Manual Override: Calculating ALT likelihood for {len(alt_sequence)}bp...")
            alt_ll = _get_total_log_likelihood(alt_sequence)
            logger.info(f"ALT Likelihood: {alt_ll}")

            delta_score = alt_ll - ref_ll
            
            # --- Confidence Calculation ---
            # This remains a heuristic, but it's based on the magnitude of the scores.
            confidence = 1 - (1 / (1 + 2 * abs(delta_score)))
            confidence = min(0.99, confidence + 0.1)

            logger.success(f"📊 Zeta Score (Δ): {delta_score}, Confidence: {confidence:.2f}")
            return {
                "zeta_score": float(delta_score),
                "confidence": float(confidence),
                "status": "success",
                "reference_likelihood": float(ref_ll),
                "alternate_likelihood": float(alt_ll)
            }
        else: # Testing tier
            print(f"🧬 Scoring sequences (Testing Tier): REF={ref_sequence[:30]}... ALT={alt_sequence[:30]}...")
            # Simple heuristic: count differences
            min_len = min(len(ref_sequence), len(alt_sequence))
            differences = sum(1 for i in range(min_len) if ref_sequence[i] != alt_sequence[i])
            len_diff = abs(len(ref_sequence) - len(alt_sequence))
            total_changes = differences + len_diff
            mock_delta = -0.1 * total_changes
            
            # Add mock confidence for testing tier
            confidence = random.uniform(0.75, 0.95)

            print(f"📊 Zeta Score (Heuristic): {mock_delta}, Confidence: {confidence:.2f}")
            return {
                "zeta_score": mock_delta,
                "confidence": confidence,
                "status": "success_testing_tier",
                "reference_likelihood": -100.0,
                "alternate_likelihood": -100.0 + mock_delta
            }

    def generate_sequence(self, prompt: str, **kwargs) -> str:
        """
        Internal method to generate a DNA sequence completion for a given prompt.
        """
        if self.tier == "production":
            print(f"🔥 Forge task received (Production Tier): Generating sequence for prompt: '{prompt}'")
            output = self.model.generate([prompt], **kwargs)
            
            if output and output.sequences and output.sequences[0]:
                completion = output.sequences[0]
                print(f"✅ Generated sequence: {completion}")
                return completion
            else:
                print("❌ Generation failed to produce output.")
                raise RuntimeError("Generation failed to produce a valid sequence.")
        else: # Testing tier
            print(f"🔥 Forge task received (Testing Tier): Generating mock sequence for prompt: '{prompt}'")
            n_tokens = kwargs.get("n_tokens", 50)
            mock_sequence = "".join(["GATTACA"[i % 7] for i in range(n_tokens)])
            print(f"✅ Generated mock sequence: {mock_sequence}")
            return mock_sequence

    def embed_sequence(self, sequence: str) -> list[float]:
        """
        Internal method to generate an embedding for a given DNA sequence.
        This now correctly uses the forward pass with embedding extraction.
        """
        import torch

        if self.tier == "production":
            print(f"🧠 Generating embedding (Production Tier): SEQ={sequence[:30]}...")

            # Tokenize the input sequence
            input_ids = torch.tensor(
                self.model.tokenizer.tokenize(sequence),
                dtype=torch.int,
            ).unsqueeze(0).to('cuda:0')

            # Layer to extract embeddings from. Intermediate layers are best.
            # From the 40B model's 50 layers, we choose a later one.
            layer_name = 'blocks.40.mlp.l3'

            # Get the embeddings via the forward pass
            _logits, embeddings = self.model(
                input_ids, 
                return_embeddings=True, 
                layer_names=[layer_name]
            )
            
            # Pool the embeddings across the sequence length and convert to a list
            embedding_tensor = embeddings[layer_name].mean(dim=1).squeeze()
            print(f"✅ Embedding generated with shape: {embedding_tensor.shape}")
            return embedding_tensor.tolist()
            
        else: # Testing tier
            print(f"🧠 Generating mock embedding (Testing Tier): SEQ={sequence[:30]}...")
            return [0.1] * 8192 # Matching the hidden_size of the 40B model

    # --- ARCHITECTURE REVERT: RESTORE THE ORIGINAL API STRUCTURE ---
    # The previous refactoring was a failure. We are reverting to the stable architecture
    # of defining the ASGI app within the class, which allows direct access to `self`.
    # The root cause of the CORS error is fixed here by adding the middleware
    # directly to this application.
    @modal.asgi_app()
    def api(self):
        from fastapi import FastAPI
        from fastapi.responses import JSONResponse
        from fastapi.middleware.cors import CORSMiddleware

        app = FastAPI(
            title="Zeta Oracle V2 API",
            description="The primary interface for all Zeta Oracle capabilities.",
        )

        origins = [
            "http://localhost:5173", # Local development frontend
            "http://localhost:5174",
            "https://crispro.io",      # Production domain
        ]

        # --- THE CORRECT FIX: Apply CORS to the *actual* running app ---
        app.add_middleware(
            CORSMiddleware,
            allow_origins=origins,
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"],
        )

        @app.post("/invoke")
        def invoke(item: dict):
            action = item.get("action")
            params = item.get("params", {})

            if not action:
                return JSONResponse(content={"status": "error", "message": "No 'action' specified."}, status_code=400)

            try:
                if action == "score":
                    ref_seq = params.get("reference_sequence", "").strip().upper()
                    alt_seq = params.get("alternate_sequence", "").strip().upper()
                    
                    if not ref_seq or not alt_seq:
                         return JSONResponse(content={"status": "error", "message": "Missing sequences for 'score' action"}, status_code=400)

                    result = self.score_variant(ref_seq, alt_seq)
                    
                    if result["zeta_score"] < -0.1:
                        interpretation = "Disruptive (Likely Pathogenic)"
                    else:
                        interpretation = "Uncertain Significance"
                    
                    response_data = {**result, "interpretation": interpretation}
                    return JSONResponse(content=response_data, status_code=200)

                elif action == "generate":
                    prompt = params.get("prompt")
                    if not prompt:
                        return JSONResponse(content={"error": "No prompt for 'generate' action"}, status_code=400)
                    gen_params = params.get("gen_params", {})
                    completion = self.generate_sequence(prompt, **gen_params)
                    return JSONResponse(content={"status": "success", "completion": completion}, status_code=200)

                elif action == "embed":
                    gene = params.get("gene_symbol")
                    hgvsp = params.get("variant")
                    
                    if not gene or not hgvsp:
                        return JSONResponse(content={"error": "Missing params for 'embed' action"}, status_code=400)
                    
                    try:
                        canonical_seq = self._get_canonical_sequence(gene)
                        mutated_seq = self._apply_hgvs_mutation(canonical_seq, hgvsp)
                    except (FileNotFoundError, ValueError) as e:
                        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=400)

                    embedding = self.embed_sequence(mutated_seq)
                    return JSONResponse(content={"status": "success", "embedding": embedding}, status_code=200)

                else:
                    return JSONResponse(content={"status": "error", "message": f"Invalid action '{action}'."}, status_code=400)

            except Exception as e:
                import traceback
                # --- DEFINITIVE FIX: Use the correct 'logging' object, not 'logger' ---
                logging.error(f"💥 ZETA ORACLE CRITICAL ERROR: {e}\n{traceback.format_exc()}")
                return JSONResponse(content={"status": "error", "message": f"Internal server error: {str(e)}"}, status_code=500)
        
        @app.get("/health")
        def health_check():
            return {"status": "healthy", "service": "zeta-oracle"}
        
        # The score_variants endpoint can be restored here if needed in the future.
        
        @app.post("/judge_interaction", response_model=JudgingResponse)
        async def judge_interaction(request: JudgingRequest):
            """
            Phase III of the Triumvirate Protocol. This is the new home for the
            Interface Confidence Score (ICS) calculation. It judges the quality
            of a predicted protein-protein interaction from AlphaFold data.
            """
            from bio.PDB import PDBParser
            from io import StringIO
            import numpy as np

            logging.info(f"--- ⚖️ JUDGEMENT INITIATED for job {request.job_id} ---")
            
            try:
                plddt = request.alphafold_results_json.get('plddt', [])
                pae = request.alphafold_results_json.get('pae', [])

                if not plddt or not pae:
                    raise ValueError("Missing pLDDT or pAE data in request.")

                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("complex", StringIO(request.pdb_string))
                
                model = structure[0]
                chain_a = model['A']
                chain_b = model['B']
                
                interface_contacts = 0
                
                # --- DOCTRINE: INTERFACE CONFIDENCE SCORE (ICS) ---
                PLDDT_THRESHOLD = 70.0
                PAE_THRESHOLD = 10.0
                DISTANCE_THRESHOLD = 10.0

                pae_matrix = np.array(pae)
                
                len_chain_a = len(list(chain_a.get_residues()))

                for res_a in chain_a:
                    if plddt[res_a.id[1] - 1] < PLDDT_THRESHOLD:
                        continue
                    
                    for res_b in chain_b:
                        if plddt[len_chain_a + res_b.id[1] - 1] < PLDDT_THRESHOLD:
                            continue
                        
                        # Check distance
                        distance = res_a['CA'] - res_b['CA']
                        if distance < DISTANCE_THRESHOLD:
                            # Check pAE
                            pae_val = pae_matrix[res_a.id[1] - 1, len_chain_a + res_b.id[1] - 1]
                            if pae_val < PAE_THRESHOLD:
                                interface_contacts += 1
                
                logging.info(f"--- ✅ JUDGEMENT COMPLETE for job {request.job_id}. ICS: {interface_contacts} ---")
                return JudgingResponse(job_id=request.job_id, interface_confidence_score=interface_contacts)

            except Exception as e:
                logging.error(f"--- ❌ JUDGEMENT FAILED for job {request.job_id}: {e} ---")
                return JudgingResponse(job_id=request.job_id, interface_confidence_score=-1.0, error=str(e))

        @app.post("/validate_inhibitors", response_model=InhibitorValidationResponse)
        async def validate_inhibitors(request: InhibitorValidationRequest):
            """
            Phase III of the Ambush Protocol. This is the doctrinally correct, dedicated
            endpoint for validating inhibitor stability.
            """
            logging.info(f"--- 💀 CONFIRMATION KILL (DEDICATED): Validating {len(request.candidate_inhibitors)} candidates... ---")
            
            validated_candidates = []
            
            for inhibitor_seq in request.candidate_inhibitors:
                bound_sequence = request.target_dna_sequence + inhibitor_seq
                
                # Note: The 'self' from the class is available here.
                score_resp = self.score_variant(
                    ref_sequence=request.target_dna_sequence,
                    alt_sequence=bound_sequence
                )

                stability_score = score_resp.get("zeta_score", -9999.0)
                
                # --- DOCTRINE REFINEMENT V3: CORRECTED THRESHOLDS ---
                # My previous attempt to apply these thresholds failed. This is the correct logic.
                if stability_score > -0.15:
                    commentary = "High Stability. Priority Candidate."
                elif stability_score > -0.35:
                    commentary = "Moderate Stability. Worthy of further investigation."
                else:
                    commentary = "Low Stability. Likely a disruptive interaction."
                
                validated_candidates.append(
                    ValidatedCandidate(
                        inhibitor_sequence=inhibitor_seq,
                        stability_score=stability_score,
                        commentary=commentary
                    )
                )

            # Rank candidates by stability score (higher is better)
            validated_candidates.sort(key=lambda x: x.stability_score, reverse=True)
            
            logging.info(f"--- ✅ CONFIRMATION KILL COMPLETE (DEDICATED) ---")
            return InhibitorValidationResponse(ranked_candidates=validated_candidates)

        return app

# The local_entrypoint for testing remains correct.
@app.local_entrypoint()
def main():
    """
    A local test function to verify the oracle is working.
    """
    print("--- ⚔️ LOCAL ORACLE TEST ⚔️ ---")
    oracle = ZetaOracle()
    
    # Test with sample sequences
    ref_sequence = "ATCGATCGATCGATCG"
    alt_sequence = "ATCGATCGATCGTTCG"
    
    print(f"🧪 Testing Zeta Oracle with sequences:")
    print(f"  REF: {ref_sequence}")
    print(f"  ALT: {alt_sequence}")
    
    try:
        # We must call the internal method for local tests as it's not a web endpoint.
        result = oracle.score_variant(ref_sequence, alt_sequence)
        print("\n--- ✅ Test Scoring Result ---")
        print(f"  Result: {result}")
        print("---------------------------")
    except Exception as e:
        print(f"\n--- ❌ TEST FAILED ---")
        print(f"  Error: {e}")
        print("---------------------") 