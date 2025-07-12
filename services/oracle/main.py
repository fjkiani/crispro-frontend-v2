# FORCE REDEPLOY: v1.6 (Adding debug check for BRAF)
# FORCE REDEPLOY: v1.8 (Internalizing sequence utils)
import modal
import os
import re

# --- Image Definition ---
# This is the battle-tested image definition from the proven Zeta Oracle and Evo Service.
# It correctly handles all dependencies and build complexities for Evo2.
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11"
    )
    .apt_install(
        ["build-essential", "cmake", "ninja-build",
            "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"]
    )
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
    })
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install(
        "fastapi",
        "uvicorn[standard]",
        "loguru",
        "pydantic",
        "numpy==1.22.0",
        "scikit-learn==1.3.2",
        "func-timeout",
        "python-dotenv",
        "pandas"
    )
)

# --- App Definition ---
app = modal.App("zeta-oracle-v2", image=evo2_image)

# --- Volume for Caching ---
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="H100:2", volumes={mount_path: volume}, scaledown_window=600, timeout=1200)
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
        Load the Evo2 model into memory when the container starts.
        The tier is controlled by the ORACLE_TIER environment variable.
        """
        from dotenv import load_dotenv
        load_dotenv()
        
        self.tier = os.getenv("ORACLE_TIER", "production")
        self.model = None

        print(f"--- CrisPRO Zeta Oracle ---")
        print(f"  MODE: TIER-{self.tier.upper()}")
        print(f"---------------------------")
        
        if self.tier == "production":
            from evo2 import Evo2
            print("🚀 Loading Evo2 40B model for PRODUCTION tier...")
            # Use the 40B base model for maximum accuracy.
            self.model = Evo2('evo2_40b')
            print("🎉 Zeta Oracle model loaded successfully! WE ARE OPERATIONAL. 💥")
        else:
            print(f"✅ Zeta Oracle loaded in TESTING tier. Using fast, low-cost heuristics.")
            print(f"   To enable the full model, set ORACLE_TIER=production in your environment.")


    def score_variant(self, ref_sequence: str, alt_sequence: str) -> dict:
        """
        Internal method to score the likelihood difference between reference and alternate sequences.
        This is the core weapon of the Oracle.
        """
        if self.tier == "production":
            print(f"🧬 Scoring sequences (Production Tier): REF={ref_sequence[:30]}... ALT={alt_sequence[:30]}...")
            log_likelihoods = self.model.score_sequences([ref_sequence, alt_sequence])
            delta_score = float(log_likelihoods[1] - log_likelihoods[0])
            
            print(f"📊 Zeta Score (Δ): {delta_score}")
            return {
                "zeta_score": delta_score,
                "status": "success",
                "reference_likelihood": float(log_likelihoods[0]),
                "alternate_likelihood": float(log_likelihoods[1])
            }
        else: # Testing tier
            print(f"🧬 Scoring sequences (Testing Tier): REF={ref_sequence[:30]}... ALT={alt_sequence[:30]}...")
            # Simple heuristic: count differences
            min_len = min(len(ref_sequence), len(alt_sequence))
            differences = sum(1 for i in range(min_len) if ref_sequence[i] != alt_sequence[i])
            len_diff = abs(len(ref_sequence) - len(alt_sequence))
            total_changes = differences + len_diff
            mock_delta = -0.1 * total_changes
            print(f"📊 Zeta Score (Heuristic): {mock_delta}")
            return {
                "zeta_score": mock_delta,
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

    @modal.asgi_app()
    def api(self):
        """A single, unified endpoint to handle all Oracle tasks (scoring and generation)."""
        from fastapi import FastAPI
        from fastapi.responses import JSONResponse
        
        app = FastAPI()

        @app.post("/invoke")
        def invoke(item: dict):
            action = item.get("action")
            params = item.get("params", {})

            if not action:
                return JSONResponse(content={"status": "error", "message": "No 'action' specified. Must be 'score', 'generate', or 'embed'."}, status_code=400)

            try:
                if action == "score":
                    ref_seq = params.get("reference_sequence", "").strip().upper()
                    alt_seq = params.get("alternate_sequence", "").strip().upper()
                    
                    if not ref_seq or not alt_seq:
                         return JSONResponse(content={"status": "error", "message": "Missing reference_sequence or alternate_sequence for 'score' action"}, status_code=400)

                    result = self.score_variant(ref_seq, alt_seq)
                    
                    if result["zeta_score"] < -0.1:
                        interpretation = "Disruptive (Likely Pathogenic)"
                    elif result["zeta_score"] > 0.1:
                        interpretation = "Tolerated (Likely Benign)"
                    else:
                        interpretation = "Uncertain Significance"
                    
                    response_data = {
                        "zeta_score": result["zeta_score"],
                        "status": "success",
                        "interpretation": interpretation,
                        "reference_likelihood": result["reference_likelihood"],
                        "alternate_likelihood": result["alternate_likelihood"],
                    }
                    return JSONResponse(content=response_data, status_code=200)

                elif action == "generate":
                    prompt = params.get("prompt")
                    if not prompt:
                        return JSONResponse(content={"error": "No prompt provided for 'generate' action"}, status_code=400)
                    
                    gen_params = params.get("gen_params", {})
                    completion = self.generate_sequence(prompt, **gen_params)
                    return JSONResponse(content={"status": "success", "completion": completion}, status_code=200)

                elif action == "embed":
                    gene = params.get("gene_symbol")
                    hgvsp = params.get("variant")
                    
                    if not gene or not hgvsp:
                        return JSONResponse(content={"error": "Missing gene_symbol or variant for 'embed' action"}, status_code=400)
                    
                    try:
                        # Generate the correct sequence using internal methods
                        canonical_seq = self._get_canonical_sequence(gene)
                        mutated_seq = self._apply_hgvs_mutation(canonical_seq, hgvsp)
                    except (FileNotFoundError, ValueError) as e:
                        print(f"SEQUENCE GENERATION ERROR for {gene} {hgvsp}: {e}")
                        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=400)

                    embedding = self.embed_sequence(mutated_seq)
                    return JSONResponse(content={"status": "success", "embedding": embedding}, status_code=200)

                else:
                    return JSONResponse(content={"status": "error", "message": f"Invalid action '{action}' specified."}, status_code=400)

            except Exception as e:
                print(f"💥 ZETA ORACLE CRITICAL ERROR: {e}")
                import traceback
                print(f"📋 Full traceback: {traceback.format_exc()}")
                return JSONResponse(
                    content={"status": "error", "message": f"Internal server error: {str(e)}"},
                    status_code=500
                )
        
        @app.get("/health")
        def health_check():
            """A simple health check endpoint to confirm the service is running."""
            return {"status": "healthy", "service": "zeta-oracle"}
        
        return app

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