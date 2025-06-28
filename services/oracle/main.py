import modal
import os

# --- Image Definition ---
# Use the same battle-tested image from the working generative model
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
    # The numpy build process has a hardcoded path for the `ar` archiver.
    # We create a symlink from the system's `ar` to the path that numpy expects.
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar",
    )
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "numpy==1.22.0", "func-timeout", "python-dotenv")
)

# --- App Definition ---
app = modal.App("zeta-oracle", image=evo2_image)

# --- Volume for Caching ---
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="H100:2", volumes={mount_path: volume}, scaledown_window=600, timeout=1200)
class ZetaOracle:
    @modal.enter()
    def load_model(self):
        """
        Load the Evo2 model into memory when the container starts.
        The tier is controlled by the ORACLE_TIER environment variable.
        """
        from dotenv import load_dotenv
        load_dotenv()
        
        self.tier = os.getenv("ORACLE_TIER", "testing")
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
            print("✅ Zeta Oracle loaded in TESTING tier. Using fast, low-cost heuristics.")
            print("   To enable the full model, set ORACLE_TIER=production in your environment.")


    def _score_variant(self, ref_sequence: str, alt_sequence: str) -> dict:
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

    def _generate_sequence(self, prompt: str, **kwargs) -> str:
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
                return JSONResponse(content={"status": "error", "message": "No 'action' specified. Must be 'score' or 'generate'."}, status_code=400)

            try:
                if action == "score":
                    ref_seq = params.get("reference_sequence", "").strip().upper()
                    alt_seq = params.get("alternate_sequence", "").strip().upper()
                    
                    if not ref_seq or not alt_seq:
                         return JSONResponse(content={"status": "error", "message": "Missing reference_sequence or alternate_sequence for 'score' action"}, status_code=400)

                    result = self._score_variant(ref_seq, alt_seq)
                    
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
                    completion = self._generate_sequence(prompt, **gen_params)
                    return JSONResponse(content={"status": "success", "completion": completion}, status_code=200)

                else:
                    return JSONResponse(content={"status": "error", "message": f"Invalid action '{action}'. Must be 'score' or 'generate'."}, status_code=400)

            except Exception as e:
                print(f"💥 ZETA ORACLE CRITICAL ERROR (Action: {action}): {e}")
                import traceback
                print(f"📋 Full traceback: {traceback.format_exc()}")
                return JSONResponse(
                    content={"status": "error", "message": f"Internal server error: {str(e)}"},
                    status_code=500
                )
        
        @app.get("/health")
        def health_check():
            """Health check endpoint to confirm service is alive."""
            model_name = "evo2_40b" if self.tier == "production" else "heuristic_testing_model"
            return {
                "status": "healthy", 
                "service": "zeta-oracle", 
                "version": "3.1.0-tiered",
                "tier": self.tier,
                "model": model_name,
                "message": "I am the Oracle. I score and I create. All queries flow through me."
            }
        
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
        result = oracle._score_variant(ref_sequence, alt_sequence)
        print("\n--- ✅ Test Scoring Result ---")
        print(f"  Result: {result}")
        print("---------------------------")
    except Exception as e:
        print(f"\n--- ❌ TEST FAILED ---")
        print(f"  Error: {e}")
        print("---------------------") 