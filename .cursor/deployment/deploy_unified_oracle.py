# CACHE BUSTER v1.1: Forcing a rebuild to fix the GenerationOutput attribute error.
import modal
import os

# --- Image Definition ---
# A single, consolidated image for all oracle tasks.
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11"
    )
    .apt_install([
        "build-essential", "cmake", "ninja-build", "git", "gcc", "g++",
        "libcudnn8", "libcudnn8-dev"
    ])
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
    })
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar",
    )
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "pydantic", "requests", "numpy==1.22.0")
)

# --- App Definition ---
app = modal.App("unified-oracle-service", image=evo2_image)

# --- Volume for Caching ---
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="H100:2", volumes={mount_path: volume}, timeout=900, scaledown_window=300)
class UnifiedOracle:
    @modal.enter()
    def load_model(self):
        """
        Load the single Evo2 model for all tasks.
        """
        from evo2 import Evo2
        print("✅ Loading Unified Evo2 Biological Foundation Model...")
        # The 7B model is a strong balance of power and reliability for both tasks.
        self.model = Evo2('evo2_7b_base')
        print("🎉 Unified Evo2 7B model loaded successfully.")

    # DOCTRINAL CORRECTION: These are internal helper methods, not standalone Modal Functions.
    # The @modal.method decorator is not needed and was the source of the call failures.
    def score_variant(self, ref_sequence: str, alt_sequence: str) -> dict:
        """Scores the likelihood difference between reference and alternate sequences."""
        try:
            log_likelihoods = self.model.score_sequences([ref_sequence, alt_sequence])
            delta_score = float(log_likelihoods[1] - log_likelihoods[0])
            
            print(f"📊 Delta score: {delta_score}")
            return {
                "delta_score": delta_score,
                "status": "success",
                "ref_likelihood": float(log_likelihoods[0]),
                "alt_likelihood": float(log_likelihoods[1])
            }
        except Exception as e:
            print(f"❌ Error scoring variant: {e}")
            return {"status": "error", "message": str(e)}

    def generate_variant(self, prompt_sequence: str, gen_params: dict) -> dict:
        """Generates a new sequence based on a prompt and generation parameters."""
        try:
            temperature = gen_params.get("temperature", 0.7)
            top_p = gen_params.get("top_p", 0.95)
            n_tokens = gen_params.get("n_tokens", 2048) # Default to a reasonable number of tokens
            
            print(f"🧬 Generating sequence from prompt: {prompt_sequence[:30]}...")
            print(f"   - Temp: {temperature}, Top-p: {top_p}, Tokens: {n_tokens}")

            completion_obj = self.model.generate(
                [prompt_sequence], # CRITICAL FIX: Pass the prompt as a list to treat it as a single batch item.
                n_tokens=n_tokens, 
                temperature=temperature, 
                top_p=top_p
            )
            
            # DOCTRINAL CORRECTION v2: The model returns a GenerationOutput object.
            # The correct attribute is .sequences, which is a list of strings.
            # We will take the first element as the primary candidate.
            completion_str = completion_obj.sequences[0]

            print(f"✅ Generated sequence: {completion_str[:30]}...")
            return {"completion": completion_str, "status": "success"}
        except Exception as e:
            print(f"❌ Error generating variant: {e}")
            return {"completion": "", "status": "error", "message": str(e)}

    @modal.fastapi_endpoint(method="POST")
    def invoke(self, item: dict):
        """
        Unified invocation endpoint for the Oracle.
        Routes requests to the appropriate method based on the 'action' field.
        """
        from fastapi.responses import JSONResponse

        action = item.get("action")
        params = item.get("params", {})

        if not action:
            return JSONResponse({"status": "error", "message": "'action' field is required"}, status_code=400)
        
        print(f"Received action: {action}")

        if action == "score":
            ref_seq = params.get("reference_sequence")
            alt_seq = params.get("alternate_sequence")
            if not ref_seq or not alt_seq:
                return JSONResponse({"status": "error", "message": "Scoring requires 'reference_sequence' and 'alternate_sequence'"}, status_code=400)
            
            # Now a direct, regular Python method call.
            result = self.score_variant(ref_seq, alt_seq)
            if result.get("status") == "success":
                delta = result.get("delta_score", 0)
                if delta < -5: result["interpretation"] = "High Impact (Deleterious)"
                elif delta < -1: result["interpretation"] = "Moderate Impact"
                elif delta > 1: result["interpretation"] = "Tolerated"
                else: result["interpretation"] = "Benign / Uncertain"
            return result
        
        elif action == "generate":
            prompt = params.get("prompt")
            gen_params = params.get("gen_params", {})
            if not prompt:
                return JSONResponse({"status": "error", "message": "Generation requires a 'prompt' sequence"}, status_code=400)
            
            # Now a direct, regular Python method call.
            return self.generate_variant(prompt, gen_params)

        else:
            return JSONResponse({"status": "error", "message": f"Unknown action: '{action}'"}, status_code=400)

    @staticmethod
    @modal.fastapi_endpoint(method="GET")
    def health():
        return {"status": "healthy", "service": "unified-oracle"}

@app.local_entrypoint()
def main():
    oracle = UnifiedOracle()
    print("--- Testing Scoring ---")
    score_res = oracle.invoke.remote({"action": "score", "params": {"reference_sequence": "ATGC", "alternate_sequence": "ATGG"}})
    print(score_res)
    
    print("\n--- Testing Generation ---")
    gen_res = oracle.invoke.remote({"action": "generate", "params": {"prompt": "ATGC", "gen_params": {"n_tokens": 10}}})
    print(gen_res) 