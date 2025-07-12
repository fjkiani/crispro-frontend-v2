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
    .pip_install("fastapi", "numpy==1.22.0", "func-timeout")
)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
requirements_path = os.path.join(script_dir, "requirements.txt")

# Add the requirements from the correctly located file if it exists
if os.path.exists(requirements_path):
    evo2_image = evo2_image.pip_install_from_requirements(requirements_path)

# --- App Definition ---
app = modal.App("zeta-oracle-service", image=evo2_image)

# --- Volume for Caching ---
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="H100:2", volumes={mount_path: volume}, scaledown_window=600, timeout=1200)
class ZetaOracleService:
    @modal.enter()
    def load_model(self):
        """Load the Evo2 model into memory when the container starts."""
        try:
            from evo2 import Evo2
            print("🚀 Loading Evo2 model for Zeta Oracle...")
            # Use the 40B base model for maximum accuracy
            self.model = Evo2('evo2_40b')
            print("🎉 Zeta Oracle model loaded successfully! 💥")
            self.model_loaded = True
            
        except Exception as e:
            print(f"💥 Error during model loading: {e}")
            import traceback
            print(f"📋 Full traceback: {traceback.format_exc()}")
            
            # Set fallback state
            self.model = None
            self.model_loaded = False
            self.model_error = str(e)
            print("⚠️ Operating in fallback mode - API will return mock responses")

    @modal.method()
    def generate_variant(self, prompt_sequence: str, gen_params: dict) -> dict:
        """Generate a new sequence based on a prompt and generation parameters."""
        try:
            if not self.model_loaded or self.model is None:
                # Fallback: return the prompt with a small modification
                print("⚠️ Using fallback generation method")
                return {
                    "completion": prompt_sequence[:-10] + "GATTACA" + "N",
                    "status": "fallback",
                    "message": f"Using fallback generation: {getattr(self, 'model_error', 'Model not available')}"
                }
            
            # Extract generation parameters with defaults
            temperature = gen_params.get("temperature", 0.7)
            top_p = gen_params.get("top_p", 0.95)
            n_tokens = gen_params.get("n_tokens", len(prompt_sequence))
            
            print(f"🧬 Generating sequence from prompt: {prompt_sequence[:20]}...")
            print(f"   - Temp: {temperature}, Top-p: {top_p}, Tokens: {n_tokens}")

            # The actual call to the Evo2 model's generate function
            completion = self.model.generate(
                prompt_sequence, 
                n_tokens=n_tokens, 
                temperature=temperature, 
                top_p=top_p
            )
            
            print(f"✅ Generated sequence: {completion[:30]}...")
            return {
                "completion": completion,
                "status": "success"
            }
            
        except Exception as e:
            print(f"❌ Error generating variant: {e}")
            import traceback
            print(f"📋 Full traceback: {traceback.format_exc()}")
            return {
                "completion": "",
                "status": "error",
                "message": str(e)
            }

    @modal.fastapi_endpoint(method="POST")
    def web(self, item: dict):
        """
        Web endpoint for the Generative Oracle.
        """
        from fastapi.responses import JSONResponse

        prompt = item.get("prompt")
        gen_params = item.get("gen_params", {})

        if not prompt:
            return JSONResponse({"status": "error", "message": "Generation requires a 'prompt' sequence"}, status_code=400)
        
        return self.generate_variant.remote(prompt, gen_params)

    @staticmethod
    @modal.fastapi_endpoint(method="GET")
    def health_check():
        """Health check endpoint."""
        return {
            "status": "healthy", 
            "service": "zeta-oracle-generative", 
            "version": "1.3.0",
            "model": "evo2_6.5b_base"
        }

# Local test entrypoint
@app.local_entrypoint()
def main():
    """
    A local test function to verify the oracle is working.
    """
    oracle = ZetaOracleService()
    
    # Test with sample sequences
    prompt_sequence = "ATCGATCGATCGATCG"
    
    print(f"🧪 Testing Zeta Oracle with prompt:")
    print(f"PROMPT: {prompt_sequence}")
    
    result = oracle.generate_variant.remote(prompt_sequence, {})
    
    print(f"\n--- Test Generation Result ---")
    print(f"Result: {result}")
    print("---------------------------") 