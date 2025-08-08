import modal
import os

# --- Image Definition ---
# Use the same battle-tested image from the discriminator model
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
    .pip_install("fastapi", "pydantic", "requests", "numpy==1.22.0")
)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
requirements_path = os.path.join(script_dir, "requirements.txt")

# Add the requirements from the correctly located file
if os.path.exists(requirements_path):
    evo2_image = evo2_image.pip_install_from_requirements(requirements_path)

# --- App Definition ---
# Create a new, separate app for this service
app = modal.App("evo2-discriminator", image=evo2_image)

# --- Volume for Caching ---
# Reuse the same volume to cache model weights
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="A100", volumes={mount_path: volume}, timeout=600)
class Evo2Discriminator:
    @modal.enter()
    def load_model(self):
        """
        Load the Evo2 model into memory when the container starts.
        """
        from evo2 import Evo2
        print("Loading Evo2 Biological Foundation Model for discrimination...")
        # Using the 7B model as it's more powerful for scoring
        self.model = Evo2('evo2_7b_base')
        print("Evo2 7B model loaded.")

    @modal.method()
    def score_sequences_real(self, ref_seq: str, alt_seq: str) -> dict:
        """Score sequences using the real Evo2 model's native score method."""
        try:
            print("Attempting to score reference sequence...")
            ref_likelihood_raw = self.model.score_sequences([ref_seq])[0]
            print(f"  Ref likelihood raw value: {ref_likelihood_raw}, type: {type(ref_likelihood_raw)}")
            
            print("Attempting to score alternate sequence...")
            alt_likelihood_raw = self.model.score_sequences([alt_seq])[0]
            print(f"  Alt likelihood raw value: {alt_likelihood_raw}, type: {type(alt_likelihood_raw)}")

            ref_likelihood = float(ref_likelihood_raw)
            alt_likelihood = float(alt_likelihood_raw)
            delta_score = alt_likelihood - ref_likelihood

            print("Scoring successful. Preparing response.")
            return {
                "ref_likelihood": ref_likelihood,
                "alt_likelihood": alt_likelihood,
                "delta_score": delta_score,
                "method": "real_evo2",
                "model": "evo2_7b_base"
            }
        except Exception as e:
            import traceback
            print(f"❌ Real model scoring failed unexpectedly.")
            print(f"  Error Type: {type(e)}")
            print(f"  Error Message: {e}")
            print("  Traceback:")
            traceback.print_exc()
            return {"error": str(e), "method": "error", "traceback": traceback.format_exc()}

    @modal.fastapi_endpoint(method="POST")
    async def web(self, item: dict):
        """
        A web endpoint to handle sequence scoring requests.
        """
        from fastapi.responses import JSONResponse

        ref_seq = item.get("reference_sequence")
        alt_seq = item.get("alternate_sequence")

        if not ref_seq or not alt_seq:
            return JSONResponse(content={"error": "Missing reference_sequence or alternate_sequence"}, status_code=400)

        # In an async endpoint, we must use `await` with the `.aio` variant of remote().
        scoring_result = await self.score_sequences_real.remote.aio(ref_seq, alt_seq)
        
        return JSONResponse(content=scoring_result)

@app.local_entrypoint()
def main():
    """
    A local test function to verify the discriminator is working.
    """
    model = Evo2Discriminator()
    ref = "ATCGATCGATCGATCGATCGATCGATCGATCG"
    alt = "ATCGATCGATCGTTCGATCGATCGATCGATCG"
    
    print(f"Calling discriminator with sequences...")
    print(f"  REF: {ref}")
    print(f"  ALT: {alt}")
    
    result = model.score_sequences_real.remote(ref, alt)
    
    print(f"\n--- Test Discrimination Result ---")
    print(result)
    print("--------------------------------") 