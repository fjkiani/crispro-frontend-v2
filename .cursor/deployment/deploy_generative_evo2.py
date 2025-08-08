import modal
import os

# --- Image Definition ---
# Use the same battle-tested image from the discriminator model
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.12"
    )
    .apt_install(
        ["build-essential", "cmake", "ninja-build",
            "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"]
    )
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
    })
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "func-timeout")
)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))
requirements_path = os.path.join(script_dir, "requirements.txt")

# Add the requirements from the correctly located file
if os.path.exists(requirements_path):
    evo2_image = evo2_image.pip_install_from_requirements(requirements_path)

# --- App Definition ---
# Create a new, separate app for this service
app = modal.App("evo2-sequence-generator-small", image=evo2_image)

# --- Volume for Caching ---
# Reuse the same volume to cache model weights
volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"

@app.cls(gpu="H100", volumes={mount_path: volume}, scaledown_window=300)
class Evo2Generator:
    @modal.enter()
    def load_model(self):
        """
        Load the Evo2 model into memory when the container starts.
        """
        from evo2 import Evo2
        print("Loading a smaller, more stable evo2 model for generation...")
        # Switch to the 1B parameter base model to ensure stability.
        self.model = Evo2('evo2_1b_base')
        print("Evo2 model loaded.")

    @modal.method()
    def generate(self, prompt: str):
        """
        Generates a DNA sequence completion for a given prompt.
        """
        print(f"Generating sequence for prompt: '{prompt}'")
        output = self.model.generate([prompt])
        
        if output and output.sequences and output.sequences[0]:
            completion = output.sequences[0]
            print(f"Generated sequence: {completion}")
            return completion
        else:
            print("Generation failed to produce output.")
            return None

    @modal.fastapi_endpoint(method="POST")
    def web(self, item: dict):
        """
        A web endpoint to handle generation requests.
        """
        from fastapi.responses import JSONResponse

        prompt = item.get("prompt")
        if not prompt:
            return JSONResponse(content={"error": "No prompt provided"}, status_code=400)

        # Use .call() for a direct, blocking invocation from within the same class.
        # This is the recommended pattern for this use case.
        # It avoids the `TypeError: object str can't be used in 'await' expression`
        # that occurs when using .remote() from an async FastAPI endpoint.
        completion = self.generate.call(prompt)
        
        return {"completion": completion}

@app.local_entrypoint()
def main():
    """
    A local test function to verify the generator is working.
    """
    model = Evo2Generator()
    prompt = "Complete the following DNA sequence: GCTTGAAATGTGTTAGGATG"
    print(f"Calling generator with prompt: {prompt}")
    
    # Call the 'generate' method remotely. This blocks until it completes.
    completion = model.generate.remote(prompt)
    
    print(f"\n--- Test Generation Result ---")
    print(f"Prompt: {prompt}")
    print(f"Completion: {completion}")
    print("----------------------------") 