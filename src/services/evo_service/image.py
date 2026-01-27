"""Modal image definition for Evo2 service.
Based on proven UnifiedOracle deployment recipe."""
import modal
import os
from pathlib import Path

# This is the battle-tested image definition from the proven UnifiedOracle deployment.
# CRITICAL: Follows deploy_unified_oracle.py configuration for PyTorch 2.3.0 + transformer_engine 1.13
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11"
    )
    .apt_install(
        [
            "build-essential", "cmake", "ninja-build",
            "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"
        ]
    )
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True"
    })
    # Fix for numpy 1.22.0 build (llvm-ar workaround)
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar"
    )
    # CRITICAL: Install wheel package before building any packages (fixes bdist_wheel error)
    .run_commands("pip install wheel setuptools")
    # CRITICAL: Pin PyTorch 2.3.0 BEFORE installing evo2 (proven working version)
    # Using cu121 for CUDA 12.4 compatibility (matching deploy_unified_oracle.py)
    .run_commands(
        "pip install --no-cache-dir torch==2.3.0 torchvision torchaudio "
        "--index-url https://download.pytorch.org/whl/cu121"
    )
    # Install Evo2 from source
    .run_commands(
        "git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && "
        "cd evo2 && pip install ."
    )
    # CRITICAL: transformer_engine>=2.0.0 requires cuDNN 9 which is not in Ubuntu repos
    # Downgrade to 1.13 to avoid cuDNN 9 requirement (uninstall both spellings)
    .run_commands("pip uninstall -y transformer-engine transformer_engine || true")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    # Install flash-attn with fallback strategy (CUDA 12.4 -> 12.1 -> generic)
    .run_commands(
        "bash -c \"pip install 'flash-attn-cu124==2.6.3' --no-build-isolation 2>/dev/null || "
        "pip install 'flash-attn-cu121==2.6.3' --no-build-isolation 2>/dev/null || "
        "pip install 'flash-attn==2.6.3' --no-build-isolation\""
    )
    # Install remaining dependencies
    .pip_install(
        "fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx"
    )
    # Add service files to /root so they can be imported at runtime
    # Path is relative to the deployment directory (evo_service/)
    .add_local_dir(
        str(Path(__file__).parent), 
        remote_path="/root"
    )
)
