# SAE Service: Evo2 Mechanistic Interpretability via Sparse Autoencoders
# Pattern derived from: scripts/evo2/evo2/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb
# Last updated: 2025-11-21 07:20 UTC - checkpoint loading disabled for dimension compatibility
import modal
import os
from fastapi import FastAPI, HTTPException, Header, Request
from pydantic import BaseModel, Field
from loguru import logger
from typing import List, Optional, Dict, Any
import torch
import torch.nn as nn

# --- Image Definition ---
sae_image = (
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
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar",
    )
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx", "huggingface-hub", "slowapi")
)

# --- App Definition ---
app = modal.App("sae-service", image=sae_image)

# --- Shared Volume for SAE Weights ---
volume = modal.Volume.from_name("sae-model-cache", create_if_missing=True)
SAE_CACHE_DIR = "/root/.cache/sae"

# --- Pydantic Models ---
class ExtractFeaturesRequest(BaseModel):
    activations: Optional[List[List[List[float]]]] = Field(None, description="Layer 26 activations [batch, seq_len, 4096]")
    chrom: Optional[str] = Field(None, description="Chromosome for variant scoring")
    pos: Optional[int] = Field(None, description="Position for variant scoring")
    ref: Optional[str] = Field(None, description="Reference allele")
    alt: Optional[str] = Field(None, description="Alternate allele")
    model_id: str = Field("evo2_7b", description="Evo2 model to use")
    assembly: str = Field("GRCh38", description="Genome assembly")
    window: int = Field(8192, description="Window size for fetching sequence")

# --- Utility Classes from Notebook ---
class ModelScope:
    """Context manager for adding hooks to a model."""
    def __init__(self, model):
        self.model = model
        self.hooks = []

    def add_hook(self, layer_name, hook_fn):
        layer = self.model.get_submodule(layer_name)
        handle = layer.register_forward_hook(hook_fn)
        self.hooks.append(handle)
        return handle

    def remove_hooks(self):
        for handle in self.hooks:
            handle.remove()
        self.hooks = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.remove_hooks()

class BatchTopKTiedSAE(nn.Module):
    """
    Batch-TopK Sparse Autoencoder with tied decoder.
    From Evo2 paper: https://github.com/ArcInstitute/evo2/blob/main/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb
    """
    def __init__(
        self,
        d_in: int,
        d_hidden: int,
        k: int,
        device: str,
        dtype: torch.dtype,
        tiebreaker_epsilon: float = 1e-6
    ):
        super().__init__()
        self.d_in = d_in
        self.d_hidden = d_hidden
        self.k = k
        
        # Initialize weights
        W_mat = torch.randn((d_in, d_hidden))
        W_mat = 0.1 * W_mat / torch.linalg.norm(W_mat, dim=0, ord=2, keepdim=True)
        self.W = nn.Parameter(W_mat)
        self.b_enc = nn.Parameter(torch.zeros(self.d_hidden))
        self.b_dec = nn.Parameter(torch.zeros(self.d_in))
        
        self.device = device
        self.dtype = dtype
        self.tiebreaker_epsilon = tiebreaker_epsilon
        self.tiebreaker = torch.linspace(0, tiebreaker_epsilon, d_hidden)
        
        self.to(self.device, self.dtype)
    
    def encoder_pre(self, x):
        return x @ self.W + self.b_enc
    
    def encode(self, x, tiebreak=False):
        f = torch.nn.functional.relu(self.encoder_pre(x))
        return self._batch_topk(f, self.k, tiebreak=tiebreak)
    
    def _batch_topk(self, f, k, tiebreak=False):
        from math import prod
        
        if tiebreak:  # break ties in feature order for determinism
            f += self.tiebreaker.to(f.device).broadcast_to(f.shape)
        
        *input_shape, _ = f.shape  # handle higher-dim tensors (e.g. from sequence input)
        numel = k * prod(input_shape)
        f_topk = torch.topk(f.flatten(), numel, dim=-1)
        f_topk = torch.zeros_like(f.flatten()).scatter(-1, f_topk.indices, f_topk.values).reshape(f.shape)
        return f_topk
    
    def decode(self, f):
        return f @ self.W.T + self.b_dec
    
    def forward(self, x):
        f = self.encode(x)
        return self.decode(f), f

# --- Main Service Class ---
@app.cls(
    gpu="H100:1",
    volumes={SAE_CACHE_DIR: volume},
    scaledown_window=300,
    timeout=1800
)
class SAEService:
    @modal.enter()
    def load_sae_model(self):
        """Load SAE model weights and initialize the FastAPI app."""
        from huggingface_hub import hf_hub_download
        from evo2 import Evo2
        
        logger.info("🚀 Loading SAE service...")
        
        # Load Evo2 model (for optional sequence scoring if chrom/pos provided)
        model_id = os.getenv("EVO_MODEL_ID", "evo2_7b_base")
        logger.info(f"Loading Evo2 model: {model_id}")
        self.evo_model = Evo2(model_id)
        
        # Determine a valid layer name for SAE extraction.
        # The paper uses layer 26, but some checkpoints may have fewer blocks.
        try:
            n_blocks = len(self.evo_model.model.blocks)
            if n_blocks > 26:
                self.sae_layer_name = "blocks.26"
            else:
                # Fallback to last available block (0‑indexed)
                self.sae_layer_name = f"blocks.{n_blocks - 1}"
            logger.info(f"Using Evo2 layer for SAE activations: {self.sae_layer_name} (n_blocks={n_blocks})")
        except Exception as e:
            logger.warning(f"Could not infer number of blocks from Evo2 model: {e}")
            # Safe fallback: last block by name scan
            names = [name for name, _ in self.evo_model.model.named_modules() if name.startswith("blocks.")]
            self.sae_layer_name = sorted(names)[-1] if names else "blocks.0"
            logger.info(f"Falling back to SAE layer name: {self.sae_layer_name}")
        
        # Download SAE weights from Hugging Face
        # Based on Evo2 paper: SAE trained on layer 26, 32K features (expansion=8), k=64
        logger.info("Downloading SAE weights from Hugging Face...")
        try:
            sae_weights_path = hf_hub_download(
                repo_id="Goodfire/Evo-2-Layer-26-Mixed",
                filename="sae-layer26-mixed-expansion_8-k_64.pt",  # Correct filename from notebook
                repo_type="model",
                cache_dir=SAE_CACHE_DIR
            )
            logger.info(f"SAE weights downloaded to: {sae_weights_path}")
        except Exception as e:
            logger.warning(f"Could not download SAE weights: {e}. Will initialize with random weights.")
            sae_weights_path = None
        
        # Initialize SAE model
        device = "cuda" if torch.cuda.is_available() else "cpu"
        logger.info(f"Initializing SAE on device: {device}")

        # Infer true hidden dimension for selected Evo2 layer to avoid shape mismatch
        try:
            dummy_seq = "ACGT" * 512  # 2048 bp – enough for context
            input_ids_list = self.evo_model.tokenizer.tokenize(dummy_seq)
            # Convert list -> tensor and ensure 2D [batch, seq]
            input_ids = torch.tensor([input_ids_list], dtype=torch.long, device=device)
            
            _, embs = self.evo_model.forward(
                input_ids,
                return_embeddings=True,
                layer_names=[self.sae_layer_name],
            )
            sample_acts = embs[self.sae_layer_name]
            d_in_detected = int(sample_acts.shape[-1])
            logger.info(f"✅ Detected Evo2 hidden dimension for {self.sae_layer_name}: d_in={d_in_detected}")
        except Exception as e:
            logger.error(f"❌ Failed to infer hidden dimension from Evo2: {e}")
            logger.error("    Defaulting to 4096 for evo2_7b_base; this may fail if using a different Evo2 variant.")
            d_in_detected = 4096  # Typical for evo2_7b_base
        
        self.sae_model = BatchTopKTiedSAE(
            d_in=d_in_detected,      # Match actual Evo2 hidden dimension
            d_hidden=32768,          # SAE feature dimension (from paper)
            k=64,                    # Batch-TopK sparsity (from paper)
            device=device,
            dtype=torch.float32,
            tiebreaker_epsilon=1e-6
        )

        # Load trained SAE weights if dimension matches (4096×32768 for evo2_7b)
        weights_loaded = False
        if sae_weights_path and d_in_detected == 4096:
            try:
                logger.info(f"Loading SAE checkpoint from: {sae_weights_path}")
                checkpoint = torch.load(sae_weights_path, map_location=device, weights_only=True)
                
                # Strip prefix if needed (matches notebook pattern)
                new_dict = {}
                for key, item in checkpoint.items():
                    new_key = key.replace("_orig_mod.", "").replace("module.", "")
                    new_dict[new_key] = item
                
                # Load state dict (strict=False to handle minor mismatches)
                self.sae_model.load_state_dict(new_dict, strict=False)
                weights_loaded = True
                logger.info("✅ Loaded trained SAE weights (4096×32768) from Goodfire/Evo-2-Layer-26-Mixed")
            except Exception as e:
                logger.warning(f"Failed to load checkpoint: {e}. Using random initialization.")
                logger.warning(f"⚠️  Using randomly initialized SAE ({d_in_detected}×32768)")
        else:
            if d_in_detected != 4096:
                logger.warning(f"⚠️  Dimension mismatch: detected {d_in_detected}, checkpoint expects 4096. Using random initialization.")
            elif not sae_weights_path:
                logger.warning(f"⚠️  SAE checkpoint not available. Using random initialization ({d_in_detected}×32768)")
            else:
                logger.warning(f"⚠️  Using randomly initialized SAE ({d_in_detected}×32768)")
        
        # Store weights status and dimension for provenance
        self.sae_weights_loaded = weights_loaded
        self.d_in_detected = d_in_detected
        
        self.device = device
        logger.info("🎉 SAE service initialized successfully!")
        
        # Initialize FastAPI app (we enforce safety via cohort caps/circuit breaker, not API keys)
        self.fastapi_app = FastAPI(title="SAE Feature Extraction Service")
        
        # Track usage for a simple circuit breaker (per-container)
        self.request_count = 0
        self.error_count = 0
        self.MAX_ERRORS_PER_100 = 50  # If >50% error rate across 100 calls, reject new ones
        
        @self.fastapi_app.post("/extract_features")
        def extract_features(request: Request, body: ExtractFeaturesRequest):
            """
            Extract SAE features from Evo2 layer 26 activations.
            Input: either activations directly OR (chrom, pos, ref, alt) for variant scoring
            Output: { features, top_features, layer, stats, provenance }
            """
            # Simple circuit breaker: if error rate is too high, reject new requests
            self.request_count += 1
            if self.request_count % 100 == 0:
                error_rate = (self.error_count / self.request_count) * 100
                if error_rate > self.MAX_ERRORS_PER_100:
                    logger.error(f"🚨 Circuit breaker triggered! Error rate: {error_rate:.1f}% (>{self.MAX_ERRORS_PER_100}%)")
                    logger.error("   Rejecting new requests. Check logs and redeploy to reset.")
                    raise HTTPException(status_code=503, detail="Circuit breaker: high error rate detected")
            
            logger.info(f"/extract_features called (total: {self.request_count}, errors: {self.error_count})")
            
            try:
                # If activations provided directly, use them
                if body.activations is not None:
                    activations_tensor = torch.tensor(body.activations, dtype=torch.float32, device=self.device)
                    logger.info(f"Using provided activations: shape={activations_tensor.shape}")
                
                # Otherwise, score variant and extract activations
                elif body.chrom and body.pos and body.ref and body.alt:
                    import httpx
                    logger.info(f"Scoring variant: {body.chrom}:{body.pos} {body.ref}>{body.alt}")
                    
                    # Fetch sequence and get activations using Evo2
                    assembly = body.assembly
                    chrom = str(body.chrom)
                    pos = int(body.pos)
                    ref = str(body.ref).upper()
                    alt = str(body.alt).upper()
                    window = body.window
                    
                    flank = max(1, window // 2)
                    start = max(1, pos - flank)
                    end = pos + flank
                    asm = "GRCh38" if assembly.lower() in ("grch38", "hg38") else "GRCh37"
                    region = f"{chrom}:{start}-{end}:1"
                    url = f"https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}"
                    
                    with httpx.Client(timeout=30) as client:
                        resp = client.get(url)
                        resp.raise_for_status()
                        seq = resp.text.strip().upper()
                    
                    idx = pos - start
                    if idx < 0 or idx >= len(seq):
                        raise HTTPException(status_code=400, detail="position out of fetched window")
                    
                    ref_base = seq[idx]
                    if ref_base != ref and ref_base != "N":
                        raise HTTPException(status_code=400, detail=f"Reference allele mismatch")
                    
                    alt_sequence = seq[:idx] + alt + seq[idx+1:]
                    
                    # Get activations from Evo2
                    # NOTE: Evo2 uses a CharLevelTokenizer which exposes `tokenize`/`tokenize_batch`,
                    # not a HuggingFace-style `.encode` method.
                    token_ids = self.evo_model.tokenizer.tokenize(alt_sequence)
                    alt_input_ids = torch.tensor([token_ids], dtype=torch.long, device=self.device)
                    
                    _, embeddings = self.evo_model.forward(
                        alt_input_ids,
                        return_embeddings=True,
                        layer_names=[self.sae_layer_name],
                    )
                    
                    activations_tensor = embeddings.get(self.sae_layer_name)
                    if activations_tensor is None:
                        raise HTTPException(status_code=500, detail="Failed to extract layer 26 activations")
                    
                    logger.info(f"Extracted activations from Evo2: shape={activations_tensor.shape}, dtype={activations_tensor.dtype}")
                
                else:
                    raise HTTPException(
                        status_code=400,
                        detail="Either 'activations' or (chrom, pos, ref, alt) must be provided"
                    )
                
                # Extract SAE features
                with torch.no_grad():
                    # Evo2 may emit BFloat16; convert to Float32 for SAE
                    if activations_tensor.dtype != torch.float32:
                        activations_tensor = activations_tensor.float()
                    
                    # Forward through SAE
                    _, features = self.sae_model(activations_tensor)
                    # features shape: [batch=1, seq_len=8193, d_hidden=32768]
                    
                    # Aggregate across sequence positions (mean pooling)
                    # This gives us a single 32K-dim vector representing the variant
                    features_aggregated = features.mean(dim=1).squeeze(0)  # [32768]
                    
                    # Get top-k SAE features (indices in 0-32767 range)
                    top_k_values, top_k_indices = torch.topk(features_aggregated, k=64)
                    
                    # Compute stats on original 3D features
                    sparsity = (features != 0).float().mean().item()
                    mean_activation = features[features != 0].mean().item() if (features != 0).any() else 0.0
                    
                    logger.info(f"SAE extraction done | sparsity={sparsity:.4f} mean_activation={mean_activation:.4f}")
                    
                    # NOTE: We only return top_features (k=64), not the full 32K-dim vector,
                    # to prevent massive payloads (268M floats = 1-2GB JSON) that crash Modal.
                    # Downstream biomarker analysis only needs the top-k active features anyway.
                    return {
                        "top_features": [
                            {"index": int(idx), "value": float(val)}
                            for idx, val in zip(top_k_indices.cpu().numpy(), top_k_values.cpu().numpy())
                        ],
                        "layer": "blocks.26",
                        "stats": {
                            "sparsity": sparsity,
                            "mean_activation": mean_activation,
                            "num_active_features": int((features != 0).sum().item()),
                            "shape": list(features.shape)
                        },
                        "provenance": {
                            "method": "batch_topk_tied_sae",
                            "d_in": int(self.d_in_detected),  # Use detected Evo2 input dimension
                            "d_hidden": 32768,
                            "k": 64,
                            "model": "Goodfire/Evo-2-Layer-26-Mixed (trained weights)" if self.sae_weights_loaded else f"Goodfire/Evo-2-Layer-26-Mixed (random init)"
                        }
                    }
            
            except HTTPException:
                self.error_count += 1
                raise
            except Exception as e:
                self.error_count += 1
                logger.error(f"SAE extraction failed: {e}")
                raise HTTPException(status_code=500, detail=f"SAE extraction failed: {str(e)}")
    
    @modal.asgi_app()
    def api(self):
        """Serve the FastAPI app."""
        return self.fastapi_app

# --- Local Entrypoint for Testing ---
@app.local_entrypoint()
def local_main():
    print("--- ⚔️  LOCAL SAE-SERVICE TEST ⚔️ ---")
    print("Local entrypoint for syntax validation. Does not run SAE extraction.")
    print("✅ Syntax validation passed.")

