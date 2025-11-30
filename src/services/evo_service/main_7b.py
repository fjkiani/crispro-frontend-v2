import modal
import os
from fastapi import FastAPI, HTTPException
from loguru import logger

# --- Image Definition (same as main) ---
evo2_image = (
    modal.Image.from_registry("nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install(["build-essential", "cmake", "ninja-build", "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"])
    .env({
        "CC": "/usr/bin/gcc",
        "CXX": "/usr/bin/g++",
        "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
        "EVO_MODEL_ID": "evo2_7b",
    })
    .run_commands("mkdir -p /tools/llvm/bin", "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar")
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .run_commands("bash -lc \"pip install 'flash-attn-cu124==2.6.3' || pip install 'flash-attn-cu121==2.6.3' || pip install 'flash-attn==2.6.3' --no-build-isolation\"")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx")
)

app = modal.App("evo-service-7b", image=evo2_image)

@ app.cls(gpu="A10G", scaledown_window=300, timeout=1800)
class Evo7BService:
    @modal.enter()
    def load_model_and_api(self):
        from evo2 import Evo2
        model_id = os.getenv("EVO_MODEL_ID", "evo2_7b")
        logger.info(f"🚀 Loading Evo2 model: {model_id} ...")
        self.model = Evo2(model_id)
        logger.info("🎉 Evo2 7B model loaded successfully!")
        self.fastapi_app = FastAPI(title="Evo2 7B Scoring Service")

        @self.fastapi_app.post("/score_delta")
        def score_delta(item: dict):
            ref = item.get("ref_sequence", ""); alt = item.get("alt_sequence", "")
            if not ref or not alt:
                raise HTTPException(status_code=400, detail="ref_sequence and alt_sequence are required")
            ll = self.model.score_sequences([ref, alt])
            return {"ref_likelihood": float(ll[0]), "alt_likelihood": float(ll[1]), "delta_score": float(ll[1]) - float(ll[0])}

        @self.fastapi_app.post("/score_variant")
        def score_variant(item: dict):
            import httpx
            asm_in = item.get("assembly", "GRCh38").lower(); chrom = str(item.get("chrom")); pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper(); alt = str(item.get("alt")).upper(); window = int(item.get("window", 8192))
            flank = max(1, window // 2); start = max(1, pos - flank); end = pos + flank
            asm = "GRCh38" if asm_in in ("grch38", "hg38") else "GRCh37"
            url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{start}-{end}:1?content-type=text/plain;coord_system_version={asm}"
            with httpx.Client(timeout=30) as client:
                r = client.get(url); r.raise_for_status(); seq = r.text.strip().upper()
            idx = pos - start
            if idx < 0 or idx >= len(seq):
                raise HTTPException(status_code=400, detail="position out of fetched window")
            ref_base = seq[idx]
            if ref_base != ref and ref_base != "N":
                raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
            ref_seq = seq; alt_seq = seq[:idx] + alt + seq[idx+1:]
            if ref_seq == alt_seq:
                raise HTTPException(status_code=400, detail=f"No change applied: sequences identical at {chrom}:{pos}")
            ll = self.model.score_sequences([ref_seq, alt_seq])
            return {"ref_likelihood": float(ll[0]), "alt_likelihood": float(ll[1]), "delta_score": float(ll[1]) - float(ll[0])}

        @self.fastapi_app.post("/score_variant_multi")
        def score_variant_multi(item: dict):
            import httpx
            asm_in = item.get("assembly", "GRCh38").lower(); chrom = str(item.get("chrom")); pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper(); alt = str(item.get("alt")).upper(); windows = item.get("windows") or [1024,2048,4096,8192]
            deltas = []
            with httpx.Client(timeout=30) as client:
                for w in windows:
                    flank = max(1, int(w)//2); start = max(1,pos-flank); end = pos+flank
                    asm = "GRCh38" if asm_in in ("grch38","hg38") else "GRCh37"
                    url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{start}-{end}:1?content-type=text/plain;coord_system_version={asm}"
                    r = client.get(url); r.raise_for_status(); seq = r.text.strip().upper()
                    idx = pos - start
                    if idx < 0 or idx >= len(seq):
                        raise HTTPException(status_code=400, detail="position out of fetched window")
                    if seq[idx] != ref and seq[idx] != "N":
                        raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{seq[idx]}' provided='{ref}' at {chrom}:{pos}")
                    ref_seq = seq; alt_seq = seq[:idx] + alt + seq[idx+1:]
                    if ref_seq == alt_seq:
                        raise HTTPException(status_code=400, detail=f"No change applied: sequences identical at {chrom}:{pos}")
                    ll = self.model.score_sequences([ref_seq, alt_seq])
                    deltas.append({"window": int(w), "delta": float(ll[1]) - float(ll[0])})
            min_entry = min(deltas, key=lambda d: d["delta"]) if deltas else {"window": None, "delta": 0.0}
            return {"deltas": deltas, "min_delta": min_entry["delta"], "window_used": min_entry["window"]}

        @self.fastapi_app.post("/score_variant_exon")
        def score_variant_exon(item: dict):
            import httpx
            asm_in = item.get("assembly", "GRCh38").lower(); chrom = str(item.get("chrom")); pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper(); alt = str(item.get("alt")).upper(); flank = int(item.get("flank", 600))
            start = max(1, pos - flank); end = pos + flank
            asm = "GRCh38" if asm_in in ("grch38", "hg38") else "GRCh37"
            url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{start}-{end}:1?content-type=text/plain;coord_system_version={asm}"
            with httpx.Client(timeout=30) as client:
                r = client.get(url); r.raise_for_status(); seq = r.text.strip().upper()
            idx = pos - start
            if idx < 0 or idx >= len(seq):
                raise HTTPException(status_code=400, detail="position out of fetched window")
            if seq[idx] != ref and seq[idx] != "N":
                raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{seq[idx]}' provided='{ref}' at {chrom}:{pos}")
            ref_seq = seq; alt_seq = seq[:idx] + alt + seq[idx+1:]
            if ref_seq == alt_seq:
                raise HTTPException(status_code=400, detail=f"No change applied: sequences identical at {chrom}:{pos}")
            ll = self.model.score_sequences([ref_seq, alt_seq])
            return {"exon_delta": float(ll[1]) - float(ll[0]), "window_used": flank * 2}

        @self.fastapi_app.post("/score_variant_with_activations")
        def score_variant_with_activations(item: dict):
            """
            Score a variant AND return layer 26 activations for SAE extraction (7B).
            Request: { assembly, chrom, pos, ref, alt, window?: int, return_activations?: bool }
            Response: { ref_likelihood, alt_likelihood, delta_score, layer_26_activations?: array, window_start, window_end }
            """
            import httpx
            import torch
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            window = int(item.get("window", 8192))
            return_activations = item.get("return_activations", True)
            logger.info(f"/score_variant_with_activations (7B) called | {chrom}:{pos} {ref}>{alt} window={window} return_activations={return_activations}")
            
            if alt == ref:
                raise HTTPException(status_code=400, detail="ref and alt must differ")
            
            flank = max(1, window // 2)
            start = max(1, pos - flank)
            end = pos + flank
            asm = "GRCh38" if assembly.lower() in ("grch38", "hg38") else "GRCh37"
            region = f"{chrom}:{start}-{end}:1"
            url = f"https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}"
            
            try:
                with httpx.Client(timeout=30) as client:
                    resp = client.get(url)
                    resp.raise_for_status()
                    seq = resp.text.strip().upper()
            except Exception as e:
                logger.error(f"Reference fetch failed: {e}")
                raise HTTPException(status_code=502, detail=f"Failed to fetch reference: {e}")
            
            idx = pos - start
            if idx < 0 or idx >= len(seq):
                raise HTTPException(status_code=400, detail="position out of fetched window")
            ref_base = seq[idx]
            if ref_base != ref and ref_base != "N":
                raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
            
            ref_sequence = seq
            alt_sequence = seq[:idx] + alt + seq[idx+1:]
            
            if return_activations:
                try:
                    # Tokenize sequences
                    ref_tokens = self.model.tokenizer.tokenize(ref_sequence)
                    alt_tokens = self.model.tokenizer.tokenize(alt_sequence)
                    if isinstance(ref_tokens, list):
                        ref_tokens = torch.tensor(ref_tokens, dtype=torch.long).unsqueeze(0)
                    if isinstance(alt_tokens, list):
                        alt_tokens = torch.tensor(alt_tokens, dtype=torch.long).unsqueeze(0)
                    
                    # Get activations from layer 26
                    with torch.no_grad():
                        _, ref_embs = self.model.forward(ref_tokens, return_embeddings=True, layer_names=["blocks.26"])
                        _, alt_embs = self.model.forward(alt_tokens, return_embeddings=True, layer_names=["blocks.26"])
                    
                    ref_acts = ref_embs.get("blocks.26", None)
                    alt_acts = alt_embs.get("blocks.26", None)
                    
                    if ref_acts is not None and alt_acts is not None:
                        # Take mean over sequence dimension if needed
                        if len(ref_acts.shape) > 1:
                            ref_acts = ref_acts.mean(dim=1).squeeze().cpu().numpy().tolist()
                        else:
                            ref_acts = ref_acts.squeeze().cpu().numpy().tolist()
                        if len(alt_acts.shape) > 1:
                            alt_acts = alt_acts.mean(dim=1).squeeze().cpu().numpy().tolist()
                        else:
                            alt_acts = alt_acts.squeeze().cpu().numpy().tolist()
                        
                        # Score sequences
                        ll = self.model.score_sequences([ref_sequence, alt_sequence])
                        delta = float(ll[1]) - float(ll[0])
                        
                        logger.info(f"/score_variant_with_activations (7B) done | delta={delta:.4f} | activations extracted")
                        return {
                            "ref_likelihood": float(ll[0]),
                            "alt_likelihood": float(ll[1]),
                            "delta_score": delta,
                            "layer_26_activations": {
                                "ref": ref_acts,
                                "alt": alt_acts
                            },
                            "window_start": start,
                            "window_end": end
                        }
                    else:
                        # Fallback to scoring without activations
                        ll = self.model.score_sequences([ref_sequence, alt_sequence])
                        delta = float(ll[1]) - float(ll[0])
                        logger.warning(f"/score_variant_with_activations (7B) | activations not available, returning scores only")
                        return {
                            "ref_likelihood": float(ll[0]),
                            "alt_likelihood": float(ll[1]),
                            "delta_score": delta,
                            "window_start": start,
                            "window_end": end
                        }
                except Exception as e:
                    logger.error(f"Activation extraction failed: {e}, falling back to scoring only")
                    ll = self.model.score_sequences([ref_sequence, alt_sequence])
                    delta = float(ll[1]) - float(ll[0])
                    return {
                        "ref_likelihood": float(ll[0]),
                        "alt_likelihood": float(ll[1]),
                        "delta_score": delta,
                        "window_start": start,
                        "window_end": end
                    }
            else:
                ll = self.model.score_sequences([ref_sequence, alt_sequence])
                delta = float(ll[1]) - float(ll[0])
                logger.info(f"/score_variant_with_activations (7B) done | delta={delta:.4f} | no activations")
                return {
                    "ref_likelihood": float(ll[0]),
                    "alt_likelihood": float(ll[1]),
                    "delta_score": delta,
                    "window_start": start,
                    "window_end": end
                }

    @modal.asgi_app()
    def api(self):
        return self.fastapi_app 