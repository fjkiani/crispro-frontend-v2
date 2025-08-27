import modal
import os
from fastapi import FastAPI, HTTPException
from loguru import logger

# Image config cloned from main.py (known-good)
evo2_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11"
    )
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
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx")
)

app = modal.App("evo-service-7b-alt", image=evo2_image)

@app.cls(gpu="H100:1", scaledown_window=300, timeout=1800)
class Evo7BServiceAlt:
    @modal.enter()
    def load_model_and_api(self):
        from evo2 import Evo2
        model_id = os.getenv("EVO_MODEL_ID", "evo2_7b")
        logger.info(f"🚀 Loading Evo2 model: {model_id} ...")
        self.model = Evo2(model_id)
        logger.info("🎉 Evo2 7B model loaded successfully!")
        self.fastapi_app = FastAPI(title="Evo2 7B Scoring Service (Alt)")

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

    @modal.asgi_app()
    def api(self):
        return self.fastapi_app 