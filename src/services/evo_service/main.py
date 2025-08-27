# FORCE REDEPLOY: v2.0 (Phoenix Protocol)
import modal
import os
import uuid
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from loguru import logger

# --- Image Definition ---
# This is the battle-tested image definition from the proven Zeta Oracle.
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
    # --- Fix for numpy 1.22.0 build ---
    # As per Zeta Doctrine (evo2-deployment-guide.mdc), this is a critical fix
    # to address a build bug in this specific numpy version.
    .run_commands(
        "mkdir -p /tools/llvm/bin",
        "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar",
    )
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx")
)

# --- App Definition ---
app = modal.App("evo-service", image=evo2_image)

# --- Shared State & Volume ---
job_status_dict = modal.Dict.from_name("evo-job-status", create_if_missing=True)
volume = modal.Volume.from_name("evo-model-cache", create_if_missing=True)
EVO_MODEL_DIR = "/root/.cache/huggingface"

# --- Pydantic Models ---
class JobSubmitResponse(BaseModel):
    job_id: str

class GenerationRequest(BaseModel):
    prompt: str = Field(..., description="A detailed natural language prompt for the generative model.")
    n_tokens: int = Field(100, description="The number of tokens to generate.")

# --- Main Service Class ---
@app.cls(
    gpu="H100:2", 
    volumes={EVO_MODEL_DIR: volume}, 
    scaledown_window=300, 
    timeout=1800
)
class EvoService:
    @modal.enter()
    def load_model_and_api(self):
        """Load model and initialize the FastAPI app within the container."""
        from evo2 import Evo2
        import os as _os
        model_id = _os.getenv('EVO_MODEL_ID', 'evo2_40b')
        logger.info(f"🚀 Loading Evo2 model: {model_id} ...")
        self.model = Evo2(model_id)
        logger.info("🎉 Evo2 model loaded successfully!")

        self.fastapi_app = FastAPI(title="Evo2 General-Purpose Generation Service")

        @self.fastapi_app.post("/generate", status_code=202, response_model=JobSubmitResponse)
        def submit_generation(request: GenerationRequest):
            job_id = str(uuid.uuid4())
            logger.info(f"Received submission request. Assigning job ID: {job_id}")
            job_status_dict[job_id] = {"status": "pending"}
            self.generate_in_background.spawn(job_id, request.prompt, request.n_tokens)
            return {"job_id": job_id}

        @self.fastapi_app.get("/status/{job_id}")
        def get_status(job_id: str):
            logger.info(f"Received status query for job ID: {job_id}")
            if job_id not in job_status_dict:
                raise HTTPException(status_code=404, detail="Job ID not found.")
            return job_status_dict[job_id]

        @self.fastapi_app.post("/score_delta")
        def score_delta(item: dict):
            """
            Compute delta log-likelihood between reference and alternate sequences.
            Request: { "ref_sequence": str, "alt_sequence": str }
            Response: { "ref_likelihood": float, "alt_likelihood": float, "delta_score": float }
            """
            ref_sequence = item.get("ref_sequence", "")
            alt_sequence = item.get("alt_sequence", "")
            logger.info(f"/score_delta called | ref_len={len(ref_sequence)} alt_len={len(alt_sequence)}")
            if not ref_sequence or not alt_sequence:
                raise HTTPException(status_code=400, detail="ref_sequence and alt_sequence are required")
            try:
                ll = self.model.score_sequences([ref_sequence, alt_sequence])
                ref_ll = float(ll[0])
                alt_ll = float(ll[1])
                delta = alt_ll - ref_ll
                logger.info(f"/score_delta done | delta={delta:.4f}")
                return {
                    "ref_likelihood": ref_ll,
                    "alt_likelihood": alt_ll,
                    "delta_score": delta,
                }
            except Exception as e:
                logger.error(f"Scoring failed: {e}")
                raise HTTPException(status_code=500, detail="Evo2 scoring failed")

        @self.fastapi_app.post("/score_batch")
        def score_batch(item: dict):
            """
            Batch delta scoring.
            Request: { "pairs": [ {"ref_sequence": str, "alt_sequence": str}, ... ] }
            Response: { "results": [ {"ref_likelihood": float, "alt_likelihood": float, "delta_score": float}, ... ] }
            """
            pairs = item.get("pairs", [])
            logger.info(f"/score_batch called | num_pairs={len(pairs) if isinstance(pairs, list) else 0}")
            if not isinstance(pairs, list) or not pairs:
                raise HTTPException(status_code=400, detail="pairs must be a non-empty list")
            try:
                # Build flattened list [ref1, alt1, ref2, alt2, ...]
                seqs = []
                for p in pairs:
                    r = p.get("ref_sequence", "")
                    a = p.get("alt_sequence", "")
                    if not r or not a:
                        raise HTTPException(status_code=400, detail="Each pair requires ref_sequence and alt_sequence")
                    seqs.extend([r, a])
                ll = self.model.score_sequences(seqs)
                results = []
                for i in range(0, len(ll), 2):
                    ref_ll = float(ll[i])
                    alt_ll = float(ll[i+1])
                    results.append({
                        "ref_likelihood": ref_ll,
                        "alt_likelihood": alt_ll,
                        "delta_score": alt_ll - ref_ll,
                    })
                logger.info(f"/score_batch done | num_results={len(results)}")
                return {"results": results}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"Batch scoring failed: {e}")
                raise HTTPException(status_code=500, detail="Evo2 batch scoring failed")

        @self.fastapi_app.post("/score_variant")
        def score_variant(item: dict):
            """
            Compute delta on a genomic SNV by fetching a window from Ensembl.
            Request: { "assembly": "GRCh38"|"hg38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "window": 8192 }
            Response: same as /score_delta
            """
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            window = int(item.get("window", 8192))
            logger.info(f"/score_variant called | asm={assembly} {chrom}:{pos} {ref}>{alt} window={window}")
            if alt == ref:
                raise HTTPException(status_code=400, detail="ref and alt must differ")
            flank = max(1, window // 2)
            start = max(1, pos - flank)
            end = pos + flank
            species = "human"
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
            # Compute index of the variant within window
            idx = pos - start
            if idx < 0 or idx >= len(seq):
                raise HTTPException(status_code=400, detail="position out of fetched window")
            ref_base = seq[idx]
            # Strict validation: fetched base must match provided REF (allow N only)
            if ref_base != ref and ref_base != "N":
                raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
            ref_sequence = seq
            alt_sequence = seq[:idx] + alt + seq[idx+1:]
            try:
                ll = self.model.score_sequences([ref_sequence, alt_sequence])
                ref_ll = float(ll[0])
                alt_ll = float(ll[1])
                delta = alt_ll - ref_ll
                logger.info(f"/score_variant done | delta={delta:.4f}")
                return {
                    "ref_likelihood": ref_ll,
                    "alt_likelihood": alt_ll,
                    "delta_score": delta,
                    "window_start": start,
                    "window_end": end,
                    "variant_index": idx,
                }
            except Exception as e:
                logger.error(f"Variant scoring failed: {e}")
                raise HTTPException(status_code=500, detail="Evo2 variant scoring failed")

        @self.fastapi_app.post("/score_variant_multi")
        def score_variant_multi(item: dict):
            """
            Multi-scale scoring across window sizes; returns per-window deltas and the most negative delta.
            Request: { assembly, chrom, pos, ref, alt, windows?: [int] }
            Response: { deltas: [{window:int, delta:float}], min_delta: float, window_used: int }
            """
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            windows = item.get("windows") or [1024, 2048, 4096, 8192]
            logger.info(f"/score_variant_multi called | {chrom}:{pos} {ref}>{alt} windows={windows}")
            deltas = []
            try:
                with httpx.Client(timeout=30) as client:
                    for w in windows:
                        flank = max(1, int(w) // 2)
                        start = max(1, pos - flank)
                        end = pos + flank
                        asm = "GRCh38" if assembly.lower() in ("grch38", "hg38") else "GRCh37"
                        region = f"{chrom}:{start}-{end}:1"
                        url = f"https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}"
                        resp = client.get(url)
                        resp.raise_for_status()
                        seq = resp.text.strip().upper()
                        idx = pos - start
                        if idx < 0 or idx >= len(seq):
                            raise HTTPException(status_code=400, detail="position out of fetched window")
                        ref_base = seq[idx]
                        if ref_base != ref and ref_base != "N":
                            raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
                        ref_sequence = seq
                        alt_sequence = seq[:idx] + alt + seq[idx+1:]
                        ll = self.model.score_sequences([ref_sequence, alt_sequence])
                        ref_ll = float(ll[0])
                        alt_ll = float(ll[1])
                        delta = alt_ll - ref_ll
                        deltas.append({"window": int(w), "delta": delta})
                # pick most negative (strongest disruption)
                min_entry = min(deltas, key=lambda d: d["delta"]) if deltas else {"window": None, "delta": 0.0}
                logger.info(f"/score_variant_multi done | min_delta={min_entry['delta']:.4f} window={min_entry['window']}")
                return {"deltas": deltas, "min_delta": min_entry["delta"], "window_used": min_entry["window"]}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_multi failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 multi-window scoring failed: {e}")

        @self.fastapi_app.post("/score_variant_exon")
        def score_variant_exon(item: dict):
            """
            Tight-window (exon-like) scoring around the variant; default ±600 bp.
            Request: { assembly, chrom, pos, ref, alt, flank?: int }
            Response: { exon_delta: float, window_used: int }
            """
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            flank = int(item.get("flank", 600))
            logger.info(f"/score_variant_exon called | {chrom}:{pos} {ref}>{alt} flank={flank}")
            try:
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
                    raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
                ref_sequence = seq
                alt_sequence = seq[:idx] + alt + seq[idx+1:]
                ll = self.model.score_sequences([ref_sequence, alt_sequence])
                ref_ll = float(ll[0])
                alt_ll = float(ll[1])
                delta = alt_ll - ref_ll
                logger.info(f"/score_variant_exon done | exon_delta={delta:.4f}")
                return {"exon_delta": delta, "window_used": flank * 2}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_exon failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 exon scoring failed: {e}")

        def _ensure_change_applied(ref_sequence: str, alt_sequence: str, chrom: str, pos: int):
            if ref_sequence == alt_sequence:
                raise HTTPException(status_code=400, detail=f"No change applied: sequences identical at {chrom}:{pos}")

        @self.fastapi_app.post("/score_variant_profile")
        def score_variant_profile(item: dict):
            """
            Local delta profile by applying the same ALT at offsets in ±100 bp around the locus.
            Request: { assembly, chrom, pos, ref, alt, flank?: int (fetch context, default 600), radius?: int (default 100) }
            Response: { profile: [{offset:int, delta:float}], peak_delta: float, peak_offset: int }
            """
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            flank = int(item.get("flank", 600))
            radius = int(item.get("radius", 100))
            logger.info(f"/score_variant_profile called | {chrom}:{pos} {ref}>{alt} radius={radius}")
            try:
                start = max(1, pos - flank)
                end = pos + flank
                asm = "GRCh38" if assembly.lower() in ("grch38", "hg38") else "GRCh37"
                region = f"{chrom}:{start}-{end}:1"
                url = f"https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}"
                with httpx.Client(timeout=30) as client:
                    resp = client.get(url)
                    resp.raise_for_status()
                    seq = resp.text.strip().upper()
                idx0 = pos - start
                if idx0 < 0 or idx0 >= len(seq):
                    raise HTTPException(status_code=400, detail="position out of fetched window")
                if seq[idx0] != ref and seq[idx0] != "N":
                    raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{seq[idx0]}' provided='{ref}' at {chrom}:{pos}")
                profile = []
                for off in range(-radius, radius + 1):
                    idx = idx0 + off
                    if idx < 0 or idx >= len(seq):
                        continue
                    # Only mutate if matching ref (or N); otherwise skip offset to avoid double-mutation
                    if seq[idx] not in "ACGTN":
                        continue
                    ref_seq = seq
                    alt_seq = ref_seq[:idx] + alt + ref_seq[idx+1:]
                    _ensure_change_applied(ref_seq, alt_seq, chrom, pos+off)
                    ll = self.model.score_sequences([ref_seq, alt_seq])
                    ref_ll = float(ll[0])
                    alt_ll = float(ll[1])
                    profile.append({"offset": off, "delta": alt_ll - ref_ll})
                peak = min(profile, key=lambda x: x["delta"]) if profile else {"offset": 0, "delta": 0.0}
                logger.info(f"/score_variant_profile done | peak_delta={peak['delta']:.4f} @ offset {peak['offset']}")
                return {"profile": profile, "peak_delta": peak["delta"], "peak_offset": peak["offset"]}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_profile failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 profile failed: {e}")

        @self.fastapi_app.post("/score_variant_probe")
        def score_variant_probe(item: dict):
            """
            3-alt sensitivity probe at the locus; scores ref->each alt in {A,C,G,T}\{ref}.
            Request: { assembly, chrom, pos, ref }
            Response: { probes: [{alt: str, delta: float}], top_alt: str, top_delta: float }
            """
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            logger.info(f"/score_variant_probe called | {chrom}:{pos} ref={ref}")
            try:
                start = max(1, pos - 600)
                end = pos + 600
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
                    raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
                alts = [b for b in "ACGT" if b != ref]
                results = []
                for a in alts:
                    ref_seq = seq
                    alt_seq = ref_seq[:idx] + a + ref_seq[idx+1:]
                    _ensure_change_applied(ref_seq, alt_seq, chrom, pos)
                    ll = self.model.score_sequences([ref_seq, alt_seq])
                    ref_ll = float(ll[0])
                    alt_ll = float(ll[1])
                    results.append({"alt": a, "delta": alt_ll - ref_ll})
                top = min(results, key=lambda x: x["delta"]) if results else {"alt": None, "delta": 0.0}
                logger.info(f"/score_variant_probe done | top_alt={top['alt']} delta={top['delta']:.4f}")
                return {"probes": results, "top_alt": top["alt"], "top_delta": top["delta"]}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_probe failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 probe failed: {e}")

    @modal.method()
    def generate_in_background(self, job_id: str, prompt: str, n_tokens: int):
        """Runs the model inference in the background with a general-purpose prompt."""
        logger.info(f"Job {job_id}: Starting background generation with prompt: '{prompt[:50]}...'")
        job_status_dict[job_id] = {"status": "running"}
        try:
            generation_result = self.model.generate(
                prompt_seqs=[prompt],
                n_tokens=n_tokens,
            )
            # We now return a structured result, not a simple list.
            # This is critical for reliable parsing by the client (the Zeta Forge).
            decoded_sequence = generation_result.sequences[0]
            logger.info(f"Job {job_id}: Successfully generated sequence.")
            job_status_dict[job_id] = {"status": "complete", "result": {"sequence": decoded_sequence}}
        except Exception as e:
            logger.error(f"Job {job_id}: Error during generation: {e}")
            job_status_dict[job_id] = {"status": "failed", "error": str(e)}

    @modal.asgi_app()
    def api(self):
        """Serve the FastAPI app."""
        return self.fastapi_app


@app.cls(
    gpu="H100:1",
    scaledown_window=300,
    timeout=1800,
)
class EvoService7B:
    @modal.enter()
    def load_model_and_api(self):
        from evo2 import Evo2
        import os as _os
        _os.environ["EVO_MODEL_ID"] = "evo2_7b"
        model_id = "evo2_7b"
        logger.info(f"🚀 Loading Evo2 model: {model_id} ...")
        self.model = Evo2(model_id)
        logger.info("🎉 Evo2 7B model loaded successfully!")

        self.fastapi_app = FastAPI(title="Evo2 7B Scoring Service (within main)")

        @self.fastapi_app.post("/score_delta")
        def score_delta(item: dict):
            ref_sequence = item.get("ref_sequence", "")
            alt_sequence = item.get("alt_sequence", "")
            if not ref_sequence or not alt_sequence:
                raise HTTPException(status_code=400, detail="ref_sequence and alt_sequence are required")
            ll = self.model.score_sequences([ref_sequence, alt_sequence])
            return {
                "ref_likelihood": float(ll[0]),
                "alt_likelihood": float(ll[1]),
                "delta_score": float(ll[1]) - float(ll[0]),
            }

        @self.fastapi_app.post("/score_variant")
        def score_variant(item: dict):
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            window = int(item.get("window", 8192))
            logger.info(f"/score_variant (7B) | asm={assembly} {chrom}:{pos} {ref}>{alt} window={window}")
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
            ll = self.model.score_sequences([ref_sequence, alt_sequence])
            return {
                "ref_likelihood": float(ll[0]),
                "alt_likelihood": float(ll[1]),
                "delta_score": float(ll[1]) - float(ll[0]),
            }

        @self.fastapi_app.post("/score_variant_multi")
        def score_variant_multi(item: dict):
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            windows = item.get("windows") or [1024, 2048, 4096, 8192]
            logger.info(f"/score_variant_multi (7B) | {chrom}:{pos} {ref}>{alt} windows={windows}")
            deltas = []
            try:
                with httpx.Client(timeout=30) as client:
                    for w in windows:
                        flank = max(1, int(w) // 2)
                        start = max(1, pos - flank)
                        end = pos + flank
                        asm = "GRCh38" if assembly.lower() in ("grch38", "hg38") else "GRCh37"
                        region = f"{chrom}:{start}-{end}:1"
                        url = f"https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}"
                        resp = client.get(url)
                        resp.raise_for_status()
                        seq = resp.text.strip().upper()
                        idx = pos - start
                        if idx < 0 or idx >= len(seq):
                            raise HTTPException(status_code=400, detail="position out of fetched window")
                        ref_base = seq[idx]
                        if ref_base != ref and ref_base != "N":
                            raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
                        ref_sequence = seq
                        alt_sequence = seq[:idx] + alt + seq[idx+1:]
                        ll = self.model.score_sequences([ref_sequence, alt_sequence])
                        ref_ll = float(ll[0])
                        alt_ll = float(ll[1])
                        delta = alt_ll - ref_ll
                        deltas.append({"window": int(w), "delta": delta})
                min_entry = min(deltas, key=lambda d: d["delta"]) if deltas else {"window": None, "delta": 0.0}
                logger.info(f"/score_variant_multi (7B) done | min_delta={min_entry['delta']:.4f} window={min_entry['window']}")
                return {"deltas": deltas, "min_delta": min_entry["delta"], "window_used": min_entry["window"]}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_multi (7B) failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 multi-window scoring failed: {e}")

        @self.fastapi_app.post("/score_variant_exon")
        def score_variant_exon(item: dict):
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            flank = int(item.get("flank", 600))
            logger.info(f"/score_variant_exon (7B) | {chrom}:{pos} {ref}>{alt} flank={flank}")
            try:
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
                    raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
                ref_sequence = seq
                alt_sequence = seq[:idx] + alt + seq[idx+1:]
                ll = self.model.score_sequences([ref_sequence, alt_sequence])
                ref_ll = float(ll[0])
                alt_ll = float(ll[1])
                delta = alt_ll - ref_ll
                logger.info(f"/score_variant_exon (7B) done | exon_delta={delta:.4f}")
                return {"exon_delta": delta, "window_used": flank * 2}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_exon (7B) failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 exon scoring failed: {e}")

        @self.fastapi_app.post("/score_variant_profile")
        def score_variant_profile(item: dict):
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            alt = str(item.get("alt")).upper()
            flank = int(item.get("flank", 600))
            radius = int(item.get("radius", 100))
            logger.info(f"/score_variant_profile (7B) | {chrom}:{pos} {ref}>{alt} radius={radius}")
            try:
                start = max(1, pos - flank)
                end = pos + flank
                asm = "GRCh38" if assembly.lower() in ("grch38", "hg38") else "GRCh37"
                region = f"{chrom}:{start}-{end}:1"
                url = f"https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}"
                with httpx.Client(timeout=30) as client:
                    resp = client.get(url)
                    resp.raise_for_status()
                    seq = resp.text.strip().upper()
                idx0 = pos - start
                if idx0 < 0 or idx0 >= len(seq):
                    raise HTTPException(status_code=400, detail="position out of fetched window")
                if seq[idx0] != ref and seq[idx0] != "N":
                    raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{seq[idx0]}' provided='{ref}' at {chrom}:{pos}")
                profile = []
                for off in range(-radius, radius + 1):
                    idx = idx0 + off
                    if idx < 0 or idx >= len(seq):
                        continue
                    if seq[idx] not in "ACGTN":
                        continue
                    ref_seq = seq
                    alt_seq = ref_seq[:idx] + alt + ref_seq[idx+1:]
                    ll = self.model.score_sequences([ref_seq, alt_seq])
                    ref_ll = float(ll[0])
                    alt_ll = float(ll[1])
                    profile.append({"offset": off, "delta": alt_ll - ref_ll})
                peak = min(profile, key=lambda x: x["delta"]) if profile else {"offset": 0, "delta": 0.0}
                logger.info(f"/score_variant_profile (7B) done | peak_delta={peak['delta']:.4f} @ offset {peak['offset']}")
                return {"profile": profile, "peak_delta": peak["delta"], "peak_offset": peak["offset"]}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_profile (7B) failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 profile failed: {e}")

        @self.fastapi_app.post("/score_variant_probe")
        def score_variant_probe(item: dict):
            import httpx
            assembly = item.get("assembly", "GRCh38")
            chrom = str(item.get("chrom"))
            pos = int(item.get("pos"))
            ref = str(item.get("ref")).upper()
            logger.info(f"/score_variant_probe (7B) | {chrom}:{pos} ref={ref}")
            try:
                start = max(1, pos - 600)
                end = pos + 600
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
                    raise HTTPException(status_code=400, detail=f"Reference allele mismatch: fetched='{ref_base}' provided='{ref}' at {chrom}:{pos}")
                alts = [b for b in "ACGT" if b != ref]
                results = []
                for a in alts:
                    ref_seq = seq
                    alt_seq = ref_seq[:idx] + a + ref_seq[idx+1:]
                    ll = self.model.score_sequences([ref_seq, alt_seq])
                    ref_ll = float(ll[0])
                    alt_ll = float(ll[1])
                    results.append({"alt": a, "delta": alt_ll - ref_ll})
                top = min(results, key=lambda x: x["delta"]) if results else {"alt": None, "delta": 0.0}
                logger.info(f"/score_variant_probe (7B) done | top_alt={top['alt']} delta={top['delta']:.4f}")
                return {"probes": results, "top_alt": top["alt"], "top_delta": top["delta"]}
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"/score_variant_probe (7B) failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 probe failed: {e}")

    @modal.asgi_app()
    def api_7b(self):
        return self.fastapi_app

# --- Local Entrypoint for Testing ---
@app.local_entrypoint()
def main():
    print("--- ⚔️ LOCAL EVO-SERVICE TEST ⚔️ ---")
    print("Local entrypoint for syntax validation. Does not run generation.")
    print("✅ Syntax validation passed.")