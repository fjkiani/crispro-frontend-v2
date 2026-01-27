"""FastAPI endpoint handlers for Evo2 service."""
import uuid
import httpx
from fastapi import FastAPI, HTTPException
from loguru import logger

from .models import JobSubmitResponse, GenerationRequest


def register_endpoints(app: FastAPI, service_instance, job_status_dict):
    """
    Register all FastAPI endpoints with the FastAPI app.
    
    Args:
        app: FastAPI application instance
        service_instance: EvoService instance (for accessing self.model)
        job_status_dict: Modal Dict for job status tracking
    """
    
    @app.post("/generate", status_code=202, response_model=JobSubmitResponse)
    def submit_generation(request: GenerationRequest):
            job_id = str(uuid.uuid4())
            logger.info(f"Received submission request. Assigning job ID: {job_id}")
            job_status_dict[job_id] = {"status": "pending"}
            service_instance.generate_in_background.spawn(job_id, request.prompt, request.n_tokens)
            return {"job_id": job_id}

    @app.get("/status/{job_id}")
    def get_status(job_id: str):
            logger.info(f"Received status query for job ID: {job_id}")
            if job_id not in job_status_dict:
                raise HTTPException(status_code=404, detail="Job ID not found.")
            return job_status_dict[job_id]

    @app.post("/score_delta")
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
                ll = service_instance.model.score_sequences([ref_sequence, alt_sequence])
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

    @app.post("/score_batch")
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
                ll = service_instance.model.score_sequences(seqs)
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

    @app.post("/score_variant")
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
                ll = service_instance.model.score_sequences([ref_sequence, alt_sequence])
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

    @app.post("/score_variant_multi")
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
                        ll = service_instance.model.score_sequences([ref_sequence, alt_sequence])
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

    @app.post("/score_variant_exon")
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
                ll = service_instance.model.score_sequences([ref_sequence, alt_sequence])
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

    @app.post("/score_variant_profile")
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
                    ll = service_instance.model.score_sequences([ref_seq, alt_seq])
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

    @app.post("/score_variant_probe")
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
                    ll = service_instance.model.score_sequences([ref_seq, alt_seq])
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

    @app.post("/score_variant_with_activations")
    def score_variant_with_activations(item: dict):
            """
            Score a variant AND return layer 26 activations for SAE extraction.
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
            logger.info(f"/score_variant_with_activations called | {chrom}:{pos} {ref}>{alt} window={window} return_activations={return_activations}")
            
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
            
            try:
                # Tokenize sequences
                ref_input_ids = torch.tensor([service_instance.model.tokenizer.encode(ref_sequence)])
                alt_input_ids = torch.tensor([service_instance.model.tokenizer.encode(alt_sequence)])
                
                # Move to GPU if available
                device = next(service_instance.model.model.parameters()).device
                ref_input_ids = ref_input_ids.to(device)
                alt_input_ids = alt_input_ids.to(device)
                
                if return_activations:
                    # Forward pass with activation extraction for both sequences
                    ref_logits, ref_embeddings = service_instance.model.forward(
                        ref_input_ids,
                        return_embeddings=True,
                        layer_names=["blocks.26"]
                    )
                    alt_logits, alt_embeddings = service_instance.model.forward(
                        alt_input_ids,
                        return_embeddings=True,
                        layer_names=["blocks.26"]
                    )
                    
                    # Extract layer 26 activations
                    ref_layer26 = ref_embeddings.get("blocks.26")
                    alt_layer26 = alt_embeddings.get("blocks.26")
                    
                    # Compute log-likelihoods from logits
                    ref_ll = service_instance.model.score_sequences([ref_sequence])[0]
                    alt_ll = service_instance.model.score_sequences([alt_sequence])[0]
                    
                    logger.info(f"/score_variant_with_activations done | delta={alt_ll - ref_ll:.4f} | activations extracted")
                    
                    return {
                        "ref_likelihood": float(ref_ll),
                        "alt_likelihood": float(alt_ll),
                        "delta_score": float(alt_ll - ref_ll),
                        "window_start": start,
                        "window_end": end,
                        "variant_index": idx,
                        "layer_26_activations": {
                            "ref": ref_layer26.cpu().numpy().tolist() if ref_layer26 is not None else None,
                            "alt": alt_layer26.cpu().numpy().tolist() if alt_layer26 is not None else None,
                            "shape": list(ref_layer26.shape) if ref_layer26 is not None else None
                        }
                    }
                else:
                    # Standard scoring without activations
                    ll = service_instance.model.score_sequences([ref_sequence, alt_sequence])
                    ref_ll = float(ll[0])
                    alt_ll = float(ll[1])
                    delta = alt_ll - ref_ll
                    logger.info(f"/score_variant_with_activations done | delta={delta:.4f} | no activations")
                    return {
                        "ref_likelihood": ref_ll,
                        "alt_likelihood": alt_ll,
                        "delta_score": delta,
                        "window_start": start,
                        "window_end": end,
                        "variant_index": idx,
                    }
            except Exception as e:
                logger.error(f"Scoring/activation extraction failed: {e}")
                raise HTTPException(status_code=500, detail=f"Evo2 scoring/activation failed: {e}")

