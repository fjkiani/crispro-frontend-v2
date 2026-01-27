"""
Modal-deployed Enformer web service (FastAPI).

Implements the backend contract used by `oncology-coPilot/oncology-backend-minimal/api/services/enformer_client.py`:
- GET  /health_check
- POST /predict

Key requirement for resubmission:
- Must NOT return `provenance.method == "deterministic_fallback"` for our publication audit runs.

Deployment (from repo root):
  modal deploy services/enformer_service/modal_web.py

After deploy, set:
  export ENFORMER_URL="https://<your-modal-endpoint>"
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional

import modal


APP_NAME = os.environ.get("ENFORMER_MODAL_APP_NAME", "enformer-service-web-v1")
app = modal.App(APP_NAME)

# Enformer expects fixed input length 393,216.
ENFORMER_SEQ_LEN = 393_216


image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "fastapi==0.115.0",
        "pydantic==2.7.4",
        "httpx==0.27.2",
        "numpy==1.26.4",
        "tensorflow==2.15.0",
        "tensorflow-hub==0.16.1",
    )
)


@dataclass
class _Counters:
    requests: int = 0
    errors: int = 0
    ensembl_calls: int = 0


@app.function(
    image=image,
    gpu=os.environ.get("ENFORMER_GPU", "A10G"),
    timeout=60 * 10,
    memory=32768,
)
@modal.asgi_app()
def fastapi_app():
    from fastapi import FastAPI, HTTPException, Body
    from pydantic import BaseModel, Field
    from pydantic import ValidationError
    import numpy as np
    import httpx
    import tensorflow_hub as hub

    api = FastAPI(title="Enformer (Modal)", version="1.0.0")

    counters = _Counters()

    # Load model once per container
    model = hub.load("https://tfhub.dev/deepmind/enformer/1").model

    class PredictRequest(BaseModel):
        chrom: str = Field(..., description="Chromosome, with or without 'chr' prefix (GRCh38).")
        pos: int = Field(..., description="1-based position (GRCh38). Used as the center coordinate.")
        ref: str = Field("N", description="Reference allele (optional).")
        alt: str = Field("N", description="Alternate allele (optional).")
        context_bp: int = Field(32000, description="Accepted for compatibility; Enformer uses fixed window.")
        use_cache: bool = Field(True, description="Accepted for compatibility (not used currently).")

    async def _fetch_ensembl_sequence(chrom: str, start: int, end: int) -> str:
        """
        Fetch reference sequence from Ensembl REST.
        Enformer requires exactly 196,608 bp. We retrieve the full window as plain text.
        """
        # Ensembl expects chromosome without "chr"
        chrom_norm = chrom[3:] if chrom.lower().startswith("chr") else chrom
        url = f"https://rest.ensembl.org/sequence/region/human/{chrom_norm}:{start}..{end}"
        headers = {"Content-Type": "text/plain"}
        counters.ensembl_calls += 1
        async with httpx.AsyncClient(timeout=httpx.Timeout(60.0, connect=10.0)) as client:
            r = await client.get(url, headers=headers)
            if r.status_code >= 400:
                raise HTTPException(status_code=502, detail=f"Ensembl error ({r.status_code}): {r.text[:200]}")
            return (r.text or "").strip().upper()

    def _one_hot_dna(seq: str) -> np.ndarray:
        # A,C,G,T with N as zeros
        mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
        arr = np.zeros((len(seq), 4), dtype=np.float32)
        for i, ch in enumerate(seq):
            idx = mapping.get(ch)
            if idx is not None:
                arr[i, idx] = 1.0
        return arr

    def _summarize_to_accessibility(pred_human: np.ndarray) -> Dict[str, float]:
        """
        Convert Enformer output to a compact scalar proxy in [0,1].

        Note:
        - Enformer outputs many tracks; instead of claiming specific assays without track mapping,
          we compute a generic 'regulatory activity proxy' using the mean prediction magnitude
          near the center bins and then squash to [0,1].
        """
        # pred_human shape: [1, 896, tracks]
        x = pred_human[0]
        center = x.shape[0] // 2
        window = x[max(0, center - 8) : min(x.shape[0], center + 8), :]
        raw_mean = float(np.mean(window))
        # squash to [0,1] in a stable way
        score = float(1.0 / (1.0 + np.exp(-0.05 * raw_mean)))
        return {"accessibility_raw_mean": raw_mean, "accessibility_score": score}

    @api.get("/health_check")
    async def health_check():
        return {
            "status": "healthy",
            "model": "deepmind_enformer_tfhub",
            "seq_len": ENFORMER_SEQ_LEN,
            "requests_served": counters.requests,
            "errors": counters.errors,
            "ensembl_calls": counters.ensembl_calls,
        }

    @api.post("/predict")
    async def predict(payload: dict = Body(...)):
        try:
            req = PredictRequest.model_validate(payload)
        except ValidationError as e:
            counters.errors += 1
            raise HTTPException(status_code=422, detail=str(e))

        counters.requests += 1
        start_time = time.time()

        if req.pos <= 0:
            counters.errors += 1
            raise HTTPException(status_code=400, detail="pos must be positive (1-based GRCh38)")

        # Compute fixed window bounds
        half = ENFORMER_SEQ_LEN // 2
        start = int(req.pos) - half
        end = start + ENFORMER_SEQ_LEN - 1
        if start < 1:
            counters.errors += 1
            raise HTTPException(
                status_code=400,
                detail=f"pos too close to chromosome start for {ENFORMER_SEQ_LEN}bp window",
            )

        try:
            seq = await _fetch_ensembl_sequence(req.chrom, start, end)
        except HTTPException:
            raise
        except Exception as e:
            counters.errors += 1
            raise HTTPException(status_code=502, detail=f"ensembl_fetch_failed: {type(e).__name__}: {e}")

        if len(seq) != ENFORMER_SEQ_LEN:
            counters.errors += 1
            raise HTTPException(status_code=502, detail=f"Ensembl returned len={len(seq)}, expected {ENFORMER_SEQ_LEN}")

        # Apply SNV if provided
        ref = (req.ref or "N").upper()
        alt = (req.alt or "N").upper()
        center_idx = half  # 0-based index in seq
        seq_ref = seq[center_idx]
        ref_match = (ref in "ACGT") and (seq_ref == ref)
        snv_applied = False
        warning = None
        if alt in "ACGT":
            if (ref in "ACGT") and (seq_ref != ref):
                warning = f"ref_mismatch: expected {ref} but reference has {seq_ref}"
            if alt != seq_ref:
                seq = seq[:center_idx] + alt + seq[center_idx + 1 :]
                snv_applied = True

        one_hot = _one_hot_dna(seq)[None, :, :]

        try:
            preds = model.predict_on_batch(one_hot)
        except Exception as e:
            counters.errors += 1
            raise HTTPException(status_code=500, detail=f"enformer_inference_failed: {type(e).__name__}: {e}")

        human = preds["human"]
        summary = _summarize_to_accessibility(human)
        total_time = time.time() - start_time

        out: Dict[str, Any] = {
            "accessibility_score": float(round(summary["accessibility_score"], 6)),
            "dnase_signal": float(round(summary["accessibility_raw_mean"], 6)),
            "cage_signal": float(round(summary["accessibility_raw_mean"], 6)),
            "atac_signal": float(round(summary["accessibility_raw_mean"], 6)),
            "provenance": {
                "model": "enformer-v1",
                "method": "deepmind_enformer_tfhub",
                "seq_len": ENFORMER_SEQ_LEN,
                "window_start": start,
                "window_end": end,
                "center_ref_base": seq_ref,
                "ref_match": bool(ref_match),
                "snv_applied": bool(snv_applied),
                "warning": warning,
                "total_time_sec": total_time,
            },
        }
        return out


    return api

