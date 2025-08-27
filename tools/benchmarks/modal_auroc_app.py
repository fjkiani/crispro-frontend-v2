import gzip
import io
import json
import os
import random
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

import modal


app = modal.App("crispro-variant-auroc")

image = (
    modal.Image.debian_slim()
    .pip_install(
        "httpx==0.27.0",
        "pandas==2.2.2",
        "scikit-learn==1.5.1",
        "tqdm==4.66.4",
    )
)


CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    label: int  # 1=pathogenic, 0=benign


def _is_snv(ref: str, alt: str) -> bool:
    return (
        len(ref) == 1
        and len(alt) == 1
        and ref in "ACGT"
        and alt in "ACGT"
        and ref != alt
    )


async def _download(url: str) -> bytes:
    import httpx

    async with httpx.AsyncClient(timeout=None) as client:
        r = await client.get(url)
        r.raise_for_status()
        return r.content


def _iter_variant_summary_tsv_bytes(data: bytes):
    import csv

    with gzip.open(io.BytesIO(data), mode="rt", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row


def _sample_variants(
    rows, n_pos: int, n_neg: int, allow_grch37: bool, include_likely: bool
) -> Tuple[List[Variant], List[Variant]]:
    positives: List[Variant] = []
    negatives: List[Variant] = []
    pathogenic_terms = ["pathogenic", "likely_pathogenic"] if include_likely else ["pathogenic"]
    benign_terms = ["benign", "likely_benign"] if include_likely else ["benign"]

    for row in rows:
        asm = row.get("Assembly", "")
        if asm not in ("GRCh38", "GRCh37"):
            continue
        if asm == "GRCh37" and not allow_grch37:
            continue
        chrom = (row.get("Chromosome") or "").strip()
        start = (row.get("Start") or "").strip()
        ref = (row.get("ReferenceAllele") or "").upper()
        alt = (row.get("AlternateAllele") or "").upper()
        csig = (row.get("ClinicalSignificance") or "").lower()
        if not chrom or not start.isdigit():
            continue
        if not _is_snv(ref, alt):
            continue
        if "conflict" in csig:
            continue
        pos = int(start)
        cs = csig.replace(" ", "_")
        is_pos = any(t in cs for t in pathogenic_terms)
        is_neg = any(t in cs for t in benign_terms) and not is_pos

        if is_pos and len(positives) < n_pos:
            positives.append(Variant(chrom=str(chrom), pos=pos, ref=ref, alt=alt, label=1))
        elif is_neg and len(negatives) < n_neg:
            negatives.append(Variant(chrom=str(chrom), pos=pos, ref=ref, alt=alt, label=0))

        if len(positives) >= n_pos and len(negatives) >= n_neg:
            break

    return positives, negatives


async def _score_min_delta_batch(client, evo_url: str, model_id: str, batch: List[Variant]) -> List[float]:
    import httpx

    scores: List[float] = []
    for v in batch:
        payload = {
            "assembly": "GRCh38",
            "chrom": str(v.chrom),
            "pos": int(v.pos),
            "ref": v.ref,
            "alt": v.alt,
            "model_id": model_id,
        }
        # Retry/backoff
        score = 0.0
        last_err: str | None = None
        for attempt in range(3):
            try:
                r = await client.post(f"{evo_url}/score_variant_multi", json=payload)
                if r.status_code < 400:
                    js = r.json() or {}
                    md = js.get("min_delta")
                    score = abs(float(md)) if md is not None else 0.0
                    break
                last_err = r.text
            except Exception as e:  # noqa: BLE001
                last_err = str(e)
            await _async_sleep(0.5 * (attempt + 1))
        scores.append(score)
    return scores


async def _async_sleep(s: float) -> None:
    import asyncio

    await asyncio.sleep(s)


@app.function(image=image, timeout=60 * 20)
async def run_auroc(
    n_pos: int = 1000,
    n_neg: int = 1000,
    include_likely: bool = True,
    allow_grch37: bool = True,
    model_id: str = "evo2_7b",
    evo_url_7b: str | None = None,
    evo_url_40b: str | None = None,
    out_path: str = "/tmp/variant_auroc_results.json",
) -> Dict:
    """
    Compute AUROC/AUPRC for Evo2 min_delta vs ClinVar labels in Modal.

    Returns JSON with metrics and sampling details. Writes JSON to out_path in the container.
    """
    import httpx
    from sklearn.metrics import roc_auc_score, average_precision_score

    evo_url = evo_url_7b or os.environ.get("EVO_URL_7B") or os.environ.get("EVO_URL_40B") or evo_url_40b
    if not evo_url:
        raise RuntimeError("EVO_URL_7B or EVO_URL_40B must be provided via args or env")

    # Download ClinVar summary (cached in memory for simplicity)
    data = await _download(CLINVAR_URL)
    rows = list(_iter_variant_summary_tsv_bytes(data))
    positives, negatives = _sample_variants(rows, n_pos, n_neg, allow_grch37, include_likely)

    variants = positives + negatives
    labels = [v.label for v in variants]

    if len(set(labels)) < 2:
        raise RuntimeError(
            f"Not enough class diversity (pos={len(positives)}, neg={len(negatives)}). Increase sample sizes or relax filters."
        )

    # Score min_delta in chunks
    scores: List[float] = []
    chunk = 64
    async with httpx.AsyncClient(timeout=httpx.Timeout(120.0, connect=10.0)) as client:
        for i in range(0, len(variants), chunk):
            batch = variants[i : i + chunk]
            scores += await _score_min_delta_batch(client, evo_url, model_id, batch)

    # Compute metrics
    auroc = float(roc_auc_score(labels, scores))
    auprc = float(average_precision_score(labels, scores))

    summary = {
        "n_pos": len(positives),
        "n_neg": len(negatives),
        "model_id": model_id,
        "evo_url": evo_url,
        "include_likely": include_likely,
        "allow_grch37": allow_grch37,
        "metrics": {"auroc": auroc, "auprc": auprc},
    }

    # Persist JSON in container
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


@app.local_entrypoint()
def main():
    """Local runner for quick testing with your local env EVO URLs."""
    import asyncio

    res = asyncio.run(
        run_auroc.remote(
            n_pos=200,
            n_neg=200,
            include_likely=True,
            allow_grch37=True,
            model_id="evo2_7b",
        )
    )
    print(json.dumps(res, indent=2))


