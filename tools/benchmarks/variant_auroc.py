#!/usr/bin/env python3
import argparse
import csv
import gzip
import io
import json
import os
import random
import sys
from dataclasses import dataclass
from typing import List, Tuple

import httpx
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import train_test_split, GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
import numpy as np
import os


CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


def download_clinvar_variant_summary(dest_path: str) -> None:
    with httpx.stream("GET", CLINVAR_URL, timeout=None) as r:
        r.raise_for_status()
        with open(dest_path, "wb") as f:
            for chunk in r.iter_bytes():
                f.write(chunk)


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    label: int  # 1=pathogenic, 0=benign
    gene: str = ""  # GeneSymbol if available


def load_variants_from_variant_summary(
    path: str,
    n_pos: int,
    n_neg: int,
    allow_grch37: bool = False,
    include_likely: bool = True,
) -> Tuple[List[Variant], List[Variant]]:
    def is_snv(ref: str, alt: str) -> bool:
        return len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and ref != alt

    # ClinVar variant_summary.txt header fields
    # We need: Chromosome, Start, ReferenceAllele, AlternateAllele, ClinicalSignificance, Assembly
    opener = gzip.open if path.endswith(".gz") else open
    positives: List[Variant] = []
    negatives: List[Variant] = []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            asm = row.get("Assembly", "")
            if asm not in ("GRCh38", "GRCh37"):
                continue
            if asm == "GRCh37" and not allow_grch37:
                continue
            chrom = row.get("Chromosome") or ""
            start = row.get("PositionVCF") or ""  # Use VCF position
            ref = (row.get("ReferenceAlleleVCF") or "").upper()  # Use VCF alleles
            alt = (row.get("AlternateAlleleVCF") or "").upper()
            gene_symbol = (row.get("GeneSymbol") or "").strip()
            csig = (row.get("ClinicalSignificance") or "").lower()
            # basic chromosome filter
            if not chrom or not start.isdigit():
                continue
            if not is_snv(ref, alt):
                continue
            pos = int(start)

            # Labeling: positives = (pathogenic[/likely]) and not conflicting
            # negatives = (benign[/likely]) and not conflicting
            if "conflict" in csig:
                continue
            cs = csig.replace(' ', '_')
            pathogenic_terms = ["pathogenic", "likely_pathogenic"] if include_likely else ["pathogenic"]
            benign_terms = ["benign", "likely_benign"] if include_likely else ["benign"]
            is_pos = any(t in cs for t in pathogenic_terms)
            is_neg = any(t in cs for t in benign_terms) and (not is_pos)

            if is_pos and len(positives) < n_pos:
                positives.append(Variant(chrom=str(chrom), pos=pos, ref=ref, alt=alt, label=1, gene=gene_symbol))
            elif is_neg and len(negatives) < n_neg:
                negatives.append(Variant(chrom=str(chrom), pos=pos, ref=ref, alt=alt, label=0, gene=gene_symbol))

            if len(positives) >= n_pos and len(negatives) >= n_neg:
                break

    return positives, negatives


async def _post_score(client: httpx.AsyncClient, url: str, payload: dict) -> httpx.Response:
    return await client.post(url, json=payload)

async def fetch_min_delta_bidir(client: httpx.AsyncClient, api_base: str, v: Variant, model_id: str) -> float:
    import asyncio
    ref = (v.ref or "").upper()
    alt = (v.alt or "").upper()
    payload = {"assembly": "GRCh38", "chrom": str(v.chrom), "pos": int(v.pos), "ref": ref, "alt": alt, "model_id": model_id}
    for attempt in range(3):
        try:
            r1 = await _post_score(client, f"{api_base}/api/evo/score_variant_multi", payload)
            if r1.status_code == 400 and "mismatch" in (r1.text or "").lower():
                return -999.0
            md1 = float((r1.json() or {}).get("min_delta") or 0.0) if r1.status_code < 400 else 0.0
            # reverse
            payload_rev = dict(payload)
            payload_rev.update({"ref": alt, "alt": ref})
            r2 = await _post_score(client, f"{api_base}/api/evo/score_variant_multi", payload_rev)
            if r2.status_code == 400 and "mismatch" in (r2.text or "").lower():
                return -999.0
            md2 = float((r2.json() or {}).get("min_delta") or 0.0) if r2.status_code < 400 else 0.0
            return float((abs(md1) + abs(md2)) / 2.0)
        except Exception:
            await asyncio.sleep(0.3 * (attempt + 1))
    return 0.0

async def fetch_best_exon_delta_bidir(client: httpx.AsyncClient, api_base: str, v: Variant, model_id: str, flanks: List[int]) -> float:
    import asyncio
    ref = (v.ref or "").upper()
    alt = (v.alt or "").upper()
    best_abs = 0.0
    for flank in flanks:
        payload = {"assembly": "GRCh38", "chrom": str(v.chrom), "pos": int(v.pos), "ref": ref, "alt": alt, "flank": int(flank), "model_id": model_id}
        try:
            r1 = await _post_score(client, f"{api_base}/api/evo/score_variant_exon", payload)
            if r1.status_code == 400 and "mismatch" in (r1.text or "").lower():
                return -999.0
            ed1 = float((r1.json() or {}).get("exon_delta") or 0.0) if r1.status_code < 400 else 0.0
            payload_rev = dict(payload)
            payload_rev.update({"ref": alt, "alt": ref})
            r2 = await _post_score(client, f"{api_base}/api/evo/score_variant_exon", payload_rev)
            if r2.status_code == 400 and "mismatch" in (r2.text or "").lower():
                return -999.0
            ed2 = float((r2.json() or {}).get("exon_delta") or 0.0) if r2.status_code < 400 else 0.0
            cur = max(abs(ed1), abs(ed2))
            if cur > best_abs:
                best_abs = cur
        except Exception:
            await asyncio.sleep(0.2)
            continue
    return best_abs

async def fetch_features(client: httpx.AsyncClient, api_base: str, v: Variant) -> dict:
    """Extract richer features for supervised head using 7B (+optional 40B)."""
    md7 = await fetch_min_delta_bidir(client, api_base, v, "evo2_7b")
    if md7 == -999.0:
        return {"skip": True}
    # Light mode via env to keep small runs fast
    light = os.getenv("FUSION_LIGHT", "1") == "1"
    flanks7 = [4096] if light else [4096, 8192, 12500, 16384]
    ex7 = await fetch_best_exon_delta_bidir(client, api_base, v, "evo2_7b", flanks=flanks7)
    if ex7 == -999.0:
        return {"skip": True}
    # Optional 40B
    md40, ex40 = None, None
    if os.getenv("FUSION_USE_40B", "0") == "1" and not light:
        try:
            md40 = await fetch_min_delta_bidir(client, api_base, v, "evo2_40b")
            if md40 == -999.0:
                md40 = None
            ex40 = await fetch_best_exon_delta_bidir(client, api_base, v, "evo2_40b", flanks=[4096, 8192])
            if ex40 == -999.0:
                ex40 = None
        except Exception:
            md40, ex40 = None, None

    # AlphaMissense lookup (optional)
    am_score = None
    try:
        # Lazy import to avoid hard dependency
        try:
            from tools.alphamissense_client import AlphaMissenseClient  # type: ignore
        except Exception:
            from src.tools.alphamissense_client import AlphaMissenseClient  # type: ignore
        global _AM_CLIENT
        if '_AM_CLIENT' not in globals() or _AM_CLIENT is None:
            am_path = os.getenv("ALPHAMISSENSE_DATA", "/data/AlphaMissense_hg38.parquet")
            _AM_CLIENT = AlphaMissenseClient(data_path=am_path)
        chrom_str = str(v.chrom)
        if not chrom_str.startswith("chr"):
            chrom_str = f"chr{chrom_str}"
        lookup_key = f"{chrom_str}:{int(v.pos)}:{str(v.ref)}:{str(v.alt)}"
        am = _AM_CLIENT.get_score(lookup_key)
        am_score = float(am.get("score")) if am and isinstance(am.get("score", None), (int, float)) else None
    except Exception:
        am_score = None
    # HTTP fallback via Modal fusion-engine if no local parquet
    if am_score is None:
        try:
            fusion_url = os.getenv("FUSION_AM_URL", "").rstrip("/")
            if fusion_url:
                payload = {
                    "protein_sequence": "M" * 128,
                    "variants": [
                        {
                            "variant_id": "v1",
                            "hgvs": "NA",
                            "alphamissense_variant_str": lookup_key,
                        }
                    ],
                }
                r = await client.post(f"{fusion_url}/score_variants", json=payload)
                if r.status_code < 400:
                    js = r.json() or {}
                    arr = js.get("scored_variants") or []
                    if arr:
                        am_score = float(arr[0].get("alphamissense_score"))
        except Exception:
            pass
    return {
        "skip": False,
        "min7": float(md7),
        "exon7": float(ex7),
        "min40": float(md40) if md40 is not None else None,
        "exon40": float(ex40) if ex40 is not None else None,
        "am": float(am_score) if am_score is not None else None,
    }

def per_gene_z_normalize(values: List[float], genes: List[str]) -> List[float]:
    arr = np.asarray(values, dtype=float)
    genes_arr = np.asarray(genes)
    out = np.copy(arr)
    for g in np.unique(genes_arr):
        idx = np.where(genes_arr == g)[0]
        if idx.size >= 2:
            mu = float(np.mean(arr[idx]))
            sd = float(np.std(arr[idx]) + 1e-9)
            out[idx] = (arr[idx] - mu) / sd
    return out.tolist()

async def build_supervised_dataset(api_base: str, variants: List[Variant], require_am: bool = False) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    X_rows: List[List[float]] = []
    y: List[int] = []
    genes: List[str] = []
    async with httpx.AsyncClient(timeout=60.0) as client:
        for v in variants:
            feats = await fetch_features(client, api_base, v)
            if feats.get("skip"):
                continue
            min7 = feats.get("min7") or 0.0
            exon7 = feats.get("exon7") or 0.0
            min40 = feats.get("min40") if feats.get("min40") is not None else 0.0
            exon40 = feats.get("exon40") if feats.get("exon40") is not None else 0.0
            am_val = feats.get("am", None)
            if require_am and (am_val is None):
                # Skip rows without AlphaMissense coverage when required
                continue
            am = am_val if am_val is not None else 0.0
            comb7 = max(min7, exon7)
            comb40 = max(min40, exon40)
            avg_comb = comb7 if comb40 == 0.0 else (0.5 * (comb7 + comb40))
            # Features: 7B/40B combos + AlphaMissense
            X_rows.append([min7, exon7, min40, exon40, comb7, comb40, avg_comb, am])
            y.append(v.label)
            genes.append(v.gene or str(v.chrom))
    if not X_rows:
        return np.zeros((0, 8), dtype=float), np.zeros((0,), dtype=int), []
    X = np.asarray(X_rows, dtype=float)
    z1 = per_gene_z_normalize(X[:, 0].tolist(), genes)
    z2 = per_gene_z_normalize(X[:, 1].tolist(), genes)
    X = np.concatenate([X, np.asarray(z1)[:, None], np.asarray(z2)[:, None]], axis=1)
    return X, np.asarray(y, dtype=int), genes

def run_group_kfold_supervised(X: np.ndarray, y: np.ndarray, groups: List[str]) -> Tuple[float, float]:
    if X.shape[0] < 20 or len(set(y.tolist())) < 2:
        return None, None
    gkf = GroupKFold(n_splits=min(5, max(2, len(set(groups)))))
    y_true_all: List[int] = []
    y_score_all: List[float] = []
    for train_idx, test_idx in gkf.split(X, y, groups=groups):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        base = LogisticRegression(max_iter=300, solver="liblinear")
        # Calibrate with internal CV on train
        clf = CalibratedClassifierCV(base, method="isotonic", cv=3)
        clf.fit(X_train, y_train)
        proba = clf.predict_proba(X_test)[:, 1]
        y_true_all.extend(y_test.tolist())
        y_score_all.extend(proba.tolist())
    try:
        auroc = roc_auc_score(y_true_all, y_score_all)
        auprc = average_precision_score(y_true_all, y_score_all)
        return auroc, auprc
    except Exception:
        return None, None

# --------- AM-only helper for fused baseline ---------
async def fetch_am_score_only(http_client: httpx.AsyncClient, v: Variant) -> float | None:
    # Try local parquet first
    am_score: float | None = None
    try:
        try:
            from tools.alphamissense_client import AlphaMissenseClient  # type: ignore
        except Exception:
            from src.tools.alphamissense_client import AlphaMissenseClient  # type: ignore
        global _AM_CLIENT_ONLY
        if '_AM_CLIENT_ONLY' not in globals() or _AM_CLIENT_ONLY is None:
            am_path = os.getenv("ALPHAMISSENSE_DATA", "/data/AlphaMissense_hg38.parquet")
            _AM_CLIENT_ONLY = AlphaMissenseClient(data_path=am_path)
        chrom_str = str(v.chrom)
        if not chrom_str.startswith("chr"):
            chrom_str = f"chr{chrom_str}"
        key = f"{chrom_str}:{int(v.pos)}:{str(v.ref)}:{str(v.alt)}"
        am = _AM_CLIENT_ONLY.get_score(key)
        if am and isinstance(am.get("score", None), (int, float)):
            am_score = float(am["score"])
    except Exception:
        am_score = None

    if am_score is None:
        try:
            fusion_url = os.getenv("FUSION_AM_URL", "").rstrip("/")
            if fusion_url:
                chrom_str = str(v.chrom)
                if not chrom_str.startswith("chr"):
                    chrom_str = f"chr{chrom_str}"
                key = f"{chrom_str}:{int(v.pos)}:{str(v.ref)}:{str(v.alt)}"
                payload = {
                    "protein_sequence": "M" * 64,
                    "variants": [{"variant_id": "v1", "hgvs": "NA", "alphamissense_variant_str": key}],
                }
                r = await http_client.post(f"{fusion_url}/score_variants", json=payload)
                if r.status_code < 400:
                    js = r.json() or {}
                    arr = js.get("scored_variants") or []
                    if arr and isinstance(arr[0].get("alphamissense_score", None), (int, float)):
                        am_score = float(arr[0]["alphamissense_score"])
        except Exception:
            am_score = None
    return am_score

async def score_variants(api_base: str, variants: List[Variant], model_id: str, fast: bool = False) -> Tuple[List[Variant], List[float]]:
    # Simple sequential scoring to avoid event loop complexity for smoke tests
    valid_variants: List[Variant] = []
    valid_scores: List[float] = []
    async with httpx.AsyncClient(timeout=60.0) as client:
        for v in variants:
            min_delta_abs = await fetch_min_delta_bidir(client, api_base, v, model_id)
            if min_delta_abs == -999.0:
                # Skip ref/alt mismatches
                continue
            if fast:
                valid_variants.append(v)
                valid_scores.append(float(min_delta_abs))
                continue
            # Probe a small set of exon windows
            best_exon_abs = await fetch_best_exon_delta_bidir(client, api_base, v, model_id, flanks=[2048, 4096])
            if best_exon_abs == -999.0:
                continue
            # Combine: take the stronger signal
            combined = max(float(min_delta_abs), float(best_exon_abs))
            valid_variants.append(v)
            valid_scores.append(combined)
    return valid_variants, valid_scores


def main():
    parser = argparse.ArgumentParser(description="Compute AUROC/AUPRC for Evo-based variant impact vs ClinVar")
    parser.add_argument("--clinvar_path", type=str, default="variant_summary.txt.gz", help="Path to ClinVar variant_summary.txt[.gz]")
    parser.add_argument("--download", action="store_true", help="Download ClinVar variant_summary to --clinvar_path if missing")
    parser.add_argument("--n_pos", type=int, default=100, help="Number of positive (pathogenic) SNVs")
    parser.add_argument("--n_neg", type=int, default=100, help="Number of negative (benign) SNVs")
    parser.add_argument("--api_base", type=str, default="http://127.0.0.1:8000", help="Backend API base URL")
    parser.add_argument("--model_id", type=str, default="evo2_7b", help="Evo2 model id")
    parser.add_argument("--out", type=str, default="variant_auroc_results.json", help="Where to write summary JSON")
    parser.add_argument("--allow_grch37", action="store_true", help="Include GRCh37 rows if GRCh38 is sparse")
    parser.add_argument("--no_likely", action="store_true", help="Exclude likely_pathogenic/likely_benign from labels")
    parser.add_argument("--fast", action="store_true", help="Fast mode: min_delta only; skip exon and supervised")
    parser.add_argument("--require_am", action="store_true", help="Require AlphaMissense coverage for supervised dataset")
    args = parser.parse_args()

    if not os.path.exists(args.clinvar_path):
        if args.download:
            print(f"Downloading ClinVar variant_summary to {args.clinvar_path} ...", file=sys.stderr)
            download_clinvar_variant_summary(args.clinvar_path)
        else:
            print(f"ClinVar file not found: {args.clinvar_path}. Pass --download to fetch it.", file=sys.stderr)
            sys.exit(1)

    print("Loading ClinVar and sampling variants...", file=sys.stderr)
    pos, neg = load_variants_from_variant_summary(
        args.clinvar_path,
        args.n_pos,
        args.n_neg,
        allow_grch37=args.allow_grch37,
        include_likely=(not args.no_likely),
    )
    if len(pos) < args.n_pos or len(neg) < args.n_neg:
        print(f"Warning: sampled fewer variants than requested (pos={len(pos)}, neg={len(neg)})", file=sys.stderr)

    variants = pos + neg
    initial_labels = [v.label for v in variants]

    if len(set(initial_labels)) < 2:
        print("Not enough class diversity in sampled variants; try increasing --n_pos/--n_neg or relaxing filters.", file=sys.stderr)
        sys.exit(1)

    import asyncio
    print("Scoring variants via Evo proxy...", file=sys.stderr)
    valid_variants, scores = asyncio.run(score_variants(args.api_base, variants, args.model_id, fast=args.fast))
    
    # Get labels for valid variants only
    labels = [v.label for v in valid_variants]
    
    if len(valid_variants) == 0:
        print("No valid variants after filtering reference mismatches.", file=sys.stderr)
        sys.exit(1)
    
    if len(set(labels)) < 2:
        print("Not enough class diversity in valid variants after filtering.", file=sys.stderr)
        sys.exit(1)

    try:
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
    except Exception as e:
        print(f"Failed to compute metrics: {e}", file=sys.stderr)
        sys.exit(1)

    # Supervised head with isotonic calibration (skip in fast mode)
    sup_auroc, sup_auprc = None, None
    if not args.fast:
        print("Building supervised dataset...", file=sys.stderr)
        import asyncio as _asyncio
        X_all, y_all, genes = _asyncio.run(build_supervised_dataset(args.api_base, valid_variants, require_am=args.require_am))
        sup_auroc, sup_auprc = run_group_kfold_supervised(X_all, y_all, genes)

    # Fused AM-only baseline on covered subset for quick comparison
    fused_baseline = None
    try:
        import asyncio as _asyncio
        async def _am_fused(http_client: httpx.AsyncClient):
            y_cov: List[int] = []
            s_cov: List[float] = []
            for v in valid_variants:
                am = await fetch_am_score_only(http_client, v)
                if am is None or am in (-998.0, -999.0):
                    continue
                # Simple fusion: max(|Evo2|, AM)
                idx = valid_variants.index(v)
                evo = float(scores[idx])
                fused = max(evo, float(am))
                y_cov.append(v.label)
                s_cov.append(fused)
            if len(set(y_cov)) >= 2 and len(y_cov) >= 10:
                return {
                    "auroc": roc_auc_score(y_cov, s_cov),
                    "auprc": average_precision_score(y_cov, s_cov),
                    "n": len(y_cov),
                }
            return None
        async def _runner():
            async with httpx.AsyncClient(timeout=30.0) as _hc:
                return await _am_fused(_hc)
        fused_baseline = _asyncio.run(_runner())
    except Exception:
        fused_baseline = None

    summary = {
        "n_pos_initial": len(pos),
        "n_neg_initial": len(neg),
        "n_pos_valid": sum(1 for v in valid_variants if v.label == 1),
        "n_neg_valid": sum(1 for v in valid_variants if v.label == 0),
        "n_filtered": len(variants) - len(valid_variants),
        "model_id": args.model_id,
        "api_base": args.api_base,
        "metrics": {
            "auroc": auroc,
            "auprc": auprc,
        },
        "metrics_supervised_groupcv": {
            "auroc": sup_auroc,
            "auprc": sup_auprc,
        },
        "metrics_fused_baseline": fused_baseline,
    }
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()


