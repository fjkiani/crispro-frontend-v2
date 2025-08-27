#!/usr/bin/env python3
import argparse
import csv
import json
from typing import List, Dict
import hashlib
import os

import httpx
import time
from sklearn.metrics import precision_recall_fscore_support


"""
Guidance Tier Evaluation

Input CSV columns (header required):
  disease,therapy,drug_or_class,gene,hgvs_p,chrom,pos,ref,alt,build,on_label_truth

For each row, calls /api/guidance/chemo and compares:
  predicted Tier I (Yes) vs on_label_truth (Yes if on-label/guideline) as a proxy.

Outputs JSON with precision/recall/F1 and per-row results.
"""


def parse_bool(s: str) -> bool:
    return str(s).strip().lower() in {"1", "true", "yes", "y"}


def main():
    ap = argparse.ArgumentParser(description="Evaluate Guidance Tier I vs on-label truth")
    ap.add_argument("--csv", required=True, help="Input CSV with columns incl. on_label_truth")
    ap.add_argument("--api_base", default="http://127.0.0.1:8000")
    ap.add_argument("--out", default="guidance_tier_eval_results.json")
    ap.add_argument("--cache_path", default="tools/benchmarks/.guidance_eval_cache.json", help="On-disk cache file path")
    ap.add_argument("--no_cache", action="store_true", help="Disable on-disk caching")
    ap.add_argument(
        "--mode",
        choices=["chemo", "radonc", "synthetic"],
        default="chemo",
        help="Which guidance endpoint to call (default: chemo)",
    )
    args = ap.parse_args()

    rows: List[Dict] = []
    with open(args.csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    y_true = []
    y_pred = []
    detailed = []

    # Optional on-disk cache to avoid repeat API calls
    cache_path = args.cache_path
    use_cache = not args.no_cache
    cache: Dict[str, Dict] = {}
    if use_cache and os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as cf:
                cache = json.load(cf) or {}
        except Exception:
            cache = {}

    with httpx.Client(timeout=180.0) as client:
        for r in rows:
            disease = r.get("disease") or ""
            therapy = r.get("therapy") or r.get("drug_or_class") or ""
            gene = r.get("gene") or ""
            hgvs_p = r.get("hgvs_p") or ""
            chrom = r.get("chrom") or None
            pos = r.get("pos") or None
            ref = r.get("ref") or None
            alt = r.get("alt") or None
            build = r.get("build") or None
            truth = parse_bool(r.get("on_label_truth", "false"))

            payload = {"mutations": [], "api_base": args.api_base}
            # Select endpoint/mode. Route platinum/PARP to synthetic if using chemo mode.
            if args.mode == "synthetic" or therapy.lower() in {"platinum", "parp inhibitor"}:
                endpoint = f"{args.api_base}/api/guidance/synthetic_lethality"
                payload.update({"disease": disease})
            elif args.mode == "radonc":
                endpoint = f"{args.api_base}/api/guidance/radonc"
                payload.update({"disease": disease, "options": {"adaptive": True, "ensemble": True}})
            else:
                endpoint = f"{args.api_base}/api/guidance/chemo"
                payload.update({
                    "disease": disease,
                    "drug_or_class": therapy,
                    "options": {"adaptive": True, "ensemble": True},
                })
            if gene and hgvs_p:
                mut = {"gene": gene, "hgvs_p": hgvs_p}
                if chrom and pos and ref and alt:
                    try:
                        mut.update({"chrom": str(chrom), "pos": int(pos), "ref": str(ref), "alt": str(alt)})
                    except Exception:
                        pass
                if build:
                    mut["build"] = build
                payload["mutations"].append(mut)

            # Retry/backoff loop with cache
            pred_tier_i = False
            details = {}
            last_err = None
            js = None

            # Build cache key from endpoint + payload
            cache_key_src = json.dumps({"endpoint": endpoint, "payload": payload}, sort_keys=True)
            cache_key = hashlib.sha256(cache_key_src.encode("utf-8")).hexdigest()

            if use_cache and cache_key in cache:
                js = cache.get(cache_key)
            else:
                for attempt in range(3):
                    try:
                        resp = client.post(endpoint, json=payload)
                        if resp.status_code < 400:
                            js = resp.json() or {}
                            if use_cache:
                                cache[cache_key] = js
                            break
                        else:
                            last_err = resp.text
                    except Exception as e:
                        last_err = str(e)
                    time.sleep(0.8)

            if js is None:
                js = {}

            # Parse tier; handle nested synthetic_lethality shape where tier is under guidance
            tier = js.get("tier")
            if tier is None and isinstance(js.get("guidance"), dict):
                tier = (js.get("guidance") or {}).get("tier")
                # Prefer nested details when available
                if isinstance(js.get("guidance"), dict):
                    details = js.get("guidance")
            if not details:
                details = js if js else ({"error": last_err} if last_err else {})

            pred_tier_i = (tier == "I")

            y_true.append(1 if truth else 0)
            y_pred.append(1 if pred_tier_i else 0)
            detailed.append({"input": r, "prediction": details})

    precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average="binary", zero_division=0)

    out = {
        "metrics": {"precision": precision, "recall": recall, "f1": f1},
        "n": len(rows),
        "api_base": args.api_base,
        "results": detailed,
    }
    # Persist cache if enabled
    if use_cache:
        try:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            with open(cache_path, "w", encoding="utf-8") as cf:
                json.dump(cache, cf)
        except Exception:
            pass

    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()


