#!/usr/bin/env python3
import json
import time
import argparse
import httpx
from sklearn.metrics import roc_auc_score, average_precision_score

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--micro_json", required=True, help="Path to AM-covered micro JSON {variants:[]} from Modal sampler")
    ap.add_argument("--api_base", default="http://127.0.0.1:8000")
    ap.add_argument("--fusion_url", required=True, help="Fusion-engine URL for AlphaMissense lookups")
    args = ap.parse_args()

    with open(args.micro_json, "r") as f:
        js = json.load(f)
    variants = js.get("variants", [])

    y, evo_s, am_s, fused_s = [], [], [], []

    with httpx.Client(timeout=120.0) as client:
        for v in variants:
            chrom = str(v["chrom"]).replace("chr", "")
            pos = int(v["pos"])
            ref = str(v["ref"]).upper()
            alt = str(v["alt"]).upper()
            # Evo score with retries
            evo = 0.0
            for attempt in range(3):
                try:
                    r = client.post(f"{args.api_base}/api/evo/score_variant_multi", json={
                        "assembly": "GRCh38", "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "model_id": "evo2_7b"
                    })
                    if r.status_code < 400:
                        evo = abs(float((r.json() or {}).get("min_delta") or 0.0))
                        break
                except Exception:
                    pass
                time.sleep(1.5 * (attempt + 1))

            # AM score via fusion-engine
            am = 0.0
            try:
                payload = {
                    "protein_sequence": "M" * 64,
                    "variants": [{"variant_id": "v1", "hgvs": "NA", "alphamissense_variant_str": f"{v['chrom']}:{pos}:{ref}:{alt}"}],
                }
                ra = client.post(f"{args.fusion_url.rstrip('/')}/score_variants", json=payload)
                if ra.status_code < 400:
                    arr = (ra.json() or {}).get("scored_variants") or []
                    if arr and isinstance(arr[0].get("alphamissense_score", None), (int, float)):
                        am = float(arr[0]["alphamissense_score"])
            except Exception:
                pass

            fused = max(evo, am)
            y.append(int(v["label"]))
            evo_s.append(evo)
            am_s.append(am)
            fused_s.append(fused)

    def pr(name, s):
        if len(set(y)) >= 2:
            print(json.dumps({name: {"auroc": roc_auc_score(y, s), "auprc": average_precision_score(y, s), "n": len(y)}}, indent=2))
        else:
            print(json.dumps({name: {"error": "insufficient class variety", "n": len(y)}}, indent=2))

    pr("crisrPRO_only", evo_s)
    pr("am_only", am_s)
    pr("fused_max", fused_s)

if __name__ == "__main__":
    main()


