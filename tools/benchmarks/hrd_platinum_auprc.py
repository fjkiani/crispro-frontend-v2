#!/usr/bin/env python3
import argparse
import csv
import json
from typing import List, Dict

import httpx
from sklearn.metrics import average_precision_score, roc_auc_score

"""
HRD/Platinum AUPRC Benchmark (skeleton)

Input CSV columns (header required):
  disease,gene,hgvs_p,chrom,pos,ref,alt,build,outcome_platinum (1=respond,0=non-respond)

For each row:
  - Compute damage+dependency via guidance/synthetic_lethality or direct insights
  - Map to platinum/PARP score proxy (here: reuse /api/guidance/synthetic_lethality and score=1 if suggested_therapy includes platinum, else 0.5 fallback)
  - Compare proxy scores vs outcome labels; report AUPRC/AUROC
"""


def parse_int(s: str) -> int:
    try:
        return int(str(s).strip())
    except Exception:
        return 0


def main():
    ap = argparse.ArgumentParser(description="Benchmark HRD/platinum AUPRC using synthetic lethality guidance")
    ap.add_argument("--csv", required=True, help="Input cohort CSV with outcomes")
    ap.add_argument("--api_base", default="http://127.0.0.1:8000")
    ap.add_argument("--out", default="hrd_platinum_auprc_results.json")
    ap.add_argument("--use_evo", action="store_true", help="Use Evo2 /api/evo/score_variant_multi for S-only scoring instead of guidance")
    args = ap.parse_args()

    rows: List[Dict] = []
    with open(args.csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    y_true = []
    y_score = []
    detailed = []

    def evo_score(client: httpx.Client, base: str, r: Dict) -> float:
        chrom = r.get("chrom")
        pos = r.get("pos")
        ref = r.get("ref")
        alt = r.get("alt")
        build = r.get("build") or "GRCh38"
        # Require complete coords
        if not (chrom and pos and ref and alt):
            return 0.5
        try:
            payload = {"assembly": build, "chrom": str(chrom), "pos": int(pos), "ref": str(ref), "alt": str(alt)}
            resp = client.post(f"{base}/api/evo/score_variant_multi", json=payload)
            if resp.status_code >= 400:
                return 0.5
            js = resp.json() or {}
            md = js.get("min_delta")
            if md is None:
                return 0.5
            # Map |min_delta| to a [0.5,0.9] proxy score (cheap S-only heuristic)
            mag = abs(float(md))
            lift = min(0.4, mag * 50.0)  # 0.02 -> +1.0; cap at +0.4
            return max(0.0, min(0.9, 0.5 + lift))
        except Exception:
            return 0.5

    with httpx.Client(timeout=60.0) as client:
        for r in rows:
            disease = r.get("disease") or ""
            gene = r.get("gene") or ""
            hgvs_p = r.get("hgvs_p") or ""
            chrom = r.get("chrom") or None
            pos = r.get("pos") or None
            ref = r.get("ref") or None
            alt = r.get("alt") or None
            build = r.get("build") or None
            label = parse_int(r.get("outcome_platinum", 0))

            payload = {
                "disease": disease,
                "mutations": [],
                "api_base": args.api_base,
            }
            mut = {}
            if gene: mut["gene"] = gene
            if hgvs_p: mut["hgvs_p"] = hgvs_p
            if chrom and pos and ref and alt:
                try:
                    mut.update({"chrom": str(chrom), "pos": int(pos), "ref": str(ref), "alt": str(alt)})
                except Exception:
                    pass
            if build:
                mut["build"] = build
            if mut:
                payload["mutations"].append(mut)

            score = 0.5
            details = {}
            if args.use_evo:
                score = evo_score(client, args.api_base, r)
                details = {"method": "evo_min_delta_proxy"}
            else:
                try:
                    resp = client.post(f"{args.api_base}/api/guidance/synthetic_lethality", json=payload)
                    if resp.status_code < 400:
                        js = resp.json() or {}
                        details = js
                        therapy = (js.get("suggested_therapy") or "").lower()
                        if "platinum" in therapy or "parp" in therapy:
                            score = 0.8
                        # if guidance present with tier I, bump further
                        g = js.get("guidance") or {}
                        if g.get("tier") == "I":
                            score = max(score, 0.9)
                except Exception as e:
                    details = {"error": str(e)}

            y_true.append(label)
            y_score.append(score)
            detailed.append({"input": r, "prediction": details, "score": score})

    try:
        auprc = average_precision_score(y_true, y_score)
        auroc = roc_auc_score(y_true, y_score)
    except Exception:
        auprc = 0.0
        auroc = 0.0

    out = {
        "metrics": {"auprc": auprc, "auroc": auroc},
        "n": len(rows),
        "api_base": args.api_base,
        "results": detailed,
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()


                except Exception as e:
                    details = {"error": str(e)}

            y_true.append(label)
            y_score.append(score)
            detailed.append({"input": r, "prediction": details, "score": score})

    try:
        auprc = average_precision_score(y_true, y_score)
        auroc = roc_auc_score(y_true, y_score)
    except Exception:
        auprc = 0.0
        auroc = 0.0

    out = {
        "metrics": {"auprc": auprc, "auroc": auroc},
        "n": len(rows),
        "api_base": args.api_base,
        "results": detailed,
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()




                except Exception as e:
                    details = {"error": str(e)}

            y_true.append(label)
            y_score.append(score)
            detailed.append({"input": r, "prediction": details, "score": score})

    try:
        auprc = average_precision_score(y_true, y_score)
        auroc = roc_auc_score(y_true, y_score)
    except Exception:
        auprc = 0.0
        auroc = 0.0

    out = {
        "metrics": {"auprc": auprc, "auroc": auroc},
        "n": len(rows),
        "api_base": args.api_base,
        "results": detailed,
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()


                except Exception as e:
                    details = {"error": str(e)}

            y_true.append(label)
            y_score.append(score)
            detailed.append({"input": r, "prediction": details, "score": score})

    try:
        auprc = average_precision_score(y_true, y_score)
        auroc = roc_auc_score(y_true, y_score)
    except Exception:
        auprc = 0.0
        auroc = 0.0

    out = {
        "metrics": {"auprc": auprc, "auroc": auroc},
        "n": len(rows),
        "api_base": args.api_base,
        "results": detailed,
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()



