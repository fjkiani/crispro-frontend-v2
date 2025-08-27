#!/usr/bin/env python3
import argparse
import csv


def parse_bool(s: str) -> str:
    v = str(s).strip().lower()
    if v in {"1", "true", "yes", "y"}:
        return "true"
    if v in {"0", "false", "no", "n"}:
        return "false"
    return ""


def main():
    ap = argparse.ArgumentParser(description="ETL clean guidance labels CSV for tier evaluation")
    ap.add_argument("--inp", required=True, help="Input CSV path")
    ap.add_argument("--out", required=True, help="Output cleaned CSV path")
    ap.add_argument("--require_variant", action="store_true", help="Drop rows without gene+hgvs_p")
    args = ap.parse_args()

    with open(args.inp, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    cleaned = []
    for r in rows:
        disease = (r.get("disease") or "").strip()
        therapy = (r.get("therapy") or r.get("drug_or_class") or "").strip()
        gene = (r.get("gene") or "").strip()
        hgvs_p = (r.get("hgvs_p") or "").strip()
        chrom = (r.get("chrom") or "").strip()
        pos = (r.get("pos") or "").strip()
        ref = (r.get("ref") or "").strip().upper()
        alt = (r.get("alt") or "").strip().upper()
        build = (r.get("build") or "hg38").strip()
        on_label_truth = parse_bool(r.get("on_label_truth"))

        # minimal required
        if not disease or not therapy:
            continue
        if args.require_variant and not (gene and hgvs_p):
            continue

        out_row = {
            "disease": disease,
            "therapy": therapy,
            "drug_or_class": therapy,
            "gene": gene,
            "hgvs_p": hgvs_p,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "build": build,
            "on_label_truth": on_label_truth,
        }
        cleaned.append(out_row)

    fieldnames = [
        "disease","therapy","drug_or_class","gene","hgvs_p","chrom","pos","ref","alt","build","on_label_truth"
    ]
    with open(args.out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in cleaned:
            w.writerow(r)

    print(f"Wrote {len(cleaned)} rows to {args.out}")


if __name__ == "__main__":
    main()


