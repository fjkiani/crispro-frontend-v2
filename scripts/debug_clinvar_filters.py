#!/usr/bin/env python3
import csv
import gzip
import io
import os
import sys
import urllib.request
from collections import Counter

CLINVAR_URL = os.environ.get(
    "CLINVAR_SUMMARY_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
)
TARGET_GENES = {"KRAS","NRAS","BRAF","TP53","PIK3CA","FGFR3","MAP2K1","NF1"}
ASSEMBLY = "GRCh38"


def http_get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "crispro-benchmark/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def main():
    buf = http_get(CLINVAR_URL)
    with gzip.GzipFile(fileobj=io.BytesIO(buf)) as gz:
        text = gz.read().decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text), delimiter='\t')

    total = 0
    cols = reader.fieldnames or []
    print("Columns:", ", ".join(cols))

    c_assembly = Counter()
    c_type = Counter()
    c_review = Counter()
    c_genes = Counter()

    kept_after_assembly = 0
    kept_after_genes = 0
    kept_after_type = 0
    kept_after_review = 0

    examples = []

    for row in reader:
        total += 1
        asm = (row.get("Assembly", "") or "").strip()
        c_assembly[asm] += 1
        gene = (row.get("GeneSymbol", "") or "").upper()
        c_genes[gene] += 1
        typ = (row.get("Type", "") or "").strip().lower()
        c_type[typ] += 1
        review = (row.get("ReviewStatus", "") or "").lower()
        c_review[review] += 1

        if asm != ASSEMBLY:
            continue
        kept_after_assembly += 1
        if gene not in TARGET_GENES:
            continue
        kept_after_genes += 1
        if not (typ in {"snv","snp","single nucleotide variant"} or "single nucleotide variant" in typ):
            continue
        kept_after_type += 1
        if ("conflicting" in review) or ("no assertion" in review):
            continue
        if not any(h in review for h in [
            "reviewed by expert panel",
            "criteria provided, multiple submitters, no conflicts",
            "practice guideline",
        ]):
            continue
        kept_after_review += 1
        if len(examples) < 5:
            examples.append({k: row.get(k) for k in [
                "GeneSymbol","Chromosome","Start","ChromosomeStart","ReferenceAllele","AlternateAllele","ClinicalSignificance","ReviewStatus","VariationID","RCVaccession","Protein_change","Type","Assembly"
            ]})

    print(f"Total rows: {total}")
    print(f"After Assembly={ASSEMBLY}: {kept_after_assembly}")
    print(f"After Gene filter: {kept_after_genes}")
    print(f"After Type (SNV): {kept_after_type}")
    print(f"After Review (criteria provided, no conflict): {kept_after_review}")

    print("Top assemblies:", c_assembly.most_common(5))
    print("Top types:", c_type.most_common(10))
    print("Top reviews:", c_review.most_common(10))
    print("Top target genes present:", [(g, c_genes[g]) for g in TARGET_GENES])

    print("\nExamples passing all filters (5 max):")
    for ex in examples:
        print(ex)

if __name__ == "__main__":
    main() 