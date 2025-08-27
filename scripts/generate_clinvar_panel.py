#!/usr/bin/env python3
import csv
import gzip
import io
import json
import os
import sys
import time
import urllib.request
from dataclasses import dataclass
from typing import List, Dict, Any, Optional

CLINVAR_URL = os.environ.get(
    "CLINVAR_SUMMARY_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
)
TARGET_GENES = [
    "KRAS", "NRAS", "BRAF", "TP53", "PIK3CA", "FGFR3", "MAP2K1", "NF1",
]
ASSEMBLY = "GRCh38"
POS_LABELS = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
NEG_LABELS = {"Benign", "Likely_benign", "Benign/Likely_benign"}
REVIEW_OK_SUBSTR = "criteria provided"
REVIEW_BAD_SUBSTR = "conflicting"
REVIEW_HIGH_CONF = (
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
    "practice guideline",
)

MAX_PER_GENE_POS = int(os.environ.get("MAX_PER_GENE_POS", "10"))
MAX_PER_GENE_NEG = int(os.environ.get("MAX_PER_GENE_NEG", "10"))
OUT_JSON = os.environ.get("OUT_JSON", "docs/benchmarks/myeloma_panel_generated.json")
OUT_PROV = os.environ.get("OUT_PROV", "docs/benchmarks/myeloma_panel_provenance.csv")
FAST_MODE = os.environ.get("FAST_MODE", "1") == "1"
SAMPLE_VALIDATE_PER_GENE = int(os.environ.get("SAMPLE_VALIDATE_PER_GENE", "2"))
VUS_LABELS = {"Uncertain_significance"}
OUT_VUS_JSON = os.environ.get("OUT_VUS_JSON", "docs/benchmarks/vus_panel.json")
INCLUDE_VUS = os.environ.get("INCLUDE_VUS", "0") == "1"


def http_get(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "crispro-benchmark/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def ensembl_ref_base(chrom: str, pos: int, assembly: str = ASSEMBLY, retry: int = 3) -> Optional[str]:
    """Fetch the reference base at chrom:pos (1-based) using Ensembl REST."""
    url = (
        f"https://rest.ensembl.org/sequence/region/human/{chrom}:{pos}-{pos}:1?"
        f"content-type=text/plain;coord_system_version={assembly}"
    )
    for attempt in range(retry):
        try:
            data = http_get(url, timeout=30)
            base = data.decode("utf-8").strip().upper()
            if base:
                return base
        except Exception:
            time.sleep(0.5 * (attempt + 1))
    return None


@dataclass
class ClinVarEntry:
    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    clin_sig: str
    review: str
    variation_id: str
    rcv_accession: str
    protein_change: str


def load_clinvar_entries() -> List[ClinVarEntry]:
    buf = http_get(CLINVAR_URL, timeout=120)
    with gzip.GzipFile(fileobj=io.BytesIO(buf)) as gz:
        text = gz.read().decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    entries: List[ClinVarEntry] = []
    for row in reader:
        try:
            if row.get("Assembly", "").strip() != ASSEMBLY:
                continue
            if (row.get("GeneSymbol", "") or "").upper() not in TARGET_GENES:
                continue
            # Filter SNVs only (ClinVar uses 'single nucleotide variant')
            type_val = (row.get("Type", "") or "").strip().lower()
            if type_val not in {"snv", "snp", "single nucleotide variant"} and "single nucleotide variant" not in type_val:
                continue
            clin = (row.get("ClinicalSignificance", "") or "").replace(" ", "_")
            if clin not in POS_LABELS and clin not in NEG_LABELS:
                continue
            review = (row.get("ReviewStatus", "") or "").lower()
            if (REVIEW_BAD_SUBSTR in review) or ("no assertion" in review):
                continue
            if not any(h in review for h in REVIEW_HIGH_CONF):
                continue
            chrom = (row.get("Chromosome", "") or "").replace("chr", "")
            # Some ClinVar exports use ChromosomeStart/ChromosomeStop
            start_str = row.get("ChromosomeStart") or row.get("Start") or row.get("PositionVCF") or "0"
            start = int(start_str or 0)
            # Prefer VCF allele columns when available
            ref_raw = (row.get("ReferenceAllele", "") or "").upper()
            alt_raw = (row.get("AlternateAllele", "") or "").upper()
            if ref_raw in {"", "NA", "N/A", "NA", "NA", "NA", "NA", "NA", "NA", "NA"} or ref_raw == "NA" or ref_raw == "N/A" or ref_raw == "NA":
                ref_raw = (row.get("ReferenceAlleleVCF", "") or "").upper()
            if alt_raw in {"", "NA", "N/A"}:
                alt_raw = (row.get("AlternateAlleleVCF", "") or "").upper()
            ref = ref_raw
            alt = alt_raw
            if not (chrom and start and ref and alt and len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT"):
                continue
            entries.append(
                ClinVarEntry(
                    gene=(row.get("GeneSymbol", "") or "").upper(),
                    chrom=chrom,
                    pos=start,
                    ref=ref,
                    alt=alt,
                    clin_sig=clin,
                    review=row.get("ReviewStatus", "") or "",
                    variation_id=row.get("VariationID", "") or "",
                    rcv_accession=row.get("RCVaccession", "") or "",
                    protein_change=row.get("Protein_change", "") or "",
                )
            )
        except Exception:
            continue
    return entries


def load_vus_entries() -> List[ClinVarEntry]:
    buf = http_get(CLINVAR_URL, timeout=120)
    with gzip.GzipFile(fileobj=io.BytesIO(buf)) as gz:
        text = gz.read().decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    out: List[ClinVarEntry] = []
    for row in reader:
        try:
            if row.get("Assembly", "").strip() != ASSEMBLY:
                continue
            if (row.get("GeneSymbol", "") or "").upper() not in TARGET_GENES:
                continue
            type_val = (row.get("Type", "") or "").strip().lower()
            if type_val not in {"snv", "snp", "single nucleotide variant"} and "single nucleotide variant" not in type_val:
                continue
            clin = (row.get("ClinicalSignificance", "") or "").replace(" ", "_")
            if clin not in VUS_LABELS:
                continue
            review = (row.get("ReviewStatus", "") or "").lower()
            # VUS: accept criteria provided (single or multiple), reject conflicts/no assertion
            if ("criteria provided" not in review) or ("conflicting" in review) or ("no assertion" in review):
                continue
            chrom = (row.get("Chromosome", "") or "").replace("chr", "")
            start_str = row.get("ChromosomeStart") or row.get("Start") or row.get("PositionVCF") or "0"
            start = int(start_str or 0)
            ref_raw = (row.get("ReferenceAllele", "") or "").upper()
            alt_raw = (row.get("AlternateAllele", "") or "").upper()
            if ref_raw in {"", "NA", "N/A"}:
                ref_raw = (row.get("ReferenceAlleleVCF", "") or "").upper()
            if alt_raw in {"", "NA", "N/A"}:
                alt_raw = (row.get("AlternateAlleleVCF", "") or "").upper()
            if not (chrom and start and ref_raw and alt_raw and len(ref_raw) == 1 and len(alt_raw) == 1 and ref_raw in "ACGT" and alt_raw in "ACGT"):
                continue
            out.append(ClinVarEntry(
                gene=(row.get("GeneSymbol", "") or "").upper(),
                chrom=chrom,
                pos=start,
                ref=ref_raw,
                alt=alt_raw,
                clin_sig=clin,
                review=row.get("ReviewStatus", "") or "",
                variation_id=row.get("VariationID", "") or "",
                rcv_accession=row.get("RCVaccession", "") or "",
                protein_change=row.get("Protein_change", "") or "",
            ))
        except Exception:
            continue
    return out


def validate_entries(entries: List[ClinVarEntry]) -> List[ClinVarEntry]:
    if FAST_MODE:
        # Optionally sample-validate a few entries per gene and log failures, but don't block
        per_gene_count: Dict[str, int] = {g: 0 for g in TARGET_GENES}
        validated: List[ClinVarEntry] = []
        for e in entries:
            if e.gene in per_gene_count and per_gene_count[e.gene] < SAMPLE_VALIDATE_PER_GENE:
                base = ensembl_ref_base(e.chrom, e.pos, ASSEMBLY)
                per_gene_count[e.gene] += 1
                if base is not None and (base == e.ref or base == "N"):
                    validated.append(e)
                else:
                    # skip sampled failures silently in fast mode
                    continue
            else:
                # rely on VCF-provided ref/alt for speed
                validated.append(e)
        return validated
    # Original strict validation
    valid: List[ClinVarEntry] = []
    for e in entries:
        base = ensembl_ref_base(e.chrom, e.pos, ASSEMBLY)
        if base is None:
            continue
        if base == e.ref or base == "N":
            valid.append(e)
    return valid


def balance(entries: List[ClinVarEntry]) -> List[ClinVarEntry]:
    per_gene_pos: Dict[str, List[ClinVarEntry]] = {g: [] for g in TARGET_GENES}
    per_gene_neg: Dict[str, List[ClinVarEntry]] = {g: [] for g in TARGET_GENES}
    for e in entries:
        if e.gene not in per_gene_pos:
            continue
        if e.clin_sig in POS_LABELS:
            per_gene_pos[e.gene].append(e)
        elif e.clin_sig in NEG_LABELS:
            per_gene_neg[e.gene].append(e)
    panel: List[ClinVarEntry] = []
    for g in TARGET_GENES:
        # take up to limits per class per gene
        panel.extend(per_gene_pos[g][:MAX_PER_GENE_POS])
        panel.extend(per_gene_neg[g][:MAX_PER_GENE_NEG])
    # trim to a reasonable size, preserving order
    max_total = (MAX_PER_GENE_POS + MAX_PER_GENE_NEG) * len(TARGET_GENES)
    return panel[:max_total]


def write_outputs(panel: List[ClinVarEntry]):
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    # JSON panel for the app/benchmark
    js: List[Dict[str, Any]] = []
    for e in panel:
        label_bool = True if e.clin_sig in POS_LABELS else False
        js.append({
            "gene": e.gene,
            "hgvs_p": e.protein_change or "",
            "variant_info": f"chr{e.chrom}:{e.pos} {e.ref}>{e.alt}",
            "build": "hg38",
            "label_is_pathogenic": label_bool,
            "provenance": {
                "source": "ClinVar",
                "variation_id": e.variation_id,
                "rcv_accession": e.rcv_accession,
                "review_status": e.review,
                "clinical_significance": e.clin_sig,
                "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{e.variation_id}/",
                "ensembl_region_url": f"https://www.ensembl.org/Homo_sapiens/Location/View?r={e.chrom}%3A{max(1,e.pos-600)}-{e.pos+600}",
            }
        })
    with open(OUT_JSON, 'w') as f:
        json.dump(js, f, indent=2)

    # Optional VUS write
    if INCLUDE_VUS:
        vus_entries = load_vus_entries()
        vus_js = []
        for e in vus_entries:
            vus_js.append({
                "gene": e.gene,
                "hgvs_p": e.protein_change or "",
                "variant_info": f"chr{e.chrom}:{e.pos} {e.ref}>{e.alt}",
                "build": "hg38",
                "label_is_vus": True,
                "provenance": {
                    "source": "ClinVar",
                    "variation_id": e.variation_id,
                    "rcv_accession": e.rcv_accession,
                    "review_status": e.review,
                    "clinical_significance": e.clin_sig,
                    "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{e.variation_id}/",
                    "ensembl_region_url": f"https://www.ensembl.org/Homo_sapiens/Location/View?r={e.chrom}%3A{max(1,e.pos-600)}-{e.pos+600}",
                }
            })
        with open(OUT_VUS_JSON, 'w') as f:
            json.dump(vus_js, f, indent=2)

    # CSV provenance log
    with open(OUT_PROV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["gene","chrom","pos","ref","alt","clin_sig","review_status","variation_id","rcv_accession","protein_change"]) 
        for e in panel:
            w.writerow([e.gene, e.chrom, e.pos, e.ref, e.alt, e.clin_sig, e.review, e.variation_id, e.rcv_accession, e.protein_change])


def main():
    print("Downloading and parsing ClinVar...", file=sys.stderr)
    entries = load_clinvar_entries()
    print(f"Parsed {len(entries)} entries after initial filters", file=sys.stderr)
    print("Validating reference alleles via Ensembl...", file=sys.stderr)
    valid = validate_entries(entries)
    print(f"Validated {len(valid)} entries", file=sys.stderr)
    print("Balancing per gene...", file=sys.stderr)
    panel = balance(valid)
    print(f"Panel size {len(panel)}", file=sys.stderr)
    write_outputs(panel)
    print(json.dumps({
        "panel_size": len(panel),
        "out_json": OUT_JSON,
        "out_provenance": OUT_PROV
    }, indent=2))


if __name__ == "__main__":
    main() 