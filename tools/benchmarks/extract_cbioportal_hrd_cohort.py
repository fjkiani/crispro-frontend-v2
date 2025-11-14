#!/usr/bin/env python3
import argparse
import csv
import json
from typing import List, Dict
import httpx
import os


"""
Extract a small BRCA1/2 + platinum-exposed cohort from cBioPortal (e.g., TCGA OV) and
emit a CSV suitable for hrd_platinum_auprc.py.

Note: cBioPortal clinical response fields vary by study; we approximate outcome by drug exposure containing
"platinum" (carboplatin/cisplatin/oxaliplatin) and fallback to presence of exposure.
"""


CBIO_BASE = "https://www.cbioportal.org/api"


def _headers() -> Dict[str, str]:
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    token = os.getenv("CBIO_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers


def _get_samples(study_id: str) -> List[Dict]:
    with httpx.Client(timeout=60.0, headers=_headers()) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/samples")
        r.raise_for_status()
        return r.json() or []


def _post_json(url: str, data: Dict) -> Dict:
    with httpx.Client(timeout=120.0, headers=_headers()) as client:
        r = client.post(url, json=data)
        r.raise_for_status()
        try:
            return r.json() or {}
        except Exception:
            # Return text payload for diagnostics
            return {"raw": r.text}


def _get_json(url: str, params: Dict) -> Dict:
    with httpx.Client(timeout=120.0, headers=_headers()) as client:
        r = client.get(url, params=params)
        r.raise_for_status()
        try:
            return r.json() or {}
        except Exception:
            return {"raw": r.text}


def fetch_mutations_for_genes(study_id: str, genes: List[str], sample_ids: List[str]) -> Dict[str, List[Dict]]:
    # Use unauthenticated GET endpoint for the default mutations profile with the all-samples list
    profile_id = f"{study_id}_mutations"
    sample_list_id = f"{study_id}_all"
    params = {
        "sampleListId": sample_list_id,
        # Pass as repeated query params by using a list; httpx will encode hugoGeneSymbols=BRCA1&hugoGeneSymbols=BRCA2
        "hugoGeneSymbols": genes,
        "projection": "DETAILED",
        "pageSize": 10000,
    }
    js = _get_json(f"{CBIO_BASE}/molecular-profiles/{profile_id}/mutations", params)
    muts = js if isinstance(js, list) else js.get("result", js.get("items", []))
    by_sample: Dict[str, List[Dict]] = {}
    for m in muts:
        sid = m.get("sampleId") or m.get("sampleId")
        if sid:
            by_sample.setdefault(sid, []).append(m)
    return by_sample


def fetch_clinical_samples(study_id: str, sample_ids: List[str]) -> Dict[str, Dict]:
    rows: Dict[str, Dict] = {}
    # Prefer POST /clinical-data/fetch; chunk to avoid payload issues
    chunk_size = 200
    any_success = False
    for i in range(0, len(sample_ids), chunk_size):
        chunk = sample_ids[i : i + chunk_size]
        data = {
            "entityIds": chunk,
            "entityType": "SAMPLE",
            "projection": "DETAILED",
            "attributeIds": [
                "DRUG_NAME",
                "TREATMENT_TYPE",
                "THERAPY_NAME",
                "CLINICAL_TREATMENT_TYPE",
                "PHARMACEUTICAL_TX_TYPE",
            ],
        }
        js = _post_json(f"{CBIO_BASE}/clinical-data/fetch", data)
        def _ingest(row: Dict):
            sid = row.get("entityId") if isinstance(row, dict) else None
            if not sid:
                return
            existing = rows.setdefault(sid, {})
            # Merge fields shallowly
            existing.update(row)
            # Aggregate all text for robust exposure labeling later
            import json as _json
            agg = existing.get("__combined_text", "")
            existing["__combined_text"] = (agg + "\n" + _json.dumps(row)).strip()
        if isinstance(js, list):
            for row in js:
                _ingest(row)
                any_success = True
        elif isinstance(js, dict) and isinstance(js.get("result"), list):
            for row in js["result"]:
                _ingest(row)
                any_success = True
    return rows


def fetch_clinical_patients(study_id: str, patient_ids: List[str]) -> Dict[str, Dict]:
    rows: Dict[str, Dict] = {}
    # PATIENT-level clinical often carries treatment/exposure
    chunk_size = 200
    for i in range(0, len(patient_ids), chunk_size):
        chunk = patient_ids[i : i + chunk_size]
        data = {
            "entityIds": chunk,
            "entityType": "PATIENT",
            "projection": "DETAILED",
            "attributeIds": [
                "DRUG_NAME",
                "TREATMENT_TYPE",
                "THERAPY_NAME",
                "CLINICAL_TREATMENT_TYPE",
                "PHARMACEUTICAL_TX_TYPE",
                "TREATMENTS",
                "CHEMOTHERAPY",
                "TREATMENT_SUMMARY",
            ],
        }
        js = _post_json(f"{CBIO_BASE}/clinical-data/fetch", data)
        def _ingest(row: Dict):
            pid = row.get("entityId") if isinstance(row, dict) else None
            if not pid:
                return
            existing = rows.setdefault(pid, {})
            existing.update(row)
            import json as _json
            agg = existing.get("__combined_text", "")
            existing["__combined_text"] = (agg + "\n" + _json.dumps(row)).strip()
        if isinstance(js, list):
            for row in js:
                _ingest(row)
        elif isinstance(js, dict) and isinstance(js.get("result"), list):
            for row in js["result"]:
                _ingest(row)
    return rows


def label_platinum_exposure(clin_row: Dict) -> int:
    # Heuristic: search aggregated clinical text for platinum-related drugs
    combined = (clin_row or {}).get("__combined_text")
    txt = (combined if isinstance(combined, str) and combined else json.dumps(clin_row or {})).lower()
    return 1 if any(k in txt for k in ["carboplatin", "cisplatin", "oxaliplatin", "platinum"]) else 0


def main():
    ap = argparse.ArgumentParser(description="Extract HRD/platinum cohort from cBioPortal study")
    ap.add_argument("--study", default="ov_tcga", help="cBioPortal study ID (e.g., ov_tcga)")
    ap.add_argument("--out", default="tools/benchmarks/data/hrd_cohort_cbioportal.csv")
    args = ap.parse_args()

    samples = _get_samples(args.study)
    sample_ids = [s.get("sampleId") for s in samples if s.get("sampleId")]
    # Build sample->patient map when available
    samp_to_patient = {}
    for s in samples:
        sid = s.get("sampleId")
        pid = s.get("patientId") or s.get("patient_id") or s.get("patient")
        if sid and pid:
            samp_to_patient[sid] = str(pid)
    if not sample_ids:
        raise SystemExit("No samples returned")

    muts = fetch_mutations_for_genes(args.study, ["BRCA1", "BRCA2"], sample_ids)
    # Prefer PATIENT-level clinical when mapping exists; else fall back to SAMPLE-level
    patient_ids = list({p for p in samp_to_patient.values() if p})
    clin_patient = fetch_clinical_patients(args.study, patient_ids) if patient_ids else {}
    clin_sample = fetch_clinical_samples(args.study, sample_ids) if not clin_patient else {}

    out_rows: List[Dict] = []
    for sid in sample_ids:
        # Choose clinical row: patient-level if mapped, else sample-level
        pid = samp_to_patient.get(sid)
        row = (clin_patient.get(pid) if pid else None) or clin_sample.get(sid, {})
        # Use case ID when available
        disease = (row.get("cancerType") or row.get("CANCER_TYPE") or "ovarian cancer")
        # Collect gene/hgvs if present in mutations
        ms = muts.get(sid, [])
        gene = ""
        hgvs = ""
        chrom = pos = ref = alt = ""
        if ms:
            m = ms[0]
            gene = (m.get("gene") or {}).get("hugoGeneSymbol") or ""
            hgvs = m.get("proteinChange") or ""
            # Best-effort coordinates if provided by API
            chrom = str(m.get("chromosome") or m.get("chr") or "")
            try:
                pos = str(int(m.get("startPosition") or m.get("startPosition") or 0)) if (m.get("startPosition") or 0) else ""
            except Exception:
                pos = ""
            ref = str(m.get("referenceAllele") or "").upper()
            alt = str(m.get("variantAllele") or "").upper()
        outcome = label_platinum_exposure(row)
        out_rows.append({
            "disease": disease,
            "gene": gene,
            "hgvs_p": hgvs,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "build": "GRCh38",
            "outcome_platinum": outcome,
        })

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "disease","gene","hgvs_p","chrom","pos","ref","alt","build","outcome_platinum"
        ])
        w.writeheader()
        w.writerows(out_rows[:200])

    print(f"Wrote {len(out_rows[:200])} rows to {args.out}")


if __name__ == "__main__":
    main()



        chunk = patient_ids[i : i + chunk_size]
        data = {
            "entityIds": chunk,
            "entityType": "PATIENT",
            "projection": "DETAILED",
            "attributeIds": [
                "DRUG_NAME",
                "TREATMENT_TYPE",
                "THERAPY_NAME",
                "CLINICAL_TREATMENT_TYPE",
                "PHARMACEUTICAL_TX_TYPE",
                "TREATMENTS",
                "CHEMOTHERAPY",
                "TREATMENT_SUMMARY",
            ],
        }
        js = _post_json(f"{CBIO_BASE}/clinical-data/fetch", data)
        def _ingest(row: Dict):
            pid = row.get("entityId") if isinstance(row, dict) else None
            if not pid:
                return
            existing = rows.setdefault(pid, {})
            existing.update(row)
            import json as _json
            agg = existing.get("__combined_text", "")
            existing["__combined_text"] = (agg + "\n" + _json.dumps(row)).strip()
        if isinstance(js, list):
            for row in js:
                _ingest(row)
        elif isinstance(js, dict) and isinstance(js.get("result"), list):
            for row in js["result"]:
                _ingest(row)
    return rows


def label_platinum_exposure(clin_row: Dict) -> int:
    # Heuristic: search aggregated clinical text for platinum-related drugs
    combined = (clin_row or {}).get("__combined_text")
    txt = (combined if isinstance(combined, str) and combined else json.dumps(clin_row or {})).lower()
    return 1 if any(k in txt for k in ["carboplatin", "cisplatin", "oxaliplatin", "platinum"]) else 0


def main():
    ap = argparse.ArgumentParser(description="Extract HRD/platinum cohort from cBioPortal study")
    ap.add_argument("--study", default="ov_tcga", help="cBioPortal study ID (e.g., ov_tcga)")
    ap.add_argument("--out", default="tools/benchmarks/data/hrd_cohort_cbioportal.csv")
    args = ap.parse_args()

    samples = _get_samples(args.study)
    sample_ids = [s.get("sampleId") for s in samples if s.get("sampleId")]
    # Build sample->patient map when available
    samp_to_patient = {}
    for s in samples:
        sid = s.get("sampleId")
        pid = s.get("patientId") or s.get("patient_id") or s.get("patient")
        if sid and pid:
            samp_to_patient[sid] = str(pid)
    if not sample_ids:
        raise SystemExit("No samples returned")

    muts = fetch_mutations_for_genes(args.study, ["BRCA1", "BRCA2"], sample_ids)
    # Prefer PATIENT-level clinical when mapping exists; else fall back to SAMPLE-level
    patient_ids = list({p for p in samp_to_patient.values() if p})
    clin_patient = fetch_clinical_patients(args.study, patient_ids) if patient_ids else {}
    clin_sample = fetch_clinical_samples(args.study, sample_ids) if not clin_patient else {}

    out_rows: List[Dict] = []
    for sid in sample_ids:
        # Choose clinical row: patient-level if mapped, else sample-level
        pid = samp_to_patient.get(sid)
        row = (clin_patient.get(pid) if pid else None) or clin_sample.get(sid, {})
        # Use case ID when available
        disease = (row.get("cancerType") or row.get("CANCER_TYPE") or "ovarian cancer")
        # Collect gene/hgvs if present in mutations
        ms = muts.get(sid, [])
        gene = ""
        hgvs = ""
        chrom = pos = ref = alt = ""
        if ms:
            m = ms[0]
            gene = (m.get("gene") or {}).get("hugoGeneSymbol") or ""
            hgvs = m.get("proteinChange") or ""
            # Best-effort coordinates if provided by API
            chrom = str(m.get("chromosome") or m.get("chr") or "")
            try:
                pos = str(int(m.get("startPosition") or m.get("startPosition") or 0)) if (m.get("startPosition") or 0) else ""
            except Exception:
                pos = ""
            ref = str(m.get("referenceAllele") or "").upper()
            alt = str(m.get("variantAllele") or "").upper()
        outcome = label_platinum_exposure(row)
        out_rows.append({
            "disease": disease,
            "gene": gene,
            "hgvs_p": hgvs,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "build": "GRCh38",
            "outcome_platinum": outcome,
        })

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "disease","gene","hgvs_p","chrom","pos","ref","alt","build","outcome_platinum"
        ])
        w.writeheader()
        w.writerows(out_rows[:200])

    print(f"Wrote {len(out_rows[:200])} rows to {args.out}")


if __name__ == "__main__":
    main()


