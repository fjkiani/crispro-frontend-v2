#!/usr/bin/env python3
"""
Fetch public Multiple Myeloma cohorts from cBioPortal and stage them under:
  data/validation/mm_cohort/cbioportal_mm/

Outputs (all JSON):
- study_index.json: list of MM-like studies discovered (ids + basic metadata)
- gene_counts_by_study.json: per-study counts for key genes (SNV/indel)
- merged_patients.json: merged patient-level clinical dict (best-effort)
- merged_mutations.json: merged mutations (best-effort, can be large)
- merged_cna.json: merged discrete CNA calls (best-effort; only if supported)

Notes:
- This is intended to increase power for rare markers (e.g., PSMB5/CRBN),
  and to enable CNV-aware markers when discrete CNA profiles are available.
- Requires network access to https://www.cbioportal.org/api
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Set

from scripts.data_acquisition.utils.cbioportal_client import CBioportalClient


OUT_DIR = Path("data/validation/mm_cohort/cbioportal_mm")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Genes that matter for MM resistance mission
KEY_GENES: List[str] = [
    "PSMB5",
    "CRBN",
    "IKZF1",
    "IKZF3",
    "CUL4A",
    "DDB1",
    "GSPT1",
    "TP53",
    "DIS3",
    "KRAS",
    "NRAS",
    "BRAF",
    "NFE2L2",
    "KEAP1",
    "XBP1",
    "ERN1",   # IRE1
    "ABCB1",
    "CD38",
    "TNFRSF17",  # BCMA
]


def is_mm_study(study: Dict[str, Any]) -> bool:
    """
    Heuristic filter: keep studies whose id or name mentions myeloma.
    """
    sid = (study.get("studyId") or "").lower()
    name = (study.get("name") or "").lower()
    cancer_type = (study.get("cancerTypeId") or "").lower()
    txt = " ".join([sid, name, cancer_type])
    return any(k in txt for k in ["myeloma", "mmrf", "commpass", "multiple myeloma", "plasma cell"])


def extract_gene_set_from_mutations(muts: List[Dict[str, Any]]) -> Set[str]:
    genes: Set[str] = set()
    for m in muts:
        g = (m.get("gene") or {}).get("hugoGeneSymbol") if isinstance(m.get("gene"), dict) else m.get("gene")
        if not g:
            # some responses use "hugoGeneSymbol"
            g = m.get("hugoGeneSymbol")
        if g:
            genes.add(str(g).upper())
    return genes


def main() -> None:
    client = CBioportalClient()
    try:
        studies = client.list_studies()
        mm_studies = [s for s in studies if is_mm_study(s)]

        (OUT_DIR / "study_index.json").write_text(json.dumps(mm_studies, indent=2))
        print(f"Found {len(mm_studies)} MM-like studies. Wrote: {OUT_DIR/'study_index.json'}")

        gene_counts_by_study: Dict[str, Dict[str, Any]] = {}
        merged_mutations: List[Dict[str, Any]] = []
        merged_patients: Dict[str, Dict[str, Any]] = {}
        merged_cna: List[Dict[str, Any]] = []

        for s in mm_studies:
            study_id = s.get("studyId")
            if not study_id:
                continue
            print(f"\n=== Study: {study_id} ===")

            # patients + clinical
            patients = client.get_study_patients(study_id)
            clinical = client.get_clinical_data(study_id, entity_type="PATIENT")

            # sample ids (for mutations/cna)
            samples = client.get_study_samples(study_id)
            sample_ids = [x.get("sampleId") for x in samples if x.get("sampleId")]

            profile_id = client.get_mutation_profile_id(study_id)
            muts: List[Dict[str, Any]] = []
            if profile_id and sample_ids:
                # cBioPortal may accept large sampleIds; if not, split.
                batch = 1000
                for i in range(0, len(sample_ids), batch):
                    muts.extend(client.get_mutations_for_samples(profile_id, sample_ids[i : i + batch]))

            # best-effort CNA pulls
            cna_profile = client.get_discrete_cna_profile_id(study_id)
            cna_calls: List[Dict[str, Any]] = []
            if cna_profile and sample_ids:
                try:
                    batch = 1000
                    for i in range(0, len(sample_ids), batch):
                        cna_calls.extend(client.get_discrete_copy_number_for_samples(cna_profile, sample_ids[i : i + batch]))
                except Exception as e:
                    print(f"  ⚠️ CNA fetch skipped for {study_id}: {e}")

            genes_present = extract_gene_set_from_mutations(muts)
            counts = {g: 1 if g in genes_present else 0 for g in KEY_GENES}

            gene_counts_by_study[study_id] = {
                "n_patients": len(patients),
                "n_samples": len(sample_ids),
                "n_mutations_rows": len(muts),
                "mutation_gene_presence_any": counts,
                "mutation_genes_seen": sorted(list(genes_present)),
                "cna_profile_id": cna_profile,
                "n_cna_rows": len(cna_calls),
            }

            # merge outputs (best-effort; these can get large)
            for p in patients:
                pid = p.get("patientId")
                if not pid:
                    continue
                merged_patients[f"{study_id}::{pid}"] = {
                    "study_id": study_id,
                    "patient_id": pid,
                    "clinical": clinical.get(pid, {}),
                }
            for m in muts:
                m["__study_id"] = study_id
                merged_mutations.append(m)
            for c in cna_calls:
                c["__study_id"] = study_id
                merged_cna.append(c)

            print(f"  patients={len(patients)} samples={len(sample_ids)} muts={len(muts)} cna={len(cna_calls)}")

        (OUT_DIR / "gene_counts_by_study.json").write_text(json.dumps(gene_counts_by_study, indent=2))
        (OUT_DIR / "merged_patients.json").write_text(json.dumps(merged_patients, indent=2))
        (OUT_DIR / "merged_mutations.json").write_text(json.dumps(merged_mutations, indent=2))
        (OUT_DIR / "merged_cna.json").write_text(json.dumps(merged_cna, indent=2))

        print("\nDone:")
        print(f"- {OUT_DIR/'gene_counts_by_study.json'}")
        print(f"- {OUT_DIR/'merged_patients.json'}")
        print(f"- {OUT_DIR/'merged_mutations.json'}")
        print(f"- {OUT_DIR/'merged_cna.json'}")

    finally:
        client.close()


if __name__ == "__main__":
    main()




















