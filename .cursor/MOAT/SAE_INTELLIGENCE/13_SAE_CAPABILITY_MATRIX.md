# SAE Capability Matrix (Proxy)

**Purpose**: Clearly define which clinical questions the Proxy SAE system can answer reliably.
**Status**: Jan 2026 Audit Corrected.

## 1. High Capability (Validated ✅)

| Capability | Integration Method | Validation Source | Clinical Value |
| :--- | :--- | :--- | :--- |
| **Transcriptomic Risk** | MFAP4 Expression | GSE63885 | **HIGH**: Predicts Platinum Resistance (AUROC 0.763). |
| **Variant Impact** | Pathway Scores | ClinVar / COSMIC | **HIGH**: Identifies driver mutations. |
| **Functional Annotation** | Insights Bundle | UniProt | **HIGH**: Quantifies protein-level loss. |
| **Pathway Analysis** | Gene Aggregation | KEGG / Reactome | **HIGH**: Maps mutations to vulnerabilities (DDR, MAPK). |
| **Immunogenicity** | TMB/MSI | FDA Criteria | **HIGH**: Determines IO eligibility (TMB >= 20). |
| **Drug Prediction** | Mechanism Vector | Literature | **HIGH**: Ranks drugs by pathway fit (e.g., PARP for DDR). |

## 2. Moderate Capability (Partial ⚠️)

| Capability | Integration Method | Validation Source | Limitations |
| :--- | :--- | :--- | :--- |
| **Trial Matching** | Mechanism Fit | Internal | Limited by trial availability for rare profiles. |
| **Metastasis Prediction** | Repair Trends | Logic-based | Model is heuristic, not ML-based. |

## 3. Low Capability (Do Not Rely ❌)

| Capability | Integration Method | Limitations |
| :--- | :--- | :--- |
| **Nutritional Therapy** | Pathway Alignment | Evidence is weak/mixed. |
| **Precise Efficacy %** | Heuristic | Estimates are categorical (High/Med/Low), not precise probabilities. |
