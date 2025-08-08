# CommandCenter API Documentation

This document specifies the key API endpoints exposed by the `CommandCenter` service.

**Base URL:** `https://alpha-command-center.modal.run` (Note: This may change based on deployment.)

---

## `/workflow/assess_threat`

This is the primary endpoint for the "Intelligence Fusion" workflow. It provides a comprehensive, multi-source assessment of a single genetic variant.

*   **Method:** `POST`
*   **Summary:** Run Hardened Threat Assessment
*   **Description:** Fuses data from the internal `ThreatMatrix` database and external services like Ensembl VEP to generate a complete "Threat Dossier" for a given gene.

### Request Body

The endpoint expects a JSON object with the following structure:

```json
{
  "gene_symbol": "string",
  "protein_change": "string (optional)"
}
```

*   `gene_symbol` (string, required): The official symbol of the gene to assess (e.g., "ASXL1", "BRAF").
*   `protein_change` (string, optional): The specific protein change to assess (e.g., "p.Gly646fs"). If omitted, the system will retrieve a representative high-impact variant for the specified gene from the database.

### Success Response (200 OK)

On success, the endpoint returns a `ThreatDossier` object, which is a rich JSON containing the fused intelligence.

*   **Content-Type:** `application/json`

#### `ThreatDossier` Object Structure

```json
{
  "gene_symbol": "string",
  "protein_change": "string",
  "variant_profile": {
    "mutation_id": "string",
    "gene_symbol": "string",
    "mutation_cds": "string",
    "mutation_aa": "string",
    "mutation_description": "string",
    "genome_position": "string",
    "primary_site": "string",
    "primary_histology": "string",
    "pubmed_pmid": "integer"
  },
  "vep_dossier": {
    "most_severe_consequence": "string",
    "clinvar_summary": {},
    "prediction_summary": {},
    "cancer_associations": {}
  },
  "clinical_trials": [
    {
      "trial_id": "integer",
      "gene_id": "integer",
      "nct_id": "string",
      "phase": "string",
      "therapies": "string",
      "title": "string",
      "recruitment_status": "string"
    }
  ],
  "efficacy_evidence": [
    {
      "evidence_id": "integer",
      "gene_id": "integer",
      "molecular_profile": "string",
      "indication": "string",
      "response_type": "string",
      "therapy_name": "string",
      "approval_status": "string",
      "evidence_type": "string",
      "efficacy_description": "string",
      "reference_link": "string"
    }
  ],
  "literature_summary": {
    "id": "integer",
    "pubmed_id": "integer",
    "prompt": "string",
    "summary": "string",
    "model_used": "string",
    "analyzed_at": "string"
  }
}
```

### Error Responses

*   **400 Bad Request:** If the `gene_symbol` is missing.
*   **422 Unprocessable Entity:** If the request body is not a valid JSON.
*   **500 Internal Server Error:** If the `CommandCenter` encounters an unhandled error during workflow execution (e.g., cannot connect to the database). The response body will typically contain `{"error": "Error message"}`. 