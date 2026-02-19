## Tumor Board Data Provenance (receipt-only)

### Cross-check rule (MANDATORY before “Verified”)

- **Verified**: Any “live response receipt” MUST include:
  - a raw minimized nested JSON excerpt copied from the response, and
  - a path→value table derived from that same excerpt.
  
If the path→value table contradicts the raw excerpt for the same request identity, the section is **FAIL** and must be regenerated.

### 1) Endpoint contract (backend)

- **Verified**: Bundle handler signature exists at `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py`:

```127:135:oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py
@router.post("/bundle")
async def ayesha_analysis_bundle(
    level: Optional[str] = Query("all", description="Level to analyze: 'l1', 'l2', 'l3', or 'all'"),
    scenario_id: Optional[str] = Query(None, description="L2 scenario ID"),
    l3_scenario_id: Optional[str] = Query(None, description="L3 scenario ID"),
    include_synthetic_lethality: bool = Query(True, description="Include SL panel"),
    ctdna_status_override: Optional[str] = Query(None, description="Force ctDNA status (e.g. RISING) for validation"),
):
    """Canonical end-to-end bundle for Ayesha Therapy Fit."""
```

- **Verified**: Bundle response shape keys are constructed in the same handler (selected excerpts):

```166:178:oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py
            out[lk.upper()] = {
                "is_preview": bool(lk != "l1"),
                "effective_assembly": _infer_effective_assembly(profile, mutations_used, tumor_context_used),
                "inputs_used": {"mutations": mutations_used, "tumor_context": tumor_context_used},
                "efficacy": {
                    "drugs": efficacy.get("drugs", []),
                    **_normalize_pathway_scores_with_meta(efficacy.get("pathway_scores", {})),
                    "provenance": efficacy.get("provenance", {}),
                },
                "resistance_gate": resistance_gate,
                "synthetic_lethality": sl_payload,
                "completeness": efficacy.get("completeness"),
            }
```

```192:202:oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_therapy_fit.py
    return {
        "contract_version": "v2.0",
        "patient_id": "AK",
        "requested_levels": [k.upper() for k in levels],
        "generated_at": datetime.now().isoformat(),
        "levels": out,
        "tests_needed": determine_required_tests(build_profile_for_level("l1")),
        "patient_context": patient_context,
        "synthetic_lethality": l1_data.get("synthetic_lethality")
    }
```

---

### 2) One live response receipt (minimized excerpt)

- **Verified**: Request identity (single run):
  - **method**: `POST`
  - **path**: `/api/ayesha/therapy-fit/bundle?level=l1&include_synthetic_lethality=true`
  - **body**: `{}`
  - **status_code**: `200`
  - **generated_at** (from response): `2026-02-17T02:37:33.434830`

- **Verified**: Raw minimized nested excerpt (copied from that response):

```json
{
  "contract_version": "v2.0",
  "patient_id": "AK",
  "requested_levels": [
    "L1"
  ],
  "generated_at": "2026-02-17T02:37:33.434830",
  "levels": {
    "L1": {
      "efficacy": {
        "drugs": [
          {
            "name": "paclitaxel",
            "efficacy_score": 0.32,
            "confidence": 0.6,
            "evidence_tier": "insufficient",
            "clinical_band": "D (Low/VUS)",
            "rationale": [
              {
                "type": "sequence",
                "value": 0.0001,
                "percentile": 0.8
              }
            ]
          }
        ]
      },
      "synthetic_lethality": {
        "synthetic_lethality_detected": false,
        "provenance": {
          "status": "failed_true_scoring_no_variants_sent"
        }
      }
    }
  }
}
```

- **Verified**: Path→value table derived from the same run (same request identity + same values):

```json
{
  "levels.L1.efficacy.drugs[0].name": "paclitaxel",
  "levels.L1.efficacy.drugs[0].efficacy_score": 0.32,
  "levels.L1.efficacy.drugs[0].rationale[0]": {
    "percentile": 0.8,
    "type": "sequence",
    "value": 0.0001
  },
  "levels.L1.synthetic_lethality.synthetic_lethality_detected": false
}
```

---

### 2b) One live response receipt (COMPREHENSIVE mode, citations-enabled)

- **Verified**: Request identity (single run):
  - **method**: `POST`
  - **path**: `/api/ayesha/therapy-fit/bundle?level=l1&include_synthetic_lethality=true&efficacy_mode=comprehensive`
  - **body**: `{}`
  - **status_code**: `200`
  - **generated_at** (from response): `2026-02-17T02:47:06.158972`

- **Verified**: Raw minimized nested excerpt (copied from that response):

```json
{"contract_version":"v2.0","generated_at":"2026-02-17T02:47:06.158972","patient_id":"AK","requested_levels":["L1"],"levels":{"L1":{"efficacy":{"drugs":[{"name":"adavosertib","efficacy_score":0.64,"confidence":0.6,"evidence_tier":"consider","clinical_band":"B (Good)","citations_count":0,"citations":[],"rationale":[{"type":"sequence","value":0.0001,"percentile":0.8}],"badges":["PathwayAligned","SL-Detected"]}]},"synthetic_lethality":{"synthetic_lethality_detected":true,"provenance":{"status":"ok","degraded_from_true_scoring":true,"degraded_reason":"failed_true_scoring_no_variants_sent"}}}}}
```

- **Verified**: Path→value table derived from the same run (same request identity + same values):

```json
{
  "levels.L1.efficacy.drugs[0].name": "adavosertib",
  "levels.L1.efficacy.drugs[0].efficacy_score": 0.64,
  "levels.L1.efficacy.drugs[0].citations_count": 0,
  "levels.L1.efficacy.drugs[0].citations": [],
  "levels.L1.efficacy.drugs[0].rationale[0]": {
    "type": "sequence",
    "value": 0.0001,
    "percentile": 0.8
  },
  "levels.L1.synthetic_lethality.synthetic_lethality_detected": true,
  "levels.L1.synthetic_lethality.provenance.status": "ok",
  "levels.L1.synthetic_lethality.provenance.degraded_from_true_scoring": true,
  "levels.L1.synthetic_lethality.provenance.degraded_reason": "failed_true_scoring_no_variants_sent"
}
```

---

### 3) Frontend fetch receipt (required)

- **Verified**: The Tumor Board page fetches the bundle at `oncology-coPilot/oncology-frontend/src/hooks/useTumorBoardBundle.js`:

```21:51:oncology-coPilot/oncology-frontend/src/hooks/useTumorBoardBundle.js
async function fetchTumorBoardBundle({
  level = 'l1',
  scenarioId = null,
  l3ScenarioId = null,
  includeSyntheticLethality = true,
  ctdnaStatusOverride = null,
} = {}) {
  const params = new URLSearchParams();
  params.set('level', String(level || 'l1'));
  if (scenarioId) params.set('scenario_id', String(scenarioId));
  if (l3ScenarioId) params.set('l3_scenario_id', String(l3ScenarioId));
  params.set('include_synthetic_lethality', includeSyntheticLethality ? 'true' : 'false');
  if (ctdnaStatusOverride) params.set('ctdna_status_override', String(ctdnaStatusOverride));

  const res = await fetch(`${API_ROOT}/api/ayesha/therapy-fit/bundle?${params.toString()}`, {
    method: 'POST',
    headers,
    body: JSON.stringify({}),
  });

  return res.json();
}
```

