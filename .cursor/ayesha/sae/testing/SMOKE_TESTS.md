# ⚔️ SAE Smoke Tests - End-to-End Validation

## **🎯 Objective**

Validate that SAE features are correctly extracted from real data, returned via API, and displayed in the frontend for Ayesha's use case.

---

## **🧪 Test Suite**

### **Test 1: BRAF V600E - Known Hotspot** ✅ READY TO RUN

**Clinical Context**: BRAF V600E is a canonical hotspot in melanoma/colorectal cancer. Should trigger:
- ✅ `exon_disruption` (hotspot floor applied)
- ✅ `hotspot_mutation` (from hotspot_calibration)
- ✅ `essentiality_signal` (if BRAF is essential)

**Backend Test**:
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal

# Start backend
venv/bin/uvicorn api.main:app --host 127.0.0.1 --port 8000 &

# Wait 3 seconds for startup
sleep 3

# Run test
curl -X POST http://127.0.0.1:8000/api/clinical_genomics/analyze_variant \
  -H 'Content-Type: application/json' \
  -d '{
    "mutations": [
      {
        "gene": "BRAF",
        "chrom": "7",
        "pos": 140753336,
        "ref": "T",
        "alt": "A",
        "build": "GRCh38",
        "hgvs_p": "V600E"
      }
    ],
    "profile": "richer",
    "options": {
      "include_sae_features": true
    }
  }' | python3 -m json.tool > /tmp/braf_sae_test.json

# Extract SAE features
cat /tmp/braf_sae_test.json | jq '.efficacy.sae_features'
```

**Expected Output**:
```json
{
  "features": [
    {
      "id": "exon_disruption",
      "name": "Exon Disruption",
      "activation": 0.85,  // High due to hotspot floor
      "impact": "positive",
      "explanation": "Variant significantly disrupts exon structure",
      "provenance": "hotspot_calibration",
      "threshold": 0.5
    },
    {
      "id": "hotspot_mutation",
      "name": "Known Hotspot",
      "activation": 0.95,  // Very high for V600E
      "impact": "positive",
      "explanation": "Variant matches known pathogenic hotspot pattern",
      "provenance": "hotspot_calibration",
      "threshold": 0.5
    },
    {
      "id": "essentiality_signal",
      "name": "Gene Essentiality",
      "activation": 0.35,  // BRAF moderately essential
      "impact": "negative",
      "explanation": "Variant affects moderately essential gene",
      "provenance": "evo2_essentiality_endpoint",
      "threshold": 0.7
    }
  ],
  "boosting_features": ["exon_disruption", "hotspot_mutation"],
  "limiting_features": ["essentiality_signal"],
  "overall_impact": 0.48,
  "provenance": {
    "method": "real_data_transformation",
    "data_sources": ["hotspot_calibration", "evo2_essentiality_endpoint"],
    "feature_count": 3,
    "boosting_count": 2,
    "limiting_count": 1,
    "gene": "BRAF"
  }
}
```

**Acceptance Criteria**:
- ✅ At least 2 boosting features present
- ✅ `exon_disruption` activation ≥ 0.8 (hotspot floor)
- ✅ `hotspot_mutation` activation ≥ 0.9 (canonical hotspot)
- ✅ `overall_impact` > 0.3 (net positive boost)

---

### **Test 2: BRCA2 Pathogenic - Ayesha's Case** ✅ READY TO RUN

**Clinical Context**: BRCA2 pathogenic variant in ovarian cancer. Should trigger:
- ✅ `exon_disruption` (if large delta)
- ✅ `DNA_repair_capacity` (via toxicity pathway overlap with platinum agents)
- ✅ `essentiality_signal` (BRCA2 is highly essential)

**Backend Test**:
```bash
curl -X POST http://127.0.0.1:8000/api/clinical_genomics/analyze_variant \
  -H 'Content-Type: application/json' \
  -d '{
    "mutations": [
      {
        "gene": "BRCA2",
        "chrom": "13",
        "pos": 32936732,
        "ref": "C",
        "alt": "T",
        "build": "GRCh38"
      }
    ],
    "profile": "richer",
    "options": {
      "include_sae_features": true,
      "disease": "ovarian_cancer",
      "candidate_drugs": ["PARP_inhibitor", "platinum_agent"]
    }
  }' | python3 -m json.tool > /tmp/brca2_sae_test.json

cat /tmp/brca2_sae_test.json | jq '.efficacy.sae_features'
```

**Expected Output**:
```json
{
  "features": [
    {
      "id": "exon_disruption",
      "activation": 0.75,
      "impact": "positive",
      "explanation": "Variant moderately disrupts exon structure"
    },
    {
      "id": "DNA_repair_capacity",
      "activation": 0.85,
      "impact": "positive",
      "explanation": "DNA repair pathway burden detected (score: 0.85)",
      "provenance": "toxicity_pathway_mapping"
    },
    {
      "id": "essentiality_signal",
      "activation": 0.92,
      "impact": "positive",
      "explanation": "Variant affects essential gene (dependency: 0.92)",
      "provenance": "evo2_essentiality_endpoint"
    }
  ],
  "boosting_features": ["exon_disruption", "DNA_repair_capacity", "essentiality_signal"],
  "limiting_features": [],
  "overall_impact": 0.84,
  "provenance": {
    "gene": "BRCA2"
  }
}
```

**Acceptance Criteria**:
- ✅ `DNA_repair_capacity` present (key for PARP inhibitor rationale)
- ✅ `essentiality_signal` ≥ 0.9 (BRCA2 is critical)
- ✅ `overall_impact` > 0.7 (strong confidence boost)
- ✅ No limiting features (clean profile)

---

### **Test 3: VUS (Unknown Variant) - Minimal Features** ✅ READY TO RUN

**Clinical Context**: Variant of Uncertain Significance with no hotspot, low Evo2 delta, no cohort data.

**Backend Test**:
```bash
curl -X POST http://127.0.0.1:8000/api/clinical_genomics/analyze_variant \
  -H 'Content-Type: application/json' \
  -d '{
    "mutations": [
      {
        "gene": "TP53",
        "chrom": "17",
        "pos": 7577548,
        "ref": "C",
        "alt": "T",
        "build": "GRCh38"
      }
    ],
    "profile": "baseline",
    "options": {
      "include_sae_features": true
    }
  }' | python3 -m json.tool > /tmp/vus_sae_test.json

cat /tmp/vus_sae_test.json | jq '.efficacy.sae_features'
```

**Expected Output**:
```json
{
  "features": [
    {
      "id": "exon_disruption",
      "activation": 0.15,
      "impact": "negative",
      "explanation": "Variant minimally disrupts exon structure"
    },
    {
      "id": "cohort_overlap",
      "activation": 0.0,
      "impact": "negative",
      "explanation": "No cohort validation available (limited real-world data)"
    }
  ],
  "boosting_features": [],
  "limiting_features": ["exon_disruption", "cohort_overlap"],
  "overall_impact": -0.25,
  "provenance": {
    "gene": "TP53"
  }
}
```

**Acceptance Criteria**:
- ✅ ≤2 features (minimal signal)
- ✅ All features are limiting (no boosting)
- ✅ `overall_impact` < 0 (confidence penalty)
- ✅ Frontend displays "N/A" for missing features gracefully

---

## **🎨 Frontend Smoke Tests**

### **Test 4: SAEFeaturesCard Rendering** ✅ READY TO RUN

**Steps**:
1. Start frontend: `cd oncology-coPilot/oncology-frontend && npm run dev`
2. Navigate to `http://localhost:5173/clinical-genomics`
3. Click "Mechanistic Evidence" tab
4. Enter BRAF V600E variant:
   - Gene: BRAF
   - Chromosome: 7
   - Position: 140753336
   - Ref: T
   - Alt: A
5. Select "Richer S" profile
6. Click "🧬 Run Deep Analysis"

**Expected Visual Output**:
- ✅ SAE card appears **after EvidenceBand, before EfficacyCard**
- ✅ Overall Impact badge: "Overall Impact: +48%" (green chip)
- ✅ Boosting Features section:
  - Green CheckCircle icon
  - "Exon Disruption (0.88)" chip
  - "Known Hotspot (0.95)" chip
- ✅ Limiting Features section (if any):
  - Orange Warning icon
  - "Gene Essentiality (0.35)" chip
- ✅ Feature details accordion:
  - Expands to show explanation + provenance
  - "Variant significantly disrupts exon structure"
  - "Provenance: hotspot_calibration"
- ✅ RUO disclaimer at bottom

---

### **Test 5: CoPilot "Explain features?" Action** ✅ READY TO RUN

**Steps** (continue from Test 4):
1. Scroll to CoPilot Quick Actions (above SAE card)
2. Verify "Explain features?" chip appears (green, bold)
3. Click "Explain features?"

**Expected Behavior**:
- ✅ CoPilot drawer opens on right side
- ✅ Question pre-populated:
  - "Explain why confidence is boosted by these SAE features: Boosting (exon_disruption, hotspot_mutation), Limiting (essentiality_signal). What do these features mean for BRAF V600E in melanoma?"
- ✅ Context attached:
  - `sae_features` object
  - `variant` (BRAF V600E)
  - `disease` (melanoma)
- ✅ CoPilot responds with mechanistic explanation

---

## **🚨 Error Handling Tests**

### **Test 6: No SAE Features Available**

**Backend**: Return empty `sae_features` (no features extracted)

**Expected Frontend Behavior**:
- ✅ SAE card shows empty state:
  - "No SAE features available for this analysis."
  - "This may occur when variant data is incomplete or analysis profile is 'baseline'."
  - Suggestion: "Try 'Richer S' or 'Fusion' profile for enhanced feature extraction."
- ✅ "Explain features?" chip does NOT appear

---

### **Test 7: Backend Error (500)**

**Backend**: Force 500 error in SAE service

**Expected Frontend Behavior**:
- ✅ SAE card shows error state (no crash)
- ✅ Other cards (EfficacyCard, ToxicityRiskCard) still render
- ✅ User can retry analysis

---

## **📊 Smoke Test Checklist Summary**

| Test | Backend | Frontend | CoPilot | Status |
|---|---|---|---|---|
| Test 1: BRAF V600E | ✅ | ✅ | ✅ | 🔄 READY |
| Test 2: BRCA2 (Ayesha) | ✅ | ✅ | ✅ | 🔄 READY |
| Test 3: VUS (minimal) | ✅ | ✅ | N/A | 🔄 READY |
| Test 4: Card rendering | N/A | ✅ | N/A | 🔄 READY |
| Test 5: CoPilot action | N/A | ✅ | ✅ | 🔄 READY |
| Test 6: Empty state | ✅ | ✅ | N/A | 🔄 READY |
| Test 7: Error handling | ✅ | ✅ | N/A | 🔄 READY |

---

## **⚡ Quick Test Command (All-in-One)**

```bash
#!/bin/bash
# Run all SAE smoke tests

cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal

# Start backend
venv/bin/uvicorn api.main:app --host 127.0.0.1 --port 8000 > /tmp/backend.log 2>&1 &
BACKEND_PID=$!
sleep 5

echo "⚔️ Running SAE Smoke Tests..."

# Test 1: BRAF V600E
echo "Test 1: BRAF V600E (hotspot)"
curl -sS -X POST http://127.0.0.1:8000/api/clinical_genomics/analyze_variant \
  -H 'Content-Type: application/json' \
  -d '{"mutations":[{"gene":"BRAF","chrom":"7","pos":140753336,"ref":"T","alt":"A","build":"GRCh38"}],"profile":"richer","options":{"include_sae_features":true}}' \
  | jq '.efficacy.sae_features | {feature_count: .features | length, boosting_count: .boosting_features | length, overall_impact}'

# Test 2: BRCA2 (Ayesha's case)
echo "Test 2: BRCA2 (Ayesha)"
curl -sS -X POST http://127.0.0.1:8000/api/clinical_genomics/analyze_variant \
  -H 'Content-Type: application/json' \
  -d '{"mutations":[{"gene":"BRCA2","chrom":"13","pos":32936732,"ref":"C","alt":"T","build":"GRCh38"}],"profile":"richer","options":{"include_sae_features":true}}' \
  | jq '.efficacy.sae_features | {feature_count: .features | length, boosting_count: .boosting_features | length, overall_impact}'

# Test 3: VUS (minimal features)
echo "Test 3: TP53 VUS (minimal)"
curl -sS -X POST http://127.0.0.1:8000/api/clinical_genomics/analyze_variant \
  -H 'Content-Type: application/json' \
  -d '{"mutations":[{"gene":"TP53","chrom":"17","pos":7577548,"ref":"C","alt":"T","build":"GRCh38"}],"profile":"baseline","options":{"include_sae_features":true}}' \
  | jq '.efficacy.sae_features | {feature_count: .features | length, limiting_count: .limiting_features | length, overall_impact}'

# Cleanup
kill $BACKEND_PID
echo "✅ SAE smoke tests complete!"
```

**Save as**: `.cursor/ayesha/sae/testing/run_smoke_tests.sh`

---

⚔️💀 **SAE SMOKE TESTS READY - EXECUTE TO VALIDATE** 💀⚔️

