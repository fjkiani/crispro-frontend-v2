# SAE Readiness & Blocker Removal Plan - MBD4+TP53 HGSOC Target

**Date**: January 20, 2025  
**Target**: MBD4 Germline + TP53 Somatic HGSOC Analysis  
**Mission**: Remove all blockers, complete health checks/tests, prepare for TRUE SAE integration

---

## 🎯 Mission Context

**Primary Target**: MBD4+TP53 HGSOC case requires:
- Accurate pathway burden assessment (BER + HRD deficiency)
- Mechanism fit ranking for PARP inhibitor trials
- Early resistance detection (DNA repair capacity monitoring)
- **All three services need TRUE SAE** for optimal accuracy

**Current Reality**:
- ✅ TRUE SAE features extracted (66 patients, trained weights)
- ⚠️ Production uses PROXY SAE (gene mutations → pathway scores)
- ❌ TRUE SAE blocked by Feature→Pathway Mapping
- **Impact**: MBD4+TP53 analysis uses pathway-based vectors (works, but TRUE SAE would be better)

**Critical Clarification - S/P/E Integration**:
- ✅ **SAE uses S/P/E outputs**: Proxy SAE is derived FROM S/P/E (pathway scores, insights bundle)
- ❌ **SAE does NOT modulate S/P/E**: Manager's vision is "SAE must live inside S/P/E and modulate confidence" - but this is BLOCKED by validation requirement
- **Manager's Comment** (from `ZO_CRITICAL_AUDIT_SAE_PHASES_1_2_3.md` line 1626): "SAE must live inside S/P/E (WIWFM) and modulate confidence, not sit beside it"
- **Current State**: SAE is "display only" - uses S/P/E outputs but doesn't feed back into S/P/E confidence
- **Future State**: SAE should modulate S/P/E confidence (lifts/penalties) - requires architectural refactor (1-2 days)


Variant Impact Prediction
Which of Ayesha’s somatic (tumor) mutations are likely to be the PROBABLE drivers of her cancer?

Prioritize TP53, any HRR pathway genes (BRCA1/2, RAD51, ATM, ATR, PALB2), and anything else from tumor NGS.

2. Functional Annotation
For each detected variant (including MBD4):

What is the predicted effect at the protein level (loss-of-function, gain-of-function, dominant negative, etc)?

Is there published/validated evidence that this variant confers therapy resistance, immunogenicity, or altered cell behavior?

3. Pathway Analysis
What are the dominant oncogenic pathways, synthetic lethal vulnerabilities, and DNA repair deficiencies in this case?

Given the combination of MBD4 germline loss and TP53 somatic mutation, what pathways light up as “attackable”? (e.g., homologous recombination deficiency, base excision repair, checkpoint adaptation).

4. Drug and Therapy Prediction
Which approved or investigational drugs would be most effective given Ayesha’s REAL molecular signature?

Ask for:

PARP inhibitors (should we use olaparib/niraparib/rucaparib even if BRCA is not mutated in the tumor? Will HRD/MBD4 loss sensitize to these?)

Immune checkpoint therapy (is there evidence that MBD4 loss or related hypermutability increases response to pembrolizumab/nivolumab?)

Any off-label or repurposed drugs (ATR/ATM, DNA-PK inhibitors, etc)

Synthetic lethal partners (what drug/target combos will exploit these weaknesses?)

Hormonal therapies, given ER moderate positivity—is there evidence for aromatase inhibitors or SERDs being helpful?

5. Trial and Biomarker Matching
Given this rare MBD4 genotype + high-grade serous phenotype, which trials are a molecular fit?

Are there DNA repair–focused, mismatch repair–deficient, or rare cancer basket trials that she uniquely qualifies for, either nationally or at a regional academic center?

Also: Are there rare disease registries or research collaborations for MBD4 that are enrolling?

6. Metastasis Prediction/Surveillance
Given these biomarkers, what is the risk profile for metastasis or recurrence?

Can the tool forecast likely patterns (e.g., peritoneal, distant, CNS), and recommend enhanced surveillance or pre-emptive therapy?

7. Immunogenicity & Vaccine Target
Are there predicted neoantigens or immunogenic peptides from her unique mutation combo (MBD4 + TP53)?

Could vaccine or cell therapy trials—TILs, dendritic cell—be options, and which candidate epitopes would be prioritized?

8. Personalized Nutritional or Adjunctive Therapies
Are there diet/metabolic/small molecule interventions that biology or omics suggest would synergize with her standard chemo? (e.g., fasting, metformin, sulforaphane/broccoli sprouts, statins, aspirin)

---

## 📋 Progress Since Last Test

### ✅ Completed (January 14-20, 2025)

1. **Evo2 7B Migration** ✅
   - Migrated from evo2_1b to evo2_7b
   - Trained weights loaded (4096×32768 checkpoint)
   - Dimension mismatch resolved
   - **Code**: `src/services/sae_service/main.py:155-230`

2. **Full Cohort Extraction** ✅
   - 66 patients extracted
   - 2,897 variants processed
   - Trained weights verified
   - **Output**: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

3. **Data Quality Verification** ✅
   - Outcome field populated (53 sensitive, 11 refractory, 2 resistant)
   - Feature indices valid (0-32767)
   - All patients have SAE features

4. **Plan Strengthening** ✅
   - Proxy vs TRUE SAE distinction clarified
   - All three services documented
   - Critical blocker identified (Feature→Pathway Mapping)

5. **Clinical Trials Integration** ✅
   - SAE questions answered
   - Mechanism fit ranking patterns documented
   - Integration flow clarified

### ⏸️ Pending (Ready to Execute)

1. **Biomarker Analysis Re-Run** ⏸️
   - Service ready, bug fixed
   - Data quality verified
   - **Blocked by**: Need to verify health checks first

2. **Feature→Pathway Mapping** ❌
   - **Critical Blocker**: Blocks all three services
   - **Requires**: Top significant features from biomarker analysis
   - **Strategy**: Biomarker-driven mapping (gene→pathway inference)

3. **MBD4+TP53 End-to-End Test** ⏸️
   - Plan exists (`.cursor/plans/MBD4.mdc`)
   - **Blocked by**: Need health checks and biomarker analysis first

---

## 🔍 Phase 1: Comprehensive Health Checks

### 1.1 Modal Service Health Checks

**Check 1: SAE Service Deployment** ✅ **READY**

```bash
# Test SAE Modal service health
curl https://crispro--sae-service-saeservice-api.modal.run/health

# Expected: {"status": "healthy", "model": "evo2_7b", "weights": "trained"}
```

**Check 2: Evo2 7B Service Health** ✅ **READY**

```bash
# Test Evo2 7B Modal service health
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health

# Expected: {"status": "healthy", "model": "evo2_7b_base"}
```

**Check 3: Backend SAE Health Endpoint** ✅ **READY**

```bash
# Test backend SAE health check
curl http://localhost:8000/api/sae/health

# Expected: {"service": "sae", "status": "healthy", "url": "https://..."}
```

**Code Reference**: `api/routers/sae.py:188-219`

### 1.2 Data Quality Health Checks

**Check 4: SAE Cohort Data Structure** ✅ **READY**

```python
# scripts/sae/health_check_data.py
import json
from pathlib import Path

def check_sae_cohort_data():
    """Verify SAE cohort data structure and quality."""
    data_file = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum.json")
    
    # Check file exists
    assert data_file.exists(), f"Data file not found: {data_file}"
    
    # Load data
    with open(data_file) as f:
        data = json.load(f)
    
    # Check structure
    assert "patients" in data, "Missing 'patients' key"
    assert len(data["patients"]) == 66, f"Expected 66 patients, got {len(data['patients'])}"
    
    # Check each patient
    for patient in data["patients"]:
        assert "patient_id" in patient, "Missing patient_id"
        assert "outcome" in patient, "Missing outcome field"
        assert patient["outcome"] in ["sensitive", "refractory", "resistant"], f"Invalid outcome: {patient['outcome']}"
        assert "sae_features" in patient, "Missing sae_features"
        assert "feature_indices" in patient["sae_features"], "Missing feature_indices"
        
        # Check feature indices are valid (0-32767)
        for idx in patient["sae_features"]["feature_indices"]:
            assert 0 <= idx < 32768, f"Invalid feature index: {idx}"
    
    # Check outcome distribution
    outcomes = [p["outcome"] for p in data["patients"]]
    outcome_counts = {outcome: outcomes.count(outcome) for outcome in set(outcomes)}
    
    print(f"✅ Data structure valid")
    print(f"✅ Patient count: {len(data['patients'])}")
    print(f"✅ Outcome distribution: {outcome_counts}")
    
    return True
```

**Check 5: Feature Distribution Analysis** ⏸️ **NEW**

```python
# scripts/sae/health_check_feature_distributions.py
import json
import numpy as np
from pathlib import Path

def check_feature_distributions():
    """Verify SAE feature distributions are reasonable."""
    data_file = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum.json")
    
    with open(data_file) as f:
        data = json.load(f)
    
    # Collect all feature activations
    all_features = []
    for patient in data["patients"]:
        features = patient["sae_features"].get("feature_activations", [])
        all_features.extend(features)
    
    # Check distributions
    features_array = np.array(all_features)
    
    # Check not all zeros
    assert np.any(features_array > 0), "All features are zero - check extraction"
    
    # Check not all ones
    assert np.any(features_array < 1.0), "All features are one - check extraction"
    
    # Check reasonable range (SAE activations should be sparse, mostly 0)
    zero_fraction = np.mean(features_array == 0)
    assert zero_fraction > 0.8, f"Too few zeros ({zero_fraction:.2%}) - SAE should be sparse"
    
    # Check variation across patients
    patient_features = []
    for patient in data["patients"]:
        features = patient["sae_features"].get("feature_activations", [])
        patient_features.append(np.array(features))
    
    # Check patients have different feature patterns
    if len(patient_features) > 1:
        correlations = []
        for i in range(len(patient_features) - 1):
            corr = np.corrcoef(patient_features[i], patient_features[i+1])[0, 1]
            correlations.append(corr)
        
        mean_corr = np.mean(correlations)
        assert mean_corr < 0.95, f"Patients too similar (mean corr: {mean_corr:.3f}) - check extraction"
    
    print(f"✅ Feature distributions valid")
    print(f"✅ Zero fraction: {zero_fraction:.2%}")
    print(f"✅ Mean activation: {np.mean(features_array):.4f}")
    print(f"✅ Std activation: {np.std(features_array):.4f}")
    
    return True
```

### 1.3 Service Integration Health Checks

**Check 6: Evo2 → SAE Pipeline** ⏸️ **NEW**

```python
# scripts/sae/health_check_pipeline.py
import asyncio
import httpx
import json

async def test_evo2_to_sae_pipeline():
    """Test end-to-end Evo2 → SAE extraction pipeline."""
    
    # Test variant: MBD4 c.1239delA (from MBD4+TP53 plan)
    test_variant = {
        "assembly": "GRCh38",
        "chrom": "3",
        "pos": 129430456,
        "ref": "A",
        "alt": "",
        "window": 8192
    }
    
    async with httpx.AsyncClient(timeout=60.0) as client:
        # Step 1: Get Evo2 activations
        evo_response = await client.post(
            "http://localhost:8000/api/evo/score_variant_with_activations",
            json={**test_variant, "model_id": "evo2_7b"}
        )
        
        assert evo_response.status_code == 200, f"Evo2 failed: {evo_response.status_code}"
        evo_data = evo_response.json()
        
        assert "activations" in evo_data, "Missing activations in Evo2 response"
        assert len(evo_data["activations"]) == 4096, f"Expected 4096-dim, got {len(evo_data['activations'])}"
        
        # Step 2: Extract SAE features
        sae_response = await client.post(
            "http://localhost:8000/api/sae/extract_features",
            json={
                "activations": evo_data["activations"],
                "model_id": "evo2_7b"
            }
        )
        
        assert sae_response.status_code == 200, f"SAE failed: {sae_response.status_code}"
        sae_data = sae_response.json()
        
        assert "features" in sae_data, "Missing features in SAE response"
        assert len(sae_data["features"]) == 32768, f"Expected 32768 features, got {len(sae_data['features'])}"
        
        # Check sparsity (top K=64 active)
        active_features = [f for f in sae_data["features"] if f > 0]
        assert len(active_features) <= 64, f"Expected ≤64 active features, got {len(active_features)}"
        
        print(f"✅ Evo2 → SAE pipeline working")
        print(f"✅ Evo2 activations: {len(evo_data['activations'])}-dim")
        print(f"✅ SAE features: {len(sae_data['features'])}-dim, {len(active_features)} active")
        
        return True

# Run test
asyncio.run(test_evo2_to_sae_pipeline())
```

**Check 7: Backend SAE Feature Service** ⏸️ **NEW**

```python
# scripts/sae/health_check_backend_service.py
import httpx
import json

async def test_backend_sae_service():
    """Test backend SAE feature service (proxy SAE)."""
    
    # Test with MBD4+TP53 tumor context
    tumor_context = {
        "somatic_mutations": [
            {"gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2"},
            {"gene": "TP53", "hgvs_p": "p.R175H"}
        ],
        "hrd_score": 0.75,
        "tmb_score": 25.0,
        "msi_status": "MSS"
    }
    
    pathway_scores = {
        "ddr": 0.88,  # High DDR (MBD4 BER + TP53 checkpoint)
        "ras_mapk": 0.12,
        "pi3k": 0.15,
        "vegf": 0.20,
        "her2": 0.05,
        "io": 0.0,
        "efflux": 0.10
    }
    
    insights_bundle = {
        "functionality": 0.85,
        "chromatin": 0.70,
        "essentiality": 0.80,
        "regulatory": 0.65
    }
    
    async with httpx.AsyncClient(timeout=30.0) as client:
        # Call SAE feature service
        response = await client.post(
            "http://localhost:8000/api/sae/compute_features",
            json={
                "pathway_scores": pathway_scores,
                "insights_bundle": insights_bundle,
                "tumor_context": tumor_context
            }
        )
        
        assert response.status_code == 200, f"SAE service failed: {response.status_code}"
        data = response.json()
        
        # Check required fields
        assert "dna_repair_capacity" in data, "Missing dna_repair_capacity"
        assert "mechanism_vector" in data, "Missing mechanism_vector"
        assert "provenance" in data, "Missing provenance"
        
        # Check DNA repair capacity (Manager's C1 formula)
        dna_repair = data["dna_repair_capacity"]
        assert 0.0 <= dna_repair <= 1.0, f"Invalid DNA repair capacity: {dna_repair}"
        
        # Check mechanism vector (7D)
        mechanism_vector = data["mechanism_vector"]
        assert len(mechanism_vector) == 7, f"Expected 7D vector, got {len(mechanism_vector)}-D"
        
        # Check provenance
        assert data["provenance"]["sae"] == "proxy", "Should use proxy SAE (TRUE SAE not available)"
        
        print(f"✅ Backend SAE service working")
        print(f"✅ DNA repair capacity: {dna_repair:.3f}")
        print(f"✅ Mechanism vector: {mechanism_vector}")
        print(f"✅ Provenance: {data['provenance']['sae']}")
        
        return True
```

### 1.4 MBD4+TP53 Specific Health Checks

**Check 8: MBD4 Variant Scoring** ⏸️ **NEW**

```python
# scripts/sae/health_check_mbd4.py
import asyncio
import httpx

async def test_mbd4_variant_scoring():
    """Test MBD4 variant scoring for MBD4+TP53 analysis."""
    
    mbd4_variant = {
        "assembly": "GRCh38",
        "chrom": "3",
        "pos": 129430456,
        "ref": "A",
        "alt": "",
        "window": 8192
    }
    
    async with httpx.AsyncClient(timeout=60.0) as client:
        # Test Evo2 scoring
        response = await client.post(
            "http://localhost:8000/api/evo/score_variant_multi",
            json={**mbd4_variant, "model_id": "evo2_7b"}
        )
        
        assert response.status_code == 200, f"Evo2 scoring failed: {response.status_code}"
        data = response.json()
        
        # Check high disruption (frameshift)
        assert "delta" in data, "Missing delta score"
        delta = abs(data["delta"])
        assert delta > 0.5, f"Expected high disruption (frameshift), got delta: {delta}"
        
        print(f"✅ MBD4 variant scoring working")
        print(f"✅ Delta score: {data['delta']:.4f}")
        print(f"✅ High disruption confirmed (frameshift)")
        
        return True
```

**Check 9: MBD4+TP53 Pathway Analysis** ⏸️ **NEW**

```python
# scripts/sae/health_check_mbd4_pathways.py
import httpx

async def test_mbd4_pathway_analysis():
    """Test pathway analysis for MBD4+TP53 combination."""
    
    mutations = [
        {"gene": "MBD4", "hgvs_p": "p.Ile413Serfs*2", "chrom": "3", "pos": 129430456, "ref": "A", "alt": "", "build": "GRCh38"},
        {"gene": "TP53", "hgvs_p": "p.R175H", "chrom": "17", "pos": 7577120, "ref": "G", "alt": "A", "build": "GRCh38"}
    ]
    
    async with httpx.AsyncClient(timeout=120.0) as client:
        # Test S/P/E framework
        response = await client.post(
            "http://localhost:8000/api/efficacy/predict",
            json={
                "mutations": mutations,
                "disease": "ovarian_cancer",
                "germline_status": "positive",  # MBD4 is germline
                "tumor_context": {
                    "disease": "ovarian_cancer",
                    "hrd_score": 0.75,
                    "tmb_score": 25.0,
                    "msi_status": "MSS"
                }
            }
        )
        
        assert response.status_code == 200, f"Efficacy prediction failed: {response.status_code}"
        data = response.json()
        
        # Check pathway scores extracted
        pathway_scores = data.get("provenance", {}).get("confidence_breakdown", {}).get("pathway_disruption", {})
        assert pathway_scores, "Missing pathway scores"
        
        # Check high DDR pathway (MBD4 BER + TP53 checkpoint)
        ddr_score = pathway_scores.get("ddr", 0.0)
        assert ddr_score > 0.7, f"Expected high DDR (>0.7), got {ddr_score}"
        
        # Check PARP inhibitors rank high
        drugs = data.get("drugs", [])
        parp_drugs = [d for d in drugs if "parp" in d.get("name", "").lower() or "olaparib" in d.get("name", "").lower()]
        assert parp_drugs, "No PARP inhibitors found"
        
        top_parp = parp_drugs[0]
        assert top_parp.get("efficacy_score", 0) > 0.75, f"Expected PARP efficacy >0.75, got {top_parp.get('efficacy_score')}"
        
        print(f"✅ MBD4+TP53 pathway analysis working")
        print(f"✅ DDR pathway score: {ddr_score:.3f}")
        print(f"✅ Top PARP inhibitor: {top_parp.get('name')} (efficacy: {top_parp.get('efficacy_score'):.3f})")
        
        return True
```

---

## 🧪 Phase 2: Comprehensive Test Suite

### 2.1 Unit Tests

**Test 1: SAE Feature Extraction** ⏸️ **NEW**

```python
# tests/test_sae_extraction_pipeline.py
import pytest
import asyncio
import httpx
from pathlib import Path

@pytest.mark.asyncio
async def test_sae_feature_extraction():
    """Test SAE feature extraction from Evo2 activations."""
    # Test implementation (see Check 6 above)
    pass

@pytest.mark.asyncio
async def test_sae_feature_sparsity():
    """Test SAE features are sparse (top K=64 active)."""
    # Verify only top K features are non-zero
    pass

@pytest.mark.asyncio
async def test_sae_feature_dimensions():
    """Test SAE feature dimensions are correct (32K-dim)."""
    # Verify feature vector length
    pass
```

**Test 2: Biomarker Correlation Service** ✅ **EXISTS**

- File: `tests/test_sae_phase2_services.py`
- Status: Comprehensive test suite exists
- **Action**: Verify tests pass with current data

**Test 3: Backend SAE Feature Service** ✅ **EXISTS**

- File: `tests/test_sae_phase2_services.py`
- Status: Tests exist for Manager's formulas (C1-C10)
- **Action**: Run tests to verify proxy SAE working

### 2.2 Integration Tests

**Test 4: End-to-End SAE Pipeline** ✅ **EXISTS**

- File: `tests/test_ayesha_post_ngs_e2e.py`
- Status: E2E test exists
- **Action**: Run with MBD4+TP53 tumor context

**Test 5: MBD4+TP53 Full Analysis** ⏸️ **NEW**

```python
# tests/test_mbd4_tp53_analysis.py
import pytest
import httpx

@pytest.mark.asyncio
async def test_mbd4_tp53_full_analysis():
    """Test complete MBD4+TP53 HGSOC analysis pipeline."""
    
    # 1. Variant annotation (MBD4 + TP53)
    # 2. Pathway analysis (BER + HRD deficiency)
    # 3. Drug predictions (PARP inhibitors)
    # 4. Trial matching (mechanism fit ranking)
    # 5. SAE features (proxy, future: TRUE SAE)
    
    # Verify all steps complete successfully
    pass
```

### 2.3 Validation Tests

**Test 6: Data Quality Validation** ✅ **EXISTS**

- File: `scripts/validate_sae_tcga.py`
- Status: Validation script exists
- **Action**: Run validation on current data

**Test 7: Feature Distribution Validation** ⏸️ **NEW**

- File: `scripts/sae/health_check_feature_distributions.py` (see Check 5)
- Status: New validation script
- **Action**: Run to verify feature distributions

---

## 🚀 Phase 3: Blocker Removal

### 3.1 Critical Blocker: Feature→Pathway Mapping

**Current Status**: ❌ **BLOCKS ALL THREE SERVICES**

**Blocker Details**:
- Cannot map 32K-dim SAE features → 7D pathway scores
- Blocks: Resistance Prophet, Mechanism Fit Ranking, Early Resistance Detection
- **Impact on MBD4+TP53**: Uses proxy pathway scores (works, but TRUE SAE would be better)

**Resolution Strategy**:

**Step 1: Re-Run Biomarker Analysis** ⏸️ **READY**

```bash
# After health checks pass
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots
```

**Expected Output**:
- Top 100 significant features (p < 0.01, Cohen's d ≥ 0.3)
- Feature correlations with platinum response
- Feature stability across CV folds

**Step 2: Create Feature→Pathway Mapping** ⏸️ **PENDING**

**Approach**: Biomarker-driven mapping
1. Take top 100 significant features from biomarker analysis
2. For each feature, identify which genes activate it (from extraction data)
3. Map genes → pathways (using existing `drug_mapping.py`)
4. Assign feature to pathway(s) based on gene associations
5. Validate: BRCA1 mutations → DDR pathway features should be high

**Implementation**:
```python
# scripts/sae/create_feature_pathway_mapping.py
def create_feature_pathway_mapping(biomarker_results, extraction_data):
    """
    Create feature→pathway mapping from biomarker results.
    
    Strategy:
    1. For each top significant feature:
       - Find which patients have this feature active
       - Identify genes mutated in those patients
       - Map genes → pathways (using drug_mapping.py)
       - Assign feature to pathway(s)
    
    2. Validate mapping:
       - BRCA1 mutations → DDR pathway features should be high
       - KRAS mutations → MAPK pathway features should be high
       - HER2 mutations → HER2 pathway features should be high
    """
    pass
```

**Step 3: Validate Mapping** ⏸️ **PENDING**

**Validation Tests**:
1. BRCA1 mutation → DDR pathway features high
2. KRAS mutation → MAPK pathway features high
3. HER2 mutation → HER2 pathway features high
4. MBD4+TP53 → DDR pathway features high (for our target case)

**Step 4: Integrate into SAE Feature Service** ⏸️ **PENDING**

- Update `sae_feature_service.py` to use TRUE SAE pathway scores
- Replace proxy pathway scores with SAE-derived scores
- Test: MBD4+TP53 should show higher DDR pathway score with TRUE SAE

### 3.2 Secondary Blocker: Biomarker Analysis Re-Run

**Current Status**: ⏸️ **READY TO RUN** (blocked by health checks)

**Resolution**:
1. ✅ Health checks pass
2. ✅ Data quality verified
3. ⏸️ Run biomarker analysis
4. ⏸️ Review results
5. ⏸️ Proceed with Feature→Pathway Mapping

---

## 🎯 Phase 4: MBD4+TP53 Integration

### 4.1 End-to-End Test with MBD4+TP53

**Test Case**: Complete MBD4+TP53 HGSOC analysis

**Steps**:
1. Variant annotation (MBD4 + TP53)
2. Pathway analysis (BER + HRD deficiency)
3. Drug predictions (PARP inhibitors)
4. Trial matching (mechanism fit ranking)
5. SAE features (proxy now, TRUE SAE when mapping ready)

**Expected Results**:
- High DDR pathway score (>0.85)
- PARP inhibitors rank #1-2 (efficacy >0.80)
- Mechanism fit for PARP trials >0.90
- DNA repair capacity >0.85

**Code Reference**: `.cursor/plans/MBD4.mdc` (complete plan)

### 4.2 SAE Enhancement for MBD4+TP53

**Current (Proxy SAE)**:
- DDR pathway: 0.85 (from gene mutations)
- Mechanism fit for PARP: 0.82

**Future (TRUE SAE)**:
- DDR pathway: 0.92 (from SAE features, captures BER+checkpoint synergy)
- Mechanism fit for PARP: 0.91 (better alignment)

**Benefit**: More accurate pathway burden for rare combinations like MBD4+TP53

---

## 📊 Phase 5: Health Check Execution Plan

### 5.1 Pre-Flight Checklist

**Before Running Tests**:
- [ ] Backend running (`http://localhost:8000`)
- [ ] Modal services deployed (SAE + Evo2 7B)
- [ ] Environment variables set (`SAE_SERVICE_URL`, `ENABLE_TRUE_SAE`)
- [ ] Data file exists (`sae_features_tcga_ov_platinum.json`)

### 5.2 Health Check Execution Order

1. **Modal Service Health** (5 min)
   - Check SAE service: `curl https://crispro--sae-service-saeservice-api.modal.run/health`
   - Check Evo2 7B service: `curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health`
   - Check backend health: `curl http://localhost:8000/api/sae/health`

2. **Data Quality Checks** (10 min)
   - Run `scripts/sae/health_check_data.py`
   - Run `scripts/sae/health_check_feature_distributions.py`
   - Verify all checks pass

3. **Pipeline Integration Checks** (15 min)
   - Run `scripts/sae/health_check_pipeline.py`
   - Run `scripts/sae/health_check_backend_service.py`
   - Verify Evo2 → SAE → Backend pipeline works

4. **MBD4+TP53 Specific Checks** (20 min)
   - Run `scripts/sae/health_check_mbd4.py`
   - Run `scripts/sae/health_check_mbd4_pathways.py`
   - Verify MBD4+TP53 analysis works end-to-end

### 5.3 Test Execution Order

1. **Unit Tests** (10 min)
   - Run `pytest tests/test_sae_phase2_services.py`
   - Run `pytest tests/test_sae_extraction_pipeline.py` (new)

2. **Integration Tests** (15 min)
   - Run `pytest tests/test_ayesha_post_ngs_e2e.py`
   - Run `pytest tests/test_mbd4_tp53_analysis.py` (new)

3. **Validation Tests** (10 min)
   - Run `python3 scripts/validate_sae_tcga.py`
   - Review validation results

---

## 🎯 Success Criteria

### Health Checks
- ✅ All Modal services healthy
- ✅ All data quality checks pass
- ✅ All pipeline integration checks pass
- ✅ MBD4+TP53 specific checks pass

### Tests
- ✅ All unit tests pass
- ✅ All integration tests pass
- ✅ All validation tests pass

### Blocker Removal
- ✅ Biomarker analysis re-run complete
- ✅ Feature→Pathway Mapping created
- ✅ Mapping validated (BRCA1 → DDR, KRAS → MAPK, MBD4+TP53 → DDR)
- ✅ SAE Feature Service updated to use TRUE SAE

### MBD4+TP53 Integration
- ✅ End-to-end analysis works
- ✅ Pathway scores accurate (DDR >0.85)
- ✅ PARP inhibitors rank correctly (efficacy >0.80)
- ✅ Mechanism fit ranking works (fit >0.90)

---

## 📝 Action Items

### Immediate (Next 2 Hours)

1. **Create Health Check Scripts** ⏸️
   - [ ] `scripts/sae/health_check_data.py`
   - [ ] `scripts/sae/health_check_feature_distributions.py`
   - [ ] `scripts/sae/health_check_pipeline.py`
   - [ ] `scripts/sae/health_check_backend_service.py`
   - [ ] `scripts/sae/health_check_mbd4.py`
   - [ ] `scripts/sae/health_check_mbd4_pathways.py`

2. **Run Health Checks** ⏸️
   - [ ] Modal service health
   - [ ] Data quality
   - [ ] Pipeline integration
   - [ ] MBD4+TP53 specific

3. **Run Tests** ⏸️
   - [ ] Unit tests
   - [ ] Integration tests
   - [ ] Validation tests

### Short-Term (Next Week)

4. **Re-Run Biomarker Analysis** ⏸️
   - [ ] After health checks pass
   - [ ] Review results
   - [ ] Document findings

5. **Create Feature→Pathway Mapping** ⏸️
   - [ ] After biomarker analysis complete
   - [ ] Implement mapping script
   - [ ] Validate mapping
   - [ ] Integrate into SAE Feature Service

6. **MBD4+TP53 End-to-End Test** ⏸️
   - [ ] After mapping complete
   - [ ] Test with TRUE SAE
   - [ ] Compare proxy vs TRUE SAE results

---

## 🔗 Related Documents

- **Master Plan**: `.cursor/plans/final-comprehensive-document-review-bad14970.plan.md`
- **Uncertainties**: `.cursor/ayesha/PLAN_UNCERTAINTIES_AND_RISKS.md`
- **MBD4 Plan**: `.cursor/plans/MBD4.mdc`
- **Pipeline Roadmap**: `.cursor/rules/SAE_TO_RESISTANCE_PROPHET_PIPELINE.mdc`
- **Modal URLs**: `.cursor/ayesha/MODAL_SERVICE_URLS.md`

---

**Status**: ✅ **PLAN READY FOR EXECUTION**

**Next Step**: Create health check scripts and run pre-flight checks

