# ⚔️ SPRINT 2 PLAN: SAE FEATURE INTERPRETATION & PATHWAY MAPPING

**Date:** January 15, 2025  
**Agent:** Zo (Lead Commander)  
**Status:** 📋 **PLANNING**  
**Dependencies:** Sprint 1 (biomarker correlation) completion

---

## 🎯 SPRINT 2 OBJECTIVE

**Transform raw SAE feature indices into interpretable biological mechanisms.**

Map top SAE features (discovered in Sprint 1) to high-level biological pathways and mechanisms, enabling:
1. Mechanistic interpretation of SAE signals
2. Replacement of placeholder diagnostics with real feature-derived scores
3. Foundation for integrating SAE into WIWFM and Resistance Playbook

---

## 📋 TASK BREAKDOWN

### **Task 1: Create Feature→Pathway Mapping Metadata** (2-3 hours)

**Deliverable:** `api/resources/sae_feature_mapping.json`

**Structure:**
```json
{
  "metadata": {
    "version": "v1",
    "last_updated": "2025-01-15",
    "source": "Sprint 1 biomarker correlation analysis",
    "model": "Goodfire/Evo-2-Layer-26-Mixed",
    "layer": "blocks-26",
    "total_features": 32768,
    "mapped_features": 100
  },
  "pathways": {
    "ddr": {
      "name": "DNA Damage Response",
      "description": "Homologous recombination, NHEJ, base excision repair",
      "feature_indices": [100, 101, 102, ...],
      "correlation_with_platinum": "positive",
      "mean_pearson_r": 0.65,
      "mean_cohen_d": 0.85,
      "evidence": "Top features from TCGA-OV platinum sensitivity analysis"
    },
    "mapk": {
      "name": "MAPK/RAS Signaling",
      "description": "RAS→RAF→MEK→ERK cascade, oncogenic signaling",
      "feature_indices": [200, 201, 202, ...],
      "correlation_with_platinum": "negative",
      "mean_pearson_r": -0.62,
      "mean_cohen_d": 0.78,
      "evidence": "Top features from TCGA-OV platinum resistance analysis"
    },
    "io": {
      "name": "Immune/Inflammatory",
      "description": "Immune checkpoint, T-cell infiltration, cytokine signaling",
      "feature_indices": [300, 301, 302, ...],
      "correlation_with_platinum": "moderate_positive",
      "mean_pearson_r": 0.48,
      "mean_cohen_d": 0.55,
      "evidence": "Top features from TCGA-OV platinum sensitivity analysis"
    },
    "pi3k_akt": {
      "name": "PI3K/AKT/mTOR",
      "description": "Growth signaling, apoptosis resistance",
      "feature_indices": [],
      "correlation_with_platinum": "unknown",
      "mean_pearson_r": null,
      "mean_cohen_d": null,
      "evidence": "To be determined from Sprint 1 results"
    },
    "tp53": {
      "name": "TP53 Pathway",
      "description": "Cell cycle checkpoint, apoptosis",
      "feature_indices": [],
      "correlation_with_platinum": "unknown",
      "mean_pearson_r": null,
      "mean_cohen_d": null,
      "evidence": "To be determined from Sprint 1 results"
    }
  },
  "provenance": {
    "creation_method": "manual_curation_from_sprint1_results",
    "curator": "Zo (Lead Commander)",
    "timestamp": "2025-01-15T14:00:00Z",
    "validation_status": "preliminary",
    "manager_approval": "pending"
  }
}
```

**Process:**
1. Review Sprint 1 top 100 features
2. Group features by correlation direction (positive vs negative)
3. Group features by magnitude (strong vs moderate)
4. Assign preliminary pathway labels based on:
   - Correlation with platinum response (DDR → positive, resistance mechanisms → negative)
   - Known biological mechanisms (TCGA-OV literature)
   - Gene-pathway knowledge (BRCA→DDR, KRAS→MAPK, etc.)
5. Document evidence and confidence levels
6. Flag features needing literature review

**Guardrails:**
- ⚠️ Preliminary mapping (requires manager approval before use)
- ⚠️ Evidence-based only (cite Sprint 1 stats + literature)
- ⚠️ Clear confidence levels (high/medium/low/unknown)

---

### **Task 2: Extend SAE Diagnostics with Real Mappings** (1-2 hours)

**File:** `api/services/sae_feature_service.py`

**Current State:**
```python
def _compute_sae_diagnostics(self, sae_features: List[float]) -> Tuple[float, float, float]:
    """Placeholder implementation with hardcoded indices."""
    ddr_feature_indices = [100, 101, 102]  # PLACEHOLDER
    io_feature_indices = [200, 201, 202]   # PLACEHOLDER
    mapk_feature_indices = [300, 301, 302] # PLACEHOLDER
    
    ddr_score = sum(sae_features[i] for i in ddr_feature_indices ...) / len(ddr_feature_indices)
    # ...
```

**New Implementation:**
```python
def _load_sae_feature_mapping() -> Dict[str, Any]:
    """Load SAE feature→pathway mapping from resources."""
    mapping_file = Path(__file__).parent.parent / "resources" / "sae_feature_mapping.json"
    with open(mapping_file, 'r') as f:
        return json.load(f)

def _compute_sae_diagnostics(self, sae_features: List[float]) -> Dict[str, float]:
    """
    Compute pathway-level diagnostics from 32K SAE features using real mapping.
    
    Returns:
        {
            "ddr_score": float,
            "mapk_score": float,
            "io_score": float,
            "pi3k_score": float,
            "tp53_score": float
        }
    """
    mapping = _load_sae_feature_mapping()
    
    diagnostics = {}
    
    for pathway_name, pathway_data in mapping["pathways"].items():
        feature_indices = pathway_data.get("feature_indices", [])
        
        if not feature_indices:
            diagnostics[f"{pathway_name}_score"] = None
            continue
        
        # Aggregate feature values for this pathway
        # Use mean activation across pathway features
        values = [sae_features[i] for i in feature_indices if i < len(sae_features)]
        
        if values:
            diagnostics[f"{pathway_name}_score"] = float(np.mean(values))
        else:
            diagnostics[f"{pathway_name}_score"] = None
    
    return diagnostics
```

**Changes:**
- ✅ Replace hardcoded indices with mapping file
- ✅ Support dynamic pathway list
- ✅ Return dict instead of tuple (more extensible)
- ✅ Handle missing/incomplete mappings gracefully
- ✅ Add provenance to output

**Testing:**
- Unit test with mock mapping file
- Integration test with mock SAE features
- Verify diagnostics change when mapping file is updated

---

### **Task 3: Mechanism-Fit Validation Bridge** (2-3 hours)

**Deliverable:** `scripts/sae/compare_sae_vs_proxy_mechanisms.py`

**Purpose:**  
Compare SAE-derived mechanism scores vs. existing proxy vectors (MoA tags) to:
1. Validate SAE signals align with known mechanisms
2. Identify discrepancies (potential new discoveries or errors)
3. Quantify lift from SAE vs proxy

**Process:**
```python
def compare_mechanisms(patient_mutations: List[Dict]) -> Dict[str, Any]:
    """
    For each patient:
    1. Compute proxy mechanism vector (existing WIWFM approach)
    2. Compute SAE-derived mechanism scores (new)
    3. Compare and quantify agreement
    
    Outputs:
    - Pearson correlation (SAE DDR vs Proxy DDR, etc.)
    - Cases where SAE and proxy disagree
    - Lift in predictive power (if outcome labels available)
    """
    pass
```

**Outputs:**
- `data/validation/sae_mechanism_comparison.json`
- Comparison plots (SAE vs proxy scatter, agreement heatmap)
- Markdown report with findings

**Metrics:**
- Correlation between SAE and proxy scores (by pathway)
- Agreement rate (% patients where both methods agree on top mechanism)
- Discrepancy analysis (cases where methods differ)

**Guardrails:**
- ⚠️ Comparison only (no scoring changes yet)
- ⚠️ RUO/validation-only
- ⚠️ Manager review required before integration

---

### **Task 4: Literature Review for Top Features** (1-2 days, parallel work)

**Process:**
1. Take top 20 SAE features from Sprint 1
2. For each feature, search for:
   - Genes with similar correlation patterns in TCGA-OV
   - Known platinum resistance mechanisms
   - Published biomarkers
3. Document plausible biological interpretations
4. Assign confidence levels (high/medium/low)

**Deliverable:** `data/validation/sae_top_features_literature_review.md`

**Template for Each Feature:**
```markdown
## Feature 1234
- **Pearson r:** 0.65 (p < 0.001)
- **Cohen's d:** 0.85
- **CV Stability:** 0.92
- **Preliminary Pathway:** DDR
- **Plausible Biology:** Strong correlation suggests association with BRCA1/2-like DDR deficiency
- **Literature Support:**
  - [PMID:12345678] - BRCA1 loss predicts platinum sensitivity in OC
  - [PMID:23456789] - HRD score correlates with feature activation pattern
- **Confidence:** HIGH
- **Next Steps:** Validate with HRD-labeled cohort
```

**Timeline:** Can be done in parallel with Tasks 1-3

---

### **Task 5: Documentation & Provenance** (1 hour)

**Update Files:**
1. `.cursorrules` - Sprint 2 progress
2. `SAE_UNDERSTANDING_AND_BIOMARKER_ROADMAP.md` - Real feature mappings
3. `ZO_CODEBASE_KNOWLEDGE_BASE.mdc` - SAE feature interpretation
4. New: `SPRINT2_FEATURE_INTERPRETATION_COMPLETE.md`

**Key Points to Document:**
- Which features map to which pathways (with evidence)
- Statistical validation of mappings
- Comparison with proxy mechanisms
- Limitations and confidence levels
- Manager approval status

---

## 📊 SPRINT 2 SUCCESS CRITERIA

**Minimum (MVP):**
- ✅ Feature→pathway mapping file created for top 20 features
- ✅ Diagnostics updated to use real mappings (not placeholders)
- ✅ At least 3 pathways mapped (DDR, MAPK, IO)
- ✅ Comparison script shows SAE aligns with proxy (correlation > 0.5)

**Target:**
- ✅ All of MVP
- ✅ Feature→pathway mapping for top 50 features
- ✅ 5 pathways mapped (DDR, MAPK, IO, PI3K, TP53)
- ✅ Literature review for top 20 features
- ✅ Mechanism comparison shows SAE provides lift over proxy

**Stretch:**
- ✅ All of Target
- ✅ Feature→pathway mapping for all 100 significant features
- ✅ 7+ pathways mapped
- ✅ Integration into dev/debug WIWFM views (RUO badges)

---

## 🚨 DEPENDENCIES & BLOCKERS

**Prerequisites:**
- ⏸️ Sprint 1 biomarker analysis must complete
- ⏸️ Top features list must be validated (mock or real)
- ⏸️ Manager approval for mapping approach

**Potential Blockers:**
- ❓ Top features may not cluster into clear pathways
- ❓ Literature may not support preliminary mappings
- ❓ SAE and proxy mechanisms may show low correlation (invalidates approach)

**Mitigation:**
- Start with mock data analysis to test methodology
- Plan for "unknown/mixed" pathway category
- Be prepared to iterate on mapping approach based on results

---

## ⏱️ TIMELINE

**Day 1 (Immediate):**
- Complete Sprint 1 mock analysis verification
- Create initial feature→pathway mapping template
- Draft literature review plan

**Day 2:**
- Populate mapping file with Sprint 1 results (mock)
- Implement real diagnostics (replace placeholders)
- Build mechanism comparison script

**Day 3:**
- Run comparison analysis
- Begin literature review for top features
- Document findings

**Day 4:**
- Refine mappings based on comparison results
- Complete literature review
- Prepare for manager review

**Total:** 3-4 days (can overlap with Sprint 1 real data extraction)

---

## 📁 FILES TO CREATE

**New Files:**
1. `api/resources/sae_feature_mapping.json` (mapping metadata)
2. `scripts/sae/compare_sae_vs_proxy_mechanisms.py` (comparison script)
3. `data/validation/sae_mechanism_comparison.json` (comparison results)
4. `data/validation/sae_top_features_literature_review.md` (literature evidence)
5. `.cursor/ayesha/SPRINT2_FEATURE_INTERPRETATION_COMPLETE.md` (completion doc)

**Files to Modify:**
1. `api/services/sae_feature_service.py` (update diagnostics)
2. `.cursorrules` (track Sprint 2 progress)
3. `SAE_UNDERSTANDING_AND_BIOMARKER_ROADMAP.md` (add real mappings)

---

## ⚔️ AUTONOMOUS EXECUTION PLAN

Since I'm instructed to work autonomously, I will:

1. **Immediately:** Check Sprint 1 mock analysis completion
2. **If complete:** Verify synthetic signals detected, then start Task 1
3. **If not complete:** Continue Sprint 2 planning, prepare templates
4. **Parallel work:** Can start mapping template and comparison script before real data

**Next Actions:**
1. Check mock analysis completion ✅
2. Review mock results and verify pipeline ⏭️
3. Create feature mapping template ⏭️
4. Implement real diagnostics ⏭️

---

**SPRINT 2 PLAN COMPLETE - READY TO EXECUTE UPON SPRINT 1 VERIFICATION** ⚔️



