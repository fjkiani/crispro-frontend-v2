# ❓ CRITICAL EXECUTION QUESTIONS & DECISIONS

**Status:** ✅ **ALL DECISIONS APPROVED - EXECUTION READY**

**Audience:** Lower-level agents building the Food Validator

This document contains all approved decisions. DO NOT deviate from these decisions without commander approval.

---

## **📋 DECISIONS SUMMARY (FOR QUICK REFERENCE)**

| # | Question | Approved Decision | Build Priority |
|---|----------|-------------------|----------------|
| **Q1** | Evo2 Approach | Variant Scoring Proxy (with neutral S fallback) | Phase 1 |
| **Q2** | Gene Sequences | Ensembl REST + cache | Phase 1 |
| **Q3** | SAE for Supplements | Supplement-specific rules | Phase 1 |
| **Q4** | S/P/E Integration | New simplified service | Phase 2 |
| **Q5** | Target Extraction | Hardcoded first (KB), ChEMBL/LLM later | Phase 1 → Phase 3 |
| **Q6** | Treatment Line Logic | Supplement rules JSON | Phase 1 |

---

## **🔬 Q1: EVO2 API APPROACH** ✅ **DECIDED**

### **✅ APPROVED: Variant Scoring Proxy with Neutral Fallback**

**Implementation:**

```python
# Primary Path: Variant scoring at promoter regions
async def score_compound_plausibility(compound, target_gene, disease_context):
    """
    Score compound effect using promoter variant proxy.
    """
    
    # Step 1: Get gene TSS (transcription start site)
    gene_info = await fetch_ensembl_gene_info(target_gene)
    tss = gene_info['start']  # GRCh38 coordinate
    chrom = gene_info['seq_region_name']
    
    # Step 2: Create synthetic variant at promoter (TSS - 500bp)
    # Rationale: Promoter variants modulate expression (proxy for compound effect)
    promoter_pos = tss - 500
    
    # Step 3: Get reference allele at that position
    ref_allele = await fetch_reference_allele(chrom, promoter_pos)
    alt_allele = "G" if ref_allele != "G" else "A"  # Pick different allele
    
    # Step 4: Score variant
    response = await http_client.post(
        f"{API_BASE}/api/evo/score_variant_multi",
        json={
            "assembly": "GRCh38",
            "chrom": chrom,
            "pos": promoter_pos,
            "ref": ref_allele,
            "alt": alt_allele,
            "model_id": "evo2_1b"
        }
    )
    
    # Step 5: Extract delta (min_delta from response)
    result = response.json()
    delta = abs(result.get('min_delta', 0))
    
    # Step 6: Classify plausibility
    if delta > 0.5:
        plausibility = "HIGH"
    elif delta > 0.2:
        plausibility = "MODERATE"
    else:
        plausibility = "LOW"
    
    return {
        "gene": target_gene,
        "delta": delta,
        "plausibility": plausibility,
        "proxy_method": "promoter_variant",
        "coordinates": {"chrom": chrom, "pos": promoter_pos}
    }
```

**Fallback Path: Neutral S**
```python
# If binding sites unknown OR variant scoring fails:
def get_neutral_sequence_score():
    """Return neutral S score when Evo2 unavailable."""
    return {
        "sequence_score": 0.5,
        "method": "neutral_fallback",
        "reason": "Binding sites unknown - using P/E/SAE only"
    }
```

**Files to Create:**
1. `evo2_food_plausibility.py` - Implements promoter variant proxy
2. `binding_sites.json` (optional) - Override promoter default with known sites

**Agent Instructions:**
- Implement promoter variant proxy as primary
- Add neutral fallback for when TSS unavailable
- Record method in provenance
- Test with Test 6 script before using in production

**Implementation Notes (dependencies):**
- Reference allele: fetch via Ensembl REST 1-bp window at GRCh38 (ensure 1-based indexing). Persist `{chrom,pos,ref,source}` in provenance.
- `/api/evo/score_variant_multi` payload: `{ assembly: "GRCh38", chrom: "<string>", pos: <int>, ref: "A|C|G|T", alt: "A|C|G|T", model_id: "evo2_1b" }`.
- Expected response fields: prefer `min_delta`; if absent, select max abs across available deltas. Default to 0.0 on failure.

---

## **🧬 Q2: GENE SEQUENCE FETCHING** ✅ **DECIDED**

### **✅ APPROVED: Ensembl REST + In-Memory Cache**

**Implementation:**

```python
class GeneSequenceFetcher:
    """Fetch gene sequences from Ensembl with caching."""
    
    def __init__(self):
        self.cache = {}  # In-memory cache {gene: sequence}
        self.http_client = httpx.AsyncClient(timeout=30.0)
    
    async def fetch_sequence(self, gene_symbol: str) -> Optional[str]:
        """
        Fetch gene sequence from Ensembl REST API.
        
        Returns full gene sequence (may be large - 10kb to 2Mb).
        Cache to avoid repeated API calls.
        """
        
        # Check cache first
        if gene_symbol in self.cache:
            return self.cache[gene_symbol]
        
        try:
            # Step 1: Lookup gene ID
            lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json;expand=0"
            
            lookup_resp = await self.http_client.get(lookup_url, timeout=20.0)
            
            if lookup_resp.status_code != 200:
                return None
            
            gene_data = lookup_resp.json()
            gene_id = gene_data.get('id')
            
            if not gene_id:
                return None
            
            # Step 2: Get sequence
            seq_url = f"https://rest.ensembl.org/sequence/id/{gene_id}?content-type=application/json"
            seq_resp = await self.http_client.get(seq_url, timeout=20.0)
            
            if seq_resp.status_code != 200:
                return None
            
            seq_data = seq_resp.json()
            sequence = seq_data.get('seq', '')
            
            # Cache and return
            if sequence:
                self.cache[gene_symbol] = sequence
                return sequence
            
            return None
        
        except Exception as e:
            print(f"⚠️ Error fetching {gene_symbol}: {e}")
            return None
```

**Agent Instructions:**
- Use this exact pattern (proven to work in `evo.py`)
- Cache sequences in-memory (avoid repeated calls)
- Timeout: 20s per API call
- Return None on failure (graceful degradation)

---

## **⚔️ Q3: SAE FOR SUPPLEMENTS** ✅ **DECIDED**

### **✅ APPROVED: Supplement-Specific Rules (Config-Driven)**

**Files to Create:**
1. `supplement_treatment_rules.json` - Supplement rules data
2. `food_treatment_line_service.py` - Service to apply rules

**Rule Structure:**

```json
{
  "supplement_rules": {
    "NAC": {
      "mechanism": "oxidative_stress_recovery",
      "high_appropriateness_contexts": [
        "post_platinum",
        "post_chemotherapy",
        "oxidative_stress_high"
      ],
      "default_scores": {
        "line_appropriateness": 1.0,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.95
      },
      "biomarker_gates": {
        "chemotherapy_history": "platinum-based"
      }
    },
    "Vitamin D": {
      "mechanism": "dna_repair_support",
      "high_appropriateness_contexts": [
        "hrd_positive",
        "dna_repair_deficient",
        "brca1_brca2_mutated"
      ],
      "default_scores": {
        "line_appropriateness": 0.9,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.85
      },
      "biomarker_gates": {
        "HRD": "POSITIVE"
      }
    },
    "Omega-3": {
      "mechanism": "anti_inflammatory",
      "high_appropriateness_contexts": [
        "post_chemotherapy",
        "chronic_inflammation",
        "post_immunotherapy"
      ],
      "default_scores": {
        "line_appropriateness": 0.85,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.80
      }
    },
    "Curcumin": {
      "mechanism": "nfkb_inhibition",
      "high_appropriateness_contexts": [
        "chronic_inflammation",
        "nfkb_active"
      ],
      "default_scores": {
        "line_appropriateness": 0.7,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.75
      }
    }
  },
  "default_supplement": {
    "line_appropriateness": 0.6,
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.6
  }
}
```

**Service Implementation:**

```python
def compute_food_treatment_line_features(
    compound: str,
    disease_context: Dict[str, Any],
    treatment_history: Optional[Dict[str, Any]]
) -> Dict[str, float]:
    """
    Compute SAE features for dietary supplements.
    
    Returns:
        {
            "line_appropriateness": 0.9,
            "cross_resistance": 0.0,
            "sequencing_fitness": 0.85
        }
    """
    
    # Load rules
    rules = load_supplement_rules()  # From supplement_treatment_rules.json
    
    # Get rule for this compound (or default)
    compound_rule = rules.get(compound, rules['default_supplement'])
    
    # Start with default scores
    scores = compound_rule['default_scores'].copy()
    
    # Apply biomarker gates (boost if context matches)
    biomarkers = disease_context.get('biomarkers', {})
    biomarker_gates = compound_rule.get('biomarker_gates', {})
    
    for key, expected_value in biomarker_gates.items():
        if biomarkers.get(key) == expected_value:
            # Boost line appropriateness if biomarker matches
            scores['line_appropriateness'] = min(scores['line_appropriateness'] + 0.1, 1.0)
    
    # Apply treatment history gates
    if treatment_history:
        prior_therapies = treatment_history.get('prior_therapies', [])
        contexts = compound_rule.get('high_appropriateness_contexts', [])
        
        # Check if any prior therapy matches high-appropriateness context
        if 'post_platinum' in contexts and any('platin' in t.lower() for t in prior_therapies):
            scores['line_appropriateness'] = min(scores['line_appropriateness'] + 0.1, 1.0)
        
        if 'post_chemotherapy' in contexts and len(prior_therapies) > 0:
            scores['line_appropriateness'] = min(scores['line_appropriateness'] + 0.05, 1.0)
    
    return scores
```

**Agent Instructions:**
- Create `supplement_treatment_rules.json` with rules above
- Implement `compute_food_treatment_line_features()` in `food_treatment_line_service.py`
- Add biomarker + treatment history gating logic
- Return neutral defaults for unknown compounds

---

## **📊 Q4: S/P/E INTEGRATION** ✅ **DECIDED**

### **✅ APPROVED: New Simplified Service (Purpose-Built)**

**Implementation:**

```python
class FoodSPEIntegrationService:
    """
    Aggregate S/P/E + SAE for food validation.
    
    Simpler than drug efficacy_orchestrator (compound-focused, not variant-focused).
    """
    
    async def compute_spe_score(
        self,
        compound: str,
        targets: List[str],
        pathways: List[str],
        disease_context: Dict[str, Any],
        evidence_grade: str,  # From LLM: "STRONG"/"MODERATE"/"WEAK"
        treatment_history: Optional[Dict] = None,
        evo2_enabled: bool = True
    ) -> Dict[str, Any]:
        """
        Compute S/P/E + SAE score.
        
        Formula:
        - S (Sequence): Evo2 plausibility (0.4 weight) OR neutral 0.5 if disabled
        - P (Pathway): Alignment score (0.3 weight)
        - E (Evidence): Literature grade (0.3 weight)
        - SAE: Treatment line features (boosts confidence)
        
        Returns comprehensive scoring with confidence modulation.
        """
        
        # [1] SEQUENCE (S)
        if evo2_enabled:
            evo2_service = get_evo2_food_plausibility_service()
            evo2_result = await evo2_service.compute_biological_plausibility(
                compound=compound,
                targets=targets,
                disease_context=disease_context,
                pathways=pathways
            )
            sequence_score = evo2_result['overall_plausibility']
        else:
            sequence_score = 0.5  # Neutral
            evo2_result = {"method": "disabled"}
        
        # [2] PATHWAY (P)
        pathway_score = self._compute_pathway_alignment(
            compound_pathways=pathways,
            disease_pathways=disease_context.get('pathways_disrupted', [])
        )
        
        # [3] EVIDENCE (E)
        evidence_score = self._convert_evidence_grade(evidence_grade)
        
        # [4] SAE
        sae_features = None
        if treatment_history:
            sae_features = compute_food_treatment_line_features(
                compound=compound,
                disease_context=disease_context,
                treatment_history=treatment_history
            )
        
        # [5] Aggregate S/P/E
        overall_score = (
            sequence_score * 0.4 +
            pathway_score * 0.3 +
            evidence_score * 0.3
        )
        
        # [6] Confidence modulation
        confidence = self._compute_confidence(
            sequence_score=sequence_score,
            pathway_score=pathway_score,
            evidence_score=evidence_score,
            evo2_result=evo2_result if evo2_enabled else None,
            sae_features=sae_features,
            disease_context=disease_context
        )
        
        # [7] Classify verdict
        verdict = self._classify_verdict(overall_score, confidence)
        
        return {
            "overall_score": round(overall_score, 3),
            "confidence": round(confidence, 3),
            "verdict": verdict,
            "spe_breakdown": {
                "sequence": round(sequence_score, 3),
                "pathway": round(pathway_score, 3),
                "evidence": round(evidence_score, 3)
            },
            "sae_features": sae_features or {},
            "evo2_analysis": evo2_result if evo2_enabled else {"enabled": False}
        }
    
    def _compute_pathway_alignment(
        self, 
        compound_pathways: List[str],
        disease_pathways: List[str]
    ) -> float:
        """
        Score pathway alignment (compound vs disease).
        
        Simple keyword matching:
        - Aligned pathways: score 1.0 each
        - Misaligned: score 0.2 each
        - Average across all pathways
        """
        if not compound_pathways:
            return 0.5  # Neutral if unknown
        
        aligned_count = 0
        for comp_path in compound_pathways:
            for disease_path in disease_pathways:
                # Simple keyword overlap
                comp_words = set(comp_path.lower().split())
                disease_words = set(disease_path.lower().split())
                
                if comp_words & disease_words:  # Intersection
                    aligned_count += 1
                    break
        
        alignment_ratio = aligned_count / len(compound_pathways)
        
        # Score: aligned=1.0, misaligned=0.2
        score = alignment_ratio * 1.0 + (1 - alignment_ratio) * 0.2
        
        return score
    
    def _convert_evidence_grade(self, grade: str) -> float:
        """Convert literature grade to 0-1 score."""
        mapping = {
            "STRONG": 0.9,
            "MODERATE": 0.6,
            "WEAK": 0.3,
            "INSUFFICIENT": 0.1
        }
        return mapping.get(grade, 0.5)
    
    def _compute_confidence(
        self,
        sequence_score: float,
        pathway_score: float,
        evidence_score: float,
        evo2_result: Optional[Dict],
        sae_features: Optional[Dict],
        disease_context: Dict
    ) -> float:
        """
        Compute confidence with multi-stage modulation.
        
        Formula:
        - Base = (S+P+E)/3
        - Evo2 boost: +0.05 per HIGH plausibility target
        - SAE boost: (line_app + seq_fit) × 0.05
        - Biomarker boost: +0.05 if matches
        - Final = min(base + boosts, 0.95)
        """
        
        base = (sequence_score + pathway_score + evidence_score) / 3.0
        
        # Evo2 boost
        evo2_boost = 0.0
        if evo2_result and 'target_analysis' in evo2_result:
            high_plausibility = [
                t for t in evo2_result['target_analysis']
                if t.get('plausibility') == 'HIGH'
            ]
            evo2_boost = min(len(high_plausibility) * 0.05, 0.15)
        
        # SAE boost
        sae_boost = 0.0
        if sae_features:
            line_app = sae_features.get('line_appropriateness', 0)
            seq_fit = sae_features.get('sequencing_fitness', 0)
            sae_boost = (line_app + seq_fit) * 0.05
        
        # Biomarker boost
        biomarker_boost = 0.0
        biomarkers = disease_context.get('biomarkers', {})
        if biomarkers.get('HRD') == 'POSITIVE' and 'DNA repair' in str(evo2_result):
            biomarker_boost += 0.05
        if biomarkers.get('TMB', 0) >= 10:
            biomarker_boost += 0.03
        
        final = min(base + evo2_boost + sae_boost + biomarker_boost, 0.95)
        
        return final
    
    def _classify_verdict(self, score: float, confidence: float) -> str:
        """
        Classify verdict based on score AND confidence.
        
        Thresholds:
        - SUPPORTED: score ≥0.65 AND confidence ≥0.70
        - WEAK_SUPPORT: score ≥0.45 AND confidence ≥0.50
        - NOT_SUPPORTED: otherwise
        """
        
        if score >= 0.65 and confidence >= 0.70:
            return "SUPPORTED"
        elif score >= 0.45 and confidence >= 0.50:
            return "WEAK_SUPPORT"
        else:
            return "NOT_SUPPORTED"
```

**Agent Instructions:**
- Create `food_spe_integration.py` with class above
- DO NOT reuse `efficacy_orchestrator` (too complex for compounds)
- Keep formula simple: 0.4×S + 0.3×P + 0.3×E
- Apply all confidence boosts (Evo2, SAE, biomarker)

**Evidence Grade Provider & Caching:**
- Source: `llm_literature_service` with fallback grade `INSUFFICIENT` if provider unavailable.
- Cache: in-memory TTL 10 minutes per (compound, disease) key; 1 retry with exponential backoff.

**Pathway Alignment Config (optional):**
- If present, consult `disease_pathways.json` to reduce false positives; otherwise use the keyword-overlap heuristic.
- Precedence: curated config > heuristic.

---

## **🔬 Q5: TARGET EXTRACTION PRIORITY** ✅ **DECIDED**

### **✅ APPROVED: Phased Approach (KB First, ChEMBL/LLM Later)**

**Phase 1 (MVP):**
```python
def extract_targets_phase1(compound: str) -> Dict[str, Any]:
    """
    Phase 1: Use food_targets.json only.
    
    Fast, deterministic, works for 6 compounds.
    """
    
    # Load KB
    with open('food_targets.json') as f:
        food_targets = json.load(f)
    
    # Find compound
    for item in food_targets['compounds']:
        if compound.lower() in item['compound'].lower():
            return {
                "compound": item['compound'],
                "targets": item['targets'],
                "pathways": item.get('pathways', []),
                "source": "knowledge_base",
                "confidence": 0.95
            }
    
    # Not found
    return {
        "compound": compound,
        "targets": [],
        "source": "not_found",
        "confidence": 0.0,
        "error": f"Compound '{compound}' not in knowledge base"
    }
```

**Phase 3 (Enhancement):**
```python
async def extract_targets_phase3(compound: str) -> Dict[str, Any]:
    """
    Phase 3: Full extraction chain.
    
    Priority: KB → ChEMBL → LLM literature
    """
    
    # 1. Check KB (fast)
    kb_result = extract_targets_phase1(compound)
    if kb_result['targets']:
        return kb_result
    
    # 2. Try ChEMBL (medium speed)
    chembl_result = await query_chembl(compound)
    if chembl_result['targets']:
        return chembl_result
    
    # 3. Extract from literature (slow)
    llm_result = await extract_from_literature(compound)
    return llm_result
```

**Agent Instructions:**
- Build Phase 1 first (hardcoded KB)
- Test with 6 compounds (Vitamin D, Curcumin, NAC, Omega-3, Green Tea, Folate)
- Add Phase 3 AFTER Phase 1+2 work

---

## **⚔️ Q6: TREATMENT LINE LOGIC** ✅ **DECIDED**

### **✅ APPROVED: Implement Supplement Rules (Same as Q3)**

**See Q3 above for complete implementation.**

**Agent Instructions:**
- Use `supplement_treatment_rules.json` (same file as Q3)
- Function: `compute_food_treatment_line_features()`
- Apply biomarker + treatment history gates
- Return neutral defaults for unknown supplements

---

## **🎯 IMPLEMENTATION ORDER FOR AGENTS**

### **Step 1: Data Files (30 min)**
Create these JSON files FIRST:
- [ ] `.cursor/ayesha/hypothesis_validator/data/food_targets.json`
- [ ] `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json`

### **Step 2: Services (3 hours)**
Build in this order:
- [ ] `api/services/food_treatment_line_service.py` (simplest, no dependencies)
- [ ] `api/services/evo2_food_plausibility.py` (depends on Ensembl, Evo2 API)
- [ ] `api/services/compound_target_extraction.py` (depends on food_targets.json)
- [ ] `api/services/food_spe_integration.py` (depends on all above)

### **Step 3: Endpoint (1 hour)**
- [ ] `api/routers/hypothesis_validator.py` - Add `POST /api/hypothesis/validate_food_complete`

#### **API Contract: `POST /api/hypothesis/validate_food_complete`**

Request (JSON):
```json
{
  "compound": "Vitamin D",
  "disease_context": {
    "disease": "ovarian_cancer_hgs",
    "mutations": [{ "gene": "TP53", "hgvs_p": "R248Q" }],
    "biomarkers": { "HRD": "POSITIVE", "TMB": 8.2 },
    "pathways_disrupted": ["DNA repair", "Cell cycle"]
  },
  "treatment_history": {
    "current_line": "L3",
    "prior_therapies": ["carboplatin", "paclitaxel"]
  },
  "use_evo2": false,
  "force_phase": null
}
```

Response (JSON):
```json
{
  "overall_score": 0.615,
  "confidence": 0.82,
  "verdict": "WEAK_SUPPORT",
  "spe_breakdown": { "sequence": 0.45, "pathway": 0.85, "evidence": 0.60 },
  "sae_features": { "line_appropriateness": 0.9, "cross_resistance": 0.0, "sequencing_fitness": 0.85 },
  "evo2_analysis": { "enabled": false },
  "provenance": {
    "run_id": "<uuid>",
    "profile": { "use_evo2": false },
    "sources": ["food_targets.json", "llm_literature_service"],
    "timestamp": "<ISO8601>"
  }
}
```

Notes:
- All coordinates MUST be GRCh38 if/when variant proxy is used.
- `use_evo2` defaults to false until promoted from experimental.
- Response scores are 0–1. Verdict thresholds: SUPPORTED (≥0.65 and conf ≥0.70), WEAK_SUPPORT (≥0.45 and conf ≥0.50).

### **Step 4: Testing (2 hours)**
- [ ] Run all 10 PoC tests
- [ ] Validate with known compounds
- [ ] Make GO/NO-GO decision

### **Step 5: Frontend (1 hour - if GO)**
- [ ] Update `FoodValidatorAB.jsx` with new endpoint
- [ ] Add Evo2 toggle (default: OFF initially, ON after validation)
- [ ] Display S/P/E breakdown + SAE chips

---

## **✅ FINAL AGENT CHECKLIST**

**Before Building:**
- [ ] Read this entire document
- [ ] Read `testing/PROOF_OF_CONCEPT_TESTS.md`
- [ ] Understand approved decisions (Q1-Q6)

**During Build:**
- [ ] Follow implementation order (Step 1 → Step 5)
- [ ] Test each component before proceeding
- [ ] Use exact code patterns provided above
- [ ] Record all deviations in provenance
- [ ] Validate endpoint contract matches implemented router and returned schema

**After Build:**
- [ ] Run all 10 PoC tests
- [ ] Report GO/NO-GO decision
- [ ] If GO: Proceed with full build
- [ ] If NO-GO: Report failure mode + recommendation

---

## **❓ OPEN STRATEGIC QUESTIONS FOR MANAGER REVIEW**

**Status:** ⚠️ **PENDING COMMANDER DECISION** - Must resolve before finalizing execution plan

---

### **Q7: PHASED APPROACH APPROVAL** ⚠️ **CRITICAL**

**Context:**
Manager assessed Evo2 variant proxy approach as 60% confidence ("scientifically questionable but novel"), while P/E/SAE alone is 90% confidence ("definitely work").

**Manager's Recommendation:**
Build P/E/SAE first (guaranteed value), then add Evo2 as experimental toggle (potentially novel differentiation).

**Concrete Impact for Agent X:**

**If APPROVED (Phased Approach):**
- **Phase 1 (4-5 hours):** Build P/E/SAE services + endpoint with `use_evo2=false` default
- **Phase 2 (2-3 hours):** Add Evo2 toggle as experimental feature
- **Phase 3 (1-2 hours):** Validate if Evo2 provides meaningful differentiation
- **Result:** Guaranteed launch (P/E/SAE), with optional Evo2 enhancement

**If REJECTED (Original Plan):**
- **Test 6 FIRST (30 min):** Validate variant proxy works before building anything
- **Then build (7-8 hours):** Full S/P/E/SAE with Evo2 integrated from day 1
- **Result:** All-or-nothing - if Test 6 fails, back to square one

**Agent X Action:**
- [ ] **APPROVED:** Follow Steps 1→5 with `use_evo2=false` default, add toggle in Phase 2
- [ ] **REJECTED:** Run Test 6 first, build based on results

**Pending Commander Decision:** _________________

---

### **Q8: DEFINITION OF "EVO2 WORKS"**

**Context:**
Manager says: "If Evo2 works → promote to default, if not → keep experimental"

**Agent X needs to know:** What are the pass/fail criteria for promoting Evo2 from experimental → default?

**Concrete Success Criteria:**

**Option A: Technical Success (Can it run without errors?)**
- ✅ Test 6 passes: Promoter variants return non-zero deltas (delta > 0.2)
- ✅ Endpoints execute without timeouts or crashes
- ✅ Scores are stable (same input = same output ±5%)
- **Agent X Action:** If this is the criteria, just run Test 6 and verify stability

**Option B: Biological Correlation (Does it match reality?)**
- ✅ Test 10 passes: NAC (known beneficial) scores HIGH, placebo scores LOW
- ✅ Vitamin D (DNA repair) scores higher for HRD+ than HRD- patients
- ✅ False positive rate < 10% (no junk compounds scoring HIGH)
- **Agent X Action:** If this is the criteria, run Test 10 validation suite with known good/bad compounds

**Option C: Business Value (Does it justify the complexity?)**
- ✅ Evo2 S component changes ranking (not just decorative)
- ✅ Adds differentiation vs "PubMed + P/E/SAE"
- ✅ Users value the Evo2 insights
- **Agent X Action:** If this is the criteria, build A/B comparison (with/without Evo2) and measure delta

**Agent X Decision Tree:**
```
If A only → Just verify Test 6 passes, promote immediately
If B only → Run Test 10 validation, need biological proof
If C only → Need user testing and A/B comparison
If A+B → Technical + biological validation required
If A+B+C → Full validation suite before promotion
```

**Pending Commander Decision:**
- [ ] Option A only (technical bar is low)
- [ ] Option B only (biology must be correct)
- [ ] Option C only (business value required)
- [ ] Option A + B (technical + biological)
- [ ] Option A + B + C (full validation)
- [ ] Other: _________________

---

### **Q9: MVP VALUE PROPOSITION**

**Context:**
Manager says: "Still better than PubMed (integrated P/E/SAE scoring, treatment line intelligence)"

**Agent X needs to know:** Can we ship the MVP without Evo2 and still claim meaningful differentiation?

**Comparison Matrix:**

| Feature | PubMed Search | Our P/E/SAE MVP | Our P/E/SAE + Evo2 |
|---------|---------------|-----------------|---------------------|
| Literature search | ✅ | ✅ | ✅ |
| Pathway alignment | ❌ | ✅ | ✅ |
| Treatment line intelligence | ❌ | ✅ | ✅ |
| SAE confidence scoring | ❌ | ✅ | ✅ |
| Biomarker gating | ❌ | ✅ | ✅ |
| Integrated S/P/E scoring | ❌ | ✅ (P/E only, S=0.5) | ✅ (Full S/P/E) |
| Biological plausibility | ❌ | ❌ | ✅ |

**Agent X Messaging Impact:**

**If P/E/SAE alone is sufficient:**
- **Marketing:** "AI-powered supplement validation with pathway alignment and treatment line intelligence"
- **Differentiation:** Integrated P/E/SAE + biomarker gating + SAE confidence (vs. PubMed manual search)
- **Agent X Action:** Build MVP endpoint, launch with `use_evo2=false` default

**If we need Evo2:**
- **Marketing:** "Evo2-powered biological plausibility scoring for supplements"
- **Differentiation:** Sequence-level impact prediction (unique capability)
- **Agent X Action:** Must run Test 6 first, ensure Evo2 works before launch

**If hybrid (launch without Evo2, add later):**
- **Marketing Phase 1:** Focus on P/E/SAE + treatment line intelligence
- **Marketing Phase 2:** Add "New: Evo2 biological plausibility" as enhancement
- **Agent X Action:** Build MVP with toggle, promote Evo2 when validated

**Pending Commander Decision:**
- [ ] P/E/SAE alone is sufficient (ship MVP now, Evo2 optional)
- [ ] We need Evo2 for compelling differentiation (don't ship without it)
- [ ] Hybrid approach (ship MVP, add Evo2 as Phase 2 enhancement)

---

### **Q10: TESTING PRIORITY & EXECUTION ORDER**

**Context:**
Manager suggests: "Phase 3 Validation: Test if Evo2 correlates with biology (1-2 hours)"
But Test 6 (variant proxy) could be run FIRST to validate before building anything.

**Agent X needs to know:** Should I test Evo2 first or build P/E/SAE first?

**Timeline Comparison:**

**Option A: Test-First (Validate Evo2 before building):**
```
Hour 0-0.5: Run Test 6 (variant proxy)
  ├─ PASS → Hour 0.5-8.5: Build full S/P/E/SAE with Evo2
  └─ FAIL → Hour 0.5-5.5: Build P/E/SAE only

Total Time: 5.5-8.5 hours
Risk: Waste 30 min if Test 6 fails
Benefit: Don't build wrong thing
```

**Option B: Build-First (Manager's Recommendation):**
```
Hour 0-5: Build P/E/SAE MVP (guaranteed value)
Hour 5-6: Test with Vitamin D, NAC (validate P/E/SAE works)
Hour 6-8: Add Evo2 toggle (optional enhancement)
Hour 8-9: Run Test 6 + Test 10 (validate Evo2)

Total Time: 9 hours
Risk: None (P/E/SAE guaranteed to work)
Benefit: Guaranteed shippable MVP at Hour 6
```

**Option C: Hybrid (Parallel execution):**
```
Agent X: Hour 0-5: Build P/E/SAE services
Agent Y: Hour 0-0.5: Run Test 6 in parallel
Hour 5: Merge results
  ├─ Test 6 PASS → Hour 5-7: Add Evo2 to existing services
  └─ Test 6 FAIL → Hour 5-6: Ship P/E/SAE MVP

Total Time: 6-7 hours
Risk: Need 2 agents
Benefit: Best of both worlds
```

**Agent X Decision Impact:**

**If Option A (test-first):**
- **Action:** Run `PROOF_OF_CONCEPT_TESTS.md` Test 6 NOW
- **Then:** Build based on Test 6 results
- **Blocker:** Can't start building until Test 6 completes

**If Option B (build-first):**
- **Action:** Start building `food_treatment_line_service.py` NOW
- **Ignore:** Test 6 until P/E/SAE MVP complete
- **Guarantee:** Will have shippable MVP regardless of Evo2

**If Option C (hybrid):**
- **Action:** Agent X starts building, Agent Y runs Test 6
- **Merge:** Results at Hour 5

**Pending Commander Decision:**
- [ ] Option A (test-first) - Run Test 6, build based on results
- [ ] Option B (build-first) - Build P/E/SAE MVP, add Evo2 later
- [ ] Option C (hybrid) - Parallel execution with 2 agents

---

### **Q11: BUSINESS RISK TOLERANCE**

**Manager's Assessment:**
- Variant proxy: 60% confidence, "scientifically questionable but novel"
- P/E/SAE only: 90% confidence, "definitely work"

**The Question:** Is "scientifically questionable but novel" acceptable if it provides differentiation?

**Options:**
- [ ] YES - Novel differentiation is worth scientific risk (ship as "experimental/novel approach")
- [ ] NO - Must prioritize scientific rigor over novelty (only ship if scientifically validated)
- [ ] CONDITIONAL - Acceptable if we label clearly as "experimental" and have fallback

**Pending Commander Decision:** _________________

---

### **Q12: DOCTRINE UPDATE REQUIREMENT**

**Current Doctrine:** Assumes Evo2 will work and is core differentiator.

**Manager's Recommendation:** Changes this assumption (phased approach, Evo2 optional).

**The Question:** Should we update doctrine documents to reflect manager's conservative approach?

**Files to Update:**
- `EVO2_FOOD_VALIDATOR_BUILD_DOCTRINE.md` - Change from "Evo2 core" to "P/E/SAE core, Evo2 experimental"
- `EXECUTION_DECISIONS.md` - Add phased approach as approved path
- `PROOF_OF_CONCEPT_TESTS.md` - Update test priority to reflect build-first vs test-first

**Pending Commander Decision:**
- [ ] YES - Update all doctrine files to reflect phased approach
- [ ] NO - Keep doctrine as-is, let agents make tactical decisions
- [ ] PARTIAL - Update execution decisions, keep vision documents unchanged

---

### **Q13: AGENT EXECUTION PATH**

**Based on manager's recommendation, here are two execution paths:**

**Path A (Conservative - Manager's Recommendation):**
1. Build P/E/SAE services first (4-5 hours)
2. Build endpoint with Evo2 disabled (use_evo2=false default)
3. Test with known compounds (Vitamin D, NAC)
4. Then add Evo2 toggle (2-3 hours)
5. Validate Evo2 correlation (1-2 hours)

**Path B (Original - Doctrine Assumes Evo2):**
1. Run Test 6 first (30 min) - validate variant proxy
2. If Test 6 passes: Build full S/P/E/SAE with Evo2 (7-8 hours)
3. If Test 6 fails: Build P/E/SAE only (4-5 hours)

**Pending Commander Decision:**
- [ ] Path A (build P/E/SAE first, Evo2 optional)
- [ ] Path B (test Evo2 first, build based on results)
- [ ] Hybrid: Run Test 6 in parallel while building P/E/SAE, decision based on both

---

### **Q14: COMMITMENT TO ORIGINAL VISION**

**Original Doctrine:** "Evo2 Biological Plausibility" positioned as core innovation.

**Manager's Assessment:** Suggests Evo2 may not work reliably.

**The Question:** Should we pivot messaging or double down?

**Options:**
- [ ] **PIVOT:** Shift messaging to "P/E/SAE + Treatment Line Intelligence" as primary value, Evo2 as "experimental AI-powered sequence analysis"
- [ ] **DOUBLE DOWN:** Commit more R&D time to make Evo2 work (extend timeline, test more approaches)
- [ ] **HYBRID:** Keep Evo2 in vision, but de-risk with phased delivery (MVP without Evo2, Evo2 as Phase 2)

**Pending Commander Decision:** _________________

---

## **📋 MANAGER REVIEW CHECKLIST**

**Manager, please review these questions and provide decisions:**

- [ ] Q7: Phased approach approved or rejected?
- [ ] Q8: Definition of "Evo2 works" (A/B/C or combination)?
- [ ] Q9: Is P/E/SAE alone sufficient for MVP value prop?
- [ ] Q10: Testing priority (test-first vs build-first)?
- [ ] Q11: Risk tolerance (accept "scientifically questionable but novel"?)
- [ ] Q12: Update doctrine files to reflect phased approach?
- [ ] Q13: Which execution path for agents (A, B, or hybrid)?
- [ ] Q14: Pivot messaging, double down, or hybrid approach?

**Once all questions answered, agents will proceed with approved path.** ⚔️

---

**Status: ✅ COMMANDER DECISIONS FINALIZED - AGENT JR MAY PROCEED**

---

## **🎯 COMMANDER'S FINAL DECISIONS (Nov 2, 2025)**

**Zo's Strategic Analysis Completed. Commander Alpha approved phased approach.**

### **Q7: Phased Approach** ✅ **APPROVED**
- Build P/E/SAE MVP first (4-5 hours guaranteed value)
- Add Evo2 toggle as experimental feature (2-3 hours if Phase 1 validates)
- Rationale: 90% confidence on P/E/SAE vs 60% on Evo2 variant proxy

### **Q8: "Evo2 Works" Definition** ✅ **DECIDED: Option A + B**
- Technical success (delta > 0.2, stable, no crashes) AND
- Biological correlation (Test 10 passes, known good/bad compounds separate)
- Business value (Option C) deferred to post-launch user feedback

### **Q9: MVP Value Proposition** ✅ **DECIDED: Hybrid Approach**
- Launch MVP without Evo2 (P/E/SAE + treatment line intelligence sufficient)
- Add Evo2 as Phase 2 enhancement after validation
- Marketing: Focus on integrated P/E/SAE + biomarker gating + SAE confidence

### **Q10: Testing Priority** ✅ **DECIDED: Build-First (Option B)**
- Build P/E/SAE MVP first (guaranteed value)
- Test with known compounds (Vitamin D, NAC for Ayesha's case)
- Add Evo2 toggle later
- Validate Evo2 correlation as separate phase

### **Q11: Risk Tolerance** ✅ **DECIDED: CONDITIONAL**
- Acceptable if labeled "experimental" AND have P/E/SAE fallback
- Phase 1 MVP must work without Evo2
- Evo2 is enhancement, not requirement

### **Q12: Doctrine Update** ✅ **DECIDED: PARTIAL**
- Update execution decisions (this file) ✅
- Create AGENT_JR_EXECUTION_PLAN.md with phased approach ✅
- Keep vision documents (EVO2_FOOD_VALIDATOR_BUILD_DOCTRINE.md) as aspirational roadmap

### **Q13: Agent Execution Path** ✅ **DECIDED: Path A (Conservative)**
1. Build P/E/SAE services first (4-5 hours)
2. Build endpoint with `use_evo2=false` default
3. Test with Ayesha's case (Vitamin D, NAC, HRD+, L3 post-platinum)
4. Then add Evo2 toggle (2-3 hours)
5. Validate Evo2 correlation (1-2 hours)

### **Q14: Commitment to Vision** ✅ **DECIDED: HYBRID**
- Keep Evo2 in vision (future enhancement)
- De-risk with phased delivery (MVP without Evo2)
- Evo2 promoted from experimental → default only if validation passes (A+B criteria)

---

## **⚔️ AGENT JR INSTRUCTIONS**

**SEE:** `.cursor/ayesha/hypothesis_validator/AGENT_JR_EXECUTION_PLAN.md`

**Path:** Build Phase 1 P/E/SAE MVP → Test → Validate → Add Phase 2 Evo2 toggle if successful

**Timeline:** 4-5 hours Phase 1, +2-3 hours Phase 2 (optional based on validation)

**Focus:** Ayesha's ovarian cancer case (HRD+, TP53 mutant, L3 post-platinum)
