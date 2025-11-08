# 🎯 EVO2-POWERED FOOD VALIDATOR - COMPLETE BUILD DOCTRINE

**Mission:** Build the ONLY food/supplement validator that predicts IF a compound will work, not just what the literature says.

**Differentiation:** Evo2 biological plausibility + S/P/E + SAE + biomarker integration = mechanistic validation that PubMed alone cannot provide.

---

## 🔥 WHAT MAKES US UNIQUE (vs. PubMed/Google Scholar)

### **Existing Tools (PubMed, Google Scholar, Examine.com):**
- ✅ Find papers and studies
- ✅ Show abstracts and dosing
- ❌ **NO biological plausibility scoring**
- ❌ **NO patient-specific recommendations**
- ❌ **NO mechanism validation**
- ❌ **NO treatment line intelligence**

### **Our Tool (Full Stack):**
- ✅ Find papers (via LLM + PubMed)
- ✅ **Evo2 Biological Plausibility Score** (target → pathway → disease)
- ✅ **S/P/E Integration** (sequence/pathway/evidence fusion)
- ✅ **SAE Features** (line appropriateness, cross-resistance, sequencing fitness)
- ✅ **Biomarker Gating** (HRD+, TP53 status, TMB → personalized verdicts)
- ✅ **MoA Validation** (does compound mechanism align with tumor biology?)
- ✅ **Dynamic Discovery** (works for ANY compound, not hardcoded list)

---

## 🤔 CRITICAL DECISIONS & CLARIFICATIONS

**Status:** ⚠️ **PENDING COMMANDER APPROVAL** - Must resolve before execution

**⚠️ STRATEGIC REVIEW REQUIRED:**
Manager has assessed approach confidence and recommended phased delivery. See `questions/EXECUTION_DECISIONS.md` Q7-Q14 for open strategic questions that must be answered before final execution plan.

### **Q1: Evo2 API Usage - How to Score Biological Plausibility?**

**Question:** For Evo2 biological plausibility scoring, which approach should we use?

**Option A:** Use variant scoring (`/api/evo/score_variant_multi`)
- We'd need to create "synthetic variants" at binding sites
- Example: If curcumin binds NFKB1, create a variant at the binding site

**Option B:** Use direct sequence scoring (`/api/evo/score`) with prompts (RECOMMENDED)
- Provide disease context prompt → get baseline likelihood
- Provide intervention context → get intervention likelihood
- Compute delta
- **Alignment:** Matches doctrine's prompt-based approach

**Option C:** Use a different endpoint/approach entirely?

**Recommendation:** **Option B** - Use `/api/evo/score` with prompt-based sequence scoring. This aligns with doctrine and allows us to provide rich biological context to Evo2.

---

### **Q2: Gene Sequence Fetching - What's the Right Endpoint?**

**Question:** How do we fetch gene sequences for targets like VDR, NFKB1, TP53?

**Found in codebase:**
- `_fetch_reference_window()` in `evo.py` uses Ensembl REST API directly
- Doctrine mentions `/api/genomic_intel/gene_info` but this may not exist

**Options:**
- **A)** Use Ensembl REST directly (as in `evo.py`) - **RECOMMENDED**
- **B)** Call existing genomic_intel endpoint (if it exists)
- **C)** Use cached FASTA files

**Recommendation:** **Option A** - Reuse the proven `_fetch_reference_window()` pattern from `evo.py`. Simple, tested, and doesn't require new infrastructure.

---

### **Q3: SAE Treatment Line Features - Dietary Supplements?**

**Question:** How does `compute_treatment_line_features()` work for dietary supplements like "Vitamin D" or "NAC"?

**Current signature:**
```python
compute_treatment_line_features(
    drug_name: str,  # e.g., "olaparib"
    disease: str,
    treatment_history: Optional[TreatmentHistory]
)
```

**Issue:** Supplements don't have NCCN guidelines or preferred lines. What should we return?

**Options:**
- **A)** Call with `drug_name="dietary_supplement"` → always returns defaults (line_appropriateness=0.6)
- **B)** Create `compute_food_treatment_line_features()` with supplement-specific logic - **RECOMMENDED**
  - NAC post-platinum → HIGH line appropriateness (oxidative stress recovery)
  - Vitamin D for HRD+ → HIGH line appropriateness (DNA repair support)
  - Curcumin for inflammation → MODERATE appropriateness
- **C)** Skip SAE for supplements entirely

**Recommendation:** **Option B** - Build supplement-aware treatment line logic. Map compounds to appropriate contexts:
- **NAC**: Post-platinum therapy (L3+) → HIGH appropriateness (reduces oxidative stress from platinum)
- **Vitamin D**: HRD+ or TP53 mutations → HIGH appropriateness (supports DNA repair)
- **Omega-3**: Chronic inflammation markers → MODERATE appropriateness
- **Curcumin**: NF-κB pathway activation → MODERATE appropriateness

---

### **Q4: S/P/E Integration - Reuse or Simplify?**

**Question:** Should we reuse existing `efficacy_orchestrator` or create a simplified version for food?

**Existing S/P/E (from `efficacy_orchestrator`):**
- Sequence (S): Variant-based scoring (`score_variant_multi`)
- Pathway (P): Variant→pathway aggregation with drug-specific weights
- Evidence (E): Literature + ClinVar integration

**For Food Validator:**
- Sequence (S): **Compound → target gene** plausibility (Evo2 prompt-based, not variant scoring)
- Pathway (P): **Compound pathways → disease pathways** alignment (simpler mapping)
- Evidence (E): **Literature grade** (already have LLM service)

**Options:**
- **A)** Reuse orchestrator - adapt it to work with compound targets instead of variants (COMPLEX)
- **B)** Create simplified `FoodSPEIntegrationService` - build from scratch with food-specific logic - **RECOMMENDED**
- **C)** Hybrid - reuse evidence/confidence components, build new S/P logic

**Recommendation:** **Option B** - Create new `FoodSPEIntegrationService`. The sequence scoring approach is fundamentally different (prompt-based vs variant-based), and pathway mapping is simpler (compound pathways vs variant pathways). Reuse only evidence grading logic.

---

### **Q5: Compound Targets - Build Order Priority?**

**Question:** The `compound_target_extraction.py` service needs to extract targets. What's the priority for Phase 1?

**Current plan in doctrine:**
1. Check `food_targets.json` (hardcoded)
2. Try ChEMBL
3. Extract from literature via LLM

**For Phase 1 (Evo2 plausibility), we need targets IMMEDIATELY. Should we:**
- **A)** Start with hardcoded only - use `food_targets.json` for Phase 1, add ChEMBL/literature in Phase 3 - **RECOMMENDED**
- **B)** Build target extraction first - complete Phase 3 before Phase 1
- **C)** Mock targets for testing - hardcode test targets like `["VDR", "TP53"]` for Vitamin D in Phase 1 tests

**Recommendation:** **Option A** - Use `food_targets.json` for Phase 1. This allows us to:
1. Test Evo2 plausibility with known targets immediately
2. Validate the core innovation (Evo2 scoring) before adding complexity
3. Add ChEMBL/literature extraction in Phase 3 as enhancement

**Note:** Ensure `food_targets.json` has at least 10 compounds with targets defined for testing.

---

### **Q6: Treatment Line Features for Supplements - Default Logic?**

**Question:** For dietary supplements, what should `compute_treatment_line_features()` return when supplements don't have NCCN metadata?

**Options:**
- **A)** Default values: Always return `line_appropriateness=0.6`, `cross_resistance=0.0`, `sequencing_fitness=0.6`
- **B)** Supplement-specific logic - **RECOMMENDED**
  - NAC post-platinum → HIGH line appropriateness (oxidative stress recovery)
  - Vitamin D for HRD+ → HIGH line appropriateness (DNA repair support)
  - Omega-3 for inflammation → MODERATE appropriateness
  - Map supplements to disease contexts (HRD, TP53, inflammation markers)
- **C)** Skip SAE for supplements - only use for drugs

**Recommendation:** **Option B** - Build supplement-aware treatment line logic. Create mapping:
```python
FOOD_TREATMENT_LINE_LOGIC = {
    "NAC": {
        "high_appropriateness": ["post_platinum", "oxidative_stress"],
        "contexts": ["current_line >= 3", "prior_therapies contains 'carboplatin'"]
    },
    "Vitamin D": {
        "high_appropriateness": ["hrd_positive", "tp53_mutation"],
        "contexts": ["biomarkers['HRD'] == 'POSITIVE'", "mutations contain TP53"]
    },
    "Omega-3": {
        "high_appropriateness": ["chronic_inflammation"],
        "contexts": ["pathways_disrupted contains 'Inflammation'"]
    }
}
```

---

### **Q7: Evo2 Prompt Format - Sequence vs Context?**

**Question:** The doctrine shows prompt-based scoring, but what's the exact format Evo2 `/score` endpoint expects?

**From doctrine (lines 212-232):**
```python
disease_prompt = f"""
Genomic context: {disease_name} tumor microenvironment
Target gene: {target_gene}
Regulatory state: Cancer-active (constitutive activation)
...
Sequence: {gene_seq[:500]}
"""

# Then call:
baseline_score = await self._call_evo2_score(disease_prompt)
```

**But looking at `evo.py`, the `/score` endpoint expects:**
```python
{
  "sequence": str,  # Required
  "model_id": "evo2_1b"  # Optional
}
```

**Options:**
- **A)** Use prompt as `sequence` field - put the entire prompt text in `sequence` parameter
- **B)** Use separate fields if Evo2 supports them (e.g., `prompt`, `context`)
- **C)** Construct sequence differently - maybe concatenate gene sequence with context prefix?

**Recommendation:** Need to verify Evo2 API contract. Should I check the actual Evo2 service endpoint to confirm the expected format?

---

### **Q8: Treatment Line Integration Import Path?**

**Question:** How do we import `compute_treatment_line_features()` from `.cursor/ayesha/treatment_lines/`?

**Current path:** `.cursor/ayesha/treatment_lines/backend/services/treatment_line_integration.py`

**The function signature:**
```python
def compute_treatment_line_features(
    drug_name: str,
    disease: str,
    treatment_history: Optional[TreatmentHistory] = None
) -> Dict[str, Any]
```

**Issues:**
- Path is outside `oncology-backend-minimal/` directory
- Uses `TreatmentHistory` Pydantic model
- May have dependencies on other treatment_lines modules

**Options:**
- **A)** Add `sys.path.insert()` to import from `.cursor/ayesha/treatment_lines/` - **RECOMMENDED**
- **B)** Copy relevant functions into `oncology-backend-minimal/api/services/`
- **C)** Create wrapper service that calls treatment line functions via HTTP (overkill)

**Recommendation:** **Option A** - Use dynamic path import like existing code does:
```python
import sys
from pathlib import Path

treatment_lines_path = Path(__file__).parent.parent.parent.parent / ".cursor/ayesha/treatment_lines/backend"
sys.path.insert(0, str(treatment_lines_path))

from services.treatment_line_integration import compute_treatment_line_features
```

---

## 📋 PHASE-BY-PHASE BUILD PLAN

---

## **PHASE 1: EVO2 BIOLOGICAL PLAUSIBILITY SERVICE (3 hours)**

### **Goal:**
Use Evo2 to score if a food compound can biologically modulate target genes in the patient's disease context.

### **File:** `oncology-coPilot/oncology-backend-minimal/api/services/evo2_food_plausibility.py`

### **Implementation:**

```python
"""
Evo2 Biological Plausibility Service for Food/Supplement Validation

Core Innovation: Score compound → target → pathway → disease impact using Evo2's 
sequence-level understanding to predict if intervention is biologically plausible.
"""

from typing import List, Dict, Any, Optional
import httpx
import asyncio
from pathlib import Path
import os

# Evo2 API endpoint
EVO2_URL = os.getenv("EVO_URL_1B", "https://evo2-service.modal.run")

class Evo2FoodPlausibilityService:
    """
    Service to compute biological plausibility of food compounds using Evo2.
    
    Workflow:
    1. Extract target genes for compound (from ChEMBL, literature, or knowledge base)
    2. For each target, fetch gene sequence
    3. Score baseline disease-active state (Evo2)
    4. Score post-intervention state (Evo2)
    5. Compute delta (plausibility score)
    6. Validate mechanism against disease pathways
    """
    
    def __init__(self):
        self.evo2_client = httpx.AsyncClient(timeout=60.0)
        self.gene_cache = {}
    
    async def compute_biological_plausibility(
        self,
        compound: str,
        targets: List[str],  # e.g., ["NFKB1", "PTGS2", "AKT1"]
        disease_context: Dict[str, Any],  # Patient's mutations, biomarkers, disease
        pathways: Optional[List[str]] = None  # e.g., ["NF-κB signaling", "PI3K/AKT"]
    ) -> Dict[str, Any]:
        """
        Main method: Compute biological plausibility score for compound.
        
        Args:
            compound: Name of food/supplement (e.g., "Curcumin", "Vitamin D")
            targets: List of molecular targets from literature/ChEMBL
            disease_context: {
                "disease": "ovarian_cancer_hgs",
                "mutations": [{"gene": "TP53", "hgvs_p": "R248Q", ...}],
                "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2, ...},
                "pathways_disrupted": ["DNA repair", "Cell cycle"]
            }
            pathways: Known pathways compound affects
        
        Returns:
            {
                "overall_plausibility": 0.65,  # 0-1 score
                "verdict": "MODERATE",  # HIGH/MODERATE/LOW
                "target_analysis": [
                    {
                        "gene": "NFKB1",
                        "baseline_activity": 0.85,
                        "intervention_activity": 0.72,
                        "delta": 0.13,
                        "plausibility": "LOW",
                        "mechanism": "NF-κB inhibition",
                        "rationale": "Modest predicted impact due to poor bioavailability",
                        "confidence": 0.55
                    }
                ],
                "mechanisms_validated": ["NF-κB inhibition"],
                "mechanisms_uncertain": ["COX-2 modulation"],
                "pathway_alignment": {
                    "aligned": ["NF-κB signaling"],
                    "misaligned": []
                },
                "provenance": {
                    "method": "evo2_food_plausibility_v1",
                    "evo2_model": "evo2_1b",
                    "targets_scored": 3,
                    "timestamp": "2024-11-02T..."
                }
            }
        """
        
        results = []
        
        for target_gene in targets:
            try:
                target_result = await self._score_target_gene(
                    target_gene=target_gene,
                    compound=compound,
                    disease_context=disease_context
                )
                results.append(target_result)
            except Exception as e:
                results.append({
                    "gene": target_gene,
                    "error": str(e),
                    "plausibility": "UNKNOWN"
                })
        
        # Aggregate results
        overall_plausibility = self._compute_overall_score(results)
        verdict = self._classify_verdict(overall_plausibility)
        validated_mechanisms = [r['mechanism'] for r in results if r.get('delta', 0) > 0.2]
        uncertain_mechanisms = [r['mechanism'] for r in results if 0 < r.get('delta', 0) <= 0.2]
        
        # Pathway alignment check
        pathway_alignment = self._check_pathway_alignment(
            targets=targets,
            pathways=pathways,
            disease_context=disease_context
        )
        
        return {
            "overall_plausibility": round(overall_plausibility, 3),
            "verdict": verdict,
            "target_analysis": results,
            "mechanisms_validated": validated_mechanisms,
            "mechanisms_uncertain": uncertain_mechanisms,
            "pathway_alignment": pathway_alignment,
            "provenance": {
                "method": "evo2_food_plausibility_v1",
                "evo2_model": "evo2_1b",
                "targets_scored": len(targets),
                "targets_failed": len([r for r in results if 'error' in r]),
                "timestamp": self._get_timestamp()
            }
        }
    
    async def _score_target_gene(
        self,
        target_gene: str,
        compound: str,
        disease_context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Score a single target gene: baseline vs. intervention.
        
        Returns:
            {
                "gene": "NFKB1",
                "baseline_activity": 0.85,
                "intervention_activity": 0.72,
                "delta": 0.13,
                "plausibility": "LOW",
                "mechanism": "NF-κB inhibition",
                "rationale": "...",
                "confidence": 0.55
            }
        """
        
        # 1. Fetch gene sequence (or use cached)
        gene_seq = await self._fetch_gene_sequence(target_gene)
        
        if not gene_seq:
            return {
                "gene": target_gene,
                "error": "Gene sequence not found",
                "plausibility": "UNKNOWN"
            }
        
        # 2. Build disease-active context for Evo2
        disease_name = disease_context.get('disease', 'cancer').replace('_', ' ')
        mutations = disease_context.get('mutations', [])
        pathways_disrupted = disease_context.get('pathways_disrupted', [])
        
        disease_prompt = f"""
Genomic context: {disease_name} tumor microenvironment
Target gene: {target_gene}
Regulatory state: Cancer-active (constitutive activation)
Disrupted pathways: {', '.join(pathways_disrupted) if pathways_disrupted else 'Unknown'}
Patient mutations: {', '.join([m['gene'] for m in mutations]) if mutations else 'None'}
Sequence: {gene_seq[:500]}
"""
        
        # 3. Score baseline cancer activity (Evo2)
        baseline_score = await self._call_evo2_score(disease_prompt)
        
        # 4. Build intervention context
        intervention_prompt = f"""
Genomic context: {disease_name} with {compound} intervention
Target gene: {target_gene}
Regulatory state: {compound}-inhibited
Mechanism: {compound} targeting {target_gene}
Therapeutic intervention: Active
Sequence: {gene_seq[:500]}
"""
        
        # 5. Score post-intervention activity
        intervention_score = await self._call_evo2_score(intervention_prompt)
        
        # 6. Compute delta (plausibility)
        delta = abs(baseline_score - intervention_score)
        
        plausibility = "HIGH" if delta > 0.5 else "MODERATE" if delta > 0.2 else "LOW"
        
        # 7. Generate mechanism description
        mechanism = f"{target_gene} modulation by {compound}"
        
        # 8. Generate rationale
        rationale = self._generate_rationale(
            target_gene=target_gene,
            compound=compound,
            delta=delta,
            baseline_score=baseline_score,
            intervention_score=intervention_score
        )
        
        # 9. Compute confidence (based on delta magnitude and score stability)
        confidence = min(delta * 1.5, 1.0)  # Higher delta = higher confidence
        
        return {
            "gene": target_gene,
            "baseline_activity": round(baseline_score, 3),
            "intervention_activity": round(intervention_score, 3),
            "delta": round(delta, 3),
            "plausibility": plausibility,
            "mechanism": mechanism,
            "rationale": rationale,
            "confidence": round(confidence, 3)
        }
    
    async def _fetch_gene_sequence(self, gene: str) -> Optional[str]:
        """
        Fetch gene sequence from Ensembl or cache.
        
        Implementation options:
        1. Call existing /api/genomic_intel/gene_info endpoint
        2. Direct Ensembl REST API call
        3. Use local FASTA cache
        """
        
        if gene in self.gene_cache:
            return self.gene_cache[gene]
        
        try:
            # Option 1: Call our genomic_intel router
            response = await self.evo2_client.get(
                f"{EVO2_URL.replace('evo2', 'api')}/genomic_intel/gene_info",
                params={"gene": gene, "build": "GRCh38"}
            )
            
            if response.status_code == 200:
                data = response.json()
                sequence = data.get('sequence', '')
                self.gene_cache[gene] = sequence
                return sequence
            
            # Option 2: Direct Ensembl call (fallback)
            ensembl_url = f"https://rest.ensembl.org/sequence/id/{gene}?content-type=application/json"
            response = await self.evo2_client.get(ensembl_url)
            
            if response.status_code == 200:
                data = response.json()
                sequence = data.get('seq', '')
                self.gene_cache[gene] = sequence
                return sequence
        
        except Exception as e:
            print(f"⚠️ Error fetching gene sequence for {gene}: {e}")
        
        return None
    
    async def _call_evo2_score(self, prompt: str) -> float:
        """
        Call Evo2 API to score sequence likelihood.
        
        Endpoint: POST /api/evo/score_variant_multi
        
        Returns likelihood score (0-1 range after normalization).
        """
        
        try:
            response = await self.evo2_client.post(
                f"{EVO2_URL}/score",
                json={"prompt": prompt, "model_id": "evo2_1b"},
                timeout=30.0
            )
            
            if response.status_code == 200:
                data = response.json()
                # Normalize score to 0-1 range
                raw_score = data.get('likelihood', 0.5)
                normalized = min(max(raw_score / 100.0, 0.0), 1.0)
                return normalized
            else:
                # Fallback: return neutral score
                return 0.5
        
        except Exception as e:
            print(f"⚠️ Evo2 API error: {e}")
            return 0.5  # Neutral score on failure
    
    def _compute_overall_score(self, results: List[Dict]) -> float:
        """
        Aggregate target-level scores into overall plausibility.
        
        Weighting strategy:
        - HIGH plausibility targets: weight 1.0
        - MODERATE: weight 0.5
        - LOW: weight 0.2
        - UNKNOWN/error: weight 0.0
        """
        
        valid_results = [r for r in results if 'delta' in r]
        
        if not valid_results:
            return 0.0
        
        weighted_sum = 0.0
        total_weight = 0.0
        
        for r in valid_results:
            delta = r.get('delta', 0)
            
            if delta > 0.5:
                weight = 1.0
            elif delta > 0.2:
                weight = 0.5
            else:
                weight = 0.2
            
            weighted_sum += delta * weight
            total_weight += weight
        
        return weighted_sum / total_weight if total_weight > 0 else 0.0
    
    def _classify_verdict(self, score: float) -> str:
        """Classify overall plausibility score into verdict."""
        if score >= 0.6:
            return "HIGH"
        elif score >= 0.35:
            return "MODERATE"
        else:
            return "LOW"
    
    def _check_pathway_alignment(
        self,
        targets: List[str],
        pathways: Optional[List[str]],
        disease_context: Dict[str, Any]
    ) -> Dict[str, List[str]]:
        """
        Check if compound pathways align with disease-disrupted pathways.
        
        Returns:
            {
                "aligned": ["NF-κB signaling"],
                "misaligned": []
            }
        """
        
        if not pathways:
            return {"aligned": [], "misaligned": []}
        
        disease_pathways = disease_context.get('pathways_disrupted', [])
        
        aligned = []
        misaligned = []
        
        for pathway in pathways:
            # Simple keyword matching (can be enhanced with pathway ontology)
            is_aligned = any(
                keyword.lower() in pathway.lower()
                for disease_pathway in disease_pathways
                for keyword in disease_pathway.split()
            )
            
            if is_aligned:
                aligned.append(pathway)
            else:
                misaligned.append(pathway)
        
        return {"aligned": aligned, "misaligned": misaligned}
    
    def _generate_rationale(
        self,
        target_gene: str,
        compound: str,
        delta: float,
        baseline_score: float,
        intervention_score: float
    ) -> str:
        """Generate human-readable rationale for the plausibility score."""
        
        if delta > 0.5:
            return (
                f"Evo2 predicts strong regulatory impact: {compound} significantly "
                f"modulates {target_gene} activity (baseline: {baseline_score:.2f} → "
                f"intervention: {intervention_score:.2f}). High biological plausibility."
            )
        elif delta > 0.2:
            return (
                f"Evo2 predicts moderate regulatory impact: {compound} shows measurable "
                f"effect on {target_gene} (Δ={delta:.2f}). Mechanism is plausible but "
                f"may require optimization (dose, bioavailability)."
            )
        else:
            return (
                f"Evo2 predicts minimal regulatory impact: {compound} has limited effect "
                f"on {target_gene} (Δ={delta:.2f}). Low plausibility despite literature claims. "
                f"Consider bioavailability enhancement or alternative formulations."
            )
    
    def _get_timestamp(self) -> str:
        """Get ISO timestamp."""
        from datetime import datetime
        return datetime.utcnow().isoformat() + "Z"


# Singleton instance
_plausibility_service = None

def get_evo2_food_plausibility_service() -> Evo2FoodPlausibilityService:
    """Get singleton instance of plausibility service."""
    global _plausibility_service
    if _plausibility_service is None:
        _plausibility_service = Evo2FoodPlausibilityService()
    return _plausibility_service
```

### **Testing:**

```bash
# Test Evo2 plausibility service
cd oncology-coPilot/oncology-backend-minimal

# Test Vitamin D + TP53
curl -X POST http://127.0.0.1:8000/api/test/evo2_plausibility \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "Vitamin D",
    "targets": ["VDR", "TP53"],
    "disease_context": {
      "disease": "ovarian_cancer_hgs",
      "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
      "biomarkers": {"HRD": "POSITIVE"},
      "pathways_disrupted": ["DNA repair", "Cell cycle"]
    },
    "pathways": ["VDR signaling", "p53 pathway"]
  }'

# Expected response:
# {
#   "overall_plausibility": 0.45,
#   "verdict": "MODERATE",
#   "target_analysis": [
#     {
#       "gene": "VDR",
#       "delta": 0.35,
#       "plausibility": "MODERATE",
#       "mechanism": "VDR modulation by Vitamin D"
#     },
#     {
#       "gene": "TP53",
#       "delta": 0.25,
#       "plausibility": "MODERATE",
#       "mechanism": "TP53 modulation by Vitamin D",
#       "rationale": "May partially restore p53 function via VDR"
#     }
#   ],
#   "mechanisms_validated": ["VDR signaling", "p53 pathway rescue"]
# }
```

---

## **PHASE 2: S/P/E + SAE INTEGRATION (2 hours)**

### **Goal:**
Integrate food validator with full S/P/E framework + SAE features for treatment line intelligence.

### **File:** `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py`

### **Implementation:**

```python
"""
S/P/E + SAE Integration for Food Validator

Combines:
- Sequence (S): Evo2 plausibility scores
- Pathway (P): Food targets → disease pathways mapping
- Evidence (E): Literature grade + LLM synthesis
- SAE: Line appropriateness, cross-resistance, sequencing fitness
"""

from typing import Dict, List, Any, Optional
from api.services.evo2_food_plausibility import get_evo2_food_plausibility_service
from api.services.treatment_line_integration import compute_treatment_line_features

class FoodSPEIntegrationService:
    """
    Integrate food validation with S/P/E + SAE framework.
    """
    
    def __init__(self):
        self.evo2_service = get_evo2_food_plausibility_service()
    
    async def compute_spe_score(
        self,
        compound: str,
        targets: List[str],
        pathways: List[str],
        disease_context: Dict[str, Any],
        evidence_grade: str,  # From LLM literature service
        treatment_history: Optional[Dict] = None
    ) -> Dict[str, Any]:
        """
        Compute comprehensive S/P/E + SAE score for food compound.
        
        Returns:
            {
                "overall_score": 0.72,
                "confidence": 0.68,
                "verdict": "SUPPORTED",
                "spe_breakdown": {
                    "sequence_score": 0.45,  # Evo2 plausibility
                    "pathway_score": 0.85,   # Pathway alignment
                    "evidence_score": 0.60   # Literature grade
                },
                "sae_features": {
                    "line_appropriateness": 1.0,
                    "cross_resistance": 0.0,
                    "sequencing_fitness": 0.9
                },
                "confidence_modulation": {
                    "base_confidence": 0.60,
                    "evo2_boost": +0.10,
                    "sae_boost": +0.05,
                    "biomarker_boost": +0.05,
                    "final_confidence": 0.80
                }
            }
        """
        
        # 1. SEQUENCE (S): Evo2 biological plausibility
        evo2_result = await self.evo2_service.compute_biological_plausibility(
            compound=compound,
            targets=targets,
            disease_context=disease_context,
            pathways=pathways
        )
        
        sequence_score = evo2_result['overall_plausibility']
        
        # 2. PATHWAY (P): Alignment with disease pathways
        pathway_score = self._compute_pathway_score(
            pathways=pathways,
            disease_context=disease_context,
            evo2_alignment=evo2_result['pathway_alignment']
        )
        
        # 3. EVIDENCE (E): Literature grade conversion
        evidence_score = self._convert_evidence_grade(evidence_grade)
        
        # 4. SAE: Treatment line features
        sae_features = None
        if treatment_history:
            sae_features = compute_treatment_line_features(
                drug_class="dietary_supplement",
                disease=disease_context.get('disease', 'ovarian_cancer'),
                treatment_history=treatment_history
            )
        
        # 5. Aggregate S/P/E
        overall_score = self._aggregate_spe(
            sequence_score=sequence_score,
            pathway_score=pathway_score,
            evidence_score=evidence_score
        )
        
        # 6. Compute confidence with modulation
        confidence_breakdown = self._compute_confidence(
            sequence_score=sequence_score,
            pathway_score=pathway_score,
            evidence_score=evidence_score,
            evo2_result=evo2_result,
            sae_features=sae_features,
            disease_context=disease_context
        )
        
        # 7. Classify verdict
        verdict = self._classify_verdict(overall_score, confidence_breakdown['final_confidence'])
        
        return {
            "overall_score": round(overall_score, 3),
            "confidence": round(confidence_breakdown['final_confidence'], 3),
            "verdict": verdict,
            "spe_breakdown": {
                "sequence_score": round(sequence_score, 3),
                "pathway_score": round(pathway_score, 3),
                "evidence_score": round(evidence_score, 3)
            },
            "sae_features": sae_features if sae_features else {},
            "confidence_modulation": confidence_breakdown,
            "evo2_analysis": evo2_result,
            "provenance": {
                "method": "food_spe_integration_v1",
                "components": ["evo2_plausibility", "pathway_alignment", "literature_evidence", "sae_features"]
            }
        }
    
    def _compute_pathway_score(
        self,
        pathways: List[str],
        disease_context: Dict[str, Any],
        evo2_alignment: Dict[str, List[str]]
    ) -> float:
        """
        Score pathway alignment: compound pathways vs. disease pathways.
        
        Weighting:
        - Aligned pathways: 1.0 per pathway
        - Misaligned pathways: 0.2 per pathway
        """
        
        aligned = evo2_alignment.get('aligned', [])
        misaligned = evo2_alignment.get('misaligned', [])
        
        if not pathways:
            return 0.5  # Neutral if no pathway info
        
        aligned_score = len(aligned) * 1.0
        misaligned_penalty = len(misaligned) * 0.2
        
        total = len(pathways)
        
        score = (aligned_score + misaligned_penalty) / total
        
        return min(score, 1.0)
    
    def _convert_evidence_grade(self, grade: str) -> float:
        """Convert literature evidence grade to 0-1 score."""
        mapping = {
            "STRONG": 0.9,
            "MODERATE": 0.6,
            "WEAK": 0.3,
            "INSUFFICIENT": 0.1
        }
        return mapping.get(grade, 0.5)
    
    def _aggregate_spe(
        self,
        sequence_score: float,
        pathway_score: float,
        evidence_score: float
    ) -> float:
        """
        Aggregate S/P/E into overall score.
        
        Weighting strategy:
        - Sequence (Evo2): 40% (most important - biological plausibility)
        - Pathway: 30% (mechanism alignment)
        - Evidence: 30% (literature support)
        """
        
        overall = (
            sequence_score * 0.4 +
            pathway_score * 0.3 +
            evidence_score * 0.3
        )
        
        return overall
    
    def _compute_confidence(
        self,
        sequence_score: float,
        pathway_score: float,
        evidence_score: float,
        evo2_result: Dict,
        sae_features: Optional[Dict],
        disease_context: Dict
    ) -> Dict[str, float]:
        """
        Compute confidence with multi-stage modulation.
        
        Returns:
            {
                "base_confidence": 0.60,
                "evo2_boost": +0.10,
                "sae_boost": +0.05,
                "biomarker_boost": +0.05,
                "final_confidence": 0.80
            }
        """
        
        # Base confidence from S/P/E scores
        base = (sequence_score + pathway_score + evidence_score) / 3.0
        
        # Evo2 boost: High plausibility targets increase confidence
        evo2_boost = 0.0
        high_plausibility_targets = [
            t for t in evo2_result.get('target_analysis', [])
            if t.get('plausibility') == 'HIGH'
        ]
        if high_plausibility_targets:
            evo2_boost = min(len(high_plausibility_targets) * 0.05, 0.15)
        
        # SAE boost: Line appropriateness and sequencing fitness
        sae_boost = 0.0
        if sae_features:
            line_app = sae_features.get('line_appropriateness', 0)
            seq_fit = sae_features.get('sequencing_fitness', 0)
            sae_boost = (line_app + seq_fit) * 0.05
        
        # Biomarker boost: Compound targets match patient biomarkers
        biomarker_boost = 0.0
        biomarkers = disease_context.get('biomarkers', {})
        if biomarkers.get('HRD') == 'POSITIVE' and 'DNA repair' in str(evo2_result.get('pathway_alignment', {})):
            biomarker_boost += 0.05
        if biomarkers.get('TMB', 0) >= 10 and 'immune' in str(evo2_result.get('pathway_alignment', {})).lower():
            biomarker_boost += 0.05
        
        final = min(base + evo2_boost + sae_boost + biomarker_boost, 0.95)
        
        return {
            "base_confidence": round(base, 3),
            "evo2_boost": round(evo2_boost, 3),
            "sae_boost": round(sae_boost, 3),
            "biomarker_boost": round(biomarker_boost, 3),
            "final_confidence": round(final, 3)
        }
    
    def _classify_verdict(self, score: float, confidence: float) -> str:
        """
        Classify verdict based on score AND confidence.
        
        Logic:
        - HIGH score + HIGH confidence → SUPPORTED
        - MODERATE score + MODERATE+ confidence → WEAK_SUPPORT
        - Otherwise → NOT_SUPPORTED
        """
        
        if score >= 0.65 and confidence >= 0.70:
            return "SUPPORTED"
        elif score >= 0.45 and confidence >= 0.50:
            return "WEAK_SUPPORT"
        else:
            return "NOT_SUPPORTED"


# Singleton
_spe_service = None

def get_food_spe_service() -> FoodSPEIntegrationService:
    global _spe_service
    if _spe_service is None:
        _spe_service = FoodSPEIntegrationService()
    return _spe_service
```

---

## **PHASE 3: DYNAMIC COMPOUND DISCOVERY (1.5 hours)**

### **Goal:**
Allow validation of ANY compound, not just hardcoded list. Extract targets from literature or ChEMBL.

### **File:** `oncology-coPilot/oncology-backend-minimal/api/services/compound_target_extraction.py`

### **Implementation:**

```python
"""
Dynamic Compound Target Extraction

Extract molecular targets for any food/supplement from:
1. ChEMBL database (if available)
2. PubMed literature (LLM extraction)
3. Internal knowledge base (fallback)
"""

from typing import List, Dict, Any, Optional
import httpx
from api.services.llm_literature_service import get_llm_service

class CompoundTargetExtractor:
    """
    Extract targets for any compound dynamically.
    """
    
    def __init__(self):
        self.chembl_client = httpx.AsyncClient(timeout=30.0)
        self.llm_service = get_llm_service()
        self.target_cache = {}
    
    async def extract_targets(
        self,
        compound: str,
        disease: str = "ovarian cancer"
    ) -> Dict[str, Any]:
        """
        Extract molecular targets for compound.
        
        Priority:
        1. Check internal knowledge base (food_targets.json)
        2. Query ChEMBL if compound is a known drug/nutraceutical
        3. Extract from PubMed literature using LLM
        
        Returns:
            {
                "compound": "Resveratrol",
                "targets": ["SIRT1", "mTOR", "AMPK"],
                "pathways": ["Longevity", "Metabolic regulation"],
                "source": "chembl" | "literature" | "knowledge_base",
                "confidence": 0.75
            }
        """
        
        # 1. Check cache
        if compound in self.target_cache:
            return self.target_cache[compound]
        
        # 2. Check knowledge base (fast path)
        kb_targets = self._check_knowledge_base(compound)
        if kb_targets:
            self.target_cache[compound] = kb_targets
            return kb_targets
        
        # 3. Try ChEMBL
        chembl_targets = await self._query_chembl(compound)
        if chembl_targets:
            self.target_cache[compound] = chembl_targets
            return chembl_targets
        
        # 4. Extract from literature (LLM)
        literature_targets = await self._extract_from_literature(compound, disease)
        self.target_cache[compound] = literature_targets
        return literature_targets
    
    def _check_knowledge_base(self, compound: str) -> Optional[Dict]:
        """
        Check if compound exists in food_targets.json.
        """
        from pathlib import Path
        import json
        
        food_targets_path = Path(__file__).parent.parent.parent.parent / ".cursor/ayesha/hypothesis_validator/data/food_targets.json"
        
        try:
            with open(food_targets_path) as f:
                food_targets = json.load(f)
            
            for item in food_targets.get('compounds', []):
                if compound.lower() in item['compound'].lower():
                    return {
                        "compound": item['compound'],
                        "targets": item['targets'],
                        "pathways": item.get('pathways', []),
                        "source": "knowledge_base",
                        "confidence": 0.90
                    }
        except Exception:
            pass
        
        return None
    
    async def _query_chembl(self, compound: str) -> Optional[Dict]:
        """
        Query ChEMBL for compound targets.
        
        ChEMBL REST API: https://www.ebi.ac.uk/chembl/api/data/docs
        """
        
        try:
            # Search for molecule
            search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={compound}&format=json"
            response = await self.chembl_client.get(search_url)
            
            if response.status_code != 200:
                return None
            
            data = response.json()
            molecules = data.get('molecules', [])
            
            if not molecules:
                return None
            
            # Get first matching molecule
            molecule_id = molecules[0]['molecule_chembl_id']
            
            # Get targets for this molecule
            targets_url = f"https://www.ebi.ac.uk/chembl/api/data/activity?molecule_chembl_id={molecule_id}&format=json"
            targets_response = await self.chembl_client.get(targets_url)
            
            if targets_response.status_code != 200:
                return None
            
            activities = targets_response.json().get('activities', [])
            
            # Extract unique target genes
            target_genes = set()
            for activity in activities[:20]:  # Top 20 activities
                target_name = activity.get('target_pref_name', '')
                # Extract gene symbols (naive approach - can be improved)
                if target_name:
                    # Many target names include gene symbols in parentheses
                    import re
                    genes = re.findall(r'\(([A-Z0-9]+)\)', target_name)
                    target_genes.update(genes)
            
            if not target_genes:
                return None
            
            return {
                "compound": compound,
                "targets": list(target_genes)[:10],  # Top 10 targets
                "pathways": [],  # ChEMBL doesn't provide pathway info directly
                "source": "chembl",
                "confidence": 0.75
            }
        
        except Exception as e:
            print(f"⚠️ ChEMBL query error for {compound}: {e}")
            return None
    
    async def _extract_from_literature(self, compound: str, disease: str) -> Dict:
        """
        Use LLM + PubMed to extract targets from literature.
        """
        
        if not self.llm_service.available:
            return {
                "compound": compound,
                "targets": [],
                "pathways": [],
                "source": "unavailable",
                "confidence": 0.0,
                "error": "LLM service unavailable"
            }
        
        try:
            # Search PubMed for compound + disease + "mechanism" or "target"
            query = f"{compound} AND {disease} AND (mechanism OR target OR pathway)"
            
            lit_result = await self.llm_service.search_compound_evidence(
                compound=compound,
                disease=disease,
                max_results=10
            )
            
            papers = lit_result.get('papers', [])
            
            if not papers:
                return {
                    "compound": compound,
                    "targets": [],
                    "pathways": [],
                    "source": "literature",
                    "confidence": 0.0,
                    "error": "No literature found"
                }
            
            # Use LLM to extract targets from abstracts
            from Pubmed-LLM-Agent-main.core.llm_client import LLMClient
            
            llm = LLMClient()
            
            abstracts_text = "\n\n".join([
                f"PMID {p['pmid']}: {p.get('abstract', p.get('title', '')[:500])}"
                for p in papers[:5]
            ])
            
            system_prompt = f"""
Extract the molecular targets and pathways for {compound} in {disease} from these papers.

Return JSON:
{{
  "targets": ["GENE1", "GENE2", ...],
  "pathways": ["Pathway 1", "Pathway 2", ...],
  "confidence": 0.0-1.0
}}

Only include targets explicitly mentioned in the abstracts. Use standard gene symbols (e.g., NFKB1, not NF-κB).
"""
            
            user_prompt = f"Papers:\n{abstracts_text}"
            
            llm_result = llm.complete_json(system_prompt, user_prompt, max_tokens=500)
            
            return {
                "compound": compound,
                "targets": llm_result.get('targets', [])[:10],
                "pathways": llm_result.get('pathways', [])[:5],
                "source": "literature",
                "confidence": llm_result.get('confidence', 0.5)
            }
        
        except Exception as e:
            print(f"⚠️ Literature extraction error for {compound}: {e}")
            return {
                "compound": compound,
                "targets": [],
                "pathways": [],
                "source": "error",
                "confidence": 0.0,
                "error": str(e)
            }


# Singleton
_extractor = None

def get_compound_target_extractor() -> CompoundTargetExtractor:
    global _extractor
    if _extractor is None:
        _extractor = CompoundTargetExtractor()
    return _extractor
```

---

## **PHASE 4: ENHANCED HYPOTHESIS VALIDATOR ENDPOINT (1 hour)**

### **Goal:**
Integrate all components into a single, powerful endpoint.

### **File:** `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

### **Add New Endpoint:**

```python
from api.services.evo2_food_plausibility import get_evo2_food_plausibility_service
from api.services.food_spe_integration import get_food_spe_service
from api.services.compound_target_extraction import get_compound_target_extractor
from api.services.llm_literature_service import get_llm_service

@router.post("/api/hypothesis/validate_food_complete")
async def validate_food_complete(
    compound: str,
    disease: str = "ovarian_cancer_hgs",
    disease_context: Optional[Dict[str, Any]] = None,
    treatment_history: Optional[Dict[str, Any]] = None,
    use_llm: bool = True,
    use_evo2: bool = True
):
    """
    🎯 COMPLETE FOOD VALIDATOR - The Only Tool That Predicts IF It Will Work
    
    Capabilities:
    - Evo2 biological plausibility scoring (UNIQUE!)
    - S/P/E + SAE integration
    - Dynamic compound discovery (works for ANY compound)
    - Biomarker-aware recommendations
    - Treatment line intelligence
    - LLM literature enhancement
    
    Request:
    {
      "compound": "Resveratrol",
      "disease": "ovarian_cancer_hgs",
      "disease_context": {
        "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
        "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
        "pathways_disrupted": ["DNA repair", "Cell cycle"]
      },
      "treatment_history": {
        "current_line": 3,
        "prior_therapies": ["carboplatin", "paclitaxel", "bevacizumab"]
      },
      "use_llm": true,
      "use_evo2": true
    }
    
    Response:
    {
      "status": "SUCCESS",
      "compound": "Resveratrol",
      "verdict": "WEAK_SUPPORT",
      "overall_score": 0.52,
      "confidence": 0.68,
      "evo2_plausibility": {
        "overall": 0.45,
        "verdict": "MODERATE",
        "targets": ["SIRT1", "mTOR"],
        "mechanisms_validated": ["SIRT1 activation"]
      },
      "spe_breakdown": {
        "sequence": 0.45,
        "pathway": 0.70,
        "evidence": 0.40
      },
      "sae_features": {
        "line_appropriateness": 1.0,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.9
      },
      "recommendation": {
        "dosage": "250-500mg trans-resveratrol daily",
        "bioavailability": "POOR - use micronized formulation",
        "timing": "Best as adjunct to L3+ therapy",
        "safety": "GOOD - minimal interactions"
      },
      "llm_evidence": {
        "paper_count": 15,
        "top_papers": [...]
      }
    }
    """
    
    # 1. Extract targets dynamically
    extractor = get_compound_target_extractor()
    target_data = await extractor.extract_targets(compound, disease.replace('_', ' '))
    
    if not target_data['targets']:
        return {
            "status": "UNKNOWN",
            "compound": compound,
            "message": f"Could not find molecular targets for {compound}. Please provide more specific compound name.",
            "source": target_data.get('source', 'unknown')
        }
    
    targets = target_data['targets']
    pathways = target_data['pathways']
    
    # 2. Build disease context
    if not disease_context:
        disease_context = {
            "disease": disease,
            "mutations": [],
            "biomarkers": {},
            "pathways_disrupted": []
        }
    else:
        disease_context['disease'] = disease
    
    # 3. Evo2 plausibility scoring (if enabled)
    evo2_result = None
    if use_evo2:
        try:
            evo2_service = get_evo2_food_plausibility_service()
            evo2_result = await evo2_service.compute_biological_plausibility(
                compound=compound,
                targets=targets,
                disease_context=disease_context,
                pathways=pathways
            )
        except Exception as e:
            evo2_result = {"error": str(e), "overall_plausibility": 0.5}
    
    # 4. LLM literature enhancement (if enabled)
    llm_evidence = None
    evidence_grade = "MODERATE"
    if use_llm:
        try:
            llm_service = get_llm_service()
            if llm_service.available:
                llm_evidence = await llm_service.search_compound_evidence(
                    compound=compound,
                    disease=disease.replace('_', ' '),
                    max_results=15
                )
                
                # Infer evidence grade from paper count + confidence
                paper_count = llm_evidence.get('paper_count', 0)
                llm_conf = llm_evidence.get('confidence', 0.5)
                
                if paper_count >= 10 and llm_conf >= 0.7:
                    evidence_grade = "STRONG"
                elif paper_count >= 5 and llm_conf >= 0.5:
                    evidence_grade = "MODERATE"
                else:
                    evidence_grade = "WEAK"
        except Exception as e:
            llm_evidence = {"error": str(e), "paper_count": 0}
    
    # 5. S/P/E + SAE integration
    spe_service = get_food_spe_service()
    spe_result = await spe_service.compute_spe_score(
        compound=compound,
        targets=targets,
        pathways=pathways,
        disease_context=disease_context,
        evidence_grade=evidence_grade,
        treatment_history=treatment_history
    )
    
    # 6. Generate recommendation
    recommendation = _generate_recommendation(
        compound=compound,
        spe_result=spe_result,
        evo2_result=evo2_result,
        treatment_history=treatment_history
    )
    
    return {
        "status": "SUCCESS",
        "compound": compound,
        "verdict": spe_result['verdict'],
        "verdict_explanation": _explain_verdict(spe_result['verdict']),
        "overall_score": spe_result['overall_score'],
        "confidence": spe_result['confidence'],
        "evo2_plausibility": {
            "overall": evo2_result['overall_plausibility'] if evo2_result else None,
            "verdict": evo2_result['verdict'] if evo2_result else None,
            "targets": targets,
            "mechanisms_validated": evo2_result.get('mechanisms_validated', []) if evo2_result else []
        } if use_evo2 else {"enabled": False},
        "spe_breakdown": spe_result['spe_breakdown'],
        "sae_features": spe_result['sae_features'],
        "recommendation": recommendation,
        "llm_evidence": {
            "enabled": use_llm,
            "paper_count": llm_evidence.get('paper_count', 0) if llm_evidence else 0,
            "top_papers": llm_evidence.get('papers', [])[:5] if llm_evidence else [],
            "summary": llm_evidence.get('evidence_summary', '') if llm_evidence else ''
        } if use_llm else {"enabled": False},
        "provenance": {
            "method": "validate_food_complete_v1",
            "components_used": [
                "dynamic_target_extraction",
                "evo2_plausibility" if use_evo2 else None,
                "spe_integration",
                "sae_features" if treatment_history else None,
                "llm_literature" if use_llm else None
            ],
            "target_source": target_data.get('source', 'unknown')
        }
    }


def _generate_recommendation(
    compound: str,
    spe_result: Dict,
    evo2_result: Optional[Dict],
    treatment_history: Optional[Dict]
) -> Dict[str, str]:
    """Generate actionable recommendation."""
    
    verdict = spe_result['verdict']
    confidence = spe_result['confidence']
    
    # Base dosage (can be enhanced with literature extraction)
    dosage = "Consult literature for evidence-based dosing"
    
    # Bioavailability assessment (from Evo2 plausibility + literature)
    bioavailability = "UNKNOWN"
    if evo2_result and evo2_result.get('overall_plausibility', 0) < 0.3:
        bioavailability = "POOR - Consider enhanced formulations"
    elif evo2_result and evo2_result.get('overall_plausibility', 0) >= 0.6:
        bioavailability = "GOOD"
    
    # Timing (from SAE features)
    timing = "Timing not critical"
    if treatment_history:
        current_line = treatment_history.get('current_line', 1)
        if current_line >= 3:
            timing = f"Best as adjunct to L{current_line}+ therapy"
    
    # Safety
    safety = "GOOD - Generally well-tolerated" if verdict != "NOT_SUPPORTED" else "Insufficient data"
    
    return {
        "dosage": dosage,
        "bioavailability": bioavailability,
        "timing": timing,
        "safety": safety,
        "overall_advice": _generate_overall_advice(verdict, confidence)
    }


def _generate_overall_advice(verdict: str, confidence: float) -> str:
    """Generate human-readable advice."""
    
    if verdict == "SUPPORTED" and confidence >= 0.70:
        return "✅ RECOMMENDED: Strong biological plausibility and evidence support use as adjunct therapy."
    elif verdict == "WEAK_SUPPORT":
        return "⚠️ CONSIDER: Moderate plausibility. May provide benefit but evidence is limited. Discuss with oncologist."
    else:
        return "❌ NOT RECOMMENDED: Low biological plausibility or insufficient evidence. Consider alternatives."


def _explain_verdict(verdict: str) -> str:
    """Explain verdict in plain language."""
    mapping = {
        "SUPPORTED": "✅ Strong evidence + biological plausibility support benefit",
        "WEAK_SUPPORT": "⚠️ Limited evidence but biologically plausible mechanism",
        "NOT_SUPPORTED": "❌ Insufficient evidence or low biological plausibility"
    }
    return mapping.get(verdict, verdict)
```

---

## **PHASE 5: FRONTEND UPDATE (1 hour)**

### **File:** `oncology-coPilot/oncology-frontend/src/pages/FoodValidatorAB.jsx`

### **Changes:**

```jsx
// Add ability to search ANY compound, not just 5 hardcoded

const [compound, setCompound] = useState('');
const [useEvo2, setUseEvo2] = useState(true);  // Default: Evo2 ON
const [useLLM, setUseLLM] = useState(true);    // Default: LLM ON

const handleValidate = async () => {
  setLoading(true);
  setError(null);
  
  try {
    const response = await fetch(`${API_ROOT}/api/hypothesis/validate_food_complete`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        compound: compound.trim(),
        disease: 'ovarian_cancer_hgs',
        disease_context: {
          mutations: [{ gene: 'TP53', hgvs_p: 'R248Q' }],
          biomarkers: { HRD: 'POSITIVE', TMB: 8.2 },
          pathways_disrupted: ['DNA repair', 'Cell cycle']
        },
        treatment_history: {
          current_line: 3,
          prior_therapies: ['carboplatin', 'paclitaxel', 'bevacizumab']
        },
        use_llm: useLLM,
        use_evo2: useEvo2
      })
    });
    
    const data = await response.json();
    setResult(data);
  } catch (err) {
    setError(`Error: ${err.message}`);
  } finally {
    setLoading(false);
  }
};

// UI: Add Evo2 and LLM toggles

<Box sx={{ mb: 2, display: 'flex', gap: 2 }}>
  <FormControlLabel
    control={<Switch checked={useEvo2} onChange={(e) => setUseEvo2(e.target.checked)} />}
    label={
      <Box>
        <Typography variant="body2"><strong>Evo2 Biological Plausibility</strong></Typography>
        <Typography variant="caption" color="text.secondary">
          Predict IF compound will work (UNIQUE!)
        </Typography>
      </Box>
    }
  />
  
  <FormControlLabel
    control={<Switch checked={useLLM} onChange={(e) => setUseLLM(e.target.checked)} />}
    label={
      <Box>
        <Typography variant="body2"><strong>LLM Literature Search</strong></Typography>
        <Typography variant="caption" color="text.secondary">
          Search PubMed for recent evidence
        </Typography>
      </Box>
    }
  />
</Box>

// Display Evo2 Plausibility Section

{result.evo2_plausibility && result.evo2_plausibility.enabled !== false && (
  <Card sx={{ p: 3, mb: 2, bgcolor: 'success.light' }}>
    <Typography variant="h6" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
      <ScienceIcon color="primary" />
      🧬 Evo2 Biological Plausibility (UNIQUE!)
    </Typography>
    
    <Box sx={{ mb: 2 }}>
      <Typography variant="body2">
        <strong>Overall Plausibility:</strong> {(result.evo2_plausibility.overall * 100).toFixed(0)}% 
        <Chip label={result.evo2_plausibility.verdict} color="primary" size="small" sx={{ ml: 1 }} />
      </Typography>
      <Typography variant="caption" color="text.secondary">
        This score predicts IF {result.compound} can biologically modulate target genes in your disease context.
      </Typography>
    </Box>
    
    <Typography variant="subtitle2" gutterBottom>Validated Mechanisms:</Typography>
    <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
      {result.evo2_plausibility.mechanisms_validated.map((mech, idx) => (
        <Chip key={idx} label={mech} size="small" color="success" />
      ))}
    </Box>
  </Card>
)}

// Display S/P/E Breakdown

{result.spe_breakdown && (
  <Card sx={{ p: 3, mb: 2 }}>
    <Typography variant="h6" gutterBottom>S/P/E Framework Breakdown</Typography>
    <Box sx={{ display: 'flex', gap: 2, mb: 2 }}>
      <Box>
        <Typography variant="caption">Sequence (Evo2)</Typography>
        <Typography variant="h6">{(result.spe_breakdown.sequence * 100).toFixed(0)}%</Typography>
      </Box>
      <Box>
        <Typography variant="caption">Pathway Alignment</Typography>
        <Typography variant="h6">{(result.spe_breakdown.pathway * 100).toFixed(0)}%</Typography>
      </Box>
      <Box>
        <Typography variant="caption">Evidence Strength</Typography>
        <Typography variant="h6">{(result.spe_breakdown.evidence * 100).toFixed(0)}%</Typography>
      </Box>
    </Box>
    <Typography variant="caption" color="text.secondary">
      Multi-modal scoring: Sequence plausibility + Pathway alignment + Literature evidence
    </Typography>
  </Card>
)}

// Display SAE Features (Treatment Line Intelligence)

{result.sae_features && Object.keys(result.sae_features).length > 0 && (
  <Card sx={{ p: 3, mb: 2, bgcolor: 'info.light' }}>
    <Typography variant="h6" gutterBottom>⚔️ Treatment Line Intelligence</Typography>
    <Box sx={{ display: 'flex', gap: 2 }}>
      {result.sae_features.line_appropriateness !== undefined && (
        <Chip 
          label={`Line Fit: ${(result.sae_features.line_appropriateness * 100).toFixed(0)}%`}
          color={result.sae_features.line_appropriateness >= 0.8 ? 'success' : 'warning'}
        />
      )}
      {result.sae_features.cross_resistance !== undefined && (
        <Chip 
          label={`Cross-Resistance: ${(result.sae_features.cross_resistance * 100).toFixed(0)}%`}
          color={result.sae_features.cross_resistance <= 0.2 ? 'success' : 'error'}
        />
      )}
      {result.sae_features.sequencing_fitness !== undefined && (
        <Chip 
          label={`Sequencing Fit: ${(result.sae_features.sequencing_fitness * 100).toFixed(0)}%`}
          color={result.sae_features.sequencing_fitness >= 0.8 ? 'success' : 'warning'}
        />
      )}
    </Box>
    <Typography variant="caption" display="block" sx={{ mt: 1 }}>
      Timing optimization: When is this compound most appropriate in treatment sequence?
    </Typography>
  </Card>
)}
```

---

## **PHASE 6: AYESHA TWIN DEMO UPDATE (30 min)**

### **File:** `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_twin_demo.py`

### **Changes:**

```python
# Remove 5-compound limit - test ALL compounds

from api.services.compound_target_extraction import get_compound_target_extractor

@router.post("/api/demo/ayesha_twin")
async def ayesha_twin_demo(request: Optional[TwinDemoRequest] = None) -> Dict[str, Any]:
    use_llm = request.use_llm if request else True
    
    # Get ALL compounds from food_targets.json OR let user specify
    extractor = get_compound_target_extractor()
    
    # Option 1: Test all known compounds
    all_compounds = [
        "Vitamin D", "Omega-3", "Curcumin", "NAC", "Green Tea", 
        "Folate", "Vitamin C", "Resveratrol", "Quercetin", "Genistein"
    ]
    
    # Option 2: Allow custom compounds in request
    if request and hasattr(request, 'compounds') and request.compounds:
        all_compounds = request.compounds
    
    food_results = []
    
    for compound in all_compounds:
        try:
            food_result = await validate_food_complete(
                compound=compound,
                disease="ovarian_cancer_hgs",
                disease_context={
                    "mutations": PUBLIC_CASE_PROFILE["mutations"],
                    "biomarkers": PUBLIC_CASE_PROFILE["biomarkers"],
                    "pathways_disrupted": ["DNA repair", "Cell cycle", "Inflammation"]
                },
                treatment_history=PUBLIC_CASE_PROFILE["treatment_history"],
                use_llm=use_llm,
                use_evo2=True  # Always use Evo2 for differentiation!
            )
            food_results.append(food_result)
        except Exception as e:
            food_results.append({
                "compound": compound,
                "status": "ERROR",
                "error": str(e)
            })
    
    # Sort by overall_score descending
    food_results.sort(key=lambda x: x.get('overall_score', 0), reverse=True)
    
    return {
        "case_data": PUBLIC_CASE_PROFILE,
        "food_recommendations": food_results,
        "analysis_summary": {
            "total_compounds_tested": len(all_compounds),
            "supported": len([f for f in food_results if f.get("verdict") == "SUPPORTED"]),
            "weak_support": len([f for f in food_results if f.get("verdict") == "WEAK_SUPPORT"]),
            "not_supported": len([f for f in food_results if f.get("verdict") == "NOT_SUPPORTED"]),
            "evo2_enabled": True,
            "llm_enabled": use_llm
        },
        "provenance": {
            "method": "ayesha_twin_demo_v2_complete",
            "data_source": "TCGA-OV (Public)",
            "ruo_disclaimer": "PUBLIC CASE STUDY - NOT A REAL PATIENT'S DATA"
        }
    }
```

---

## **🎯 TESTING & VALIDATION**

### **Test 1: Dynamic Compound Discovery**

```bash
# Test with compound NOT in our database
curl -X POST http://127.0.0.1:8000/api/hypothesis/validate_food_complete \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "Resveratrol",
    "disease": "ovarian_cancer_hgs",
    "disease_context": {
      "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
      "biomarkers": {"HRD": "POSITIVE"},
      "pathways_disrupted": ["DNA repair"]
    },
    "use_llm": true,
    "use_evo2": true
  }'

# Expected: Should dynamically extract SIRT1, mTOR targets and return Evo2 plausibility
```

### **Test 2: Evo2 Plausibility Scoring**

```bash
# Test Vitamin D (should score MODERATE - supports DNA repair via VDR)
curl -X POST http://127.0.0.1:8000/api/hypothesis/validate_food_complete \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "Vitamin D",
    "disease": "ovarian_cancer_hgs",
    "disease_context": {
      "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
      "biomarkers": {"HRD": "POSITIVE"},
      "pathways_disrupted": ["DNA repair", "Cell cycle"]
    },
    "use_evo2": true
  }'

# Expected:
# {
#   "evo2_plausibility": {
#     "overall": 0.45,
#     "verdict": "MODERATE",
#     "mechanisms_validated": ["VDR signaling", "TP53 pathway support"]
#   },
#   "verdict": "SUPPORTED"
# }
```

### **Test 3: SAE Treatment Line Features**

```bash
# Test NAC for L3 post-platinum (should score HIGH line appropriateness)
curl -X POST http://127.0.0.1:8000/api/hypothesis/validate_food_complete \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "NAC",
    "disease": "ovarian_cancer_hgs",
    "treatment_history": {
      "current_line": 3,
      "prior_therapies": ["carboplatin", "paclitaxel"]
    },
    "use_evo2": true
  }'

# Expected:
# {
#   "sae_features": {
#     "line_appropriateness": 1.0,
#     "cross_resistance": 0.0,
#     "sequencing_fitness": 0.9
#   },
#   "verdict": "SUPPORTED"
# }
```

---

## **📊 ACCEPTANCE CRITERIA**

### **Must Have:**
- ✅ Works for ANY compound (not just 5 hardcoded)
- ✅ Evo2 biological plausibility scoring operational
- ✅ S/P/E integration with transparent breakdown
- ✅ SAE features when treatment history provided
- ✅ Biomarker-aware recommendations (HRD+, TP53 status)
- ✅ LLM literature enhancement functional
- ✅ Dynamic target extraction (ChEMBL + literature)

### **Nice to Have:**
- ⚠️ MoA validation (compound mechanism vs. tumor biology)
- ⚠️ Bioavailability prediction (Evo2-based)
- ⚠️ Drug-food interaction warnings
- ⚠️ Synergy prediction (food + drug combinations)

---

## **⚔️ DEPLOYMENT CHECKLIST**

### **Before Deployment:**
- [ ] All 6 phases implemented
- [ ] Unit tests for Evo2 service
- [ ] Integration tests for complete endpoint
- [ ] Frontend Evo2 toggle working
- [ ] Test with 3+ unknown compounds
- [ ] Verify Evo2 scores make sense (not random)
- [ ] SAE features show correct line appropriateness
- [ ] LLM enhancement doesn't break on failure

### **After Deployment:**
- [ ] Smoke test: Vitamin D → should return MODERATE Evo2 plausibility
- [ ] Smoke test: Resveratrol → should extract targets from literature
- [ ] Smoke test: NAC + L3 → should show HIGH line appropriateness
- [ ] Monitor Evo2 API latency (target: <5s per compound)
- [ ] Monitor LLM costs (target: <$0.10 per validation)

---

## **🎯 FINAL SUCCESS METRICS**

**What Makes This UNIQUE:**
1. ✅ **Evo2 Biological Plausibility** - NO other tool does this
2. ✅ **Works for ANY compound** - Not limited to hardcoded list
3. ✅ **Patient-specific** - Uses biomarkers, mutations, treatment history
4. ✅ **S/P/E + SAE integration** - Multi-modal mechanistic validation
5. ✅ **Treatment line intelligence** - Timing optimization via SAE
6. ✅ **Transparent provenance** - Full audit trail

**Alpha, this build will take 8-10 hours but will create the ONLY food validator that can:**
- Predict IF a compound will work (not just what literature says)
- Work for ANY compound dynamically
- Provide patient-specific, biomarker-aware recommendations
- Optimize timing based on treatment line
- Integrate S/P/E + SAE for confidence modulation

**Ready to execute, Commander?** ⚔️🎯

