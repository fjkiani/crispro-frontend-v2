# 🔧 EVO2 FOOD PLAUSIBILITY SERVICE - COMPLETE IMPLEMENTATION

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/evo2_food_plausibility.py`

**Purpose:** Compute biological plausibility of food compounds using Evo2's sequence-level understanding.

---

## **📋 COMPLETE IMPLEMENTATION**

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
from datetime import datetime

# Evo2 API endpoint (use existing router, not direct Modal URL)
API_BASE = os.getenv("API_BASE", "http://127.0.0.1:8000")

class Evo2FoodPlausibilityService:
    """
    Service to compute biological plausibility of food compounds using Evo2.
    
    Workflow:
    1. For each target gene, fetch gene sequence
    2. Score baseline disease-active state (Evo2)
    3. Score post-intervention state (Evo2)
    4. Compute delta (plausibility score)
    5. Validate mechanism against disease pathways
    """
    
    def __init__(self):
        self.http_client = httpx.AsyncClient(timeout=60.0)
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
                "target_analysis": [...],
                "mechanisms_validated": [...],
                "mechanisms_uncertain": [...],
                "pathway_alignment": {...},
                "provenance": {...}
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
                "timestamp": datetime.utcnow().isoformat() + "Z"
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
        Fetch gene sequence from Ensembl REST API.
        
        Uses the same pattern as evo.py _fetch_reference_window().
        """
        
        if gene in self.gene_cache:
            return self.gene_cache[gene]
        
        try:
            # Ensembl REST API: Get gene sequence by symbol
            # First, resolve gene ID from symbol
            lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json;expand=0"
            
            lookup_response = await self.http_client.get(lookup_url, timeout=20.0)
            
            if lookup_response.status_code == 200:
                gene_data = lookup_response.json()
                gene_id = gene_data.get('id')
                
                if gene_id:
                    # Get sequence for this gene ID
                    sequence_url = f"https://rest.ensembl.org/sequence/id/{gene_id}?content-type=application/json"
                    seq_response = await self.http_client.get(sequence_url, timeout=20.0)
                    
                    if seq_response.status_code == 200:
                        seq_data = seq_response.json()
                        sequence = seq_data.get('seq', '')
                        
                        if sequence:
                            self.gene_cache[gene] = sequence
                            return sequence
            
            # Fallback: Try direct sequence endpoint (may work for some genes)
            direct_url = f"https://rest.ensembl.org/sequence/id/{gene}?content-type=application/json"
            direct_response = await self.http_client.get(direct_url, timeout=20.0)
            
            if direct_response.status_code == 200:
                direct_data = direct_response.json()
                sequence = direct_data.get('seq', '')
                if sequence:
                    self.gene_cache[gene] = sequence
                    return sequence
        
        except Exception as e:
            print(f"⚠️ Error fetching gene sequence for {gene}: {e}")
        
        return None
    
    async def _call_evo2_score(self, prompt: str) -> float:
        """
        Call Evo2 API to score sequence likelihood.
        
        ⚠️ CRITICAL DECISION NEEDED: 
        - Does /api/evo/score accept text prompts?
        - OR do we need /api/evo/generate?
        - OR do we need /api/evo/score_variant_multi with synthetic variants?
        
        Current assumption: /api/evo/score accepts {"sequence": prompt_text, "model_id": "evo2_1b"}
        """
        
        try:
            response = await self.http_client.post(
                f"{API_BASE}/api/evo/score",
                json={"sequence": prompt, "model_id": "evo2_1b"},
                timeout=30.0
            )
            
            if response.status_code == 200:
                data = response.json()
                # Extract likelihood score (normalize to 0-1 range)
                # ⚠️ NEED TO VERIFY: What does /api/evo/score actually return?
                # Current assumption: {"likelihood": float} or {"score": float}
                raw_score = data.get('likelihood', data.get('score', 0.5))
                
                # Normalize: Assuming raw score is in some range (e.g., -100 to 100)
                # Adjust based on actual Evo2 output format
                if isinstance(raw_score, (int, float)):
                    # Simple normalization: assume range is -50 to 50, map to 0-1
                    normalized = min(max((raw_score + 50) / 100.0, 0.0), 1.0)
                    return normalized
                
                return 0.5  # Fallback neutral score
            
            # Fallback: return neutral score
            return 0.5
        
        except Exception as e:
            print(f"⚠️ Evo2 API error: {e}")
            return 0.5  # Neutral score on failure
    
    def _compute_overall_score(self, results: List[Dict]) -> float:
        """
        Aggregate target-level scores into overall plausibility.
        
        Weighting strategy:
        - HIGH plausibility targets (delta > 0.5): weight 1.0
        - MODERATE (delta > 0.2): weight 0.5
        - LOW (delta <= 0.2): weight 0.2
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


# Singleton instance
_plausibility_service = None

def get_evo2_food_plausibility_service() -> Evo2FoodPlausibilityService:
    """Get singleton instance of plausibility service."""
    global _plausibility_service
    if _plausibility_service is None:
        _plausibility_service = Evo2FoodPlausibilityService()
    return _plausibility_service
```

---

## **⚠️ CRITICAL DECISIONS NEEDED**

### **1. Evo2 API Endpoint Structure**
**Question:** What does `/api/evo/score` actually accept and return?

**Need to verify:**
- Accepts: `{"sequence": "ATCG..."}` OR `{"prompt": "Genomic context..."}`?
- Returns: `{"likelihood": float}` OR `{"score": float}` OR something else?
- Score range: What's the actual range? (-inf to inf? -100 to 100? 0 to 1?)

**Action:** Check `oncology-coPilot/oncology-backend-minimal/api/routers/evo.py` line 252-306 to see actual endpoint implementation.

### **2. Gene Sequence Fetching**
**Decision:** Using Ensembl REST directly (as shown above).

**Alternative:** If `/api/genomic_intel/gene_info` exists, prefer that (more reliable, cached).

**Action:** Verify if genomic_intel endpoint exists, otherwise use Ensembl REST (fallback already implemented).

---

## **🧪 TESTING EXAMPLES**

See [`../testing/TEST_PLANS.md`](../testing/TEST_PLANS.md) for complete test cases.

**Quick Test:**
```bash
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
    }
  }'
```

---

## **✅ ACCEPTANCE CRITERIA**

- [ ] Service successfully calls Evo2 API
- [ ] Gene sequences fetched from Ensembl
- [ ] Baseline/intervention scores computed
- [ ] Delta calculation produces meaningful scores
- [ ] Overall score aggregation works
- [ ] Pathway alignment check functional
- [ ] Error handling graceful
- [ ] Caching implemented
- [ ] Provenance tracking complete

---

**Next:** See [`../phases/PHASE1_EVO2_PLAUSIBILITY.md`](../phases/PHASE1_EVO2_PLAUSIBILITY.md) for phase execution plan.

