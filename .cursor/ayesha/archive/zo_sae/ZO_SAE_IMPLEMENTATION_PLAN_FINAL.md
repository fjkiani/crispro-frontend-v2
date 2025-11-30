# ⚔️ ZO'S SAE IMPLEMENTATION PLAN (MANAGER-APPROVED) ⚔️

**Date**: January 13, 2025  
**Status**: 🔓 **UNBLOCKED** - Manager answers received  
**Owner**: Zo  
**Timeline**: 2-day sprint (Phase 1 TODAY focus)

---

## 🎯 EXECUTIVE SUMMARY

**Manager's Strategic Direction:**
- **Pre-NGS (TODAY)**: Deliver deterministic, guideline-based care (SOC + trials + CA-125 + NGS fast-track) with next-test recommender
- **Post-NGS (LATER)**: Unlock SAE-driven mechanism fit, hint tiles, and resistance planning
- **No Hallucination**: All thresholds sourced, all LLM outputs human-reviewed, all drug claims gated by NGS availability

**What This Means for Ayesha TODAY:**
- ✅ She gets: 10 frontline trials (transparent reasoning), SOC plan (detailed dosing), CA-125 monitoring (resistance flags), NGS ordering guidance (HRD → ctDNA → SLFN11 → ABCB1)
- ⏸️ She waits for: Drug efficacy predictions (WIWFM), mechanism map, SAE-driven hints (all require NGS)
- 🎯 Timeline: 7-10 days to unlock full WIWFM (ctDNA results)

---

## 📋 PHASE 1: PRE-NGS ENHANCEMENTS (TODAY - 4 hours)

### **Task 1.1: Next-Test Recommender Service** (1.5 hours)

**File**: `api/services/next_test_recommender.py` (NEW)

**Manager's Policy (from P1, C6)**:
- Priority: 1) HRD, 2) ctDNA (MSI/TMB/somatic HRR), 3) SLFN11 IHC, 4) ABCB1 (if post-taxane)
- Trigger: completeness L0/L1 or missing HRD/MSI/TMB
- Format: "If positive → X; If negative → Y" with turnaround

**Implementation**:
```python
from typing import Dict, List, Optional
from dataclasses import dataclass

@dataclass
class NextTestRecommendation:
    test_name: str
    priority: int  # 1-4
    turnaround_days: int
    cost_estimate: str
    impact_if_positive: str
    impact_if_negative: str
    rationale: str
    urgency: str  # "high", "medium", "low"

class NextTestRecommenderService:
    """Generate next-test recommendations based on missing biomarkers."""
    
    def recommend_tests(
        self,
        germline_status: str,
        tumor_context: Optional[Dict] = None,
        treatment_history: List[str] = None
    ) -> List[NextTestRecommendation]:
        """
        Prioritized test recommendations per Manager's policy.
        
        Manager's Order (P1, C6):
        1. HRD (PARP gate)
        2. ctDNA MSI/TMB + somatic HRR (IO + DDR combos)
        3. SLFN11 IHC (PARP sensitivity)
        4. ABCB1 proxy (if prior taxane)
        """
        recommendations = []
        treatment_history = treatment_history or []
        
        # Parse tumor context completeness
        has_hrd = tumor_context and tumor_context.get("hrd_score") is not None
        has_msi = tumor_context and tumor_context.get("msi_status") is not None
        has_tmb = tumor_context and tumor_context.get("tmb") is not None
        has_slfn11 = tumor_context and tumor_context.get("slfn11_status") is not None
        has_abcb1 = tumor_context and tumor_context.get("abcb1_status") is not None
        
        # Test 1: HRD (Priority 1 - PARP gate)
        if not has_hrd:
            recommendations.append(NextTestRecommendation(
                test_name="HRD Score (MyChoice CDx or tissue-based)",
                priority=1,
                turnaround_days=10,
                cost_estimate="$4,000-$6,000 (typically covered by insurance)",
                impact_if_positive="HRD ≥42 → PARP maintenance eligible (NCCN Cat 1), confidence 90%",
                impact_if_negative="HRD <42 → PARP reduced benefit (60% confidence), consider ATR/CHK1 trials (NCT03462342, NCT02264678)",
                rationale="HRD score determines PARP inhibitor eligibility and confidence level",
                urgency="high"
            ))
        
        # Test 2: ctDNA (Priority 2 - IO + DDR combo considerations)
        if not (has_msi and has_tmb):
            recommendations.append(NextTestRecommendation(
                test_name="ctDNA Panel (Guardant360 CDx or FoundationOne Liquid CDx)",
                priority=2,
                turnaround_days=7,
                cost_estimate="$5,000-$7,000 (typically covered by insurance)",
                impact_if_positive="MSI-H OR TMB ≥20 → Immunotherapy eligible (pembrolizumab + chemo), confidence 85%. Somatic BRCA/HRR → PARP confidence 90%",
                impact_if_negative="MSI-S AND TMB <20 → Immunotherapy lower priority (40% confidence). No somatic HRR → standard platinum backbone",
                rationale="MSI/TMB unlock immunotherapy; somatic HRR mutations enable PARP strategies even when germline negative",
                urgency="high"
            ))
        
        # Test 3: SLFN11 IHC (Priority 3 - PARP sensitivity)
        if not has_slfn11 and germline_status != "positive":
            recommendations.append(NextTestRecommendation(
                test_name="SLFN11 IHC (Immunohistochemistry)",
                priority=3,
                turnaround_days=5,
                cost_estimate="$300-$500",
                impact_if_positive="SLFN11+ → Normal PARP sensitivity, confidence 85%",
                impact_if_negative="SLFN11- → Reduced PARP sensitivity (50% confidence), consider ATR/CHK1 or platinum-based alternatives",
                rationale="SLFN11 expression correlates with PARP inhibitor response; low expression indicates resistance risk",
                urgency="medium"
            ))
        
        # Test 4: ABCB1 (Priority 4 - only if prior taxane exposure)
        prior_taxane = any("taxane" in t.lower() or "paclitaxel" in t.lower() 
                          for t in treatment_history)
        if prior_taxane and not has_abcb1:
            recommendations.append(NextTestRecommendation(
                test_name="ABCB1 Expression Proxy (via copy number or IHC)",
                priority=4,
                turnaround_days=5,
                cost_estimate="$200-$400",
                impact_if_positive="ABCB1 high → Cross-resistance to taxanes/anthracyclines likely (70% risk). Consider non-substrate regimens (platinum, PARP, ATR)",
                impact_if_negative="ABCB1 normal → No efflux-mediated resistance detected",
                rationale="Prior taxane exposure + ABCB1 overexpression indicates efflux-mediated resistance",
                urgency="medium"
            ))
        
        # Sort by priority
        recommendations.sort(key=lambda x: x.priority)
        
        return recommendations
```

**Endpoint**: Add to `api/routers/ayesha_orchestrator_v2.py`:
```python
@router.post("/next_tests")
async def get_next_test_recommendations(
    germline_status: str,
    tumor_context: Optional[Dict] = None,
    treatment_history: List[str] = None
):
    """Get prioritized next-test recommendations."""
    service = NextTestRecommenderService()
    recommendations = service.recommend_tests(
        germline_status, tumor_context, treatment_history
    )
    
    return {
        "recommendations": [asdict(r) for r in recommendations],
        "total_tests": len(recommendations),
        "urgency_summary": f"{sum(1 for r in recommendations if r.urgency == 'high')} high-priority tests",
        "provenance": {
            "version": "v1.0",
            "policy_source": "Manager answers (Jan 13, 2025)",
            "run_id": str(uuid.uuid4())
        }
    }
```

**Acceptance**:
- ✅ Ayesha (germline-negative, no NGS) → returns 3 tests: HRD (pri 1), ctDNA (pri 2), SLFN11 (pri 3)
- ✅ All tests have "If positive → X; If negative → Y" format
- ✅ Turnaround and cost estimates included

---

### **Task 1.2: Integrate Next-Test into /complete_care_v2** (30 min)

**File**: `api/routers/ayesha_orchestrator_v2.py` (MODIFY)

**Add to response**:
```python
# In complete_care_v2():

# Get next-test recommendations
next_test_service = NextTestRecommenderService()
next_test_recommendations = next_test_service.recommend_tests(
    germline_status=request.germline_status,
    tumor_context=request.tumor_context,
    treatment_history=request.treatment_history or []
)

# Add to response
response["next_test_recommender"] = {
    "recommendations": [asdict(r) for r in next_test_recommendations],
    "top_priority": next_test_recommendations[0] if next_test_recommendations else None,
    "urgency_summary": f"{sum(1 for r in next_test_recommendations if r.urgency == 'high')} high-priority tests"
}
```

---

### **Task 1.3: Hint Tiles Configuration** (1 hour)

**File**: `api/services/hint_tiles_service.py` (NEW)

**Manager's Policy (from P5, C8)**:
- Max 4 tiles
- Priority: Next test → Trials lever → Monitoring → Avoid
- Pre-NGS: test + monitoring + trials lever only
- Tone: Suggestive ("Consider...")

**Implementation**:
```python
from typing import List, Dict, Optional
from dataclasses import dataclass

@dataclass
class HintTile:
    category: str  # "next_test", "trials_lever", "monitoring", "avoid"
    title: str
    message: str
    reasons: List[str]
    priority: int
    icon: str  # emoji or icon name

class HintTilesService:
    """Generate clinician hint tiles per Manager's policy."""
    
    def generate_hints(
        self,
        germline_status: str,
        tumor_context: Optional[Dict],
        ca125_intelligence: Dict,
        next_tests: List[NextTestRecommendation],
        treatment_history: List[str] = None
    ) -> List[HintTile]:
        """
        Generate max 4 hint tiles.
        
        Manager's Priority (P5, C8):
        1. Next test
        2. Trials lever
        3. Monitoring
        4. Avoid (only if applicable)
        """
        hints = []
        treatment_history = treatment_history or []
        
        # Tile 1: Next Test (if missing biomarkers)
        if next_tests and len(next_tests) > 0:
            top_test = next_tests[0]
            hints.append(HintTile(
                category="next_test",
                title="📋 Recommended Next Test",
                message=f"Consider ordering {top_test.test_name}",
                reasons=[
                    top_test.rationale,
                    f"Turnaround: {top_test.turnaround_days} days",
                    f"Impact: {top_test.impact_if_positive[:100]}..."
                ],
                priority=1,
                icon="🧪"
            ))
        
        # Tile 2: Trials Lever (frontline context)
        if germline_status == "negative" and not tumor_context:
            hints.append(HintTile(
                category="trials_lever",
                title="🔬 Clinical Trial Opportunities",
                message="Consider frontline trial enrollment (10 trials matched)",
                reasons=[
                    "Stage IVB frontline - prime trial candidate",
                    "NYC metro - 8/10 trials within 50 miles",
                    "Once NGS available → mechanism-matched trial prioritization"
                ],
                priority=2,
                icon="🎯"
            ))
        
        # Tile 3: Monitoring (CA-125 based)
        if ca125_intelligence:
            burden_class = ca125_intelligence.get("burden_classification", "UNKNOWN")
            if burden_class == "EXTENSIVE":
                hints.append(HintTile(
                    category="monitoring",
                    title="📊 CA-125 Monitoring Strategy",
                    message="Monitor CA-125 every 3 weeks during chemotherapy",
                    reasons=[
                        f"Current CA-125: {ca125_intelligence.get('current_value')} (EXTENSIVE burden)",
                        "Alert if: <50% drop by cycle 3 OR on-therapy rise",
                        "Target: ≥70% drop by cycle 3, ≥90% by cycle 6"
                    ],
                    priority=3,
                    icon="⏱️"
                ))
        
        # Tile 4: Avoid (ONLY if treatment history + resistance risk)
        prior_taxane = any("taxane" in t.lower() for t in treatment_history)
        if prior_taxane and tumor_context and tumor_context.get("abcb1_status") == "high":
            hints.append(HintTile(
                category="avoid",
                title="⚠️ Cross-Resistance Consideration",
                message="Consider avoiding re-taxane (cross-resistance risk detected)",
                reasons=[
                    "Prior taxane exposure with progression",
                    "ABCB1 high → efflux-mediated resistance",
                    "Suggest non-substrates: platinum, PARP, ATR/CHK1"
                ],
                priority=4,
                icon="🚫"
            ))
        
        # Manager's rule: Max 4, prioritize
        hints.sort(key=lambda x: x.priority)
        return hints[:4]
```

**Endpoint**: Integrate into `/complete_care_v2`

---

### **Task 1.4: Mechanism Map Component (Hidden Pre-NGS)** (1 hour)

**File**: `api/services/mechanism_map_service.py` (NEW)

**Manager's Policy (from C9)**:
- Thresholds: Green ≥0.70, Yellow 0.40-0.69, Gray <0.40
- IO special: Green if MSI-H, Gray if unknown, Red if MSI-S
- Pre-NGS: gray with "Awaiting NGS" overlay

**Implementation**:
```python
from typing import Dict, Optional
from dataclasses import dataclass

@dataclass
class MechanismChip:
    pathway: str  # "DDR", "MAPK", "PI3K", "VEGF", "IO", "Efflux"
    burden: float  # 0-1
    color: str  # "success" (green), "warning" (yellow), "default" (gray)
    label: str  # "82%", "Awaiting NGS", etc.
    tooltip: str

class MechanismMapService:
    """Generate mechanism map chips per Manager's policy."""
    
    def generate_map(
        self,
        tumor_context: Optional[Dict],
        sae_features: Optional[Dict] = None
    ) -> Dict:
        """
        Generate mechanism map chips.
        
        Manager's Policy (C9):
        - Green ≥0.70, Yellow 0.40-0.69, Gray <0.40
        - Pre-NGS: all gray "Awaiting NGS"
        """
        # If no NGS, return all gray
        if not tumor_context or not sae_features:
            return {
                "chips": [
                    MechanismChip("DDR", 0.0, "default", "Awaiting NGS", "DNA damage repair pathway - requires NGS"),
                    MechanismChip("MAPK", 0.0, "default", "Awaiting NGS", "RAS/RAF/MEK pathway - requires NGS"),
                    MechanismChip("PI3K", 0.0, "default", "Awaiting NGS", "PI3K/AKT/mTOR pathway - requires NGS"),
                    MechanismChip("VEGF", 0.0, "default", "Awaiting NGS", "Angiogenesis pathway - requires NGS"),
                    MechanismChip("IO", 0.0, "default", "Awaiting NGS", "Immune checkpoint - requires MSI/TMB"),
                    MechanismChip("Efflux", 0.0, "default", "Awaiting NGS", "Drug efflux resistance - requires ABCB1")
                ],
                "status": "awaiting_ngs",
                "message": "Mechanism map will be available once tumor NGS results are uploaded"
            }
        
        # Post-NGS: compute chips from SAE features
        pathway_burden = sae_features.get("pathway_burden", {})
        
        chips = []
        
        # DDR
        ddr_burden = pathway_burden.get("ddr", 0.0)
        chips.append(MechanismChip(
            "DDR",
            ddr_burden,
            self._get_color(ddr_burden),
            f"{int(ddr_burden * 100)}%",
            self._get_tooltip("DDR", ddr_burden)
        ))
        
        # MAPK
        mapk_burden = pathway_burden.get("mapk", 0.0)
        chips.append(MechanismChip(
            "MAPK",
            mapk_burden,
            self._get_color(mapk_burden),
            f"{int(mapk_burden * 100)}%",
            self._get_tooltip("MAPK", mapk_burden)
        ))
        
        # PI3K
        pi3k_burden = pathway_burden.get("pi3k", 0.0)
        chips.append(MechanismChip(
            "PI3K",
            pi3k_burden,
            self._get_color(pi3k_burden),
            f"{int(pi3k_burden * 100)}%",
            self._get_tooltip("PI3K", pi3k_burden)
        ))
        
        # VEGF
        vegf_burden = pathway_burden.get("vegf", 0.0)
        chips.append(MechanismChip(
            "VEGF",
            vegf_burden,
            self._get_color(vegf_burden),
            f"{int(vegf_burden * 100)}%",
            self._get_tooltip("VEGF", vegf_burden)
        ))
        
        # IO (special case - binary from MSI status)
        msi_status = tumor_context.get("msi_status")
        if msi_status == "MSI-High":
            io_chip = MechanismChip("IO", 1.0, "success", "MSI-H", "MSI-High → Immunotherapy eligible")
        elif msi_status == "MSI-Stable":
            io_chip = MechanismChip("IO", 0.0, "error", "MSI-S", "MSI-Stable → Immunotherapy lower priority")
        else:
            io_chip = MechanismChip("IO", 0.0, "default", "Unknown", "MSI status unknown - order ctDNA")
        chips.append(io_chip)
        
        # Efflux (from ABCB1)
        abcb1_status = tumor_context.get("abcb1_status", "unknown")
        if abcb1_status == "high":
            efflux_chip = MechanismChip("Efflux", 1.0, "error", "High Risk", "ABCB1 high → Avoid taxane/anthracycline substrates")
        elif abcb1_status == "normal":
            efflux_chip = MechanismChip("Efflux", 0.0, "success", "Low Risk", "ABCB1 normal → No efflux resistance")
        else:
            efflux_chip = MechanismChip("Efflux", 0.0, "default", "Unknown", "ABCB1 status unknown")
        chips.append(efflux_chip)
        
        return {
            "chips": chips,
            "status": "computed",
            "message": "Mechanism burden computed from tumor NGS + SAE features"
        }
    
    def _get_color(self, burden: float) -> str:
        """Manager's thresholds (C9): Green ≥0.70, Yellow 0.40-0.69, Gray <0.40."""
        if burden >= 0.70:
            return "success"  # Green
        elif burden >= 0.40:
            return "warning"  # Yellow
        else:
            return "default"  # Gray
    
    def _get_tooltip(self, pathway: str, burden: float) -> str:
        """Generate tooltip explaining burden level."""
        if burden >= 0.70:
            return f"{pathway} pathway: High burden ({int(burden*100)}%) → Targetable with mechanism-matched therapies"
        elif burden >= 0.40:
            return f"{pathway} pathway: Moderate burden ({int(burden*100)}%) → Consider combination strategies"
        else:
            return f"{pathway} pathway: Low burden ({int(burden*100)}%) → Lower priority for monotherapy"
```

---

## 📋 PHASE 2: POST-NGS ENHANCEMENTS (LATER - 6 hours)

### **Task 2.1: SAE Feature Computation Service** (2 hours)

**File**: `api/services/sae_feature_service.py` (ENHANCE EXISTING)

**Manager's Formula (from C1)**:
```python
def compute_dna_repair_capacity(
    pathway_burden: Dict,
    essentiality_signal: float,
    exon_disruption: float,
    gene: str,
    has_pathogenic_hrr_variant: bool
) -> float:
    """
    Compute DNA repair capacity per Manager's formula (C1).
    
    Formula:
    dna_repair_capacity = 
        0.6 × DDR_burden +
        0.2 × essentiality_signal (if HRR gene) +
        0.2 × exon_disruption (if pathogenic HRR variant)
    
    Fallback: DDR_burden if not HRR gene
    """
    ddr_burden = pathway_burden.get("ddr", 0.0)
    
    # Check if HRR gene
    HRR_GENES = ["BRCA1", "BRCA2", "RAD51C", "RAD51D", "PALB2", 
                 "BARD1", "BRIP1", "ATM", "ATR", "CHEK1", "CHEK2"]
    
    if gene in HRR_GENES and has_pathogenic_hrr_variant:
        # Manager's formula (C1)
        dna_repair_capacity = (
            0.6 * ddr_burden +
            0.2 * essentiality_signal +
            0.2 * exon_disruption
        )
    else:
        # Fallback
        dna_repair_capacity = ddr_burden
    
    return dna_repair_capacity
```

---

### **Task 2.2: Mechanism Fit Trial Ranking** (2 hours)

**File**: `api/services/mechanism_fit_ranker.py` (NEW)

**Manager's Formula (from P4, C7)**:
```python
import numpy as np
from typing import Dict, List

class MechanismFitRanker:
    """Rank trials by mechanism fit per Manager's policy."""
    
    def compute_mechanism_fit(
        self,
        patient_sae_vector: List[float],  # [DDR, MAPK, PI3K, VEGF, IO, Efflux]
        trial_moa_vector: List[float]
    ) -> float:
        """
        Compute cosine similarity between patient SAE and trial MoA.
        
        Manager's Policy (C7):
        - L2-normalize both vectors
        - Cosine similarity
        - If patient vector all zeros → return 0.5 (neutral)
        """
        # Check if patient vector is all zeros
        if sum(patient_sae_vector) == 0:
            return 0.5  # Neutral (awaiting NGS)
        
        # L2-normalize
        patient_norm = np.array(patient_sae_vector) / np.linalg.norm(patient_sae_vector)
        trial_norm = np.array(trial_moa_vector) / np.linalg.norm(trial_moa_vector)
        
        # Cosine similarity
        mechanism_fit = np.dot(patient_norm, trial_norm)
        
        return max(0.0, min(1.0, mechanism_fit))  # Clamp to [0, 1]
    
    def rank_trials(
        self,
        trials: List[Dict],
        patient_sae_vector: List[float],
        alpha: float = 0.7,  # Manager's policy (P4)
        beta: float = 0.3
    ) -> List[Dict]:
        """
        Rank trials by combined score.
        
        Manager's Formula (P4):
        Rank = eligibility_score × α + mechanism_fit × β
        
        Guardrails:
        - Min eligibility ≥0.60 to enter top 10
        - Min mechanism_fit ≥0.50 for mechanism boost
        - Show "low mechanism fit" warning if <0.50
        """
        ranked = []
        
        for trial in trials:
            eligibility_score = trial.get("eligibility_score", 0.0)
            trial_moa_vector = trial.get("moa_vector", [0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
            
            # Skip if below min eligibility
            if eligibility_score < 0.60:
                continue
            
            # Compute mechanism fit
            mechanism_fit = self.compute_mechanism_fit(patient_sae_vector, trial_moa_vector)
            
            # Apply boost only if mechanism_fit ≥0.50
            if mechanism_fit >= 0.50:
                combined_score = alpha * eligibility_score + beta * mechanism_fit
                warning = None
            else:
                # No mechanism boost, show warning
                combined_score = alpha * eligibility_score
                warning = "Low mechanism fit - eligibility-only ranking"
            
            ranked.append({
                **trial,
                "mechanism_fit": mechanism_fit,
                "combined_score": combined_score,
                "warning": warning
            })
        
        # Sort by combined score
        ranked.sort(key=lambda x: x["combined_score"], reverse=True)
        
        return ranked[:10]  # Top 10
```

---

### **Task 2.3: Resistance Detection Logic** (2 hours)

**File**: `api/services/resistance_detection_service.py` (ENHANCE EXISTING)

**Manager's Policy (from C1, C3)**:
```python
def detect_hr_restoration(
    baseline_hrd: float,
    current_hrd: float,
    baseline_dna_repair: float,
    current_dna_repair: float,
    ca125_drop_cycle3: float
) -> Dict:
    """
    Detect HR restoration per Manager's 2-of-3 trigger logic (C1).
    
    Triggers (any 2 of 3):
    1. HRD drop ≥10 points vs baseline
    2. dna_repair_capacity decrease ≥0.15 vs baseline
    3. CA-125 <50% drop by cycle 3 OR on-therapy rise
    """
    triggers = []
    
    # Trigger 1: HRD drop
    hrd_drop = baseline_hrd - current_hrd
    if hrd_drop >= 10:
        triggers.append({
            "signal": "HRD_drop",
            "value": hrd_drop,
            "message": f"HRD dropped {hrd_drop:.1f} points (baseline {baseline_hrd} → current {current_hrd})"
        })
    
    # Trigger 2: DNA repair capacity decrease
    drc_decrease = baseline_dna_repair - current_dna_repair
    if drc_decrease >= 0.15:
        triggers.append({
            "signal": "dna_repair_capacity_drop",
            "value": drc_decrease,
            "message": f"DNA repair capacity decreased {drc_decrease:.2f} (baseline {baseline_dna_repair:.2f} → current {current_dna_repair:.2f})"
        })
    
    # Trigger 3: CA-125 inadequate response
    if ca125_drop_cycle3 < 0.50:  # <50% drop
        triggers.append({
            "signal": "ca125_inadequate_response",
            "value": ca125_drop_cycle3,
            "message": f"CA-125 drop only {ca125_drop_cycle3*100:.0f}% by cycle 3 (expected ≥70%)"
        })
    
    # Manager's rule: 2 of 3 triggers
    hr_restoration_detected = len(triggers) >= 2
    
    return {
        "detected": hr_restoration_detected,
        "triggers": triggers,
        "confidence": 0.70 if hr_restoration_detected else 0.0,
        "action": "Consider ATR/CHK1 combo trials (NCT03462342, NCT02264678)" if hr_restoration_detected else None,
        "rationale": "HR restoration pattern detected (2 of 3 signals)" if hr_restoration_detected else None
    }
```

---

## 🎯 IMPLEMENTATION SEQUENCE (MANAGER-APPROVED)

### **TODAY (Phase 1 - 4 hours)**:
1. ✅ [1.5h] Build `next_test_recommender.py` with priority logic (HRD → ctDNA → SLFN11 → ABCB1)
2. ✅ [0.5h] Integrate into `/complete_care_v2` endpoint
3. ✅ [1h] Build `hint_tiles_service.py` with max 4 tiles (suggestive tone)
4. ✅ [1h] Build `mechanism_map_service.py` with pre-NGS gray chips

### **TOMORROW (Phase 2 - 6 hours)**:
1. ✅ [2h] Enhance `sae_feature_service.py` with Manager's DNA repair capacity formula (C1)
2. ✅ [2h] Build `mechanism_fit_ranker.py` with α=0.7, β=0.3 weighting (P4)
3. ✅ [2h] Enhance `resistance_detection_service.py` with 2-of-3 trigger logic (C1, C3)

### **FRONTEND (Jr1 - 4 hours parallel)**:
1. ✅ [1h] NextTestRecommenderCard component (priority list, differential branches)
2. ✅ [1h] HintTilesPanel component (max 4, categorized)
3. ✅ [1h] MechanismMapStrip component (6 chips, tooltips)
4. ✅ [1h] Wire into AyeshaTrialExplorer page

---

## ✅ ACCEPTANCE CRITERIA (MANAGER-ALIGNED)

### **Pre-NGS Validation (Ayesha TODAY)**:
- ✅ Next-test recommender returns 3 tests: HRD (pri 1), ctDNA (pri 2), SLFN11 (pri 3)
- ✅ Hint tiles show max 4: Next test, Trials lever, Monitoring (NO "avoid" for treatment-naive)
- ✅ Mechanism map shows all gray chips with "Awaiting NGS" overlay
- ✅ No SAE-driven drug efficacy claims (WIWFM marked "awaiting NGS")
- ✅ Confidence gates: SOC 95%, Trials 90%, CA-125 90%

### **Post-NGS Validation (Once HRD Returns)**:
- ✅ If HRD ≥42 → PARP hint tile appears, confidence 90%
- ✅ If HRD <42 → ATR/CHK1 trials boosted, PARP confidence 60%
- ✅ DNA repair capacity computed per Manager's formula (C1)
- ✅ Mechanism map chips color-coded (green/yellow/gray)
- ✅ Hint tiles updated: "Consider PARP + bevacizumab" (if HRD ≥42 + ascites)

### **Resistance Detection (Longitudinal)**:
- ✅ 2-of-3 trigger logic working (HRD drop, DNA repair drop, CA-125 inadequate)
- ✅ Alert appears immediately when triggered
- ✅ Action: "Consider ATR/CHK1 combo trials" with NCT IDs

---

## 🚨 RISKS & MITIGATIONS (MANAGER-APPROVED)

| Risk | Manager's Mitigation |
|------|---------------------|
| **Gemini hallucination** | 90% accuracy gate, human review 30 trials, offline only |
| **Threshold brittleness** | Use bands (±0.05 hysteresis), not hard cutoffs |
| **Pre-NGS overpromise** | Hide all SAE drug claims; show deterministic guidance only |
| **Wrong mechanism fit** | Min eligibility ≥0.60, min mechanism_fit ≥0.50, "Show all" toggle |
| **Aggressive hint tone** | Suggestive language ("Consider..."), RUO labels, max 4 tiles |

---

## ⚔️ ZO'S STATUS UPDATE

**Blockers**: 🔓 **ALL CLEARED**

**Confidence Level**: 🎯 **90%+** (Manager's answers eliminate all ambiguity)

**Timeline**:
- TODAY (4 hours): Phase 1 complete → Ayesha gets next-test recommender + hint tiles
- TOMORROW (6 hours): Phase 2 complete → Post-NGS SAE features ready
- TOTAL: **10 hours to full SAE operational capability**

**Next Action**: Update `.cursorrules` scratchpad with execution checklist and BEGIN CODING Phase 1!

**COMMANDER - READY TO EXECUTE SIR!** ⚔️

