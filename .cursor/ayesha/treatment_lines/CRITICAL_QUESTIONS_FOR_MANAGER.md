# ⚔️ CRITICAL QUESTIONS FOR MANAGER - TREATMENT LINE INTEGRATION

**Agent**: Zo  
**Mission**: Full treatment line integration (P0→P1→P2→P3)  
**Timeline**: 8-12 hours  
**Status**: Awaiting clarification before execution

---

## 🎯 STRATEGIC QUESTIONS


## ✅ MANAGER ANSWERS (APPROVED FOR EXECUTION)

**Q1 (Disease Panels)**: A — Ovarian + Breast only. Focus on Ovarian L2 + HER2 playbook first.  
**Q2 (Treatment History)**: A — User manual input for P0 (structured fields).  
**Q3 (Cross-Resistance)**: A now; add C (genomic resistance markers) in P1.  
**Q4 (NCCN Enforcement)**: C (display only, RUO); add D (user filter) in P1.  
**Q5 (Prior Outcomes)**: C for P0; add A + D (PFS duration + genomic resistance) in P1.  
**Q6 (SAE Names)**: C — Hybrid labels: "Treatment Line Fit", "Resistance Risk", "Sequencing Score" with plain-language tooltips.  
**Q7 (Line Definition)**: A + D — Strict NCCN definition with user override.  
**Q8 (Confidence Penalty)**: A — Linear: `confidence -= cross_resistance_risk × 0.2` (cap -20%); add D in P1.  
**Q9 (Frontend Form)**: B for P0 (structured dropdowns, multi-select); C in P1 (timeline builder).  
**Q10 (Testing)**: D — All six smoke tests (ovarian L1/L2, breast L1/L2/L3, edge L4).

---

## 🤔 FOLLOW-UP CLARIFICATION NEEDED

### **Cross-Resistance MoA Mapping Strategy**

**Context**: Manager approved "A - Simple MoA overlap" for Q3, but I need to clarify implementation approach.

**Example Scenario**:
- AK receives "carboplatin + paclitaxel" (1st line)
- Now considering "olaparib" (PARP inhibitor, 2nd line)
- Both affect DNA repair → cross-resistance risk

**Three Implementation Options**:

#### **Option A: Broad MoA Categories**
```python
DRUG_MOA_MAP = {
    "carboplatin": "DNA_crosslinking",
    "paclitaxel": "microtubule_disruption",
    "olaparib": "DNA_repair_blockade"
}

# Check if MoAs in same category
if "DNA" in prior_moa and "DNA" in candidate_moa:
    cross_resistance_risk = 0.4
```
**Pros**: Simple, fast  
**Cons**: Too broad, false positives

---

#### **Option B: MoA Tag System**
```python
DRUG_MOA_TAGS = {
    "carboplatin": ["DNA_DAMAGE", "DNA_REPAIR_PATHWAY", "platinum_agent"],
    "olaparib": ["DNA_REPAIR_BLOCKADE", "DNA_REPAIR_PATHWAY", "PARP_inhibitor"],
    "trastuzumab": ["HER2_BLOCKADE", "antibody"]
}

# Calculate overlap
prior_tags = get_moa_tags("carboplatin")
candidate_tags = get_moa_tags("olaparib")
overlap = len(set(prior_tags) & set(candidate_tags))
cross_resistance_risk = min(overlap * 0.2, 1.0)  # 1 tag = 0.2, 2 tags = 0.4
```
**Pros**: Scalable, flexible  
**Cons**: Requires tag maintenance

---

#### **Option C: Drug-Specific Cross-Resistance Map** ⭐ (RECOMMENDED)
```python
CROSS_RESISTANCE_MAP = {
    # Ovarian cancer
    "PARP_inhibitor": {
        "cross_resistant_with": ["platinum_agent"],
        "risk_level": 0.4,
        "rationale": "DNA repair pathway overlap"
    },
    "platinum_agent": {
        "cross_resistant_with": ["PARP_inhibitor"],
        "risk_level": 0.4,
        "rationale": "DNA repair pathway overlap"
    },
    
    # Breast HER2+ cancer
    "T-DXd": {
        "cross_resistant_with": ["trastuzumab", "pertuzumab"],
        "risk_level": 0.3,
        "rationale": "HER2-targeted mechanism overlap"
    },
    "tucatinib": {
        "cross_resistant_with": ["trastuzumab", "T-DXd", "lapatinib"],
        "risk_level": 0.2,
        "rationale": "HER2 TKI resistance after prior HER2 blockade"
    }
}

# Usage
def get_cross_resistance_risk(prior_drug: str, candidate_drug: str) -> float:
    candidate_class = get_drug_class(candidate_drug)  # e.g., "PARP_inhibitor"
    prior_class = get_drug_class(prior_drug)  # e.g., "platinum_agent"
    
    cross_res_info = CROSS_RESISTANCE_MAP.get(candidate_class, {})
    cross_resistant_with = cross_res_info.get("cross_resistant_with", [])
    
    if prior_class in cross_resistant_with:
        return cross_res_info.get("risk_level", 0.4)
    
    return 0.0  # No known cross-resistance
```
**Pros**: Clinically accurate, explicit, easy to validate  
**Cons**: Manual curation per drug

---

### **My Recommendation: Option C for P0**

**Rationale**:
1. **Clinical Accuracy**: Matches real oncology practice (doctors know PARP/platinum cross-resistance)
2. **Explainability**: Clear rationale for each cross-resistance relationship
3. **Fast to Build**: Only need 10 drugs (5 ovarian + 5 breast)
4. **Easy to Test**: Deterministic, no ambiguity
5. **Upgrade Path**: Can migrate to Option B (tags) in P1 when we add more drugs

**Implementation for Ovarian L2 Case**:
```python
# Ovarian L2 case: Prior platinum, considering PARP
prior = "carboplatin"  # maps to "platinum_agent"
candidate = "olaparib"  # maps to "PARP_inhibitor"

cross_resistance_risk = get_cross_resistance_risk(prior, candidate)
# Returns: 0.4 (from CROSS_RESISTANCE_MAP)

# SAE Feature
{
    "id": "cross_resistance_risk",
    "activation": 0.4,
    "impact": "negative",
    "explanation": "Moderate cross-resistance risk with prior platinum (DNA repair pathway overlap)",
    "rationale": "Both platinum and PARP target DNA repair mechanisms"
}
```

---

## ❓ MANAGER DECISION REQUIRED

**Question**: Should I proceed with **Option C (drug-specific map)** for P0?

**If YES**: I will start Phase 1 immediately (building ovarian + breast panels with cross-resistance map)

**If NO**: Please specify:
- Use Option A (broad categories)?
- Use Option B (tag system)?
- Use different approach?

---

⚔️💀 **AWAITING FINAL CONFIRMATION TO BEGIN EXECUTION!** 💀⚔️

### **Q1: Disease Panels Priority** 🤔

**Context**: We need to build disease-specific drug panels with treatment line metadata.

**Question**: Which disease panels should I prioritize?

**Options**:
- **A) Ovarian + Breast only** (Ovarian L2 case + Dr. Lustberg's work)
- **B) Ovarian + Breast + Lung** (top 3 solid tumors)
- **C) Ovarian + Breast + Multiple Myeloma** (leverage existing MM work)
- **D) Other priority?**

**My Recommendation**: **A) Ovarian + Breast only** (ship fast, validate with Ovarian L2 case, expand later)

---

### **Q2: Treatment History Input Source** 🤔

**Context**: We need prior therapy data to compute cross-resistance and line appropriateness.

**Question**: Where does treatment history come from?

**Options**:
- **A) User manual input** (doctor types in prior therapies)
- **B) EHR integration** (parse from patient record - future)
- **C) Inferred from progression timeline** (e.g., if variant seen at relapse, assume prior platinum)
- **D) Clinical trial data** (extract from trial seeding)

**My Recommendation**: **A) User manual input for P0** (simple form: "Prior therapies: platinum + taxane, 8 months PFS, partial response")

---

### **Q3: Cross-Resistance Modeling Depth** 🤔

**Context**: Platinum and PARP have DNA repair MoA overlap → cross-resistance risk.

**Question**: How deep should cross-resistance modeling go?

**Options**:
- **A) Simple MoA overlap** (if prior MoA matches drug MoA → flag risk)
- **B) Pathway-based overlap** (if prior pathway weights overlap ≥50% → flag risk)
- **C) Genomic resistance markers** (if BRCA reversion mutation → mark PARP-resistant)
- **D) Time-dependent decay** (if prior therapy >12 months ago, reduce risk)

**My Recommendation**: **A) Simple MoA overlap for P0** (fast, interpretable), **C) Genomic resistance markers in P1** (leverage existing variant analysis)

---

### **Q4: NCCN Category Enforcement** 🤔

**Context**: NCCN categories define evidence strength (1 = uniform consensus, 3 = disagreement).

**Question**: Should we enforce NCCN categories or just surface them?

**Options**:
- **A) Enforce (hard filter)** - Hide drugs with Category 3 or no NCCN support
- **B) Soft filter (rank lower)** - Show all drugs, rank NCCN Category 1 higher
- **C) Display only (no filtering)** - Show NCCN category as badge, let doctor decide
- **D) User preference** - Toggle "Show only NCCN-supported" checkbox

**My Recommendation**: **C) Display only for RUO** (we're research-use, not prescriptive), **D) User preference in production** (let doctors control filtering)

---

### **Q5: Prior Therapy Outcomes Impact** 🤔

**Context**: If a patient had 8 months PFS on platinum (partial response), should this affect PARP confidence?

**Question**: Should we factor prior therapy outcomes into confidence?

**Options**:
- **A) Yes - PFS duration matters** (if prior PFS <6mo → reduce confidence in similar MoA)
- **B) Yes - Response type matters** (if prior = complete response → increase confidence in similar pathway)
- **C) No - Too complex for P0** (just track cross-resistance MoA overlap)
- **D) Conditional - Only for genomic resistance** (if progression mutation detected → flag resistance)

**My Recommendation**: **C) No for P0** (keep simple), **A) + D) in P1** (add PFS duration + genomic resistance when we have that data)

---

### **Q6: SAE Feature Naming** 🤔

**Context**: New SAE features need doctor-friendly names.

**Question**: What should we call these features in the UI?

**Current Technical Names**:
- `line_appropriateness`
- `cross_resistance_risk`
- `sequencing_fitness`

**UI Name Options**:
- **A) Medical terminology**: "Line Appropriateness", "Cross-Resistance Risk", "Sequencing Fitness"
- **B) Plain language**: "Right for This Stage?", "Prior Therapy Concerns", "Best Next Step Score"
- **C) Hybrid**: "Treatment Line Fit", "Resistance Risk", "Sequencing Score"
- **D) Emoji-enhanced**: "✅ Line Fit", "⚠️ Resistance Risk", "🎯 Sequencing Score"

**My Recommendation**: **C) Hybrid** (professional but clear), **B) Plain language in tooltips** (explain in hover text)

---

### **Q7: Treatment Line Definition** 🤔

**Context**: What counts as a "line of therapy"?

**Question**: How do we define treatment lines?

**Options**:
- **A) Strict NCCN definition** (first-line = upfront, second-line = first progression, etc.)
- **B) Consolidation/maintenance separate** (e.g., PARP maintenance after platinum = still line 1)
- **C) Intent-based** (curative vs palliative affects line counting)
- **D) User-defined** (doctor specifies line number)

**My Recommendation**: **A) Strict NCCN for consistency**, **D) User override** (let doctor correct if system guesses wrong)

---

### **Q8: Confidence Penalty for Cross-Resistance** 🤔

**Context**: If cross-resistance risk = 0.4, how much should we penalize confidence?

**Question**: What's the penalty formula?

**Options**:
- **A) Linear penalty**: `confidence -= cross_resistance_risk * 0.2` (max -20% confidence)
- **B) Exponential penalty**: `confidence *= (1 - cross_resistance_risk * 0.3)` (proportional reduction)
- **C) Threshold-based**: Only penalize if risk ≥0.5 (ignore minor overlaps)
- **D) Evidence-modulated**: Penalty reduced if strong evidence overrides resistance concern

**My Recommendation**: **A) Linear penalty for P0** (simple, transparent), **D) Evidence-modulated in P1** (RCT evidence can override resistance concerns)

---

### **Q9: Frontend Treatment History Form** 🤔

**Context**: Doctor needs to input prior therapies.

**Question**: What's the UI for treatment history input?

**Options**:
- **A) Free text** (doctor types "platinum + taxane, 8mo PFS, PR")
- **B) Structured form** (dropdowns for drug, PFS, response type)
- **C) Timeline builder** (visual timeline with drag-drop therapy cards)
- **D) EHR paste** (paste structured data from EHR, we parse)

**My Recommendation**: **B) Structured form for P0** (fast to build, data quality), **C) Timeline builder in P1** (better UX for complex cases)

---

### **Q10: Testing Strategy** 🤔

**Context**: We need to validate treatment line logic end-to-end.

**Question**: What test cases should I prioritize?

**Options**:
- **A) Ovarian L2 case** (ovarian, BRCA2, 2nd line post-platinum)
- **B) Dr. Lustberg's case** (breast, HER2+, 3rd line post-T-DXd)
- **C) Edge cases** (1st line vs 4th line, conflicting prior therapies)
- **D) All of the above**

**My Recommendation**: **D) All of the above** - write 6 comprehensive smoke tests:
1. Ovarian 1st line (platinum preference)
2. Ovarian 2nd line post-platinum (PARP appropriateness + cross-resistance)
3. Breast HER2+ 1st line (T/P/taxane preference)
4. Breast HER2+ 2nd line post-trastuzumab (T-DXd appropriateness)
5. Breast HER2+ 3rd line post-T-DXd (tucatinib sequencing)
6. Edge case: 4th line (all options exhausted, show trials)

---

## 📋 MANAGER RESPONSE TEMPLATE

**Please answer each question with A/B/C/D or provide guidance:**

```
Q1 (Disease Panels): A
Q2 (Treatment History): A
Q3 (Cross-Resistance): A + C (P1)
Q4 (NCCN Enforcement): C + D (P1)
Q5 (Prior Outcomes): C (P0) → A + D (P1)
Q6 (SAE Names): C (Hybrid)
Q7 (Line Definition): A + D
Q8 (Confidence Penalty): A (P0) → D (P1)
Q9 (Frontend Form): B (P0) → C (P1)
Q10 (Testing): D (All 6 tests)
```

## ✅ AGENT X IMPLEMENTATION ANSWERS (FINAL)

### Scope and Strategy (P0)
- Disease panels: Ovarian + Breast (HER2+) only
- Treatment history: Structured user input (current_line, prior_therapies[], optional outcomes)
- Cross‑resistance: Option C (drug‑specific map) for P0; add genomic markers in P1
- NCCN: Display only (RUO); add user filter in P1
- Prior outcomes: Skip in P0; add PFS duration + genomic resistance in P1
- SAE labels: "Treatment Line Fit", "Resistance Risk", "Sequencing Score" (tooltips in plain language)
- Line definition: Strict NCCN with user override
- Confidence penalty: Linear −0.2 × cross_resistance_risk (cap −0.2); evidence‑modulated in P1
- Frontend: Structured form now; timeline builder in P1
- Testing: Implement all 6 smoke tests listed

### Drug Classes (for mapping and sequencing)
- Ovarian: `platinum_agent`, `PARP_inhibitor`, `bevacizumab_combo`, `topotecan`, `trial_flag`
- Breast HER2+: `TP_taxane_combo` (trastuzumab+pertuzumab+taxane), `T-DXd`, `tucatinib_combo`, `neratinib_combo`, `trial_flag`

### get_drug_class(drug_name) examples
```python
DRUG_CLASS_MAP = {
  # Ovarian
  "carboplatin": "platinum_agent",
  "cisplatin": "platinum_agent",
  "olaparib": "PARP_inhibitor",
  "niraparib": "PARP_inhibitor",
  "rucaparib": "PARP_inhibitor",
  "bevacizumab": "bevacizumab_combo",  # treat combos that include bevacizumab as this class
  "topotecan": "topotecan",

  # Breast HER2+
  "trastuzumab+pertuzumab+taxane": "TP_taxane_combo",
  "trastuzumab deruxtecan": "T-DXd",
  "tdxd": "T-DXd",
  "tucatinib+trastuzumab+capecitabine": "tucatinib_combo",
  "neratinib+capecitabine": "neratinib_combo",
}
```

### Cross‑resistance map (P0)
```python
CROSS_RESISTANCE_MAP = {
  # Ovarian
  "PARP_inhibitor": {"cross_resistant_with": ["platinum_agent"], "risk_level": 0.4,
                      "rationale": "DNA repair pathway overlap"},
  "platinum_agent": {"cross_resistant_with": ["PARP_inhibitor"], "risk_level": 0.4,
                      "rationale": "DNA repair pathway overlap"},

  # Breast HER2+
  "T-DXd": {"cross_resistant_with": ["trastuzumab", "pertuzumab"], "risk_level": 0.3,
             "rationale": "HER2-targeted mechanism overlap"},
  "tucatinib_combo": {"cross_resistant_with": ["trastuzumab", "T-DXd", "lapatinib"], "risk_level": 0.2,
                       "rationale": "HER2 TKI resistance after prior HER2 blockade"},
}
```

### Schema and file placement
- Panels with treatment line metadata: `oncology-backend-minimal/api/services/pathway/panel_config.py`
- Cross‑resistance map: `oncology-backend-minimal/api/services/pathway/cross_resistance_map.py`
- Request schema (add treatment history): `oncology-backend-minimal/api/services/efficacy_orchestrator/models.py`
- SAE integration: `oncology-backend-minimal/api/services/sae_service.py`
- Confidence integration: `oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py`

### Treatment history model (P0)
```python
class TreatmentHistory(BaseModel):
  current_line: int
  prior_therapies: List[str] = []  # raw drug or regimen strings mapped via DRUG_CLASS_MAP
  outcomes: Optional[List[Dict[str, Any]]] = None  # optional in P0
```

### Provenance additions
```json
provenance.treatment_line = {
  "current_line": 2,
  "prior_therapies": ["carboplatin"],
  "line_appropriateness": 1.0,
  "cross_resistance_risk": 0.4,
  "sequencing_fitness": 0.85,
  "nccn_category": "1"
}
```

### Defaults and guardrails
- If treatment history missing: `line_appropriateness = 0.6`, `cross_resistance_risk = 0.0`, show "insufficient history" chip
- If drug not in class map: fall back to pathway MoA tags in P1; P0 treat as neutral

### Acceptance (P0)
- Ovarian L2 post‑platinum: PARP → `line_fit=1.0`, `cross_res≈0.4`, `seq_fitness≥0.8`; confidence higher than same case L1
- Breast HER2+ L3 post‑T‑DXd: tucatinib_combo → `line_fit=1.0`, `cross_res≈0.2`; confidence appropriate for line
- EvidenceBand shows treatment line breakdown + NCCN badge; SAE shows 3 new features

### Immediate execution checklist (P0)
1) Build ovarian + breast panels with `treatment_line_metadata`
2) Add `cross_resistance_map.py` and wire into SAE + confidence
3) Extend `EfficacyRequest` with `treatment_history`
4) Update EvidenceBand + SAE card rendering for new fields
5) Implement 6 smoke tests and record before/after deltas

---
---

## 🚀 EXECUTION PLAN (After Answers)

**Phase 0: Setup** (30 min)
- Create folder structure
- Write schemas and types
- Set up test fixtures

**Phase 1: Backend Drug Panels** (2-3h)
- Build ovarian panel (5 drugs)
- Build breast panel (5 drugs)
- Add treatment_line_metadata
- Write panel loader tests

**Phase 2: SAE Treatment Line Features** (3-4h)
- Add 3 new SAE features
- Integrate into sae_service.py
- Write feature computation tests

**Phase 3: Confidence Integration** (1-2h)
- Update orchestrator with treatment history
- Modulate confidence with line fitness
- Add provenance tracking

**Phase 4: Frontend Display** (2-3h)
- Treatment history input form
- SAE features display
- CoPilot actions

**Phase 5: Testing & Documentation** (1-2h)
- Write 6 smoke tests
- Update treatment_line_plan.mdc
- Generate final completion report

---

⚔️💀 **AWAITING MANAGER'S ANSWERS TO PROCEED!** 💀⚔️

