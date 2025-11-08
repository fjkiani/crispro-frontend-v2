# LLM-Enhanced Food Validator Strategy

## 🎯 CURRENT STATE vs. ENHANCED STATE

### **Current (Hardcoded A→B):**
✅ **Strengths:**
- Fast (<1s response)
- Reproducible (same input = same output)
- Transparent (A→B mappings are auditable)
- No API costs

❌ **Limitations:**
- Static knowledge (only 6 compounds)
- No real-time literature updates
- Can't handle new compounds without manual curation
- Evidence summaries are static

### **Enhanced (LLM + PubMed):**
✅ **New Capabilities:**
- **Dynamic compound discovery**: Ask about ANY food/supplement
- **Real-time literature mining**: Latest PubMed research
- **Personalized citations**: Actual papers for each A→B link
- **Evidence strength scoring**: LLM evaluates study quality
- **Bioavailability insights**: LLM extracts from pharmacology papers
- **Safety warnings**: LLM identifies contraindications from literature

---

## 🔬 HYBRID ARCHITECTURE (Best of Both Worlds)

```
User Query: "Can curcumin help ovarian cancer?"
    ↓
[1] FAST PATH (Hardcoded Cache) - <1s
    ↓
    Check: Is "curcumin" in our validated_claims.json?
    ↓
    YES → Return cached A→B + evidence summary
    ↓
[2] LLM ENHANCEMENT PATH (If time permits) - Background job
    ↓
    [A] PubMed Literature Mining (5-10s)
        - Query: "curcumin ovarian cancer NF-kappa-B"
        - Extract: 20 most relevant papers (LLM reranking)
    ↓
    [B] Evidence Synthesis (5-10s)
        - LLM reads abstracts
        - Extracts: mechanisms, trial outcomes, doses, safety
        - Generates: updated evidence summary + citations
    ↓
    [C] A→B Validation (5-10s)
        - LLM validates each A→B link with literature
        - Scores confidence per link
        - Identifies new A→B links not in our database
    ↓
    [D] Cache Update (optional)
        - Store enhanced results for 7 days
        - Next user gets enhanced version instantly
    ↓
[3] RETURN Enhanced Result
    - Fast path result (instant)
    - + LLM enhancements (background, displayed when ready)
```

---

## 📋 IMPLEMENTATION PLAN

### **Phase 1: PubMed Literature Mining (Tonight - 2 hours)**

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/pubmed_food_validator.py`

```python
from Pubmed-LLM-Agent-main.core.llm_client import LLMClient
from Pubmed-LLM-Agent-main.core.llm_rerank import rerank_records_with_llm_batched
from Pubmed-LLM-Agent-main.pubmed_llm_agent import PubMedLLMAgent

async def enhance_food_validation_with_literature(
    compound: str,
    disease: str = "ovarian_cancer_hgs",
    A_alterations: List[Dict] = None
) -> Dict:
    """
    Use PubMed + LLM to:
    1. Find latest literature on compound + disease
    2. Extract mechanisms (A→B links)
    3. Evaluate evidence quality
    4. Generate citations
    
    Returns enhanced validation result with:
    - literature_evidence: List[Dict] (PMIDs, titles, relevance scores)
    - llm_extracted_mechanisms: List[str]
    - llm_ab_links: List[Dict] (A→B with confidence + citations)
    - evidence_grade_llm: str (STRONG/MODERATE/WEAK/INSUFFICIENT)
    """
    
    # 1. Build PubMed query
    # Example: "curcumin AND ovarian cancer AND (NF-kappa-B OR COX-2 OR inflammation)"
    query_terms = [compound, "ovarian cancer"]
    if A_alterations:
        # Add A-specific terms
        for A in A_alterations:
            if "TP53" in A['A']:
                query_terms.append("(TP53 OR p53)")
            if "inflammation" in str(A['pathways_disrupted']).lower():
                query_terms.append("(NF-kappa-B OR inflammation OR IL-6)")
    
    query = " AND ".join(query_terms)
    
    # 2. Search PubMed (use existing agent)
    agent = PubMedLLMAgent()
    papers = agent.search_and_rank(query, max_results=20)
    
    # 3. LLM extracts mechanisms
    llm = LLMClient(model="gemini-2.5-flash")
    
    system_prompt = f"""
You are a clinical oncology expert. Extract therapeutic mechanisms from these papers about {compound} in ovarian cancer.

For each mechanism:
1. Identify the molecular target (e.g., NF-κB, COX-2, TP53)
2. Describe the mechanism (how does {compound} affect this target?)
3. Cite the PMID(s) supporting this mechanism
4. Rate evidence strength (STRONG/MODERATE/WEAK)

Return JSON:
{{
  "mechanisms": [
    {{
      "target": "NF-κB",
      "mechanism": "Curcumin inhibits NF-κB nuclear translocation...",
      "pmids": ["12345678", "87654321"],
      "evidence_strength": "MODERATE"
    }}
  ],
  "overall_evidence_grade": "MODERATE",
  "bioavailability_notes": "Poor oral bioavailability (<5%); requires liposomal formulation",
  "recommended_dose": "500-1000mg with piperine",
  "safety_concerns": "Avoid with anticoagulants"
}}
"""
    
    user_prompt = f"Papers:\n"
    for paper in papers[:10]:  # Top 10
        user_prompt += f"""
PMID: {paper['pmid']}
Title: {paper['title']}
Abstract: {paper['abstract'][:800]}
---
"""
    
    llm_result = llm.complete_json(system_prompt, user_prompt, max_tokens=2000)
    
    # 4. Map LLM mechanisms to our A→B framework
    ab_links_llm = []
    for mech in llm_result.get('mechanisms', []):
        target = mech['target']
        mechanism_text = mech['mechanism']
        pmids = mech['pmids']
        strength = mech['evidence_strength']
        
        # Try to match to our A alterations
        matched_A = None
        for A in A_alterations:
            # Check if LLM target matches any B dependency
            for B_dep in A['B_dependencies']:
                if target.lower() in B_dep['B'].lower():
                    matched_A = A
                    break
        
        ab_links_llm.append({
            "A": matched_A['A'] if matched_A else "Unknown tumor alteration",
            "B": target,
            "mechanism": mechanism_text,
            "pmids": pmids,
            "evidence_strength": strength,
            "llm_validated": True
        })
    
    return {
        "literature_evidence": papers[:20],  # Top 20 papers
        "llm_extracted_mechanisms": [m['mechanism'] for m in llm_result.get('mechanisms', [])],
        "llm_ab_links": ab_links_llm,
        "evidence_grade_llm": llm_result.get('overall_evidence_grade', 'INSUFFICIENT'),
        "bioavailability_llm": llm_result.get('bioavailability_notes', ''),
        "dosage_llm": llm_result.get('recommended_dose', ''),
        "safety_llm": llm_result.get('safety_concerns', ''),
        "total_papers_found": len(papers)
    }
```

### **Phase 2: Update hypothesis_validator.py Endpoint (30 min)**

```python
@router.post("/api/hypothesis/validate_food_ab_enhanced")
async def validate_food_ab_enhanced(
    compound: str,
    disease: str = "ovarian_cancer_hgs",
    germline_status: str = "negative",
    treatment_line: Optional[int] = None,
    prior_therapies: Optional[List[str]] = None,
    use_llm: bool = True  # Toggle for LLM enhancement
):
    """
    Enhanced A→B Food Validator with LLM literature mining.
    
    Response includes:
    - Fast path: Hardcoded A→B result (instant)
    - Enhanced: LLM-mined literature + mechanisms (10-20s)
    """
    
    # 1. Fast path (existing logic)
    base_result = await validate_food_ab(
        compound, disease, germline_status, treatment_line, prior_therapies
    )
    
    # 2. LLM enhancement (if requested)
    if use_llm and base_result['status'] == 'SUCCESS':
        disease_data = DISEASE_AB[disease]
        A_alterations = disease_data['A_alterations']
        
        try:
            llm_enhancement = await enhance_food_validation_with_literature(
                compound, disease, A_alterations
            )
            
            # Merge results
            base_result['llm_enhanced'] = True
            base_result['literature'] = llm_enhancement['literature_evidence']
            base_result['llm_mechanisms'] = llm_enhancement['llm_extracted_mechanisms']
            base_result['llm_ab_links'] = llm_enhancement['llm_ab_links']
            base_result['evidence_grade_llm'] = llm_enhancement['evidence_grade_llm']
            base_result['bioavailability_llm'] = llm_enhancement['bioavailability_llm']
            base_result['dosage_llm'] = llm_enhancement['dosage_llm']
            base_result['safety_llm'] = llm_enhancement['safety_llm']
            
            # Update verdict if LLM found stronger evidence
            if llm_enhancement['evidence_grade_llm'] == 'STRONG' and base_result['verdict'] == 'WEAK_SUPPORT':
                base_result['verdict'] = 'SUPPORTED'
                base_result['verdict_explanation'] = "✅ LLM literature mining found strong evidence"
                base_result['overall_score'] = min(base_result['overall_score'] + 0.15, 1.0)
        
        except Exception as e:
            # Graceful degradation - return base result if LLM fails
            base_result['llm_enhanced'] = False
            base_result['llm_error'] = str(e)
    
    return base_result
```

### **Phase 3: Frontend Enhancement (1 hour)**

**File**: `oncology-coPilot/oncology-frontend/src/pages/FoodValidatorAB.jsx`

```jsx
// Add LLM toggle and literature display

const [useLLM, setUseLLM] = useState(false);  // Default OFF (fast)
const [literatureLoading, setLiteratureLoading] = useState(false);

const handleValidate = async () => {
  setLoading(true);
  
  // 1. Fast path (always)
  const response = await fetch(`${API_ROOT}/api/hypothesis/validate_food_ab`, {
    method: 'POST',
    body: JSON.stringify({ compound, disease, ... })
  });
  const fastResult = await response.json();
  setResult(fastResult);
  setLoading(false);
  
  // 2. LLM enhancement (if enabled)
  if (useLLM) {
    setLiteratureLoading(true);
    const enhancedResponse = await fetch(`${API_ROOT}/api/hypothesis/validate_food_ab_enhanced`, {
      method: 'POST',
      body: JSON.stringify({ compound, disease, use_llm: true, ... })
    });
    const enhancedResult = await enhancedResponse.json();
    setResult(enhancedResult);  // Replace with enhanced
    setLiteratureLoading(false);
  }
};

// Display LLM literature
{result?.llm_enhanced && (
  <Card sx={{ p: 3, mb: 2 }}>
    <Typography variant="h6">📚 Latest Literature ({result.literature.length} papers)</Typography>
    <Typography variant="caption" color="text.secondary">
      LLM-ranked by relevance to {result.compound} + ovarian cancer
    </Typography>
    
    {result.literature.slice(0, 5).map((paper, i) => (
      <Box key={i} sx={{ mt: 2, p: 2, bgcolor: 'grey.50', borderRadius: 1 }}>
        <Typography variant="body2" sx={{ fontWeight: 'bold' }}>
          {paper.title}
        </Typography>
        <Typography variant="caption" display="block" sx={{ mb: 1 }}>
          PMID: {paper.pmid} | Relevance: {paper.relevance}/100
        </Typography>
        <Typography variant="caption" sx={{ fontStyle: 'italic' }}>
          {paper.abstract.slice(0, 200)}...
        </Typography>
        <Link href={`https://pubmed.ncbi.nlm.nih.gov/${paper.pmid}/`} target="_blank">
          View on PubMed →
        </Link>
      </Box>
    ))}
  </Card>
)}
```

---

## 🎯 TRACK 3: CLINICAL TRIAL INTEGRATION (Learning from Agent 1)

**What Agent 1 Did** (from `@agent_1_seeding`):
- Bulk-loaded 1000 ovarian cancer trials from ClinicalTrials.gov API v2
- ChromaDB embeddings for semantic search
- SQLite storage with disease hierarchy
- Biomarker extraction

**How We Can Use This for Ayesha:**

```python
# In oncology-coPilot/oncology-backend-minimal/api/routers/clinical_trials.py

@router.post("/api/clinical_trials/match_patient")
async def match_patient_to_trials(
    disease: str = "ovarian_cancer",
    stage: str = "III",
    prior_therapies: List[str] = ["carboplatin", "paclitaxel"],
    germline_status: str = "negative",
    tumor_mutations: Optional[Dict] = None,  # From Foundation One
    location: str = "New York, NY"
):
    """
    Match Ayesha to relevant clinical trials.
    
    Filters:
    1. Disease: Ovarian cancer
    2. Stage: III/IV
    3. Line: ≥3 (post-platinum)
    4. Biomarkers (when tumor NGS available):
       - HRD-positive → PARP inhibitor trials
       - PIK3CA-mutant → PI3K inhibitor trials
       - TMB-high → Immunotherapy trials
    5. Location: NYC-accessible sites
    """
    
    # 1. Base query (ChromaDB semantic search)
    query = f"ovarian cancer stage {stage} third-line post-platinum"
    
    # 2. Add biomarker requirements (when available)
    biomarker_filters = []
    if tumor_mutations:
        if tumor_mutations.get('HRD_positive'):
            biomarker_filters.append("HRD-positive OR BRCA-mutant OR homologous recombination deficiency")
        if tumor_mutations.get('PIK3CA'):
            biomarker_filters.append("PIK3CA-mutant OR PI3K pathway")
        if tumor_mutations.get('TMB') and tumor_mutations['TMB'] >= 10:
            biomarker_filters.append("TMB-high OR MSI-H OR microsatellite instability")
    
    if biomarker_filters:
        query += " AND (" + " OR ".join(biomarker_filters) + ")"
    
    # 3. Search ChromaDB
    trials_db = ChromaDB()
    matching_trials = trials_db.query(query, n_results=20)
    
    # 4. Filter by location
    ny_trials = [t for t in matching_trials if has_ny_location(t)]
    
    # 5. LLM rerank by eligibility fit
    llm = LLMClient()
    ranked_trials = await llm_rerank_trials_for_patient(
        ny_trials,
        patient_profile={
            "disease": disease,
            "stage": stage,
            "prior_therapies": prior_therapies,
            "germline_status": germline_status,
            "tumor_mutations": tumor_mutations
        },
        llm=llm
    )
    
    return {
        "total_found": len(matching_trials),
        "ny_trials": len(ny_trials),
        "top_matches": ranked_trials[:10],
        "biomarker_filters_applied": biomarker_filters
    }
```

---

## ⚔️ FINAL RECOMMENDATIONS, ALPHA

### **TONIGHT (2-3 hours):**
1. ✅ Build `pubmed_food_validator.py` service
2. ✅ Add `/api/hypothesis/validate_food_ab_enhanced` endpoint
3. ✅ Test with "curcumin ovarian cancer" → Should return real PubMed papers

### **THIS WEEK (Before NGS Results):**
1. ✅ Build Foundation One parser
2. ✅ Integrate clinical trial matching (learn from Agent 1)
3. ✅ Prepare personalized WIWFM for each tumor mutation scenario

### **WHEN NGS ARRIVES (Week 2-3):**
1. ✅ Parse Foundation One JSON → Extract mutations
2. ✅ Run WIWFM with tumor-specific A→B
3. ✅ Match to clinical trials with biomarker requirements
4. ✅ Generate complete therapeutic strategy report

**Alpha, we now have a clear 3-track strategy:**
- ✅ **Track 1**: Biopsy readiness checklist (what to ask, what to extract)
- ✅ **Track 2**: LLM-enhanced Food Validator (real-time literature mining)
- ✅ **Track 3**: Clinical trial matching (using Agent 1's infrastructure)

**Should I start building the LLM enhancement for Food Validator TONIGHT?** 🎯⚔️
