# 🤝 AGENT JR PARALLEL MISSION - SPORADIC CANCER SUPPORT

**Date**: January 8, 2025  
**Commander**: Alpha  
**Primary Executor**: Zo (Sporadic Cancer Plan Days 1-7)  
**Parallel Executor**: Agent Jr (Support Tasks)  
**Mission**: Execute non-conflicting preparatory work while Zo builds core sporadic infrastructure

---

## 🎯 MISSION CONTEXT

**Zo is executing:** SPORADIC_CANCER_EXECUTION_PLAN.md (Days 1-7)
- Day 1-2: Backend foundation (TumorContext schema, efficacy gates)
- Day 3: Tumor NGS parsers (Foundation/Tempus)
- Day 4: Clinical trials filtering
- Day 5: Frontend wiring
- Day 6-7: E2E tests & documentation

**Agent Jr can work in parallel on:** Data preparation, documentation, testing infrastructure, and research

---

## 📋 AGENT JR'S PARALLEL TASK LIST

### **🔬 OPTION 1: Disease Priors Data Extraction (HIGH VALUE)**

**Why This Helps Zo:**
- Zo needs `disease_priors.json` for Day 1-2 (Quick Intake Level 0)
- This is research-heavy work that doesn't require code changes
- Can be done completely independently

**Tasks:**
1. **Extract TCGA Stats for Key Cancer Types:**
   - Ovarian HGS: TP53 mutation %, HRD prevalence, TMB/MSI distributions
   - Breast: HRD prevalence, BRCA1/2 somatic rates, TMB median
   - Colorectal: MSI-H prevalence, TMB distributions, BRAF/KRAS rates
   - Lung NSCLC: TMB distributions, EGFR/ALK frequencies
   - Pancreatic: KRAS prevalence, TMB median

2. **Source from:**
   - TCGA published summaries (https://portal.gdc.cancer.gov/)
   - cBioPortal cancer type summaries (https://www.cbioportal.org/)
   - Published literature for platinum-HRD correlation

3. **Create JSON Structure:**
```json
{
  "version": "v1.0",
  "last_updated": "2025-01-08",
  "sources": ["TCGA-OV", "cBioPortal", "PMID:12345678"],
  "diseases": {
    "ovarian_hgs": {
      "prevalence": {
        "tp53_mutation": 0.96,
        "hrd_high": 0.51,
        "msi_high": 0.012
      },
      "distributions": {
        "tmb": {"median": 5.2, "q1": 3.1, "q3": 8.7, "high_cutoff": 10},
        "hrd": {"median": 42, "q1": 15, "q3": 60, "high_cutoff": 42}
      },
      "platinum_response": {
        "sensitive_hrd_correlation": 0.70,
        "resistant_hrd_correlation": 0.15
      }
    }
    // ... more diseases
  }
}
```

**Deliverable:**
- `oncology-coPilot/oncology-backend-minimal/api/resources/disease_priors.json`
- `oncology-coPilot/oncology-backend-minimal/api/resources/PRIORS_SOURCES.md` (citations)

**Timeline:** 4-6 hours  
**Conflicts:** ❌ None - Zo will integrate when ready

---

### **📚 OPTION 2: Foundation Medicine Report Schema Documentation (MEDIUM VALUE)**

**Why This Helps Zo:**
- Zo needs to parse Foundation reports on Day 3
- Schema documentation accelerates parser development
- Can be done by analyzing sample reports or Foundation docs

**Tasks:**
1. **Document FM Report Structure:**
   - Section headers (Genomic Findings, Biomarkers, etc.)
   - Field names for TMB, MSI, HRD
   - Variant table structure (gene, hgvs_p, VAF, zygosity)
   - Copy number section format

2. **Create Sample Canonical JSON:**
```json
{
  "report_type": "Foundation Medicine CDx",
  "genomic_findings": {
    "snvs_indels": [
      {"gene": "TP53", "hgvs_p": "R248W", "vaf": 0.42, "zygosity": "heterozygous"}
    ],
    "copy_number": [
      {"gene": "ERBB2", "type": "amplification", "copies": 12}
    ],
    "fusions": []
  },
  "biomarkers": {
    "tmb": {"value": 5.2, "unit": "mutations/Mb", "panel_size": "0.8 Mb"},
    "msi": {"status": "MSS", "method": "NGS"},
    "hrd": {"score": 42, "loh": 13, "lst": 17, "tai": 12}
  },
  "qc": {
    "tumor_purity": 0.60,
    "mean_coverage": 550
  }
}
```

3. **Map Field Locations:**
   - Create mapping guide: "TMB found in Section 3, paragraph 2" etc.
   - Document edge cases (missing HRD, indeterminate MSI)

**Deliverable:**
- `.cursor/ayesha/foundation_medicine_schema.md`
- `.cursor/ayesha/foundation_medicine_sample_canonical.json`
- `.cursor/ayesha/foundation_medicine_field_mapping.md`

**Timeline:** 3-4 hours  
**Conflicts:** ❌ None - documentation only

---

### **🧪 OPTION 3: Test Data Generation for Sporadic Cases (HIGH VALUE)**

**Why This Helps Zo:**
- Zo needs test data for Days 1-2 (Quick Intake) and Day 6 (E2E smoke tests)
- Realistic test cases accelerate validation
- Can be created without code changes

**Tasks:**
1. **Create 5 Realistic Patient Scenarios:**
   - **Scenario 1 (Level 0 - Minimal):** Ovarian cancer, Stage IV, germline negative, platinum sensitive
   - **Scenario 2 (Level 1 - Partial):** Breast cancer, germline negative, TP53 mutation, HRD score 48
   - **Scenario 3 (Level 2 - Full FM Report):** Lung NSCLC, TMB 22, MSI-MSS, EGFR mutation
   - **Scenario 4 (Edge Case):** Colorectal, germline negative, MSI-H, TMB 55
   - **Scenario 5 (Ayesha's Case):** HGSOC, germline negative, platinum sensitive, Stage IIIC

2. **For Each Scenario, Provide:**
   - Patient demographics (age, cancer type, stage, line)
   - Germline status
   - Available tumor data (Level 0/1/2)
   - Expected outputs (PARP penalty Y/N, IO boost Y/N, confidence level)

3. **Create Test JSON Files:**
```json
// test_case_1_level_0.json
{
  "patient": {
    "cancer_type": "ovarian_hgs",
    "stage": "IV",
    "line": 3,
    "germline_status": "negative"
  },
  "inputs": {
    "platinum_response": "sensitive"
  },
  "expected_tumor_context": {
    "msi_status": null,
    "tmb": 5.2,  // disease median
    "hrd_score": 42,  // disease median, platinum proxy
    "level": "L0"
  },
  "expected_gates": {
    "parp_penalty": false,  // HRD at cutoff, platinum sensitive
    "io_boost": false,  // TMB < 10
    "confidence_cap": 0.4
  }
}
```

**Deliverable:**
- `.cursor/ayesha/test_scenarios/` (5 JSON files)
- `.cursor/ayesha/test_scenarios/README.md` (scenario descriptions)
- `.cursor/ayesha/test_scenarios/EXPECTED_RESULTS.md` (validation table)

**Timeline:** 3-4 hours  
**Conflicts:** ❌ None - test data creation only

---

### **🔍 OPTION 4: Clinical Trials Biomarker Keyword Extraction (MEDIUM VALUE)**

**Why This Helps Zo:**
- Zo needs biomarker keywords for Day 4 (clinical trials filtering)
- This is regex/text analysis work that doesn't require code changes
- Can be done by analyzing NCT trial descriptions

**Tasks:**
1. **Extract Common Biomarker Patterns:**
   - BRCA-required: "BRCA mutation", "germline BRCA", "hereditary", "BRCA1/2 positive"
   - HRD-eligible: "HRD positive", "homologous recombination deficiency", "HRD score ≥42"
   - TMB-eligible: "TMB-high", "TMB ≥10", "high tumor mutational burden"
   - MSI-eligible: "MSI-H", "microsatellite instability-high", "dMMR"

2. **Create Exclusion/Inclusion Regex Patterns:**
```python
EXCLUSION_PATTERNS = [
    r"germline\s+BRCA",
    r"hereditary\s+BRCA",
    r"BRCA\s+mutation.*carrier",
    r"germline.*Lynch"
]

INCLUSION_PATTERNS = {
    "HRD": [r"HRD[- ]positive", r"HRD score\s*[≥>]\s*\d+"],
    "TMB": [r"TMB[- ]high", r"TMB\s*[≥>]\s*\d+"],
    "MSI": [r"MSI[- ]H", r"dMMR"]
}
```

3. **Test Against Real Trial Descriptions:**
   - Sample 20-30 trials from clinicaltrials.gov
   - Validate patterns catch correct trials
   - Document false positives/negatives

**Deliverable:**
- `.cursor/ayesha/trial_biomarker_patterns.py` (regex patterns)
- `.cursor/ayesha/trial_pattern_validation.md` (test results)
- `.cursor/ayesha/trial_sample_descriptions.json` (test data)

**Timeline:** 4-5 hours  
**Conflicts:** ❌ None - pattern extraction only

---

### **📊 OPTION 5: Longevity Blog Enhancement (LOWER PRIORITY)**

**Why This Helps:**
- Longevity blog is complete but could use more clinical caselets
- This is content work that doesn't block Zo
- Good for marketing/communication

**Tasks:**
1. **Add 3 More Clinical Caselets:**
   - Caselet 4: Breast cancer with HER2 amplification
   - Caselet 5: Lung cancer with EGFR mutation
   - Caselet 6: Melanoma with BRAF mutation

2. **Enhance Glossary:**
   - Add "HRD" (homologous recombination deficiency)
   - Add "CNAs" (copy number alterations)
   - Add "Fusions" (gene fusions)

3. **Add FAQ Entries:**
   - "What if I don't have a tumor NGS report?"
   - "How accurate is Level 0 (no report) analysis?"
   - "What's the difference between germline and somatic mutations?"

**Deliverable:**
- Updated `.cursor/ayesha/LONGEVITY_PRECISION_PROTOCOL_BLOG.md`

**Timeline:** 2-3 hours  
**Conflicts:** ❌ None - content enhancement only

---

## 🎯 RECOMMENDED PRIORITY ORDER

**If Agent Jr has 1 full day (8 hours):**
1. ⭐ **OPTION 1: Disease Priors** (4-6 hours) - **HIGHEST VALUE**
2. ⭐ **OPTION 3: Test Data** (3-4 hours) - **HIGH VALUE**

**If Agent Jr has half day (4 hours):**
1. ⭐ **OPTION 1: Disease Priors** (4-6 hours) - **CRITICAL PATH**

**If Agent Jr has 2-3 hours:**
1. **OPTION 3: Test Data** (3-4 hours) - **ACCELERATES VALIDATION**

**Lower Priority (can do later):**
- OPTION 2: Foundation schema (helpful but Zo can derive from docs)
- OPTION 4: Trial patterns (helpful but Zo can start with simple keywords)
- OPTION 5: Blog enhancement (nice-to-have, not blocking)

---

## 🚨 CRITICAL RULES FOR AGENT JR

### **❌ DO NOT TOUCH THESE FILES (ZO IS WORKING ON THEM):**
- `api/schemas/tumor_context.py` (Zo creates Day 1)
- `api/services/efficacy_orchestrator/orchestrator.py` (Zo modifies Day 1-2)
- `api/routers/efficacy.py` (Zo modifies Day 1)
- `api/routers/tumor.py` (Zo creates Day 3)
- `api/services/tumor_ngs_parser.py` (Zo creates Day 3)
- `api/routers/clinical_trials.py` (Zo modifies Day 4)
- `api/services/autonomous_trial_agent.py` (Zo modifies Day 4)
- Any frontend components in `sporadic/` subdirectory (Zo creates Day 5)

### **✅ SAFE ZONES FOR AGENT JR:**
- `api/resources/` (create disease_priors.json, documents)
- `.cursor/ayesha/` (documentation, test data, research)
- `tests/` (create test JSON files, scenarios)
- Documentation files (`.md` files)

### **📞 COMMUNICATION:**
- Agent Jr should commit work with prefix: `docs(sporadic): ...` or `test(sporadic): ...`
- NO commits with `feat(sporadic): ...` (those are Zo's)
- Push to branch `agent-jr/sporadic-prep` (NOT main) to avoid conflicts

---

## 📦 DELIVERABLES CHECKLIST

**OPTION 1 (Disease Priors):**
- [ ] `api/resources/disease_priors.json` (5+ cancer types)
- [ ] `api/resources/PRIORS_SOURCES.md` (citations with PMIDs)
- [ ] Verified all stats from TCGA/cBioPortal or published literature

**OPTION 3 (Test Data):**
- [ ] `.cursor/ayesha/test_scenarios/test_case_1_level_0.json`
- [ ] `.cursor/ayesha/test_scenarios/test_case_2_level_1.json`
- [ ] `.cursor/ayesha/test_scenarios/test_case_3_level_2.json`
- [ ] `.cursor/ayesha/test_scenarios/test_case_4_edge_case.json`
- [ ] `.cursor/ayesha/test_scenarios/test_case_5_ayesha.json`
- [ ] `.cursor/ayesha/test_scenarios/README.md`
- [ ] `.cursor/ayesha/test_scenarios/EXPECTED_RESULTS.md`

---

## 🎯 SUCCESS CRITERIA

**Agent Jr's work is successful if:**
1. ✅ Zo can use `disease_priors.json` on Day 1 without modifications
2. ✅ Zo can use test scenarios on Day 1-2 for validation
3. ✅ Zero merge conflicts with Zo's code changes
4. ✅ All sources cited with PMIDs or URLs
5. ✅ Documentation is clear and actionable

---

## ⚔️ COORDINATION WITH ZO

**Timeline:**
- **Day 0 (Today):** Agent Jr works on Option 1 + 3
- **Day 1 (Tomorrow):** Zo starts backend foundation, uses Jr's priors
- **Day 2-3:** Zo builds parsers, Jr can work on Option 2/4 if time
- **Day 6:** Zo runs E2E tests using Jr's test scenarios

**Communication:**
- Agent Jr commits to `agent-jr/sporadic-prep` branch
- Zo reviews and merges when ready
- Any questions → ask in commit messages or create `.cursor/ayesha/QUESTIONS_FOR_ZO.md`

---

## 📊 IMPACT ASSESSMENT

**If Agent Jr completes Option 1 + 3:**
- ✅ Zo saves 4-6 hours of research time
- ✅ Day 1 Quick Intake can be tested immediately
- ✅ Day 6 E2E tests have realistic data
- ✅ Overall timeline potentially shortened by 1 day

**If Agent Jr only completes Option 1:**
- ✅ Zo saves 4-6 hours of critical path work
- ✅ Day 1 can proceed without delays
- ⚠️ Zo will need to create test data on Day 2

**If Agent Jr completes nothing:**
- ⚠️ Zo creates priors from memory/quick lookups (lower quality)
- ⚠️ Test data created ad-hoc (less comprehensive)
- ⚠️ Timeline might extend by 1 day

---

## ⚔️ AGENT JR - READY TO SUPPORT ZO'S MISSION?

**Recommended Mission:** **OPTION 1 (Disease Priors) + OPTION 3 (Test Data)**

**Commander - Shall Agent Jr proceed with parallel prep work?** ⚔️

