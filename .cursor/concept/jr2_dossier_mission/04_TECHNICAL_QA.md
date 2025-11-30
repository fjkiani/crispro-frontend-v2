# ⚔️ JR2 TECHNICAL Q&A - ALL QUESTIONS ANSWERED ⚔️

**Date**: January 13, 2025  
**From**: Zo (Lead Commander)  
**To**: JR2 (Dossier Sidekick)  
**Status**: ✅ **ALL QUESTIONS ANSWERED - NO BLOCKERS**

---

## 📋 **DATA SOURCES & FORMATS**

### **Q1: Patient Profile Structure - None vs "UNKNOWN"**
✅ **ANSWER**: Use "UNKNOWN" string (not None)

```python
patient_profile = {
    "her2_status": "UNKNOWN",  # ✅ Use this (string)
    "hrd_score": None,          # ❌ Don't use this for biomarkers
    "brca_status": "WILDTYPE",  # ✅ Known value
}
```

**Rationale**: Frontend displays "UNKNOWN" as "⚠️ PENDING - Order test". None would cause display errors.

---

### **Q2: Trial Data Sources - Exact AstraDB Schema**
✅ **ANSWER**: Full AstraDB document structure provided

See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for complete schema.

**Key Fields**:
- `nct_id`, `title`, `status`, `phase`
- `eligibility_text` (truncated to 7500 bytes - need to scrape full page)
- `locations_data` (array of location objects)
- `$vector` (768-dim embedding)

---

### **Q3: Trial Scraping - Tool & Rate Limits**
✅ **ANSWER**: Use Diffbot (already integrated, better extraction)

**Why Diffbot**: Already set up in codebase, better HTML parsing, handles JavaScript-rendered content.

**Existing Integration**:
- Endpoint: `POST /api/evidence/extract` (at `api/routers/evidence/extraction.py`)
- Config: `DIFFBOT_TOKEN` from environment
- Returns: `{title, author, date, site_name, text, html, tags}`

**Code**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for scraping function using Diffbot.

**Rate Limits**: Handled by Diffbot service (no manual rate limiting needed).

---

## 🔧 **IMPLEMENTATION DETAILS**

### **Q4: Drug Mechanism Database - Exact List**
✅ **ANSWER**: 20 drugs pre-populated

See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for complete `DRUG_MECHANISM_DB` with all 20 drugs.

**Structure**: `drug_name → {class, target, moa, layman, breakthrough}`

**If drug not in database**: Flag for Zo review with `[UNKNOWN DRUG - NEEDS RESEARCH]`

---

### **Q5: Eligibility Matching Logic**
✅ **ANSWER**: Complete matching functions provided in `CLIENT_DOSSIER_DOCTRINE.mdc`

**Functions**:
- `assess_disease_match()` - Disease/stage matching
- `assess_treatment_line_match()` - Treatment line matching
- `assess_biomarker_match()` - Biomarker gates with exact logic

**Biomarker Gate Logic**:
- `UNKNOWN` → `⚠️ PENDING` status with action
- Match → `✅ PASS` status
- Mismatch → `❌ FAIL` status with explanation

---

### **Q6: Evidence Synthesis - Integration Strategy**
✅ **ANSWER**: Use existing `EnhancedEvidenceService` + fallback to `DRUG_MECHANISM_DB`

**Code**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for integration example.

**Fallback Strategy**:
1. Query `EnhancedEvidenceService`
2. If no evidence → Use `DRUG_MECHANISM_DB`
3. If not in DB → Flag `[NEEDS EVIDENCE RESEARCH]`

---

## 📊 **OUTPUT FORMATS & QUALITY**

### **Q7: File Naming Convention**
✅ **ANSWER**: `dossier_{nct_id}_{patient_id}_{timestamp}.md`

**Examples**:
- `dossier_NCT06819007_ayesha_001_20250113T220000.md`
- `dossier_NCT03705156_ayesha_001_20250113T230000.md`

**Storage**:
- Generated: `.cursor/ayesha/dossiers/{nct_id}/`
- Approved: `.cursor/ayesha/dossiers/approved/{nct_id}/`
- Rejected: `.cursor/ayesha/dossiers/rejected/{nct_id}/`

---

### **Q8: How to Flag Uncertain Claims**
✅ **ANSWER**: Use confidence flags in markdown

**Flags**:
- `[INFERRED]` - Medium confidence (inferred from data)
- `[NEEDS VERIFICATION]` - Low confidence (needs Zo review)
- `[MISSING]` - Missing data (critical gap)
- `[CONFLICT]` - Conflicting data (needs resolution)

**Examples**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md)

---

### **Q9: Error Handling - Drug Mechanism Fallback**
✅ **ANSWER**: 3-tier fallback strategy

1. **Tier 1**: Check `DRUG_MECHANISM_DB`
2. **Tier 2**: Check `EnhancedEvidenceService`
3. **Tier 3**: Flag for Zo review

**Code**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md)

---

## 🔌 **INTEGRATION POINTS**

### **Q10: Router Location & Database Storage**
✅ **ANSWER**: Create new router at `api/routers/dossiers.py`

**Endpoints**:
- `POST /api/dossiers/generate` - Generate dossier
- `GET /api/dossiers/{id}` - Get dossier
- `POST /api/dossiers/{id}/approve` - Zo review/approve

**Database**: Store in AstraDB collection `clinical_dossiers`

**Code**: See [09_API_SPECIFICATIONS.md](./09_API_SPECIFICATIONS.md)

---

### **Q11: Caching Strategy**
✅ **ANSWER**: File-based caching (24hr TTL)

**Location**: `.cursor/ayesha/cache/trial_{nct_id}.json`

**TTL**: 24 hours for scraped trial data

**Code**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md)

---

### **Q12: Dependencies**
✅ **ANSWER**: Reuse existing utilities

**Already Exists**:
- Logging: `api.utils.logger`
- Database: `api.services.database_connections`
- Evidence: `api.services.enhanced_evidence_service`

**Add to requirements.txt**:
- `beautifulsoup4==4.12.2`
- `lxml==4.9.3`

---

## 🎯 **PRIORITY & SCOPE**

### **Q13: MVP vs Full Implementation**
✅ **ANSWER**: All 10 sections required

**Rollout Plan**:
- Phase 1: JR1 (Trial Seeding)
- Phase 2: JR2 (Dossier Generation - all 10 sections)
- Phase 3: ZO (Manual Review & Approval)

**Build incrementally**: One section at a time, test, then next section.

---

### **Q14: Sample Data**
✅ **ANSWER**: Use NCT06819007 + Ayesha profile

**Sample Trial**: NCT06819007 (already verified by Zo)

**Edge Cases Covered**:
- Missing biomarkers → `UNKNOWN` → `⚠️ PENDING`
- Incomplete trial data → Merge `astradb_record` + `scraped_full`

---

### **Q15: Timeline & Iteration**
✅ **ANSWER**: Submit complete dossier (all 10 sections) for Zo review

**Process**:
1. JR2 generates complete dossier
2. Submit to Zo for review
3. Zo reviews → Returns `edits_required` list
4. JR2 iterates based on feedback

**Expected**: 1-2 rounds of feedback typical

---

## ⚔️ **FINAL CHECKLIST**

**All Questions Answered** ✅:
- ✅ Q1: Use "UNKNOWN" string (not None)
- ✅ Q2: AstraDB schema provided
- ✅ Q3: Use Diffbot (already integrated, better extraction)
- ✅ Q4: 20 drugs pre-populated in DRUG_MECHANISM_DB
- ✅ Q5: Eligibility matching logic provided
- ✅ Q6: Use EnhancedEvidenceService + fallback
- ✅ Q7: File naming convention specified
- ✅ Q8: Confidence flags defined
- ✅ Q9: 3-tier fallback for drug mechanisms
- ✅ Q10: Router location specified
- ✅ Q11: File-based caching (24hr TTL)
- ✅ Q12: Reuse existing utilities; Diffbot already integrated
- ✅ Q13: All 10 sections required
- ✅ Q14: Sample data provided
- ✅ Q15: Review process defined

**NO MORE BLOCKERS - JR2 CAN START BUILDING!** 🔥⚔️

---

**Next Steps**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for code examples

