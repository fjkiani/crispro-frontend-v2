# Clinical Trials Integration Points

**Date:** January 28, 2025  
**Status:** ✅ **ACTIVE** - Integration points documented  
**Location:** `.cursor/MOAT/CLINICAL_TRIALS/06_INTEGRATION_POINTS.md`  
**See Also:** [00_MISSION.mdc](00_MISSION.mdc) for mission overview, [01_SYSTEM_ARCHITECTURE.md](01_SYSTEM_ARCHITECTURE.md) for architecture

---

## 🔗 Integration with Other MOAT Modules

### **Module 01: Data Extraction**
- **Input:** Patient mutations, disease, biomarkers
- **Output:** Structured mutation data for trial matching
- **Integration:** Trial matching uses extracted mutations for query generation

### **Module 02: Biomarker Intelligence**
- **Input:** TMB, MSI status, HRD score
- **Output:** IO eligibility, biomarker boosts
- **Integration:** Trial matching uses biomarkers for filtering and boosting

### **Module 04: Drug Efficacy (S/P/E Framework)**
- **Input:** Pathway disruption scores, top-ranked drugs
- **Output:** Mechanism vector (7D), auto-inferred interventions
- **Integration:** Mechanism vector used for mechanism fit ranking

### **Module 05: Trial Matching Agent**
- **Input:** Patient profile, mechanism vector
- **Output:** Matched trials with mechanism fit scores
- **Integration:** Core trial matching capability

### **Module 07: Care Plan Agent**
- **Input:** Trial matches from Module 05
- **Output:** Care plan with trial recommendations
- **Integration:** Trial matches included in care plan

### **Module 08: Monitoring Agent**
- **Input:** Trial enrollment status, patient updates
- **Output:** Trial status updates, new trial alerts
- **Integration:** Monitors trial enrollment and updates

---

## 📋 API Endpoints

### **Basic Search**
- `POST /api/trials/search` - Semantic search via AstraDB
- `POST /api/trials/search-optimized` - Graph-optimized search
- `POST /api/trials/agent/search` - Autonomous agent search

### **Advanced Queries**
- `POST /api/trials/advanced-query` - Multi-criteria queries with mechanism fit

### **Integration Endpoints**
- `POST /api/complete_care/v2` - Complete care plan (includes trial matching)
- `POST /api/ayesha/complete_care_v2` - Ayesha-specific care plan

---

## 🔗 Related Files

**System Architecture:**
- [01_SYSTEM_ARCHITECTURE.md](01_SYSTEM_ARCHITECTURE.md) - Current architecture

**Trial Matching Agent:**
- [05_TRIAL_MATCHING_AGENT.mdc](05_TRIAL_MATCHING_AGENT.mdc) - Agent specification

**Core Deliverables:**
- `.cursor/MOAT/CORE_DELIVERABLES/` - Mechanism-Based Trial Matching (core contribution)

---

*Document Owner: Zo*  
*Last Updated: January 28, 2025*  
*Status: ✅ ACTIVE - Integration points documented*



