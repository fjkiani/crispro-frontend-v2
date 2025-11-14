# ⚔️ AYESHA DEMO PLAYBOOK - COMPLETE CONSOLIDATION ⚔️

**Date**: January 13, 2025  
**Status**: ✅ **CONSOLIDATED SOURCE OF TRUTH**  
**Consolidated From**:
- `AYESHA_DEMO_SCRIPT.md` - Clinical trial explorer demo script
- `AYESHA_DEMO_WORKFLOW_COMPLETE.md` - Complete 8-step workflow
- `DEMO_EXECUTION_MASTER_PLAN.md` - Complete execution playbook
- `AYESHA_DEMO_READY_STATUS.md` - Demo readiness status
- `V2_DEMO_COMPREHENSIVE_PLAN.md` - V2 multi-modal demo plan
- `V2_DEMO_FINAL_EXECUTION_PLAN.md` - V2 final execution plan
- `METASTASIS_DEMO_V2_PLAN.md` - Metastasis demo plan

---

## 📚 TABLE OF CONTENTS

1. [Executive Summary](#executive-summary)
2. [Quick Start Guide](#quick-start-guide)
3. [Demo Scripts](#demo-scripts)
4. [Validation & Testing](#validation--testing)
5. [Troubleshooting](#troubleshooting)
6. [References to Archived Files](#references-to-archived-files)

---

## 🎯 EXECUTIVE SUMMARY

This master document consolidates all demo scripts, workflows, and execution plans for Ayesha's clinical care demonstrations.

**Key Demos**:
- **Clinical Trial Explorer**: 5-7 minute demo showing trial matching + SOC recommendations
- **Sporadic Cancer Workflow**: 8-10 minute demo showing tumor-centric analysis
- **V2 Multi-Modal Demo**: Complete Oracle→Forge→Gauntlet→Dossier workflow

**Demo Status**: ✅ **100% READY** - All scripts, workflows, and test data complete

---

## 🚀 QUICK START GUIDE

### **Setup (2 Minutes)**

**Step 1: Start Backend** (30 seconds)
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

**Step 2: Start Frontend** (30 seconds)
```bash
cd oncology-coPilot/oncology-frontend
npm run dev
```

**Step 3: Verify Health** (1 minute)
```bash
curl http://127.0.0.1:8000/healthz
curl http://127.0.0.1:8000/api/ayesha/trials/health
```

**Expected Output**: `{"status": "ok"}`

---

## 🎬 DEMO SCRIPTS

### **Demo 1: Clinical Trial Explorer** (5-7 minutes)

**Purpose**: Show oncologist complete Ayesha care plan system

**Key Messages**:
1. **Transparency**: Every recommendation has clear reasoning
2. **Speed**: Complete care plan in 3-5 minutes (vs 2-4 weeks manual)
3. **Confidence**: Deterministic confidence gates (not black-box AI)
4. **Early Detection**: CA-125 monitoring catches resistance 3-6 weeks early

**Flow**:
1. **Profile Summary** (30 sec) - Stage IVB, CA-125 2842, germline-negative
2. **SOC Recommendation** (1 min) - Carboplatin + Paclitaxel + Bevacizumab (95% confidence)
3. **CA-125 Intelligence** (1 min) - Forecast, resistance signals, monitoring strategy
4. **Top 10 Trials** (2 min) - Ranked with transparent reasoning
5. **Eligibility Checklists** (30 sec) - Hard/soft criteria split
6. **NGS Fast-Track** (30 sec) - Test ordering checklist
7. **Provenance** (30 sec) - Run ID, profile, confidence gates

**Full Script**: See `AYESHA_DEMO_SCRIPT.md` (290 lines)

---

### **Demo 2: Sporadic Cancer Workflow** (8-10 minutes)

**Purpose**: Show tumor-centric analysis for 85-90% of cancer patients

**Key Messages**:
1. **Majority Focus**: 85-90% of cancers are sporadic (non-hereditary)
2. **Tumor Genomics**: Analysis based on tumor mutations, not just germline
3. **Confidence Caps**: Transparent confidence limits (L0: 40%, L1: 60%, L2: none)
4. **PARP Rescue**: HRD ≥42 rescues PARP penalty (somatic HRD)

**Flow**:
1. **Germline Status** (30 sec) - Negative banner, explain majority
2. **Quick Intake** (1 min) - Fill form, generate Level 0 estimates
3. **Efficacy Prediction** (2 min) - Show PARP penalty + IO boost
4. **Provenance Card** (1 min) - Show gate explanations
5. **NGS Upload** (2 min) - Upload Level 2 report, show rescue
6. **Re-Analysis** (1 min) - Show improved confidence (L2)
7. **WIWFM Integration** (1 min) - Show drug ranking with sporadic context
8. **Export** (30 sec) - Copy to clipboard, provider report

**Full Script**: See `AYESHA_DEMO_WORKFLOW_COMPLETE.md` (840 lines)

---

### **Demo 3: V2 Multi-Modal AI Showcase** (10-15 minutes)

**Purpose**: Show complete Oracle→Forge→Gauntlet→Dossier workflow

**Key Messages**:
1. **Multi-Modal AI**: Sequence + Structure + Evidence integration
2. **Structural Validation**: Prevents "wet noodle" proteins (pLDDT ≥70)
3. **FDA-Ready**: Complete IND package generation
4. **Explainable**: SAE features show mechanistic interpretation

**Flow**:
1. **Oracle Analysis** (2 min) - Variant scoring with Evo2
2. **Forge Generation** (2 min) - Therapeutic protein design
3. **Gauntlet Validation** (3 min) - Structural validation (pLDDT scoring)
4. **Dossier Generation** (3 min) - IND package creation
5. **SAE Interpretation** (2 min) - Explainable features

**Full Script**: See `V2_DEMO_COMPREHENSIVE_PLAN.md` (528 lines)

---

## 🧪 VALIDATION & TESTING

### **Automated Validation Suite**

**File**: `test_data/DEMO_VALIDATION_SUITE.py`

**Tests** (6 total):
1. ✅ Backend health check
2. ✅ Quick Intake (Level 0)
3. ✅ Efficacy prediction (Level 0)
4. ✅ NGS upload (Level 2)
5. ✅ Efficacy prediction (Level 2 - PARP rescue)
6. ✅ IO boost validation

**Run Command**:
```bash
cd .cursor/ayesha/test_data
python DEMO_VALIDATION_SUITE.py
```

**Expected Output**: 6/6 tests passing (green)

---

### **Test Data Files**

1. **`test_data/ayesha_level0_intake.json`** - Level 0 Quick Intake
2. **`test_data/ayesha_tumor_ngs.json`** - Level 2 NGS Report

**Usage**: Load into Quick Intake form or NGS upload component

---

## ⚠️ TROUBLESHOOTING

### **Backend Won't Start**
- Check Python version: `python3 --version` (need 3.8+)
- Check dependencies: `pip install -r requirements.txt`
- Check port 8000 not in use: `lsof -i :8000`

### **Frontend Won't Start**
- Check Node version: `node --version` (need 16+)
- Check dependencies: `npm install`
- Check port 5173 not in use: `lsof -i :5173`

### **API Call Fails**
- Check backend is running: `curl http://localhost:8000/health`
- Check CORS: Backend should allow `http://localhost:5173`
- Check network tab: Look for 500 errors, check response body

### **Trials List Empty**
- **Expected**: If AstraDB not seeded or no matching trials
- **Not a Bug**: System is working correctly, just no matches
- **Solution**: Seed more frontline ovarian trials or adjust filters

---

## 📖 REFERENCES TO ARCHIVED FILES

All original demo files have been preserved in `archive/demo_iterations/`:

- **Demo Script**: `archive/demo_iterations/AYESHA_DEMO_SCRIPT.md`
- **Demo Workflow**: `archive/demo_iterations/AYESHA_DEMO_WORKFLOW_COMPLETE.md`
- **Execution Plan**: `archive/demo_iterations/DEMO_EXECUTION_MASTER_PLAN.md`
- **Demo Status**: `archive/demo_iterations/AYESHA_DEMO_READY_STATUS.md`
- **V2 Plans**: `archive/demo_iterations/V2_DEMO_*.md`
- **Metastasis Demo**: `archive/demo_iterations/METASTASIS_DEMO_V2_PLAN.md`

---

## ⚔️ DOCTRINE STATUS: ACTIVE

**LAST UPDATED:** January 13, 2025  
**APPLIES TO:** All Ayesha demo scripts, workflows, and execution plans  
**ENFORCEMENT:** Mandatory for all demo preparations and presentations

**This master document represents the complete consolidation of all demo playbooks. Every script, workflow, and execution plan is preserved and organized for maximum clarity and actionability.**

