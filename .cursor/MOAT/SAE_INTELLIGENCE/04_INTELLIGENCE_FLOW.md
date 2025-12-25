# Strategic Intelligence Flow: Test → Signals → Patterns → Actions

**Date:** January 28, 2025  
**Status:** ✅ **ACTIVE** - Intelligence flow documented  
**Location:** `.cursor/MOAT/SAE_INTELLIGENCE/04_INTELLIGENCE_FLOW.md`  
**See Also:** [00_MISSION.mdc](00_MISSION.mdc) for mission overview, [03_GENERALS_BATTLE_MAP.mdc](03_GENERALS_BATTLE_MAP.mdc) for 6 Pillars framework

---

## 🎯 THE STRATEGIC INTELLIGENCE FLOW

### **How Tests → Signals → Patterns → Actions**

```
1. TEST UPLOADED
   ↓
2. BIOMARKER INTELLIGENCE (Pillar 1: Tumor Burden)
   - CA-125: 2,842 → 1,200 → 800 → 1,500 (rising after nadir)
   - Pattern: Resistance emerging
   ↓
3. TRIGGER SYSTEM DETECTS
   - Condition: CA-125 > baseline * 1.25
   - Severity: HIGH
   ↓
4. AUTOMATED ACTIONS
   - Order ctDNA sequencing (Next Test Recommender)
   - Run resistance analysis (Resistance Prediction)
   - Re-rank drugs (Drug Efficacy)
   - Search salvage trials (Trial Matching)
   - Notify oncologist
   ↓
5. CTDNA RESULTS (Pillar 2: Genomic Evolution)
   - New TP53 mutation detected
   - DDR pathway burden: 0.88 → 0.65 (dropping)
   - Pattern: Clonal evolution, resistance mechanism
   ↓
6. UPDATED ACTIONS
   - Drug Efficacy: Re-rank (avoid current drug, prioritize alternatives)
   - Trial Matching: Search PARP+ATR salvage trials
   - Care Plan: Update regimen
   - Monitoring: Increase frequency
```

---

## 📊 EXAMPLE FLOWS

### **Example 1: CA-125 Rising After Nadir**

```
1. TEST: CA-125 = 1,500 (was 800 at nadir)
   ↓
2. PILLAR: Tumor Burden (Pillar 1)
   ↓
3. SIGNAL: Rising after nadir (2x consecutive rise)
   ↓
4. PATTERN: Resistance emerging (doubling time <90 days)
   ↓
5. TRIGGER: resistance_detected (HIGH severity)
   ↓
6. ACTIONS:
   - Order ctDNA sequencing (Next Test Recommender)
   - Run resistance analysis (Resistance Prediction)
   - Re-rank drugs (Drug Efficacy - avoid current, prioritize alternatives)
   - Search salvage trials (Trial Matching - mechanism-aligned)
   - Notify oncologist
   - Escalate to tumor board if no response in 24h
```

---

### **Example 2: New TP53 Mutation in ctDNA**

```
1. TEST: ctDNA shows new TP53 mutation (not in baseline)
   ↓
2. PILLAR: Genomic Evolution (Pillar 2)
   ↓
3. SIGNAL: New mutation detected (clonal evolution)
   ↓
4. PATTERN: Resistance mechanism emerging (TP53 → checkpoint bypass)
   ↓
5. TRIGGER: clonal_evolution (HIGH severity)
   ↓
6. ACTIONS:
   - Update resistance prediction (Resistance Prediction)
   - Re-rank drugs (Drug Efficacy - account for new mutation)
   - Search trials targeting TP53 (Trial Matching)
   - Update care plan (Care Plan - new regimen)
   - Increase monitoring frequency (Monitoring)
```

---

### **Example 3: TMB-High Detected**

```
1. TEST: TMB = 15.2 mut/Mb
   ↓
2. PILLAR: Immune Status (Pillar 3)
   ↓
3. SIGNAL: TMB ≥10 (TMB-High)
   ↓
4. PATTERN: IO eligible (new treatment option unlocked)
   ↓
5. TRIGGER: tmb_high_detected (MEDIUM severity)
   ↓
6. ACTIONS:
   - Update IO eligibility (Biomarker Intelligence)
   - Re-rank drugs (Drug Efficacy - boost checkpoint inhibitors)
   - Match IO trials (Trial Matching - mechanism fit for IO pathway)
   - Notify oncologist
```

---

## 🔗 EXISTING CONNECTIONS

### **1. Biomarker Intelligence → Trigger System → Actions**
- ✅ CA-125 rise → resistance_detected trigger → ctDNA order + trial search
- ✅ TMB ≥10 → tmb_high_detected trigger → IO eligibility + IO trials
- ✅ MSI-H → msi_high_detected trigger → IO eligibility + Lynch screening

### **2. Resistance Prediction → Drug Efficacy → Trial Matching**
- ✅ DIS3 detected → alternatives (carfilzomib, daratumumab) → drug_efficacy handoff
- ✅ NF1 detected → alternatives (olaparib, trametinib) → drug_efficacy handoff
- ✅ Mechanism vector extracted from drug efficacy → trial matching with mechanism fit

### **3. Next Test Recommender → State Management → Pattern Recognition**
- ✅ HRD missing → recommend HRD test → track HRD score over time
- ✅ ctDNA missing → recommend ctDNA → track mutations over time
- ✅ Historical tracking enables pattern recognition (pathway changes, resistance signals)

### **4. Trigger System → State Management → Historical Analysis**
- ✅ CA-125 trends tracked → baseline comparison → resistance detection
- ✅ Pathway burden tracked: DDR 0.88 → 0.73 → 0.65 → resistance signal
- ✅ Mechanism vector changes tracked → pathway escape detection

---

## 🎯 PATTERN RECOGNITION

### **Resistance Emerging Pattern**

**Conditions:**
- CA-125 rising after nadir
- Pathway burden dropping
- New mutations in ctDNA

**Confidence:** HIGH

**Actions:**
- Order ctDNA sequencing
- Search salvage trials
- Escalate to tumor board

---

### **Molecular Progression Pattern**

**Conditions:**
- Rising VAF in ctDNA
- New driver mutations

**Confidence:** MEDIUM

**Actions:**
- Update resistance prediction
- Re-rank drugs

---

### **IO Eligible Pattern**

**Conditions:**
- TMB ≥10 OR MSI-High

**Confidence:** HIGH

**Actions:**
- Boost checkpoint inhibitors
- Match IO trials

---

## 🔗 Related Files

**Strategic Framework:**
- [03_GENERALS_BATTLE_MAP.mdc](03_GENERALS_BATTLE_MAP.mdc) - 6 Pillars framework

**SAE Capabilities:**
- [05_SAE_CAPABILITIES.md](05_SAE_CAPABILITIES.md) - SAE capabilities mapped to 6 Pillars

**Strategic Vision:**
- [06_STRATEGIC_VISION.md](06_STRATEGIC_VISION.md) - Strategic vision, next steps

---

*Document Owner: Zo*  
*Last Updated: January 28, 2025*  
*Status: ✅ ACTIVE - Intelligence flow documented*



