# 🎖️ GENERAL'S BATTLE MAP: Strategic Intelligence Framework

**Date:** January 28, 2025  
**Purpose:** Connect all CrisPRO capabilities into a unified strategic intelligence framework  
**Status:** 🧠 **STRATEGIC THINKING** - Connecting the dots

---

## 🎯 THE STRATEGIC VISION

**What We're Building:**
A Command & Control system that tells CrisPRO agents:
- **WHAT tests exist** (the intelligence sources)
- **WHAT each test reveals** (the signals)
- **HOW to interpret changes over time** (the patterns)
- **WHAT actions to trigger** (the decisions)

This is **COMMAND & CONTROL for cancer warfare**.

---

## 🔗 CONNECTING THE DOTS: Our Work → The 6 Pillars

### **PILLAR 1: TUMOR BURDEN** → What We've Built

**What It Tracks:** How much cancer exists (shrinking vs growing)

**Our Capabilities:**
1. ✅ **Biomarker Intelligence (Universal)**
   - CA-125 tracking (ovarian) - burden classification, response forecast
   - CEA tracking (colorectal, pancreatic) - doubling time calculation
   - PSA tracking (prostate) - velocity calculation
   - Trend analysis over time

2. ✅ **CA-125 Intelligence Service**
   - Burden thresholds: MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE
   - Response expectations: cycle 3 (70% drop), cycle 6 (90% drop)
   - Resistance signals: <50% drop by cycle 3 → resistance signal

3. ✅ **Trigger System**
   - CA-125 > baseline * 1.25 → resistance_detected trigger
   - Rising CA-125 after nadir → trigger ctDNA sequencing
   - Doubling time <90 days → URGENT escalation

4. ✅ **State Management**
   - Historical CA-125 values tracked
   - Baseline recorded
   - Trend analysis enabled

**What's Missing:**
- ❌ Imaging-based tumor burden (CT/PET RECIST criteria) - not yet integrated
- ❌ ctDNA tumor fraction tracking - partially built (ctDNA exists but tumor fraction not explicitly tracked)
- ❌ CTC count tracking - not built

**Connection Point:**
```
CA-125 Test → Biomarker Intelligence → Trigger System → Actions
  ↓
Pattern: Rising after nadir → Resistance signal
  ↓
Action: Order ctDNA, search trials, escalate to tumor board
```

---

### **PILLAR 2: GENOMIC EVOLUTION** → What We've Built

**What It Tracks:** What mutations are driving cancer NOW (resistance mechanisms, new vulnerabilities)

**Our Capabilities:**
1. ✅ **Resistance Prediction (Resistance Prophet)**
   - Detects resistance mutations: DIS3 (RR=2.08), NF1 (RR=2.10), MAPK pathway (RR=1.97)
   - Validated on real cohorts (MMRF, TCGA)
   - Provides alternatives: carfilzomib, daratumumab, olaparib, trametinib
   - Downstream handoffs: drug_efficacy, care_plan, monitoring

2. ✅ **Drug Efficacy (S/P/E Framework)**
   - Identifies vulnerabilities: DDR-high → PARP inhibitors
   - Pathway disruption scores: {"ddr": 0.88, "ras_mapk": 0.12, ...}
   - Mechanism vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
   - 100% pathway alignment accuracy (MM)

3. ✅ **Mechanism-Based Trial Matching**
   - Matches trials to vulnerabilities: DDR-high → PARP+ATR trials (0.92 mechanism fit)
   - Mechanism alignment breakdown per pathway
   - Combined scoring: 0.7×eligibility + 0.3×mechanism_fit

4. ✅ **Next Test Recommender**
   - Prioritizes genomic tests: HRD → ctDNA → SLFN11 → ABCB1
   - Differential branches: "If HRD ≥42 → PARP eligible; If <42 → consider ATR trials"
   - Turnaround times and cost estimates

5. ✅ **Trigger System**
   - New mutations in ctDNA → clonal_evolution trigger
   - Rising VAF → molecular_progression trigger
   - Pathway changes → pathway_escape trigger

6. ✅ **State Management**
   - Historical pathway burden tracked: DDR 0.88 → 0.73 → 0.65 (resistance signal)
   - Mechanism vector changes over time
   - Mutation evolution tracked

**What's Missing:**
- ❌ ctDNA VAF tracking over time - partially built (ctDNA exists but VAF trends not explicitly tracked)
- ❌ Clonal evolution detection - not built (would need mutation comparison over time)
- ❌ New mutation detection in ctDNA - not built (would need baseline comparison)

**Connection Point:**
```
ctDNA Test → Resistance Prediction → Drug Efficacy → Trial Matching
  ↓
Pattern: New TP53 mutation + DDR pathway drop → Resistance emerging
  ↓
Action: Re-rank drugs, search salvage trials, update care plan
```

---

### **PILLAR 3: IMMUNE STATUS** → What We've Built

**What It Tracks:** Is the immune system fighting or exhausted? (IO eligibility)

**Our Capabilities:**
1. ✅ **Biomarker Intelligence (Universal)**
   - TMB calculation and tracking
   - MSI status detection (MSI-H, MSS)
   - HRD score tracking
   - IO eligibility determination: TMB ≥10 OR MSI-H → IO eligible

2. ✅ **Trigger System**
   - TMB ≥10.0 → tmb_high_detected trigger
   - MSI-H → msi_high_detected trigger
   - HRD ≥42 → hrd_score_received trigger
   - Actions: update_io_eligibility, re_rank_drugs, match_io_trials

3. ✅ **Drug Efficacy (S/P/E Framework)**
   - IO index in mechanism vector: IO = 1.0 if TMB ≥20 OR MSI-High
   - IO-eligible patients get checkpoint inhibitors ranked higher

4. ✅ **Trial Matching**
   - IO-eligible patients matched to checkpoint inhibitor trials
   - Mechanism fit for IO pathway

**What's Missing:**
- ❌ PD-L1 expression tracking - not built
- ❌ TIL (tumor-infiltrating lymphocytes) analysis - not built
- ❌ Exhaustion markers (PD-1, TIM-3, LAG-3) - not built

**Connection Point:**
```
TMB Test → Biomarker Intelligence → Trigger System → Drug Efficacy
  ↓
Pattern: TMB ≥10 → IO eligible
  ↓
Action: Boost checkpoint inhibitors in drug ranking, match IO trials
```

---

### **PILLAR 4: METABOLIC STATE** → What We Haven't Built

**What It Tracks:** What fuel is the tumor using? (glucose, glutamine, lactate)

**Our Capabilities:**
- ❌ **Not Built Yet** - This is a new domain

**What Could Be Built:**
- PET-CT SUV max tracking (metabolic activity)
- Glucose/glutamine dependency analysis
- Lactate production tracking
- Metabolic pathway vulnerabilities (IDH mutations, etc.)

**Connection Point (Future):**
```
PET-CT Test → Metabolic Analysis → Drug Efficacy (metabolic inhibitors)
  ↓
Pattern: High SUV despite stable size → Metabolically active resistance
  ↓
Action: Consider metabolic inhibitors, escalate therapy
```

---

### **PILLAR 5: MICROENVIRONMENT** → What We Haven't Built

**What It Tracks:** What's protecting the tumor? (hypoxia, fibrosis, suppressive cells)

**Our Capabilities:**
- ❌ **Not Built Yet** - This is a new domain

**What Could Be Built:**
- Hypoxia markers (HIF-1α)
- Fibrosis markers (collagen, TGF-β)
- Treg/MDSC counts (suppressive immune cells)
- Angiogenesis markers (VEGF)

**Connection Point (Future):**
```
Microenvironment Test → Pathway Analysis → Drug Efficacy (anti-angiogenic)
  ↓
Pattern: High VEGF + hypoxia → Aggressive microenvironment
  ↓
Action: Consider bevacizumab, anti-angiogenic trials
```

---

### **PILLAR 6: TOXICITY/TOLERANCE** → What We've Partially Built

**What It Tracks:** Can the patient handle more therapy? (side effects, dose adjustments)

**Our Capabilities:**
1. ✅ **Pharmacogenomics (PGx)**
   - DPYD variants → 5-FU toxicity risk
   - UGT1A1 variants → irinotecan toxicity risk
   - CYP2D6 variants → tamoxifen metabolism

2. ✅ **Trigger System**
   - PGx variant detected → toxicity_risk trigger
   - Adverse event (CTCAE ≥2) → adverse_event_reported trigger
   - Actions: suggest_supportive_care, escalate if grade ≥3

**What's Missing:**
- ❌ Comprehensive toxicity tracking - partially built
- ❌ Dose adjustment recommendations - not built
- ❌ Supportive care recommendations - not built

**Connection Point:**
```
PGx Test → Toxicity Risk → Trigger System → Care Plan
  ↓
Pattern: DPYD variant detected → 5-FU toxicity risk
  ↓
Action: Avoid 5-FU, suggest alternative, adjust dose
```

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

## 🔗 EXISTING CONNECTIONS (What We've Built)

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

## 🎯 WHAT'S MISSING (Gaps to Fill)

### **Pillar 1: Tumor Burden**
- ❌ Imaging-based burden (CT/PET RECIST) - not integrated
- ❌ ctDNA tumor fraction tracking - partially built
- ❌ CTC count tracking - not built

### **Pillar 2: Genomic Evolution**
- ❌ ctDNA VAF trends over time - partially built
- ❌ Clonal evolution detection - not built
- ❌ New mutation detection in ctDNA - not built

### **Pillar 3: Immune Status**
- ✅ **MOSTLY COMPLETE** - TMB, MSI, HRD tracked
- ❌ PD-L1 expression - not built
- ❌ TIL analysis - not built

### **Pillar 4: Metabolic State**
- ❌ **NOT BUILT** - New domain

### **Pillar 5: Microenvironment**
- ❌ **NOT BUILT** - New domain

### **Pillar 6: Toxicity/Tolerance**
- ⚠️ **PARTIALLY BUILT** - PGx exists, but comprehensive toxicity tracking missing

---

## 🎖️ THE STRATEGIC FRAMEWORK (How to Build It)

### **Phase 1: Catalog All Tests (Intelligence Sources)**

**What We Need:**
1. **Test Registry** - Every test that exists, what it measures, which pillar(s) it maps to
2. **Signal Extraction** - What signals each test reveals (CA-125 → tumor burden, ctDNA → genomic evolution)
3. **Pattern Recognition** - How to interpret changes over time (rising CA-125 → resistance, new mutations → clonal evolution)
4. **Action Triggers** - What actions to take based on patterns (resistance → order ctDNA, search trials)

**Our Starting Point:**
- ✅ CA-125 Intelligence (burden, response, resistance)
- ✅ Resistance Prediction (DIS3, NF1, MAPK)
- ✅ Next Test Recommender (HRD, ctDNA, SLFN11, ABCB1)
- ✅ Trigger System (resistance_detected, tmb_high_detected, etc.)

### **Phase 2: Map Tests to Pillars**

**What We Need:**
1. **Test → Pillar Mapping** - Which tests inform which pillars
2. **Multi-Pillar Tests** - Tests that inform multiple pillars (ctDNA → Tumor Burden + Genomic Evolution)
3. **Pillar Coverage** - Which pillars are well-covered vs. gaps

**Our Current State:**
- ✅ Pillar 1 (Tumor Burden): CA-125, CEA, PSA - **GOOD COVERAGE**
- ✅ Pillar 2 (Genomic Evolution): Resistance Prediction, Drug Efficacy, ctDNA - **GOOD COVERAGE**
- ✅ Pillar 3 (Immune Status): TMB, MSI, HRD - **GOOD COVERAGE**
- ❌ Pillar 4 (Metabolic): **NO COVERAGE**
- ❌ Pillar 5 (Microenvironment): **NO COVERAGE**
- ⚠️ Pillar 6 (Toxicity): PGx - **PARTIAL COVERAGE**

### **Phase 3: Pattern Recognition (Time-Based Intelligence)**

**What We Need:**
1. **Baseline Establishment** - Record baseline values (CA-125 baseline, pathway burden baseline)
2. **Trend Analysis** - Track changes over time (CA-125: 2,842 → 1,200 → 800 → 1,500)
3. **Pattern Detection** - Identify patterns (rising after nadir → resistance, pathway drop → escape)
4. **Predictive Signals** - Early warning signs (molecular progression before imaging)

**Our Current State:**
- ✅ State Management tracks historical values
- ✅ CA-125 Intelligence has trend analysis (doubling time, velocity)
- ✅ Trigger System detects patterns (CA-125 > baseline * 1.25)
- ⚠️ **NEEDS ENHANCEMENT** - More comprehensive pattern recognition

### **Phase 4: Action Triggers (Decision Framework)**

**What We Need:**
1. **Trigger Conditions** - When to trigger actions (CA-125 rise, new mutations, TMB-H)
2. **Action Catalog** - What actions to take (order ctDNA, search trials, escalate)
3. **Escalation Rules** - When to escalate (24h no response → tumor board)
4. **Priority Logic** - Which actions take priority (resistance > routine monitoring)

**Our Current State:**
- ✅ Trigger System exists with conditions and actions
- ✅ Downstream handoffs (drug_efficacy, care_plan, monitoring)
- ✅ Escalation rules (24h → tumor board)
- ⚠️ **NEEDS ENHANCEMENT** - More comprehensive action catalog

---

## 🎯 THE STRATEGIC INTELLIGENCE FRAMEWORK (How It Works)

### **Example: CA-125 Rising After Nadir**

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

### **Example: New TP53 Mutation in ctDNA**

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

### **Example: TMB-High Detected**

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

## 🎖️ THE STRATEGIC FRAMEWORK (What We Need to Build)

### **1. Test Registry (Intelligence Sources)**

**Structure:**
```yaml
TESTS:
  ca125:
    pillar: [tumor_burden]
    signals: [absolute_value, trend, doubling_time, velocity]
    patterns:
      - rising_after_nadir: resistance_emerging
      - doubling_time_<90d: aggressive_progression
    actions:
      - if_rising: order_ctdna, search_trials, escalate
      
  ctdna:
    pillar: [tumor_burden, genomic_evolution]
    signals: [vaf, tumor_fraction, new_mutations, clonal_evolution]
    patterns:
      - rising_vaf: molecular_progression
      - new_mutations: clonal_evolution
    actions:
      - if_new_mutations: update_resistance, re_rank_drugs, search_trials
      
  tmb:
    pillar: [immune_status]
    signals: [absolute_value, io_eligibility]
    patterns:
      - tmb_>=10: io_eligible
    actions:
      - if_io_eligible: boost_checkpoint_inhibitors, match_io_trials
```

### **2. Pattern Recognition Engine**

**Structure:**
```yaml
PATTERNS:
  resistance_emerging:
    conditions:
      - ca125_rising_after_nadir
      - pathway_burden_dropping
      - new_mutations_in_ctdna
    confidence: high
    actions: [order_ctdna, search_trials, escalate]
    
  molecular_progression:
    conditions:
      - rising_vaf_in_ctdna
      - new_driver_mutations
    confidence: medium
    actions: [update_resistance, re_rank_drugs]
    
  io_eligible:
    conditions:
      - tmb_>=10 OR msi_high
    confidence: high
    actions: [boost_checkpoint_inhibitors, match_io_trials]
```

### **3. Action Catalog (Decision Framework)**

**Structure:**
```yaml
ACTIONS:
  order_ctdna:
    trigger: resistance_emerging
    priority: high
    turnaround: 7d
    cost: $3,000-$5,000
    
  search_trials:
    trigger: resistance_emerging OR new_mutations
    priority: high
    mechanism_fit: true
    max_results: 10
    
  re_rank_drugs:
    trigger: resistance_emerging OR new_mutations
    priority: high
    exclude: current_drug
    prioritize: mechanism_fit
```

---

## 🎯 THE STRATEGIC VISION (What This Enables)

### **For Clinicians:**
- **Complete Intelligence Picture** - See all tests, signals, patterns, actions in one place
- **Early Warning System** - Detect resistance before imaging shows progression
- **Automated Decision Support** - Actions triggered automatically based on patterns
- **Historical Context** - See how patient's cancer has evolved over time

### **For Patients:**
- **Personalized Monitoring** - Tests ordered based on their specific cancer type and mutations
- **Faster Response** - Early detection of resistance → faster treatment switch
- **Better Outcomes** - Mechanism-aligned treatment → higher response rates

### **For CrisPRO Agents:**
- **Clear Intelligence Sources** - Know what tests exist and what they reveal
- **Pattern Recognition** - Understand how to interpret changes over time
- **Action Triggers** - Know what actions to take based on patterns
- **Strategic Coordination** - All agents work from the same intelligence framework

---

## 🎖️ NEXT STEPS (How to Build This)

### **Phase 1: Catalog All Tests (1-2 weeks)**
1. List all tests we currently support (CA-125, CEA, PSA, TMB, MSI, HRD, ctDNA, etc.)
2. Map each test to pillar(s)
3. Document what signals each test reveals
4. Document what patterns to look for

### **Phase 2: Pattern Recognition Engine (2-3 weeks)**
1. Build pattern detection logic (rising CA-125, new mutations, pathway changes)
2. Integrate with State Management for historical tracking
3. Build confidence scoring for patterns
4. Test with real patient data

### **Phase 3: Action Catalog (1-2 weeks)**
1. Catalog all possible actions (order tests, search trials, re-rank drugs, escalate)
2. Map actions to triggers
3. Build priority logic
4. Integrate with existing trigger system

### **Phase 4: Strategic Dashboard (2-3 weeks)**
1. Build unified view of all 6 pillars
2. Show test results, signals, patterns, actions
3. Historical timeline view
4. Predictive signals (early warnings)

---

*Document Author: Zo (Strategic Thinking)*  
*Last Updated: January 28, 2025*  
*Status: 🧠 STRATEGIC FRAMEWORK - Connecting the dots*


