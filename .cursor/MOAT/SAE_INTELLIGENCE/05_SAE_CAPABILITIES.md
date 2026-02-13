# SAE Capabilities Mapped to 6 Pillars of Cancer Intelligence

**Date:** January 28, 2025  
**Status:** ✅ **ACTIVE** - SAE capabilities mapped to 6 Pillars  
**Location:** `.cursor/MOAT/SAE_INTELLIGENCE/05_SAE_CAPABILITIES.md`  
**See Also:** [00_MISSION.mdc](00_MISSION.mdc) for mission overview, [03_GENERALS_BATTLE_MAP.mdc](03_GENERALS_BATTLE_MAP.mdc) for 6 Pillars framework

---

## 🔗 SAE CAPABILITIES → 6 PILLARS MAPPING

### **PILLAR 1: TUMOR BURDEN** → SAE Contribution

**What SAE Provides:**
- ✅ **CA-125 Intelligence Integration**
  - SAE resistance detection uses CA-125 trends (2-of-3 triggers)
  - CA-125 inadequate response (<50% drop by Cycle 3) → resistance trigger
  - Baseline CA-125 tracked for resistance monitoring

**What's Missing:**
- ❌ Imaging-based tumor burden (CT/PET RECIST) - not integrated
- ❌ ctDNA tumor fraction tracking - partially built
- ❌ CTC count tracking - not built

**Connection Point:**
```
CA-125 Test → SAE Resistance Detection → Trigger System → Actions
  ↓
Pattern: Rising after nadir → Resistance signal
  ↓
Action: Order ctDNA, search trials, escalate to tumor board
```

---

### **PILLAR 2: GENOMIC EVOLUTION** → SAE Contribution ⭐ **PRIMARY PILLAR**

**What SAE Provides:**
1. ✅ **DNA Repair Capacity** (Manager's C1 Formula)
   - Computes DNA repair capacity from mutations
   - Formula: `0.6×DDR_pathway + 0.2×HRR_essentiality + 0.2×exon_disruption`
   - Example: BRCA1 biallelic → DNA repair = 0.82 (HIGH)

2. ✅ **Pathway Burden** (7D Mechanism Vector)
   - Computes 7D mechanism vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
   - Source: **Proxy Logic** (Validated Genes)

3. ✅ **Hotspot Detection** (COSMIC Database)
   - Detects KRAS/BRAF/NRAS hotspots
   - Maps to pathways (KRAS → MAPK, BRAF → MAPK)

4. ✅ **Transcriptomic Risk (L1)**
   - **MFAP4 Biomarker**: High expression -> EMT Phenotype.
   - **AUROC**: 0.763 (Gold Standard).

5. ✅ **Resistance Detection** (2-of-3 Triggers)
   - HRD score drop ≥10 points
   - MFAP4 Spike (New)
   - CA-125 inadequate response

**What's Missing:**
- ❌ ctDNA VAF tracking over time - partially built
- ❌ Clonal evolution detection - not built
- ❌ New mutation detection in ctDNA - not built

**Connection Point:**
```
ctDNA Test → SAE Feature Service → Mechanism Vector → Trial Matching
  ↓
Pattern: New TP53 mutation + DDR pathway drop → Resistance emerging
  ↓
Action: Re-rank drugs, search salvage trials, update care plan
```

---

### **PILLAR 3: IMMUNE STATUS** → SAE Contribution

**What SAE Provides:**
- ✅ **IO Pathway in Mechanism Vector**
  - IO index in 7D mechanism vector: IO = 1.0 if TMB ≥20 OR MSI-High
  - IO-eligible patients get checkpoint inhibitors ranked higher
  - Mechanism fit for IO pathway in trial matching

**What's Missing:**
- ❌ PD-L1 expression tracking - not built
- ❌ TIL (tumor-infiltrating lymphocytes) analysis - not built
- ❌ Exhaustion markers (PD-1, TIM-3, LAG-3) - not built

**Connection Point:**
```
TMB Test → SAE Mechanism Vector (IO index) → Drug Efficacy → Trial Matching
  ↓
Pattern: TMB ≥10 → IO eligible
  ↓
Action: Boost checkpoint inhibitors in drug ranking, match IO trials
```

---

### **PILLAR 4: METABOLIC STATE** → SAE Contribution

**What SAE Provides:**
- ❌ **Not Built Yet** - This is a new domain

**What Could Be Built:**
- PET-CT SUV max tracking (metabolic activity)
- Glucose/glutamine dependency analysis
- Lactate production tracking
- Metabolic pathway vulnerabilities (IDH mutations, etc.)

**Connection Point (Future):**
```
PET-CT Test → Metabolic Analysis → SAE Metabolic Features → Drug Efficacy
  ↓
Pattern: High SUV despite stable size → Metabolically active resistance
  ↓
Action: Consider metabolic inhibitors, escalate therapy
```

---

### **PILLAR 5: MICROENVIRONMENT** → SAE Contribution

**What SAE Provides:**
- ❌ **Not Built Yet** - This is a new domain

**What Could Be Built:**
- Hypoxia markers (HIF-1α)
- Fibrosis markers (collagen, TGF-β)
- Treg/MDSC counts (suppressive immune cells)
- Angiogenesis markers (VEGF) - partially in mechanism vector

**Connection Point (Future):**
```
Microenvironment Test → SAE Microenvironment Features → Drug Efficacy
  ↓
Pattern: High VEGF + hypoxia → Aggressive microenvironment
  ↓
Action: Consider bevacizumab, anti-angiogenic trials
```

---

### **PILLAR 6: TOXICITY/TOLERANCE** → SAE Contribution

**What SAE Provides:**
- ⚠️ **Indirect Contribution**
  - SAE doesn't directly compute toxicity, but:
  - Resistance detection can trigger dose adjustments
  - Pathway burden changes can indicate tolerance issues

**What's Missing:**
- ❌ Comprehensive toxicity tracking - partially built
- ❌ Dose adjustment recommendations - not built
- ❌ Supportive care recommendations - not built

**Connection Point:**
```
PGx Test → Toxicity Risk → SAE Resistance Detection → Care Plan
  ↓
Pattern: DPYD variant detected → 5-FU toxicity risk
  ↓
Action: Avoid 5-FU, suggest alternative, adjust dose
```

---

## 📊 SAE FEATURES BY PILLAR

| SAE Feature | Pillar(s) | Status | Impact |
|-------------|-----------|--------|--------|
| **DNA Repair Capacity** | Pillar 2 (Genomic Evolution) | ✅ Operational | HIGH - PARP sensitivity prediction |
| **Pathway Burden (7D)** | Pillar 2 (Genomic Evolution) | ✅ Operational | HIGH - Mechanism-based trial matching |
| **Hotspot Detection** | Pillar 2 (Genomic Evolution) | ✅ Operational | MEDIUM - Targeted therapy hints |
| **Resistance Detection** | Pillar 1 (Tumor Burden) + Pillar 2 (Genomic Evolution) | ✅ Operational | HIGH - Early resistance detection |
| **IO Pathway Index** | Pillar 3 (Immune Status) | ✅ Operational | MEDIUM - IO eligibility in mechanism vector |
| **VEGF Pathway** | Pillar 5 (Microenvironment) | ⚠️ Partial | LOW - In mechanism vector, not fully utilized |

---

## 🎯 STRATEGIC VALUE

### **SAE as Intelligence Source**

**For Pillar 2 (Genomic Evolution):**
- ✅ **Primary Intelligence Source** - SAE is the core capability
- ✅ **DNA Repair Capacity** - Predicts PARP sensitivity
- ✅ **Pathway Burden** - Guides mechanism-based matching
- ✅ **Resistance Detection** - Early warning system

**For Other Pillars:**
- ⚠️ **Supporting Role** - SAE contributes but not primary
- ⚠️ **Integration Points** - SAE features used by other systems
- ⚠️ **Future Expansion** - Can add metabolic/microenvironment features

---

## 🔗 Related Files

**SAE System Debrief:**
- [01_SAE_SYSTEM_DEBRIEF.mdc](01_SAE_SYSTEM_DEBRIEF.mdc) - Complete debrief

**Strategic Framework:**
- [03_GENERALS_BATTLE_MAP.mdc](03_GENERALS_BATTLE_MAP.mdc) - 6 Pillars framework

**Intelligence Flow:**
- [04_INTELLIGENCE_FLOW.md](04_INTELLIGENCE_FLOW.md) - Test → Signals → Patterns → Actions

---

*Document Owner: Zo*  
*Last Updated: January 28, 2025*  
*Status: ✅ ACTIVE - SAE capabilities mapped to 6 Pillars*


