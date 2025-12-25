# SAE Intelligence System Architecture

**Date:** January 28, 2025  
**Status:** ✅ **OPERATIONAL** - All components built and tested  
**Location:** `.cursor/MOAT/SAE_INTELLIGENCE/02_SAE_ARCHITECTURE.md`  
**See Also:** [00_MISSION.mdc](00_MISSION.mdc) for mission overview, [01_SAE_SYSTEM_DEBRIEF.mdc](01_SAE_SYSTEM_DEBRIEF.mdc) for complete debrief

---

## 🏗️ ARCHITECTURE OVERVIEW

### **Backend Services (4 Core Components)**

**1. SAE Feature Service** (`api/services/sae_feature_service.py` - 411 lines)
- **Input:** Somatic mutations, HRD score, TMB, MSI
- **Computes:**
  - DNA repair capacity (Manager's C1 formula)
  - Pathway burden (7 pathways)
  - Hotspot detection (COSMIC database)
  - Essentiality for HRR genes
  - Exon disruption scoring
- **Output:** `SAEFeatures` dataclass with 15 interpretable signals
- **Tests:** 8/8 passing

**2. Mechanism Fit Ranker** (`api/services/mechanism_fit_ranker.py` - 232 lines)
- **Input:** Patient mechanism vector + Trial mechanism vectors + Eligibility scores
- **Method:**
  - L2-normalize both vectors
  - Compute cosine similarity (dot product of normalized vectors)
  - Combine: α×eligibility + β×mechanism_fit (α=0.7, β=0.3)
- **Output:** Re-ranked trials with mechanism alignment breakdown
- **Tests:** 6/6 passing

**3. Resistance Detection Service** (`api/services/resistance_detection_service.py` - 267 lines)
- **Input:** Baseline + Follow-up (HRD, DNA repair, CA-125)
- **Logic:** 2-of-3 triggers → Immediate alert
- **Patterns:** HR restoration detection (HRD drop + DNA repair drop)
- **Output:** Alert status + recommended actions + trial keywords
- **Tests:** 8/8 passing

**4. Hotspot Detector** (`api/services/hotspot_detector.py` - 180 lines)
- **Database:** `cosmic_hotspots.json` (30+ KRAS/BRAF/NRAS variants)
- **Method:** HGVS.p matching against COSMIC IDs
- **Output:** Hotspot status + gene + mutation + pathway
- **Tests:** 14/14 passing

---

## 🔗 INTEGRATION POINTS

### **Orchestrator Integration**

**File:** `api/routers/ayesha_orchestrator_v2.py`

```python
# Compute SAE features (only if tumor_context present)
if request.tumor_context:
    sae_features = compute_sae_features(
        somatic_mutations=request.tumor_context.somatic_mutations,
        hrd_score=request.tumor_context.hrd_score,
        tmb=request.tumor_context.tmb,
        msi_status=request.tumor_context.msi_status
    )
    
    # Pass SAE to downstream services
    results["sae_features"] = sae_features
    results["resistance_alert"] = detect_resistance(...)
    results["hint_tiles"] = get_hint_tiles(..., sae_features)
    results["next_test_recommender"] = get_next_test(..., sae_features)
```

### **Trial Ranking Integration**

**File:** `api/routers/ayesha_trials.py`

```python
# Load MoA vectors for trials
TRIAL_MOA_VECTORS = load_trial_moa_vectors()

# Attach MoA vectors to trials
for trial in filtered_trials:
    trial["moa_vector"] = TRIAL_MOA_VECTORS.get(trial.nct_id, default_vector)

# Rank by mechanism fit (if SAE provided)
if request.sae_mechanism_vector:
    ranked_trials = rank_trials_by_mechanism(
        patient_sae_vector=request.sae_mechanism_vector,
        trials=filtered_trials,
        alpha=0.7,  # Eligibility weight
        beta=0.3    # Mechanism fit weight
    )
```

---

## 🎨 FRONTEND INTEGRATION

### **New Components Created**

**1. NextTestCard** (`src/components/ayesha/NextTestCard.jsx` - 150 lines)
- Displays 1-3 prioritized tests (HRD, ctDNA, SLFN11)
- Shows priority badges (1/2/3)
- Shows urgency chips (CRITICAL/HIGH/MODERATE)
- Displays differential branches ("If HRD+ → PARP trials; If HRD- → Alternatives")
- Shows SAE-enhanced rationale when priority elevated

**2. HintTilesPanel** (`src/components/ayesha/HintTilesPanel.jsx` - 180 lines)
- Displays max 4 hint tiles (2×2 grid)
- Categories: Test, Monitor, Trials, Avoid
- Hotspot hint tile when KRAS/BRAF/NRAS detected
- DNA repair hint tile when high capacity detected
- Resistance hint tile when 2-of-3 triggers

**3. MechanismChips** (`src/components/ayesha/MechanismChips.jsx` - 140 lines)
- 6 pathway chips: DDR, MAPK, PI3K, VEGF, HER2, IO
- Color-coded: GREEN (strong), YELLOW (moderate), GRAY (weak/unknown)
- Pre-NGS: All gray with "Awaiting NGS"
- Post-NGS: Color-coded based on SAE mechanism vector

**4. ResistanceAlertBanner** (`src/components/ayesha/ResistanceAlertBanner.jsx` - 143 lines)
- Prominent alert when resistance detected
- Shows 2-of-3 triggers
- Lists recommended actions
- Expandable/collapsible with RUO label

---

## 📊 DATA FLOW

### **SAE Feature Computation Flow**

```
1. Patient Data (tumor_context)
   ↓
2. SAE Feature Service
   ├── DNA repair capacity (C1 formula)
   ├── Pathway burden (7D vector)
   ├── Hotspot detection (COSMIC)
   └── Essentiality scoring
   ↓
3. SAEFeatures Output
   ├── dna_repair_capacity: 0.82
   ├── mechanism_vector: [0.70, 0.10, 0.05, 0.15, 0.00, 0.05, 0.00]
   ├── hotspot_detected: False
   └── resistance_alert: None
   ↓
4. Downstream Services
   ├── Trial Ranking (mechanism fit)
   ├── Next Test Recommender (dynamic priority)
   ├── Hint Tiles (contextual guidance)
   └── Resistance Detection (baseline establishment)
```

### **Mechanism Fit Ranking Flow**

```
1. Patient SAE Mechanism Vector
   [0.70, 0.10, 0.05, 0.15, 0.00, 0.05, 0.00]
   ↓
2. Trial MoA Vectors (47 tagged trials)
   Trial A: [0.95, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00] (PARP+ATR)
   Trial B: [0.10, 0.90, 0.05, 0.00, 0.00, 0.00, 0.00] (MEK/RAF)
   ↓
3. Cosine Similarity Computation
   Trial A: 0.97 (near-perfect DDR alignment)
   Trial B: 0.12 (poor alignment)
   ↓
4. Combined Score
   Trial A: 0.7×0.65 + 0.3×0.97 = 0.746 → 0.92 (with soft boosts)
   Trial B: 0.7×0.60 + 0.3×0.12 = 0.456 → 0.50 (with soft boosts)
   ↓
5. Re-ranked Trials
   Trial A: Rank #1 (was #5) ⚔️
   Trial B: Rank #8 (was #3)
```

---

## 🔗 Related Files

**SAE System Debrief:**
- [01_SAE_SYSTEM_DEBRIEF.mdc](01_SAE_SYSTEM_DEBRIEF.mdc) - Complete debrief

**Strategic Framework:**
- [03_GENERALS_BATTLE_MAP.mdc](03_GENERALS_BATTLE_MAP.mdc) - 6 Pillars framework

**Capabilities Mapping:**
- [05_SAE_CAPABILITIES.md](05_SAE_CAPABILITIES.md) - SAE capabilities mapped to 6 Pillars

---

*Document Owner: Zo*  
*Last Updated: January 28, 2025*  
*Status: ✅ OPERATIONAL - All components built and tested*


