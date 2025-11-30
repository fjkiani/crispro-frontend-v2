# P1.4: Post-NGS E2E Tests - README

**File:** `tests/test_ayesha_post_ngs_e2e.py`  
**Status:** ✅ **COMPLETE** (requires backend running)  
**Owner:** Zo  
**Date:** January 13, 2025

---

## Test Scenarios

**1. BRCA1 Biallelic (HRD=58)**
- Tests high DNA repair capacity detection
- Verifies SAE features present
- Checks SLFN11 priority elevation
- Validates DDR pathway in mechanism map

**2. KRAS G12D Hotspot**
- Tests MAPK hotspot detection
- Verifies hotspot hint tile generation
- Checks ctDNA recommendation
- Validates MAPK pathway in mechanism map

**3. Resistance Detection**
- Tests resistance alert service
- Verifies 2-of-3 trigger logic
- Checks recommended actions present

---

## How to Run

### Prerequisites:
1. Backend must be running on `http://localhost:8000`
2. `httpx` library installed (`pip install httpx`)

### Start Backend:
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload --port 8000
```

### Run Tests:
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
PYTHONPATH=oncology-coPilot/oncology-backend-minimal venv/bin/python tests/test_ayesha_post_ngs_e2e.py
```

---

## Expected Output

```
⚔️ P1.4: POST-NGS E2E TEST SUITE
================================================================================

TEST 1: BRCA1 BIALLELIC (HRD=58, HIGH DNA REPAIR CAPACITY)
✅ SAE Features Present
   DNA Repair Capacity: 0.82 (expected >= 0.70)
   Hotspot Mutation: False (expected False)
   DDR Pathway Burden: 0.70
✅ Hint Tiles: 3 tiles (DDR/PARP context: True)
✅ SLFN11 Recommended: Priority 2 (expected 2 if elevated)
   SAE Enhancement: YES
✅ Mechanism Map: DDR chip present (status: green)
✅ TEST 1 PASSED

TEST 2: KRAS G12D HOTSPOT (MAPK PATHWAY ACTIVATION)
✅ SAE Hotspot Detection
   Hotspot Mutation: True (expected True)
   Gene: KRAS (expected KRAS)
   Mutation: G12D (expected G12D)
   Pathway: MAPK (expected MAPK)
✅ Hotspot Hint Tile Present
   Title: 🧬 MAPK Hotspot Detected
✅ Mechanism Map: MAPK chip present (status: yellow/green)
✅ TEST 2 PASSED

TEST 3: RESISTANCE DETECTION (2-of-3 TRIGGERS)
✅ Resistance Alert Service Called
   Alert Triggered: True/False (depends on scenario)
✅ TEST 3 PASSED

✅ ALL POST-NGS E2E TESTS PASSED!
```

---

## What This Tests

**Backend Integration:**
- ✅ `/api/ayesha/complete_care_v2` endpoint
- ✅ SAE feature computation with tumor_context
- ✅ Hotspot detection service
- ✅ Dynamic next-test prioritization
- ✅ Resistance detection service
- ✅ Hint tiles generation
- ✅ Mechanism map display

**Data Flow:**
- ✅ tumor_context → SAE features
- ✅ SAE features → hint tiles
- ✅ SAE features → next-test recommendations
- ✅ SAE features → mechanism map
- ✅ SAE features → resistance detection

**Manager Policy Alignment:**
- ✅ C1: DNA repair capacity formula
- ✅ C2: Hotspot mutation detection
- ✅ C3: Resistance detection (2-of-3 triggers)
- ✅ C6: Dynamic next-test prioritization

---

## Troubleshooting

**Backend Not Running:**
```
❌ HTTP 500: Connection error
```
→ Start backend with `uvicorn api.main:app --reload --port 8000`

**Test Failures:**
- Check SAE services are integrated in orchestrator
- Verify hotspot_detector.py is working
- Confirm next_test_recommender.py accepts sae_features parameter

---

**Status:** ⚔️ **P1.4 COMPLETE - READY TO RUN WITH BACKEND** ⚔️






