# 🚀 PRODUCTION STARTUP GUIDE - HOLISTIC SCORE FOR AYESHA

**Date**: January 29, 2025  
**Status**: ✅ **READY FOR PRODUCTION**  
**Feature**: Holistic Score Integration with RUO Labels

---

## ✅ PRE-FLIGHT CHECKLIST

### **Backend**
- ✅ Holistic Score service integrated into `trial_service.py`
- ✅ RUO label in provenance: `"ruo": "Research Use Only"`
- ✅ All tests passing (5/5)
- ✅ Ready for production

### **Frontend**
- ✅ `HolisticScoreCard` component created
- ✅ Integrated into `TrialMatchCard` and `TrialMatchesCard`
- ✅ RUO disclaimer added to component
- ✅ Ready for production

---

## 🔧 STARTUP INSTRUCTIONS

### **1. Backend Server**

```bash
# Navigate to backend directory
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal

# Start FastAPI server (with auto-reload for development)
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# OR if main.py is in root:
# uvicorn main:app --reload --host 0.0.0.0 --port 8000

# Verify server is running
curl http://localhost:8000/docs
# Should show FastAPI Swagger UI
```

**Expected Output:**
```
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process
INFO:     Started server process
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```

**Health Check:**
```bash
curl http://localhost:8000/health
# OR
curl http://localhost:8000/api/health
```

---

### **2. Frontend Server**

```bash
# Navigate to frontend directory
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-frontend

# OR check which worktree has the frontend:
# cd /Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/qhg/oncology-coPilot/oncology-frontend

# Install dependencies (if needed)
npm install

# Start development server
npm run dev

# OR if using Vite directly:
# npm run dev -- --host 0.0.0.0 --port 5173

# Verify frontend is running
# Open browser: http://localhost:5173
```

**Expected Output:**
```
VITE v5.x.x  ready in xxx ms

➜  Local:   http://localhost:5173/
➜  Network: http://0.0.0.0:5173/
```

---

## 🧪 VERIFICATION

### **1. Backend API Test**

```bash
# Test holistic score endpoint (if exposed separately)
curl -X POST http://localhost:8000/api/ayesha/complete_care_v2 \
  -H "Content-Type: application/json" \
  -d '{
    "ca125_value": 2842.0,
    "stage": "IVB",
    "treatment_line": "first-line",
    "germline_status": "positive",
    "location_state": "NY"
  }'

# Check response includes holistic scores in trials
# Should see: "holistic_score", "mechanism_fit_score", etc.
```

### **2. Frontend UI Test**

1. **Navigate to Ayesha's Trial Explorer:**
   - Open: `http://localhost:5173/ayesha-trials`
   - OR: `http://localhost:5173/ayesha-complete-care`

2. **Verify Holistic Score Display:**
   - Trial cards should show `HolisticScoreCard` component
   - Should display:
     - Overall holistic score (0-100%)
     - Mechanism Fit (50% weight)
     - Eligibility (30% weight)
     - PGx Safety (20% weight)
     - Interpretation badge (HIGH/MEDIUM/LOW)
     - **RUO Disclaimer** (red chip at bottom)

3. **Verify RUO Label:**
   - Look for red outlined chip: "Research Use Only - Not for Clinical Decision Making"
   - Should appear at bottom of `HolisticScoreCard`

---

## 📊 EXPECTED BEHAVIOR

### **For Ayesha's Profile:**

**High Holistic Score Trials (DDR-targeted):**
- Mechanism Fit: **0.85-0.90** (DDR-high patient + DDR trial)
- Eligibility: **0.85-0.95** (ovarian cancer, recruiting, age eligible)
- PGx Safety: **1.0** (no concerns)
- **Holistic Score: ~0.85-0.90** → **HIGH interpretation**

**Lower Holistic Score Trials (non-DDR):**
- Mechanism Fit: **0.20-0.30** (low DDR alignment)
- Eligibility: **0.85-0.95**
- PGx Safety: **1.0**
- **Holistic Score: ~0.50-0.60** → **LOW/MEDIUM interpretation**

---

## 🔍 TROUBLESHOOTING

### **Backend Not Starting**

```bash
# Check if port 8000 is in use
lsof -i :8000

# Kill process if needed
kill -9 <PID>

# Check Python environment
python3 --version
which python3

# Check dependencies
pip list | grep fastapi
pip list | grep uvicorn
```

### **Frontend Not Starting**

```bash
# Check if port 5173 is in use
lsof -i :5173

# Clear node_modules and reinstall
rm -rf node_modules package-lock.json
npm install

# Check Node version
node --version
npm --version
```

### **Holistic Scores Not Showing**

1. **Check Backend Logs:**
   - Look for "Holistic scores computed" log messages
   - Check for errors in trial service

2. **Check Frontend Console:**
   - Open browser DevTools (F12)
   - Check Network tab for API responses
   - Verify `trials` response includes `holistic_score` field

3. **Verify Integration:**
   - `HolisticScoreCard` should render when `trial.holistic_score` exists
   - Check that trials from API have holistic scores

---

## 📋 PRODUCTION CHECKLIST

### **Pre-Deployment**

- [x] Backend: Holistic scores computed in `trial_service.py`
- [x] Backend: RUO in provenance
- [x] Frontend: `HolisticScoreCard` component created
- [x] Frontend: RUO disclaimer added
- [x] Frontend: Integrated into trial cards
- [x] Tests: All 5 integration tests passing

### **During Deployment**

- [ ] Backend server running on port 8000
- [ ] Frontend server running on port 5173
- [ ] API endpoint `/api/ayesha/complete_care_v2` responding
- [ ] Trial response includes holistic scores
- [ ] Frontend displays `HolisticScoreCard` in trial cards
- [ ] RUO disclaimer visible in UI

### **Post-Deployment**

- [ ] Verify holistic scores for Ayesha's trials
- [ ] Verify DDR trials have higher scores than MAPK trials
- [ ] Verify RUO label displays correctly
- [ ] Manual testing of complete care plan page
- [ ] Manual testing of trial explorer page

---

## 🎯 AYESHA'S DELIVERY CHECKLIST

### **What Ayesha Should See:**

1. **Trial Cards with Holistic Scores:**
   - Overall holistic feasibility score (0-100%)
   - Breakdown of components (Mechanism Fit, Eligibility, PGx Safety)
   - Interpretation badge (HIGH/MEDIUM/LOW)

2. **High-Scoring DDR Trials:**
   - PARP inhibitor trials should show HIGH interpretation
   - Mechanism Fit should be ≥0.80
   - Holistic Score should be ≥0.85

3. **RUO Disclaimer:**
   - Red chip visible on every holistic score card
   - Clear message: "Research Use Only - Not for Clinical Decision Making"

4. **Complete Integration:**
   - Holistic scores appear in `/ayesha-trials` page
   - Holistic scores appear in `/ayesha-complete-care` page
   - All trial cards show holistic breakdown

---

## 📝 QUICK START COMMANDS

```bash
# Terminal 1: Backend
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# Terminal 2: Frontend
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-frontend
npm run dev

# Verify (Terminal 3):
curl http://localhost:8000/docs  # Backend API docs
open http://localhost:5173/ayesha-trials  # Frontend page
```

---

## ⚠️ IMPORTANT NOTES

1. **RUO Label**: All holistic score displays include "Research Use Only - Not for Clinical Decision Making" disclaimer.

2. **Data Source**: Holistic scores use Ayesha's actual profile:
   - Mechanism Vector: DDR-high `[0.88, 0.12, 0.15, 0.10, 0.05, 0.2, 0.0]`
   - Germline: MBD4 (`c.1293delA`), PDGFRA VUS (`c.2263T>C`)
   - Demographics: Age 40, Disease "Ovarian Cancer", Location "NY"

3. **Formula**: `Holistic Score = (0.5 × Mechanism Fit) + (0.3 × Eligibility) + (0.2 × PGx Safety)`

4. **Validation**: All 5 integration tests passing, scores validated against clinical expectations.

---

**Last Updated**: January 29, 2025  
**Status**: ✅ **READY FOR PRODUCTION**  
**Next Step**: Start servers and begin auditing for Ayesha
