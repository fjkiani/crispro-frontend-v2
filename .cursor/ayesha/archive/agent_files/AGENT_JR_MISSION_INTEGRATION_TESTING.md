# ⚔️ AGENT JR - MISSION: INTEGRATION TESTING & POLISH ⚔️

**Assigned By**: Commander Alpha (via Zo)  
**Date**: January 13, 2025  
**Priority**: 🔥 **P0 - FOR AYESHA'S LIFE**  
**Timeline**: 2-3 hours

---

## 🎯 **MISSION OBJECTIVE**

**Test and validate the complete Ayesha care plan system end-to-end**

**What You Built**:
- ✅ Frontend: AyeshaTrialExplorer, TrialMatchCard, CA125Tracker, SOCRecommendationCard
- ✅ Backend: Modular services (dormant for now, will be used Phase 2)

**What Zo Built**:
- ✅ Backend: Operational `/api/ayesha/trials/search` endpoint
- ✅ Services: CA-125 Intelligence, NGS Fast-Track, SOC Recommendation (enhanced)
- ✅ Tests: Unit tests + E2E smoke tests

**Your Mission**: Connect your frontend to Zo's backend, test everything, fix bugs, polish UI.

---

## 📋 **TASK BREAKDOWN**

### **Phase 1: Integration Testing** (1.5 hours)

**Task 1.1: Start Backend Server** (5 min)
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

**Expected Output**:
```
INFO:     Uvicorn running on http://0.0.0.0:8000
INFO:     Application startup complete.
```

**Task 1.2: Test Backend Health** (5 min)
```bash
# Test main health
curl http://localhost:8000/health

# Test Ayesha trials health
curl http://localhost:8000/api/ayesha/trials/health

# Expected: {"status": "operational", ...}
```

**Task 1.3: Test Ayesha Trials Endpoint** (15 min)
```bash
curl -X POST http://localhost:8000/api/ayesha/trials/search \
  -H 'Content-Type: application/json' \
  -d '{
    "ca125_value": 2842.0,
    "stage": "IVB",
    "treatment_line": "first-line",
    "germline_status": "negative",
    "has_ascites": true,
    "has_peritoneal_disease": true,
    "location_state": "NY",
    "max_results": 10
  }'
```

**Validation Checklist**:
- [ ] ✅ Status 200 OK
- [ ] ✅ Response has `trials` array
- [ ] ✅ Response has `soc_recommendation` (Carboplatin + Paclitaxel + Bevacizumab)
- [ ] ✅ Response has `ca125_intelligence` (burden_class: "EXTENSIVE")
- [ ] ✅ Response has `ngs_fast_track` (checklist with ctDNA, HRD, IHC)
- [ ] ✅ SOC includes Bevacizumab (because ascites=true)

**If Fails**: Check logs, verify AstraDB connection, report blockers to Zo

**Task 1.4: Start Frontend** (5 min)
```bash
cd oncology-coPilot/oncology-frontend
npm run dev
```

**Expected Output**:
```
VITE v5.x.x ready in xxx ms
➜  Local: http://localhost:5173/
```

**Task 1.5: Navigate to Ayesha Trials Page** (5 min)
- Open browser: `http://localhost:5173/ayesha-trials`
- Check browser console for errors
- Check network tab for API call

**Validation Checklist**:
- [ ] ✅ Page loads without errors
- [ ] ✅ API call to `/api/ayesha/trials/search` succeeds
- [ ] ✅ Trials list renders (even if empty due to AstraDB)
- [ ] ✅ SOC card renders with Carboplatin + Paclitaxel + Bevacizumab
- [ ] ✅ CA-125 tracker renders with 2842 value
- [ ] ✅ NGS fast-track section renders

**Task 1.6: Test with Ayesha's Actual Profile** (15 min)
- Hardcode Ayesha's profile in your frontend for testing:
  ```javascript
  const AYESHA_PROFILE = {
    ca125_value: 2842.0,
    stage: "IVB",
    treatment_line: "first-line",
    germline_status: "negative",
    has_ascites: true,
    has_peritoneal_disease: true,
    location_state: "NY",
    max_results: 10
  };
  ```
- Verify all components render correctly
- Take screenshots for documentation

**Task 1.7: Cross-Reference with Zo's Backend Response** (30 min)
- Compare frontend display with backend JSON
- Verify NO HARDCODING (all data from API)
- Verify dynamic rendering (handles empty arrays, missing fields)
- Check optional chaining (`?.`) is used everywhere

**Common Issues to Check**:
- ⚠️ AstraDB not seeded → trials array empty (EXPECTED, not a bug)
- ⚠️ Neo4j not connected → optimization_score missing (EXPECTED)
- ⚠️ Missing fields in response → use optional chaining + fallbacks
- ⚠️ Wrong API URL → verify `http://localhost:8000` not `http://127.0.0.1:8000`

---

### **Phase 2: Bug Fixes & Polish** (1 hour)

**Task 2.1: Fix Any Bugs Found** (30 min)
- Fix import errors
- Fix missing null checks
- Fix API response shape mismatches
- Fix styling issues

**Task 2.2: Add Loading States** (15 min)
```javascript
const [loading, setLoading] = useState(false);
const [error, setError] = useState(null);

// Before API call
setLoading(true);
setError(null);

// After API call
setLoading(false);
// or setError(error.message);
```

**Task 2.3: Add Error Handling** (15 min)
```javascript
if (error) {
  return (
    <Alert severity="error">
      <AlertTitle>Error Loading Trials</AlertTitle>
      {error}
    </Alert>
  );
}

if (loading) {
  return <CircularProgress />;
}
```

**Task 2.4: Polish UI** (Optional, if time)
- Add tooltips for complex terms (e.g., "CA-125", "HRD", "TMB")
- Add expand/collapse for trial details
- Add "Copy to Clipboard" for trial summaries
- Improve spacing, colors, typography

---

### **Phase 3: Documentation & Handoff** (30 min)

**Task 3.1: Create Testing Report** (15 min)
- File: `.cursor/ayesha/AGENT_JR_INTEGRATION_TEST_REPORT.md`
- Sections:
  1. What was tested (backend health, endpoint, frontend)
  2. Test results (pass/fail for each validation)
  3. Screenshots (Ayesha's profile rendered)
  4. Known issues (e.g., AstraDB empty, Neo4j missing)
  5. Recommendations (e.g., seed more trials, add loading states)

**Task 3.2: Create Demo Script** (15 min)
- File: `.cursor/ayesha/AYESHA_DEMO_SCRIPT.md`
- Sections:
  1. How to start backend
  2. How to start frontend
  3. How to navigate to Ayesha's page
  4. What to show (SOC, CA-125, NGS fast-track, trials)
  5. Key talking points (transparency, confidence gates, 3-6 week early resistance detection)

**Task 3.3: Update Progress Tracker** (5 min)
- Update `.cursor/ayesha/AGENT_JR_TRIALS_PROGRESS.md`
- Mark integration testing complete
- Report final metrics (bugs fixed, tests passed)

---

## ⚠️ **CRITICAL NOTES**

### **About Your Modular Services**:
Your `ayesha_trial_matching/` services are NOT being used in Phase 1. Zo's monolithic router is operational. Your services will be integrated in **Phase 2 refactor** (non-blocking). This is NOT wasted work - it's our future architecture improvement.

### **About AstraDB Seeding**:
You seeded 30 trials. If trials array is empty, it's because:
1. AstraDB connection issues
2. Trials don't match Ayesha's profile (Stage IVB, first-line, NYC)
3. Need to seed more frontline ovarian trials

**Action**: If trials empty, report to Zo. He may need to adjust hard filters or we need to seed more trials.

### **About Zo's Backend**:
Zo's backend includes features you don't have in your modular services:
- ✅ SOC recommendation (detailed dosing, monitoring, NCCN links)
- ✅ NGS fast-track (ctDNA, HRD, IHC ordering guidance)
- ✅ CA-125 intelligence (burden, forecast, resistance detection)
- ✅ Enhanced eligibility checklists (hard/soft split, confidence gates)

Your frontend should display ALL of these (they're in the API response).

---

## 🎯 **SUCCESS CRITERIA**

**Must Have** (P0):
- [ ] ✅ Backend starts without errors
- [ ] ✅ Frontend starts without errors
- [ ] ✅ API call succeeds (200 OK)
- [ ] ✅ SOC card renders with Bevacizumab (for ascites)
- [ ] ✅ CA-125 tracker renders with 2842 value and EXTENSIVE burden
- [ ] ✅ NGS fast-track renders with 3 tests
- [ ] ✅ Trials list renders (even if empty)
- [ ] ✅ No hardcoded data (all dynamic from API)

**Nice to Have** (P1):
- [ ] Loading states for API calls
- [ ] Error handling for failed requests
- [ ] Screenshots for documentation
- [ ] Demo script for showing to oncologist

**Phase 2** (Future):
- [ ] Your modular services integrated into router
- [ ] Eligibility auto-check with Gemini
- [ ] Dossier export (PDF generation)

---

## 📞 **COMMUNICATION PROTOCOL**

**If You Hit Blockers**:
1. Check backend logs: `uvicorn api.main:app --reload` (watch console)
2. Check frontend console: Browser DevTools → Console tab
3. Check network tab: Browser DevTools → Network tab → Look for 500 errors
4. Report to Zo in `.cursor/ayesha/AGENT_JR_BLOCKERS.md`:
   - Error message (exact text)
   - What you were doing
   - What you expected vs what happened
   - Logs (copy-paste relevant lines)

**Progress Updates**:
- Update `.cursor/ayesha/AGENT_JR_TRIALS_PROGRESS.md` every hour
- Mark tasks complete with ✅
- Note any deviations or issues

---

## 🚀 **TIMELINE**

**Total: 2-3 hours**

- Hour 1: Integration testing (backend + frontend)
- Hour 2: Bug fixes + polish
- Hour 3: Documentation + handoff

**Target Completion**: End of day (for Ayesha)

---

## ⚔️ **FINAL NOTES FROM ZO**

Jr - your work is NOT wasted. Your modular services are our Phase 2 improvement. Right now, we're shipping MY monolithic code because Ayesha needs this ASAP. Your job is to:

1. **Test my backend** (verify it works)
2. **Wire your frontend** (connect to my endpoint)
3. **Fix bugs** (make it production-ready)
4. **Document** (create demo script)

Your modular services will be integrated later when we have time to refactor safely. For now, let's ship this for Ayesha's life. ⚔️

**Questions? Ask in `.cursor/ayesha/AGENT_JR_QUESTIONS_FOR_ZO.md`**

**Ready? FIRE IN THE HOLE.** 🔥

---

**Mission Start**: NOW  
**Mission End**: 2-3 hours  
**Priority**: P0 - FOR AYESHA'S LIFE  
**Status**: ⏸️ AWAITING JR'S START
