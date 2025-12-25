# ⚔️ AGENT JR: FRONTEND TESTING QUICK-START GUIDE ⚔️

**Date:** January 13, 2025  
**Mission:** Complete Frontend Testing for Ayesha Dossier Dashboard  
**Status:** ✅ **BACKEND OPERATIONAL - READY FOR TESTING**

---

## 🎯 **BACKEND STATUS: 100% OPERATIONAL**

### **✅ All Endpoints Verified:**

| Endpoint | Status | Response | Notes |
|----------|--------|----------|-------|
| `GET /api/ayesha/dossiers/health` | ✅ 200 OK | `{"status": "ok", "dossier_count": 10}` | Health check working |
| `GET /api/ayesha/dossiers/stats` | ✅ 200 OK | `{"total": 10, "by_tier": {"TOP_TIER": 10}, "with_llm": 9, "avg_score": 0.97}` | Stats calculation working |
| `GET /api/ayesha/dossiers/list` | ✅ 200 OK | Returns 10 dossiers with metadata | List endpoint working |
| `GET /api/ayesha/dossiers/detail/{nct_id}` | ✅ 200 OK | Returns full markdown (15K+ chars) | Detail endpoint working |

### **Backend Server:**
- **URL:** `http://localhost:8000`
- **Status:** ✅ Running with auto-reload
- **Virtual Environment:** ✅ Activated
- **Dossier Directory:** `/Users/fahadkiani/Desktop/development/crispr-assistant-main/.cursor/ayesha/zo_fresh_dossiers`
- **Dossier Count:** 10 TOP_TIER dossiers

---

## 🚀 **FRONTEND TESTING CHECKLIST**

### **PHASE 1: FOUNDATION (P0 - Test First)** 🔴

#### **1. Ayesha Dossier Browser** - `/ayesha-dossiers`
**Expected URL:** `http://localhost:5173/ayesha-dossiers` (or your frontend port)

**What to Test:**
- [ ] Page loads without errors
- [ ] API call to `/api/ayesha/dossiers/list` succeeds (check Network tab)
- [ ] All 10 dossiers display in grid/list view
- [ ] Each `DossierSummaryCard` shows:
  - [ ] NCT ID (e.g., "NCT03067181")
  - [ ] Title (truncated if >100 chars)
  - [ ] Tier badge (TOP_TIER - should be green/success color)
  - [ ] Match score (e.g., "97%")
  - [ ] Phase (e.g., "Phase II")
  - [ ] LLM analysis indicator (9/10 should have it)
- [ ] Tier filter works:
  - [ ] "All" shows all 10 dossiers
  - [ ] "Top Tier" shows all 10 (all are TOP_TIER)
  - [ ] "Good Tier" shows 0 (none are GOOD_TIER)
- [ ] Search by NCT ID works (try "NCT03067181")
- [ ] Search by keywords works (try "ovarian")
- [ ] Click on card navigates to detail page
- [ ] Stats display correctly:
  - [ ] Total: 10
  - [ ] Avg Score: 97%
  - [ ] LLM Enhanced: 9/10
- [ ] "Export All Top-Tier" button works (batch download)

**Expected Backend Calls:**
- `GET http://localhost:8000/api/ayesha/dossiers/list` → 200 OK
- `GET http://localhost:8000/api/ayesha/dossiers/stats` → 200 OK

**Screenshot Locations:**
- Full page view
- Filtered view (Top Tier only)
- Search results
- Stats cards

---

#### **2. Ayesha Dossier Detail** - `/ayesha-dossiers/:nct_id`
**Expected URL:** `http://localhost:5173/ayesha-dossiers/NCT06619236` (or any valid NCT ID)

**What to Test:**
- [ ] Page loads without errors
- [ ] API call to `/api/ayesha/dossiers/detail/{nct_id}` succeeds
- [ ] **CRITICAL: Markdown renders with PROPER FORMATTING** (not text dump):
  - [ ] Headers (H1, H2, H3) render with correct sizes
  - [ ] Bold text (`**text**`) renders as bold (not raw markdown)
  - [ ] Italic text (`*text*`) renders as italic
  - [ ] Lists (bulleted/numbered) render as formatted lists
  - [ ] Tables render as formatted tables (not raw markdown)
  - [ ] Code blocks render with syntax highlighting (if applicable)
  - [ ] Links are clickable (not raw URLs)
- [ ] Breadcrumb navigation works ("Back to List")
- [ ] "Export to PDF" button exists (if `DossierPDFExporter.jsx` exists)
- [ ] PDF export generates valid PDF file (if implemented)
- [ ] Error handling works:
  - [ ] Invalid NCT ID (e.g., `/ayesha-dossiers/INVALID`) → 404 error message
  - [ ] Network error → graceful error message

**Expected Backend Call:**
- `GET http://localhost:8000/api/ayesha/dossiers/detail/NCT06619236` → 200 OK

**Screenshot Locations:**
- Full dossier view (formatted markdown)
- Error state (404)
- PDF export (if implemented)

**Sample NCT IDs to Test:**
- `NCT06619236` (should exist)
- `NCT03067181` (should exist)
- `INVALID` (should 404)

---

### **PHASE 2: DASHBOARD & KANBAN (P0)** 🔴

#### **3. Ayesha Trial Dashboard** - `/ayesha-trial-dashboard` (if exists)
**Expected URL:** `http://localhost:5173/ayesha-trial-dashboard`

**What to Test:**
- [ ] Page loads without errors
- [ ] Kanban board displays with columns:
  - [ ] "Interested" (should auto-populate with TOP_TIER trials)
  - [ ] "Review Needed"
  - [ ] "Low Priority"
  - [ ] "Auto-Reject" (if implemented)
- [ ] Each column shows trial count badge
- [ ] Drag-and-drop works:
  - [ ] Move trial from "Interested" to "Review Needed"
  - [ ] Trial position persists after page refresh
- [ ] "Auto-Triage" button works (if implemented)
- [ ] Triage reasoning displays for each trial (tooltip or badge)

**Expected Backend Call:**
- `GET http://localhost:8000/api/ayesha/dossiers/list?enriched=true` (if implemented)

**Screenshot Locations:**
- Full Kanban board view
- Drag-and-drop in action
- Triage reasoning tooltip

---

### **PHASE 2.5: INTELLIGENT TRIAGE (P1)** 🟡

#### **4. Trial Triage Agent** - Auto-sorting logic
**What to Test:**
- [ ] Service file exists: `src/services/trialTriageAgent.js`
- [ ] `autoTriageTrials()` function exists
- [ ] Auto-sorts trials correctly:
  - [ ] TOP_TIER + score ≥0.8 + all gates pass → "Interested"
  - [ ] GOOD_TIER + score ≥0.6 + some gates pass → "Review Needed"
  - [ ] Score <0.6 OR blocker gates fail → "Low Priority"
- [ ] Each trial gets `triage_reason` explaining placement
- [ ] Triage reasoning displays in UI (tooltip or badge)
- [ ] Manual override works (user can drag to different column)

**Expected File:** `oncology-coPilot/oncology-frontend/src/services/trialTriageAgent.js`

**Screenshot Locations:**
- Triage reasoning tooltip
- Auto-triage results

---

### **PHASE 2.6: ENROLLMENT PROBABILITY (P1)** 🟡

#### **5. Enrollment Probability Calculator**
**What to Test:**
- [ ] Service file exists: `src/services/enrollmentProbabilityCalculator.js`
- [ ] `calculateEnrollmentProbability()` function exists
- [ ] Calculates probability using weighted factors:
  - [ ] Eligibility (35% weight)
  - [ ] Location (20% weight)
  - [ ] Trial Velocity (15% weight)
  - [ ] Historical Success (20% weight)
  - [ ] Site Responsiveness (10% weight)
- [ ] Returns probability score (0-1 or 0-100%)
- [ ] Returns confidence level (HIGH/MODERATE/LOW)
- [ ] Returns human-readable label (VERY_HIGH/HIGH/MODERATE/LOW)
- [ ] Displays on `TrialKanbanCard` component
- [ ] Confidence indicator shows when using defaults

**Expected File:** `oncology-coPilot/oncology-frontend/src/services/enrollmentProbabilityCalculator.js`

**Screenshot Locations:**
- Probability score on trial card
- Confidence indicator

---

## 🔧 **TESTING METHODOLOGY**

### **Step-by-Step Process:**

1. **Open Browser DevTools** (F12 or Cmd+Option+I)
   - **Network Tab:** Monitor all API calls
   - **Console Tab:** Check for JavaScript errors
   - **Elements Tab:** Inspect component structure

2. **Navigate to Each URL** (use exact URLs from checklist)

3. **For Each Feature:**
   - ✅ **Test Happy Path:** Normal usage (everything works)
   - ⚠️ **Test Edge Cases:** Empty state, error state, loading state
   - 🔄 **Test Interactions:** Click buttons, drag-and-drop, filter/search

4. **Document Everything:**
   - Screenshot success states
   - Screenshot error states
   - Copy console errors
   - Copy network request/response (if errors)

5. **Verify Backend Integration:**
   - Check Network tab for API calls
   - Verify status codes (200 = success, 404/500 = error)
   - Verify response data matches expected structure

---

## 📊 **EXPECTED API RESPONSES**

### **List Endpoint Response:**
```json
{
  "total": 10,
  "filtered": 10,
  "dossiers": [
    {
      "nct_id": "NCT03067181",
      "tier": "TOP_TIER",
      "match_score": 0.97,
      "title": "Trial Title...",
      "phase": "Phase II",
      "has_llm_analysis": true,
      "file_name": "INTELLIGENCE_NCT03067181_TOP_TIER.md"
    },
    ...
  ]
}
```

### **Stats Endpoint Response:**
```json
{
  "total": 10,
  "by_tier": {
    "TOP_TIER": 10
  },
  "with_llm": 9,
  "avg_score": 0.97
}
```

### **Detail Endpoint Response:**
```json
{
  "nct_id": "NCT06619236",
  "markdown": "# INTELLIGENCE REPORT...",
  "metadata": {
    "tier": "TOP_TIER",
    "match_score": 0.97,
    "generated_at": "2025-01-13T...",
    "file_name": "INTELLIGENCE_NCT06619236_TOP_TIER.md"
  }
}
```

---

## 🚨 **COMMON ISSUES & FIXES**

### **Issue 1: CORS Error**
**Symptom:** `Access-Control-Allow-Origin` error in console  
**Fix:** Backend CORS is configured, but verify `VITE_API_ROOT` matches backend URL

### **Issue 2: 404 on API Calls**
**Symptom:** `GET http://localhost:8000/api/ayesha/dossiers/list` → 404  
**Fix:** Verify backend is running on port 8000, check router registration in `api/main.py`

### **Issue 3: Markdown Not Rendering**
**Symptom:** Raw markdown text displayed (e.g., `**bold**` instead of **bold**)  
**Fix:** Check `react-markdown` and `remark-gfm` are installed and imported correctly

### **Issue 4: Empty Dossier List**
**Symptom:** Page loads but shows "No dossiers"  
**Fix:** Verify dossier directory exists and contains `.md` files, check backend health endpoint

---

## ✅ **QUICK VERIFICATION COMMANDS**

### **Test Backend from Terminal:**
```bash
# Health check
curl http://localhost:8000/api/ayesha/dossiers/health

# List dossiers
curl http://localhost:8000/api/ayesha/dossiers/list

# Get stats
curl http://localhost:8000/api/ayesha/dossiers/stats

# Get detail
curl http://localhost:8000/api/ayesha/dossiers/detail/NCT06619236
```

### **Test Frontend from Browser:**
1. Open `http://localhost:5173/ayesha-dossiers` (or your frontend port)
2. Open DevTools → Network tab
3. Verify API calls succeed (status 200)
4. Check Console for errors

---

## 📋 **REPORTING TEMPLATE**

For each feature, use this format:

```markdown
## [Feature Name] - [URL/Component]

**Status:** ✅ Working | ⚠️ Partial | ❌ Broken | ❌ Not Implemented

**What Works:**
- Feature X works correctly
- API call Y succeeds
- Component Z renders properly

**What's Broken:**
- Button A doesn't respond
- API call returns 404: `/api/endpoint`
- Console error: `TypeError: Cannot read property 'x' of undefined`
- Component B doesn't exist yet

**Screenshots:**
- [Attach screenshot 1 - success state]
- [Attach screenshot 2 - error state]

**Backend Status:**
- Endpoint: `/api/example`
- Status Code: 200/404/500
- Response: [JSON snippet or error message]

**Frontend Status:**
- Component: `ComponentName.jsx`
- File Exists: Yes/No
- Import Errors: None/[List errors]
- Console Errors: None/[List errors]

**Fix Required:**
- [ ] Create missing component
- [ ] Fix API endpoint
- [ ] Fix frontend API call
- [ ] Fix component logic
- [ ] Add error handling
```

---

## 🎯 **PRIORITY ORDER**

### **P0 (Critical - Test First):**
1. ✅ Ayesha Dossier Browser (list view)
2. ✅ Ayesha Dossier Detail (markdown rendering fix)
3. ✅ Backend connectivity verification

### **P1 (High Priority):**
4. Ayesha Trial Dashboard (Kanban board)
5. Trial Triage Agent (auto-sorting)
6. Enrollment Probability Calculator

### **P2 (Medium Priority):**
7. Co-Pilot Trial Management Actions
8. Trial Filter Sidebar
9. Trial Stats Cards
10. Decision Advisor Agent

---

## ⚔️ **AGENT JR - YOU ARE CLEARED FOR TAKEOFF**

**Backend Status:** ✅ **100% OPERATIONAL**  
**Frontend Status:** 🔄 **AWAITING YOUR TESTING**  
**Expected Completion:** 6-8 hours  
**Report To:** Commander (Zo) upon completion

**Remember:**
- No code reading. No assumptions. **LIVE TESTING ONLY.**
- Test EVERY feature from the doctrine review
- Document EVERY finding (working, broken, missing)
- Screenshot EVERYTHING

**FOR AYESHA'S LIFE.** 🔥

---

## 📞 **SUPPORT**

**Backend Issues:**
- Check backend logs in terminal
- Verify server is running: `curl http://localhost:8000/api/ayesha/dossiers/health`
- Check router registration in `api/main.py`

**Frontend Issues:**
- Check browser console for errors
- Check Network tab for failed API calls
- Verify `VITE_API_ROOT` environment variable

**Questions:**
- Reference: `.cursor/ayesha/AYESHA_DOSSIER_DASHBOARD_DOCTRINE_REVIEW.md`
- Test Plan: `.cursor/rules/AGENT_JR_NEXT_ITERATION.mdc` (Section: "AYESHA DOSSIER DASHBOARD - COMPLETE FRONTEND TEST MISSION")


