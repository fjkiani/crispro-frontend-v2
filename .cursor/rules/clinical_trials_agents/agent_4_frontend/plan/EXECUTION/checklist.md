# AGENT 4: EXECUTION CHECKLIST

## **Pre-Flight Checks**

### **Backend Status:**
- [ ] Agent 1 complete (1000 trials in database with `locations_data`)
- [ ] Agent 2 endpoint `/api/trials/refresh_status` operational
- [ ] Test endpoint: `curl -X POST http://localhost:8000/api/trials/refresh_status -H "Content-Type: application/json" -d '{"nct_ids": ["NCT12345"]}'`

### **Frontend Environment:**
- [ ] Frontend dev server running (`npm start`)
- [ ] MUI components available
- [ ] Existing `Research.jsx` page accessible
- [ ] No console errors on page load

### **Database Verification:**
```bash
# Check trials with locations_data
sqlite3 backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials WHERE locations_data IS NOT NULL;"
```

---

## **Implementation Checklist**

### **Phase 1: Foundation**
- [ ] Module 7: `hooks/useTrialRefresh.js` created
- [ ] Module 7: Hook tested with mock endpoint
- [ ] Module 3: `components/research/LocationCard.jsx` created
- [ ] Module 3: LocationCard renders with test data
- [ ] Module 1: `components/research/TrialFilters.jsx` created
- [ ] Module 1: Filters render and handle state changes

### **Phase 2: Integration Components**
- [ ] Module 2: `components/research/RefreshStatusButton.jsx` created
- [ ] Module 2: Refresh button calls hook correctly
- [ ] Module 4: `ResultsDisplay.jsx` modified with LocationCard
- [ ] Module 4: Locations display in trial cards

### **Phase 3: Main Integration**
- [ ] Module 6: `utils/exportTrialsPDF.js` created
- [ ] Module 6: PDF export opens print dialog
- [ ] Module 5: `Research.jsx` modified with all components
- [ ] Module 5: Filters apply to search results
- [ ] Module 5: Refresh button updates results
- [ ] Module 5: PDF export works

### **Phase 4: Testing**
- [ ] Module 8: Test files created
- [ ] All unit tests pass
- [ ] Integration tests pass
- [ ] E2E tests pass (if applicable)

---

## **Verification Steps**

### **1. Filter Functionality:**
```bash
# Manual test:
1. Navigate to /research
2. Search for "ovarian cancer"
3. Select "Gynecologic Oncology" filter
4. Verify only matching trials show
5. Select "NY" state filter
6. Verify only NY trials show
```

### **2. Refresh Status:**
```bash
# Manual test:
1. Search for trials
2. Click "Refresh Live Status" button
3. Verify button shows loading state
4. Verify trial status updates
5. Verify "Live status refreshed" badge appears
```

### **3. Location Display:**
```bash
# Manual test:
1. Search for trials
2. Verify locations display below trial title
3. Verify contact info shows (phone/email)
4. Verify status badge shows (RECRUITING/NOT_YET_RECRUITING)
```

### **4. PDF Export:**
```bash
# Manual test:
1. Search for trials
2. Click "Export PDF" button
3. Verify print dialog opens
4. Verify PDF contains trial summaries
5. Verify locations included in PDF
```

---

## **Acceptance Criteria Verification**

### **Must Have:**
- [ ] Filters work (Disease/Phase/State)
- [ ] Refresh button updates trial status
- [ ] Locations display with contact info
- [ ] PDF export generates summary
- [ ] All components integrate into Research.jsx
- [ ] No breaking changes to existing search

### **Nice to Have:**
- [ ] Filters persist in URL params
- [ ] Loading states for all async operations
- [ ] Error messages user-friendly
- [ ] Responsive design (mobile-friendly)

---

## **Post-Implementation**

### **Documentation:**
- [ ] Update component docs
- [ ] Add usage examples
- [ ] Document API endpoint requirements

### **Cleanup:**
- [ ] Remove console.logs
- [ ] Remove unused imports
- [ ] Verify no lint errors

---

## **Rollback Plan**

If critical issues:
1. Git commit current state
2. Revert `Research.jsx` to previous version
3. Keep new components (they're independent)
4. Test existing search still works
5. Re-enable features one by one

---

**STATUS: READY FOR EXECUTION** ✅

