# AGENT 4: STEP-BY-STEP IMPLEMENTATION GUIDE

## **Build Order & Dependencies**

---

## **Phase 1: Foundation (1.5 hours)**

### **Step 1: Create Folder Structure (5 min)**
```bash
cd oncology-coPilot/oncology-frontend/src
mkdir -p components/research utils hooks
```

### **Step 2: Module 7 - Refresh Hook (15 min)**
**Why First:** Foundation for Module 2 (Refresh Button)

**File:** `hooks/useTrialRefresh.js`
- Create hook with API call logic
- Test with mock endpoint

**Dependencies:** None

### **Step 3: Module 3 - Location Card (30 min)**
**Why Next:** Standalone component, no dependencies

**File:** `components/research/LocationCard.jsx`
- Build location display component
- Test with mock location data

**Dependencies:** None

### **Step 4: Module 1 - Trial Filters (45 min)**
**Why Next:** Standalone component, needed for integration

**File:** `components/research/TrialFilters.jsx`
- Build filter component
- Test filter state changes

**Dependencies:** None

---

## **Phase 2: Integration Components (1 hour)**

### **Step 5: Module 2 - Refresh Button (30 min)**
**Why Now:** Uses Module 7 hook

**File:** `components/research/RefreshStatusButton.jsx`
- Import and use `useTrialRefresh` hook
- Build button UI
- Test refresh functionality

**Dependencies:** Module 7

### **Step 6: Module 4 - Enhanced ResultsDisplay (30 min)**
**Why Now:** Uses Module 3 LocationCard

**File:** `components/research/ResultsDisplay.jsx` (modify)
- Import LocationCard
- Add location display section
- Add live refresh indicator
- Test with existing trial data

**Dependencies:** Module 3

---

## **Phase 3: Main Integration (1 hour)**

### **Step 7: Module 6 - PDF Export (30 min)**
**Why Now:** Standalone utility, needed for integration

**File:** `utils/exportTrialsPDF.js`
- Build HTML template
- Test PDF generation

**Dependencies:** None

### **Step 8: Module 5 - Research.jsx Integration (30 min)**
**Why Last:** Wires all components together

**File:** `pages/Research.jsx` (modify)
- Import all new components
- Add filter state and logic
- Add refresh handler
- Add PDF export button
- Test complete workflow

**Dependencies:** All modules (1, 2, 3, 4, 6, 7)

---

## **Phase 4: Testing (30 min)**

### **Step 9: Module 8 - Test Suite (30 min)**
**File:** `tests/agent_4_frontend/`

**Tests:**
- `test_trial_filters.test.jsx` - Filter component tests
- `test_refresh_status.test.jsx` - Refresh button tests
- `test_location_display.test.jsx` - Location card tests
- `test_pdf_export.test.jsx` - PDF export tests

**Run:**
```bash
npm test -- tests/agent_4_frontend/
```

---

## **Execution Checklist**

### **Pre-Flight:**
- [ ] Review existing `Research.jsx` and `ResultsDisplay.jsx`
- [ ] Verify backend endpoints available (`/api/trials/refresh_status`)
- [ ] Check database has `locations_data` field populated
- [ ] Frontend dev server running

### **Build Order:**
1. [ ] Module 7: Refresh Hook (15 min)
2. [ ] Module 3: Location Card (30 min)
3. [ ] Module 1: Trial Filters (45 min)
4. [ ] Module 2: Refresh Button (30 min)
5. [ ] Module 4: Enhanced Results (30 min)
6. [ ] Module 6: PDF Export (30 min)
7. [ ] Module 5: Integration (30 min)
8. [ ] Module 8: Tests (30 min)

**Total: ~3 hours**

---

## **Testing Strategy**

### **Unit Tests (Per Component):**
- Component renders
- Props handling
- Event handlers
- Edge cases

### **Integration Tests:**
- Filters + Search workflow
- Refresh + Filter workflow
- Location display in results

### **E2E Tests:**
- Complete user flow: Search → Filter → Refresh → Export

---

## **Rollback Plan**

If issues occur:
1. Revert `Research.jsx` changes
2. Keep new components (don't break existing search)
3. Re-enable features incrementally

---

**READY TO BUILD** 🚀

