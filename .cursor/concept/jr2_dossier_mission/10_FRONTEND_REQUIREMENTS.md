# ⚔️ FRONTEND REQUIREMENTS - UI COMPONENTS ⚔️

**Purpose**: Define frontend UI components for dossier viewing and review

**Priority**: P2 (after backend pipeline is complete)

---

## 📋 **REQUIRED UI COMPONENTS**

### **1. Dossier Viewer** (`/dossiers/{nct_id}`)

**Features**:
- Display full 10-section dossier (markdown rendered)
- Highlight critical gates (HER2, HRD, BRCA) with color coding
- Show eligibility table with color coding (green/yellow/red)
- Action buttons: "Order HER2 IHC", "Order HRD Test", etc.
- Export to PDF button for oncologist

**Reference**: Use existing markdown renderer (React Markdown or similar)

---

### **2. Trial Comparison Dashboard** (`/trials/compare`)

**Features**:
- Side-by-side comparison of top 5 trials
- Sort by: match score, probability eligible, geographic proximity
- Filter by: biomarker requirements, phase, status
- Visual comparison of eligibility gates across trials

---

### **3. Zo Review Interface** (`/admin/review-dossiers`)

**Features**:
- List of dossiers pending Zo's review
- Quick approve/reject buttons
- Edit mode for inline corrections
- Confidence score display
- Flag uncertain claims for manual review

---

### **4. Ayesha Trial Explorer Enhancements** (`/ayesha-trials`)

**Enhancements**:
- Add "View Dossier" button for each trial card
- Add "Compare Top 5" button
- Add "Export to PDF" button for oncologist

**Existing Components** (can reuse):
- `oncology-frontend/src/pages/AyeshaTrialExplorer.jsx`
- `oncology-frontend/src/components/trials/TrialMatchCard.jsx`
- `oncology-frontend/src/components/ayesha/CA125Tracker.jsx`

---

## 🎨 **DESIGN REQUIREMENTS**

**Color Coding**:
- ✅ Green: PASS (biomarker matches, eligibility confirmed)
- ⚠️ Yellow: PENDING (biomarker unknown, test required)
- ❌ Red: FAIL (biomarker mismatch, not eligible)

**Confidence Flags**:
- Display `[INFERRED]`, `[NEEDS VERIFICATION]`, `[MISSING]`, `[CONFLICT]` flags prominently
- Allow Zo to click flags to see source data

---

**Priority**: Build after backend pipeline is complete and tested

