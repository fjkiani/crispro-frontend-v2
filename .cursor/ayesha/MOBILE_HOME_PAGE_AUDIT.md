# 📱 MOBILE HOME PAGE AUDIT - AK'S FIRST LOAD EXPERIENCE

**Date**: January 29, 2025  
**Issue**: AK sees "book-like" legacy dashboard on `/home` instead of patient-relevant content  
**Priority**: **P1 - CRITICAL** - First impression is broken

---

## 🚨 PROBLEM IDENTIFIED

### **What AK Sees on `/home`:**

**Current Route**: `/home` → `Home.jsx` → `DisplayInfo.jsx`

**DisplayInfo Component Shows:**
- ❌ **Legacy metrics dashboard** (Specialist Appointments, Treatment Progress, Folders, Screenings)
- ❌ **Not patient-relevant** - These are generic admin/doctor metrics
- ❌ **Looks like a "book"** - Static grid of metric cards with no actionable content
- ❌ **Not mobile-optimized** - Uses Tailwind grid classes that may not work well on mobile

**Root Cause:**
1. `/home` route (App.jsx line 117) renders `<Home />` which only renders `<DisplayInfo />`
2. `/` route (App.jsx line 111) also renders `<Home />` 
3. **No patient-specific redirect** - Patients land on generic dashboard instead of `/ayesha-trials`
4. **No PatientRoute wrapper** on `/home` or `/` routes
5. **No mobile-first design** - DisplayInfo uses desktop-first grid layouts

---

## 📊 CURRENT ROUTING STATE

### **App.jsx Routes (Lines 101-194):**

```jsx
{/* Root route shows the home page */}
<Route path="/" element={<Home />} />
<Route path="/home" element={<Home />} />

{/* Patient routes - NO REDIRECT FROM /home */}
<Route path="/ayesha-trials" element={<AyeshaTrialExplorer />} />
<Route path="/ayesha-complete-care" element={<AyeshaCompleteCare />} />
```

**Issues:**
- ❌ `/` and `/home` both render legacy `Home` component
- ❌ No patient role check or redirect logic
- ❌ Patients land on wrong page after login

### **DisplayInfo.jsx Component:**

```jsx
const DisplayInfo = () => {
  // Shows generic metrics cards:
  // - Specialist Appointments Pending
  // - Treatment Progress Update  
  // - Total Folders
  // - Total Screenings
  // - Completed/Pending/Overdue Screenings
  
  return (
    <div className="flex flex-wrap gap-[26px]">
      <div className="mt-7 grid w-full gap-4 sm:grid-cols-2 sm:gap-6 lg:grid-cols-2">
        {metricsData.slice(0, 2).map((metric) => (
          <MetricsCard key={metric.title} {...metric} />
        ))}
      </div>
      <div className="mt-[9px] grid w-full gap-4 sm:grid-cols-2 sm:gap-6 lg:grid-cols-4">
        {metricsData.slice(2).map((metric) => (
          <MetricsCard key={metric.title} {...metric} />
        ))}
      </div>
    </div>
  );
};
```

**Mobile Issues:**
- ❌ Uses `sm:grid-cols-2 lg:grid-cols-4` - Desktop-first breakpoints
- ❌ Large gaps (`gap-[26px]`) take up mobile screen space
- ❌ Cards not optimized for mobile touch targets
- ❌ No mobile navigation or quick actions

---

## ✅ WHAT AK SHOULD SEE (Per CLINICAL_MASTER_FRONTEND_AUDIT.md)

### **Expected Patient Landing Page:**

**After login, AK should be redirected to `/ayesha-trials`** which shows:

1. ✅ **Profile Summary** - Stage IVB, CA-125 2,842, germline-positive (MBD4), awaiting NGS
2. ✅ **SOC Recommendation** - Carboplatin + Paclitaxel + Bevacizumab
3. ✅ **CA-125 Tracker** - Burden classification, response forecast
4. ✅ **Resistance Alert Banner** - If resistance detected
5. ✅ **Next Test Recommender** - HRD test, ctDNA panel priorities
6. ✅ **Hint Tiles** - Max 4 actionable hints
7. ✅ **Mechanism Map** - 6 pathway chips (DDR, MAPK, PI3K, VEGF, HER2, IO)
8. ✅ **Top 10 Clinical Trials** - **WITH HOLISTIC SCORES** (just integrated)
9. ✅ **Resistance Playbook** - Alternative strategies
10. ✅ **SAE Features** - DNA repair capacity, pathway burden

**Mobile-Optimized:**
- ✅ Stack components vertically (single column on mobile)
- ✅ Collapsible sections for dense information
- ✅ Touch-friendly buttons and cards
- ✅ Clear navigation to other pages (Complete Care, Profile)

---

## 🎯 FIXES REQUIRED

### **Fix 1: Redirect Patients from `/home` and `/` to `/ayesha-trials`**

**Option A: Smart Redirect in Home Component**

```jsx
// Home.jsx
import { useNavigate } from 'react-router-dom';
import { useAuth } from '../context/AuthContext';
import { usePatient } from '../context/PatientContext';

const Home = () => {
  const navigate = useNavigate();
  const { profile } = useAuth();
  const { hasProfile } = usePatient();

  useEffect(() => {
    // Redirect patients to their dashboard
    if (profile?.role === 'patient') {
      if (hasProfile) {
        navigate('/ayesha-trials', { replace: true });
      } else {
        navigate('/patient/onboarding', { replace: true });
      }
      return;
    }
    
    // Doctors/other roles see DisplayInfo
    // Keep existing DisplayInfo logic
  }, [profile, hasProfile, navigate]);

  // Show loading while redirecting
  return <CircularProgress />;
};
```

**Option B: Protected Route Wrapper**

```jsx
// App.jsx - Replace existing /home and / routes
<Route 
  path="/" 
  element={
    <PatientRoute>
      <Navigate to="/ayesha-trials" replace />
    </PatientRoute>
  } 
/>
<Route 
  path="/home" 
  element={
    <PatientRoute>
      <Navigate to="/ayesha-trials" replace />
    </PatientRoute>
  } 
/>
```

**Recommendation**: **Option A** - Keeps legacy DisplayInfo for non-patient roles, redirects patients automatically

---

### **Fix 2: Mobile-First AyeshaTrialExplorer**

**Current Issues:**
- Components may not be mobile-optimized
- Need to audit responsive breakpoints
- Need collapsible sections for dense information

**Action Items:**
1. [ ] Audit `AyeshaTrialExplorer.jsx` mobile responsiveness
2. [ ] Verify all cards stack vertically on mobile
3. [ ] Add mobile navigation (bottom nav or hamburger menu)
4. [ ] Test touch targets (buttons, cards) are at least 44x44px

---

### **Fix 3: Remove or Hide Legacy Routes for Patients**

**Per CLINICAL_MASTER_FRONTEND_AUDIT.md:**
> "MOAT navigation" with only relevant pages, not all legacy/experimental ones

**Action Items:**
1. [ ] Hide `/home` route from patient navigation
2. [ ] Ensure sidebar/navigation only shows patient-relevant pages:
   - `/ayesha-trials` - Trial Explorer
   - `/ayesha-complete-care` - Complete Care Plan
   - `/patient/profile` - Patient Profile
3. [ ] Remove or gate legacy routes (agent-dashboard, genomic-analysis, etc.) from patient access

---

## 📋 COMPONENT AUDIT (Per CLINICAL_MASTER_FRONTEND_AUDIT.md)

### **Expected Components on `/ayesha-trials`:**

| Component | Status | Mobile-Ready | Notes |
|-----------|--------|--------------|-------|
| **SOCRecommendationCard** | ✅ EXISTS | ⚠️ Needs audit | Should stack vertically on mobile |
| **CA125Tracker** | ✅ EXISTS | ⚠️ Needs audit | May need collapsible sections |
| **ResistanceAlertBanner** | ✅ EXISTS | ⚠️ Needs audit | Should be prominent but not overwhelming |
| **NextTestCard** | ✅ EXISTS | ⚠️ Needs audit | Cards should be touch-friendly |
| **HintTilesPanel** | ✅ EXISTS | ⚠️ Needs audit | 4 tiles should wrap on mobile |
| **MechanismChips** | ✅ EXISTS | ⚠️ Needs audit | 6 chips should wrap to multiple rows |
| **TrialMatchCard** | ✅ EXISTS | ⚠️ Needs audit | **NOW WITH HOLISTIC SCORES** - verify display |
| **ResistancePlaybook** | ✅ EXISTS | ⚠️ Needs audit | Accordion should work on mobile |
| **AyeshaSAEFeaturesCard** | ✅ EXISTS | ⚠️ Needs audit | Dense data - needs mobile optimization |
| **HolisticScoreCard** | ✅ **NEW** | ⚠️ **NEEDS MOBILE AUDIT** | Just created - verify mobile layout |

---

## 🚨 CRITICAL ISSUES

### **Issue 1: Wrong Landing Page** 🔴 **CRITICAL**

**Problem**: Patients land on `/home` which shows generic metrics dashboard (DisplayInfo)

**Impact**: 
- ❌ Poor first impression ("looks like a book")
- ❌ No actionable content for patients
- ❌ Confusing navigation

**Fix**: Redirect patients to `/ayesha-trials` automatically

---

### **Issue 2: Missing PatientRoute Wrapper** 🔴 **CRITICAL**

**Problem**: `/home` and `/` routes not wrapped with PatientRoute

**Impact**:
- ❌ No role-based access control
- ❌ No redirect to onboarding if profile missing
- ❌ Legacy dashboard visible to all users

**Fix**: Add PatientRoute wrapper or smart redirect logic

---

### **Issue 3: Mobile-First Not Implemented** 🟡 **HIGH PRIORITY**

**Problem**: Components not optimized for mobile viewport

**Impact**:
- ❌ Poor mobile UX ("long wait" - likely rendering/performance issues)
- ❌ Components may overflow or be hard to interact with
- ❌ Not aligned with "mobile-first" requirement

**Fix**: Audit and optimize all components for mobile

---

### **Issue 4: Legacy Routes Visible** 🟡 **MEDIUM PRIORITY**

**Problem**: Sidebar/navigation shows all routes including legacy/experimental ones

**Impact**:
- ❌ Confusing navigation for patients
- ❌ Not "MOAT navigation" as requested

**Fix**: Filter navigation based on user role, hide irrelevant routes

---

## ✅ ALIGNMENT WITH CLINICAL_MASTER_FRONTEND_AUDIT.md

### **Expected vs Actual:**

| Requirement | Expected | Actual | Status |
|-------------|----------|--------|--------|
| **Patient Landing** | `/ayesha-trials` with all capabilities | `/home` with legacy DisplayInfo | ❌ **BROKEN** |
| **Mobile-First** | Responsive, touch-friendly | Desktop-first, not optimized | ❌ **NOT IMPLEMENTED** |
| **MOAT Navigation** | Only relevant pages | All legacy routes visible | ❌ **NOT IMPLEMENTED** |
| **Holistic Scores** | Displayed in trial cards | ✅ **JUST ADDED** | ✅ **COMPLETE** |
| **All Components** | Mobile-optimized | Needs audit | ⚠️ **NEEDS AUDIT** |

---

## 🔧 IMMEDIATE FIXES (Priority Order)

### **Fix 1: Smart Redirect for Patients (30 minutes)**

**File**: `src/pages/Home.jsx`

**Change**:
- Add auth/patient context checks
- Redirect patients to `/ayesha-trials`
- Keep DisplayInfo for non-patient roles

---

### **Fix 2: Mobile Navigation Audit (2 hours)**

**Files**:
- `src/pages/AyeshaTrialExplorer.jsx`
- All component files in audit

**Changes**:
- Verify all components use responsive breakpoints correctly
- Ensure mobile-first breakpoints (start with mobile, add desktop)
- Test touch targets (44x44px minimum)
- Add collapsible sections where needed

---

### **Fix 3: PatientRoute Integration (30 minutes)**

**File**: `src/App.jsx`

**Changes**:
- Wrap `/ayesha-trials` and `/ayesha-complete-care` with PatientRoute (if not already)
- Ensure redirect logic works for unauthenticated patients

---

### **Fix 4: Navigation Filtering (1 hour)**

**File**: `src/components/Sidebar.jsx` (or wherever navigation is defined)

**Changes**:
- Filter routes based on user role
- Show only patient-relevant pages:
  - Trial Explorer
  - Complete Care
  - Profile
- Hide legacy/experimental routes

---

## 📱 MOBILE-FIRST CHECKLIST

### **AyeshaTrialExplorer Mobile Audit:**

- [ ] All components stack vertically on mobile (single column)
- [ ] Touch targets are ≥44x44px
- [ ] Text is readable (minimum 16px font size)
- [ ] No horizontal scrolling
- [ ] Cards have proper spacing on mobile
- [ ] HolisticScoreCard displays correctly on mobile
- [ ] Trial cards are scrollable (not all on one page)
- [ ] Navigation is accessible (bottom nav or hamburger menu)
- [ ] Loading states are mobile-friendly
- [ ] Error states are mobile-friendly

---

## 🎯 EXPECTED MOBILE EXPERIENCE

### **AK's First Load (After Fixes):**

1. **Login** → Redirect to `/ayesha-trials`
2. **Loading State** → Mobile-optimized spinner/ skeleton
3. **Profile Summary** → Compact card, single column
4. **SOC Recommendation** → Prominent card, touch-friendly
5. **CA-125 Tracker** → Collapsible if needed, clear visualization
6. **Trials List** → Scrollable cards with holistic scores visible
7. **Bottom Navigation** → Quick access to Complete Care, Profile

**No "book-like" static dashboard** - Everything is actionable and relevant to AK's care.

---

## 📋 FILES TO MODIFY

1. **`src/pages/Home.jsx`** - Add patient redirect logic
2. **`src/pages/AyeshaTrialExplorer.jsx`** - Mobile optimization audit
3. **`src/components/trials/HolisticScoreCard.jsx`** - Mobile layout audit
4. **`src/components/trials/TrialMatchCard.jsx`** - Mobile layout audit
5. **`src/App.jsx`** - Route protection (if needed)
6. **`src/components/Sidebar.jsx`** - Navigation filtering (if exists)

---

**Last Updated**: January 29, 2025  
**Status**: 🚨 **CRITICAL ISSUES IDENTIFIED - FIXES REQUIRED**  
**Priority**: **P1 - Fix patient redirect immediately**
