# Persona-Based Access Control System

**Date:** January 2025  
**Status:** ✅ **IMPLEMENTATION PLAN**  
**Purpose:** Define access control for 3 personas: Patient, Oncologist, R&D/Researcher

---

## 🎯 Persona Definitions

### 1. **Patient**
- **Purpose:** Patients accessing their own medical information and care plans
- **Access Level:** Limited to personal data and patient-facing features
- **Database Role:** `patient` (maps to `role` field in `user_profiles`)

### 2. **Oncologist** 
- **Purpose:** Clinical oncologists managing patient care
- **Access Level:** Clinical tools, patient management, treatment planning
- **Database Role:** `clinician` or `oncologist` (maps to `role` field in `user_profiles`)

### 3. **R&D / Researcher**
- **Purpose:** Researchers and R&D teams using advanced research tools
- **Access Level:** Full research capabilities, experimental features, data analysis
- **Database Role:** `researcher` (maps to `role` field in `user_profiles`)

---

## 📋 Access Matrix

### **Patient Persona**

#### ✅ **Allowed Pages/Features:**
- `/patient/profile` - View own profile
- `/patient/onboarding` - Complete onboarding
- `/ayesha-complete-care` - View own care plan
- `/ayesha-trials` - View matching trials
- `/patient/tasks` - View assigned tasks
- Patient-specific medical records
- Basic drug efficacy predictions (own data only)

#### ❌ **Restricted Pages/Features:**
- `/admin/*` - Admin dashboard
- `/tools` - Research tools
- `/crispr-designer` - CRISPR design tools
- `/orchestrator` - Orchestrator dashboard
- `/research` - Research portal
- `/clinical-genomics` - Clinical genomics command center
- `/agent-dashboard` - Agent dashboard
- `/doctor-dashboard` - Doctor dashboard
- `/universal-dossiers` - Universal dossier browser
- `/universal-trial-intelligence` - Universal trial intelligence
- `/metastasis` - Metastasis dashboard
- `/synthetic-lethality` - Synthetic lethality analyzer
- `/dosing-guidance` - Dosing guidance (read-only, own data)
- `/validate` - Hypothesis validator
- `/myeloma-digital-twin` - Myeloma digital twin
- `/sporadic-cancer` - Sporadic cancer page
- `/radonc-co-pilot` - RadOnc co-pilot
- `/clinical-dossier-test` - Clinical dossier test
- `/orchestrator` - Orchestrator dashboard
- `/protein-synthesis` - Protein synthesis
- `/structure-predictor` - Structure predictor
- `/demo-summarizer` - Demo summarizer
- `/campaigns/*` - Campaign runner
- `/runx-conquest/*` - Runx conquest
- `/q2c-test` - Q2C test
- `/phase3-demo` - Phase 3 demo
- `/copilot-smoke-test` - CoPilot smoke test
- `/copilot-gap-analysis` - CoPilot gap analysis

---

### **Oncologist Persona**

#### ✅ **Allowed Pages/Features:**
- `/patient/profile` - View patient profiles
- `/patient/onboarding` - Patient onboarding
- `/ayesha-complete-care` - Complete care plans
- `/ayesha-trials` - Trial matching
- `/patient/tasks` - Patient task management
- `/medical-records` - Medical records (all patients)
- `/medical-records/:id` - Individual record details
- `/medical-records/:patientId/research` - Patient research
- `/medical-records/:patientId/tasks` - Patient tasks
- `/research` - Research portal
- `/clinical-genomics` - Clinical genomics command center
- `/validate` - Hypothesis validator
- `/myeloma-digital-twin` - Myeloma digital twin
- `/metastasis` - Metastasis dashboard
- `/synthetic-lethality` - Synthetic lethality analyzer
- `/dosing-guidance` - Dosing guidance (full access)
- `/threat-assessor` - Threat assessor
- `/radonc-co-pilot` - RadOnc co-pilot
- `/universal-dossiers` - Universal dossier browser
- `/universal-trial-intelligence` - Universal trial intelligence
- `/doctor-dashboard` - Doctor dashboard
- `/workload-dashboard` - Follow-up task board
- `/screening-schedules` - Screening schedules
- `/outreach` - Outreach dashboard
- `/mutation-explorer` - Mutation explorer
- `/agent-dashboard` - Agent dashboard
- `/agents` - Agents page
- `/profile` - Own profile
- `/home` - Home page
- `/dashboard` - Dashboard

#### ❌ **Restricted Pages/Features:**
- `/admin/*` - Admin dashboard (unless admin role)
- `/tools` - Research tools (unless researcher role)
- `/crispr-designer` - CRISPR design tools
- `/orchestrator` - Orchestrator dashboard
- `/protein-synthesis` - Protein synthesis
- `/structure-predictor` - Structure predictor
- `/demo-summarizer` - Demo summarizer
- `/campaigns/*` - Campaign runner
- `/runx-conquest/*` - Runx conquest
- `/q2c-test` - Q2C test
- `/phase3-demo` - Phase 3 demo
- `/copilot-smoke-test` - CoPilot smoke test
- `/copilot-gap-analysis` - CoPilot gap analysis
- `/agent-studio` - Agent studio
- `/agent-demo/:agentId` - Agent demo (unless DEV mode)

---

### **R&D / Researcher Persona**

#### ✅ **Allowed Pages/Features:**
- **ALL Pages** - Full access to all features
- `/tools` - Research tools (Armory)
- `/crispr-designer` - CRISPR designer
- `/orchestrator` - Orchestrator dashboard
- `/research` - Research portal
- `/clinical-genomics` - Clinical genomics command center
- `/validate` - Hypothesis validator
- `/myeloma-digital-twin` - Myeloma digital twin
- `/metastasis` - Metastasis dashboard
- `/synthetic-lethality` - Synthetic lethality analyzer
- `/dosing-guidance` - Dosing guidance (with validation tab)
- `/threat-assessor` - Threat assessor
- `/sporadic-cancer` - Sporadic cancer page
- `/radonc-co-pilot` - RadOnc co-pilot
- `/universal-dossiers` - Universal dossier browser
- `/universal-trial-intelligence` - Universal trial intelligence
- `/doctor-dashboard` - Doctor dashboard
- `/workload-dashboard` - Follow-up task board
- `/screening-schedules` - Screening schedules
- `/outreach` - Outreach dashboard
- `/mutation-explorer` - Mutation explorer
- `/agent-dashboard` - Agent dashboard
- `/agents` - Agents page
- `/agent-studio` - Agent studio
- `/agent-demo/:agentId` - Agent demo
- `/protein-synthesis` - Protein synthesis
- `/structure-predictor` - Structure predictor
- `/demo-summarizer` - Demo summarizer
- `/campaigns/*` - Campaign runner
- `/runx-conquest/*` - Runx conquest
- `/clinical-dossier-test` - Clinical dossier test
- `/profile` - Own profile
- `/home` - Home page
- `/dashboard` - Dashboard
- `/admin/*` - Admin dashboard (if admin role)

#### ❌ **Restricted Pages/Features:**
- None (full access)

---

## 🔧 Implementation

### **1. Database Schema Update**

Update `user_profiles` table to support persona:
```sql
-- Add persona field (if not exists)
ALTER TABLE public.user_profiles 
ADD COLUMN IF NOT EXISTS persona VARCHAR(50) 
CHECK (persona IN ('patient', 'oncologist', 'researcher')) 
DEFAULT 'researcher';

-- Map existing roles to personas
UPDATE public.user_profiles 
SET persona = CASE 
  WHEN role = 'clinician' THEN 'oncologist'
  WHEN role = 'researcher' THEN 'researcher'
  WHEN role = 'admin' THEN 'researcher'
  ELSE 'patient'
END
WHERE persona IS NULL;
```

### **2. Frontend Persona Context**

Create `src/context/PersonaContext.jsx`:
- Provides current user's persona
- Provides `hasAccess(page)` function
- Provides `canAccessFeature(feature)` function

### **3. Persona-Based Route Guards**

Create `src/components/auth/PersonaRoute.jsx`:
- Wraps routes with persona-based access control
- Redirects unauthorized users to appropriate page

### **4. Access Control Hook**

Create `src/hooks/usePersonaAccess.js`:
- Returns access control functions
- Checks persona against access matrix

---

## 📝 Usage Examples

### **In Route Definitions:**
```jsx
<Route 
  path="/tools" 
  element={
    <PersonaRoute allowedPersonas={['researcher']}>
      <Armory />
    </PersonaRoute>
  } 
/>

<Route 
  path="/patient/profile" 
  element={
    <PersonaRoute allowedPersonas={['patient', 'oncologist']}>
      <PatientProfile />
    </PersonaRoute>
  } 
/>
```

### **In Components:**
```jsx
import { usePersonaAccess } from '../hooks/usePersonaAccess';

function MyComponent() {
  const { canAccess, persona } = usePersonaAccess();
  
  if (!canAccess('crispr-designer')) {
    return <AccessDenied />;
  }
  
  return <CrisprDesigner />;
}
```

---

## 🔐 Security Notes

1. **Backend Validation:** Always validate persona on backend endpoints
2. **Default Deny:** Unknown personas default to most restrictive access
3. **Audit Logging:** Log all access attempts for security auditing
4. **Role Mapping:** Map database roles to personas consistently
5. **Feature Flags:** Combine persona access with tier-based feature flags

---

## 📊 Migration Path

1. **Phase 1:** Add persona field to database
2. **Phase 2:** Create PersonaContext and hooks
3. **Phase 3:** Add PersonaRoute component
4. **Phase 4:** Update App.jsx routes with PersonaRoute
5. **Phase 5:** Add backend persona validation
6. **Phase 6:** Update UI to show persona-specific features

---

**Last Updated:** January 2025  
**Status:** Ready for Implementation

