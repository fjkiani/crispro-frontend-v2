# ⚔️ SAAS TRANSFORMATION MASTER PLAN

**Commander:** Alpha  
**Architect:** Zo  
**Mission:** Transform CrisPRO platform into production SaaS with user management, feature flags, and data persistence  
**Timeline:** 2-3 weeks  
**Status:** 🎯 **READY FOR EXECUTION** (See `AUDIT_REPORT.md` for current state)

---

## 🎯 EXECUTIVE SUMMARY

**Current State:** Partial SaaS infrastructure exists - Supabase is integrated, but **no authentication, user management, or tier system** exists.

**Target State:** Production SaaS with:
- ✅ Multi-tenant user management (auth, roles, permissions)
- ✅ Feature flag system (free vs. paid tiers)
- ✅ Data persistence (user sessions, analyses, results)
- ✅ Usage tracking & quotas (rate limiting, credits)
- ✅ Billing integration (Stripe for subscriptions)
- ✅ Admin dashboard (manage users, features, quotas)

**Business Model:**
```
FREE TIER (Research Institutions):
- 10 variant analyses/month
- 5 drug efficacy queries/month
- 3 food validator queries/month
- Basic insights (no SAE features)
- Community support

PRO TIER ($499/month - Individual Researcher):
- 100 analyses/month
- Unlimited drug/food queries
- Full SAE features
- Clinical trials matching
- Priority support
- Export to PDF/CSV

ENTERPRISE TIER ($5,000/month - Academic Medical Center):
- Unlimited usage
- Custom integrations
- Dedicated Neo4j graph
- White-label options
- SLA guarantees
- 24/7 support
```

---

## 📊 CURRENT STATE AUDIT

**See:** `AUDIT_REPORT.md` for complete audit findings

### **✅ What Exists:**
- Supabase infrastructure (REST API client)
- Session management router (`/api/sessions`)
- Analysis history with Supabase
- Basic logging/analytics

### **❌ What's Missing:**
- Authentication system (no JWT, no login/signup)
- User management (no profiles, no tiers)
- Feature flags (only env-based, not user-based)
- Quotas & usage tracking
- Protected endpoints
- Billing integration

---

## 🏗️ COMPONENT STATUS

| Component | Status | Priority | Location | Notes |
|-----------|--------|----------|----------|-------|
| **1. Auth & User Management** | 🔴 Not Started | P0 | `components/1_auth/` | Supabase exists, needs Auth SDK |
| **2. Feature Flags** | 🔴 Not Started | P0 | `components/2_feature_flags/` | Env flags exist, need user-based |
| **3. Quotas & Usage Tracking** | 🔴 Not Started | P0 | `components/3_quotas/` | Logging exists, need quotas |
| **4. Session Persistence** | 🟡 Partial | P1 | `components/4_sessions/` | Sessions router exists, needs user linking |
| **5. Admin Dashboard** | 🔴 Not Started | P1 | `components/5_admin/` | Analytics dashboard exists |
| **6. Billing Integration** | 🔴 Not Started | P2 | `components/6_billing/` | No Stripe integration |

---

## 📋 EXECUTION PLAN

### **PHASE 1: FOUNDATION (Days 1-4) - Auth & Database**

#### **Day 1: Supabase Auth Setup** ✅
- [ ] Verify Supabase Auth is enabled in project
- [ ] Get JWT secret for token verification
- [ ] Run SaaS database schema (see `schemas/database_schema.sql`)
- [ ] Create indexes
- [ ] Test authentication (signup, login, logout)
- [ ] Configure email templates
- [ ] **Deliverable:** Supabase Auth + SaaS schema deployed

#### **Day 2: Backend Auth Integration** ✅
- [ ] Install dependencies (`pyjwt`, `python-multipart`, `supabase`)
- [ ] Extend `database_connections.py` with PostgreSQL support (if needed)
- [ ] Create `api/middleware/auth_middleware.py`:
  - JWT verification using Supabase JWT secret
  - `verify_token()` dependency
  - `get_current_user()` dependency
- [ ] Create `api/services/auth_service.py`:
  - `get_user_profile()` - Get user metadata
  - `update_user_profile()` - Update profile
- [ ] Create `api/routers/auth.py`:
  - `POST /api/auth/signup` - User registration (via Supabase Auth)
  - `POST /api/auth/login` - User login (via Supabase Auth)
  - `POST /api/auth/logout` - User logout
  - `GET /api/auth/profile` - Get user profile
  - `PUT /api/auth/profile` - Update profile
  - `POST /api/auth/refresh` - Refresh token
- [ ] **Deliverable:** Working auth endpoints

#### **Day 3: Frontend Auth** ✅
- [ ] Verify `@supabase/supabase-js` is installed (already installed)
- [ ] Create `src/context/AuthContext.jsx`:
  - Auth state management
  - Sign in/up/out functions using Supabase Auth
  - Session management
- [ ] Create `src/pages/auth/Login.jsx`:
  - Login form using Supabase Auth
  - Error handling
  - Redirect after login
- [ ] Create `src/pages/auth/Signup.jsx`:
  - Signup form using Supabase Auth
  - Profile creation in `user_profiles` table
  - Email confirmation flow
- [ ] Create `src/components/auth/ProtectedRoute.jsx`:
  - Route protection
  - Redirect to login if not authenticated
- [ ] Update `src/App.jsx`:
  - Add AuthProvider wrapper
  - Add protected routes
  - Add login/signup routes
- [ ] **Deliverable:** Working auth UI

#### **Day 4: Protected Endpoints** ✅
- [ ] Add auth middleware to existing endpoints:
  - `POST /api/insights/*` → require auth
  - `POST /api/efficacy/*` → require auth
  - `POST /api/design/*` → require auth
  - `POST /api/sessions/*` → require auth (update existing router)
- [ ] Update `sessions.py` router to use `Depends(get_current_user)`
- [ ] Update `AnalysisHistoryContext` to use authenticated user ID
- [ ] Test with Postman/curl
- [ ] **Deliverable:** All endpoints require authentication

---

### **PHASE 2: FEATURE FLAGS & QUOTAS (Days 5-7)**

#### **Day 5: Feature Flag System** ✅
- [ ] Create `api/services/feature_flag_service.py`:
  - `get_user_features()` - Get enabled features for user
  - `has_feature()` - Check if user has feature
  - Tier-based feature mapping
- [ ] Create `api/middleware/feature_flag_middleware.py`:
  - `require_feature()` dependency
- [ ] Add feature checks to premium endpoints:
  - `POST /api/efficacy/predict` → require `sae_features` for Pro+
  - `POST /api/datasets/extract_and_benchmark` → require `cohort_lab` for Enterprise
  - `POST /api/design/generate_guide_rna` → require `crispr_design` for Enterprise
- [ ] Test tier-based access
- [ ] **Deliverable:** Feature flag system operational

#### **Day 6: Quota System** ✅
- [ ] Create `api/services/quota_service.py`:
  - `get_user_quotas()` - Get quota usage
  - `check_quota()` - Check if user has quota
  - `increment_usage()` - Increment usage counter
  - `reset_quotas_if_needed()` - Monthly reset
- [ ] Create `api/middleware/quota_middleware.py`:
  - `check_quota()` dependency
- [ ] Add quota checks to all endpoints:
  - `POST /api/insights/*` → check `variant_analyses` quota
  - `POST /api/efficacy/*` → check `drug_queries` quota
  - `POST /api/trials/*` → check `clinical_trials` quota
- [ ] Test quota limits (free tier)
- [ ] **Deliverable:** Quota enforcement working

#### **Day 7: Usage Tracking** ✅
- [ ] Create `api/services/usage_tracking_service.py`:
  - Log usage to `usage_logs` table
- [ ] Add usage logging to all endpoints
- [ ] Create `src/pages/UsageDashboard.jsx`:
  - Display quota usage
  - Show remaining quotas
  - Upgrade prompts
- [ ] Test quota warnings
- [ ] **Deliverable:** Usage dashboard visible

---

### **PHASE 3: SESSION PERSISTENCE (Days 8-9)**

#### **Day 8: Backend Session Storage** ✅
- [ ] Update `sessions.py` router to link sessions to authenticated users
- [ ] Create `saved_analyses` table (if not exists)
- [ ] Implement save/load analysis with user linking
- [ ] Test session persistence
- [ ] **Deliverable:** Sessions API working with user auth

#### **Day 9: Frontend Session UI** ✅
- [ ] Create `src/pages/MyAnalyses.jsx`:
  - List saved analyses for current user
  - Load/resume analyses
  - Delete analyses
- [ ] Update `AnalysisHistoryContext` to use authenticated user
- [ ] Add save/load buttons to analysis pages
- [ ] Add session history
- [ ] Test cross-page resume
- [ ] **Deliverable:** Session UI complete

---

### **PHASE 4: ADMIN & BILLING (Days 10-14)**

#### **Day 10-11: Admin Dashboard** ✅
- [ ] Create admin authentication (admin role check)
- [ ] Create `api/routers/admin.py`:
  - User management endpoints
  - Usage analytics endpoints
  - Feature flag management endpoints
  - Quota override endpoints
- [ ] Create `src/pages/admin/` pages:
  - User management
  - Usage analytics
  - Feature flag management
- [ ] **Deliverable:** Admin panel functional

#### **Day 12-13: Billing Integration** ✅
- [ ] Create Stripe account and get API keys
- [ ] Create `api/services/stripe_service.py`:
  - Subscription creation
  - Webhook handling
  - Payment processing
- [ ] Create `api/routers/billing.py`:
  - Create checkout session
  - Webhook handler
  - Subscription management
- [ ] Create `src/pages/Billing.jsx`:
  - Pricing display
  - Checkout flow
  - Subscription management
- [ ] Test subscription lifecycle
- [ ] **Deliverable:** Billing system operational

#### **Day 14: Testing & Polish** ✅
- [ ] End-to-end testing:
  - Free tier signup → usage → quota hit
  - Upgrade flow
  - Premium features
  - Session persistence
- [ ] Security audit
- [ ] Performance testing
- [ ] **Deliverable:** Production-ready SaaS

---

## 🔧 TECHNICAL STACK

### **Backend:**
- **Database:** PostgreSQL (Supabase) + Redis (caching)
- **Auth:** Supabase Auth (JWT)
- **API:** FastAPI with middleware
- **Billing:** Stripe

### **Frontend:**
- **Auth:** Supabase JS client (already installed)
- **State:** React Context (AuthContext)
- **UI:** Material-UI components

### **Infrastructure:**
- **Hosting:** Existing backend (FastAPI)
- **Database:** Supabase (managed PostgreSQL)
- **Cache:** Redis (existing)

---

## 📊 COMPONENT DETAILS

### **Component 1: Auth & User Management**
**Location:** `components/1_auth/`

**Files:**
- `api/middleware/auth_middleware.py` - JWT verification (NEW)
- `api/services/auth_service.py` - Auth operations (NEW)
- `api/routers/auth.py` - Auth endpoints (NEW)
- `src/context/AuthContext.jsx` - Frontend auth (NEW)
- `src/pages/auth/Login.jsx` - Login page (NEW)
- `src/pages/auth/Signup.jsx` - Signup page (NEW)

**Status:** 🔴 Not Started

**Dependencies:**
- Supabase Auth SDK (`supabase` Python package)
- JWT secret from Supabase

---

### **Component 2: Feature Flags**
**Location:** `components/2_feature_flags/`

**Files:**
- `api/services/feature_flag_service.py` - Feature flag logic (NEW)
- `api/middleware/feature_flag_middleware.py` - Feature checks (NEW)
- Database schema: `features` and `user_feature_flags` tables

**Status:** 🔴 Not Started

**Note:** Existing env-based flags in `config.py` should remain for system-level flags. User-based flags are for tier gating.

---

### **Component 3: Quotas & Usage**
**Location:** `components/3_quotas/`

**Files:**
- `api/services/quota_service.py` - Quota management (NEW)
- `api/middleware/quota_middleware.py` - Quota checks (NEW)
- `api/services/usage_tracking_service.py` - Usage logging (NEW)
- `src/pages/UsageDashboard.jsx` - Usage UI (NEW)

**Status:** 🔴 Not Started

---

### **Component 4: Session Persistence**
**Location:** `components/4_sessions/`

**Files:**
- `api/routers/sessions.py` - Session endpoints (EXISTS - needs update)
- `src/context/AnalysisHistoryContext.jsx` - Analysis history (EXISTS - needs update)
- `src/pages/MyAnalyses.jsx` - Saved analyses UI (NEW)

**Status:** 🟡 Partial (exists but needs user linking)

---

### **Component 5: Admin Dashboard**
**Location:** `components/5_admin/`

**Files:**
- `api/routers/admin.py` - Admin endpoints (NEW)
- `src/pages/admin/` - Admin UI pages (NEW)
- `src/components/analytics/SupabaseDashboard.jsx` - Analytics (EXISTS)

**Status:** 🔴 Not Started

---

### **Component 6: Billing Integration**
**Location:** `components/6_billing/`

**Files:**
- `api/routers/billing.py` - Billing endpoints (NEW)
- `api/services/stripe_service.py` - Stripe integration (NEW)
- `src/pages/Billing.jsx` - Billing UI (NEW)

**Status:** 🔴 Not Started

---

## 🗄️ DATABASE SCHEMA

**Location:** `schemas/database_schema.sql`

**Tables:**
1. `auth.users` (managed by Supabase Auth)
2. `public.user_profiles` (NEW)
3. `public.user_subscriptions` (NEW)
4. `public.user_quotas` (NEW)
5. `public.user_feature_flags` (NEW)
6. `public.features` (NEW)
7. `public.user_sessions` (EXISTS - needs update)
8. `public.saved_analyses` (NEW - different from `analysis_history`)
9. `public.usage_logs` (NEW)

**Existing Tables (Keep):**
- `mdt_runs`, `mdt_events`, `mdt_run_variants` (logging)
- `analysis_history` (frontend analysis history)
- `session_items` (session items)

**See:** `schemas/database_schema.sql` for full schema

---

## 🧪 TESTING PLAN

**Location:** `tests/`

**Test Files:**
- `tests/test_auth.py` - Auth endpoints
- `tests/test_feature_flags.py` - Feature flag system
- `tests/test_quotas.py` - Quota enforcement
- `tests/test_sessions.py` - Session persistence
- `tests/test_admin.py` - Admin endpoints
- `tests/test_billing.py` - Billing integration

---

## 📝 NOTES & KNOWN ISSUES

### **Current State (From Audit):**
- ✅ Supabase infrastructure exists (`supabase_service.py`)
- ✅ Session router exists (`sessions.py`) but uses optional `x_user_id`
- ✅ Analysis history exists but uses `user_id: null`
- ❌ No authentication system
- ❌ No user profiles or tiers
- ❌ No feature flags (user-based)
- ❌ No quotas

### **Integration Points:**
1. **Sessions Router:** Update to use `Depends(get_current_user)` instead of optional header
2. **Analysis History:** Update `AnalysisHistoryService` to use authenticated user ID
3. **Supabase Service:** Keep existing service, add Auth SDK for authentication
4. **Database:** Run SaaS schema, keep existing tables for logging

### **Manager Questions (See `AUDIT_REPORT.md`):**
1. Supabase Auth setup status?
2. Existing data migration strategy?
3. Feature flag migration strategy?
4. Authentication requirements?

---

## 🚀 QUICK START

### **For Agent Execution:**

1. **Review Audit:** Read `AUDIT_REPORT.md` for current state
2. **Answer Manager Questions:** Get answers to questions in audit
3. **Start Component 1 (Auth):**
   ```bash
   cd components/1_auth/
   # Follow README.md in that folder
   ```
4. **Then Component 2 (Feature Flags):**
   ```bash
   cd components/2_feature_flags/
   # Follow README.md in that folder
   ```
5. **Continue sequentially through components**

---

## 📊 PROGRESS TRACKING

**Last Updated:** 2024-12-XX

**Overall Progress:** 0% (0/6 components complete)

**Current Phase:** Phase 1 - Foundation

**Next Milestone:** Auth system operational

**Audit Status:** ✅ Complete (see `AUDIT_REPORT.md`)

---

**This master plan is the single source of truth for SaaS transformation. All component work should update this document.**
