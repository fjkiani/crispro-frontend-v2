<!-- 08c68ec1-d889-4a20-8bc9-17363ed77e9f 7ca913c2-a983-4d60-9b1b-15d649151f9d -->
# SaaS Transformation - Deep Audit & Compliance Plan

## Executive Summary

**Current State (Code-Validated):**

- ✅ Authentication system: COMPLETE (auth_middleware.py, auth_service.py, auth router)
- ✅ Admin system: COMPLETE (admin_middleware.py, admin_service.py, admin router)
- ✅ Session management: EXISTS (sessions.py with optional auth)
- ⚠️ Security: CORS allows all origins (*), no encryption enforcement
- ❌ HIPAA/PII compliance: NOT IMPLEMENTED (only plans exist)
- ❌ Quota system: NOT IMPLEMENTED (no quota_service.py found)
- ❌ Feature flags: NOT IMPLEMENTED (no feature_flag_service.py found)
- ⚠️ Admin promotion: Manual SQL only (no UI)

**Critical Gaps Identified:**

1. No HIPAA/PII protection middleware
2. No encryption at rest/in transit enforcement
3. No audit logging for PHI access
4. No MFA implementation
5. CORS allows all origins (security risk)
6. RLS policies defined but not verified as active
7. No quota enforcement
8. No feature flag system
9. Admin promotion requires manual SQL

---

## Part 1: Code-Validated Current State

### What Actually Exists (Code Inspection)

#### Authentication System ✅ COMPLETE

**Files Found:**

- `api/middleware/auth_middleware.py` - JWT verification, get_current_user(), get_optional_user()
- `api/services/auth_service.py` - Signup/login/profile operations
- `api/routers/auth.py` - All auth endpoints (signup, login, logout, profile, refresh)
- `src/context/AuthContext.jsx` - Frontend auth state
- `src/pages/auth/Login.jsx` - Login page
- `src/pages/auth/Signup.jsx` - Signup page

**Code Evidence:**

- JWT verification uses SUPABASE_JWT_SECRET
- Auth endpoints functional (signup, login, profile)
- Optional auth pattern used: `get_optional_user()` allows anonymous access
- Sessions router uses optional auth (line 77: `user: Optional[Dict] = Depends(get_optional_user)`)

#### Admin System ✅ COMPLETE

**Files Found:**

- `api/middleware/admin_middleware.py` - require_admin(), require_admin_or_self()
- `api/services/admin_service.py` - User management, analytics, activity tracking
- `api/routers/admin.py` - Admin endpoints (users, analytics, activity)
- `src/pages/admin/Dashboard.jsx` - Basic admin dashboard
- `src/pages/admin/Users.jsx` - User list page

**Code Evidence:**

- Admin middleware checks user_profiles.role == 'admin' (line 39-40)
- Admin service has get_users(), get_user_details(), update_user(), get_analytics_overview()
- Admin endpoints protected with `Depends(require_admin)`
- Admin promotion: Manual SQL update required (no UI endpoint)

#### Session Management ✅ EXISTS

**Files Found:**

- `api/routers/sessions.py` - Session CRUD operations
- Uses optional auth: supports both authenticated and anonymous sessions
- Tables: `user_sessions`, `session_items` (referenced in code)

#### Database Schema ✅ DEFINED

**Files Found:**

- `.cursor/rules/saas_transformation/schemas/database_schema.sql`
- Tables defined: user_profiles, user_subscriptions, user_quotas, user_feature_flags, features, user_sessions, saved_analyses, usage_logs
- RLS policies defined but NOT VERIFIED as active in Supabase

### What's Missing (Code Inspection)

#### Quota System ❌ NOT FOUND

**Expected Files:**

- `api/services/quota_service.py` - NOT FOUND
- `api/middleware/quota_middleware.py` - NOT FOUND
- No quota enforcement in endpoints

**Evidence:**

- `grep -r "quota_service\|check_quota"` returns no matches
- Endpoints don't check quotas before processing

#### Feature Flag System ❌ NOT FOUND

**Expected Files:**

- `api/services/feature_flag_service.py` - NOT FOUND
- `api/middleware/feature_flag_middleware.py` - NOT FOUND
- Only env-based flags exist in config.py

**Evidence:**

- `grep -r "feature_flag_service\|has_feature"` returns no matches
- No tier-based feature gating

#### HIPAA/PII Compliance ❌ NOT IMPLEMENTED

**Expected Files:**

- `api/middleware/hipaa_pii.py` - NOT FOUND
- `api/audit/writer.py` - NOT FOUND
- `api/security/secrets.py` - NOT FOUND
- `api/utils/data_sanitizer.py` - NOT FOUND

**Evidence:**

- Only plans exist: `docs/hipaa_compliance_plan.md`, `.cursor/rules/remaining_categories/archive/agent_build_hipaa_foundation.mdc`
- No PHI/PII redaction middleware
- No audit logging for PHI access
- No encryption enforcement

#### Security Issues ⚠️ FOUND

**Code Evidence:**

- `api/main.py` line 141: `allow_origins=["*"]` - CORS allows all origins
- No HTTPS enforcement
- No security headers middleware
- No rate limiting per user

---

## Part 2: Critical Gaps Analysis

### Gap 1: HIPAA/PII Compliance (CRITICAL)

**Current State:**

- No PHI/PII detection or redaction
- No audit logging for PHI access
- No encryption at rest/in transit enforcement
- No data classification (PHI vs NON_PHI)
- No retention policies
- No MFA

**Required Implementation:**

1. PHI/PII Detection Middleware

   - File: `api/middleware/hipaa_pii.py`
   - Detect: names, emails, phone, MRN, DOB, address, SSN, genomic data
   - Redact from logs only (don't mutate requests)
   - Configurable via env: `HIPAA_PHI_FIELDS`

2. Audit Trail System

   - File: `api/audit/writer.py`
   - Append-only audit log with SHA-256 hash chaining
   - Log all PHI access: user_id, action, resource_type, timestamp
   - Table: `hipaa_audit_log` (from schema plan)

3. Encryption Enforcement

   - File: `api/middleware/security_headers.py`
   - Enforce HTTPS in production
   - Set HSTS, secure cookies, strict CORS
   - At-rest encryption: Supabase handles (verify)

4. Data Classification

   - Add `data_classification` field to all tables storing PHI
   - Values: 'PHI', 'NON_PHI', 'DERIVED'
   - Auto-classify genomic data as PHI

5. MFA Implementation

   - Integrate Supabase Auth MFA
   - Require MFA for admin users
   - Require MFA for PHI access

### Gap 2: Admin User Management (HIGH)

**Current State:**

- Admin promotion: Manual SQL only
- No admin creation endpoint
- No admin role assignment UI
- No admin audit logging

**Required Implementation:**

1. Admin Promotion Endpoint

   - File: `api/routers/admin.py` - Add `POST /api/admin/users/{user_id}/promote`
   - Requires existing admin role
   - Logs promotion action in audit log
   - Updates user_profiles.role = 'admin'

2. Admin Creation UI

   - File: `src/pages/admin/Users.jsx` - Add "Promote to Admin" button
   - Confirmation dialog
   - Success/error feedback

3. Admin Audit Logging

   - Log all admin actions: user updates, promotions, suspensions
   - Table: `admin_audit_log` (new table)
   - Fields: admin_user_id, action, target_user_id, changes, timestamp

4. Admin Role Hierarchy

   - Super Admin: Can promote/demote admins
   - Admin: Can manage users, view analytics
   - Consider: `is_super_admin` flag in user_profiles

### Gap 3: Quota System (HIGH)

**Current State:**

- Quota tables exist in schema
- No quota service implementation
- No quota enforcement

**Required Implementation:**

1. Quota Service

   - File: `api/services/quota_service.py`
   - Methods: get_user_quotas(), check_quota(), increment_usage(), reset_quotas_if_needed()
   - Tier-based limits: free (10/5/3), pro (100/unlimited/50), enterprise (unlimited)

2. Quota Middleware

   - File: `api/middleware/quota_middleware.py`
   - Dependency: `check_quota(quota_type)`
   - Returns 429 if quota exceeded
   - Headers: X-Quota-Limit, X-Quota-Used, X-Quota-Reset

3. Endpoint Integration

   - Add quota checks to: `/api/efficacy/predict`, `/api/insights/*`, `/api/design/*`
   - Map endpoints to quota types: variant_analyses, drug_queries, food_queries

### Gap 4: Feature Flag System (MEDIUM)

**Current State:**

- Feature tables exist in schema
- Only env-based flags exist
- No tier-based feature gating

**Required Implementation:**

1. Feature Flag Service

   - File: `api/services/feature_flag_service.py`
   - Methods: get_user_features(), has_feature(), enable_feature(), disable_feature()
   - Tier mapping: free → basic features, pro → SAE features, enterprise → all features

2. Feature Flag Middleware

   - File: `api/middleware/feature_flag_middleware.py`
   - Dependency: `require_feature(feature_name)`
   - Returns 403 if feature not available

3. Endpoint Integration

   - Add feature checks to premium endpoints
   - `/api/efficacy/predict` → require 'sae_features' for Pro+
   - `/api/datasets/extract_and_benchmark` → require 'cohort_lab' for Enterprise

### Gap 5: Security Hardening (CRITICAL)

**Current State:**

- CORS allows all origins (*)
- No HTTPS enforcement
- No security headers
- No rate limiting

**Required Implementation:**

1. CORS Configuration

   - File: `api/main.py` - Update CORS middleware
   - Allow specific origins: `allow_origins=[os.getenv("ALLOWED_ORIGINS", "http://localhost:3000")]`
   - Environment-based: dev vs prod

2. Security Headers Middleware

   - File: `api/middleware/security_headers.py`
   - Set: HSTS, X-Content-Type-Options, Referrer-Policy, X-Frame-Options
   - Enforce HTTPS in production

3. Rate Limiting

   - File: `api/middleware/rate_limit.py`
   - Per-user rate limits (Redis-based)
   - Tier-based limits: free (10/min), pro (100/min), enterprise (unlimited)

---

## Part 3: HIPAA/PII Compliance Roadmap

### Phase 1: Foundation (Week 1-2)

**Deliverables:**

1. PHI/PII Detection Middleware

   - Detect common PHI patterns (regex + heuristics)
   - Redact from logs only
   - Configurable via env

2. Structured Audit Logging

   - JSON logging with fixed schema
   - Append-only audit sink
   - SHA-256 hash chaining

3. Security Headers

   - HSTS, secure cookies, strict CORS
   - HTTPS enforcement

**Files to Create:**

- `api/middleware/hipaa_pii.py`
- `api/audit/writer.py`
- `api/middleware/security_headers.py`
- `api/utils/logging.py` (structured JSON logger)

### Phase 2: Data Protection (Week 3-4)

**Deliverables:**

1. Data Classification

   - Add `data_classification` field to tables
   - Auto-classify genomic data as PHI
   - Classification helper functions

2. Encryption Enforcement

   - Verify Supabase encryption at rest
   - Enforce TLS 1.3 for all connections
   - Key management interface (Vault/KMS shim)

3. Access Control

   - RBAC enhancements
   - PHI access logging
   - Role-based PHI access restrictions

**Files to Create:**

- `api/security/secrets.py`
- `api/utils/data_sanitizer.py`
- Database migration: Add data_classification columns

### Phase 3: Advanced Features (Week 5-6)

**Deliverables:**

1. MFA Implementation

   - Supabase Auth MFA integration
   - Require MFA for admin users
   - Require MFA for PHI access

2. Data Retention

   - Retention policy helpers
   - Automated deletion (cron job)
   - Data Subject Request (DSR) utilities

3. Monitoring & Alerts

   - SIEM integration hooks
   - Alert thresholds (configurable)
   - Incident response runbook

**Files to Create:**

- `api/security/mfa.py`
- `api/utils/retention.py`
- `docs/hipaa/incident_response.md`

---

## Part 4: Admin Management Architecture

### Admin User Creation Flow

**Option 1: Manual Promotion (Current)**

```
1. User signs up normally
2. Admin runs SQL: UPDATE user_profiles SET role = 'admin' WHERE email = '...'
3. User logs in → admin access granted
```

**Option 2: Admin Promotion Endpoint (Recommended)**

```
1. Existing admin calls: POST /api/admin/users/{user_id}/promote
2. Backend checks: requester is admin
3. Backend updates: user_profiles.role = 'admin'
4. Backend logs: admin_audit_log entry
5. Returns: success confirmation
```

**Option 3: Super Admin Creation (First Admin)**

```
1. Seed script: python scripts/create_first_admin.py
2. Creates user via Supabase Auth
3. Sets role = 'admin' in user_profiles
4. Returns: credentials (one-time use)
```

### Admin Capabilities Matrix

| Capability | Admin | Super Admin |

|------------|-------|-------------|

| View users | ✅ | ✅ |

| Update users | ✅ | ✅ |

| Suspend users | ✅ | ✅ |

| Promote to admin | ❌ | ✅ |

| Demote admin | ❌ | ✅ |

| View audit logs | ✅ | ✅ |

| Manage feature flags | ✅ | ✅ |

| Override quotas | ✅ | ✅ |

| Delete users | ❌ | ✅ |

### Admin UI Components Needed

1. **User Management Page** (`src/pages/admin/Users.jsx`)

   - User list table (EXISTS - needs enhancement)
   - User detail modal/page
   - Promote to admin button (NEW)
   - Suspend/activate buttons (EXISTS)

2. **Admin Settings Page** (`src/pages/admin/Settings.jsx`) - NEW

   - Admin role management
   - Super admin designation
   - Admin audit log viewer

3. **Admin Audit Log** (`src/pages/admin/AuditLog.jsx`) - NEW

   - All admin actions
   - Filterable by admin, action, date
   - Export to CSV

---

## Part 5: Implementation Priority

### P0 (Critical - Security)

1. CORS configuration fix (allow_origins)
2. Security headers middleware
3. PHI/PII detection middleware
4. Audit logging foundation

### P1 (High - Core Features)

1. Quota system implementation
2. Feature flag system implementation
3. Admin promotion endpoint
4. RLS policy verification

### P2 (Medium - Compliance)

1. MFA implementation
2. Data classification
3. Encryption enforcement
4. Retention policies

### P3 (Low - Enhancements)

1. Admin UI enhancements
2. Analytics charts
3. Export functionality
4. Advanced admin features

---

## Part 6: Developer Implementation Guide

### Step 1: Security Hardening (Day 1)

1. Fix CORS in `api/main.py`
2. Create `api/middleware/security_headers.py`
3. Add security headers to FastAPI app
4. Test: Verify headers in response

### Step 2: HIPAA Foundation (Days 2-3)

1. Create `api/middleware/hipaa_pii.py`
2. Create `api/audit/writer.py`
3. Create `api/utils/logging.py`
4. Register middleware in `api/main.py` (when HIPAA_MODE=true)
5. Test: Verify no PHI in logs

### Step 3: Quota System (Days 4-5)

1. Create `api/services/quota_service.py`
2. Create `api/middleware/quota_middleware.py`
3. Add quota checks to endpoints
4. Test: Verify quota enforcement

### Step 4: Feature Flags (Day 6)

1. Create `api/services/feature_flag_service.py`
2. Create `api/middleware/feature_flag_middleware.py`
3. Add feature checks to premium endpoints
4. Test: Verify tier-based access

### Step 5: Admin Enhancements (Day 7)

1. Add admin promotion endpoint
2. Create admin audit log table
3. Add promote button to admin UI
4. Test: Verify admin promotion flow

---

## Validation Checklist

### Code Inspection Validation

- [x] Read actual code files (not just docs)
- [x] Traced execution paths (auth flow, admin flow)
- [x] Mapped integration points (auth → sessions → admin)
- [x] Documented actual state vs planned
- [x] Identified real gaps from code inspection
- [x] Added code references for each finding

### Security Validation

- [ ] CORS fixed (specific origins)
- [ ] Security headers present
- [ ] HTTPS enforced in production
- [ ] Rate limiting implemented
- [ ] PHI/PII redaction working
- [ ] Audit logging functional

### Compliance Validation

- [ ] HIPAA middleware active
- [ ] Audit trail hash chain valid
- [ ] MFA required for admin
- [ ] Data classification complete
- [ ] Retention policies enforced
- [ ] RLS policies active in Supabase

---

## Next Steps

1. **Immediate:** Fix CORS configuration (5 min)
2. **Week 1:** Implement HIPAA foundation (PHI detection, audit logging)
3. **Week 2:** Implement quota and feature flag systems
4. **Week 3:** Admin enhancements and security hardening
5. **Week 4:** MFA and advanced compliance features

**This plan is code-validated and ready for implementation.**

### To-dos

- [ ] Generate Exploit Vector Analysis (Offensive Disclosure).
- [ ] Generate Mitigation Pseudocode (Defensive Compliance).