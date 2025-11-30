# ⚔️ SAAS TRANSFORMATION - IMPLEMENTATION DECISIONS

**Date:** 2024-12-XX  
**Decision Maker:** Zo (with manager guidance)

---

## 📋 MANAGER GUIDANCE

### **Confirmed:**
- ✅ **Authentication:** Use Supabase Auth
- ✅ **Data Migration:** Use best judgment (decided to preserve existing data)
- ✅ **Everything Else:** Use best judgment

### **Uncertain:**
- ❓ **Supabase Auth Status:** Not sure if enabled (will check and guide setup)

---

## 🎯 DECISIONS MADE

### **1. Supabase Auth Setup**
**Decision:** Create setup guide and verification script. If not enabled, guide through enabling it.

**Implementation:**
- Create `scripts/verify_supabase_auth.py` to check Auth status
- Create setup guide in Component 1 README
- If Auth not enabled, provide step-by-step enable instructions

---

### **2. Existing Data Migration Strategy**
**Decision:** **Preserve existing data, create new schema, link gradually**

**Rationale:**
- Existing `analysis_history` table is working - don't break it
- Existing `user_sessions` and `session_items` are working - keep them
- Create new SaaS tables (`user_profiles`, `user_quotas`, etc.)
- Link existing data to new user system when users authenticate
- Provide migration utilities for linking anonymous data to users

**Migration Plan:**
1. **Phase 1:** Create new SaaS schema tables (don't touch existing)
2. **Phase 2:** When user signs up, offer to link existing anonymous data
3. **Phase 3:** Create migration utility for linking `analysis_history` to `saved_analyses`
4. **Phase 4:** Gradually migrate sessions to use authenticated user IDs

**Tables to Keep:**
- `analysis_history` - Keep as-is (frontend analysis history)
- `user_sessions` - Keep as-is, update to use authenticated user IDs
- `session_items` - Keep as-is, update to use authenticated user IDs
- `mdt_runs`, `mdt_events`, etc. - Keep for logging

**New Tables to Create:**
- `user_profiles` - User metadata
- `user_subscriptions` - Subscription management
- `user_quotas` - Quota tracking
- `user_feature_flags` - Feature flags
- `features` - Feature registry
- `saved_analyses` - New structure (can link to `analysis_history`)
- `usage_logs` - Usage tracking

---

### **3. Feature Flags Strategy**
**Decision:** **Hybrid approach - Keep env flags for system config, add user flags for tier gating**

**Rationale:**
- Environment flags (`DISABLE_EVO2`, `DISABLE_LITERATURE`) are system-level configs
- User flags (`sae_features`, `clinical_trials`) are tier-based features
- Both serve different purposes, keep both

**Implementation:**
- Keep existing env flags in `config.py` for system-level control
- Add user-based feature flags in database for tier gating
- Combine both in middleware (system flag OR user flag = enabled)

---

### **4. Authentication Implementation**
**Decision:** **Use Supabase Auth SDK for both backend and frontend**

**Backend:**
- Use `supabase` Python package for Auth operations
- Use JWT verification with Supabase JWT secret
- Create middleware for token verification

**Frontend:**
- Use existing `@supabase/supabase-js` (already installed)
- Create `AuthContext` using Supabase Auth methods
- Use Supabase Auth UI components if needed

---

### **5. Session Management**
**Decision:** **Update existing sessions router to use authenticated users**

**Rationale:**
- Sessions router already exists and works
- Just need to update to use `Depends(get_current_user)` instead of optional header
- Link existing anonymous sessions to users when they sign up

**Implementation:**
- Update `sessions.py` to require authentication
- Add migration utility to link anonymous sessions to user IDs
- Update `AnalysisHistoryContext` to use authenticated user

---

### **6. Analysis History**
**Decision:** **Keep existing `analysis_history`, create new `saved_analyses` for new structure**

**Rationale:**
- `analysis_history` is working and used by frontend
- Create new `saved_analyses` table with better structure for SaaS
- Link `analysis_history` entries to `saved_analyses` when user authenticates
- Gradually migrate frontend to use `saved_analyses`

**Implementation:**
- Keep `analysis_history` table as-is
- Create `saved_analyses` table per SaaS schema
- When user authenticates, offer to link `analysis_history` to `saved_analyses`
- Update `AnalysisHistoryService` to use authenticated user ID

---

### **7. Database Schema Migration**
**Decision:** **Incremental migration - create new tables, keep existing**

**Strategy:**
1. Run SaaS schema SQL to create new tables
2. Don't drop existing tables
3. Add foreign keys where appropriate
4. Create migration utilities for data linking

---

### **8. Protected Endpoints**
**Decision:** **Add auth middleware gradually, starting with premium endpoints**

**Strategy:**
- Phase 1: Add auth to premium endpoints (efficacy, design, clinical trials)
- Phase 2: Add auth to all endpoints
- Phase 3: Add quota checks
- Phase 4: Add feature flag checks

---

### **9. Email Confirmation**
**Decision:** **Require email confirmation for security**

**Rationale:**
- Standard SaaS practice
- Prevents fake accounts
- Supabase Auth handles this automatically

**Implementation:**
- Configure Supabase Auth to require email confirmation
- Show pending confirmation message in frontend
- Allow resend confirmation email

---

### **10. OAuth Providers**
**Decision:** **Start with email/password only, add OAuth later if needed**

**Rationale:**
- Faster initial implementation
- Can add Google/GitHub OAuth in Phase 2
- Email/password is sufficient for MVP

**Implementation:**
- Use Supabase Auth email/password signup/login
- Leave OAuth configuration for future enhancement

---

## 🚀 EXECUTION PLAN

### **Phase 1: Auth Setup (Day 1)**
1. Create Supabase Auth verification script
2. Check if Auth is enabled
3. If not, guide through enabling
4. Get JWT secret
5. Run SaaS schema (new tables only)

### **Phase 2: Auth Implementation (Days 2-3)**
1. Backend auth middleware
2. Backend auth endpoints
3. Frontend auth context
4. Frontend login/signup pages
5. Protected routes

### **Phase 3: Integration (Day 4)**
1. Update sessions router to use auth
2. Update analysis history to use auth
3. Add auth to premium endpoints
4. Test end-to-end

---

## ✅ ACCEPTANCE CRITERIA

- [ ] Supabase Auth is enabled and working
- [ ] Users can sign up with email/password
- [ ] Users can log in and receive JWT token
- [ ] Protected endpoints require authentication
- [ ] Existing data is preserved
- [ ] New SaaS tables are created
- [ ] Migration utilities are available for linking data

---

**These decisions guide the implementation. All components should follow these principles.**








