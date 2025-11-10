# 🔧 Supabase Setup Guide for Auth

**This guide walks you through setting up Supabase for authentication.**

---

## 📋 STEP 1: Verify Supabase Project

1. **Go to Supabase Dashboard:** https://app.supabase.com
2. **Select your project** (or create a new one)
3. **Verify project URL:**
   - Go to Settings → API
   - Copy "Project URL" (should be `https://xxxxx.supabase.co`)
   - This is your `SUPABASE_URL`

---

## 📋 STEP 2: Enable Supabase Auth

1. **Go to Authentication → Settings**
2. **Enable "Email Auth":**
   - Toggle "Enable Email Signup" to ON
   - Toggle "Enable Email Signin" to ON
3. **Configure Email Templates (optional):**
   - Go to Authentication → Email Templates
   - Customize signup/login emails if desired
4. **Configure Redirect URLs (optional):**
   - Go to Authentication → URL Configuration
   - Add your frontend URL (e.g., `http://localhost:5173`)

---

## 📋 STEP 3: Get API Keys

1. **Go to Settings → API**
2. **Copy the following:**
   - **Project URL** → `SUPABASE_URL`
   - **anon public** key → `SUPABASE_ANON_KEY`
   - **service_role** key → `SUPABASE_KEY` (optional, for admin operations)
   - **JWT Secret** → `SUPABASE_JWT_SECRET` ⚠️ **REQUIRED FOR AUTH**

3. **Add to `.env` file:**
   ```bash
   SUPABASE_URL=https://xxxxx.supabase.co
   SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...
   SUPABASE_JWT_SECRET=your-jwt-secret-here  # ⚠️ CRITICAL
   ```

---

## 📋 STEP 4: Run Database Schema

1. **Go to Supabase Dashboard → SQL Editor**
2. **Open `.cursor/rules/saas_transformation/schemas/database_schema.sql`**
3. **Copy entire SQL script**
4. **Paste into Supabase SQL Editor**
5. **Click "Run"**
6. **Verify tables created:**
   - `user_profiles`
   - `user_subscriptions`
   - `user_quotas`
   - `user_feature_flags`
   - `features`
   - `saved_analyses`
   - `usage_logs`

---

## 📋 STEP 5: Verify Setup

Run the verification script:

```bash
cd oncology-coPilot/oncology-backend-minimal
python3 scripts/verify_supabase_auth.py
```

**Expected output:**
```
✅ Supabase URL: https://xxxxx.supabase.co...
✅ Supabase Key: eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...
✅ Supabase Auth is accessible
```

---

## 📋 STEP 6: Test Auth Endpoints

Once setup is complete, test the endpoints:

```bash
cd oncology-coPilot/oncology-backend-minimal
bash scripts/test_auth_endpoints.sh
```

---

## 🔍 TROUBLESHOOTING

### **Issue: "Supabase Auth not accessible"**
- **Solution:** Check that Email Auth is enabled in Authentication → Settings
- **Solution:** Verify SUPABASE_URL and SUPABASE_ANON_KEY are correct

### **Issue: "JWT verification failed"**
- **Solution:** Ensure SUPABASE_JWT_SECRET is set correctly
- **Solution:** JWT Secret must match the one in Supabase Dashboard → Settings → API

### **Issue: "Table does not exist"**
- **Solution:** Run database schema SQL in Supabase SQL Editor
- **Solution:** Check table names match exactly (case-sensitive)

### **Issue: "Email confirmation required"**
- **Solution:** This is expected - Supabase Auth requires email confirmation by default
- **Solution:** Check spam folder for confirmation email
- **Solution:** Or disable email confirmation in Authentication → Settings (for testing only)

---

## ✅ CHECKLIST

- [ ] Supabase project created
- [ ] Email Auth enabled
- [ ] SUPABASE_URL added to .env
- [ ] SUPABASE_ANON_KEY added to .env
- [ ] SUPABASE_JWT_SECRET added to .env
- [ ] Database schema run in Supabase SQL Editor
- [ ] Verification script passes
- [ ] Test endpoints work

---

**Once all steps are complete, authentication will be fully functional!**




