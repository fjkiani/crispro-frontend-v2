# ⚡ Quick Setup Guide - Auth Component

**5-minute setup guide for Component 1 (Auth)**

---

## 🚀 FASTEST PATH

### **1. Install Dependencies (1 minute)**
```bash
cd oncology-coPilot/oncology-backend-minimal
bash scripts/setup_auth_dependencies.sh
```

Or manually:
```bash
pip install supabase PyJWT python-multipart
```

### **2. Create .env Template (30 seconds)**
```bash
bash scripts/create_env_template.sh
```

### **3. Supabase Setup (3-5 minutes)**
Follow `SUPABASE_SETUP_GUIDE.md`:
1. Enable Email Auth in Supabase Dashboard
2. Get JWT Secret from Settings → API
3. Add to `.env`:
   ```
   SUPABASE_JWT_SECRET=your-actual-jwt-secret
   ```
4. Run database schema in Supabase SQL Editor

### **4. Verify (30 seconds)**
```bash
python3 scripts/verify_supabase_auth.py
```

### **5. Test (1 minute)**
```bash
bash scripts/test_auth_endpoints.sh
```

---

## ✅ DONE!

**Total time: ~5-10 minutes**

**After setup, authentication is fully functional!**

---

## 📋 WHAT'S ALREADY DONE (No Setup Required)

✅ All backend code (auth middleware, service, router)  
✅ All frontend code (AuthContext, login/signup pages)  
✅ Integration with existing endpoints  
✅ Test scripts  
✅ Documentation  

**Nothing else needs to be done outside of Supabase setup!**







