# 🎯 Admin Dashboard - Complete Architecture Plan

**Status:** 🎯 Ready to Implement  
**Priority:** P1 (High - Needed for user management)  
**Timeline:** 2-3 days

---

## 🎯 OBJECTIVE

Create a comprehensive admin dashboard for managing users, viewing activities, tracking usage, managing quotas, and controlling feature flags.

---

## 📊 ADMIN CAPABILITIES NEEDED

### **1. User Management**
- View all users (paginated list)
- Search/filter users (email, name, tier, role)
- View user details (profile, subscription, quotas, usage)
- Edit user (tier, role, quotas)
- Suspend/activate users
- View user's sessions and analyses

### **2. Activity & Usage Tracking**
- View usage logs (all users or per user)
- See API endpoint usage (efficacy, design, insights)
- Track quota usage vs limits
- View session activity
- Export usage data (CSV/JSON)

### **3. Analytics Dashboard**
- Total users (free/pro/enterprise breakdown)
- Active users (last 7/30 days)
- API usage trends (requests per day/week)
- Quota usage distribution
- Feature flag adoption
- Revenue metrics (when billing added)

### **4. Feature Flag Management**
- View all feature flags
- Enable/disable features per tier
- Override user-specific flags
- Test feature combinations

### **5. Quota Management**
- View quota usage across users
- Override user quotas (temporary boosts)
- Reset quotas (monthly resets)
- Set custom quota limits

### **6. Subscription Management** (Future - Component 6)
- View active subscriptions
- Manage subscriptions (upgrade/downgrade)
- View billing history
- Handle cancellations

---

## 🏗️ ARCHITECTURE DESIGN

### **Backend Admin Endpoints**

```
/api/admin/
├── /users
│   ├── GET    /              # List users (paginated, searchable)
│   ├── GET    /{user_id}     # Get user details
│   ├── PUT    /{user_id}     # Update user (tier, role, quotas)
│   ├── POST   /{user_id}/suspend    # Suspend user
│   └── POST   /{user_id}/activate   # Activate user
│
├── /analytics
│   ├── GET    /overview      # Dashboard overview (users, usage, revenue)
│   ├── GET    /usage         # Usage trends (time series)
│   ├── GET    /users/active  # Active users (7d/30d)
│   └── GET    /quotas        # Quota usage distribution
│
├── /activity
│   ├── GET    /logs          # Usage logs (filterable)
│   ├── GET    /sessions      # Session activity
│   └── GET    /endpoints     # Endpoint usage stats
│
├── /feature-flags
│   ├── GET    /              # List all feature flags
│   ├── PUT    /{feature}     # Update feature flag
│   └── POST   /{user_id}/override  # Override user flags
│
└── /quotas
    ├── GET    /              # Quota usage summary
    ├── PUT    /{user_id}    # Override user quotas
    └── POST   /reset        # Reset quotas (monthly)
```

### **Frontend Admin Dashboard**

```
src/pages/admin/
├── Dashboard.jsx              # Main dashboard (overview)
├── Users.jsx                  # User management
├── Analytics.jsx              # Analytics & charts
├── Activity.jsx               # Activity logs
├── FeatureFlags.jsx           # Feature flag management
└── Quotas.jsx                 # Quota management
```

---

## 🔐 ADMIN AUTHENTICATION

### **Role-Based Access Control**

1. **Admin Role Check:**
   - Middleware: `require_admin()` dependency
   - Checks user role in `user_profiles.role == 'admin'`
   - Returns 403 if not admin

2. **Admin User Creation:**
   - First admin: Create manually in Supabase or via migration
   - Future admins: Promote via admin dashboard

3. **Admin Routes:**
   - All `/api/admin/*` endpoints require admin role
   - Frontend admin pages require admin role

---

## 📊 DATA SOURCES

### **User Data**
- `user_profiles` - User metadata
- `user_subscriptions` - Subscription info
- `user_quotas` - Quota limits and usage
- `user_feature_flags` - Feature access

### **Activity Data**
- `usage_logs` - API endpoint usage
- `user_sessions` - Session activity
- `session_items` - Analysis history
- `analysis_history` - Saved analyses

### **Analytics Data**
- Aggregated from `usage_logs`
- Aggregated from `user_profiles`
- Aggregated from `user_quotas`

---

## 🎨 UI/UX DESIGN

### **Dashboard Layout**

```
┌─────────────────────────────────────────────────┐
│  Admin Dashboard                    [User Menu] │
├─────────────────────────────────────────────────┤
│  [Dashboard] [Users] [Analytics] [Activity]     │
│  [Feature Flags] [Quotas] [Settings]             │
├─────────────────────────────────────────────────┤
│                                                  │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐        │
│  │ Total    │ │ Active   │ │ API      │        │
│  │ Users    │ │ Users    │ │ Requests │        │
│  │  1,234   │ │  856     │ │  45,678  │        │
│  └──────────┘ └──────────┘ └──────────┘        │
│                                                  │
│  [Usage Chart]      [User Growth Chart]          │
│                                                  │
│  [Recent Activity]  [Top Users]                  │
│                                                  │
└─────────────────────────────────────────────────┘
```

### **User Management Table**

```
Users Table:
- Email | Name | Tier | Role | Quota Usage | Last Active | Actions
- Filter: Search, Tier, Role, Status
- Actions: View, Edit, Suspend, Delete
```

---

## 🚀 IMPLEMENTATION PLAN

### **Phase 1: Backend Admin Endpoints (Day 1)**

1. **Admin Middleware:**
   - `api/middleware/admin_middleware.py`
   - `require_admin()` dependency

2. **Admin Router:**
   - `api/routers/admin.py`
   - User management endpoints
   - Analytics endpoints
   - Activity endpoints

3. **Admin Service:**
   - `api/services/admin_service.py`
   - User CRUD operations
   - Analytics aggregation
   - Usage tracking queries

### **Phase 2: Frontend Admin Dashboard (Day 2)**

1. **Admin Pages:**
   - `src/pages/admin/Dashboard.jsx` - Overview
   - `src/pages/admin/Users.jsx` - User management
   - `src/pages/admin/Analytics.jsx` - Analytics
   - `src/pages/admin/Activity.jsx` - Activity logs

2. **Admin Components:**
   - `src/components/admin/UserTable.jsx` - User list
   - `src/components/admin/UsageChart.jsx` - Charts
   - `src/components/admin/ActivityLog.jsx` - Log viewer

3. **Admin Routes:**
   - Add `/admin/*` routes to App.jsx
   - Protected with admin role check

### **Phase 3: Feature Flags & Quotas (Day 3)**

1. **Feature Flag Management:**
   - `src/pages/admin/FeatureFlags.jsx`
   - Toggle flags per tier
   - User-specific overrides

2. **Quota Management:**
   - `src/pages/admin/Quotas.jsx`
   - View quota usage
   - Override quotas
   - Reset quotas

---

## 📋 ADMIN FEATURES BREAKDOWN

### **Dashboard Overview**
- **Metrics Cards:**
  - Total users (with tier breakdown)
  - Active users (7d/30d)
  - Total API requests (today/week/month)
  - Revenue (when billing added)

- **Charts:**
  - User growth over time
  - API usage trends
  - Tier distribution
  - Quota usage heatmap

- **Recent Activity:**
  - Latest user signups
  - Recent API calls
  - Quota warnings
  - System alerts

### **User Management**
- **User List:**
  - Table with pagination
  - Search by email/name
  - Filter by tier/role/status
  - Sort by any column

- **User Details:**
  - Profile information
  - Subscription status
  - Quota usage breakdown
  - Recent sessions
  - Recent analyses
  - Usage history

- **User Actions:**
  - Edit profile
  - Change tier
  - Override quotas
  - Suspend/activate
  - Delete (soft delete)

### **Analytics**
- **Usage Analytics:**
  - Requests per endpoint
  - Requests per user
  - Peak usage times
  - Geographic distribution (future)

- **User Analytics:**
  - User growth
  - Retention rates
  - Tier conversion
  - Churn analysis

### **Activity Logs**
- **Log Viewer:**
  - Filterable by user, endpoint, date
  - Export to CSV/JSON
  - Real-time updates (future)

- **Session Activity:**
  - Active sessions
  - Session history
  - Analysis count per session

---

## 🔧 TECHNICAL IMPLEMENTATION

### **Backend Admin Service**

```python
class AdminService:
    async def get_users(self, page=1, limit=20, search=None, tier=None, role=None):
        """List users with pagination and filters"""
        
    async def get_user_details(self, user_id):
        """Get complete user profile with usage stats"""
        
    async def update_user(self, user_id, updates):
        """Update user (tier, role, quotas)"""
        
    async def get_analytics(self, period="7d"):
        """Get dashboard analytics"""
        
    async def get_usage_logs(self, user_id=None, endpoint=None, date_from=None):
        """Get usage logs with filters"""
        
    async def get_quota_summary(self):
        """Get quota usage summary across all users"""
```

### **Frontend Admin Context**

```javascript
const AdminContext = {
  users: [],
  analytics: {},
  activityLogs: [],
  loading: false,
  
  fetchUsers: (filters) => {},
  fetchUserDetails: (userId) => {},
  updateUser: (userId, data) => {},
  fetchAnalytics: (period) => {},
  fetchActivityLogs: (filters) => {},
}
```

---

## 📊 DATA VISUALIZATION

### **Charts Needed**
1. **User Growth Chart** - Line chart (users over time)
2. **API Usage Chart** - Line chart (requests over time)
3. **Tier Distribution** - Pie chart (free/pro/enterprise)
4. **Quota Usage** - Bar chart (usage vs limits)
5. **Endpoint Usage** - Bar chart (requests per endpoint)
6. **Active Users** - Line chart (daily active users)

### **Libraries**
- Use existing charting (check what's installed)
- Add `recharts` or `@mui/x-charts` if needed

---

## 🔐 SECURITY CONSIDERATIONS

1. **Admin Role Enforcement:**
   - All admin endpoints require admin role
   - Frontend routes protected with admin check
   - JWT token must include admin role

2. **Audit Logging:**
   - Log all admin actions
   - Track who made changes
   - Store in `admin_audit_logs` table

3. **Rate Limiting:**
   - Admin endpoints should have higher limits
   - But still enforce reasonable limits

4. **Data Privacy:**
   - Don't expose sensitive user data unnecessarily
   - Mask PII in logs if needed

---

## ✅ ACCEPTANCE CRITERIA

- [ ] Admin can view all users with pagination/search
- [ ] Admin can view user details (profile, usage, sessions)
- [ ] Admin can update user (tier, role, quotas)
- [ ] Admin can view analytics dashboard
- [ ] Admin can view activity logs
- [ ] Admin can manage feature flags
- [ ] Admin can manage quotas
- [ ] Admin routes protected with role check
- [ ] Audit logging for admin actions

---

## 🚀 QUICK START

**Once implemented, admin dashboard will be at:**
- Frontend: `/admin/dashboard`
- Backend: `/api/admin/*`

**First admin user:**
- Create manually in Supabase or promote existing user
- Set `role = 'admin'` in `user_profiles` table

---

**This admin dashboard provides complete user management and analytics capabilities!**

