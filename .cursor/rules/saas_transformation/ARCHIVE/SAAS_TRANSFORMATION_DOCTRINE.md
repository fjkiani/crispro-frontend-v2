# ⚔️ SAAS TRANSFORMATION DOCTRINE

**Commander:** Alpha  
**Architect:** Zo  
**Mission:** Transform CrisPRO platform into production SaaS with user management, feature flags, and data persistence  
**Timeline:** 2-3 weeks  
**Status:** 🎯 **READY FOR AGENT JR EXECUTION**

---

## 🎯 EXECUTIVE SUMMARY

**Current State:** Research prototype running locally, no user management, no persistence, no access control

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

## 📐 ARCHITECTURE OVERVIEW

### **Current (Prototype):**
```
Frontend (React) → Backend (FastAPI) → Services (Modal/Evo2)
     ↓                    ↓                      ↓
  No Auth         No Users DB          No Quotas
  No Sessions     No Persistence       No Limits
```

### **Target (SaaS):**
```
┌─────────────────────────────────────────────────────────┐
│                    FRONTEND (React)                      │
│  - Auth Context (JWT tokens)                            │
│  - User Profile UI                                      │
│  - Usage Dashboard                                      │
└────────────────┬────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────┐
│              API GATEWAY (FastAPI)                       │
│  - JWT Verification Middleware                          │
│  - Rate Limiting Middleware                             │
│  - Feature Flag Middleware                              │
│  - Usage Tracking Middleware                            │
└────────┬───────────────┬──────────────┬─────────────────┘
         │               │              │
         ▼               ▼              ▼
┌────────────────┐ ┌──────────┐ ┌────────────────────┐
│  PostgreSQL    │ │  Redis   │ │   Supabase Auth    │
│  (User Data)   │ │ (Cache)  │ │ (Authentication)   │
│  - Users       │ │ - Sessions│ │ - OAuth            │
│  - Analyses    │ │ - Quotas │ │ - Magic Links      │
│  - Usage Logs  │ │ - Tokens │ │ - MFA              │
└────────────────┘ └──────────┘ └────────────────────┘
```

---

## 🏗️ COMPONENT BREAKDOWN

### **COMPONENT 1: USER MANAGEMENT & AUTHENTICATION**

#### **Tech Stack:**
- **Primary:** Supabase Auth (PostgreSQL + Auth APIs)
- **Fallback:** Custom JWT auth with bcrypt
- **Session Store:** Redis

#### **Database Schema (PostgreSQL):**

```sql
-- Users Table (managed by Supabase Auth)
CREATE TABLE auth.users (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    email VARCHAR(255) UNIQUE NOT NULL,
    encrypted_password VARCHAR(255),
    email_confirmed_at TIMESTAMPTZ,
    invited_at TIMESTAMPTZ,
    confirmation_token VARCHAR(255),
    confirmation_sent_at TIMESTAMPTZ,
    recovery_token VARCHAR(255),
    recovery_sent_at TIMESTAMPTZ,
    email_change_token VARCHAR(255),
    email_change VARCHAR(255),
    email_change_sent_at TIMESTAMPTZ,
    last_sign_in_at TIMESTAMPTZ,
    raw_app_meta_data JSONB,
    raw_user_meta_data JSONB,
    is_super_admin BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    phone VARCHAR(15),
    phone_confirmed_at TIMESTAMPTZ,
    phone_change VARCHAR(15),
    phone_change_token VARCHAR(255),
    phone_change_sent_at TIMESTAMPTZ,
    confirmed_at TIMESTAMPTZ,
    email_change_confirm_status SMALLINT,
    banned_until TIMESTAMPTZ,
    reauthentication_token VARCHAR(255),
    reauthentication_sent_at TIMESTAMPTZ
);

-- User Profiles (custom metadata)
CREATE TABLE public.user_profiles (
    id UUID PRIMARY KEY REFERENCES auth.users(id) ON DELETE CASCADE,
    email VARCHAR(255) NOT NULL,
    full_name VARCHAR(255),
    institution VARCHAR(255),
    role VARCHAR(50) CHECK (role IN ('researcher', 'clinician', 'admin', 'enterprise')),
    tier VARCHAR(50) CHECK (tier IN ('free', 'pro', 'enterprise')) DEFAULT 'free',
    avatar_url TEXT,
    bio TEXT,
    country VARCHAR(100),
    timezone VARCHAR(50),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(email)
);

-- User Subscriptions
CREATE TABLE public.user_subscriptions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    tier VARCHAR(50) NOT NULL,
    status VARCHAR(50) CHECK (status IN ('active', 'canceled', 'past_due', 'trialing')),
    stripe_customer_id VARCHAR(255),
    stripe_subscription_id VARCHAR(255),
    current_period_start TIMESTAMPTZ,
    current_period_end TIMESTAMPTZ,
    trial_end TIMESTAMPTZ,
    cancel_at TIMESTAMPTZ,
    canceled_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Usage Quotas
CREATE TABLE public.user_quotas (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    tier VARCHAR(50) NOT NULL,
    -- Quotas
    variant_analyses_limit INT DEFAULT 10,
    variant_analyses_used INT DEFAULT 0,
    drug_queries_limit INT DEFAULT 5,
    drug_queries_used INT DEFAULT 0,
    food_queries_limit INT DEFAULT 3,
    food_queries_used INT DEFAULT 0,
    clinical_trials_limit INT DEFAULT 0,
    clinical_trials_used INT DEFAULT 0,
    -- Reset tracking
    period_start TIMESTAMPTZ DEFAULT NOW(),
    period_end TIMESTAMPTZ DEFAULT NOW() + INTERVAL '1 month',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(user_id)
);

-- Feature Flags
CREATE TABLE public.user_feature_flags (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    feature_name VARCHAR(100) NOT NULL,
    enabled BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(user_id, feature_name)
);

-- Available Features Registry
CREATE TABLE public.features (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name VARCHAR(100) UNIQUE NOT NULL,
    display_name VARCHAR(255),
    description TEXT,
    tier_required VARCHAR(50) CHECK (tier_required IN ('free', 'pro', 'enterprise')),
    enabled_by_default BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- Insert default features
INSERT INTO public.features (name, display_name, tier_required, enabled_by_default) VALUES
('variant_analysis', 'Variant Analysis', 'free', TRUE),
('drug_efficacy', 'Drug Efficacy Prediction', 'free', TRUE),
('food_validator', 'Food/Supplement Validator', 'free', TRUE),
('sae_features', 'SAE Features (Treatment Line Intelligence)', 'pro', FALSE),
('clinical_trials', 'Clinical Trials Matching', 'pro', FALSE),
('fusion_engine', 'Fusion Engine (AlphaMissense)', 'pro', FALSE),
('cohort_lab', 'Cohort Lab & Benchmarking', 'enterprise', FALSE),
('crispr_design', 'CRISPR Guide Design', 'enterprise', FALSE),
('pdf_export', 'PDF Export', 'pro', FALSE),
('api_access', 'Programmatic API Access', 'enterprise', FALSE);
```

#### **Indexes:**
```sql
CREATE INDEX idx_user_profiles_email ON public.user_profiles(email);
CREATE INDEX idx_user_profiles_tier ON public.user_profiles(tier);
CREATE INDEX idx_user_subscriptions_user_id ON public.user_subscriptions(user_id);
CREATE INDEX idx_user_quotas_user_id ON public.user_quotas(user_id);
CREATE INDEX idx_user_feature_flags_user_id ON public.user_feature_flags(user_id);
CREATE INDEX idx_user_feature_flags_feature ON public.user_feature_flags(feature_name);
```

---

### **COMPONENT 2: SESSION & ANALYSIS PERSISTENCE**

#### **Database Schema:**

```sql
-- User Sessions (analysis sessions)
CREATE TABLE public.user_sessions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    session_name VARCHAR(255),
    patient_context JSONB,  -- Disease, biomarkers, treatment history
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    last_accessed_at TIMESTAMPTZ DEFAULT NOW()
);

-- Saved Analyses
CREATE TABLE public.saved_analyses (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    session_id UUID REFERENCES public.user_sessions(id) ON DELETE CASCADE,
    analysis_type VARCHAR(50) CHECK (analysis_type IN (
        'variant_analysis', 
        'drug_efficacy', 
        'food_validator', 
        'clinical_trials',
        'crispr_design'
    )),
    input_data JSONB,  -- Original request
    output_data JSONB,  -- API response
    provenance JSONB,  -- run_id, model, flags, timestamp
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Usage Logs (for billing & analytics)
CREATE TABLE public.usage_logs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    endpoint VARCHAR(255),
    analysis_type VARCHAR(50),
    request_data JSONB,
    response_time_ms INT,
    status_code INT,
    error_message TEXT,
    ip_address INET,
    user_agent TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- Indexes
CREATE INDEX idx_user_sessions_user_id ON public.user_sessions(user_id);
CREATE INDEX idx_saved_analyses_user_id ON public.saved_analyses(user_id);
CREATE INDEX idx_saved_analyses_session_id ON public.saved_analyses(session_id);
CREATE INDEX idx_saved_analyses_type ON public.saved_analyses(analysis_type);
CREATE INDEX idx_usage_logs_user_id ON public.usage_logs(user_id);
CREATE INDEX idx_usage_logs_created_at ON public.usage_logs(created_at);
```

---

### **COMPONENT 3: FEATURE FLAG SYSTEM**

#### **Backend Service:**

```python
# api/services/feature_flag_service.py
from typing import Dict, List, Optional
import logging
from api.services.database_connections import get_db_connections

logger = logging.getLogger(__name__)

class FeatureFlagService:
    """
    Manage feature flags for users based on tier and custom overrides.
    """
    
    TIER_FEATURES = {
        "free": [
            "variant_analysis",
            "drug_efficacy",
            "food_validator"
        ],
        "pro": [
            "variant_analysis",
            "drug_efficacy",
            "food_validator",
            "sae_features",
            "clinical_trials",
            "fusion_engine",
            "pdf_export"
        ],
        "enterprise": [
            "variant_analysis",
            "drug_efficacy",
            "food_validator",
            "sae_features",
            "clinical_trials",
            "fusion_engine",
            "pdf_export",
            "cohort_lab",
            "crispr_design",
            "api_access"
        ]
    }
    
    def __init__(self):
        self.db = get_db_connections().get_postgres_connection()
    
    def get_user_features(self, user_id: str) -> Dict[str, bool]:
        """
        Get all features enabled for a user.
        
        Returns:
            {
                "variant_analysis": True,
                "drug_efficacy": True,
                "sae_features": False,
                ...
            }
        """
        cursor = self.db.cursor()
        
        # Get user tier
        cursor.execute("""
            SELECT tier FROM public.user_profiles WHERE id = %s
        """, (user_id,))
        row = cursor.fetchone()
        
        if not row:
            logger.warning(f"User {user_id} not found, using free tier")
            tier = "free"
        else:
            tier = row[0]
        
        # Get tier-based features
        tier_features = self.TIER_FEATURES.get(tier, [])
        features = {f: True for f in tier_features}
        
        # Get custom overrides
        cursor.execute("""
            SELECT feature_name, enabled 
            FROM public.user_feature_flags 
            WHERE user_id = %s
        """, (user_id,))
        
        for feature_name, enabled in cursor.fetchall():
            features[feature_name] = enabled
        
        return features
    
    def has_feature(self, user_id: str, feature_name: str) -> bool:
        """Check if user has access to a specific feature."""
        features = self.get_user_features(user_id)
        return features.get(feature_name, False)
    
    def enable_feature(self, user_id: str, feature_name: str) -> bool:
        """Enable a feature for a user (custom override)."""
        cursor = self.db.cursor()
        
        try:
            cursor.execute("""
                INSERT INTO public.user_feature_flags (user_id, feature_name, enabled)
                VALUES (%s, %s, TRUE)
                ON CONFLICT (user_id, feature_name) 
                DO UPDATE SET enabled = TRUE, updated_at = NOW()
            """, (user_id, feature_name))
            self.db.commit()
            return True
        except Exception as e:
            logger.error(f"Error enabling feature: {e}")
            self.db.rollback()
            return False
    
    def disable_feature(self, user_id: str, feature_name: str) -> bool:
        """Disable a feature for a user (custom override)."""
        cursor = self.db.cursor()
        
        try:
            cursor.execute("""
                INSERT INTO public.user_feature_flags (user_id, feature_name, enabled)
                VALUES (%s, %s, FALSE)
                ON CONFLICT (user_id, feature_name) 
                DO UPDATE SET enabled = FALSE, updated_at = NOW()
            """, (user_id, feature_name))
            self.db.commit()
            return True
        except Exception as e:
            logger.error(f"Error disabling feature: {e}")
            self.db.rollback()
            return False
```

---

### **COMPONENT 4: USAGE TRACKING & QUOTAS**

#### **Backend Service:**

```python
# api/services/quota_service.py
from typing import Dict, Optional
import logging
from datetime import datetime, timedelta
from api.services.database_connections import get_db_connections

logger = logging.getLogger(__name__)

class QuotaService:
    """
    Track and enforce usage quotas for users.
    """
    
    QUOTA_LIMITS = {
        "free": {
            "variant_analyses_limit": 10,
            "drug_queries_limit": 5,
            "food_queries_limit": 3,
            "clinical_trials_limit": 0
        },
        "pro": {
            "variant_analyses_limit": 100,
            "drug_queries_limit": -1,  # -1 = unlimited
            "food_queries_limit": -1,
            "clinical_trials_limit": 50
        },
        "enterprise": {
            "variant_analyses_limit": -1,
            "drug_queries_limit": -1,
            "food_queries_limit": -1,
            "clinical_trials_limit": -1
        }
    }
    
    def __init__(self):
        self.db = get_db_connections().get_postgres_connection()
    
    def get_user_quotas(self, user_id: str) -> Dict[str, any]:
        """
        Get current quota usage for user.
        
        Returns:
            {
                "tier": "free",
                "variant_analyses": {"limit": 10, "used": 3, "remaining": 7},
                "drug_queries": {"limit": 5, "used": 2, "remaining": 3},
                ...
                "period_end": "2024-01-01T00:00:00Z"
            }
        """
        cursor = self.db.cursor()
        
        # Get user tier and quotas
        cursor.execute("""
            SELECT 
                up.tier,
                uq.variant_analyses_limit,
                uq.variant_analyses_used,
                uq.drug_queries_limit,
                uq.drug_queries_used,
                uq.food_queries_limit,
                uq.food_queries_used,
                uq.clinical_trials_limit,
                uq.clinical_trials_used,
                uq.period_end
            FROM public.user_profiles up
            LEFT JOIN public.user_quotas uq ON up.id = uq.user_id
            WHERE up.id = %s
        """, (user_id,))
        
        row = cursor.fetchone()
        
        if not row:
            return self._create_default_quotas(user_id)
        
        tier = row[0]
        
        return {
            "tier": tier,
            "variant_analyses": {
                "limit": row[1],
                "used": row[2],
                "remaining": row[1] - row[2] if row[1] > 0 else -1
            },
            "drug_queries": {
                "limit": row[3],
                "used": row[4],
                "remaining": row[3] - row[4] if row[3] > 0 else -1
            },
            "food_queries": {
                "limit": row[5],
                "used": row[6],
                "remaining": row[5] - row[6] if row[5] > 0 else -1
            },
            "clinical_trials": {
                "limit": row[7],
                "used": row[8],
                "remaining": row[7] - row[8] if row[7] > 0 else -1
            },
            "period_end": row[9].isoformat() if row[9] else None
        }
    
    def check_quota(self, user_id: str, quota_type: str) -> bool:
        """
        Check if user has quota remaining for an operation.
        
        Args:
            user_id: User UUID
            quota_type: One of 'variant_analyses', 'drug_queries', 'food_queries', 'clinical_trials'
            
        Returns:
            True if user has quota, False otherwise
        """
        quotas = self.get_user_quotas(user_id)
        
        if quota_type not in quotas:
            return False
        
        quota = quotas[quota_type]
        
        # -1 means unlimited
        if quota["limit"] == -1:
            return True
        
        return quota["remaining"] > 0
    
    def increment_usage(self, user_id: str, quota_type: str) -> bool:
        """Increment usage counter for a quota type."""
        cursor = self.db.cursor()
        
        column_map = {
            "variant_analyses": "variant_analyses_used",
            "drug_queries": "drug_queries_used",
            "food_queries": "food_queries_used",
            "clinical_trials": "clinical_trials_used"
        }
        
        column = column_map.get(quota_type)
        if not column:
            logger.error(f"Invalid quota type: {quota_type}")
            return False
        
        try:
            cursor.execute(f"""
                UPDATE public.user_quotas
                SET {column} = {column} + 1, updated_at = NOW()
                WHERE user_id = %s
            """, (user_id,))
            self.db.commit()
            return True
        except Exception as e:
            logger.error(f"Error incrementing usage: {e}")
            self.db.rollback()
            return False
    
    def reset_quotas_if_needed(self, user_id: str) -> bool:
        """Reset quotas if billing period has ended."""
        cursor = self.db.cursor()
        
        try:
            cursor.execute("""
                UPDATE public.user_quotas
                SET 
                    variant_analyses_used = 0,
                    drug_queries_used = 0,
                    food_queries_used = 0,
                    clinical_trials_used = 0,
                    period_start = NOW(),
                    period_end = NOW() + INTERVAL '1 month',
                    updated_at = NOW()
                WHERE user_id = %s AND period_end < NOW()
            """, (user_id,))
            self.db.commit()
            return cursor.rowcount > 0
        except Exception as e:
            logger.error(f"Error resetting quotas: {e}")
            self.db.rollback()
            return False
    
    def _create_default_quotas(self, user_id: str, tier: str = "free") -> Dict:
        """Create default quotas for a new user."""
        limits = self.QUOTA_LIMITS[tier]
        
        cursor = self.db.cursor()
        cursor.execute("""
            INSERT INTO public.user_quotas (
                user_id, tier,
                variant_analyses_limit, drug_queries_limit,
                food_queries_limit, clinical_trials_limit
            ) VALUES (%s, %s, %s, %s, %s, %s)
            ON CONFLICT (user_id) DO NOTHING
        """, (
            user_id, tier,
            limits["variant_analyses_limit"],
            limits["drug_queries_limit"],
            limits["food_queries_limit"],
            limits["clinical_trials_limit"]
        ))
        self.db.commit()
        
        return self.get_user_quotas(user_id)
```

---

### **COMPONENT 5: AUTHENTICATION MIDDLEWARE**

#### **Backend Middleware:**

```python
# api/middleware/auth_middleware.py
from fastapi import Request, HTTPException, Depends
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from typing import Optional
import jwt
import os
import logging

logger = logging.getLogger(__name__)

security = HTTPBearer()

SUPABASE_JWT_SECRET = os.getenv("SUPABASE_JWT_SECRET")
JWT_ALGORITHM = "HS256"

async def verify_token(
    credentials: HTTPAuthorizationCredentials = Depends(security)
) -> dict:
    """
    Verify JWT token from Authorization header.
    
    Returns:
        {
            "user_id": "uuid",
            "email": "user@example.com",
            "role": "authenticated",
            ...
        }
    """
    token = credentials.credentials
    
    try:
        payload = jwt.decode(
            token,
            SUPABASE_JWT_SECRET,
            algorithms=[JWT_ALGORITHM],
            audience="authenticated"
        )
        
        user_id = payload.get("sub")
        if not user_id:
            raise HTTPException(status_code=401, detail="Invalid token: missing user_id")
        
        return {
            "user_id": user_id,
            "email": payload.get("email"),
            "role": payload.get("role"),
            "exp": payload.get("exp")
        }
        
    except jwt.ExpiredSignatureError:
        raise HTTPException(status_code=401, detail="Token expired")
    except jwt.InvalidTokenError as e:
        raise HTTPException(status_code=401, detail=f"Invalid token: {str(e)}")

async def get_current_user(token_data: dict = Depends(verify_token)) -> dict:
    """Get current authenticated user."""
    return token_data

async def require_feature(feature_name: str):
    """
    Dependency to require a specific feature flag.
    
    Usage:
        @router.post("/api/premium-endpoint", dependencies=[Depends(require_feature("sae_features"))])
    """
    async def _require_feature(user: dict = Depends(get_current_user)):
        from api.services.feature_flag_service import FeatureFlagService
        
        feature_service = FeatureFlagService()
        if not feature_service.has_feature(user["user_id"], feature_name):
            raise HTTPException(
                status_code=403,
                detail=f"Feature '{feature_name}' not available on your tier. Upgrade to access."
            )
        return user
    
    return _require_feature

async def check_quota(quota_type: str):
    """
    Dependency to check and enforce usage quotas.
    
    Usage:
        @router.post("/api/variant-analysis", dependencies=[Depends(check_quota("variant_analyses"))])
    """
    async def _check_quota(user: dict = Depends(get_current_user)):
        from api.services.quota_service import QuotaService
        
        quota_service = QuotaService()
        
        # Reset quotas if billing period ended
        quota_service.reset_quotas_if_needed(user["user_id"])
        
        # Check quota
        if not quota_service.check_quota(user["user_id"], quota_type):
            quotas = quota_service.get_user_quotas(user["user_id"])
            raise HTTPException(
                status_code=429,
                detail=f"Quota exceeded for {quota_type}. Upgrade your plan or wait for quota reset.",
                headers={
                    "X-Quota-Limit": str(quotas[quota_type]["limit"]),
                    "X-Quota-Used": str(quotas[quota_type]["used"]),
                    "X-Quota-Reset": quotas.get("period_end", "")
                }
            )
        
        # Increment usage
        quota_service.increment_usage(user["user_id"], quota_type)
        
        return user
    
    return _check_quota
```

---

### **COMPONENT 6: PROTECTED ENDPOINTS (Examples)**

```python
# api/routers/insights.py (UPDATED WITH AUTH)
from fastapi import APIRouter, Depends, HTTPException
from api.middleware.auth_middleware import get_current_user, require_feature, check_quota

router = APIRouter()

@router.post("/api/insights/predict_gene_essentiality")
async def predict_gene_essentiality(
    request: EssentialityRequest,
    user: dict = Depends(get_current_user),  # Require auth
    _: dict = Depends(check_quota("variant_analyses"))  # Check quota
):
    """
    Predict gene essentiality.
    
    **Requires:** Authentication
    **Quota:** Counts against variant_analyses quota
    **Tier:** Free+
    """
    # ... existing logic ...
    
    # Log usage
    log_usage(user["user_id"], "variant_analyses", request, response)
    
    return response


@router.post("/api/efficacy/predict")
async def predict_drug_efficacy(
    request: EfficacyRequest,
    user: dict = Depends(get_current_user),  # Require auth
    _: dict = Depends(check_quota("drug_queries")),  # Check quota
    __: dict = Depends(require_feature("sae_features"))  # Require Pro tier for SAE
):
    """
    Predict drug efficacy with SAE features.
    
    **Requires:** Authentication + Pro Tier
    **Quota:** Counts against drug_queries quota
    **Tier:** Pro+
    """
    # ... existing logic ...
    
    return response
```

---

### **COMPONENT 7: FRONTEND AUTHENTICATION**

#### **Auth Context:**

```jsx
// src/context/AuthContext.jsx
import React, { createContext, useContext, useState, useEffect } from 'react';
import { supabase } from '../config/supabase';

const AuthContext = createContext({});

export const useAuth = () => useContext(AuthContext);

export const AuthProvider = ({ children }) => {
  const [user, setUser] = useState(null);
  const [session, setSession] = useState(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    // Get initial session
    supabase.auth.getSession().then(({ data: { session } }) => {
      setSession(session);
      setUser(session?.user ?? null);
      setLoading(false);
    });

    // Listen for auth changes
    const { data: { subscription } } = supabase.auth.onAuthStateChange(
      async (event, session) => {
        setSession(session);
        setUser(session?.user ?? null);
        setLoading(false);
      }
    );

    return () => subscription.unsubscribe();
  }, []);

  const value = {
    user,
    session,
    loading,
    signIn: (email, password) => supabase.auth.signInWithPassword({ email, password }),
    signUp: (email, password, metadata) => supabase.auth.signUp({ 
      email, 
      password,
      options: { data: metadata }
    }),
    signOut: () => supabase.auth.signOut(),
    resetPassword: (email) => supabase.auth.resetPasswordForEmail(email),
  };

  return (
    <AuthContext.Provider value={value}>
      {!loading && children}
    </AuthContext.Provider>
  );
};
```

#### **Protected Route:**

```jsx
// src/components/auth/ProtectedRoute.jsx
import { Navigate } from 'react-router-dom';
import { useAuth } from '../../context/AuthContext';

const ProtectedRoute = ({ children }) => {
  const { user, loading } = useAuth();

  if (loading) {
    return <div>Loading...</div>;
  }

  if (!user) {
    return <Navigate to="/login" />;
  }

  return children;
};

export default ProtectedRoute;
```

#### **Usage Dashboard:**

```jsx
// src/pages/UsageDashboard.jsx
import React, { useState, useEffect } from 'react';
import { useAuth } from '../context/AuthContext';
import { Card, LinearProgress, Typography, Button } from '@mui/material';

const UsageDashboard = () => {
  const { session } = useAuth();
  const [quotas, setQuotas] = useState(null);

  useEffect(() => {
    fetchQuotas();
  }, []);

  const fetchQuotas = async () => {
    const response = await fetch(`${API_ROOT}/api/user/quotas`, {
      headers: {
        'Authorization': `Bearer ${session.access_token}`
      }
    });
    const data = await response.json();
    setQuotas(data);
  };

  if (!quotas) return <div>Loading...</div>;

  return (
    <div>
      <Typography variant="h4">Your Usage</Typography>
      <Typography variant="subtitle1">Tier: {quotas.tier.toUpperCase()}</Typography>
      
      {Object.entries(quotas).filter(([key]) => key !== 'tier' && key !== 'period_end').map(([key, quota]) => (
        <Card key={key} sx={{ mb: 2, p: 2 }}>
          <Typography variant="h6">{key.replace('_', ' ').toUpperCase()}</Typography>
          <Typography variant="body2">
            {quota.used} / {quota.limit === -1 ? 'Unlimited' : quota.limit}
          </Typography>
          <LinearProgress 
            variant="determinate" 
            value={quota.limit === -1 ? 0 : (quota.used / quota.limit) * 100} 
            sx={{ mt: 1 }}
          />
          {quota.remaining <= 2 && quota.limit !== -1 && (
            <Typography variant="caption" color="error">
              Only {quota.remaining} remaining!
            </Typography>
          )}
        </Card>
      ))}
      
      <Button variant="contained" color="primary" href="/upgrade">
        Upgrade Plan
      </Button>
    </div>
  );
};

export default UsageDashboard;
```

---

## 📋 AGENT JR EXECUTION CHECKLIST

### **PHASE 1: DATABASE & AUTH SETUP (3-4 days)**

#### **Day 1: Supabase Setup**
- [ ] Create Supabase project
- [ ] Run SQL schema (users, profiles, subscriptions, quotas, features)
- [ ] Create indexes
- [ ] Test authentication (signup, login, logout)
- [ ] Configure email templates

#### **Day 2: Backend Auth Integration**
- [ ] Install dependencies (`pyjwt`, `python-multipart`, `psycopg2`)
- [ ] Create `auth_middleware.py`
- [ ] Create `feature_flag_service.py`
- [ ] Create `quota_service.py`
- [ ] Create `/api/auth/*` endpoints (login, signup, profile)

#### **Day 3: Protected Endpoints**
- [ ] Add auth middleware to existing endpoints
- [ ] Add quota checks to endpoints
- [ ] Add feature flag checks to premium endpoints
- [ ] Test with Postman/curl

#### **Day 4: Frontend Auth**
- [ ] Install Supabase JS client
- [ ] Create `AuthContext.jsx`
- [ ] Create login/signup pages
- [ ] Create `ProtectedRoute` component
- [ ] Update App.jsx routing

### **PHASE 2: USAGE TRACKING & QUOTAS (2-3 days)**

#### **Day 5: Usage Logging**
- [ ] Create `usage_logs` table
- [ ] Create `usage_tracking_service.py`
- [ ] Add usage logging to all endpoints
- [ ] Create admin endpoint to view usage

#### **Day 6: Quota Enforcement**
- [ ] Test quota limits (free tier)
- [ ] Test quota reset (monthly)
- [ ] Create frontend quota display
- [ ] Create quota warning notifications

#### **Day 7: Upgrade Flow**
- [ ] Create pricing page
- [ ] Create checkout flow (Stripe)
- [ ] Create subscription management page
- [ ] Test tier upgrades

### **PHASE 3: FEATURE FLAGS & TIERS (2 days)**

#### **Day 8: Feature Flag System**
- [ ] Test tier-based feature access
- [ ] Create admin panel for feature flags
- [ ] Test custom feature overrides
- [ ] Document feature flag usage

#### **Day 9: Tier UI Components**
- [ ] Add "Pro" badges to premium features
- [ ] Add upgrade prompts
- [ ] Add feature comparison page
- [ ] Test feature gating

### **PHASE 4: SESSION PERSISTENCE (2 days)**

#### **Day 10: Session Storage**
- [ ] Create `user_sessions` and `saved_analyses` tables
- [ ] Create session management endpoints
- [ ] Implement save/load analysis
- [ ] Test session persistence

#### **Day 11: Frontend Session UI**
- [ ] Create "My Analyses" page
- [ ] Add save/load buttons to analysis pages
- [ ] Add session history
- [ ] Test cross-page resume

### **PHASE 5: ADMIN DASHBOARD (2-3 days)**

#### **Day 12-13: Admin Panel**
- [ ] Create admin authentication
- [ ] Create user management page
- [ ] Create usage analytics page
- [ ] Create feature flag management
- [ ] Create quota override tools

#### **Day 14: Billing Integration**
- [ ] Stripe webhook handling
- [ ] Subscription lifecycle management
- [ ] Payment failure handling
- [ ] Cancellation flow

### **PHASE 6: TESTING & POLISH (2-3 days)**

#### **Day 15: End-to-End Testing**
- [ ] Test free tier signup → usage → quota hit
- [ ] Test upgrade flow
- [ ] Test premium features
- [ ] Test session persistence

#### **Day 16-17: Production Hardening**
- [ ] Add rate limiting (Redis)
- [ ] Add monitoring (Sentry)
- [ ] Add analytics (PostHog/Mixpanel)
- [ ] Security audit
- [ ] Performance testing

---

## 🎯 SUCCESS METRICS

### **Technical:**
- ✅ 100% uptime for auth service
- ✅ < 100ms auth middleware overhead
- ✅ Zero quota bypass vulnerabilities
- ✅ 99.9% quota accuracy

### **Business:**
- ✅ Free → Pro conversion rate: 5%+
- ✅ Monthly recurring revenue: $10K+ (20 Pro users)
- ✅ Churn rate: < 10%
- ✅ Support tickets: < 5/week

### **User Experience:**
- ✅ Signup → first analysis: < 2 minutes
- ✅ Quota visibility: always visible
- ✅ Upgrade flow: < 3 clicks
- ✅ Session restore: 100% success rate

---

## 💰 PRICING JUSTIFICATION

### **FREE TIER:**
**Value:** $0/month  
**Cost to Us:** ~$5/user/month (compute + storage)  
**Strategy:** Lead generation, viral growth

### **PRO TIER:**
**Value:** $499/month  
**Cost to Us:** ~$50/user/month  
**Margin:** 90%  
**Target:** Individual researchers, small labs

### **ENTERPRISE TIER:**
**Value:** $5,000/month  
**Cost to Us:** ~$500/user/month  
**Margin:** 90%  
**Target:** Academic medical centers, pharma companies

---

## ⚔️ FINAL VERDICT

**Timeline:** 2-3 weeks  
**Complexity:** MEDIUM-HIGH  
**ROI:** VERY HIGH (enables monetization)  
**Risk:** LOW (use proven stack: Supabase + Stripe)

**This doctrine gives Agent Jr everything needed to build a production-ready SaaS.**

**Ready to execute, Commander?** ⚔️

— Zo, SaaS Architect

