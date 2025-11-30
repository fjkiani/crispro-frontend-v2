# 📊 ITERATION 5: FRONTEND ARCHITECTURE & USER EXPERIENCE

**Status**: ✅ **COMPLETE**  
**Duration**: 3-4 hours  
**Created**: January 14, 2025

---

### **5.1 TECHNOLOGY STACK**

#### **5.1.1 Core Framework**:
- **React 18.2.0**: Component-based UI library
- **React Router DOM 6.4.4**: Client-side routing
- **Vite 3**: Build tool and dev server
- **TypeScript**: Partial adoption (some `.tsx` files)

#### **5.1.2 UI Libraries**:
- **Material-UI (MUI) 6.5.0**: Primary component library
  - Components: `Box`, `Card`, `Typography`, `Button`, `Chip`, `LinearProgress`, `Accordion`, `Tabs`, `Dialog`
  - Icons: `@mui/icons-material`
  - Charts: `@mui/x-charts`
- **Tailwind CSS 3.2.4**: Utility-first CSS framework
- **Radix UI**: Headless UI components (`@radix-ui/react-select`, `@radix-ui/react-tabs`, `@radix-ui/react-slot`)
- **Class Variance Authority**: Component variant management

#### **5.1.3 State Management**:
- **React Context API**: Primary state management
  - `AuthContext`: User authentication
  - `SporadicContext`: Sporadic cancer workflow state
  - `CoPilotContext`: AI assistant state
  - `ActivityContext`: Activity tracking
  - `AnalysisHistoryContext`: Analysis history
  - `AgentContext`: Agent system state
  - `ClinicalGenomicsContext`: Clinical genomics state
- **Local Storage**: Persistent state (e.g., Kanban tasks)
- **React Hooks**: Custom hooks for API calls and data fetching

#### **5.1.4 Additional Libraries**:
- **Supabase 2.56.0**: Authentication and database
- **React Markdown 9.0.1**: Markdown rendering
- **React Joyride 2.9.3**: User onboarding/tours
- **React Toastify 11.0.5**: Toast notifications
- **Notistack 3.0.2**: Snackbar notifications
- **React Spring 10.0.1**: Animation library
- **DnD Kit**: Drag-and-drop functionality

---

### **5.2 APPLICATION STRUCTURE**

#### **5.2.1 Directory Organization**:
```
src/
├── components/          # Reusable UI components
│   ├── ayesha/         # Ayesha-specific components
│   ├── CoPilot/        # AI assistant components
│   ├── common/         # Shared components
│   ├── dashboard/      # Dashboard components
│   └── ...
├── pages/              # Page-level components (routes)
├── context/            # React Context providers
├── hooks/              # Custom React hooks
├── services/           # API clients and services
├── utils/              # Utility functions
├── config/             # Configuration files
├── constants/          # Constants and enums
└── features/           # Feature modules
```

#### **5.2.2 Entry Point** (`main.jsx`):
```jsx
import { BrowserRouter as Router } from "react-router-dom";
import { StateContextProvider } from "./context";
import App from "./App";

root.render(
  <Router>
    <StateContextProvider>
      <App />
    </StateContextProvider>
  </Router>
);
```

**Provider Hierarchy** (`App.jsx:79-189`):
```jsx
<ErrorBoundary>
  <AuthProvider>
    <AgentProvider>
      <SporadicProvider>
        <CoPilotProvider>
          <AnalysisHistoryProvider>
            <ActivityProvider>
              {/* App content */}
            </ActivityProvider>
          </AnalysisHistoryProvider>
        </CoPilotProvider>
      </SporadicProvider>
    </AgentProvider>
  </AuthProvider>
</ErrorBoundary>
```

---

### **5.3 ROUTING ARCHITECTURE**

#### **5.3.1 Route Structure** (`App.jsx:93-172`):

**Auth Routes** (Public):
- `/login` → `Login`
- `/signup` → `Signup`

**Admin Routes** (Protected):
- `/admin/dashboard` → `AdminDashboard`
- `/admin/users` → `AdminUsers`

**Main Routes**:
- `/` → `Home`
- `/dashboard` → `DoctorDashboard`
- `/profile` → `Profile`
- `/onboarding` → `Onboarding`
- `/medical-records` → `MedicalRecords`
- `/medical-records/:id` → `SingleRecordDetails`
- `/research` → `Research`
- `/mutation-explorer` → `MutationExplorer`
- `/agent-dashboard` → `AgentDashboard`
- `/agents` → `AgentsPage`
- `/agent-studio` → `AgentStudio`

**Ayesha Routes**:
- `/ayesha-complete-care` → `AyeshaCompleteCare`
- `/ayesha-trials` → `AyeshaTrialExplorer`
- `/ayesha-dossiers` → `AyeshaDossierBrowser`
- `/ayesha-dossiers/:nct_id` → `AyeshaDossierDetail`
- `/sporadic-cancer` → `SporadicCancerPage`
- `/ayesha-twin-demo` → `AyeshaTwinDemo`

**Research/Design Routes**:
- `/metastasis` → `MetastasisDashboard`
- `/crispr-designer` → `CrisprDesigner`
- `/protein-synthesis` → `ProteinSynthesis`
- `/structure-predictor` → `StructurePredictor`
- `/validate` → `HypothesisValidator`
- `/food-validator` → `FoodValidatorAB`
- `/batch-food-validator` → `BatchFoodValidator`

**Clinical Routes**:
- `/clinical-genomics` → `ClinicalGenomicsCommandCenter`
- `/myeloma-digital-twin` → `MyelomaDigitalTwin`
- `/radonc-co-pilot` → `RadOncCoPilot`
- `/threat-assessor` → `ThreatAssessor`

#### **5.3.2 Protected Routes**:
- **Pattern**: `<ProtectedRoute><Component /></ProtectedRoute>`
- **Implementation**: Checks authentication status before rendering

---

### **5.4 STATE MANAGEMENT PATTERNS**

#### **5.4.1 Context API Pattern**:

**Example: SporadicContext** (`context/SporadicContext.jsx`):
```jsx
export const SporadicProvider = ({ children }) => {
  const [germlineStatus, setGermlineStatus] = useState('unknown');
  const [tumorContext, setTumorContext] = useState(null);
  const [contextId, setContextId] = useState(null);
  const [dataLevel, setDataLevel] = useState('L0');

  const updateTumorContext = useCallback((data) => {
    if (data?.tumor_context) {
      setTumorContext(data.tumor_context);
      setContextId(data.context_id);
      // Determine data level from completeness score
      const completeness = data.tumor_context.completeness_score || 0;
      if (completeness >= 0.7) setDataLevel('L2');
      else if (completeness >= 0.3) setDataLevel('L1');
      else setDataLevel('L0');
    }
  }, []);

  const getEfficacyPayload = useCallback((basePayload) => {
    return {
      ...basePayload,
      germline_status: germlineStatus,
      tumor_context: tumorContext,
    };
  }, [germlineStatus, tumorContext]);

  return (
    <SporadicContext.Provider value={{ ...state, ...actions }}>
      {children}
    </SporadicContext.Provider>
  );
};
```

**Key Features**:
- **State Slices**: Separate state for different concerns
- **Actions**: `useCallback`-wrapped action functions
- **Computed Values**: Derived state (e.g., `hasTumorContext`, `isSporadic`)
- **Integration Helpers**: `getEfficacyPayload` injects context into API calls

#### **5.4.2 Local Storage Pattern**:
```jsx
// Load from localStorage
const [tasks, setTasks] = useState(() => {
  const savedTasks = localStorage.getItem(KANBAN_TASKS_KEY);
  return savedTasks ? JSON.parse(savedTasks) : [];
});

// Save to localStorage
useEffect(() => {
  localStorage.setItem(KANBAN_TASKS_KEY, JSON.stringify(tasks));
}, [tasks]);
```

---

### **5.5 API INTEGRATION PATTERNS**

#### **5.5.1 API Client Hook** (`hooks/useApiClient.js`):
```jsx
export default function useApiClient(modelId) {
  const API_BASE_URL = import.meta.env.VITE_API_ROOT || 'http://127.0.0.1:8023';

  const client = useMemo(() => {
    const post = async (endpoint, payload = {}, opts = {}) => {
      const controller = new AbortController();
      const timeoutMs = opts.timeoutMs || DEFAULT_TIMEOUT_MS;
      const timer = setTimeout(() => controller.abort(), timeoutMs);
      
      try {
        const body = JSON.stringify({ ...payload, model_id: modelId || 'evo2_7b' });
        const res = await fetch(`${API_BASE_URL}${endpoint}`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body,
          signal: controller.signal,
        });
        
        const json = await res.json().catch(() => ({}));
        if (!res.ok) {
          throw new Error(json?.detail || `HTTP ${res.status}`);
        }
        return json;
      } finally {
        clearTimeout(timer);
      }
    };
    return { post };
  }, [API_BASE_URL, modelId]);

  return client;
}
```

**Features**:
- **Timeout**: 10 minutes default (configurable)
- **AbortController**: Request cancellation
- **Error Handling**: Extracts `detail` from error response
- **Model ID Injection**: Automatically adds `model_id` to payload

#### **5.5.2 Advanced API Client** (`components/ClinicalGenomicsCommandCenter/utils/genomicsUtils.js`):
```jsx
export async function apiPost(path, body, { signal, useCache = true, skipRetry = false } = {}) {
  const url = `${API_BASE}${path}`;
  
  // Check cache first
  if (useCache) {
    const cacheKey = getCacheKey(path, body);
    const cached = getCached(cacheKey);
    if (cached) return cached;
  }
  
  const controller = new AbortController();
  const abortSignal = signal || controller.signal;
  const timeoutId = setTimeout(() => controller.abort(), DEFAULT_TIMEOUT);
  
  const doFetch = async (attempt = 1) => {
    try {
      const response = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
        body: JSON.stringify(body),
        signal: abortSignal
      });
      
      clearTimeout(timeoutId);
      
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.detail || `HTTP ${response.status}`);
      }
      
      const data = await response.json();
      
      // Cache successful result
      if (useCache) {
        const cacheKey = getCacheKey(path, body);
        setCache(cacheKey, data);
      }
      
      return data;
    } catch (error) {
      clearTimeout(timeoutId);
      
      // Retry logic (exponential backoff)
      if (!skipRetry && attempt < 3 && error.name !== 'AbortError') {
        const delay = attempt * 1000; // 1s, 2s
        await new Promise(resolve => setTimeout(resolve, delay));
        return doFetch(attempt + 1);
      }
      
      throw error;
    }
  };
  
  return doFetch();
}
```

**Advanced Features**:
- **Caching**: 10-minute TTL, keyed by path + body
- **Retry Logic**: Exponential backoff (2 retries)
- **Timeout**: 60 seconds default
- **Abort Support**: Signal-based cancellation

#### **5.5.3 Custom Hooks Pattern** (`hooks/useMetastasis.js`):
```jsx
export function useMetastasisAssess(params, enabled = true) {
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [lastFetch, setLastFetch] = useState(null);

  const cacheKey = JSON.stringify(params);
  const CACHE_TTL_MS = 10 * 60 * 1000; // 10 minutes

  const fetchAssessment = useCallback(async () => {
    if (!enabled || !params.mutations || params.mutations.length === 0) {
      return;
    }

    // Check cache freshness
    const now = Date.now();
    if (lastFetch && (now - lastFetch) < CACHE_TTL_MS && data) {
      return; // Use cached data
    }

    setLoading(true);
    setError(null);

    try {
      const response = await fetch(`${API_ROOT}/api/metastasis/assess`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(params)
      });

      if (!response.ok) {
        throw new Error(`Assessment failed (${response.status})`);
      }

      const result = await response.json();
      setData(result);
      setLastFetch(now);
    } catch (err) {
      setError(err.message);
      setData(null);
    } finally {
      setLoading(false);
    }
  }, [cacheKey, enabled, lastFetch, data]);

  useEffect(() => {
    fetchAssessment();
  }, [fetchAssessment]);

  return { data, loading, error, refetch: () => fetchAssessment() };
}
```

**Pattern**:
- **State**: `data`, `loading`, `error`
- **Caching**: TTL-based cache with `lastFetch` timestamp
- **Auto-fetch**: `useEffect` triggers on dependency changes
- **Manual Refetch**: `refetch` function for explicit refresh

---

### **5.6 COMPONENT PATTERNS**

#### **5.6.1 Ayesha Components** (`components/ayesha/`):

**DrugRankingPanel** (`DrugRankingPanel.jsx`):
- **Purpose**: Display ranked drug recommendations
- **Props**: `drugs[]`, `onViewDetails`
- **Features**:
  - Tier badges (supported/consider/insufficient)
  - Evidence badges (RCT, Guideline, ClinVar-Strong)
  - Efficacy score visualization (LinearProgress)
  - Expandable details (Accordion)
  - MUI-based styling

**FoodRankingPanel** (`FoodRankingPanel.jsx`):
- **Purpose**: Display food/supplement recommendations
- **Similar pattern** to DrugRankingPanel

**IntegratedConfidenceBar** (`IntegratedConfidenceBar.jsx`):
- **Purpose**: Visualize integrated confidence scores
- **Shows**: Drug + Food confidence breakdown

**CA125Tracker** (`CA125Tracker.jsx`):
- **Purpose**: Track CA-125 biomarker over time
- **Features**: Time series visualization, response forecasting

**ResistanceAlertBanner** (`ResistanceAlertBanner.jsx`):
- **Purpose**: Display resistance warnings
- **Features**: Alert styling, actionable recommendations

#### **5.6.2 CoPilot Components** (`components/CoPilot/`):

**CoPilot** (Main Component):
- **Purpose**: AI assistant interface
- **Features**: Chat interface, context awareness, action execution

**Q2CRouter** (Question-to-Capability Router):
- **Purpose**: Route user questions to appropriate capabilities
- **Pattern**: Intent classification → Capability selection → Action execution

**Evidence Components**:
- **Purpose**: Display clinical evidence
- **Features**: PubMed citations, ClinVar data, pathway information

**Action Components**:
- **Purpose**: Execute actions (e.g., run analysis, generate report)
- **Features**: Status tracking, progress indicators

#### **5.6.3 Common Components** (`components/common/`):
- **ToolRunner**: Execute tools/analyses
- **LoadingSkeleton**: Loading state placeholders
- **ErrorBoundary**: Error handling wrapper

---

### **5.7 STYLING PATTERNS**

#### **5.7.1 MUI Styling** (Primary):
```jsx
<Box sx={{ p: 4, maxWidth: 1400, mx: 'auto' }}>
  <Card sx={{ p: 3 }}>
    <Typography variant="h5" sx={{ fontWeight: 'bold' }}>
      Title
    </Typography>
    <LinearProgress variant="determinate" value={75} />
  </Card>
</Box>
```

**Features**:
- **SX Prop**: Inline styling with theme access
- **Theme System**: Consistent colors, spacing, typography
- **Responsive**: Breakpoint-based responsive design

#### **5.7.2 Tailwind CSS** (Secondary):
```jsx
<div className="sm:-8 relative flex min-h-screen flex-row bg-white p-4">
  <div className="relative mr-10 hidden sm:flex">
    <Sidebar />
  </div>
</div>
```

**Features**:
- **Utility Classes**: Rapid styling
- **Responsive**: Breakpoint prefixes (`sm:`, `md:`, `lg:`)
- **Custom Configuration**: Tailwind config for custom values

#### **5.7.3 Styled Components** (MUI):
```jsx
const StyledModal = styled(Modal)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '16px',
}));
```

---

### **5.8 KEY INSIGHTS**

#### **Frontend Architecture**:
1. **Hybrid UI Libraries**: MUI (primary) + Tailwind (utility)
2. **Context-Heavy**: Extensive use of React Context for state
3. **Hook-Based**: Custom hooks for API calls and data fetching
4. **Component Organization**: Feature-based + shared components

#### **State Management**:
1. **Context API**: Primary pattern (no Redux)
2. **Local Storage**: Persistent state (tasks, preferences)
3. **Derived State**: Computed values in context providers
4. **Integration Helpers**: Context-aware API payload builders

#### **API Integration**:
1. **Fetch API**: Native fetch (no axios)
2. **AbortController**: Request cancellation
3. **Caching**: TTL-based in-memory cache (10 minutes)
4. **Retry Logic**: Exponential backoff (2-3 attempts)
5. **Error Handling**: Extract `detail` from error responses

#### **Component Patterns**:
1. **Panel Components**: Ranking panels for drugs/foods
2. **Card Components**: MUI Card-based layouts
3. **Modal/Dialog**: MUI Dialog for overlays
4. **Accordion**: Expandable details sections
5. **Loading States**: Skeleton loaders + progress indicators

#### **Routing**:
1. **React Router v6**: Latest routing library
2. **Protected Routes**: Auth-based route protection
3. **Nested Routes**: Patient-scoped routes (`/medical-records/:patientId/research`)
4. **Route Organization**: Feature-based route grouping

---

**Status**: 🔄 **ITERATION 6 IN PROGRESS** - Clinical Systems & Workflows Deep Dive  
**Next**: Complete I6 documentation, then move to I7 (Research & Design Systems)

---
