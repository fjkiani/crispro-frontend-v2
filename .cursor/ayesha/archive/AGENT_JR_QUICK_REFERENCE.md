# ⚔️ AGENT JR - QUICK REFERENCE GUIDE ⚔️

**Purpose**: Quick lookup for common patterns, existing code references, and execution shortcuts

---

## 🔍 **EXISTING CODE REFERENCES**

### **Backend Services to Reuse**:
1. **`HybridTrialSearchService`** - `api/services/hybrid_trial_search.py`
   - Already supports `germline_status` and `tumor_context` parameters
   - Returns `List[Dict[str, Any]]` with trial data
   - **USE THIS** - don't rewrite search logic

2. **`ClinicalTrialSearchService`** - `api/services/clinical_trial_search_service.py`
   - AstraDB semantic search
   - **USE THIS** - already integrated with HybridTrialSearchService

### **Frontend Components to Reference**:
1. **`ResultsDisplay.jsx`** - `src/components/research/ResultsDisplay.jsx`
   - Shows trial card styling patterns
   - LocationCard integration
   - BiomarkerMatchBadge usage
   - **REFERENCE THIS** for styling

2. **`SporadicProvenanceCard.jsx`** - `src/components/sporadic/SporadicProvenanceCard.jsx`
   - Provenance display patterns
   - **REFERENCE THIS** for ProvenanceCard component

---

## 📋 **COMMON PATTERNS**

### **Backend Pattern - Service Module**:
```python
# api/services/ayesha_trial_matching/eligibility_filters.py
import logging
from typing import Dict, List, Any
from api.services.hybrid_trial_search import HybridTrialSearchService

logger = logging.getLogger(__name__)

class EligibilityFilters:
    """Hard eligibility filters for Ayesha's trials."""
    
    def __init__(self):
        self.search_service = HybridTrialSearchService()
    
    async def apply_hard_filters(self, profile: Dict) -> List[Dict]:
        """Apply all hard filters."""
        # Implementation
        pass
```

### **Frontend Pattern - Display Component**:
```jsx
// src/components/trials/TrialReasoningSection.jsx
import React from 'react';
import { Box, Typography, Chip } from '@mui/material';

export default function TrialReasoningSection({ reasoning }) {
  if (!reasoning?.why_eligible?.length) return null;
  
  return (
    <Box sx={{ mb: 2 }}>
      <Typography variant="h6">Why Eligible</Typography>
      {reasoning.why_eligible.map((reason, idx) => (
        <Chip key={idx} label={reason} sx={{ mr: 1, mb: 1 }} />
      ))}
    </Box>
  );
}
```

### **Frontend Pattern - API Hook**:
```jsx
// src/hooks/useAyeshaTrials.js
import { useState, useEffect } from 'react';

export default function useAyeshaTrials(profile) {
  const [trials, setTrials] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  
  useEffect(() => {
    const fetchTrials = async () => {
      setLoading(true);
      try {
        const response = await fetch('/api/ayesha/trials/search', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(profile)
        });
        const data = await response.json();
        setTrials(data.trials || []);
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };
    
    if (profile) fetchTrials();
  }, [profile]);
  
  return { trials, loading, error };
}
```

---

## 🚨 **CRITICAL REMINDERS**

### **NO HARDCODING**:
- ❌ Don't hardcode trial counts (use `trials?.length || 0`)
- ❌ Don't hardcode trial data (all from backend)
- ❌ Don't hardcode scoring thresholds (use backend values)
- ✅ Use `.map()` for ALL arrays
- ✅ Use optional chaining (`?.`)
- ✅ Provide fallbacks (`|| 'Unknown'`)

### **USE EXISTING SERVICES**:
- ✅ `HybridTrialSearchService` - already exists, just call it
- ✅ `ClinicalTrialSearchService` - already exists, don't rewrite
- ❌ Don't create new search logic from scratch

### **FOLLOW MODULAR PATTERNS**:
- ✅ Create folder structure: `ayesha_trial_matching/`
- ✅ One file per responsibility
- ✅ Test each module independently
- ✅ Build in dependency order (see modular plan)

---

## 🔧 **QUICK COMMANDS**

### **Backend Testing**:
```bash
# Start backend
cd oncology-coPilot/oncology-backend-minimal
venv/bin/uvicorn api.main:app --reload

# Test endpoint
curl -X POST http://127.0.0.1:8000/api/ayesha/trials/search \
  -H 'Content-Type: application/json' \
  -d '{
    "ca125": 2842,
    "germline_status": "NEGATIVE",
    "stage": "IVB",
    "treatment_line": 0,
    "location": "NYC"
  }'
```

### **Frontend Testing**:
```bash
# Start frontend
cd oncology-coPilot/oncology-frontend
npm start

# Navigate to
http://localhost:3000/ayesha-trials
```

### **AstraDB Seeding**:
```bash
cd oncology-coPilot/oncology-backend-minimal
venv/bin/python scripts/seed_astra_trials.py --disease ovarian --count 200
```

---

## 📁 **FILE STRUCTURE QUICK REFERENCE**

```
Backend:
api/
├── schemas/
│   └── ayesha_trials.py
├── services/
│   ├── ca125_intelligence.py
│   └── ayesha_trial_matching/
│       ├── __init__.py
│       ├── eligibility_filters.py
│       ├── scoring_engine.py
│       ├── reasoning_generator.py
│       └── match_orchestrator.py
└── routers/
    └── ayesha_trials.py

Frontend:
src/
├── components/
│   ├── ayesha/
│   │   ├── CA125Tracker.jsx
│   │   ├── AyeshaClinicalProfile.jsx
│   │   └── ProvenanceCard.jsx
│   └── trials/
│       ├── TrialMatchCard.jsx
│       ├── TrialReasoningSection.jsx
│       ├── TrialConditionalSection.jsx
│       └── TrialRedFlagsSection.jsx
├── pages/
│   └── AyeshaTrialExplorer.jsx
└── hooks/
    └── useAyeshaTrials.js
```

---

## ✅ **CHECKLIST BEFORE STARTING EACH MODULE**

- [ ] Read existing similar code (reference files above)
- [ ] Understand module dependencies (see modular plan)
- [ ] Check if service/component already exists (don't duplicate)
- [ ] Plan module structure (single responsibility)
- [ ] Write module (follow patterns above)
- [ ] Test module independently
- [ ] Update progress tracker

---

## 🎯 **SUCCESS CRITERIA PER MODULE**

Each module should:
- ✅ Have single, clear responsibility
- ✅ Handle missing data gracefully
- ✅ Use dynamic patterns (no hardcoding)
- ✅ Follow existing code patterns
- ✅ Have proper error handling
- ✅ Be testable independently

---

**Last Updated**: [Date]  
**For**: Agent Jr


