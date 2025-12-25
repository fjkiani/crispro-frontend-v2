# 🚀 Synthetic Lethality Analyzer V2 — UI/UX Enhancement Plan

**Date:** January 28, 2025  
**Status:** 📋 READY FOR AGENT EXECUTION  
**Predecessor:** `SYNTHETIC_LETHALITY_FRONTEND_PLAN.md` (✅ V1 Complete)  
**Goal:** Enhance UI/UX and integrate AI-powered explanations

---

## 📍 Current State (V1 Complete)

**Location:** `oncology-coPilot/oncology-frontend/src/components/SyntheticLethality/`

**Existing Files:**
```
SyntheticLethality/
├── SyntheticLethalityAnalyzer.jsx    ← Main page (functional)
├── index.js
├── hooks/
│   └── useSyntheticLethality.js      ← API hook (working)
└── components/
    ├── EssentialityScoreCard.jsx     ← Basic cards
    ├── PathwayDependencyDiagram.jsx  ← Static diagram
    ├── TherapyRecommendationList.jsx ← List view
    ├── MutationInputForm.jsx         ← Form inputs
    └── ClinicalDossierModal.jsx      ← Export modal
```

**What Works:**
- ✅ Mutation input and analysis
- ✅ Essentiality scoring display
- ✅ Pathway dependency visualization
- ✅ Drug recommendations
- ✅ Clinical dossier export

**What Needs Enhancement:**
- ❌ Static UI (no animations)
- ❌ Basic card design (not premium feel)
- ❌ No AI explanations
- ❌ No conversational interface
- ❌ Limited interactivity in pathway diagram

---

## 🎯 V2 Enhancement Goals

### 1. **AI-Powered Explanations** (LLM Integration)
Use `src/tools/llm_api.py` to generate:
- Natural language interpretation of results
- "Explain Like I'm a Patient" summaries
- Clinical rationale for drug recommendations
- Q&A capability about the analysis

### 2. **Premium UI/UX**
- Animated transitions and micro-interactions
- Interactive pathway diagram (click to explore)
- Glassmorphism/modern card design
- Better color hierarchy and typography
- Loading skeletons and progress indicators

### 3. **Enhanced Interactivity**
- Click pathway nodes to see drug details
- Hover effects with tooltips
- Expandable sections for deep-dive
- Side-by-side comparison view

---

## 🔧 Implementation Tasks

### Task 1: Create LLM Explanation Hook

**File:** `hooks/useLLMExplanation.js`

**Purpose:** Call LLM API to generate natural language explanations

```javascript
// hooks/useLLMExplanation.js
import { useState, useCallback } from 'react';

const API_BASE_URL = import.meta.env.VITE_API_ROOT || 'http://localhost:8000';

/**
 * Hook for AI-powered explanations using LLM
 */
export function useLLMExplanation() {
  const [loading, setLoading] = useState(false);
  const [explanation, setExplanation] = useState(null);
  const [error, setError] = useState(null);

  /**
   * Generate explanation for synthetic lethality results
   * @param {Object} results - Analysis results
   * @param {string} audienceType - "clinician" | "patient" | "researcher"
   */
  const generateExplanation = useCallback(async (results, audienceType = 'clinician') => {
    setLoading(true);
    setError(null);

    const prompt = buildExplanationPrompt(results, audienceType);

    try {
      // Call backend LLM endpoint (create if needed)
      const response = await fetch(`${API_BASE_URL}/api/llm/explain`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          prompt,
          provider: 'gemini', // or configurable
          context: 'synthetic_lethality'
        })
      });

      if (!response.ok) throw new Error('LLM API failed');
      
      const data = await response.json();
      setExplanation(data.explanation);
      return data.explanation;
    } catch (err) {
      setError(err.message);
      return null;
    } finally {
      setLoading(false);
    }
  }, []);

  /**
   * Ask a follow-up question about the analysis
   */
  const askQuestion = useCallback(async (question, context) => {
    setLoading(true);
    
    const prompt = `
Given this synthetic lethality analysis context:
${JSON.stringify(context, null, 2)}

User question: ${question}

Provide a clear, evidence-based answer suitable for a clinical oncologist.
`;

    try {
      const response = await fetch(`${API_BASE_URL}/api/llm/chat`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ prompt, provider: 'gemini' })
      });

      const data = await response.json();
      return data.response;
    } catch (err) {
      setError(err.message);
      return null;
    } finally {
      setLoading(false);
    }
  }, []);

  return {
    generateExplanation,
    askQuestion,
    explanation,
    loading,
    error
  };
}

// Helper: Build prompt based on audience
function buildExplanationPrompt(results, audienceType) {
  const { essentiality, pathway_analysis, recommended_therapies } = results;

  const baseContext = `
## Synthetic Lethality Analysis Results

### Essentiality Scores:
${essentiality?.map(e => `- ${e.gene}: ${(e.score * 100).toFixed(0)}% (${e.pathwayImpact})`).join('\n')}

### Pathway Analysis:
- Broken Pathways: ${pathway_analysis?.broken_pathways?.join(', ') || 'None'}
- Essential Backups: ${pathway_analysis?.essential_pathways?.join(', ') || 'None'}
- Double-Hit Effect: ${pathway_analysis?.double_hit_detected ? 'Yes' : 'No'}

### Top Therapies:
${recommended_therapies?.slice(0, 3).map((t, i) => `${i + 1}. ${t.drug} (${t.target}) - ${(t.confidence * 100).toFixed(0)}% confidence`).join('\n')}
`;

  const audienceInstructions = {
    clinician: `
Explain these results for a practicing oncologist. Include:
1. Clinical significance of each gene's essentiality
2. Mechanism of synthetic lethality
3. Rationale for drug recommendations
4. Key monitoring considerations
Use medical terminology appropriately.
`,
    patient: `
Explain these results for a cancer patient with no medical background. Include:
1. What the genetic mutations mean in simple terms
2. Why certain treatments might work better
3. What "synthetic lethality" means (use analogies)
4. Reassuring but honest tone
Avoid jargon. Use 8th-grade reading level.
`,
    researcher: `
Provide a detailed scientific explanation including:
1. Molecular mechanisms of pathway disruption
2. Evidence from literature supporting synthetic lethality
3. Potential resistance mechanisms to monitor
4. Suggestions for combination therapy rationale
Include relevant gene/pathway references.
`
  };

  return `${baseContext}\n\n${audienceInstructions[audienceType]}`;
}

export default useLLMExplanation;
```

---

### Task 2: Create Backend LLM Endpoint

**File:** `oncology-coPilot/oncology-backend-minimal/api/routers/llm.py`

**Purpose:** Backend endpoint that wraps `src/tools/llm_api.py`

```python
# api/routers/llm.py
from fastapi import APIRouter, HTTPException
from typing import Dict, Any
import sys
import os

# Add tools path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'src', 'tools'))

router = APIRouter(prefix="/api/llm", tags=["llm"])

@router.post("/explain")
async def explain_results(request: Dict[str, Any]):
    """Generate LLM explanation for analysis results."""
    try:
        from llm_api import query_llm
        
        prompt = request.get("prompt", "")
        provider = request.get("provider", "gemini")
        
        if not prompt:
            raise HTTPException(status_code=400, detail="prompt required")
        
        response = query_llm(prompt, provider=provider)
        
        return {"explanation": response, "provider": provider}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"LLM error: {e}")

@router.post("/chat")
async def chat_with_llm(request: Dict[str, Any]):
    """Chat endpoint for follow-up questions."""
    try:
        from llm_api import query_llm
        
        prompt = request.get("prompt", "")
        provider = request.get("provider", "gemini")
        
        response = query_llm(prompt, provider=provider)
        
        return {"response": response}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"LLM chat error: {e}")
```

**Register in `api/index.py`:**
```python
from .routers import llm
app.include_router(llm.router)
```

---

### Task 3: Create AI Explanation Panel Component

**File:** `components/AIExplanationPanel.jsx`

**Purpose:** Display AI-generated explanations with audience selector

```jsx
// components/AIExplanationPanel.jsx
import React, { useState } from 'react';
import {
  Box,
  Paper,
  Typography,
  Button,
  ToggleButton,
  ToggleButtonGroup,
  TextField,
  CircularProgress,
  Collapse,
  IconButton,
  Divider,
  Stack,
  Chip,
  Alert
} from '@mui/material';
import {
  Psychology,
  LocalHospital,
  Person,
  Science,
  Send,
  ExpandMore,
  ExpandLess,
  AutoAwesome,
  QuestionAnswer
} from '@mui/icons-material';
import { useLLMExplanation } from '../hooks/useLLMExplanation';

const AIExplanationPanel = ({ results }) => {
  const [audienceType, setAudienceType] = useState('clinician');
  const [question, setQuestion] = useState('');
  const [chatHistory, setChatHistory] = useState([]);
  const [expanded, setExpanded] = useState(true);
  
  const {
    generateExplanation,
    askQuestion,
    explanation,
    loading,
    error
  } = useLLMExplanation();

  const handleGenerateExplanation = async () => {
    await generateExplanation(results, audienceType);
  };

  const handleAskQuestion = async () => {
    if (!question.trim()) return;
    
    const answer = await askQuestion(question, results);
    if (answer) {
      setChatHistory(prev => [...prev, { q: question, a: answer }]);
      setQuestion('');
    }
  };

  return (
    <Paper 
      elevation={3} 
      sx={{ 
        p: 3, 
        borderRadius: 3,
        background: 'linear-gradient(145deg, #f5f7fa 0%, #e8ecf1 100%)',
        border: '1px solid',
        borderColor: 'primary.light'
      }}
    >
      {/* Header */}
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
        <Stack direction="row" spacing={1} alignItems="center">
          <AutoAwesome sx={{ color: 'primary.main' }} />
          <Typography variant="h6" fontWeight="bold">
            AI Clinical Interpretation
          </Typography>
          <Chip label="Powered by LLM" size="small" color="primary" variant="outlined" />
        </Stack>
        <IconButton onClick={() => setExpanded(!expanded)}>
          {expanded ? <ExpandLess /> : <ExpandMore />}
        </IconButton>
      </Box>

      <Collapse in={expanded}>
        {/* Audience Selector */}
        <Box sx={{ mb: 3 }}>
          <Typography variant="body2" color="text.secondary" gutterBottom>
            Generate explanation for:
          </Typography>
          <ToggleButtonGroup
            value={audienceType}
            exclusive
            onChange={(e, v) => v && setAudienceType(v)}
            size="small"
          >
            <ToggleButton value="clinician">
              <LocalHospital sx={{ mr: 1 }} /> Clinician
            </ToggleButton>
            <ToggleButton value="patient">
              <Person sx={{ mr: 1 }} /> Patient
            </ToggleButton>
            <ToggleButton value="researcher">
              <Science sx={{ mr: 1 }} /> Researcher
            </ToggleButton>
          </ToggleButtonGroup>
        </Box>

        {/* Generate Button */}
        <Button
          variant="contained"
          startIcon={loading ? <CircularProgress size={20} /> : <Psychology />}
          onClick={handleGenerateExplanation}
          disabled={loading || !results}
          fullWidth
          sx={{ mb: 2 }}
        >
          {loading ? 'Generating...' : 'Generate AI Explanation'}
        </Button>

        {error && <Alert severity="error" sx={{ mb: 2 }}>{error}</Alert>}

        {/* Explanation Display */}
        {explanation && (
          <Paper 
            variant="outlined" 
            sx={{ 
              p: 2, 
              mb: 3, 
              backgroundColor: 'white',
              maxHeight: 400,
              overflow: 'auto'
            }}
          >
            <Typography 
              variant="body1" 
              sx={{ 
                whiteSpace: 'pre-wrap',
                lineHeight: 1.7
              }}
            >
              {explanation}
            </Typography>
          </Paper>
        )}

        <Divider sx={{ my: 2 }} />

        {/* Q&A Section */}
        <Box>
          <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 2 }}>
            <QuestionAnswer color="secondary" />
            <Typography variant="subtitle1" fontWeight="medium">
              Ask a Question
            </Typography>
          </Stack>

          <Stack direction="row" spacing={1}>
            <TextField
              fullWidth
              size="small"
              placeholder="e.g., Why is PARP inhibitor recommended over platinum?"
              value={question}
              onChange={(e) => setQuestion(e.target.value)}
              onKeyPress={(e) => e.key === 'Enter' && handleAskQuestion()}
            />
            <Button
              variant="contained"
              color="secondary"
              onClick={handleAskQuestion}
              disabled={loading || !question.trim()}
            >
              <Send />
            </Button>
          </Stack>

          {/* Chat History */}
          {chatHistory.length > 0 && (
            <Box sx={{ mt: 2 }}>
              {chatHistory.map((chat, idx) => (
                <Box key={idx} sx={{ mb: 2 }}>
                  <Typography variant="body2" color="primary.main" fontWeight="medium">
                    Q: {chat.q}
                  </Typography>
                  <Typography variant="body2" sx={{ pl: 2, mt: 0.5 }}>
                    {chat.a}
                  </Typography>
                </Box>
              ))}
            </Box>
          )}
        </Box>
      </Collapse>
    </Paper>
  );
};

export default AIExplanationPanel;
```

---

### Task 4: Enhanced EssentialityScoreCard (Animated)

**File:** Update `components/EssentialityScoreCard.jsx`

**Enhancements:**
- Animated progress bar (count-up effect)
- Glassmorphism card design
- Hover tooltip with details
- Pulsing effect for high scores

```jsx
// Key additions to EssentialityScoreCard.jsx

import { keyframes } from '@mui/system';

// Animation keyframes
const pulseAnimation = keyframes`
  0% { box-shadow: 0 0 0 0 rgba(244, 67, 54, 0.4); }
  70% { box-shadow: 0 0 0 10px rgba(244, 67, 54, 0); }
  100% { box-shadow: 0 0 0 0 rgba(244, 67, 54, 0); }
`;

const countUpAnimation = keyframes`
  from { width: 0%; }
  to { width: var(--target-width); }
`;

// Card styling with glassmorphism
const cardStyles = {
  background: 'rgba(255, 255, 255, 0.85)',
  backdropFilter: 'blur(10px)',
  border: '1px solid rgba(255, 255, 255, 0.3)',
  borderRadius: 4,
  transition: 'all 0.3s ease',
  '&:hover': {
    transform: 'translateY(-4px)',
    boxShadow: '0 12px 24px rgba(0,0,0,0.15)'
  }
};

// For high essentiality scores, add pulsing border
const highScoreStyles = {
  animation: `${pulseAnimation} 2s infinite`,
  borderColor: 'error.main'
};
```

---

### Task 5: Interactive Pathway Diagram

**File:** Update `components/PathwayDependencyDiagram.jsx`

**Enhancements:**
- Clickable pathway nodes
- Animated arrows/connections
- Tooltip on hover with drug info
- Visual highlighting of selected pathway

```jsx
// Key additions to PathwayDependencyDiagram.jsx

// Add state for selected pathway
const [selectedPathway, setSelectedPathway] = useState(null);

// Animated connection line
const ConnectionLine = ({ active }) => (
  <Box
    sx={{
      width: 60,
      height: 4,
      background: active 
        ? 'linear-gradient(90deg, #4caf50, #2196f3)' 
        : '#e0e0e0',
      borderRadius: 2,
      position: 'relative',
      overflow: 'hidden',
      '&::after': active ? {
        content: '""',
        position: 'absolute',
        top: 0,
        left: '-100%',
        width: '100%',
        height: '100%',
        background: 'linear-gradient(90deg, transparent, rgba(255,255,255,0.8), transparent)',
        animation: 'shimmer 2s infinite'
      } : {}
    }}
  />
);

// Clickable pathway chip
const PathwayChip = ({ pathway, type, selected, onClick }) => (
  <Chip
    label={pathway}
    color={type === 'broken' ? 'error' : type === 'essential' ? 'warning' : 'success'}
    onClick={() => onClick(pathway)}
    sx={{
      cursor: 'pointer',
      transition: 'all 0.2s',
      transform: selected ? 'scale(1.1)' : 'scale(1)',
      boxShadow: selected ? '0 4px 12px rgba(0,0,0,0.2)' : 'none'
    }}
  />
);
```

---

### Task 6: Enhanced Clinical Dossier Modal

**File:** Update `components/ClinicalDossierModal.jsx`

**Enhancements:**
- AI-generated summary section
- Better PDF export formatting
- Clinician signature placeholder
- Version/timestamp tracking

---

### Task 7: Integrate AI Panel into Main Page

**File:** Update `SyntheticLethalityAnalyzer.jsx`

**Add AIExplanationPanel after results:**

```jsx
import AIExplanationPanel from './components/AIExplanationPanel';

// In the results section, add:
{results && (
  <Stack spacing={3}>
    {/* ... existing components ... */}
    
    {/* AI Explanation Panel */}
    <AIExplanationPanel results={results} />
    
    {/* ... therapy recommendations ... */}
  </Stack>
)}
```

---

## 📋 Task Checklist for Agent

| # | Task | Priority | Complexity |
|---|------|----------|------------|
| 1 | Create `useLLMExplanation.js` hook | HIGH | Medium |
| 2 | Create backend `/api/llm/explain` endpoint | HIGH | Low |
| 3 | Create `AIExplanationPanel.jsx` component | HIGH | Medium |
| 4 | Enhance `EssentialityScoreCard.jsx` with animations | MEDIUM | Low |
| 5 | Enhance `PathwayDependencyDiagram.jsx` interactivity | MEDIUM | Medium |
| 6 | Improve `ClinicalDossierModal.jsx` with AI summary | MEDIUM | Medium |
| 7 | Integrate AI panel into main page | HIGH | Low |
| 8 | Add loading skeletons for better UX | LOW | Low |
| 9 | Add dark mode support | LOW | Medium |
| 10 | Mobile responsive improvements | LOW | Medium |

---

## 🎨 Design Tokens

```javascript
// Design constants for consistency
const designTokens = {
  colors: {
    highEssentiality: '#f44336',   // Red
    moderateEssentiality: '#ff9800', // Orange
    lowEssentiality: '#4caf50',     // Green
    
    brokenPathway: '#ef5350',
    essentialPathway: '#ffa726',
    drugTarget: '#66bb6a',
    
    aiAccent: '#7c4dff',           // Purple for AI features
    glassBg: 'rgba(255,255,255,0.85)'
  },
  
  animations: {
    cardHover: 'transform 0.3s ease, box-shadow 0.3s ease',
    fadeIn: 'opacity 0.5s ease-in',
    slideUp: 'transform 0.4s ease-out'
  },
  
  shadows: {
    card: '0 4px 12px rgba(0,0,0,0.08)',
    cardHover: '0 12px 24px rgba(0,0,0,0.15)',
    aiPanel: '0 8px 32px rgba(124,77,255,0.15)'
  }
};
```

---

## 🧪 Testing Scenarios

1. **AI Explanation Test:**
   - Load Ayesha's MBD4+TP53 case
   - Generate clinician explanation → Verify medical accuracy
   - Generate patient explanation → Verify readability
   - Ask follow-up question → Verify contextual answer

2. **UI Animation Test:**
   - Verify score cards animate on load
   - Verify pathway diagram click interactions
   - Verify hover effects work smoothly

3. **Export Test:**
   - Generate dossier with AI summary
   - Export as PDF → Verify formatting
   - Copy to clipboard → Verify content

---

## 📂 Files to Create/Modify

### New Files:
```
oncology-coPilot/oncology-frontend/src/components/SyntheticLethality/
├── hooks/
│   └── useLLMExplanation.js       ← NEW
└── components/
    └── AIExplanationPanel.jsx     ← NEW

oncology-coPilot/oncology-backend-minimal/api/routers/
└── llm.py                         ← NEW
```

### Files to Modify:
```
- SyntheticLethalityAnalyzer.jsx   (add AI panel)
- EssentialityScoreCard.jsx        (add animations)
- PathwayDependencyDiagram.jsx     (add interactivity)
- ClinicalDossierModal.jsx         (add AI summary)
- api/index.py                     (register LLM router)
```

---

## 🔑 LLM API Configuration

**Required in `.env`:**
```
GEMINI_API_KEY=your_gemini_key_here
# OR
OPENAI_API_KEY=your_openai_key_here
```

**Fallback:** If no API key configured, show message and disable AI features gracefully.

---

## ✅ Success Criteria

| Criteria | Measurement |
|----------|-------------|
| AI explanations work | Can generate for all 3 audiences |
| Q&A works | Can ask 3+ follow-up questions |
| Animations smooth | 60fps, no jank |
| Interactive diagram | Click pathways to see details |
| Export includes AI | Dossier has AI summary section |
| Graceful fallback | Works without API key (AI disabled) |

---

**Agent Instructions:**

1. Start with Task 2 (backend endpoint) — it's the foundation
2. Then Task 1 (hook) — connects frontend to backend  
3. Then Task 3 (AI panel) — user-facing feature
4. Then Tasks 4-6 (UI polish) — visual enhancements
5. Finally Task 7 (integration) — bring it together

Test with Ayesha's case (MBD4 + TP53) at each step.

---

**Ready for agent execution! 🚀**




