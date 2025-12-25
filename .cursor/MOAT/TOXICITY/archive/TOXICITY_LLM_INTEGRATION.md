# 🤖 LLM Integration for Toxicity Risk Assessment

**Date:** January 28, 2025  
**Status:** ✅ **COMPLETE**

---

## 📋 Summary

Added AI-powered LLM explanations to the toxicity risk assessment system, enabling patient-friendly, clinician-focused, and researcher-level explanations of toxicity risk results.

---

## 🎯 What Was Added

### **1. New Hook: `useToxicityLLM.js`**

**Location:** `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/hooks/useToxicityLLM.js`

**Features:**
- `generateExplanation(result, audienceType)` - Generate AI explanation for toxicity risk
- `askQuestion(question, context)` - Ask follow-up questions about toxicity
- `clearExplanation()` - Clear current explanation
- Supports 3 audience types: `clinician`, `patient`, `researcher`

**API Integration:**
- Uses `/api/llm/explain` endpoint for explanations
- Uses `/api/llm/chat` endpoint for follow-up questions
- Provider: Gemini (configurable)

---

### **2. Enhanced Component: `ToxicityRiskCard.jsx`**

**Location:** `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx`

**New Features:**
1. **AI Explanation Button** - "Generate AI Explanation" button with loading state
2. **Audience Selector** - Dropdown to choose explanation style:
   - **Clinician**: Medical terminology, clinical significance, monitoring considerations
   - **Patient**: Simple language (8th-grade level), empathetic tone, no jargon
   - **Researcher**: Detailed scientific explanation, mechanisms, literature references
3. **Explanation Display** - Collapsible section showing AI-generated explanation
4. **Mitigating Foods Display** - Shows recommended foods/supplements (already existed, now visible)

**UI Components Added:**
- `AutoAwesomeIcon` - Sparkle icon for AI features
- `FormControl` + `Select` - Audience type selector
- `Collapse` - Smooth expand/collapse for explanation
- `CloseIcon` - Button to clear explanation

---

## 🔧 Technical Details

### **Prompt Engineering**

The hook builds context-aware prompts based on:

1. **Drug Information:**
   - Drug name
   - Mechanism of Action (MoA)
   - Risk level and score
   - Confidence level

2. **Risk Summary:**
   - Reason for risk assessment
   - Contributing factors (germline variants, pharmacogenes, pathway overlaps)
   - Mitigating foods/supplements

3. **Audience-Specific Instructions:**
   - **Clinician**: Medical terminology, clinical decisions, monitoring
   - **Patient**: Simple language, empathetic, no jargon
   - **Researcher**: Detailed mechanisms, literature, research gaps

### **Example Prompt (Clinician):**

```
## Toxicity Risk Assessment Results

### Drug Information:
- Drug: carboplatin
- Mechanism of Action: platinum_agent
- Risk Level: HIGH (Score: 100%)
- Confidence: 48%

### Risk Summary:
MoA overlaps toxicity pathways with germline variants...

### Contributing Factors:
- BRCA1 variant (Type: germline, Weight: 50%, Confidence: 80%)

### Mitigating Foods/Supplements:
- NAC (N-Acetyl Cysteine): 600mg twice daily (post-chemo) - Glutathione precursor, supports DNA repair enzymes
- Vitamin D3: 5000 IU daily (continuous) - Modulates DNA repair gene expression

Explain these toxicity risk results for a practicing oncologist. Include:
1. Clinical significance of the risk score...
2. Mechanism of toxicity...
3. Key risk factors identified...
4. Rationale for mitigating foods/supplements...
5. Monitoring considerations...
6. Confidence in the assessment...
```

---

## 📊 Usage Examples

### **In ToxicityRiskCard Component:**

```jsx
import { useToxicityLLM } from '../hooks/useToxicityLLM';

const { 
  generateExplanation, 
  explanation, 
  loading: llmLoading, 
  error: llmError,
  clearExplanation 
} = useToxicityLLM();

// Generate explanation
await generateExplanation(result, 'patient');

// Display explanation
{explanation && (
  <Paper>
    <Typography>{explanation}</Typography>
  </Paper>
)}
```

### **In Standalone Page:**

The `ToxicityRiskAssessment.jsx` page automatically benefits from the enhanced card component - no additional code needed.

---

## 🎨 UI/UX Features

1. **Loading States:**
   - Button shows "Generating Explanation..." while loading
   - Button disabled during generation

2. **Error Handling:**
   - Error alerts displayed if LLM call fails
   - Graceful fallback to original reason text

3. **User Control:**
   - Audience selector before generation
   - Close button to clear explanation
   - Collapsible explanation section

4. **Visual Design:**
   - Sparkle icon (✨) indicates AI features
   - Left border accent on explanation box
   - Consistent with existing Material-UI design

---

## 🔌 API Endpoints Used

### **1. `/api/llm/explain` (POST)**

**Request:**
```json
{
  "prompt": "Full prompt string with context and instructions",
  "provider": "gemini",
  "context": "toxicity_risk"
}
```

**Response:**
```json
{
  "explanation": "Generated explanation text...",
  "provider": "gemini",
  "context": "toxicity_risk"
}
```

### **2. `/api/llm/chat` (POST)**

**Request:**
```json
{
  "prompt": "User question with context",
  "provider": "gemini"
}
```

**Response:**
```json
{
  "response": "LLM answer text..."
}
```

---

## ✅ Testing Checklist

- [x] Hook created and exported
- [x] Component enhanced with LLM features
- [x] Audience selector working (clinician/patient/researcher)
- [x] Explanation generation working
- [x] Error handling implemented
- [x] Loading states implemented
- [x] UI components styled correctly
- [x] Mitigating foods display working
- [x] No linter errors

---

## 🚀 Next Steps (Optional Enhancements)

1. **Caching:** Cache explanations to avoid re-generating for same results
2. **Streaming:** Stream LLM responses for better UX
3. **History:** Save explanation history for user reference
4. **Export:** Allow exporting explanations as PDF/text
5. **Multi-language:** Support explanations in multiple languages
6. **Custom Prompts:** Allow users to customize prompt instructions

---

## 📝 Files Modified

| File | Changes |
|------|---------|
| `hooks/useToxicityLLM.js` | ✅ **NEW** - LLM hook for toxicity explanations |
| `cards/ToxicityRiskCard.jsx` | ✅ **ENHANCED** - Added LLM explanation UI |

---

## 🎯 Benefits

1. **Patient Communication:** Clear, empathetic explanations for patients
2. **Clinical Decision Support:** Detailed explanations for oncologists
3. **Research Context:** Scientific explanations for researchers
4. **Accessibility:** Multiple audience types improve understanding
5. **User Experience:** Interactive AI features enhance engagement

---

**Status:** ✅ **PRODUCTION READY**

All LLM features integrated and tested. Ready for use in toxicity risk assessment pages.

---

**Last Updated:** January 28, 2025  
**Implemented By:** Zo


