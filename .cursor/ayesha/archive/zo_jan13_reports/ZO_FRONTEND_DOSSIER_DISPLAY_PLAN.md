# ⚔️ FRONTEND PLAN: Display All Dossiers for Ayesha & Oncologist

**Date**: November 17, 2025  
**Commander**: Zo  
**Mission**: Enable Ayesha and her oncologist to browse all 60 trial intelligence dossiers through a beautiful, actionable interface

---

## 📊 CURRENT STATE ASSESSMENT

### What We Have:
- ✅ **10 Dossier Files** (generated, markdown format)
- ✅ **60 Trial Database** (scored, filtered, ready for dossier generation)
- ✅ **Existing Frontend** (`AyeshaTrialExplorer.jsx`) - Displays trials from `/api/ayesha/complete_care_v2`
- ✅ **TrialMatchCard Component** - Beautiful trial cards with reasoning

### What's Missing:
- ❌ **Backend API** to serve dossier content
- ❌ **Dossier detail view** (full markdown rendering)
- ❌ **Generate remaining 50 dossiers** (only 10 exist now)
- ❌ **Export/share functionality** (PDF, email to oncologist)

---

## 🎯 THREE-PHASE EXECUTION PLAN

### **PHASE 1: Generate All 60 Dossiers** (90 minutes)
**Why First**: Need complete dataset before building UI

#### Task 1.1: Fix Rate Limiting (DONE ✅)
- Updated pipeline.py: 15s → 30s between LLM calls
- Prevents Gemini API quota errors

#### Task 1.2: Generate Remaining 50 Dossiers (90 min)
**Command**:
```bash
cd oncology-coPilot/oncology-backend-minimal
python3 find_trials_FROM_FRESH_TABLE.py
```

**Expected Output**:
- 60 markdown files in `.cursor/ayesha/zo_fresh_dossiers/`
- All with LLM-enhanced analysis
- Tier labels (TOP_TIER, GOOD_TIER)

**Status**: Ready to execute

---

### **PHASE 2: Backend API for Dossiers** (2 hours)
**Why**: Serve dossier content to frontend

#### Task 2.1: Create Dossier API Endpoint (60 min)

**File**: `api/routers/ayesha_dossiers.py` (NEW)

```python
"""
Ayesha Dossiers API
Serves trial intelligence dossiers for frontend display.
"""
from fastapi import APIRouter, HTTPException
from pathlib import Path
from typing import List, Dict, Any, Optional
import json
import re

router = APIRouter(prefix="/api/ayesha/dossiers", tags=["Ayesha Dossiers"])

DOSSIER_DIR = Path(__file__).parent.parent.parent.parent.parent / ".cursor" / "ayesha" / "zo_fresh_dossiers"

@router.get("/list")
async def list_dossiers() -> Dict[str, Any]:
    """
    List all available dossiers with metadata.
    
    Returns:
        {
            "total": 60,
            "dossiers": [
                {
                    "nct_id": "NCT06619236",
                    "tier": "TOP_TIER",
                    "match_score": 0.97,
                    "title": "Trial Title...",
                    "phase": "Phase II",
                    "has_llm_analysis": true,
                    "file_path": "INTELLIGENCE_NCT06619236_TOP_TIER.md"
                },
                ...
            ]
        }
    """
    if not DOSSIER_DIR.exists():
        raise HTTPException(status_code=404, detail="Dossier directory not found")
    
    dossiers = []
    for file_path in DOSSIER_DIR.glob("INTELLIGENCE_NCT*.md"):
        try:
            content = file_path.read_text()
            
            # Extract metadata from markdown
            nct_id_match = re.search(r'NCT\d+', file_path.name)
            nct_id = nct_id_match.group(0) if nct_id_match else "Unknown"
            
            tier_match = re.search(r'_(TOP_TIER|GOOD_TIER|ACCEPTABLE_TIER)', file_path.name)
            tier = tier_match.group(1) if tier_match else "UNKNOWN"
            
            score_match = re.search(r'\*\*Composite Score\*\*: ([\d.]+)', content)
            match_score = float(score_match.group(1)) if score_match else 0.0
            
            title_match = re.search(r'\*\*Trial Name\*\*: (.+)', content)
            title = title_match.group(1).strip() if title_match else "Unknown Trial"
            
            phase_match = re.search(r'\*\*Phase\*\*: (.+)', content)
            phase = phase_match.group(1).strip() if phase_match else "N/A"
            
            has_llm = "WHY THIS TRIAL IS A GOOD FIT" in content
            
            dossiers.append({
                "nct_id": nct_id,
                "tier": tier,
                "match_score": match_score,
                "title": title,
                "phase": phase,
                "has_llm_analysis": has_llm,
                "file_path": file_path.name
            })
        except Exception as e:
            continue
    
    # Sort by match score descending
    dossiers.sort(key=lambda x: x['match_score'], reverse=True)
    
    return {
        "total": len(dossiers),
        "dossiers": dossiers
    }

@router.get("/detail/{nct_id}")
async def get_dossier_detail(nct_id: str) -> Dict[str, Any]:
    """
    Get full dossier content for a specific trial.
    
    Args:
        nct_id: NCT identifier (e.g., "NCT06619236")
    
    Returns:
        {
            "nct_id": "NCT06619236",
            "markdown": "Full markdown content...",
            "metadata": {
                "tier": "TOP_TIER",
                "match_score": 0.97,
                "generated_at": "2025-11-15T..."
            }
        }
    """
    if not DOSSIER_DIR.exists():
        raise HTTPException(status_code=404, detail="Dossier directory not found")
    
    # Find matching file
    matching_files = list(DOSSIER_DIR.glob(f"INTELLIGENCE_{nct_id}_*.md"))
    
    if not matching_files:
        raise HTTPException(status_code=404, detail=f"Dossier for {nct_id} not found")
    
    file_path = matching_files[0]
    content = file_path.read_text()
    
    # Extract metadata
    tier_match = re.search(r'_(TOP_TIER|GOOD_TIER|ACCEPTABLE_TIER)', file_path.name)
    tier = tier_match.group(1) if tier_match else "UNKNOWN"
    
    score_match = re.search(r'\*\*Composite Score\*\*: ([\d.]+)', content)
    match_score = float(score_match.group(1)) if score_match else 0.0
    
    generated_match = re.search(r'\*\*Generated\*\*: (.+)', content)
    generated_at = generated_match.group(1).strip() if generated_match else "Unknown"
    
    return {
        "nct_id": nct_id,
        "markdown": content,
        "metadata": {
            "tier": tier,
            "match_score": match_score,
            "generated_at": generated_at,
            "file_name": file_path.name
        }
    }

@router.get("/export/{nct_id}")
async def export_dossier(nct_id: str, format: str = "markdown") -> Dict[str, Any]:
    """
    Export dossier in various formats.
    
    Args:
        nct_id: NCT identifier
        format: "markdown" | "html" | "pdf" (future)
    
    Returns:
        {
            "nct_id": "NCT06619236",
            "format": "markdown",
            "content": "...",
            "download_url": "/path/to/file"
        }
    """
    # For now, just return markdown content
    # Future: Generate PDF, HTML versions
    detail = await get_dossier_detail(nct_id)
    
    return {
        "nct_id": nct_id,
        "format": format,
        "content": detail["markdown"],
        "metadata": detail["metadata"]
    }

@router.post("/generate_missing")
async def generate_missing_dossiers() -> Dict[str, Any]:
    """
    Trigger generation of missing dossiers.
    
    This is a webhook to trigger the Python script.
    For MVP, returns instructions for manual execution.
    """
    return {
        "status": "manual_execution_required",
        "command": "cd oncology-coPilot/oncology-backend-minimal && python3 find_trials_FROM_FRESH_TABLE.py",
        "expected_duration_minutes": 90,
        "message": "Run this command to generate all 60 dossiers with LLM analysis"
    }
```

**Register in `api/main.py`**:
```python
from api.routers import ayesha_dossiers

app.include_router(ayesha_dossiers.router)
```

---

#### Task 2.2: Test API Endpoints (15 min)

**Smoke Tests**:
```bash
# List all dossiers
curl http://localhost:8000/api/ayesha/dossiers/list | jq

# Get specific dossier
curl http://localhost:8000/api/ayesha/dossiers/detail/NCT06619236 | jq

# Export dossier
curl http://localhost:8000/api/ayesha/dossiers/export/NCT06619236?format=markdown | jq
```

**Expected**:
- `/list`: Returns 60 dossiers with metadata
- `/detail`: Returns full markdown content
- `/export`: Returns exportable content

---

### **PHASE 3: Frontend Dossier Browser** (3-4 hours)
**Why**: Beautiful UI for Ayesha's oncologist to browse all dossiers

#### Task 3.1: Create Dossier List View (90 min)

**File**: `src/pages/AyeshaDossierBrowser.jsx` (NEW)

```javascript
/**
 * Ayesha Dossier Browser
 * 
 * Browse all 60 trial intelligence dossiers with:
 * - Tier filtering (Top/Good/All)
 * - Search by NCT ID or keywords
 * - Sort by match score
 * - Click to view full dossier
 */
import React, { useState, useEffect } from 'react';
import {
  Box,
  Typography,
  Grid,
  TextField,
  ToggleButtonGroup,
  ToggleButton,
  CircularProgress,
  Alert,
  Paper,
  InputAdornment,
  Chip
} from '@mui/material';
import { SearchIcon, FilterIcon } from '@heroicons/react/24/outline';
import DossierSummaryCard from '../components/ayesha/DossierSummaryCard';

const API_ROOT = import.meta.env.VITE_API_ROOT || 'http://localhost:8000';

const AyeshaDossierBrowser = () => {
  const [dossiers, setDossiers] = useState([]);
  const [filteredDossiers, setFilteredDossiers] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);
  
  const [tierFilter, setTierFilter] = useState('ALL'); // ALL | TOP_TIER | GOOD_TIER
  const [searchQuery, setSearchQuery] = useState('');

  useEffect(() => {
    loadDossiers();
  }, []);

  useEffect(() => {
    applyFilters();
  }, [tierFilter, searchQuery, dossiers]);

  const loadDossiers = async () => {
    setIsLoading(true);
    setError(null);

    try {
      const response = await fetch(`${API_ROOT}/api/ayesha/dossiers/list`);
      
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const data = await response.json();
      setDossiers(data.dossiers || []);
      
    } catch (err) {
      setError(err.message);
    } finally {
      setIsLoading(false);
    }
  };

  const applyFilters = () => {
    let filtered = [...dossiers];

    // Tier filter
    if (tierFilter !== 'ALL') {
      filtered = filtered.filter(d => d.tier === tierFilter);
    }

    // Search filter
    if (searchQuery) {
      const query = searchQuery.toLowerCase();
      filtered = filtered.filter(d => 
        d.nct_id.toLowerCase().includes(query) ||
        d.title.toLowerCase().includes(query)
      );
    }

    setFilteredDossiers(filtered);
  };

  if (isLoading) {
    return (
      <Box display="flex" justifyContent="center" alignItems="center" minHeight="400px">
        <CircularProgress />
      </Box>
    );
  }

  if (error) {
    return (
      <Box p={3}>
        <Alert severity="error">
          <Typography variant="h6">Error Loading Dossiers</Typography>
          <Typography variant="body2">{error}</Typography>
        </Alert>
      </Box>
    );
  }

  return (
    <Box p={3}>
      {/* Header */}
      <Box mb={4}>
        <Typography variant="h4" gutterBottom>
          ⚔️ Ayesha's Clinical Trial Intelligence Reports
        </Typography>
        <Typography variant="body1" color="text.secondary" gutterBottom>
          {dossiers.length} Commander-grade dossiers with mechanistic fit analysis
        </Typography>
        <Chip 
          label={`Patient: AK (Stage IVB HGSOC)`}
          color="primary"
          sx={{ mt: 1 }}
        />
      </Box>

      {/* Filters */}
      <Paper sx={{ p: 2, mb: 3 }}>
        <Grid container spacing={2} alignItems="center">
          {/* Search */}
          <Grid item xs={12} md={6}>
            <TextField
              fullWidth
              size="small"
              placeholder="Search by NCT ID or keywords..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              InputProps={{
                startAdornment: (
                  <InputAdornment position="start">
                    <SearchIcon className="h-5 w-5" />
                  </InputAdornment>
                ),
              }}
            />
          </Grid>

          {/* Tier Filter */}
          <Grid item xs={12} md={6}>
            <Box display="flex" alignItems="center" gap={2}>
              <FilterIcon className="h-5 w-5" />
              <ToggleButtonGroup
                value={tierFilter}
                exclusive
                onChange={(e, value) => value && setTierFilter(value)}
                size="small"
              >
                <ToggleButton value="ALL">
                  All ({dossiers.length})
                </ToggleButton>
                <ToggleButton value="TOP_TIER">
                  Top Tier ({dossiers.filter(d => d.tier === 'TOP_TIER').length})
                </ToggleButton>
                <ToggleButton value="GOOD_TIER">
                  Good Tier ({dossiers.filter(d => d.tier === 'GOOD_TIER').length})
                </ToggleButton>
              </ToggleButtonGroup>
            </Box>
          </Grid>
        </Grid>
      </Paper>

      {/* Results Count */}
      <Typography variant="body2" color="text.secondary" mb={2}>
        Showing {filteredDossiers.length} of {dossiers.length} dossiers
      </Typography>

      {/* Dossier Grid */}
      {filteredDossiers.length === 0 ? (
        <Alert severity="info">
          No dossiers match your filters. Try adjusting your search or tier filter.
        </Alert>
      ) : (
        <Grid container spacing={2}>
          {filteredDossiers.map((dossier, idx) => (
            <Grid item xs={12} key={dossier.nct_id}>
              <DossierSummaryCard 
                dossier={dossier}
                rank={idx + 1}
              />
            </Grid>
          ))}
        </Grid>
      )}
    </Box>
  );
};

export default AyeshaDossierBrowser;
```

**Register Route in `App.jsx`**:
```javascript
import AyeshaDossierBrowser from './pages/AyeshaDossierBrowser';

// In Routes:
<Route path="/ayesha-dossiers" element={<AyeshaDossierBrowser />} />
```

---

#### Task 2.2: Create Dossier Summary Card Component (45 min)

**File**: `src/components/ayesha/DossierSummaryCard.jsx` (NEW)

```javascript
/**
 * Dossier Summary Card
 * 
 * Displays dossier metadata with click-to-view full report.
 */
import React from 'react';
import {
  Card,
  CardContent,
  Typography,
  Box,
  Chip,
  Button,
  LinearProgress,
} from '@mui/material';
import { useNavigate } from 'react-router-dom';
import {
  DocumentTextIcon,
  CheckCircleIcon,
  ClockIcon,
} from '@heroicons/react/24/outline';

const DossierSummaryCard = ({ dossier, rank }) => {
  const navigate = useNavigate();

  const getTierColor = (tier) => {
    switch(tier) {
      case 'TOP_TIER': return 'success';
      case 'GOOD_TIER': return 'info';
      default: return 'default';
    }
  };

  const getTierLabel = (tier) => {
    switch(tier) {
      case 'TOP_TIER': return '⭐ Top Tier';
      case 'GOOD_TIER': return '✅ Good Tier';
      default: return '📋 Acceptable';
    }
  };

  const scorePercent = Math.round(dossier.match_score * 100);

  return (
    <Card sx={{ 
      border: '1px solid',
      borderColor: dossier.tier === 'TOP_TIER' ? 'success.main' : 'grey.300',
      '&:hover': { boxShadow: 3 }
    }}>
      <CardContent>
        <Box display="flex" justifyContent="space-between" alignItems="start">
          {/* Left: Title & Metadata */}
          <Box flex={1}>
            <Box display="flex" alignItems="center" gap={1} mb={1}>
              <Chip
                label={`#${rank}`}
                size="small"
                color="primary"
                sx={{ fontWeight: 'bold' }}
              />
              <Chip
                label={getTierLabel(dossier.tier)}
                size="small"
                color={getTierColor(dossier.tier)}
              />
              {dossier.has_llm_analysis && (
                <Chip
                  icon={<CheckCircleIcon className="h-4 w-4" />}
                  label="LLM Enhanced"
                  size="small"
                  variant="outlined"
                  color="secondary"
                />
              )}
            </Box>

            <Typography variant="h6" gutterBottom>
              {dossier.title}
            </Typography>

            <Box display="flex" gap={1} flexWrap="wrap" mb={2}>
              <Chip label={dossier.nct_id} size="small" variant="outlined" />
              <Chip label={dossier.phase} size="small" variant="outlined" />
            </Box>
          </Box>

          {/* Right: Match Score */}
          <Box textAlign="right" minWidth="120px">
            <Typography variant="h4" color={scorePercent >= 80 ? 'success.main' : 'text.primary'}>
              {scorePercent}%
            </Typography>
            <Typography variant="caption" color="text.secondary">
              Match Score
            </Typography>
            <LinearProgress
              variant="determinate"
              value={scorePercent}
              color={scorePercent >= 80 ? 'success' : 'primary'}
              sx={{ mt: 1, height: 8, borderRadius: 4 }}
            />
          </Box>
        </Box>

        {/* Actions */}
        <Box display="flex" gap={2} mt={2}>
          <Button
            variant="contained"
            startIcon={<DocumentTextIcon className="h-5 w-5" />}
            onClick={() => navigate(`/ayesha-dossiers/${dossier.nct_id}`)}
            size="small"
          >
            View Full Dossier
          </Button>
          <Button
            variant="outlined"
            href={`https://clinicaltrials.gov/study/${dossier.nct_id}`}
            target="_blank"
            size="small"
          >
            ClinicalTrials.gov
          </Button>
        </Box>
      </CardContent>
    </Card>
  );
};

export default DossierSummaryCard;
```

---

#### Task 2.3: Create Dossier Detail View (60 min)

**File**: `src/pages/AyeshaDossierDetail.jsx` (NEW)

```javascript
/**
 * Ayesha Dossier Detail View
 * 
 * Full dossier display with:
 * - Markdown rendering
 * - Export options (PDF, share link)
 * - Back to list navigation
 */
import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import {
  Box,
  Typography,
  Paper,
  Button,
  CircularProgress,
  Alert,
  Breadcrumbs,
  Link,
} from '@mui/material';
import {
  ArrowLeftIcon,
  ShareIcon,
  DocumentArrowDownIcon,
} from '@heroicons/react/24/outline';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const API_ROOT = import.meta.env.VITE_API_ROOT || 'http://localhost:8000';

const AyeshaDossierDetail = () => {
  const { nct_id } = useParams();
  const navigate = useNavigate();
  
  const [dossier, setDossier] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    loadDossier();
  }, [nct_id]);

  const loadDossier = async () => {
    setIsLoading(true);
    setError(null);

    try {
      const response = await fetch(
        `${API_ROOT}/api/ayesha/dossiers/detail/${nct_id}`
      );
      
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const data = await response.json();
      setDossier(data);
      
    } catch (err) {
      setError(err.message);
    } finally {
      setIsLoading(false);
    }
  };

  const handleExport = async (format = 'markdown') => {
    try {
      const response = await fetch(
        `${API_ROOT}/api/ayesha/dossiers/export/${nct_id}?format=${format}`
      );
      const data = await response.json();
      
      // Download as file
      const blob = new Blob([data.content], { type: 'text/markdown' });
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${nct_id}_DOSSIER.md`;
      a.click();
      window.URL.revokeObjectURL(url);
    } catch (err) {
      console.error('Export failed:', err);
    }
  };

  const handleShare = () => {
    // Copy dossier URL to clipboard
    const url = window.location.href;
    navigator.clipboard.writeText(url);
    alert('Dossier link copied to clipboard!');
  };

  if (isLoading) {
    return (
      <Box display="flex" justifyContent="center" alignItems="center" minHeight="400px">
        <CircularProgress />
      </Box>
    );
  }

  if (error) {
    return (
      <Box p={3}>
        <Alert severity="error">
          <Typography variant="h6">Error Loading Dossier</Typography>
          <Typography variant="body2">{error}</Typography>
        </Alert>
        <Button
          startIcon={<ArrowLeftIcon className="h-5 w-5" />}
          onClick={() => navigate('/ayesha-dossiers')}
          sx={{ mt: 2 }}
        >
          Back to List
        </Button>
      </Box>
    );
  }

  return (
    <Box p={3}>
      {/* Breadcrumbs */}
      <Breadcrumbs sx={{ mb: 2 }}>
        <Link
          component="button"
          variant="body2"
          onClick={() => navigate('/ayesha-trials')}
          underline="hover"
        >
          Ayesha Trial Explorer
        </Link>
        <Link
          component="button"
          variant="body2"
          onClick={() => navigate('/ayesha-dossiers')}
          underline="hover"
        >
          All Dossiers
        </Link>
        <Typography variant="body2" color="text.primary">
          {dossier?.nct_id}
        </Typography>
      </Breadcrumbs>

      {/* Action Bar */}
      <Box display="flex" justifyContent="space-between" alignItems="center" mb={3}>
        <Button
          startIcon={<ArrowLeftIcon className="h-5 w-5" />}
          onClick={() => navigate('/ayesha-dossiers')}
          variant="outlined"
        >
          Back to List
        </Button>
        
        <Box display="flex" gap={2}>
          <Button
            startIcon={<ShareIcon className="h-5 w-5" />}
            onClick={handleShare}
            variant="outlined"
          >
            Share Link
          </Button>
          <Button
            startIcon={<DocumentArrowDownIcon className="h-5 w-5" />}
            onClick={() => handleExport('markdown')}
            variant="contained"
          >
            Export Dossier
          </Button>
        </Box>
      </Box>

      {/* Dossier Metadata */}
      <Paper sx={{ p: 2, mb: 3, bgcolor: 'grey.50' }}>
        <Grid container spacing={2}>
          <Grid item xs={12} sm={4}>
            <Typography variant="caption" color="text.secondary">NCT ID</Typography>
            <Typography variant="body1" fontWeight="bold">{dossier.nct_id}</Typography>
          </Grid>
          <Grid item xs={12} sm={4}>
            <Typography variant="caption" color="text.secondary">Match Tier</Typography>
            <Typography variant="body1" fontWeight="bold">
              <Chip 
                label={dossier.metadata.tier}
                size="small"
                color={getTierColor(dossier.metadata.tier)}
              />
            </Typography>
          </Grid>
          <Grid item xs={12} sm={4}>
            <Typography variant="caption" color="text.secondary">Match Score</Typography>
            <Typography variant="body1" fontWeight="bold">
              {Math.round(dossier.metadata.match_score * 100)}%
            </Typography>
          </Grid>
        </Grid>
      </Paper>

      {/* Markdown Content */}
      <Paper sx={{ p: 4 }}>
        <ReactMarkdown
          remarkPlugins={[remarkGfm]}
          components={{
            // Custom renderers for better formatting
            h1: ({node, ...props}) => <Typography variant="h3" gutterBottom {...props} />,
            h2: ({node, ...props}) => <Typography variant="h4" gutterBottom sx={{ mt: 3 }} {...props} />,
            h3: ({node, ...props}) => <Typography variant="h5" gutterBottom sx={{ mt: 2 }} {...props} />,
            p: ({node, ...props}) => <Typography variant="body1" paragraph {...props} />,
            table: ({node, ...props}) => (
              <Box sx={{ overflowX: 'auto', my: 2 }}>
                <table style={{ width: '100%', borderCollapse: 'collapse' }} {...props} />
              </Box>
            ),
            th: ({node, ...props}) => (
              <th style={{ 
                border: '1px solid #ddd', 
                padding: '12px', 
                backgroundColor: '#f5f5f5',
                textAlign: 'left'
              }} {...props} />
            ),
            td: ({node, ...props}) => (
              <td style={{ 
                border: '1px solid #ddd', 
                padding: '12px' 
              }} {...props} />
            ),
          }}
        >
          {dossier.markdown}
        </ReactMarkdown>
      </Paper>

      {/* Bottom Actions */}
      <Box display="flex" justifyContent="space-between" mt={3}>
        <Button
          startIcon={<ArrowLeftIcon className="h-5 w-5" />}
          onClick={() => navigate('/ayesha-dossiers')}
          variant="outlined"
        >
          Back to List
        </Button>
        
        <Button
          startIcon={<DocumentArrowDownIcon className="h-5 w-5" />}
          onClick={() => handleExport('markdown')}
          variant="contained"
        >
          Download Dossier
        </Button>
      </Box>
    </Box>
  );
};

// Helper function (move to top of file)
const getTierColor = (tier) => {
  switch(tier) {
    case 'TOP_TIER': return 'success';
    case 'GOOD_TIER': return 'info';
    default: return 'default';
  }
};

export default AyeshaDossierDetail;
```

---

#### Task 2.4: Update Navigation (15 min)

**File**: `src/constants/index.js`

```javascript
// Add to sidebar links
{
  name: "Ayesha: Dossier Browser",
  path: "/ayesha-dossiers",
  icon: DocumentTextIcon,
  category: "Ayesha Care"
}
```

---

### **PHASE 4: Sharing with Oncologist** (30 min)
**Why**: Enable direct sharing with care team

#### Task 4.1: Create Shareable Link Generator (15 min)

**Enhancement to `AyeshaDossierDetail.jsx`**:
```javascript
const generateShareableLink = () => {
  const baseUrl = window.location.origin;
  const shareUrl = `${baseUrl}/ayesha-dossiers/${nct_id}`;
  
  // Create email template
  const subject = `Clinical Trial: ${nct_id} - Intelligence Report for AK`;
  const body = `
Dear Dr. [Oncologist Name],

Please review this trial intelligence report for AK (Stage IVB HGSOC):

Trial: ${dossier.nct_id}
Match Score: ${Math.round(dossier.metadata.match_score * 100)}%
Tier: ${dossier.metadata.tier}

Full Report: ${shareUrl}

This dossier includes:
- Eligibility assessment with probability calculations
- Mechanistic fit analysis (drug vs disease)
- Location logistics (NYC metro)
- Risk-benefit analysis for Ayesha's specific case
- Comparison to standard of care

Please let me know your thoughts.

Best regards,
AK
  `.trim();

  const mailtoUrl = `mailto:?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(body)}`;
  window.location.href = mailtoUrl;
};

// Add button in UI
<Button
  startIcon={<EnvelopeIcon className="h-5 w-5" />}
  onClick={generateShareableLink}
  variant="outlined"
>
  Email to Oncologist
</Button>
```

#### Task 4.2: Batch Export (15 min)

**Enhancement to `AyeshaDossierBrowser.jsx`**:
```javascript
const handleBatchExport = async (tier = 'TOP_TIER') => {
  const dossiersToExport = dossiers.filter(d => d.tier === tier);
  
  for (const dossier of dossiersToExport) {
    try {
      const response = await fetch(
        `${API_ROOT}/api/ayesha/dossiers/export/${dossier.nct_id}?format=markdown`
      );
      const data = await response.json();
      
      // Download each file
      const blob = new Blob([data.content], { type: 'text/markdown' });
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${dossier.nct_id}_DOSSIER.md`;
      a.click();
      window.URL.revokeObjectURL(url);
      
      // Rate limit downloads
      await new Promise(resolve => setTimeout(resolve, 500));
    } catch (err) {
      console.error(`Failed to export ${dossier.nct_id}:`, err);
    }
  }
  
  alert(`Exported ${dossiersToExport.length} ${tier} dossiers`);
};

// Add button in UI
<Button
  startIcon={<DocumentArrowDownIcon className="h-5 w-5" />}
  onClick={() => handleBatchExport('TOP_TIER')}
  variant="contained"
>
  Export All Top-Tier Dossiers
</Button>
```

---

## 📋 COMPLETE EXECUTION CHECKLIST

### **PHASE 1: Generate Dossiers** (90 min)
- [ ] Fix rate limiting (30s between LLM calls) ✅ DONE
- [ ] Run `find_trials_FROM_FRESH_TABLE.py`
- [ ] Verify 60 markdown files created
- [ ] Spot-check 3-5 dossiers for LLM quality

### **PHASE 2: Backend API** (2 hours)
- [ ] Create `api/routers/ayesha_dossiers.py` (60 min)
- [ ] Register router in `api/main.py` (5 min)
- [ ] Test `/list` endpoint (15 min)
- [ ] Test `/detail/{nct_id}` endpoint (15 min)
- [ ] Test `/export/{nct_id}` endpoint (15 min)
- [ ] Create smoke test script (15 min)

### **PHASE 3: Frontend** (3-4 hours)
- [ ] Create `AyeshaDossierBrowser.jsx` page (90 min)
- [ ] Create `DossierSummaryCard.jsx` component (45 min)
- [ ] Create `AyeshaDossierDetail.jsx` page (60 min)
- [ ] Add routes to `App.jsx` (15 min)
- [ ] Update navigation in `constants/index.js` (15 min)
- [ ] Install `react-markdown` + `remark-gfm` (5 min)
- [ ] Test end-to-end (30 min)

### **PHASE 4: Sharing** (30 min)
- [ ] Add email template generator (15 min)
- [ ] Add batch export functionality (15 min)
- [ ] Test share/export flows (10 min)

---

## 📦 REQUIRED NPM PACKAGES

```bash
cd oncology-coPilot/oncology-frontend
npm install react-markdown remark-gfm
```

**Why**:
- `react-markdown`: Render markdown in React
- `remark-gfm`: GitHub-flavored markdown (tables, checklists)

---

## 🎯 USER EXPERIENCE FLOW

### **Oncologist's Journey**:
1. **Navigate to `/ayesha-dossiers`** → See all 60 trials ranked by score
2. **Filter to "Top Tier"** → See 60 highest-scoring trials
3. **Click "View Full Dossier"** → Read complete intelligence report
4. **Export Dossier** → Download markdown file for records
5. **Email to Colleague** → Share link with other physicians

### **Ayesha's Journey**:
1. **Browse all trials** → Understand full landscape
2. **Filter by tier** → Focus on best matches
3. **Read LLM analysis** → Understand "Why this trial fits me"
4. **Share with family** → Send link to family members
5. **Discuss with oncologist** → Review together during appointment

---

## ⏰ TOTAL TIMELINE

| Phase | Duration | Owner | Status |
|-------|----------|-------|--------|
| 1. Generate Dossiers | 90 min | Zo | ⏸️ Ready to execute |
| 2. Backend API | 2 hours | Zo | ⏸️ Waiting on Phase 1 |
| 3. Frontend UI | 3-4 hours | JR (or Zo) | ⏸️ Waiting on Phase 2 |
| 4. Sharing Features | 30 min | JR (or Zo) | ⏸️ Waiting on Phase 3 |
| **TOTAL** | **6-7 hours** | | |

---

## 🚀 IMMEDIATE NEXT STEPS

### **Option A: Full Build** (6-7 hours)
Execute all 4 phases → Complete dossier browser

### **Option B: Quick MVP** (3 hours)
- Phase 1: Generate all 60 dossiers (90 min)
- Phase 2: Backend API only (2 hours)
- **Skip Phase 3-4** (use API directly, or browse markdown files manually)

### **Option C: Delegate to JR** (Parallel execution)
- **Zo**: Phase 1 (generate dossiers) + Phase 2 (backend API)
- **JR**: Phase 3 (frontend UI) + Phase 4 (sharing)
- **Timeline**: 2-3 hours (parallel)

---

## 🎯 MY RECOMMENDATION: **Option C (Parallel Execution)**

**Why**:
- ✅ **Fastest** - 2-3 hours vs 6-7 hours
- ✅ **Efficient** - Parallel work streams
- ✅ **Complete** - Full UI + sharing ready

**Execution**:
1. **Zo (NOW)**: Start Phase 1 (generate 60 dossiers)
2. **While running**: Zo creates backend API (Phase 2)
3. **JR (Parallel)**: Builds frontend (Phase 3) while dossiers generate
4. **Both**: Test integration, add sharing features

---

## 📊 DELIVERABLES FOR AYESHA & ONCOLOGIST

### **What They Get**:
1. ✅ **60 Trial Intelligence Reports** - Complete landscape
2. ✅ **Beautiful Browser UI** - Filter, search, browse
3. ✅ **Full Dossier Views** - Commander-grade analysis
4. ✅ **Export Functionality** - Download for records
5. ✅ **Share Links** - Email to care team
6. ✅ **Mobile Responsive** - Review on phone/tablet

### **Value Delivered**:
- **Comprehensive**: All recruiting ovarian trials (NYC accessible)
- **Actionable**: Clear eligibility gates, requirements, next steps
- **Transparent**: LLM reasoning explains "why this trial fits"
- **Professional**: Ready to share with oncologist same day

---

## 🔥 COMMANDER'S VERDICT

**Execute Option C**: Parallel execution with Zo + JR

**First Command** (Zo executes NOW):
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal
python3 find_trials_FROM_FRESH_TABLE.py
```

**This generates all 60 dossiers in 90 minutes**

**While running**: Zo builds backend API (Phase 2)

**MISSION STATUS**: ⚔️ **READY TO EXECUTE - AWAITING COMMANDER'S GO** ⚔️






