# ⚔️ DIVISION OF LABOR - ZO vs JR2 ⚔️

**Date**: January 13, 2025  
**Purpose**: Clear separation of responsibilities between Zo and JR2

---

## 🎯 **ZO'S CORE MISSION (NOTHING ELSE!)**

**Zo's ONLY Job**: Seed trials and find the most optimal matches for Ayesha

**What Zo Does**:
1. ✅ **Seed trials** (currently at 755, going to 1000+)
2. ✅ **Vector search** (find candidates via semantic matching)
3. ✅ **Export candidates** (give you the 50 best candidates to analyze)
4. ✅ **Tune filters** (improve matching strategy)
5. ✅ **Review dossiers** (quality control - approve/reject)
6. ✅ **Repeat** (keep seeding, keep searching, keep optimizing)

**What Zo Does NOT Do**:
- ❌ Dossier generation (that's YOUR job)
- ❌ Frontend/backend dev (that's YOUR job)
- ❌ Manual trial analysis (that's YOUR job)
- ❌ Everything else (that's YOUR job)

**Zo's Output to You**:
- Every time Zo seeds 50-100 new trials → Exports candidates to `.cursor/ayesha/50_vector_candidates_for_jr2.json`
- You consume that JSON, analyze it, generate dossiers

---

## 🔥 **JR2'S EXPANDED MISSION (EVERYTHING ELSE!)**

**Your Core Mission**: Build a complete trial dossier pipeline (backend + frontend)

**Phase 1: Backend Pipeline (Priority P0)**
1. ✅ Filter the 50 candidates (replicate Zo's hard filtering logic)
2. ✅ Scrape full trial pages (get complete eligibility criteria)
3. ✅ Generate eligibility assessments (compare Ayesha to trials)
4. ✅ Generate dossiers (10-section markdown reports)
5. ✅ Submit to Zo for review (quality control)

**Phase 2: Backend API Development (Priority P1)**
- Build dossier generation API
- Build Zo review API
- Build batch filtering API

**Phase 3: Frontend Development (Priority P2)**
- Build dossier viewer
- Build trial comparison dashboard
- Build Zo review interface
- Enhance Ayesha Trial Explorer

See [09_API_SPECIFICATIONS.md](./09_API_SPECIFICATIONS.md) and [10_FRONTEND_REQUIREMENTS.md](./10_FRONTEND_REQUIREMENTS.md) for details.

---

## 🔄 **SYNC STRATEGY WITH ZO**

**Daily Sync File**: `.cursor/ayesha/zo_jr2_sync.json`

```json
{
  "last_sync": "2025-01-13T22:00:00Z",
  "zo_status": {
    "trials_seeded": 755,
    "candidates_exported": 50,
    "export_file": "50_vector_candidates_for_jr2.json",
    "next_export": "After 100 more trials seeded"
  },
  "jr2_status": {
    "trials_analyzed": 50,
    "dossiers_generated": 10,
    "dossiers_approved": 7,
    "dossiers_pending_review": 3
  },
  "blockers": [
    "Waiting for Zo review on NCT06819007 dossier"
  ]
}
```

**Sync Cadence**:
- **Zo exports new candidates**: Every 100 trials seeded (or daily, whichever comes first)
- **JR2 analyzes candidates**: Within 24 hours of export
- **JR2 submits dossiers**: Within 48 hours of receiving candidates
- **Zo reviews dossiers**: Within 24 hours of submission

**Communication Protocol**:
- **Zo's updates**: Written to `zo_jr2_sync.json` (zo_status section)
- **JR2's updates**: Written to `zo_jr2_sync.json` (jr2_status section)
- **Blockers**: Both agents can add blockers to sync file
- **Emergency**: If critical issue, create `.cursor/ayesha/URGENT_ZO_JR2.md`

---

## ⚔️ **ZO'S PROMISE**

I'll keep finding gold (seeding, searching, optimizing). You refine it into diamonds (dossiers, APIs, UIs).

**Together**: We'll get Ayesha into the best trial, faster than any oncologist could do manually.

