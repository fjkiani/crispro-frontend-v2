---
alwaysApply: false
description: Evidence Intelligence – Execution Doctrine (Sept 2025). Plain‑language website copy + execution plan for in‑silico evidence panels (confidence, tiers, badges, citations, provenance) with sublinks and JSON schemas (RUO).
globs: 
---


export const coPilotDetailsData: Record<string, CoPilotDetailContent> = {
  "evidence": {
    slug: "evidence",
    pageTitle: "Evidence Intelligence: Confidence, Tiers, Badges, Citations",
    heroSubtitle: "Turn raw findings into a clear evidence story: confidence, tier, badges, and citations — all with provenance (RUO).",
    vision: "Make decisions easier by showing exactly how strong the evidence is — and why — in one view you can trust and share.",

    // Website value props (plain)
    valueProps: [
      {
        audience: 'For Tumor Boards',
        icon: 'ShieldCheck',
        points: [
          'Confidence (0–1) and an evidence tier you can read.',
          'Badges that explain strength: Guideline, RCT, ClinVar‑Strong, Pathway‑Aligned.',
          'Citations and provenance on every result (RUO).'
        ]
      },
      {
        audience: 'For Researchers',
        icon: 'BookOpen',
        points: [
          'ClinVar priors, literature tabs, and cohort overlays (optional).',
          'Caching and provider fallbacks for stability.',
          'Exportable panels you can reuse in notes.'
        ]
      }
    ],

    buildsOn: "What this runs on (today vs roadmap)",
    buildsOnStackPoints: [
      "**Today:** ClinVar prior + evidence panel with confidence/tier/badges/citations.",
      "**Optional today:** Literature (PubMed/OpenAlex/Semantic Scholar) with caching/fallbacks.",
      "**Roadmap:** Guideline/on‑label cues and cohort evidence overlays by disease."
    ],

    kpis: [
      { label: 'ClinVar AUROC (total n=53,210)', value: '0.957' },
      { label: 'SpliceVarDB AUROC (n=4,950)', value: '~0.825–0.826' },
      { label: 'VUS target', value: '→ 15% (≈$2.1M saved/program)' }
    ],

    observedOutcomes: [
      'Tier promotions when ClinVar‑Strong and Pathway‑Aligned co‑occur (~10–20% Consider→Supported)',
      'Confidence +0.05–0.12 with literature on and cohort overlays aligned',
      'Fewer reworks: shareable panel with run ID speeds alignment'
    ],

    genomicInsightsOverview: "The Evidence panel shows confidence (0–1), tier (Supported/Consider/Insufficient), badges (Guideline/RCT/ClinVar‑Strong/Pathway‑Aligned), citations, and provenance. Tabs let you dive into ClinVar, Literature, and Cohort evidence without losing the simple story.",
    coreProblemIntro: "Raw results are hard to trust. We turn them into a concise evidence story you can read and share.",
    coreProblemPoints: [
      "No single place to see strength and sources.",
      "Hard to explain why a result deserves action.",
      "Difficult to reproduce what was shown last time."
    ],

    genomicUseCasesGrid: [
      { label: "Evidence overview panel", iconName: "Layout", color: "text-blue-400" },
      { label: "ClinVar prior tab", iconName: "Database", color: "text-green-400" },
      { label: "Literature tab (optional)", iconName: "BookOpen", color: "text-purple-400" },
      { label: "Cohort evidence tab (optional)", iconName: "Users", color: "text-orange-400" },
      { label: "Citations export", iconName: "FileText", color: "text-cyan-400" },
      { label: "Provenance (run ID)", iconName: "Hash", color: "text-pink-400" }
    ],

    keyCapabilities: [
      {
        title: "Evidence Overview (live)",
        technical: "Render confidence/tier/badges and citations with provenance. Tier is derived from rules; badges reflect strength sources (Guideline/RCT/ClinVar‑Strong/Pathway‑Aligned).",
        scientific: "Presents a compact strength summary with transparent sources and methods (RUO).",
        business: `
- **Clarity:** Teams align faster when strength is explicit.
- **Audit:** Sources and run IDs reduce rework.
`,
        genomicUseCasesParagraph: "Today: \n1. **Evidence panel** on therapy views. \n2. **Export citations** and show provenance."
      },
      {
        title: "ClinVar Prior (live)",
        technical: "ClinVar lookup via `/api/evidence/clinvar` with caching; summarize relevant assertions and review status.",
        scientific: "Provides a prior for the variant/gene that feeds the tier and badges.",
        business: `
- **Trust:** Recognized sources with concise summarization.
`,
        genomicUseCasesParagraph: "Today: \n1. **ClinVar tab** summarizes assertions and review status."
      },
      {
        title: "Literature (optional)",
        technical: "Provider fallback (PubMed/OpenAlex/Semantic Scholar) with Redis caching, dedupe by title/DOI, and timeouts/retries.",
        scientific: "Adds study context and MoA‑aware relevance without blocking decisions when rates limit.",
        business: `
- **Context:** Stronger story when references are available.
`,
        genomicUseCasesParagraph: "Today (optional): \n1. **Literature tab** with compact relevance and links."
      },
      {
        title: "Cohort Evidence (optional)",
        technical: "When Cohort Lab provides overlays, display a small prevalence/metric snippet tied to the gene/variant or pathway.",
        scientific: "Grounds the result in real‑world context (RUO).",
        business: `
- **Confidence:** Supplemental lift when context aligns.
`,
        genomicUseCasesParagraph: "Today (when present): \n1. **Cohort tab** with a compact card (n, metric, study)."
      }
    ],

    valuePropositionSections: [
      {
        audience: "For the Care Team",
        points: [
          "A single place to see strength and sources.",
          "Short, shareable evidence story with citations (RUO).",
          "Run IDs for repeatability and review."
        ]
      }
    ],

    conclusion: "Evidence that’s simple to read and easy to trust. Confidence. Tier. Badges. Citations. Provenance (RUO)."
  },
};


# Evidence Intelligence – Execution Doctrine (Sept 2025)

This rule describes how to deliver the Evidence panel (overview + sublinks) with transparent provenance and stable provider behavior.

## What we can deliver today (research‑grade)
- An Evidence panel with:
  - **Overview:** confidence (0–1), tier (Supported/Consider/Insufficient), badges, citations, provenance (run_id, profile).
  - **Tabs:** ClinVar, Literature (optional), Cohort (optional), Citations, Provenance.

## Numbers That Matter (Research‑Mode)
- Variant foundations: ClinVar AUROC 0.957 (n=53,210) support reliable priors; SpliceVarDB ~0.825–0.826 (n=4,950) improves regulatory context.
- Business target: VUS ~40% → ~15% coupled with evidence tiering and badges; ≈$2.1M/program saved by focusing follow‑ups.
- Ops: caching and fallback reduce latency and flakiness; exportable citations speed documentation.

## Observed Outcomes (pilot; research‑mode)
- Consider→Supported promotions ~10–20% when ClinVar‑Strong + Pathway‑Aligned and literature present.
- Confidence lifts +0.05–0.12 with aligned literature/cohort overlays.
- Faster reviews: exportable citations + run IDs reduced rework and meeting time.

## Plain English: Why this matters
- Decisions need strength and sources. We show both, clearly and consistently, with transparent provenance.

## Where the code lives (reused components)
- Evidence router: `[evidence.py](mdc:oncology-coPilot/oncology-backend-minimal/api/routers/evidence.py)`
- Efficacy orchestrator (panel consumer): `[efficacy.py](mdc:oncology-coPilot/oncology-backend-minimal/api/routers/efficacy.py)`
- Cohort Lab (optional overlays): `[datasets.py](mdc:oncology-coPilot/oncology-backend-minimal/api/routers/datasets.py)`

## How to run Evidence now
- ClinVar prior:
```bash
curl -sS "http://127.0.0.1:8000/api/evidence/clinvar?chrom=7&pos=140453136&ref=T&alt=A"
```

## Methodology (simple, transparent)
- **Tier rules (display):** Supported when meaningful literature or ClinVar‑Strong + pathway alignment; Consider for partial; Insufficient when sparse.
- **Badges:** Guideline, RCT, ClinVar‑Strong, Pathway‑Aligned.
- **Providers:** PubMed, OpenAlex, Semantic Scholar via fallback; Redis caching; dedupe by title/DOI; timeouts/retries.

## Frontend integration (tabs and contracts)
- Tabs: Overview | ClinVar | Literature | Cohort | Citations | Provenance
- Render Overview always; enable Literature/Cohort when providers/overlays present.

## Subcomponents (FE spec – keep it all here)
- EvidencePanel (container)
  - Purpose: orchestrates tabs; fetches data; renders Overview + sublinks
  - Props: evidence (overview), tabsEnabled {clinvar, literature, cohort}, provenance
- OverviewCard
  - Children: ConfidenceBar (0–1), TierLabel (+ TierRulesPopover), EvidenceBadges, CitationCount, ProvenanceMini
  - Props: confidence, tier, badges[], citationsCount, provenance
- ClinVarTab
  - Children: ClinVarSummary (review_status, assertions, accession), VariantLocusChip (chrom, pos, ref, alt)
  - Props: variant, summary {review_status, assertions[], accession}
- LiteratureTab (optional)
  - Children: ProviderStatusChip (PubMed/OpenAlex/S2), LiteratureItemCard (title, year, authors, link, provider, relevance), ExportCitationsButton
  - Props: providers[], items[], cache {ttl_hours, hit}
- CohortTab (optional)
  - Children: CohortSnippetCard (study, n, prevalence, metrics), MetricsMini (AUROC/AUPRC)
  - Props: study, metrics?, snippet? {gene|variant, n, prevalence}
- CitationsTab
  - Children: CitationTable (title, doi/url, provider, year), ExportMenu (CSV/JSON)
  - Props: items[], count
- ProvenanceTab
  - Children: RunInfo (runId, profile, timestamp), FeatureFlagsList (disableLiterature, disableFusion, profile)
  - Props: run_id, profile, flags{}

### Cross‑cutting building blocks
- ConfidenceBar(value, size, showValue?)
- TierRulesPopover(tier, rules)
- EvidenceBadges(badges[], tooltipMap, onClick?)
- ProviderStatusChip(name, enabled, healthy)
- EvidenceSkeleton / EvidenceEmpty / EvidenceError(code, retry)
- ProvenanceMini(runId, profile, timestamp)

## JSON Schemas (drop‑in)
```json
{
  "EvidencePanel": {
    "confidence": 0.52,
    "tier": "Consider",
    "badges": ["Pathway-Aligned"],
    "citations": [
      {"provider": "PubMed", "title": "BRAF pathway in MM", "year": 2021, "url": "https://..."}
    ],
    "provenance": {"run_id": "abc-123", "profile": "Baseline", "timestamp": "2025-09-01T10:12:00Z"}
  },
  "ClinVarTab": {
    "variant": {"chrom": "7", "pos": 140453136, "ref": "T", "alt": "A"},
    "summary": {"review_status": "criteria provided, multiple submitters", "assertions": ["Likely pathogenic"], "accession": "RCV..."}
  },
  "LiteratureTab": {
    "providers": ["PubMed", "OpenAlex"],
    "items": [
      {"title": "...", "year": 2023, "authors": ["..."], "url": "https://...", "source": "PubMed"}
    ],
    "cache": {"ttl_hours": 48, "hit": true}
  },
  "CohortTab": {
    "study": "TCGA-OV",
    "metrics": {"auroc": 0.50, "auprc": 0.50},
    "snippet": {"gene": "BRAF", "n": 42, "prevalence": 0.07}
  },
  "CitationsTab": {
    "count": 3,
    "items": [
      {"title": "...", "doi": "10.xxxx/yyy", "provider": "OpenAlex", "url": "https://..."}
    ]
  },
  "ProvenanceTab": {
    "run_id": "abc-123",
    "profile": "Baseline",
    "flags": {"disableLiterature": true}
  },
  "TierRules": {
    "supported": "meaningful literature or ClinVar-Strong + pathway alignment",
    "consider": "partial support or mixed signals",
    "insufficient": "sparse or conflicting evidence"
  },
  "ProviderStatus": {
    "PubMed": {"enabled": true, "healthy": true},
    "OpenAlex": {"enabled": true, "healthy": true},
    "SemanticScholar": {"enabled": false, "healthy": false}
  }
}
```

## Navigation / sublinks (site IA)
- Evidence (overview)
- ClinVar (prior)
- Literature (optional)
- Cohort (optional)
- Citations (export)
- Provenance (run details)

## Current limitations (transparent)
- Some providers may be unavailable in dev; keep panel functional with priors and badges.
- Literature coverage varies; caching reduces flakiness but not completeness.

## Success criteria
- Overview always renders with confidence/tier/badges/citations/provenance.
- Tabs load without blocking; export of citations works.
- RUO visible; tooltips explain tiers and badges.

