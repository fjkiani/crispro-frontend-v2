import { UseCaseMetrics } from '../metrics/types';
import { discriminativeMetrics, generativeMetrics, businessMetrics } from '../metrics/core-metrics';

export const multipleMyelomaUseCase: UseCaseMetrics = {
  useCaseId: 'multiple-myeloma',
  title: 'Multiple Myeloma Digital Twin',
  description: 'In-silico guidance you can read, trust, and share (research-mode). From variants to therapy ranking and trials in minutes — with confidence, evidence, and provenance.',
  metrics: {
    discriminative: discriminativeMetrics,
    generative: generativeMetrics,
    business: businessMetrics,
    validation: [
      {
        title: 'ClinVar AUROC (total n=53,210)',
        value: { value: 0.957, format: 'decimal', precision: 3 },
        description: 'Overall accuracy across coding/non-coding SNVs and non-SNVs',
        dataset: 'ClinVar',
        sampleSize: 53210,
        source: 'Internal benchmark rollup',
        category: 'validation'
      },
      {
        title: 'SpliceVarDB AUROC (n=4,950)',
        value: { value: 0.826, format: 'decimal', precision: 3 },
        description: 'Exonic/intronic splice prediction accuracy (~82.5–82.6%)',
        dataset: 'SpliceVarDB',
        sampleSize: 4950,
        source: 'SpliceVarDB benchmark',
        category: 'validation'
      },
      {
        title: 'BRCA1 Supervised (coding SNV)',
        value: { value: 0.94, format: 'decimal', precision: 2 },
        description: 'AUROC 0.94, AUPRC 0.84 — oncology benchmark',
        dataset: 'BRCA1',
        sampleSize: 3893,
        source: 'Oncology benchmarks',
        category: 'validation'
      }
    ]
  },
  whyItMatters: [
    'Reduce VUS from ~40% to ~15% (target) to unblock decisions and experiments.',
    'Explainable therapy ranking with citations speeds tumor board alignment.',
    'Provenance (run IDs) ensures repeatability and auditability.'
  ],
  delivered: [
    'Variant Insight chips; Therapy Fit table; Pathway View; Trials shortlist.',
    'Toxicity Risk chip; CRISPR Readiness (demo).',
    'Exports (JSON/CSV) and provenance across views.'
  ],
  howToRead: [
    'Confidence (0–1) is a certainty hint; Evidence Tier is Supported/Consider/Insufficient.',
    'Badges show strength (Guideline, RCT, ClinVar-Strong, Pathway-Aligned).',
    'Fusion labeled when eligible; Baseline remains deterministic.'
  ],
  specificFindings: [
    {
      title: 'MM Research Signals (Observed)',
      description: 'Specific findings from Multiple Myeloma research applications',
      metrics: [
        {
          title: 'WIWFM Confidence (BRAF V600E)',
          value: { value: 0.48, format: 'decimal', precision: 2 },
          description: 'Will-It-Work-For-Me confidence range for BRAF V600E variants',
          dataset: 'MM Research',
          sampleSize: 50,
          source: 'MM research observations',
          category: 'validation'
        },
        {
          title: 'Efficacy Score Range',
          value: { value: 0.22, format: 'decimal', precision: 2 },
          description: 'Efficacy score range (0.17-0.26) for MM variants',
          dataset: 'MM Research',
          sampleSize: 50,
          source: 'MM research observations',
          category: 'validation'
        },
        {
          title: 'Fusion Coverage Usage',
          value: { value: 100, format: 'percentage', precision: 0 },
          description: 'Fusion profile used only when AlphaMissense coverage exists',
          dataset: 'MM Research',
          sampleSize: 50,
          source: 'MM research observations',
          category: 'validation'
        }
      ]
    },
    {
      title: 'Two-Hit Hypothesis (MM)',
      description: 'Multiple Myeloma follows a two-hit model with driver and cooperating alterations',
      metrics: [
        {
          title: 'MAPK Driver Frequency',
          value: { value: 60, format: 'percentage', precision: 0 },
          description: 'Frequency of MAPK pathway activation (BRAF/NRAS/KRAS)',
          dataset: 'MM Genomics',
          sampleSize: 1000,
          source: 'MM genomic studies',
          category: 'validation'
        },
        {
          title: 'TP53/17p Cooperation',
          value: { value: 25, format: 'percentage', precision: 0 },
          description: 'Frequency of TP53/17p alterations as cooperating hits',
          dataset: 'MM Genomics',
          sampleSize: 1000,
          source: 'MM genomic studies',
          category: 'validation'
        },
        {
          title: 'MYC Amplification',
          value: { value: 15, format: 'percentage', precision: 0 },
          description: 'Frequency of MYC amplification as cooperating alteration',
          dataset: 'MM Genomics',
          sampleSize: 1000,
          source: 'MM genomic studies',
          category: 'validation'
        }
      ]
    },
    {
      title: 'Clinical Trial Shortlist Compression',
      description: 'Efficiency gains in clinical trial matching for MM patients',
      metrics: [
        {
          title: 'Trial Shortlist Reduction',
          value: { value: 85, format: 'percentage', precision: 0 },
          description: 'Reduction from 50+ to ~5-12 relevant trials',
          dataset: 'MM Clinical Trials',
          sampleSize: 100,
          source: 'MM trial matching analysis',
          category: 'business'
        },
        {
          title: 'Time to Shortlist',
          value: { value: 5, format: 'integer' },
          description: 'Minutes to generate trial shortlist',
          dataset: 'MM Clinical Trials',
          sampleSize: 100,
          source: 'MM trial matching analysis',
          category: 'business'
        }
      ]
    }
  ]
};

// Semantic alignment: Map stakeholder language → our in‑silico outputs (MM, RUO)
// Labels: live = available now; partial = scaffolded via proxies; roadmap = planned
export const mmSemanticAlignment = {
  safety: [
    { theySay: 'Off‑target editing', weSay: 'Genome‑wide off‑target profiling', status: 'partial', sourceFrom: 'CRISPR Readiness (demo), design/doctrine', notes: 'Preview only; full off‑target modeling is roadmap.' },
    { theySay: 'Cytotoxicity', weSay: 'Computational cytotoxicity prediction', status: 'live', sourceFrom: 'Toxicity Risk chip (/api/toxicity/hint)', notes: 'Plain caution + confidence + sources; RUO.' },
    { theySay: 'Immunogenicity', weSay: 'AI‑predicted immunogenicity risk', status: 'roadmap', sourceFrom: 'evidence/insights (future service)', notes: 'Planned endpoint surfaced as chip with provenance.' },
    { theySay: 'Pharmacokinetics (PK)', weSay: 'In‑silico PK modeling', status: 'roadmap', sourceFrom: 'PK model integration', notes: 'Expose as profile toggle; provenance required.' },
    { theySay: 'Pharmacodynamics (PD)', weSay: 'Computational PD profiling', status: 'partial', sourceFrom: 'Pathway impact + efficacy rationale', notes: 'PD proxy from pathway modulation and S/P signals.' },
    { theySay: 'Toxicity profiling', weSay: 'Multi‑organ toxicity prediction', status: 'roadmap', sourceFrom: 'toxicity bundle (future)', notes: 'Start with organ panels; add sources + confidence.' },
    { theySay: 'Dose‑limiting toxicity (DLT)', weSay: 'Preclinical DLT prediction', status: 'roadmap', sourceFrom: 'PK/PD + cohort overlays', notes: 'Surface conservative DLT risk with provenance.' }
  ],
  efficacy: [
    { theySay: 'Target engagement', weSay: 'Computationally validated target binding', status: 'partial', sourceFrom: 'Functionality + PathwayAligned + sequence disruption', notes: 'Binding affinity per‑se is roadmap; current proxy via S/P fit.' },
    { theySay: 'Biomarker response', weSay: 'AI‑predicted biomarker modulation', status: 'partial', sourceFrom: 'Regulatory/Chromatin chips + pathway impact', notes: 'Signals suggest modulation; pair with evidence when available.' },
    { theySay: 'Disease progression', weSay: 'Computational disease modeling', status: 'roadmap', sourceFrom: 'longitudinal simulator', notes: 'Future state; not claimed today.' },
    { theySay: 'Response rate', weSay: 'Predictive response stratification', status: 'partial', sourceFrom: 'Efficacy ranking + confidence + cohort_signals', notes: 'Stratification proxy via ranked therapies and small, auditable lifts.' },
    { theySay: 'Progression‑free survival (PFS)', weSay: 'In‑silico PFS modeling', status: 'roadmap', sourceFrom: 'simulator + cohorts', notes: 'Planned; gated by dataset availability.' },
    { theySay: 'Overall survival (OS)', weSay: 'Computational OS prediction', status: 'roadmap', sourceFrom: 'simulator + cohorts', notes: 'Planned; explicit RUO posture.' },
    { theySay: 'Objective response rate (ORR)', weSay: 'AI‑optimized ORR', status: 'roadmap', sourceFrom: 'predictor + cohorts', notes: 'Expose as endpoint with provenance when ready.' }
  ],
  moa: [
    { theySay: 'Target validation', weSay: 'Computational target validation', status: 'live', sourceFrom: 'Essentiality + truncation/frameshift gate', notes: 'Gene essentiality + Triumvirate compliance.' },
    { theySay: 'Pathway modulation', weSay: 'AI‑guided pathway engineering', status: 'live', sourceFrom: 'Pathway impact (P) + therapy PathwayAligned badges', notes: 'Explained in rationale.' },
    { theySay: 'Gene expression', weSay: 'Predictive gene expression modeling', status: 'roadmap', sourceFrom: 'expression model integration', notes: 'Display as chip with citations when available.' },
    { theySay: 'Protein‑protein interactions', weSay: 'Computational PPI prediction', status: 'roadmap', sourceFrom: 'structure/ML integration', notes: 'Surface as supportive rationale.' },
    { theySay: 'Binding affinity', weSay: 'AI‑optimized binding affinity', status: 'roadmap', sourceFrom: 'docking/ML integration', notes: 'Map into target engagement panel.' },
    { theySay: 'Selectivity', weSay: 'Computational selectivity profiling', status: 'roadmap', sourceFrom: 'off‑target/affinity differentials', notes: 'Pair with off‑target profiling outputs.' }
  ],
  trials: [
    { theySay: 'Primary endpoint', weSay: 'AI‑optimized primary endpoint design', status: 'roadmap', sourceFrom: 'trial design agent', notes: 'Template alongside shortlist.' },
    { theySay: 'Secondary endpoint', weSay: 'Computational secondary endpoint validation', status: 'roadmap', sourceFrom: 'trial design agent', notes: 'Propose with provenance.' },
    { theySay: 'Surrogate endpoint', weSay: 'AI‑validated surrogate biomarkers', status: 'partial', sourceFrom: 'biomarker chips + evidence', notes: 'Pair chips to candidate surrogates (RUO).' },
    { theySay: 'Biomarker validation', weSay: 'Computational biomarker discovery', status: 'partial', sourceFrom: 'insights + evidence', notes: 'Map to rationale with citations when active.' },
    { theySay: 'Patient stratification', weSay: 'AI‑driven patient stratification', status: 'partial', sourceFrom: 'efficacy ranking + confidence + cohort_signals', notes: 'Explain in trials one‑pager.' }
  ],
  regulatory: [
    { theySay: 'IND (Investigational New Drug)', weSay: 'AI‑enhanced IND package', status: 'partial', sourceFrom: 'Dossier export + provenance', notes: 'Clear RUO labels; live provenance.' },
    { theySay: 'Clinical trial design', weSay: 'Computational trial optimization', status: 'partial', sourceFrom: 'shortlist + rationale + future design agent', notes: 'Counts and categories are live in provenance.' },
    { theySay: 'Regulatory submission', weSay: 'AI‑validated regulatory package', status: 'roadmap', sourceFrom: 'evidence + trials design agent', notes: 'Keep conservative copy.' },
    { theySay: 'CMC', weSay: 'AI‑optimized CMC strategy', status: 'roadmap', sourceFrom: 'CMC planner', notes: 'Out of scope for MM demo; planned for platform.' }
  ],
  killerPhrases: [
    'From phenotypic screening to genotypic precision',
    'From empirical dosing to computational PK/PD modeling',
    'From biomarker discovery to AI‑validated endpoints',
    'From target identification to computational target validation',
    'From trial‑and‑error to predictive modeling',
    'From off‑target risk to genome‑wide off‑target profiling',
    'From patient heterogeneity to precision stratification'
  ],
  fearToValue: [
    { fear: 'Unknown mechanism', counter: 'Computational mechanism elucidation (pathway impact + rationale)' },
    { fear: 'Off‑target effects', counter: 'Genome‑wide off‑target profiling (demo today, roadmap full)' },
    { fear: 'Patient variability', counter: 'AI‑driven patient stratification (ranked therapies + cohort overlays)' },
    { fear: 'Dose optimization', counter: 'Computational dose modeling (PK/PD roadmap)' },
    { fear: 'Biomarker discovery', counter: 'AI‑accelerated biomarker identification (chips + evidence)' }
  ]
};

// How to use this in slides and product copy (MM)
export const mmNarrativeGuide = {
  executive: [
    'We turn noisy variants into a therapy game plan you can read and defend.',
    'Scores have receipts: every lift is small, toggleable, and provenanced (run IDs).',
    'Safety is built‑in: germline toxicity hint, with room to add PK/PD and organ panels.'
  ],
  connectOutputsToCare: [
    { label: 'Efficacy ranking', why: 'prioritize best bets; focus budget', from: 'drugs[*].{efficacy_score, confidence, tier, badges, rationale}' },
    { label: 'Pathway story', why: 'explain mechanism to boards', from: 'rationale[type=pathway], badges=PathwayAligned' },
    { label: 'Toxicity hint', why: 'plan conservatively', from: 'toxicity chip JSON' },
    { label: 'Cohorts/Trials', why: 'context and acceleration', from: 'cohort_signals, provenance.trials' }
  ],
  jsonAnchors: [
    'provenance.run_id (cite on slides)',
    'provenance.sequence_scoring.mode (Fusion/Baseline)',
    'cohort_signals (coverage_by_gene, response_rates)',
    'trials.shortlist_compression + categories'
  ]
};

// MM-specific capabilities
export const multipleMyelomaCapabilities = {
  variantInsight: {
    title: 'Variant Insight (VUS)',
    description: 'Four chips (Function, Regulatory, Essentiality, Chromatin) in plain language. Turn unknowns into readable signals with helper text and thresholds.',
    features: [
      'Function chip: Protein impact prediction',
      'Regulatory chip: Non-coding variant effects',
      'Essentiality chip: Gene essentiality scores',
      'Chromatin chip: Chromatin accessibility impact'
    ],
    whyItMatters: [
      'Reduces intake ambiguity; sets up therapy/pathway reasoning.'
    ],
    whatWeDelivered: [
      'Live chips with helpers and provenance; export-ready.'
    ]
  },
  therapyFit: {
    title: 'Therapy Fit (Chemo)',
    description: 'Ranked classes with confidence, short "why," and citations. Fusion labeled when AM coverage exists; Baseline deterministic otherwise.',
    features: [
      'PI (Proteasome Inhibitors) ranking',
      'IMiD (Immunomodulatory drugs) ranking',
      'Anti-CD38 monoclonal antibodies',
      'MAPK-aligned agents when RAS/RAF signal'
    ],
    whyItMatters: [
      'Explainable starting point for treatment planning (RUO).',
      'Shortens meetings by making rationale and sources explicit.'
    ],
    whatWeDelivered: [
      'Score, confidence, tier, badges, rationale; export + provenance.'
    ]
  },
  pathwayView: {
    title: 'Pathway View',
    description: 'Top 3 MM pathways with one-line "why" and contribution bars; links to therapy alignment.',
    features: [
      'MAPK pathway (BRAF/NRAS/KRAS)',
      'TP53/DDR pathway',
      'Proteostasis/CRBN pathway'
    ],
    whyItMatters: [
      'A fast biology story that justifies therapy choices.'
    ],
    whatWeDelivered: [
      'Stable top-3 pathways with bars, one-liners, and provenance.'
    ]
  },
  toxicityRisk: {
    title: 'Toxicity Risk (Germline)',
    description: 'Simple caution chip to plan conservatively. Confidence and sources included (RUO).',
    features: [
      'Germline variant screening',
      'Drug metabolism variants',
      'Toxicity risk scoring'
    ],
    whyItMatters: [
      'Flags potential sensitivity early; improves patient communication.'
    ],
    whatWeDelivered: [
      'Caution chip with helper, confidence, sources, provenance.'
    ]
  },
  crisprReadiness: {
    title: 'CRISPR Readiness (Demo)',
    description: 'Feasibility, access, off-target preview, delivery notes (demo). 1M-token context enables richer prompts.',
    features: [
      'On-target feasibility scoring',
      'Off-target prediction',
      'Delivery optimization',
      'Guide RNA design'
    ],
    whyItMatters: [
      'Faster, safer starts for design exploration (research-mode).'
    ],
    whatWeDelivered: [
      'Safety-gated candidates, access chip, demo off-target/delivery notes with provenance.'
    ]
  },
  clinicalTrials: {
    title: 'Clinical Trials Co-Pilot',
    description: 'Fast shortlist with Likely/Potential/Unlikely and a shareable one-pager. Synonym/biomarker-aware search and structured eligibility.',
    features: [
      'Smart trial matching',
      'Eligibility assessment',
      'One-pager export',
      'Real-time updates'
    ],
    whyItMatters: [
      'Reduces 50+ trials to ~5–12 in minutes; improves patient/board alignment.'
    ],
    whatWeDelivered: [
      'Shortlist with labels and “why”; export with run ID and sources.'
    ]
  }
};
