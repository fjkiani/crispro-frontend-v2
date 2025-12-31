import { CoPilotDetailContent } from '@/types/copilot-types';

export const toxicityData: CoPilotDetailContent = {
  slug: "toxicity-risk",
  pageTitle: "Toxicity Risk Assessment",
  heroSubtitle: "Identifies pharmacogene variants and pathway overlap risks. Recommends mitigating foods with personalized explanations.",
  vision: "Prevent drug toxicity by identifying germline risks and recommending protective foods.",

  valueProps: [
    {
      audience: 'For Clinicians',
      icon: 'Shield',
      points: [
        'Flags pharmacogene variants (DPYD, TPMT, UGT1A1, CYP2D6) that affect drug metabolism.',
        'Identifies pathway overlap between patient germline variants and drug toxic pathways.',
        'Recommends pathway-specific foods (NAC, Omega-3, CoQ10) with timing and dosage.'
      ]
    },
    {
      audience: 'For Researchers',
      icon: 'FileText',
      points: [
        'Transparent pathway-based risk assessment (DNA repair, inflammation, cardiometabolic).',
        'Mechanism-based food recommendations with LLM-enhanced explanations.',
        'Conservative approach: Better to flag potential risks than miss them.'
      ]
    }
  ],

  buildsOn: "Core Capabilities",
  buildsOnStackPoints: [
    "**Pharmacogene detection:** 20+ drug-metabolizing enzymes (DPYD, TPMT, UGT1A1, CYP2D6).",
    "**Pathway overlap:** Maps drug MoA to toxic pathways and computes overlap with patient variants.",
    "**Mitigating foods:** Pathway-specific recommendations (NAC for DNA repair, Omega-3 for inflammation, CoQ10 for cardiometabolic)."
  ],

  kpis: [
    { label: 'Pharmacogenes', value: '20+ genes' },
    { label: 'Toxicity Pathways', value: '3 pathways' },
    { label: 'Drug MoAs', value: '15+ mechanisms' },
    { label: 'Mitigating Foods', value: '9 compounds' }
  ],

  observedOutcomes: [
    {
      title: "Pharmacogene Detection",
      keyMetric: "20+ genes",
      description: "Identifies variants in DPYD, TPMT, UGT1A1, CYP2D6, and other drug-metabolizing enzymes.",
      icon: "Dna",
      color: "blue"
    },
    {
      title: "Pathway Overlap Analysis",
      keyMetric: "3 pathways",
      description: "Computes overlap between patient germline variants and drug MoA toxic pathways.",
      icon: "Activity",
      color: "green"
    },
    {
      title: "Mitigating Foods",
      keyMetric: "9 compounds",
      description: "Recommends pathway-specific foods (NAC, Vitamin D, Omega-3, CoQ10) with timing and dosage.",
      icon: "Apple",
      color: "purple"
    },
    {
      title: "LLM Explanations",
      keyMetric: "✅ Built",
      description: "Personalized explanations for why foods mitigate toxicity, with patient-friendly summaries.",
      icon: "MessageSquare",
      color: "orange"
    }
  ],

  genomicInsightsOverview: "Identifies pharmacogene variants and pathway overlap risks. Recommends mitigating foods with personalized LLM explanations.",
  coreProblemIntro: "Germline variants can cause severe drug toxicity. We identify risks and recommend protective foods.",
  coreProblemPoints: [
    "Pharmacogene variants affect drug breakdown (DPYD, TPMT, UGT1A1).",
    "Drug MoA overlaps with patient's germline pathway vulnerabilities.",
    "Need mechanism-based food recommendations."
  ],

  genomicUseCasesGrid: [
    { label: "Pharmacogene detection", iconName: "Dna", color: "text-blue-400" },
    { label: "Pathway overlap", iconName: "Activity", color: "text-green-400" },
    { label: "Mitigating foods", iconName: "Apple", color: "text-purple-400" },
    { label: "LLM explanations", iconName: "MessageSquare", color: "text-orange-400" }
  ],

  keyCapabilities: [
    {
      title: "Pharmacogene Detection",
      technical: {
        title: "Implementation",
        keyMetric: "20+ genes",
        description: "Detects variants in drug-metabolizing enzymes. High-impact genes (DPYD, TPMT, UGT1A1) get risk weight 0.4.",
        icon: "Dna",
        color: "blue",
        components: [
          {
            title: "Gene List",
            subtitle: "DPYD, TPMT, UGT1A1, CYP2D6, CYP2C19, ABCB1, SLCO1B1, and 13+ more",
            iconName: "Dna",
            color: "blue"
          },
          {
            title: "Risk Weighting",
            subtitle: "High-impact (0.4), CYP enzymes (0.3), others (0.2)",
            iconName: "Gauge",
            color: "teal"
          }
        ]
      },
      scientific: {
        title: "Foundation",
        keyMetric: "Mechanism-Based",
        description: "DPYD variants cause 5-FU toxicity. TPMT variants cause thiopurine toxicity. UGT1A1*28 causes irinotecan toxicity.",
        icon: "BookOpen",
        color: "teal",
        components: [
          {
            title: "High-Impact",
            subtitle: "DPYD (5-FU), TPMT (thiopurines), UGT1A1 (irinotecan), G6PD, NUDT15",
            iconName: "AlertTriangle",
            color: "blue"
          },
          {
            title: "CYP Enzymes",
            subtitle: "CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP3A5",
            iconName: "Activity",
            color: "teal"
          }
        ]
      },
      business: {
        title: "Value",
        keyMetric: "Risk Prevention",
        description: "Flags patients at risk before prescribing. Enables dose adjustments or alternative drugs.",
        icon: "Shield",
        color: "indigo",
        components: [
          {
            title: "Proactive Screening",
            subtitle: "Flags variants before drug administration",
            iconName: "Shield",
            color: "blue"
          },
          {
            title: "Dose Guidance",
            subtitle: "High-risk → avoid drug or reduce dose 50-90%",
            iconName: "Pill",
            color: "teal"
          }
        ]
      },
      genomicUseCasesParagraph: "**Built:** 20+ pharmacogenes detected with risk weighting. Integrated with safety service and food validation."
    },
    {
      title: "Pathway Overlap Analysis",
      technical: {
        title: "Implementation",
        keyMetric: "3 pathways",
        description: "Maps drug MoA to toxic pathways and computes overlap with patient germline variants.",
        icon: "Activity",
        color: "blue",
        components: [
          {
            title: "MoA Mapping",
            subtitle: "15+ drug MoAs → toxic pathways (platinum → DNA repair: 0.9, anthracycline → cardiometabolic: 0.9)",
            iconName: "Map",
            color: "blue"
          },
          {
            title: "Pathway Genes",
            subtitle: "DNA repair (BRCA1/2, MBD4, TP53), inflammation (TNF, IL6, NFKB), cardiometabolic (APOB, MTOR)",
            iconName: "Dna",
            color: "teal"
          }
        ]
      },
      scientific: {
        title: "Foundation",
        keyMetric: "Mechanism-Based",
        description: "Platinum agents damage DNA. Anthracyclines cause cardiotoxicity. Checkpoint inhibitors cause iRAEs.",
        icon: "Microscope",
        color: "teal",
        components: [
          {
            title: "DNA Repair",
            subtitle: "Platinum, PARP inhibitors, alkylating agents → DNA damage",
            iconName: "Shield",
            color: "blue"
          },
          {
            title: "Inflammation",
            subtitle: "Checkpoint inhibitors, CAR-T, IMiDs → Immune activation",
            iconName: "Flame",
            color: "teal"
          },
          {
            title: "Cardiometabolic",
            subtitle: "Anthracyclines, BRAF/MEK inhibitors → Cardiac toxicity",
            iconName: "Heart",
            color: "indigo"
          }
        ]
      },
      business: {
        title: "Value",
        keyMetric: "Mechanism-Based Risk",
        description: "Identifies which toxicity pathway is at risk. Enables targeted food recommendations.",
        icon: "Target",
        color: "indigo",
        components: [
          {
            title: "Pathway-Specific",
            subtitle: "Identifies DNA repair, inflammation, or cardiometabolic risk",
            iconName: "Target",
            color: "blue"
          },
          {
            title: "Food Mapping",
            subtitle: "Pathway → specific foods (NAC, Omega-3, CoQ10)",
            iconName: "Apple",
            color: "teal"
          }
        ]
      },
      genomicUseCasesParagraph: "**Built:** 3 toxic pathways, 15+ drug MoAs mapped. Overlap computation with patient variants."
    },
    {
      title: "Mitigating Foods",
      technical: {
        title: "Implementation",
        keyMetric: "9 compounds",
        description: "Maps pathway overlap to mitigating foods with timing, dosage, and LLM explanations.",
        icon: "Apple",
        color: "blue",
        components: [
          {
            title: "Pathway Mapping",
            subtitle: "DNA repair → NAC, Vitamin D, Folate; Inflammation → Omega-3, Curcumin, EGCG; Cardiometabolic → CoQ10, L-Carnitine, Magnesium",
            iconName: "Apple",
            color: "blue"
          },
          {
            title: "Timing & Dosage",
            subtitle: "Post-infusion (NAC), continuous (Vitamin D), between meals (Curcumin), with fatty meal (CoQ10)",
            iconName: "Clock",
            color: "teal"
          }
        ]
      },
      scientific: {
        title: "Foundation",
        keyMetric: "Evidence-Based",
        description: "NAC (glutathione for DNA repair), Omega-3 (resolvin for inflammation), CoQ10 (mitochondrial support for cardiotoxicity).",
        icon: "BookOpen",
        color: "teal",
        components: [
          {
            title: "DNA Repair",
            subtitle: "NAC (glutathione), Vitamin D (VDR-mediated), Folate (DNA synthesis)",
            iconName: "Shield",
            color: "blue"
          },
          {
            title: "Anti-Inflammatory",
            subtitle: "Omega-3 (NF-κB), Curcumin (COX-2), EGCG (STAT3)",
            iconName: "Flame",
            color: "teal"
          },
          {
            title: "Cardioprotective",
            subtitle: "CoQ10 (mitochondrial), L-Carnitine (fatty acids), Magnesium (QT)",
            iconName: "Heart",
            color: "indigo"
          }
        ]
      },
      business: {
        title: "Value",
        keyMetric: "THE MOAT",
        description: "Connects toxicity assessment to actionable food recommendations. LLM explanations improve adherence.",
        icon: "Target",
        color: "indigo",
        components: [
          {
            title: "Actionable",
            subtitle: "Specific foods, doses, timing for each pathway",
            iconName: "CheckCircle",
            color: "blue"
          },
          {
            title: "Personalized",
            subtitle: "LLM-generated explanations for patient's specific context",
            iconName: "MessageSquare",
            color: "teal"
          }
        ]
      },
      genomicUseCasesParagraph: "**Built:** 9 compounds with pathway mapping, timing, dosage. LLM enhancement for personalized explanations."
    }
  ],

  valuePropositionSections: [
    {
      audience: "For Clinicians",
      points: [
        "Identifies pharmacogene variants that affect drug metabolism.",
        "Computes pathway overlap between patient variants and drug toxic pathways.",
        "Recommends mitigating foods with personalized explanations."
      ]
    }
  ],

  conclusion: "Toxicity risk assessment identifies pharmacogene variants and pathway overlap risks. Recommends pathway-specific mitigating foods with personalized LLM explanations.",

  inSilicoOverview: {
    coreConcepts: [
      {
        icon: "Dna",
        title: "Pharmacogene Detection",
        description: "20+ drug-metabolizing enzymes (DPYD, TPMT, UGT1A1, CYP2D6) with risk weighting.",
        color: "blue"
      },
      {
        icon: "Activity",
        title: "Pathway Overlap",
        description: "Computes overlap between patient germline variants and drug MoA toxic pathways.",
        color: "teal"
      },
      {
        icon: "Apple",
        title: "Mitigating Foods",
        description: "Pathway-specific foods (NAC, Omega-3, CoQ10) with timing, dosage, and LLM explanations.",
        color: "purple"
      }
    ],
    valuePropositions: [
      {
        icon: "Shield",
        title: "Risk Prevention",
        description: "Identifies patients at risk before prescribing",
        metric: "20+ pharmacogenes",
        color: "blue"
      },
      {
        icon: "Activity",
        title: "Mechanism-Based",
        description: "Pathway overlap identifies which toxicity pathway is at risk",
        metric: "3 pathways",
        color: "teal"
      },
      {
        icon: "Apple",
        title: "Actionable Foods",
        description: "Pathway-specific foods with timing, dosage, and explanations",
        metric: "9 compounds",
        color: "indigo"
      }
    ],
    deliverables: [
      {
        icon: "Dna",
        title: "Pharmacogene Detection",
        description: "Variants in drug-metabolizing enzymes with risk weighting"
      },
      {
        icon: "Activity",
        title: "Pathway Overlap",
        description: "Overlap between patient variants and drug toxic pathways"
      },
      {
        icon: "Apple",
        title: "Mitigating Foods",
        description: "Pathway-specific foods with timing, dosage, and LLM explanations"
      }
    ]
  }
};
