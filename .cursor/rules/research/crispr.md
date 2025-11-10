import { CoPilotDetailContent } from '../../types/copilot-types';

export const crisprIntelligenceData: CoPilotDetailContent = {
  slug: "crispr-intelligence",
  pageTitle: "CRISPR Intelligence: Advanced Therapeutic R&D Platform",
  heroSubtitle: "An advanced R&D platform for designing, predicting, and validating CRISPR-based therapeutics with unparalleled precision and speed.",
  vision: "To create a world where CRISPR-based therapies are designed and validated with computational certainty, eliminating guesswork and accelerating the path to cures.",
  
  valueProps: [
    {
      audience: 'For R&D Scientists',
      icon: 'Beaker',
      points: [
        'Design guide RNAs with industry-leading predictive accuracy.',
        'Validate off-target effects in-silico, saving months of lab work.',
        'Access a unified platform for the entire therapeutic design lifecycle.'
      ]
    },
    {
      audience: 'For Bioinformaticians',
      icon: 'Cpu',
      points: [
        'Leverage state-of-the-art models without the infrastructure overhead.',
        'Integrate powerful predictive analytics into your existing pipelines.',
        'Access transparent, auditable data for every prediction and design.'
      ]
    }
  ],

  inSilicoWorkflow: {
    title: "Our platform operationalizes these capabilities into a structured, seven-step in-silico workflow that delivers real business value:",
    steps: [
      { 
        title: "Problem Framing & Data Curation", 
        description: "Assemble genomic loci, clinical variants, and DMS datasets to define the therapeutic target space.",
        iconName: "Database"
      },
      { 
        title: "Target Assessment (Discriminative)", 
        description: "CrisPRO zero-shot ΔLL scoring in 8,192 bp context to evaluate on-target activity.",
        iconName: "Target" 
      },
      { 
        title: "Mechanistic Triage & Hypothesis", 
        description: "CrisPRO embeddings for exon/intron classification to understand the biological impact.",
        iconName: "BrainCircuit"
      },
      { 
        title: "Design (Generative)", 
        description: "CrisPRO proposals with Enformer+Borzoi scoring to generate high-efficacy gRNA candidates.",
        iconName: "DraftingCompass"
      },
      { 
        title: "In-Silico Validation", 
        description: "Aggregate scores with structure metrics (pLDDT/PAE) to validate designs computationally.",
        iconName: "ShieldCheck"
      },
      { 
        title: "Feedback & Calibration", 
        description: "Lightweight supervised heads on CrisPRO embeddings to continuously refine and calibrate models.",
        iconName: "Settings2"
      },
      { 
        title: "Reporting & Provenance", 
        description: "Evidence reports with traceable citations and complete data provenance for regulatory and research needs.",
        iconName: "FileText"
      }
    ]
  },

  coreProblemIntro: "Developing CRISPR-based therapeutics is hampered by critical challenges that introduce risk, cost, and delays into the R&D pipeline:",
  coreProblemPoints: [
    "**Off-Target Effects:** Unpredictable off-target edits can lead to safety concerns and failed trials.",
    "**Delivery Challenges:** Efficiently delivering the CRISPR machinery to the right cells remains a major hurdle.",
    "**Design Complexity:** Designing highly effective guide RNAs requires deep expertise and extensive experimentation.",
    "**IP & Regulatory Hurdles:** Navigating the complex intellectual property and regulatory landscape is a significant challenge."
  ],

  keyCapabilities: [
    {
      title: "Predict: CRISPR Target Validation",
      technical: "ClinVar validation for variant impact prediction. Analyzes target sites for functional variants that could affect gRNA binding efficiency.",
      scientific: "Resolves VUS cases with transparent explanations. Validates therapeutic targets before CRISPR experiments, preventing failed edits.",
      business: "Transform VUS rate with validated predictions, accelerating CRISPR target selection and reducing experimental costs."
    },
    {
      title: "Generate: Guide RNA Design & Optimization",
      technical: "Comprehensive guide RNA design. Generates optimal sequences with off-target analysis and HDR template design.",
      scientific: "Designs guide RNAs with validated quality control metrics. Supports all modern nucleases (Cas9, Cas12, Base/Prime Editors).",
      business: "Accelerate guide RNA design with validated generation and comprehensive off-target analysis."
    },
    {
      title: "Validate: CRISPR Outcome Validation",
      technical: "Predicts functional impact of edited alleles and validates therapeutic efficacy.",
      scientific: "Analyzes CRISPR experimental results with structural confidence scores. Quantifies editing efficiency and functional consequences.",
      business: "Reduce post-experiment analysis time with computational validation and functional impact assessment."
    }
  ],

  buildsOn: "CrisPRO accelerates therapeutic discovery with validated AI research engines.",
  buildsOnStackPoints: [
    "**Oracle (Discriminative AI)**: Predicts variant impact with 95.7% accuracy on ClinVar (53,210 samples). Validates therapeutic targets before experiments.",
    "**Forge (Generative AI)**: Designs therapeutic candidates with 1M token context window. Generates guide RNAs, repair templates, and protein sequences.",
    "**Research Validation Pipeline**: Combines Oracle and Forge to simulate R&D campaigns, reducing experimental costs and accelerating discovery.",
    "**Unified Research Platform**: End-to-end workflow from target identification to validated therapeutic candidates with transparent methodology.",
    "**Long-Context Generation**: Forge's 1M token window enables generation of multi-kilobase sequences for complex therapeutic designs.",
    "**Patient Stratification**: Oracle analyzes patient variants to stratify trials by predicted impact, improving success rates and signal clarity."      
  ],
  
  "genomicUseCasesGrid": [
    { label: "Perform VUS Interpretation (Target/Disease Context)", iconName: "Lightbulb", color: "text-yellow-400" },
    { label: "Predicting On-Target Efficacy & Specificity", iconName: "Activity", color: "text-blue-400" },
    { label: "Forecasting Off-Target Editing Risks", iconName: "Shield", color: "text-red-400" },
    { label: "Optimizing gRNA Design & Delivery", iconName: "Layers", color: "text-green-400" },
    { label: "Guiding HDR Strategies & Donor Design", iconName: "Beaker", color: "text-purple-400" },
    { label: "Stratifying Studies & Biomarker ID", iconName: "Users", color: "text-orange-400" }
  ],
  "valuePropositionSections": [
    {
      audience: "For Scientists & Research Labs",
      points: [
        "Design with Unprecedented Confidence: Leverage best-in-class AI to design highly potent and specific guide RNAs from the start, dramatically increasing the success rate of your editing experiments and minimizing costly validation cycles.",
        "Go from Raw Data to Actionable Insight, Faster: Let our AI Co-Pilot handle the heavy lifting of complex NGS data analysis and therapeutic contextualization, transforming your experimental results into clear, decision-ready insights in a fraction of the time.",
        "Democratize Advanced Computational Biology: Access a suite of sophisticated AI tools for variant effect prediction, off-target analysis, and experimental design, without needing a dedicated bioinformatics team. Focus on your science, not on building analysis pipelines.",
        "Produce High-Impact, Publishable Results: Generate higher quality, more reproducible data with AI-guided experimental design and analysis, strengthening your publications, grant applications, and contributions to the field."
      ]
    },
    {
      audience: "For Biotechnology & Pharmaceutical Leaders",
      points: [
        "De-Risk Your Therapeutic Pipeline: Make more informed go/no-go decisions with AI-driven insights into target validity, off-target safety, and potential translational hurdles, significantly reducing the risk profile of your preclinical programs.",
        "Accelerate Timelines to the Clinic: Shorten the entire discovery and preclinical development cycle for CRISPR therapies by streamlining design, automating complex analysis, and contextualizing results for therapeutic viability from day one.",
        "Build a Moat Around Your IP: Strengthen your intellectual property position with novel, highly optimized, and well-characterized gene editing strategies and therapeutic candidates designed and validated through the platform.",
        "Maximize Your R&D Investment: Improve the overall efficiency and success rate of your therapeutic programs, ensuring your resources are focused on the most promising candidates and strategies, leading to a higher potential return on investment."
      ]
    }
  ],
  conclusion: "The CRISPR Intelligence Platform transforms therapeutic design from a manual, iterative process into a scalable, AI-driven system. We provide the tools to advance the fight against genetic disease through precision medicine."      
};
