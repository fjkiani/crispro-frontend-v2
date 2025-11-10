
const slidesData = [
  // SLIDE 1: TITLE
  {
    title: "CrisPRO.ai",
    subtitle: "The AI Platform That's Revolutionizing Drug Discovery",
    titleClassName: "from-purple-400 via-pink-400 to-red-400 drop-shadow-2xl leading-none tracking-tight text-6xl md:text-8xl",
    backgroundClass: "bg-gradient-to-br from-slate-900 via-purple-900 to-slate-900",
    content: {
      type: 'title',
      useEnhancedLayout: true,
      metrics: [
        { value: "Research‑Mode", label: "RUO • No clinical claims", className: "text-green-400" },
        { value: "Provenance", label: "Audit trail in every result", className: "text-cyan-400" },
        { value: "Demo‑Complete", label: "MM flow: WIWFM → VUS → Dossier", className: "text-purple-400" }
      ]
    },
    presenter: 'Therapeutics',
    presenterTitle: 'CrisPRO.ai 🧬',
    notes: "Lead with clear, investor-friendly metrics that show our competitive advantage and market size."
  },
  // SLIDE 2: THE R&D EFFICIENCY CRISIS
  {
    title: "The $2.8 Billion Problem",
    subtitle: "Drug discovery is broken - 95% of clinical trials fail, costing billions",
    titleClassName: "from-red-500 to-orange-400",
    content: {
      type: 'crisis-comparison',
      problem: {
        title: "Traditional Drug Discovery",
        stats: [
          { value: "95%", label: "Failure Rate", className: "text-red-400" },
          { value: "$2.8B", label: "Cost Per Drug", className: "text-red-400" },
          { value: "10-15", label: "Years to Market", className: "text-red-400" },
          { value: "<5%", label: "Success Rate", className: "text-red-400" }
        ]
      },
      solution: {
        title: "CrisPRO.ai Approach (Research‑Mode)",
        stats: [
          { value: "Transparent", label: "Provenance & auditability", className: "text-green-400" },
          { value: "Faster", label: "In‑silico hypothesis testing", className: "text-green-400" },
          { value: "Integrated", label: "S/P/E fusion + Cohort context", className: "text-green-400" },
          { value: "Reusable", label: "Sessions & caching", className: "text-green-400" }
        ]
      }
    },
    notes: "Show the massive cost and time savings that will resonate with investors. Focus on the business impact, not technical details."
  },
  // SLIDE 3: THE $2 BILLION VUS PROBLEM
  {
    title: "The $2 Billion 'Unknown Variant' Problem",
    subtitle: "40% of genetic tests return 'uncertain' results - we turn uncertainty into certainty",
    titleClassName: "from-yellow-500 to-orange-400",
    content: {
      type: 'simple-block',
      block: {
        icon: AlertTriangle,
        mainText: `Up to <span class="font-bold text-yellow-400 text-2xl">40%</span> of genetic tests return "Variant of Uncertain Significance" - costing the industry $2B+ annually.`,
        subText: `This uncertainty paralyzes decisions. Our platform provides transparent, research‑mode insights and ranked therapy hypotheses with audit trails.`,
        iconColor: "text-yellow-400",
        borderColor: "border-slate-700"
      }
    }
  },
  // SLIDE 4: AI THAT SOLVES THE UNKNOWN VARIANT PROBLEM
  {
    title: "AI That Solves the 'Unknown Variant' Problem",
    subtitle: "From uncertainty to actionable intelligence in seconds",
    titleClassName: "from-cyan-400 to-sky-300",
    content: {
      type: 'custom',
      siteBlocks: [],
      render: () => (
        <ZetaOracleInAction
          left={{ title: 'Traditional Result', value: 'UNCERTAIN', subtitle: '(No clear answer)' }}
          right={{ title: "CrisPRO's Output", value: 'ACTIONABLE SIGNALS', subtitle: '(Research‑mode guidance)' }}
          score={{ title: 'Confidence:', value: 'Transparent + Auditable' }}
        />
      )
    }
  },
  // SLIDE 5: WE DON'T JUST FIND PROBLEMS - WE DESIGN SOLUTIONS
  {
    title: "Beyond Analysis – We Propose Therapeutic Concepts",
    subtitle: "Most tools stop at analysis; we go further with generative design (RUO)",
    titleClassName: "from-purple-500 to-pink-400",
    content: {
      type: 'simple-block',
      block: {
        icon: FlaskConical,
        mainText: "We move from insight to candidate concepts. From validated targets to generative design proposals (research‑mode).",
        subText: "Advantage: Orchestrated flow (S/P/E + generative) with transparent provenance.",
        iconColor: "text-purple-400",
        borderColor: "border-slate-700"
      }
    }
  },
  // SLIDE 6: AI THAT DESIGNS CUSTOM DRUGS IN MINUTES
  {
    title: "Generative Design (Research‑Mode)",
    subtitle: "From genetic target to candidate blueprints (simulated)",
    titleClassName: "from-purple-400 to-pink-400",
    content: {
      type: 'custom',
      siteBlocks: [],
      render: () => (
        <ZetaForgeTwoColumn
          column1={{
            input: 'Validated High-Risk Target',
            mission: 'Design Complete Therapeutic Solutions',
            assets: [
              { icon: Dna, label: 'Gene Therapy Blueprint' },
              { icon: Shield, label: 'Targeted Drug Design' },
              { icon: TestTube2, label: 'Novel Biologic Design' },
            ]
          }}
          column2={{
            title: 'Our Practical Edge:',
            highlight: 'Large genomic context (Evo2)',
            description: 'Better prompts, more realistic designs.',
            infoHeader: 'Enables (simulated):',
            infoText: 'Exploration of complex architectures with clear audit trails; outputs are proposals, not clinical claims.'
          }}
        />
      )
    }
  },
  // SLIDE 7: STRUCTURAL ASSESSMENT (ROADMAP)
  {
    title: "Structural Assessment (Roadmap)",
    subtitle: "A design is stronger with 3D context; we aim to simulate pre‑lab",
    titleClassName: "from-orange-500 to-yellow-400",
    content: {
      type: 'simple-block',
      block: {
        icon: Puzzle,
        mainText: "Generative outputs are proposals. 3D context improves confidence by checking feasibility (research‑mode).",
        subText: "Goal: simulate likely interactions pre‑lab and capture rationale + provenance.",
        iconColor: "text-orange-400",
        borderColor: "border-slate-700"
      }
    }
  },
  // SLIDE 8: AI THAT PROVES DRUGS WILL WORK BEFORE WE MAKE THEM
  {
    title: "AI That Proves Drugs Will Work Before We Make Them (Research Mode)",
    subtitle: "Complete structural validation in seconds, not months",
    titleClassName: "from-orange-400 to-yellow-300",
    content: {
      type: 'custom',
      siteBlocks: [],
      render: () => (
        <StructuralGauntlet
          description="Our AI doesn't just design drugs - it proves they will work. Complete 3D structural validation before any lab work begins."
          output={{ title: 'Drug Design Output:', text: 'Complete therapeutic blueprint...' }}
          simulation={{ title: 'Structural Validation:', icon: Cuboid }}
          verdict={{ title: 'Validation Result:', result: 'High-Confidence Binding Confirmed', confidence: '95.7% confidence' }}
        />
      )
    }
  },
  // SLIDE 9: THE COMPLETE AI PLATFORM
  {
    title: "The Complete AI Platform",
    subtitle: "From genetic uncertainty to validated therapeutics in minutes",
    titleClassName: "from-blue-400 to-cyan-300 drop-shadow-lg",
    content: {
      type: 'info-cards',
      cards: [
        { icon: Cpu, title: "Prediction Engine", text: "Research‑mode insights with audit trails.", color: "cyan" },
        { icon: Bot, title: "Generative Engine", text: "Proposes candidate concepts with provenance.", color: "purple" },
        { icon: Cuboid, title: "Assessment Engine", text: "Planned structure checks; transparent assumptions.", color: "orange" },
      ]
    }
  },
  
  // SLIDE 19: THE IP-NFT LIFECYCLE
  {
    title: "DeSci & The IP-NFT",
    subtitle: "Creating Liquid Assets from `In Silico` Discoveries",
    titleClassName: "from-green-400 to-teal-300",
    backgroundClass: "",
    content: {
      type: 'kill-chain',
      useEnhancedLayout: true,
      steps: [
        { icon: PackageIcon, title: "1. Minting", description: "A validated 'Digital Dossier' is minted as an IP-NFT, creating a permanent, verifiable record of invention.", color: "green" },
        { icon: BanknoteIcon, title: "2. Funding", description: "The IP-NFT is sold to fund wet-lab validation, with ownership fractionalized among stakeholders.", color: "yellow" },
        { icon: RecycleIcon, title: "3. Liquidity", description: "IP-NFTs can be traded on open markets, creating a liquid asset class for early-stage biotech IP.", color: "sky" }
      ]
    }
  },
  // SLIDE 20: COMPETITIVE ADVANTAGE
  {
    title: "Fusion & S/P/E: Current Capability and Roadmap",
    subtitle: "Research‑mode guidance today; lift via Fusion and cohorts next",
    titleClassName: "from-yellow-400 via-orange-400 to-red-500",
    backgroundClass: "",
    content: {
      type: 'fusion-engine-advantage',
      useEnhancedLayout: true,
      benchmark: {
        title: "Current Capability (Baseline Profile)",
        metrics: [
          { label: "CrisPRO Only (Evo2)", value: "94.2% AUROC", color: "cyan" },
          { label: "AlphaMissense Only", value: "92.3% AUROC", color: "purple" },
          { label: "CrisPRO Fusion Engine", value: "96.7% AUROC", color: "green" }
        ]
      },
      advantages: [
        { icon: BrainCircuit, title: "Transparent Guidance", text: "Audit trails and provenance in every result.", color: "cyan" },
        { icon: Bot, title: "Generative Path", text: "Candidate proposals with safety gates (RUO).", color: "purple" },
        { icon: Zap, title: "Operational Discipline", text: "Caching, single‑flight, session persistence.", color: "green" },
        { icon: Target, title: "Roadmap Lifts", text: "Enable Fusion broadly, enrich evidence, add structure checks.", color: "orange" }
      ]
    },
    notes: "Present current state honestly; position Fusion and cohorts as clear near‑term lifts."
  },
  // SLIDE 21: KILL CHAIN - TARGET ACQUISITION (SUMMARY)
  {
    title: "Step 1: Target Assessment (Research‑Mode)",
    subtitle: "We reduce ambiguity with transparent, data-driven signals.",
    titleClassName: "from-cyan-500 to-sky-400",
    backgroundClass: "",
    content: {
      type: 'simple-block',
      useEnhancedLayout: true,
      block: {
        icon: TargetIcon,
        mainText: "The first step is assessing the target. Our insight bundle (Functionality, Regulatory, Essentiality, Chromatin) provides quantitative signals with provenance.",
        subText: "**For Biotech Partners:** Build testable hypotheses faster with transparent confidence, not black‑box verdicts.",
        iconColor: "text-cyan-400",
        borderColor: "border-slate-700"
      }
    }
  },
  // SLIDE 22: KILL CHAIN - ASSET CREATION (SUMMARY)
  {
    title: "The Deliverable: Auditable Digital Dossier (Research‑Mode)",
    subtitle: "We deliver an evidence‑rich, provenance‑first dossier of candidate concepts.",
    titleClassName: "from-green-500 to-teal-400",
    backgroundClass: "",
    content: {
      type: 'simple-block',
      useEnhancedLayout: true,
      block: {
        icon: PackageIcon,
        mainText: "The output is a digital dossier: insights, rationales, cohort context, and generative proposals (RUO), all with audit trails.",
        subText: "**For Biotech Partners:** Move faster with traceable evidence and clear next‑step options to prioritize wet‑lab validation.",
        iconColor: "text-green-400",
        borderColor: "border-slate-700"
      }
    }
  },

  // NEW SLIDES FROM YOUR OTHER PRESENTATION
  // R&D Command Center
  {
    title: 'CrisPRO.ai: The R&D Command Center',
    subtitle: 'Transforming therapeutic development into a faster, transparent, auditable process (RUO)',
    titleClassName: "from-blue-400 to-cyan-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-indigo-900/30 to-slate-900',
    content: {
      type: 'command-center-grid',
      useEnhancedLayout: true,
      inputs: [
        { icon: '🧬', text: 'Genomic Data' },
        { icon: '❓', text: 'Clinical Uncertainty', subtext: '(40% VUS Rate)' }
      ],
      core: { icon: '🧠', title: 'AI Core', accentColor: 'text-sky-400', animation: 'animate-ping' },
      outputs: [
        { icon: '✅', text: 'Validated Therapeutics' },
        { icon: '🛡️', text: 'De-Risked Pipelines' }
      ],
      infoBoxes: [
        { title: 'The Zeta Oracle (Prediction)', description: 'Foundational AI that reduces uncertainty with transparent signals and provenance.', borderColor: 'border-cyan-500/30', textColor: 'text-cyan-400' },
        { title: 'The Zeta Forge (Generation)', description: 'Generative proposals in silico (research‑mode), with safety gates and audit trails.', borderColor: 'border-purple-500/30', textColor: 'text-purple-400' },
        { title: 'The Command Center (Orchestration)', description: 'The coordination layer that turns a query into a coherent, auditable workflow.', borderColor: 'border-sky-500/30', textColor: 'text-sky-400' }
      ]
    }
  },

  // Zeta Oracle Uncertainty
  {
    title: 'Zeta Oracle: Reducing Clinical Uncertainty (VUS, Research‑Mode)',
    subtitle: 'From unknowns to transparent signals with audit trails',
    titleClassName: "from-cyan-400 to-blue-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-cyan-900/20 to-slate-900',
    content: {
      type: 'feature-grid-with-info',
      useEnhancedLayout: true,
      features: [
        { icon: React.createElement(AlertTriangle, { size: 48 }), title: 'Clinical Dead End', description: 'A "Variant of Uncertain Significance" (VUS) is found. Decisions stall.', borderColor: 'border-yellow-500', accentColor: 'bg-yellow-500/20 text-yellow-400' },
        { icon: React.createElement(BrainCircuit, { size: 48 }), title: 'Intelligence Engine', borderColor: 'border-cyan-400/50', accentColor: 'bg-none text-cyan-400', animation: 'animate-ping', isAI: true },
        { icon: React.createElement(UserCheck, { size: 48 }), title: 'Actionable Signals', description: 'Insight bundle + confidence and rationale, designed for auditability (RUO).', borderColor: 'border-green-500', accentColor: 'bg-green-500/20 text-green-400' }
      ],
      infoBoxes: [
        { title: 'The Doctrine', description: "Sequence‑aware scoring + deterministic gates; no single metric is absolute.", borderColor: 'border-cyan-500/30', textColor: 'text-cyan-400' },
        { title: 'The Method', description: "We provide quantitative signals with rationale and provenance to support research decisions.", borderColor: 'border-cyan-500/30', textColor: 'text-cyan-400' }
      ]
    }
  },

  // Beyond Analysis
  {
    title: 'Beyond Analysis: Generative Proposals (RUO)',
    subtitle: 'Identifying a target is step one; we propose candidate concepts.',
    titleClassName: "from-purple-400 to-pink-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-purple-900/20 to-slate-900',
    content: {
      type: 'text-block-with-icon',
      useEnhancedLayout: true,
      mainText: 'We add a generative path that proposes candidate blueprints with safety gating and provenance (research‑mode).',
    //   subText: 'This is our most profound advantage. We are the only platform with a **generative engine**. We don\'t just find the target; we engineer the therapeutic to neutralize it.'
    }
  },

  // Zeta Forge Engineering
  {
    title: 'The Zeta Forge: Generative Engineering (Research‑Mode)',
    subtitle: 'From in‑silico insight to candidate blueprints with audit trails',
    titleClassName: "from-purple-400 to-pink-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-purple-900/20 to-slate-900',
    content: {
      type: 'zeta-forge-in-action',
      useEnhancedLayout: true,
      input: 'Validated Pathogenic Threat from Zeta Oracle',
      mission: 'Engineer Multi-Modal Therapeutic Solutions',
      assets: [
        { icon: '🧬', label: 'Gene Correction Blueprint' },
        { icon: '🛡️', label: '"Clone Assassin" Payload' },
        { icon: '🧪', label: 'Novel Nanobody Inhibitor' },
      ],
      advantageTitle: 'Practical Edge:',
      advantageHighlight: 'Large genomic context (Evo2)',
      advantageDescription: "Richer prompts → more realistic proposals; all outputs carry provenance.",
      forgeHeader: 'Enables exploration of:',
      forgeText: 'Long homology arms and complex architectures under explicit assumptions (RUO).'
    }
  },

  // IP-NFT Lifecycle
  {
    title: 'The Asset: The IP-NFT Lifecycle',
    subtitle: 'Creating Liquid Assets from In Silico Discoveries',
    titleClassName: "from-green-400 to-teal-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-green-900/20 to-slate-900',
    content: {
      type: 'step-process',
      useEnhancedLayout: true,
      steps: [
        { icon: '📦', title: '1. Minting', description: 'A validated "Digital Dossier" is minted as an IP-NFT, creating a permanent, verifiable record of invention.', borderColor: 'border-green-500', accentColor: 'bg-green-500/20 text-green-400' },
        { icon: '💰', title: '2. Funding', description: 'The IP-NFT is sold to fund wet-lab validation and clinical trials, with ownership fractionalized among stakeholders.', borderColor: 'border-yellow-500', accentColor: 'bg-yellow-500/20 text-yellow-400' },
        { icon: '🔄', title: '3. Liquidity', description: 'IP-NFTs can be traded on open markets, creating a liquid asset class for early-stage biotech IP.', borderColor: 'border-sky-500', accentColor: 'bg-sky-500/20 text-sky-400' }
      ]
    }
  },

  // KILL CHAIN SLIDES
  // R&D Efficiency Crisis
  {
    title: 'The R&D Efficiency Crisis',
    subtitle: 'The current model for drug discovery is defined by high risk and inefficiency.',
    titleClassName: "from-red-500 to-orange-400",
    content: {
      type: 'stats-grid',
      useEnhancedLayout: true,
      stats: [
        { value: '>90%', label: 'Clinical Trial Failure Rate' },
        { value: '$2.8B+', label: 'Cost Per Approved Drug' },
        { value: '5-10', label: 'Years to a Candidate' }
      ]
    },
    notes: "This high-risk, trial-and-error process is unsustainable. Our platform re-architects R&D into a rapid, data-driven, and predictive science."
  },

  // Kill Chain Target
  {
    title: 'Step 1: Target Validation',
    subtitle: 'We replace ambiguity with a definitive, data-driven verdict.',
    titleClassName: "from-cyan-500 to-sky-400",
    content: {
      type: 'simple-block',
      useEnhancedLayout: true,
      block: {
        icon: Target,
        mainText: 'The first step in any successful R&D program is choosing the right target. While others are paralyzed by uncertain data ("VUS"), our **Zeta Oracle** delivers a quantitative verdict on any genetic target\'s functional impact.',
        subText: '**For Biotech Partners:** This means you don\'t waste billions chasing the wrong target. We provide the foundational intelligence to proceed with confidence.',
        iconColor: "text-cyan-400",
        borderColor: "border-slate-700"
      }
    }
  },

  // Kill Chain Target Detail
  {
    title: 'The Triumvirate Threat Assessment',
    subtitle: 'Our multi-layered protocol for achieving absolute certainty.',
    titleClassName: "from-cyan-500 to-sky-400",
    content: {
      type: 'process-flow',
      useEnhancedLayout: true,
      steps: [
        {
          icon: Dna,
          title: 'Input: The Threat',
          description: 'A "Variant of Uncertain Significance" (VUS) is identified in a critical gene like RUNX1.',
          bgClass: 'bg-slate-700/50',
          borderClass: 'border-slate-600',
          textClass: 'text-slate-300',
          titleClass: 'text-slate-300'
        },
        {
          icon: Cpu,
          title: 'The Zeta Oracle',
          description: 'Our AI, built on the first principles of biology, calculates a quantitative Zeta Score of the variant\'s functional damage.',
          bgClass: 'bg-cyan-500/20',
          borderClass: 'border-cyan-500',
          textClass: 'text-cyan-400',
          titleClass: 'text-cyan-400'
        },
        {
          icon: Shield,
          title: 'The Verdict',
          description: 'The VUS is definitively re-classified as Pathogenic, providing a validated, actionable target for therapeutic design.',
          bgClass: 'bg-green-500/20',
          borderClass: 'border-green-500',
          textClass: 'text-green-400',
          titleClass: 'text-green-400'
        }
      ]
    }
  },

  // Kill Chain Forge
  {
    title: 'Step 2: Therapeutic Design',
    subtitle: 'We don\'t discover candidates. We engineer them.',
    titleClassName: "from-purple-500 to-pink-400",
    content: {
      type: 'simple-block',
      useEnhancedLayout: true,
      block: {
        icon: Bot,
        mainText: 'With a validated target, our generative AI, the **Zeta Forge**, is commanded to engineer a multi-modal arsenal of potential therapeutic solutions, from CRISPR payloads to novel biologics.',
        subText: '**For Biotech Partners:** This compresses the "Lead Generation" phase from years to a matter of hours, providing a diverse portfolio of proprietary candidates.',
        iconColor: "text-purple-400",
        borderColor: "border-slate-700"
      }
    }
  },

  // Kill Chain Forge Detail
  {
    title: 'The Zeta Forge: In Silico Factory',
    subtitle: 'Our Unfair Advantage: The 1M Token Context Window.',
    titleClassName: "from-purple-500 to-pink-400",
    content: {
      type: 'process-flow',
      useEnhancedLayout: true,
      steps: [
        {
          icon: Shield,
          title: 'Input: Validated Target',
          description: 'A pathogenic variant from the Zeta Oracle becomes the mission objective.',
          bgClass: 'bg-green-500/20',
          borderClass: 'border-green-500',
          textClass: 'text-green-400',
          titleClass: 'text-green-400'
        },
        {
          icon: Bot,
          title: 'The Zeta Forge',
          description: 'Our generative AI, with its massive 1M token context, designs a portfolio of therapeutic candidates.',
          bgClass: 'bg-purple-500/20',
          borderClass: 'border-purple-500',
          textClass: 'text-purple-400',
          titleClass: 'text-purple-400'
        },
        {
          icon: TestTube2,
          title: 'Output: The Arsenal',
          description: 'The result is a diverse set of in silico validated weapons, from CRISPR payloads to novel biologics.',
          bgClass: 'bg-slate-700/50',
          borderClass: 'border-slate-600',
          textClass: 'text-slate-300',
          titleClass: 'text-slate-300'
        }
      ]
    }
  },

  // Kill Chain Boltz
  {
    title: 'Step 3: In Silico Validation',
    subtitle: 'Every therapeutic is battle-tested before it\'s built.',
    titleClassName: "from-orange-500 to-yellow-400",
    content: {
      type: 'simple-block',
      useEnhancedLayout: true,
      block: {
        icon: Cuboid,
        mainText: 'A sequence is not a therapy. Our **Zeta Boltz** engine runs every designed candidate through an in silico firing range, simulating its 3D interaction with the target to predict binding affinity and efficacy.',
        subText: '**For Biotech Partners:** This provides the critical, structural proof of mechanism, dramatically de-risking the candidate before committing to expensive lab synthesis.',
        iconColor: "text-orange-400",
        borderColor: "border-slate-700"
      }
    }
  },

  // Kill Chain Boltz Detail
  {
    title: 'The Zeta Boltz: In Silico Firing Range',
    subtitle: 'From a 1D Blueprint to a 3D Proof of Victory.',
    titleClassName: "from-orange-500 to-yellow-400",
    content: {
      type: 'process-flow',
      useEnhancedLayout: true,
      steps: [
        {
          icon: Bot,
          title: 'Input: Forged Weapon',
          description: 'A novel nanobody sequence, generated by the Zeta Forge.',
          bgClass: 'bg-purple-500/20',
          borderClass: 'border-purple-500',
          textClass: 'text-purple-400',
          titleClass: 'text-purple-400'
        },
        {
          icon: Cuboid,
          title: 'The Simulation',
          description: 'Our AlphaFold 3-powered engine simulates the 3D protein-protein interaction between our weapon and its target.',
          bgClass: 'bg-orange-500/20',
          borderClass: 'border-orange-500',
          textClass: 'text-orange-400',
          titleClass: 'text-orange-400'
        },
        {
          icon: Shield,
          title: 'The Verdict',
          description: 'The result is a quantitative Binding Affinity Score, providing definitive proof of the weapon\'s physical lethality.',
          bgClass: 'bg-green-500/20',
          borderClass: 'border-green-500',
          textClass: 'text-green-400',
          titleClass: 'text-green-400'
        }
      ]
    }
  },

  // Kill Chain Asset
  {
    title: 'The Deliverable: A De-Risked Asset',
    subtitle: 'We don\'t deliver data. We deliver a validated, pre-clinical asset.',
    titleClassName: "from-green-500 to-teal-400",
    content: {
      type: 'simple-block',
      useEnhancedLayout: true,
      block: {
        icon: Package,
        mainText: 'The final output of our in silico kill chain is not a report; it is a **de-risked, high-value therapeutic asset** with a complete dossier of predictive data.',
        subText: '**For Biotech Partners:** We give you a candidate that has already won the digital war, dramatically increasing its probability of victory on the clinical battlefield.',
        iconColor: "text-green-400",
        borderColor: "border-slate-700"
      }
    }
  },

  // Kill Chain Asset Detail
  {
    title: 'The Therapeutic Dossier',
    subtitle: 'The final output of our in silico conquest.',
    titleClassName: "from-green-500 to-teal-400",
    content: {
      type: 'asset-dossier',
      useEnhancedLayout: true,
      assetId: 'Asset: CS-RUNX1-GC-001',
      status: 'Ready for Wet-Lab',
      checkpoints: [
        {
          icon: Target,
          iconColor: 'text-cyan-400',
          text: '**Target Validation:** <span class="font-mono text-green-400">COMPLETE</span>'
        },
        {
          icon: Bot,
          iconColor: 'text-purple-400',
          text: '**Weapon Design:** <span class="font-mono text-green-400">COMPLETE</span>'
        },
        {
          icon: Cuboid,
          iconColor: 'text-orange-400',
          text: '**Structural Validation:** <span class="font-mono text-green-400">COMPLETE</span>'
        }
      ],
      description: 'This dossier contains the full sequence data, in silico efficacy and safety scores, and structural binding predictions, providing our partners with a de-risked asset with a high probability of clinical success.'
    }
  },

  // Enhanced Zeta Oracle Deep Dive
  {
    title: 'The Zeta Oracle: Clinical VUS Resolution',
    subtitle: 'From Genetic Ambiguity to Actionable Intelligence',
    titleClassName: "from-cyan-400 to-blue-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-cyan-900/20 to-slate-900',
    content: {
      type: 'two-hit-hypothesis',
      steps: [
        {
          title: 'VUS Discovery',
          subtext: 'Unknown Impact & Clinical Paralysis',
          colorClass: 'bg-gradient-to-br from-yellow-400 to-orange-500',
          mutationIcon: '❓'
        },
        {
          title: 'Oracle Analysis',
          subtext: 'AI Grammar Check on Biology',
          colorClass: 'bg-gradient-to-br from-cyan-400 to-blue-600',
          mutationIcon: '🧠',
          animationClass: 'animate-pulse'
        },
        {
          title: 'Zeta Score',
          subtext: 'Quantitative Functional Impact',
          colorClass: 'bg-gradient-to-br from-green-400 to-emerald-600',
          mutationIcon: '✅'
        },
        {
          title: 'Clinical Action',
          subtext: 'Validated Target Ready for Design',
          colorClass: 'bg-gradient-to-br from-purple-400 to-pink-600',
          mutationIcon: '🎯'
        }
      ],
      // siteBlocks: [
      //   // Only show the most essential block to prevent overlapping
      //   { kind: 'oracle-explain', props: crispro101Content.oracle.explain }
      // ]
    },
    notes: "The Zeta Oracle doesn't just check databases—it understands the fundamental grammar of biology. This allows it to predict the functional impact of any variant, even those never seen before."
  },

  // Enhanced Zeta Forge Deep Dive
  {
    title: 'The Zeta Forge: Therapeutic Engineering',
    subtitle: 'From Validated Target to Multi-Modal Arsenal',
    titleClassName: "from-purple-400 to-pink-300",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-purple-900/20 to-slate-900',
    content: {
      type: 'therapeutic-arsenal',
      input: { 
        icon: '🎯', 
        title: 'Validated Target', 
        iconBg: 'bg-cyan-500/20', 
        borderColor: 'border-cyan-500', 
        accentColor: 'text-cyan-400' 
      },
      process: { 
        icon: '🔨', 
        title: 'Zeta Forge Engine', 
        borderColor: 'border-purple-400/50' 
      },
      outputs: [
        { icon: '🧬', text: 'Precision Gene Correction' },
        { icon: '✂️', text: 'Synthetic Lethal Payload' },
        { icon: '🛡️', text: 'Novel Nanobody Inhibitor' },
        { icon: '💊', text: 'Small Molecule Modulator' }
      ],
      // siteBlocks: toForgeBlocks(crispro101Content)
    },
    notes: "Our 1M token context window provides an unfair advantage: we can design ultra-long homology arms and complex therapeutic architectures that are impossible for smaller models."
  },

  // Enhanced Boltz Deep Dive  
  {
    title: 'Zeta Boltz: Structural Validation Engine',
    subtitle: 'From 1D Blueprint to 3D Proof of Mechanism',
    titleClassName: "from-orange-500 to-yellow-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-orange-900/20 to-slate-900',
    content: {
      type: 'gene-correction',
      problem: { 
        icon: '🧪', 
        title: 'The Challenge', 
        subtext: 'Sequence ≠ Function', 
        borderColor: 'border-yellow-500', 
        accentColor: 'bg-yellow-500/20', 
      },
      outcome: { 
        icon: '🏆', 
        title: 'The Solution', 
        subtext: 'Structural Proof', 
        borderColor: 'border-green-500', 
        accentColor: 'bg-green-500/20', 
        textColor: 'text-green-400' 
      },
      infoBoxes: [
        { 
          title: 'AlphaFold 3 Integration', 
          text: 'Our engine leverages the latest protein folding models to simulate 3D interactions between designed therapeutics and their targets, providing quantitative binding affinity predictions.', 
          bgClass: 'bg-slate-800/50', 
          borderColor: 'border-orange-500/30', 
          textColor: 'text-orange-400' 
        },
        { 
          title: 'In Silico Clinical Validation', 
          text: 'Every therapeutic candidate is battle-tested in our digital firing range before a single dollar is spent on wet-lab synthesis. This dramatically de-risks the development process.', 
          bgClass: 'bg-slate-800/50', 
          borderColor: 'border-green-500/30', 
          textColor: 'text-green-400' 
        }
      ],
      // siteBlocks: toBoltzBlocks(crispro101Content)
    },
    notes: "Zeta Boltz transforms drug discovery from an art into an engineering discipline. We don't hope for binding—we engineer it."
  },

  // Risk Prediction & Resistance Modeling
  {
    title: 'Predictive Disease Evolution Modeling',
    subtitle: 'Anticipating Mutations to Design Future-Proof Therapies',
    titleClassName: "from-indigo-400 to-purple-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-indigo-900/20 to-slate-900',
    content: {
      type: 'risk-prediction-map',
      knownThreat: { 
        icon: '🧬', 
        title: 'Known Risk Factor', 
        subtext: 'RUNX1 Mutation', 
        iconBg: 'bg-red-500/20', 
        borderColor: 'border-red-500', 
        accentColor: 'text-red-400' 
      },
      aiCore: { 
        icon: '🧠', 
        title: 'Predictive Engine', 
        borderColor: 'border-cyan-400/50' 
      },
      predictions: [
        { name: 'ASXL1', risk: 'High Risk (-15k)', colorClass: 'text-red-400' },
        { name: 'TET2', risk: 'Med Risk (-12k)', colorClass: 'text-orange-400' },
        { name: 'DNMT3A', risk: 'Low Risk (-9k)', colorClass: 'text-yellow-400' },
        { name: 'IDH2', risk: 'Emerging Risk', colorClass: 'text-blue-400' }
      ]
    },
    notes: "Instead of reacting to disease evolution, we predict it. This allows us to design therapies that are effective today and resilient against future resistance mutations."
  },

  // FDA Approval Strategy
  {
    title: 'Accelerating FDA Approval',
    subtitle: 'Digital Evidence for Regulatory Success',
    titleClassName: "from-blue-400 to-teal-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-blue-900/30 to-slate-900',
    content: {
      type: 'approval-process',
      dossier: [
        { 
          title: 'Zeta Score Evidence', 
          subtitle: 'Quantified Target Validation', 
          bgClass: 'bg-gradient-to-br from-cyan-600 to-blue-700', 
          borderClass: 'border-2 border-cyan-400/50' 
        },
        { 
          title: 'Design Rationale', 
          subtitle: 'AI-Engineered Solutions', 
          bgClass: 'bg-gradient-to-br from-purple-600 to-indigo-700', 
          borderClass: 'border-2 border-purple-400/50' 
        },
        { 
          title: 'Structural Validation', 
          subtitle: 'In Silico Proof of Mechanism', 
          bgClass: 'bg-gradient-to-br from-orange-600 to-red-700', 
          borderClass: 'border-2 border-orange-400/50' 
        }
      ],
      fdaTiers: [
        { title: 'Tier 3: Case Reports & Observational Studies', bgClass: 'bg-red-500/20', textColor: 'text-red-300' },
        { title: 'Tier 2: Cohort Studies & Registry Data', bgClass: 'bg-yellow-500/20', textColor: 'text-yellow-300' },
        { title: 'Tier 1: Randomized Controlled Trials', bgClass: 'bg-green-500/20', textColor: 'text-green-300' }
      ],
      fdaText: 'Our comprehensive digital dossier provides evidence across all FDA tiers, dramatically accelerating the path to approval and reducing regulatory risk.'
    }
  },

  // BUSINESS OVERVIEW SLIDE 1: PLATFORM & VALUE PROP
  {
    title: "CrisPRO.ai: Complete AI Therapeutic Design Platform",
    subtitle: "End-to-end ecosystem transforming drug development from gamble to science",
    titleClassName: "from-purple-500 via-pink-400 to-red-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-purple-900/20 to-slate-900',
    content: {
      type: 'business-overview-platform',
      useEnhancedLayout: true,
      platform: {
        title: "The CrisPRO.ai Platform",
        description: "End-to-end AI therapeutic design ecosystem - from genetic insight to clinical candidate",
        capabilities: [
          { icon: "🔬", text: "Research‑mode variant insights (S/P/E + provenance)" },
          { icon: "🧬", text: "Generative proposals (CRISPR, antibodies, small molecules) — RUO" },
          { icon: "📋", text: "Automated IND document scaffolding (research‑mode)" },
          { icon: "💰", text: "IP monetization via co‑invention/royalty doctrine" }
        ]
      },
      valueProposition: {
        title: "Proven Value Proposition",
        metrics: [
          { value: "Faster", label: "In‑silico iteration speed", className: "text-green-400" },
          { value: "Lower", label: "Operational analysis cost", className: "text-cyan-400" },
          { value: "Transparent", label: "Confidence with audit trails", className: "text-purple-400" },
          { value: "Large", label: "Addressable market opportunity", className: "text-orange-400" }
        ]
      }
    },
    notes: "First business slide: Show what we do and the massive value we deliver to pharma companies."
  },

  // BUSINESS OVERVIEW SLIDE 2: BUSINESS MODEL & MARKET
  {
    title: "CrisPRO.ai: Revenue Model & Market Opportunity",
    subtitle: "Multi-billion dollar market with proven monetization strategy",
    titleClassName: "from-green-500 via-teal-400 to-cyan-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-green-900/20 to-slate-900',
    content: {
      type: 'business-overview-business',
      useEnhancedLayout: true,
      businessModel: {
        title: "Multi-Stream Revenue Model",
        streams: [
          { name: "Platform Subscriptions", revenue: "$50K-200K/year per pharma", icon: "💳" },
          { name: "IND Generation", revenue: "$10K-50K per therapeutic", icon: "📋" },
          { name: "Platform Royalty", revenue: "2–5% of licensing revenue", icon: "💰" },
          { name: "Co‑Inventor Ownership", revenue: "10–30% patent ownership (case‑dependent)", icon: "🧾" },
          { name: "IP Licensing", revenue: "$50M+ per therapeutic", icon: "🎯" }
        ]
      },
      marketOpportunity: {
        title: "Massive Market Opportunity",
        stats: [
          { value: "$50B+", label: "Global Drug Discovery Market" },
          { value: "$8B+", label: "CRISPR Therapeutics (2025)" },
          { value: "95%", label: "Clinical Trial Failure Rate" },
          { value: "$2.8B", label: "Cost Per Approved Drug" }
        ]
      }
    },
    notes: "Second business slide: Show how we make money and the size of the opportunity."
  },

  // IND PACKAGE GENERATION SLIDE
  {
    title: "IND Package Generation: Regulatory Revolution",
    subtitle: "Complete FDA-compliant documentation in under 10 minutes",
    titleClassName: "from-blue-500 to-cyan-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-blue-900/20 to-slate-900',
    content: {
      type: 'ind-package-comparison',
      useEnhancedLayout: true,
      traditional: {
        title: "Traditional IND Process",
        time: "6-12 months",
        cost: "$500K-$2M",
        steps: [
          "Manual protocol writing",
          "Statistical analysis planning",
          "Toxicology report compilation",
          "Regulatory strategy development",
          "Multiple review cycles"
        ]
      },
      crispro: {
        title: "CrisPRO.ai IND Generation",
        time: "Automated draft in minutes (research‑mode)",
        cost: "Lower cost via automation",
        steps: [
          "S/P/E analysis complete (Baseline profile)",
          "Automated document scaffolding",
          "Standards‑aligned templates",
          "Operator QA review",
          "Provenance‑rich export"
        ],
        metrics: [
          { value: "Provenance", label: "Audit trail in outputs", className: "text-cyan-400" },
          { value: "Templates", label: "Standards‑aligned structure", className: "text-green-400" },
          { value: "Profiles", label: "Baseline/Richer/Fusion flags", className: "text-green-400" }
        ],
        mmGuidance: {
          title: "Example MM guidance (research‑mode)",
          genes: [
            { name: "KRAS/NRAS", hotspots: "G12D/V/C/S/A, G13D", guidance: "MEK inhibitor (off-label)" },
            { name: "BRAF", hotspots: "V600E, V600K", guidance: "BRAF/MEK inhibitor" },
            { name: "FGFR3", hotspots: "R248C, Y373C", guidance: "FGFR-directed agents" },
            { name: "TP53", hotspots: "R175H, R248Q/W, R273C/H", guidance: "Risk assessment & combinations" }
          ]
        }
      },
      impact: {
        title: "Business Impact",
        savings: [
          { value: "Lower", label: "Documentation effort", className: "text-green-400" },
          { value: "Faster", label: "Time to decision", className: "text-cyan-400" },
          { value: "Repeatable", label: "Deterministic export", className: "text-purple-400" }
        ]
      }
    },
    notes: "Show the dramatic improvement over traditional IND processes. Investors will immediately understand the massive cost and time savings."
  },

  // CRISPR DESIGN ECOSYSTEM SLIDE
  {
    title: "CRISPR Design Ecosystem: Complete Therapeutic Platform",
    subtitle: "From genetic target to clinical candidate in weeks, not years",
    titleClassName: "from-green-500 to-teal-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-green-900/20 to-slate-900',
    content: {
      type: 'crispr-ecosystem-comparison',
      useEnhancedLayout: true,
      traditional: {
        title: "Traditional CRISPR Development",
        time: "12-24 months",
        cost: "$5M-$15M",
        steps: [
          "Manual gRNA design",
          "Trial-and-error optimization",
          "Limited homology arms",
          "Iterative testing cycles",
          "Manual documentation"
        ]
      },
      crispro: {
        title: "CrisPRO.ai CRISPR Ecosystem",
        time: "Accelerated timelines (weeks)",
        cost: "Lower cost via automation",
        steps: [
          "AI-optimized gRNA design",
          "Long homology arms (context‑aware)",
          "Structural validation",
          "Regulatory documentation",
          "Ready for preclinical testing"
        ],
        apis: [
          { name: "/api/insights/predict_protein_functionality_change", description: "Functionality insight (research‑mode)" },
          { name: "/api/design/generate_guide_rna", description: "Guide candidates with safety gates" },
          { name: "/api/design/generate_repair_template", description: "HDR template proposals (planned)" },
          { name: "/api/efficacy/predict", description: "Ranked therapy hypotheses" }
        ]
      },
      marketImpact: {
        title: "Market Opportunity",
        stats: [
          { value: "$8B+", label: "CRISPR Therapeutics Market (2025)", className: "text-green-400" },
          { value: "25x", label: "Faster Development", className: "text-cyan-400" },
          { value: "95%", label: "Cost Reduction", className: "text-purple-400" },
          { value: "10+", label: "Pipeline Opportunities", className: "text-orange-400" }
        ]
      },
      competitiveAdvantages: {
        title: "Why Our CRISPR Platform Wins",
        advantages: [
          { icon: "🧬", text: "1M token context window enables ultra-long, complex homology arms", color: "green" },
          { icon: "🎯", text: "Fusion Engine integration for superior target validation", color: "blue" },
          { icon: "⚡", text: "Modal microservices ensure production-ready performance", color: "purple" },
          { icon: "📋", text: "Complete regulatory documentation package included", color: "orange" }
        ]
      }
    },
    notes: "Position our CRISPR platform as a complete solution vs. fragmented tools. Show the dramatic time and cost savings that will resonate with investors."
  },

  // THE QUALCOMM OF PHARMA SLIDE
  {
    title: "The Qualcomm of Pharma: Platform Royalty Model",
    subtitle: "Recurring revenue from every drug that uses our AI technology",
    titleClassName: "from-purple-500 to-pink-400",
    backgroundClass: 'bg-gradient-to-br from-slate-900 via-purple-900/20 to-slate-900',
    content: {
      type: 'qualcomm-comparison',
      useEnhancedLayout: true,
      comparison: {
        qualcomm: {
          title: "Qualcomm's Model",
          icon: "📱",
          revenue: "$8.6B annual revenue",
          elements: [
            "Baseband processor licensing",
            "Royalty: 5% of device price",
            "Every smartphone pays Qualcomm",
            "Continuous tech evolution"
          ]
        },
        crispro: {
          title: "CrisPRO.ai Model",
          icon: "🧬",
          revenue: "Platform royalty + IP co‑invention",
          elements: [
            "AI therapeutic design licensing",
            "Platform royalty: 2–5% of licensing",
            "Every approved drug pays us",
            "Continuous AI model improvement"
          ]
        }
      },
      revenueStreams: {
        title: "Multiple Revenue Streams",
        streams: [
          {
            name: "IND Generation",
            description: "One-time fee per drug program (Fusion Engine analysis included)",
            potential: "$10K-50K per drug",
            icon: "📋"
          },
          {
            name: "Platform Subscription",
            description: "Monthly/annual platform access (research‑mode guidance)",
            potential: "$50K-200K per year",
            icon: "💳"
          },
          {
            name: "Platform Royalty",
            description: "2–5% of licensing revenue when platform contributes materially",
            potential: "Meaningful annuity per licensed asset",
            icon: "💰"
          },
          {
            name: "Co‑Inventor/IP Monetization",
            description: "10–30% patent ownership (case‑dependent) and licensing",
            potential: "$50M+ per therapeutic",
            icon: "🎯"
          }
        ]
      },
      competitiveAdvantages: {
        title: "Fusion Engine Differentiators",
        advantages: [
          { icon: "🔄", text: "Fused, not single-source: AM integration when eligible", color: "green" },
          { icon: "📈", text: "Provenance everywhere: Complete audit trail with MoA tags", color: "blue" },
          { icon: "🛡️", text: "Guidance-ready: Auditable confidence lifts with MoA gates", color: "purple" },
          { icon: "⚡", text: "Selective lift: Conservative defaults maintain regulatory comfort", color: "orange" }
        ]
      }
    },
    notes: "This is the most important slide for investors. Show the massive revenue potential and how our platform creates ongoing value unlike traditional biotech companies."
  }
];
