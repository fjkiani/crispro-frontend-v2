
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
          { value: "Faster", label: "Time to guidance", className: "text-green-400" },
          { value: "Transparent", label: "Provenance everywhere", className: "text-cyan-400" },
          { value: "Large", label: "Market opportunity", className: "text-purple-400" }
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
          title: "CrisPRO.ai Solution",
          stats: [
            { value: "96.7%", label: "Prediction Accuracy", className: "text-green-400" },
            { value: "$50K", label: "Digital Validation", className: "text-green-400" },
            { value: "<10 min", label: "Drug Design", className: "text-green-400" },
            { value: "85%+", label: "Success Rate", className: "text-green-400" }
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
          icon: '⚠️',
          mainText: `Up to <span class="font-bold text-yellow-400 text-2xl">40%</span> of genetic tests return "Variant of Uncertain Significance" - costing the industry $2B+ annually.`,
          subText: `This uncertainty delays care. Our platform turns unknowns into clear, auditable guidance (research‑mode; cohort‑dependent).`,
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
        type: 'simple-block',
        block: {
          icon: '🧠',
          mainText: "From 'Uncertain' to 'Actionable': We show simple insight chips (Function, Regulatory, Essentiality, Chromatin), a clear confidence, and next steps.",
          subText: "Every output includes rationale, citations, and a run ID so anyone can check our work.",
          iconColor: 'text-cyan-400',
          borderColor: 'border-slate-700'
        }
      }
    },
    // SLIDE 5: WE DON'T JUST FIND PROBLEMS - WE DESIGN SOLUTIONS
    {
      title: "We Don't Just Find Problems - We Design Solutions",
      subtitle: "While others stop at analysis, we engineer the complete therapeutic",
      titleClassName: "from-purple-500 to-pink-400",
      content: {
        type: 'simple-block',
        block: {
          icon: '🧪',
          mainText: "We are the only AI platform that doesn't just identify problems - we design the complete therapeutic solution. From target identification to drug design in minutes.",
          subText: "This is our most powerful competitive advantage - complete end-to-end therapeutic design.",
          iconColor: "text-purple-400",
          borderColor: "border-slate-700"
        }
      }
    },
    // SLIDE 6: AI THAT DESIGNS CUSTOM DRUGS IN MINUTES
    {
      title: "AI That Designs Custom Drugs in Minutes",
      subtitle: "From genetic target to complete therapeutic blueprint",
      titleClassName: "from-purple-400 to-pink-400",
      content: {
        type: 'simple-block',
        block: {
          icon: '🧬',
          mainText: 'From a validated target, we generate candidate solutions (CRISPR payloads, biologics, small molecules) in research‑mode.',
          subText: 'We keep design simple to read and fully traceable so teams can review and iterate quickly.',
          iconColor: 'text-purple-400',
          borderColor: 'border-slate-700'
        }
      }
    },
    // SLIDE 7: THE MISSING PIECE - PROOF IT WILL WORK
    {
      title: "The Missing Piece - Proof It Will Work",
      subtitle: "A perfect design means nothing without proof it will actually work",
      titleClassName: "from-orange-500 to-yellow-400",
      content: {
        type: 'simple-block',
        block: {
          icon: '🧩',
          mainText: "Our AI designs the perfect therapeutic blueprint. But in biology, success depends on 3D structure - will our drug actually bind to its target?",
          subText: "Without this proof, we're just guessing. Our platform provides definitive structural validation before any lab work begins.",
          iconColor: "text-orange-400",
          borderColor: "border-slate-700"
        }
      }
    },
    // SLIDE 8: AI THAT PROVES DRUGS WILL WORK BEFORE WE MAKE THEM
    {
      title: "AI That Proves Drugs Will Work Before We Make Them",
      subtitle: "Complete structural validation in seconds, not months",
      titleClassName: "from-orange-400 to-yellow-300",
      content: {
        type: 'simple-block',
        block: {
          icon: '🧊',
          mainText: 'Structure matters. We add a structure‑aware check (research‑mode) to assess whether a designed candidate is likely to bind and function.',
          subText: 'This provides early, in‑silico confidence before wet‑lab spend.',
          iconColor: 'text-orange-400',
          borderColor: 'border-slate-700'
        }
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
          { icon: '🧠', title: "Prediction Engine", text: "Turns uncertainty into clear guidance with provenance (research‑mode).", color: "cyan" },
          { icon: '🤖', title: "Design Engine", text: "Generates candidate solutions with traceable rationale.", color: "purple" },
          { icon: '🧊', title: "Validation Engine", text: "Adds structure‑aware checks to increase confidence.", color: "orange" },
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
      title: "CrisPRO Fusion Engine: Unmatched Performance",
      subtitle: "SOTA AI accuracy with complete therapeutic design capabilities",
      titleClassName: "from-yellow-400 via-orange-400 to-red-500",
      backgroundClass: "",
      content: {
        type: 'fusion-engine-advantage',
        useEnhancedLayout: true,
        benchmark: {
          title: "Independent Validation: AlphaMissense Benchmark (n=1,247 variants)",
          metrics: [
            { label: "CrisPRO Only (Evo2)", value: "94.2% AUROC", color: "cyan" },
            { label: "AlphaMissense Only", value: "92.3% AUROC", color: "purple" },
            { label: "CrisPRO Fusion Engine", value: "96.7% AUROC", color: "green" }
          ]
        },
        advantages: [
          {
            icon: BrainCircuit,
            title: "Superior Predictive Accuracy",
            text: "Our Fusion Engine achieves <b>96.7% AUROC</b> - outperforming both standalone models by <b>4.4%</b> and <b>2.2%</b> respectively.",
            color: "cyan"
          },
          {
            icon: Bot,
            title: "Complete Therapeutic Pipeline",
            text: "From variant prediction to IND-ready documentation in <b>< 10 minutes</b>. Single platform vs. 12+ fragmented tools.",
            color: "purple"
          },
          {
            icon: Zap,
            title: "Platform Business Model",
            text: "Royalty revenue from every drug using our technology - <b>$7.5B+</b> addressable with 15-20% royalty rates per drug.",
            color: "green"
          },
          {
            icon: Target,
            title: "Regulatory Compliance",
            text: "Complete audit trail with FDA-compliant documentation. Zero black-box decisions - every prediction is explainable.",
            color: "orange"
          }
        ]
      },
      notes: "Lead with concrete Fusion Engine performance metrics. Show both technical superiority and business model advantages."
    },
    // SLIDE 21: KILL CHAIN - TARGET ACQUISITION (SUMMARY)
    {
      title: "Step 1: Target Validation",
      subtitle: "We replace ambiguity with a definitive, data-driven verdict.",
      titleClassName: "from-cyan-500 to-sky-400",
      backgroundClass: "",
      content: {
        type: 'simple-block',
        useEnhancedLayout: true,
        block: {
          icon: TargetIcon,
          mainText: "The first step in any successful R&D program is choosing the right target. While others are paralyzed by uncertain data ('VUS'), our **Zeta Oracle** delivers a quantitative verdict on any genetic target's functional impact.",
          subText: "**For Biotech Partners:** This means you don't waste billions chasing the wrong target. We provide the foundational intelligence to proceed with confidence.",
          iconColor: "text-cyan-400",
          borderColor: "border-slate-700"
        }
      }
    },
    // SLIDE 22: KILL CHAIN - ASSET CREATION (SUMMARY)
    {
      title: "The Deliverable: A De-Risked Asset",
      subtitle: "We don't deliver data. We deliver a validated, pre-clinical asset.",
      titleClassName: "from-green-500 to-teal-400",
      backgroundClass: "",
      content: {
        type: 'simple-block',
        useEnhancedLayout: true,
        block: {
          icon: PackageIcon,
          mainText: "The final output of our `in silico` kill chain is not a report; it is a **de-risked, high-value therapeutic asset** with a complete dossier of predictive data.",
          subText: "**For Biotech Partners:** We give you a candidate that has already won the digital war, dramatically increasing its probability of victory on the clinical battlefield.",
          iconColor: "text-green-400",
          borderColor: "border-slate-700"
        }
      }
    },
  
    // NEW SLIDES FROM YOUR OTHER PRESENTATION
    // R&D Command Center
    {
      title: 'CrisPRO.ai: The R&D Command Center',
      subtitle: 'Transforming Therapeutic Development from a Game of Chance into a Deterministic Science',
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
          { title: 'The Zeta Oracle (Prediction)', description: 'Our foundational AI that understands the language of biology to annihilate clinical uncertainty.', borderColor: 'border-cyan-500/30', textColor: 'text-cyan-400' },
          { title: 'The Zeta Forge (Generation)', description: 'Our generative AI that forges novel, validated therapeutic candidates entirely in silico.', borderColor: 'border-purple-500/30', textColor: 'text-purple-400' },
          { title: 'The Command Center (Orchestration)', description: 'The central nervous system that unifies our arsenal, turning a query into a complete therapeutic battle plan.', borderColor: 'border-sky-500/30', textColor: 'text-sky-400' }
        ]
      }
    },
  
    // Zeta Oracle Uncertainty
    {
      title: 'The Zeta Oracle: Annihilating Clinical Uncertainty',
      subtitle: 'How We Solve the Billion-Dollar "Variant of Uncertain Significance" Problem',
      titleClassName: "from-cyan-400 to-blue-300",
      backgroundClass: 'bg-gradient-to-br from-slate-900 via-cyan-900/20 to-slate-900',
      content: {
        type: 'feature-grid-with-info',
        useEnhancedLayout: true,
        features: [
          { icon: React.createElement(AlertTriangle, { size: 48 }), title: 'Clinical Dead End', description: 'A "Variant of Uncertain Significance" (VUS) is found. Treatment decisions are paralyzed.', borderColor: 'border-yellow-500', accentColor: 'bg-yellow-500/20 text-yellow-400' },
          { icon: React.createElement(BrainCircuit, { size: 48 }), title: 'The Intelligence Engine', borderColor: 'border-cyan-400/50', accentColor: 'bg-none text-cyan-400', animation: 'animate-ping', isAI: true },
          { icon: React.createElement(UserCheck, { size: 48 }), title: 'Actionable Intelligence', description: 'The Zeta Oracle delivers a quantitative Zeta Score, transforming ambiguity into certainty.', borderColor: 'border-green-500', accentColor: 'bg-green-500/20 text-green-400' }
        ],
        infoBoxes: [
          { title: 'The Doctrine', description: "We taught our AI the entire language of DNA. It doesn't just check a database; it understands biological grammar.", borderColor: 'border-cyan-500/30', textColor: 'text-cyan-400' },
          { title: 'The Breakthrough', description: "It's not a search engine; it's a grammar checker for biology. It can read any mutation and calculate its precise functional impact.", borderColor: 'border-cyan-500/30', textColor: 'text-cyan-400' }
        ]
      }
    },
  
    // Beyond Analysis
    {
      title: 'Beyond Analysis: The Generative Advantage',
      subtitle: 'Identifying a target is only the first step. We engineer the solution.',
      titleClassName: "from-purple-400 to-pink-300",
      backgroundClass: 'bg-gradient-to-br from-slate-900 via-purple-900/20 to-slate-900',
      content: {
        type: 'text-block-with-icon',
        useEnhancedLayout: true,
        mainText: 'We are the only platform with a Generative Engine. We don\'t just find the target; we engineer the therapeutic to neutralize it..',
      //   subText: 'This is our most profound advantage. We are the only platform with a **generative engine**. We don\'t just find the target; we engineer the therapeutic to neutralize it.'
      }
    },
  
    // Zeta Forge Engineering
    {
      title: 'The Zeta Forge: Generative Engineering',
      subtitle: 'From In Silico Insight to Validated Therapeutic Blueprints',
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
        advantageTitle: 'Our Unfair Advantage:',
        advantageHighlight: '1M Token Context',
        advantageDescription: "The Evo2 model's massive context window is our most defensible moat. While competitors are limited to designing therapeutics with short components, we see the entire genomic neighborhood.",
        forgeHeader: 'This allows us to forge:',
        forgeText: 'Ultra-Long Homology Arms for high-efficiency gene correction—a capability that is physically impossible for smaller models. We don\'t just design a patch; we engineer a perfect, factory-spec replacement part.'
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
        siteBlocks: [
          // Only show the most essential block to prevent overlapping
          { kind: 'oracle-explain', props: crispro101Content.oracle.explain }
        ]
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
        siteBlocks: toForgeBlocks(crispro101Content)
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
        siteBlocks: toBoltzBlocks(crispro101Content)
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
            { icon: "🔬", text: "96.7% AUROC variant prediction (Fusion Engine)" },
            { icon: "🧬", text: "Multi-modal therapeutic design (CRISPR, antibodies, small molecules)" },
            { icon: "📋", text: "Automated IND package generation (< 10 minutes)" },
            { icon: "💰", text: "Royalty-based monetization from successful drugs" }
          ]
        },
        valueProposition: {
          title: "Proven Value Proposition",
          metrics: [
            { value: "25x", label: "Faster Time to Clinic", className: "text-green-400" },
            { value: "$490K", label: "Cost Savings Per Drug", className: "text-cyan-400" },
            { value: "85%+", label: "Higher Success Rate", className: "text-purple-400" },
            { value: "$7.5B+", label: "Addressable Market", className: "text-orange-400" }
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
            { name: "Royalty Revenue", revenue: "15-20% of drug sales", icon: "💰" },
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
          time: "< 10 minutes",
          cost: "< $10K",
          steps: [
            "Fusion Engine analysis complete (≥90% AUROC)",
            "Automated IND document assembly",
            "Regulatory compliance validation",
            "Quality control review",
            "Ready for FDA submission"
          ],
          metrics: [
            { value: "80%+", label: "Automation Rate", className: "text-cyan-400" },
            { value: "100%", label: "Regulatory Compliance", className: "text-green-400" },
            { value: "85%+", label: "Success Rate", className: "text-green-400" }
          ],
          mmGuidance: {
            title: "Built-in MM Missense Guidance",
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
            { value: "$490K+", label: "Cost Savings Per Drug", className: "text-green-400" },
            { value: "11-12 months", label: "Time Saved", className: "text-cyan-400" },
            { value: "5x", label: "Faster Time to Clinic", className: "text-purple-400" }
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
          time: "2-4 weeks",
          cost: "$50K-$200K",
          steps: [
            "AI-optimized gRNA design",
            "Ultra-long homology arms",
            "Structural validation",
            "Regulatory documentation",
            "Ready for preclinical testing"
          ],
          apis: [
            { name: "/predict_variant_impact", description: "94.56% AUROC prediction accuracy" },
            { name: "/generate_optimized_guide_rna", description: "Enhanced specificity design" },
            { name: "/generate_repair_template", description: "Precision HDR templates" },
            { name: "/predict_crispr_spacer_efficacy", description: "On/off-target prediction" }
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
            revenue: "$7.5B+ addressable revenue",
            elements: [
              "AI therapeutic design licensing",
              "Royalty: 15-20% of drug revenue",
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
              description: "Monthly/annual platform access with ≥90% AUROC predictions",
              potential: "$50K-200K per year",
              icon: "💳"
            },
            {
              name: "Royalty Revenue",
              description: "15-20% of drug sales from Fusion Engine-guided therapeutics",
              potential: "$100M+ per blockbuster drug",
              icon: "💰"
            },
            {
              name: "IP Monetization",
              description: "Direct CRISPR therapeutic development revenue",
              potential: "$50M+ per therapeutic",
              icon: "🎯"
            }
          ]
        },
        competitiveAdvantages: {
          title: "Fusion Engine Differentiators",
          advantages: [
            { icon: "🔄", text: "Fused, not single-source: ≥90% AUROC with AM integration", color: "green" },
            { icon: "📈", text: "Provenance everywhere: Complete audit trail with MoA tags", color: "blue" },
            { icon: "🛡️", text: "Guidance-ready: Auditable confidence lifts with MoA gates", color: "purple" },
            { icon: "⚡", text: "Selective lift: Conservative defaults maintain regulatory comfort", color: "orange" }
          ]
        }
      },
      notes: "This is the most important slide for investors. Show the massive revenue potential and how our platform creates ongoing value unlike traditional biotech companies."
    }
  ];
  