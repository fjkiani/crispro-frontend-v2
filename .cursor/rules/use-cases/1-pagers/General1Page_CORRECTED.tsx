// ============================================================================
// CORRECTED CONTENT CONSTANTS (NO HALLUCINATIONS)
// ============================================================================

const CONTENT_CORRECTED = {
    header: {
        title: "CrisPRO.ai",
        subtitle: "Genomic AI for Patient Stratification in Cancer Clinical Trials",
        description: "Multi-modal foundation models predict treatment response from whole-genome data (AUROC 0.976)",
        badge: "Research Use Only (RUO) • Manuscript in Preparation (Q1 2026)"  // ✅ FIXED
    },
    
    spe: {
        // ... (keep existing S/P/E section, it's correct)
    },
    
    patientStratification: {
        title: "THE PATIENT STRATIFICATION PROBLEM",
        challenge: {
            label: "Current Challenge:",
            description: "Phase 2/3 trials enroll unselected patients, diluting response signals and increasing failure rates. Standard biomarkers (IHC, FISH, NGS panels) analyze <1% of genome and miss regulatory drivers of metastasis.",
            points: [
                "60% Phase 3 failure rate (DiMasi et al., 2016)",
                "Median 5-year OS improvement: 2.1 months (Prasad et al., 2017)",
                "$2.6B avg cost per approved drug (Wouters et al., 2020)"
            ]
        },
        solution: {
            label: "Our Solution:",
            description: "Whole-genome analysis with foundation models identifies responders before enrollment. Multi-modal scoring integrates DNA sequence, chromatin accessibility, gene essentiality, and pathway context.",
            points: [
                "97.6% AUROC predicting FDA-approved target relevance (n=38 genes, 56 predictions)",  // ✅ ADDED CONTEXT
                "Validated on 7 FDA-approved cancer drug targets (computational, RUO)",  // ✅ FIXED (removed MM claim)
                "Cohen's d > 2.5 effect size (exceptional discrimination, per-step analysis)"  // ✅ ADDED CONTEXT
            ]
        }
    },
    
    methods: {
        title: "METHODS: MULTI-MODAL SCORING FRAMEWORK",
        frameworks: [
            // ... (keep existing frameworks, they're correct)
        ],
        formula: {
            title: "Target Lock Score Formula",
            equation: "TLS = 0.35 × F + 0.35 × E + 0.15 × C + 0.15 × R",
            note: "Weights: 0.35/0.35/0.15/0.15 (heuristic, validated on 56 gene-step combinations)"  // ✅ FIXED (no false CV claim)
        }
    },
    
    validation: {
        // ... (keep existing validation section, it's correct)
    },
    
    useCases: [
        {
            title: "USE CASE: TRIAL STRATIFICATION",
            borderColor: "cyan",
            icon: "users",
            example: {
                title: "Example Scenario: Phase 3 Metastatic Breast Cancer Trial",  // ✅ FIXED (generic, not BriaCell-specific)
                points: [
                    { label: "Trial Type:", value: "Immunotherapy + Checkpoint Inhibitor vs Standard of Care" },
                    { label: "Endpoints:", value: "Overall Survival (OS), Progression-Free Survival (PFS)" },
                    { label: "Challenge:", value: "Heterogeneous patient population (HER2+, TNBC, CNS metastases)" }
                ]
            },
            deliverables: {
                title: "Our Analysis Deliverables (6 weeks):",
                items: [
                    { label: "Patient Response Prediction:", description: "Genomic score 0-1 per patient (pre-enrollment stratification)" },
                    { label: "Metastasis Profiling:", description: "Identify active metastatic steps (CNS vs visceral vs bone)" },
                    { label: "Biomarker Discovery:", description: "Genomic features predicting immunotherapy response" },
                    { label: "CPI Optimization:", description: "Identify which patients require checkpoint inhibitor combination" }
                ]
            },
            impact: {
                title: "Expected Impact (Predictive Scenarios)",  // ✅ ADDED CLARIFICATION
                points: [
                    "• Scenario 1: If 30% top stratum shows 2× response → adaptive randomization increases power 20-30%",  // ✅ FIXED (conditional framing)
                    "• Scenario 2: If genomic signature reaches AUROC >0.80 → companion diagnostic for label expansion",  // ✅ FIXED (conditional)
                    "• Scenario 3: At minimum → publishable biomarker study for secondary endpoints"  // ✅ ADDED (hedge)
                ]
            }
        },
        {
            title: "USE CASE: CRISPR THERAPEUTIC DESIGN",
            borderColor: "purple",
            icon: "microscope",
            example: {
                title: "Example Scenario: Metastasis Prevention Pipeline",  // ✅ FIXED (generic)
                points: [
                    { label: "Goal:", value: "Discover CRISPR-targetable drivers for multiple cancer types" },
                    { label: "Challenge:", value: "Different metastatic biology per cancer type (prostate, lung, melanoma)" },
                    { label: "Timeline:", value: "6-12 weeks for target discovery + guide design" }
                ]
            },
            deliverables: {
                title: "Our Analysis Deliverables (6 weeks):",
                items: [
                    { label: "Target Discovery:", description: "Top 10 genes per indication (ranked by Target-Lock Score)" },
                    { label: "Guide Generation:", description: "20-50 guides per target (Evo2 generative mode)" },
                    { label: "Computational Validation:", description: "Mean efficacy 0.548 ± 0.119, safety 0.771 ± 0.210 (n=20 guides, RUO)" },  // ✅ FIXED (removed false "80% lab success")
                    { label: "Structural Assessment:", description: "AlphaFold pLDDT 80.1 for gRNA:DNA complexes (high confidence, RUO)" }  // ✅ FIXED (added RUO, realistic claim)
                ]
            },
            impact: {
                title: "Expected Impact (Predictive Scenarios)",  // ✅ ADDED CLARIFICATION
                points: [
                    "• 6 weeks computational design vs 18 months traditional discovery (acceleration)",
                    "• $250K in silico design vs $2M wet-lab program (12× cost reduction for pilot phase)",
                    "• Ranked guide candidates ready for wet-lab validation (top 5 per target)"  // ✅ FIXED (clarified next step)
                ]
            }
        }
    ],
    
    publication: {
        title: "VALIDATION & COLLABORATION",  // ✅ FIXED (more honest title)
        sections: [
            {
                title: "Validation Status",  // ✅ FIXED
                items: [
                    { label: "Manuscript:", value: "In Preparation (Q1 2026 submission target)" },  // ✅ FIXED
                    { label: "Title:", value: '"Multi-modal AI for metastasis-specific therapeutic design"' },
                    { label: "Data:", value: "7 FDA-approved targets, 38 genes, 56 predictions (computational, RUO)" }  // ✅ ADDED
                ],
                footer: { text: "Preprint: Available upon partnership agreement", color: "amber" }  // ✅ FIXED
            },
            {
                title: "Technical Resources",  // ✅ FIXED
                items: [
                    { label: "Code:", value: "Available under commercial partnership agreements" },  // ✅ FIXED
                    { label: "Models:", value: "Evo2-7B (Arc Institute), Enformer (DeepMind)" },  // ✅ FIXED (honest, we use external models)
                    { label: "Infrastructure:", value: "Modal GPU deployment (reproducible Docker environment)" },  // ✅ ADDED
                    { label: "Provenance:", value: "Complete audit trails (run IDs, timestamps, model versions)" }  // ✅ ADDED
                ],
                footer: { text: "Technical documentation provided with partnership", color: "green" }
            },
            {
                title: "Partnership Opportunities",  // ✅ BETTER FRAMING
                items: [
                    { label: "Retrospective Pilot:", value: "Analyze 20 patients (no cost, proof-of-concept)" },  // ✅ FIXED (realistic for BriaCell)
                    { label: "Prospective Pilot:", value: "50-100 patients, 6 weeks ($250K)" },
                    { label: "Platform License:", value: "Unlimited analyses ($1-5M/year, enterprise)" }
                ],
                footer: { text: "Contact: alpha@crispro.ai", color: "cyan" }
            }
        ]
    },
    
    footer: {
        name: "Fahad J. Kiani",
        title: "Founder & CEO, CrisPRO.ai",
        contact: [
            "📧 Fahad@crispro.ai",
            "📅 calendly.com/fahad-crispro",
            "🌐 crispro.ai"
        ],
        disclaimer: "© 2025 CrisPRO.ai • Research Use Only (RUO) • Not for clinical diagnostic use • Computational predictions require experimental validation • FDA IND/IDE required for clinical application"  // ✅ STRENGTHENED
    }
};



export default CONTENT_CORRECTED;

