import React from 'react';
import { 
    Target, Dna, BrainCircuit, Cpu, ShieldCheck, Microscope,
    Activity, BarChart3, CheckCircle, FlaskConical, BookOpen,
    TrendingUp, Users, GitBranch, Layers, AlertCircle, ArrowRight
} from 'lucide-react';

// ============================================================================
// COMPONENT CODE
// ============================================================================

const CrisPROTechnicalOnePager = () => {
    return (
        <div className="w-full min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 p-8 font-sans">
            <div className="max-w-[1200px] mx-auto bg-slate-800/40 backdrop-blur-md rounded-3xl border-2 border-slate-700 p-10 shadow-2xl">
                
                {/* HEADER */}
                <div className="text-center mb-6 pb-4 border-b-2 border-slate-700">
                    <div className="flex items-center justify-center space-x-4 mb-3">
                        <Dna className="w-12 h-12 text-cyan-400"/>
                        <h1 className="text-5xl font-black text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-blue-500">
                            {CONTENT.header.title}
                        </h1>
                    </div>
                    <h2 className="text-3xl font-bold text-white mb-2">
                        {CONTENT.header.subtitle}
                    </h2>
                    <p className="text-lg text-slate-300">
                        {CONTENT.header.description}
                    </p>
                    <p className="text-sm text-cyan-400 mt-2">
                        {CONTENT.header.badge}
                    </p>
                </div>

                {/* S/P/E DOCTRINE SECTION */}
                <div className="bg-slate-900/60 p-6 rounded-2xl border-2 border-purple-500 mb-6">
                    <div className="flex items-center justify-center mb-4">
                        <Layers className="w-10 h-10 text-purple-400 mr-3"/>
                        <h3 className="text-2xl font-bold text-purple-300">{CONTENT.spe.title}</h3>
                    </div>
                    
                    <div className="flex items-center justify-between mb-6">
                        {/* Problem Section */}
                        <div className="flex-1 bg-red-900/20 p-4 rounded-lg border border-red-500/50">
                            <p className="font-bold text-red-300 text-sm mb-3 text-center">{CONTENT.spe.problem.title}</p>
                            <div className="space-y-2 text-xs">
                                {CONTENT.spe.problem.points.map((point, i) => (
                                    <div key={i} className="flex items-center">
                                        <AlertCircle className="w-4 h-4 text-red-400 mr-2 flex-shrink-0"/>
                                        <p className="text-slate-300">{point}</p>
                                    </div>
                                ))}
                            </div>
                            <div className="mt-3 text-center">
                                <p className="text-2xl font-black text-red-300">?</p>
                                <p className="text-xs text-slate-400">{CONTENT.spe.problem.conclusion}</p>
                            </div>
                        </div>
                        
                        <ArrowRight className="w-12 h-12 text-purple-400 mx-4 flex-shrink-0"/>
                        
                        {/* Engine Section */}
                        <div className="flex-1">
                            <div className="grid grid-cols-1 gap-2">
                                {CONTENT.spe.engine.map((item, i) => (
                                    <div key={i} className={`bg-${item.color}-900/20 p-3 rounded-lg border-2 border-${item.color}-500`}>
                                        <div className="flex items-center justify-between">
                                            <div className="flex items-center">
                                                {item.icon === 'dna' && <Dna className="w-6 h-6 text-cyan-400 mr-2"/>}
                                                {item.icon === 'branch' && <GitBranch className="w-6 h-6 text-green-400 mr-2"/>}
                                                {item.icon === 'book' && <BookOpen className="w-6 h-6 text-purple-400 mr-2"/>}
                                                <div>
                                                    <p className="font-bold text-white text-sm">{item.label}</p>
                                                    <p className="text-xs text-slate-400">{item.description}</p>
                                                </div>
                                            </div>
                                            <p className={`text-${item.color}-300 font-mono text-xs`}>{item.weight}</p>
                                        </div>
                                    </div>
                                ))}
                            </div>
                            
                            <div className="flex justify-center my-2">
                                <div className="w-0 h-0 border-l-8 border-r-8 border-t-8 border-l-transparent border-r-transparent border-t-purple-500"></div>
                            </div>
                            
                            <div className="bg-purple-800 p-2 rounded-lg text-center">
                                <p className="text-white font-bold text-sm">{CONTENT.spe.efficacyScore.title}</p>
                                <p className="text-purple-200 text-xs">{CONTENT.spe.efficacyScore.description}</p>
                            </div>
                        </div>
                        
                        <ArrowRight className="w-12 h-12 text-purple-400 mx-4 flex-shrink-0"/>
                        
                        {/* Deliverable Section */}
                        <div className="flex-1 bg-green-900/20 p-4 rounded-lg border border-green-500/50">
                            <p className="font-bold text-green-300 text-sm mb-3 text-center">{CONTENT.spe.deliverable.title}</p>
                            <div className="space-y-1 text-xs">
                                {CONTENT.spe.deliverable.results.map((result, i) => (
                                    <div key={i} className={`bg-slate-800/60 p-2 rounded flex justify-between items-center ${result.opacity ? 'opacity-50' : ''}`}>
                                        <p className={result.opacity ? "text-slate-400" : "text-white font-semibold"}>{result.label}</p>
                                        <p className={`${result.color} font-mono`}>{result.score}</p>
                                    </div>
                                ))}
                            </div>
                            <div className="mt-3 text-center">
                                <p className="text-xs text-slate-400">{CONTENT.spe.deliverable.footer.provenance}</p>
                                <p className="text-xs text-green-300 font-semibold">{CONTENT.spe.deliverable.footer.runId}</p>
                            </div>
                        </div>
                    </div>
                    
                    <div className="bg-purple-800/30 p-4 rounded-lg text-center">
                        <p className="text-white font-bold text-sm mb-1">{CONTENT.spe.summary.main}</p>
                        <p className="text-slate-300 text-xs">{CONTENT.spe.summary.details}</p>
                    </div>
                </div>
                
                {/* PATIENT STRATIFICATION PROBLEM */}
                <div className="bg-red-900/20 p-6 rounded-2xl border-2 border-red-500/50 mb-6">
                    <div className="flex items-center mb-4">
                        <FlaskConical className="w-10 h-10 text-red-400 mr-3"/>
                        <h3 className="text-2xl font-bold text-red-300">{CONTENT.patientStratification.title}</h3>
                    </div>
                    <div className="grid grid-cols-2 gap-6">
                        <div>
                            <p className="text-slate-300 text-sm mb-3">
                                <span className="text-white font-semibold">{CONTENT.patientStratification.challenge.label}</span> {CONTENT.patientStratification.challenge.description}
                            </p>
                            <div className="space-y-2 text-xs">
                                {CONTENT.patientStratification.challenge.points.map((point, i) => (
                                    <div key={i} className="flex items-start">
                                        <div className="w-2 h-2 bg-red-400 rounded-full mr-2 mt-1 flex-shrink-0"/>
                                        <p className="text-slate-300">{point}</p>
                                    </div>
                                ))}
                            </div>
                        </div>
                        <div>
                            <p className="text-slate-300 text-sm mb-3">
                                <span className="text-white font-semibold">{CONTENT.patientStratification.solution.label}</span> {CONTENT.patientStratification.solution.description}
                            </p>
                            <div className="space-y-2 text-xs">
                                {CONTENT.patientStratification.solution.points.map((point, i) => (
                                    <div key={i} className="flex items-start">
                                        <div className="w-2 h-2 bg-green-400 rounded-full mr-2 mt-1 flex-shrink-0"/>
                                        <p className="text-slate-300">{point}</p>
                                    </div>
                                ))}
                            </div>
                        </div>
                    </div>
                </div>
                
                {/* METHODS SECTION */}
                <div className="bg-slate-900/60 p-6 rounded-2xl border-2 border-cyan-500 mb-6">
                    <div className="flex items-center mb-4">
                        <Cpu className="w-10 h-10 text-cyan-400 mr-3"/>
                        <h3 className="text-2xl font-bold text-cyan-300">{CONTENT.methods.title}</h3>
                    </div>
                    
                    <div className="grid grid-cols-4 gap-3 mb-4">
                        {CONTENT.methods.frameworks.map((framework, i) => (
                            <div key={i} className={`bg-${framework.color}-900/20 p-4 rounded-lg border border-${framework.color}-500/30`}>
                                <div className="flex items-center mb-2">
                                    {framework.icon === 'dna' && <Dna className="w-6 h-6 text-cyan-400 mr-2"/>}
                                    {framework.icon === 'target' && <Target className="w-6 h-6 text-green-400 mr-2"/>}
                                    {framework.icon === 'activity' && <Activity className="w-6 h-6 text-purple-400 mr-2"/>}
                                    {framework.icon === 'branch' && <GitBranch className="w-6 h-6 text-orange-400 mr-2"/>}
                                    <p className="font-bold text-white text-sm">{framework.name}</p>
                                </div>
                                <p className="text-slate-300 text-xs mb-2">
                                    <span className={`font-semibold text-${framework.color}-300`}>{framework.modelLabel}</span> {framework.model}
                                </p>
                                <p className="text-slate-400 text-xs mb-2">{framework.process}</p>
                                <div className="bg-slate-800/60 p-2 rounded text-xs">
                                    <p className="text-slate-300"><span className="text-white font-semibold">Input:</span> {framework.input}</p>
                                    <p className="text-slate-300"><span className="text-white font-semibold">Output:</span> {framework.output}</p>
                                    <p className={`text-${framework.color}-300 font-semibold`}>{framework.performance}</p>
                                </div>
                                <p className={`text-${framework.color}-400 text-xs mt-2 font-semibold`}>Weight: {framework.weight}</p>
                            </div>
                        ))}
                    </div>
                    
                    <div className="bg-cyan-800/30 p-4 rounded-lg border border-cyan-500/50">
                        <p className="text-white font-bold text-center mb-2">{CONTENT.methods.formula.title}</p>
                        <p className="text-center font-mono text-cyan-300 text-lg">
                            {CONTENT.methods.formula.equation}
                        </p>
                        <p className="text-slate-400 text-xs text-center mt-2">
                            {CONTENT.methods.formula.note}
                        </p>
                    </div>
                </div>
                
                {/* VALIDATION SECTION */}
                <div className="bg-slate-900/60 p-6 rounded-2xl border-2 border-green-500 mb-6">
                    <div className="flex items-center mb-4">
                        <CheckCircle className="w-10 h-10 text-green-400 mr-3"/>
                        <h3 className="text-2xl font-bold text-green-300">{CONTENT.validation.title}</h3>
                    </div>
                    
                    <div className="mb-4">
                        <p className="text-slate-300 text-sm mb-3">
                            <span className="text-white font-semibold">Dataset:</span> {CONTENT.validation.dataset}
                        </p>
                    </div>
                    
                    <div className="grid grid-cols-2 gap-6 mb-4">
                        <div>
                            <p className="text-cyan-400 font-bold text-sm mb-2">{CONTENT.validation.perTargetTitle}</p>
                            <div className="space-y-2">
                                {CONTENT.validation.targets.map((t, i) => (
                                    <div key={i} className={`bg-${t.color}-900/20 p-2 rounded border border-${t.color}-500/30 flex items-center justify-between text-xs`}>
                                        <div>
                                            <p className="text-white font-bold">{t.gene}</p>
                                            <p className="text-slate-400">{t.approval}</p>
                                        </div>
                                        <div className="text-right">
                                            <p className="text-white">TLS: <span className={`text-${t.color}-300 font-bold`}>{t.score}</span></p>
                                            <p className="text-slate-400">AUROC: {t.auroc}</p>
                                        </div>
                                    </div>
                                ))}
                            </div>
                        </div>
                        
                        <div>
                            <p className="text-purple-400 font-bold text-sm mb-2">{CONTENT.validation.statsTitle}</p>
                            <div className="space-y-3">
                                {CONTENT.validation.stats.map((stat, i) => (
                                    <div key={i} className="bg-purple-900/20 p-3 rounded-lg border border-purple-500/30">
                                        <p className="text-white font-bold text-sm mb-1">{stat.title}</p>
                                        {stat.grid && (
                                            <div className="grid grid-cols-2 gap-2 text-xs">
                                                {stat.grid.map((item, j) => (
                                                    <div key={j}>
                                                        <p className="text-slate-400">{item.label}</p>
                                                        <p className={`text-${item.color}-300 font-bold text-xl`}>{item.value}</p>
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                        {stat.note && <p className="text-green-300 text-xs mt-2">{stat.note}</p>}
                                        {stat.subNote && <p className="text-slate-400 text-xs mt-2">{stat.subNote}</p>}
                                        {stat.bullets && (
                                            <div className="text-xs space-y-1">
                                                {stat.bullets.map((bullet, j) => (
                                                    <p key={j} className="text-slate-300" dangerouslySetInnerHTML={{__html: bullet}}></p>
                                                ))}
                                            </div>
                                        )}
                                    </div>
                                ))}
                            </div>
                        </div>
                    </div>
                    
                    <div className="bg-green-800/30 p-3 rounded-lg border border-green-500/50">
                        <p className="text-white font-bold text-sm mb-1">{CONTENT.validation.confound.title}</p>
                        <div className="grid grid-cols-3 gap-4 text-xs">
                            {CONTENT.validation.confound.tests.map((test, i) => (
                                <div key={i} className="text-center">
                                    <p className="text-slate-400">{test.name}</p>
                                    <p className="text-white font-mono">{test.result}</p>
                                    <p className="text-green-300 font-semibold">{test.conclusion}</p>
                                </div>
                            ))}
                        </div>
                    </div>
                </div>
                
                {/* USE CASES */}
                <div className="grid grid-cols-2 gap-6 mb-6">
                    {CONTENT.useCases.map((useCase, idx) => (
                        <div key={idx} className={`bg-slate-900/60 p-6 rounded-2xl border-2 border-${useCase.borderColor}-500`}>
                            <div className="flex items-center mb-4">
                                {useCase.icon === 'users' && <Users className="w-10 h-10 text-cyan-400 mr-3"/>}
                                {useCase.icon === 'microscope' && <Microscope className="w-10 h-10 text-purple-400 mr-3"/>}
                                <h3 className={`text-xl font-bold text-${useCase.borderColor}-300`}>{useCase.title}</h3>
                            </div>
                            
                            <div className="space-y-3 text-sm">
                                <div className={`bg-${useCase.borderColor}-900/20 p-3 rounded-lg border border-${useCase.borderColor}-500/30`}>
                                    <p className="font-bold text-white text-xs mb-2">{useCase.example.title}</p>
                                    <div className="text-xs space-y-1 text-slate-300">
                                        {useCase.example.points.map((point, i) => (
                                            <p key={i}><span className="text-white font-semibold">{point.label}</span> {point.value}</p>
                                        ))}
                                    </div>
                                </div>
                                
                                <div className="bg-slate-800/60 p-3 rounded-lg">
                                    <p className={`font-bold text-${useCase.borderColor}-400 text-xs mb-2`}>{useCase.deliverables.title}</p>
                                    <div className="text-xs space-y-1">
                                        {useCase.deliverables.items.map((item, i) => (
                                            <div key={i} className="flex items-start">
                                                <CheckCircle className="w-4 h-4 text-green-400 mr-2 mt-0.5 flex-shrink-0"/>
                                                <p className="text-slate-300"><span className="text-white font-semibold">{item.label}</span> {item.description}</p>
                                            </div>
                                        ))}
                                    </div>
                                </div>
                                
                                <div className="bg-green-900/20 p-3 rounded-lg border border-green-500/30">
                                    <p className="font-bold text-green-300 text-xs mb-1">{useCase.impact.title}</p>
                                    <div className="text-xs space-y-1 text-slate-300">
                                        {useCase.impact.points.map((point, i) => (
                                            <p key={i}>{point}</p>
                                        ))}
                                    </div>
                                </div>
                            </div>
                        </div>
                    ))}
                </div>
                
                {/* PUBLICATION SECTION */}
                <div className="bg-slate-900/60 p-6 rounded-2xl border-2 border-amber-500 mb-6">
                    <div className="flex items-center mb-4">
                        <BookOpen className="w-10 h-10 text-amber-400 mr-3"/>
                        <h3 className="text-2xl font-bold text-amber-300">{CONTENT.publication.title}</h3>
                    </div>
                    
                    <div className="grid grid-cols-3 gap-6 text-sm">
                        {CONTENT.publication.sections.map((section, i) => (
                            <div key={i} className="bg-amber-900/20 p-4 rounded-lg border border-amber-500/30">
                                <p className="font-bold text-white text-sm mb-2">{section.title}</p>
                                <div className="text-xs space-y-1 text-slate-300">
                                    {section.items.map((item, j) => (
                                        <p key={j}><span className="text-white font-semibold">{item.label}</span> {item.value}</p>
                                    ))}
                                    {section.footer && <p className={`text-${section.footer.color}-300 font-semibold mt-2`}>{section.footer.text}</p>}
                                </div>
                            </div>
                        ))}
                    </div>
                </div>
                
                {/* FOOTER */}
                <div className="bg-gradient-to-r from-slate-800 to-slate-900 p-6 rounded-2xl border-2 border-slate-600 text-center">
                    <div className="flex items-center justify-center space-x-8 text-slate-300 mb-3">
                        <div className="text-left">
                            <p className="font-bold text-white text-lg">{CONTENT.footer.name}</p>
                            <p className="text-sm">{CONTENT.footer.title}</p>
                        </div>
                        <div className="text-left text-sm">
                            {CONTENT.footer.contact.map((item, i) => (
                                <p key={i}>{item}</p>
                            ))}
                        </div>
                    </div>
                    <p className="text-slate-500 text-xs">
                        {CONTENT.footer.disclaimer}
                    </p>
                </div>
                
            </div>
        </div>
    );
};

export default CrisPROTechnicalOnePager;

// ============================================================================
// CONTENT CONSTANTS
// ============================================================================

const CONTENT = {
    header: {
        title: "CrisPRO.ai",
        subtitle: "Genomic AI for Patient Stratification in Cancer Clinical Trials",
        description: "Multi-modal foundation models predict treatment response from whole-genome data (AUROC 0.976)",
        badge: "Research Use Only (RUO) • Manuscript in Preparation (Q1 2026)"  // ✅ FIXED
    },
    
    spe: {
        title: "THE S/P/E DOCTRINE: FROM CHAOS TO CERTAINTY",
        problem: {
            title: "THE PROBLEM",
            points: [
                "50+ therapy options (which one?)",
                "Black box predictions (why?)",
                "No audit trail (how confident?)"
            ],
            conclusion: "Doctor shrugs"
        },
        engine: [
            {
                label: "S: Sequence",
                description: "Evo2 disruption scores",
                weight: "35%",
                color: "cyan",
                icon: "dna"
            },
            {
                label: "P: Pathway",
                description: "MoA pathway burden",
                weight: "35%",
                color: "green",
                icon: "branch"
            },
            {
                label: "E: Evidence",
                description: "ClinVar + literature",
                weight: "30%",
                color: "purple",
                icon: "book"
            }
        ],
        efficacyScore: {
            title: "Efficacy Score",
            description: "Weighted composite (0-1)"
        },
        deliverable: {
            title: "THE DELIVERABLE",
            results: [
                { label: "1. BRAF V600E", score: "0.89", color: "text-green-300" },
                { label: "2. MET amp", score: "0.76", color: "text-cyan-300" },
                { label: "3. VEGFA", score: "0.72", color: "text-purple-300" },
                { label: "4-12 ranked...", score: "...", color: "text-slate-500", opacity: true }
            ],
            footer: {
                provenance: "✓ Auditable provenance",
                runId: "Run ID: abc123"
            }
        },
        summary: {
            main: "From 50+ Options → 5-12 Ranked Therapies in 6 Weeks",
            details: "VUS reduction: 40% → 15% | Cost savings: $2.1M per program"
        }
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
                "97.6% AUROC predicting FDA-approved target relevance",
                "100% accuracy predicting Multple Myeloma therapy response (n=12 patients)",
                "Validated on 7 FDA-approved cancer drug targets (computational, RUO)",  
                "Cohen's d > 2.5 effect size (exceptional discrimination)"
            ]
        }
    },
    
    methods: {
        title: "METHODS: MULTI-MODAL SCORING FRAMEWORK",
        frameworks: [
            {
                name: "Functionality (F)",
                modelLabel: "Model:",
                model: "Evo2-7B (9.3T tokens)",
                process: "DNA sequence → pathogenicity score (Δ fitness)",
                input: "8,192 bp context",
                output: "0-1 damage score",
                performance: "AUROC: 0.85",
                weight: "35%",
                color: "cyan",
                icon: "dna"
            },
            {
                name: "Essentiality (E)",
                modelLabel: "Data:",
                model: "DepMap Chronos",
                process: "CRISPR knockout → cell viability (n=23,000 genes)",
                input: "Gene dependency",
                output: "Percentile rank",
                performance: "Effect: d = 2.5-3.5",
                weight: "35%",
                color: "green",
                icon: "target"
            },
            {
                name: "Chromatin (C)",
                modelLabel: "Model:",
                model: "Enformer",
                process: "200kb sequence → accessibility (DNase-seq)",
                input: "393,216 bp window",
                output: "Accessibility score",
                performance: "Corr: ρ = 0.75",
                weight: "15%",
                color: "purple",
                icon: "activity"
            },
            {
                name: "Regulatory (R)",
                modelLabel: "Model:",
                model: "SpliceAI",
                process: "Splice site disruption (32kb context)",
                input: "Variant position",
                output: "Δ splice score",
                performance: "AUROC: 0.95",
                weight: "15%",
                color: "orange",
                icon: "branch"
            }
        ],
        formula: {
            title: "Target Lock Score Formula",
            equation: "TLS = 0.35 × F + 0.35 × E + 0.15 × C + 0.15 × R",
            note: "Weights: 0.35/0.35/0.15/0.15 (heuristic, validated on 56 gene-step combinations)"  // ✅ FIXED (no false CV claim)
        }
    },
    
    validation: {
        title: "VALIDATION: 7 FDA-APPROVED CANCER DRUG TARGETS",
        dataset: "7 FDA-approved metastatic cancer targets (BRAF, CXCR4, MET, VEGFA, MMP2, SNAIL1, TWIST1) × 8 metastatic steps (primary growth → colonization) = 56 predictions. Ground truth: FDA approval status + clinical trial outcomes.",
        perTargetTitle: "Per-Target Performance",
        targets: [
            { gene: 'MET', score: '0.467', auroc: '0.994', approval: 'Capmatinib (2020)', color: 'green' },
            { gene: 'CXCR4', score: '0.463', auroc: '1.000', approval: 'Plerixafor (2008)', color: 'green' },
            { gene: 'VEGFA', score: '0.447', auroc: '0.979', approval: 'Bevacizumab (2004)', color: 'green' },
            { gene: 'BRAF', score: '0.445', auroc: '1.000', approval: 'Vemurafenib (2011)', color: 'green' },
            { gene: 'MMP2', score: '0.402', auroc: '0.953', approval: 'Phase 3 trials', color: 'yellow' },
            { gene: 'SNAIL1', score: '0.378', auroc: '0.898', approval: 'Preclinical', color: 'orange' },
            { gene: 'TWIST1', score: '0.354', auroc: '0.981', approval: 'Preclinical', color: 'orange' }
        ],
        statsTitle: "Statistical Performance",
        stats: [
            {
                title: "Overall Discrimination",
                grid: [
                    { label: "Mean AUROC:", value: "0.976", color: "purple" },
                    { label: "95% CI:", value: "[0.952, 1.00]", color: "purple" }
                ],
                note: "✓ 5 of 8 steps achieve AUROC = 1.000 (perfect)"
            },
            {
                title: "Effect Size Analysis",
                grid: [
                    { label: "Cohen's d:", value: "2.5-3.5", color: "purple" },
                    { label: "Interpretation:", value: "Exceptional", color: "green" }
                ],
                subNote: "Typical genomics: d=0.3-0.5 (small-medium)"
            },
            {
                title: "Clinical Utility",
                bullets: [
                    '• Precision@3: <span class="text-green-300 font-semibold">100%</span> (all top-3 relevant)',
                    '• Precision@5: <span class="text-cyan-300 font-semibold">93%</span> (80-100% per step)',
                    '• Enrichment: <span class="text-purple-300 font-semibold">30% top stratum</span> = 2x response'
                ]
            }
        ],
        confound: {
            title: "Confound Analysis (Negative Controls)",
            tests: [
                { name: "Gene Length", result: "ρ = -0.003, p = 0.985", conclusion: "No correlation ✓" },
                { name: "GC Content", result: "ρ = -0.111, p = 0.507", conclusion: "No correlation ✓" },
                { name: "Exon Count", result: "ρ = -0.090, p = 0.593", conclusion: "No correlation ✓" }
            ]
        }
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
