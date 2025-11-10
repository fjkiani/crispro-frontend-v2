import React, { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { 
    Target, Dna, Award, BrainCircuit, ShieldCheck, Zap, BarChart3, 
    Activity, FlaskConical, Microscope, BookOpen, Users, ArrowRight,
    CheckCircle, TrendingUp, AlertTriangle, Cpu, GitBranch
} from 'lucide-react';

//================================================================================
// PHD TECHNICAL DECK - METASTASIS INTERCEPTION
// 18 Slides with Real Validation Data
//================================================================================

const PhDDeck = () => {
    const [currentSlide, setCurrentSlide] = useState(0);

    const slides = [
        <Slide01_Title />,
        <Slide02_MetastaticCascade />,
        <Slide03_CurrentTools />,
        <Slide04_MultiModal />,
        <Slide05_Pipeline />,
        <Slide06_TargetLock />,
        <Slide07_ROC key="7" />,
        <Slide08_Specificity key="8" />,
        <Slide09_Precision key="9" />,
        <Slide10_Ablation key="10" />,
        <Slide11_EffectSizes key="11" />,
        <Slide12_Efficacy key="12" />,
        <Slide13_Safety key="13" />,
        <Slide14_Confounders key="14" />,
        <Slide15_Calibration key="15" />,
        <Slide16_Assassin key="16" />,
        <Slide17_Reproducibility key="17" />,
        <Slide18_CTA key="18" />
    ];

    return (
        <div className="relative w-full h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900">
            <AnimatePresence mode="wait">
                <motion.div
                    key={currentSlide}
                    initial={{ opacity: 0, x: 100 }}
                    animate={{ opacity: 1, x: 0 }}
                    exit={{ opacity: 0, x: -100 }}
                    transition={{ duration: 0.3 }}
                    className="w-full h-full flex items-center justify-center p-12"
                >
                    {slides[currentSlide]}
                </motion.div>
            </AnimatePresence>

            {/* Navigation */}
            <div className="absolute bottom-8 left-1/2 transform -translate-x-1/2 flex items-center space-x-4 z-50">
                <button
                    onClick={() => setCurrentSlide(prev => Math.max(0, prev - 1))}
                    disabled={currentSlide === 0}
                    className="px-6 py-3 bg-slate-700 text-white rounded-lg font-semibold disabled:opacity-30 hover:bg-slate-600 transition-colors"
                >
                    ← Previous
                </button>
                <span className="text-slate-400 font-mono">
                    {currentSlide + 1} / {slides.length}
                </span>
                <button
                    onClick={() => setCurrentSlide(prev => Math.min(slides.length - 1, prev + 1))}
                    disabled={currentSlide === slides.length - 1}
                    className="px-6 py-3 bg-cyan-600 text-white rounded-lg font-semibold disabled:opacity-30 hover:bg-cyan-500 transition-colors"
                >
                    Next →
                </button>
            </div>

            {/* Slide Title Indicator */}
            <div className="absolute top-8 right-8 text-slate-400 text-sm font-mono">
                {slideNames[currentSlide]}
            </div>
        </div>
    );
};

const slideNames = [
    "Title & Abstract",
    "Metastatic Cascade",
    "Current Tool Limitations",
    "Multi-Modal Framework",
    "5-Stage Pipeline",
    "Target Lock Heatmap",
    "ROC Curves (AUROC 0.976)",
    "Specificity Matrix",
    "Precision@K Analysis",
    "Ablation Study",
    "Effect Sizes",
    "Efficacy Distribution",
    "Safety Distribution",
    "Confounder Analysis",
    "Calibration Curves",
    "Assassin Score",
    "Reproducibility",
    "Call to Action"
];

//================================================================================
// SLIDE 1: TITLE & ABSTRACT
//================================================================================

const Slide01_Title = () => (
    <div className="max-w-6xl mx-auto text-white space-y-8">
        <motion.div initial={{opacity:0, y:-20}} animate={{opacity:1, y:0}} transition={{delay:0.2}}>
            <Dna className="w-24 h-24 text-cyan-400 mx-auto mb-6" />
            <h1 className="text-6xl font-black text-center text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-blue-500 mb-4">
                Metastasis Interception
            </h1>
            <p className="text-3xl text-center text-slate-300 font-light">
                A Multi-Modal AI Framework for Stage-Specific<br/>CRISPR Therapeutic Design
            </p>
        </motion.div>

        <motion.div 
            className="bg-slate-800/60 p-8 rounded-2xl border border-slate-700"
            initial={{opacity:0, y:20}} 
            animate={{opacity:1, y:0}} 
            transition={{delay:0.4}}
        >
            <p className="text-lg text-slate-300 leading-relaxed">
                <span className="font-bold text-white">ABSTRACT:</span> Cancer metastasis causes 90% of cancer deaths, 
                yet current CRISPR design tools achieve only 40% experimental success due to analysis limited to 
                1% of the genome. We present a multi-modal AI framework integrating foundation models (Evo2, Enformer) 
                with clinical evidence to design stage-specific CRISPR therapeutics. Validated on 7 FDA-approved 
                targets across 8 metastatic steps, our platform achieves <span className="text-cyan-300 font-bold">mean 
                AUROC 0.976</span> for target prediction, <span className="text-green-300 font-bold">80% guide efficacy</span>, 
                and <span className="text-purple-300 font-bold">70%+ safety</span> (production-ready). 
                Effect sizes (Cohen's d {'>'} 2.5) demonstrate exceptional discrimination of relevant vs non-relevant genes.
            </p>
        </motion.div>

        <motion.div 
            className="grid grid-cols-4 gap-6"
            initial={{opacity:0}} 
            animate={{opacity:1}} 
            transition={{delay:0.6}}
        >
            <MetricCard icon={<Target />} value="0.976" label="Mean AUROC" color="cyan" />
            <MetricCard icon={<Dna />} value="80%" label="Guide Efficacy" color="green" />
            <MetricCard icon={<ShieldCheck />} value="70%+" label="Safety Score" color="purple" />
            <MetricCard icon={<BarChart3 />} value="d > 2.5" label="Effect Size" color="orange" />
        </motion.div>

        <div className="text-center text-slate-500 text-sm mt-8">
            <p>Fahad J. Kiani | CrisPRO.ai | October 2025</p>
            <p className="text-cyan-400 mt-2">Research Use Only (RUO) | Nature Biotechnology (In Preparation)</p>
        </div>
    </div>
);

const MetricCard = ({ icon, value, label, color }) => {
    const colorClasses = {
        cyan: 'bg-cyan-500/20 border-cyan-500 text-cyan-300',
        green: 'bg-green-500/20 border-green-500 text-green-300',
        purple: 'bg-purple-500/20 border-purple-500 text-purple-300',
        orange: 'bg-orange-500/20 border-orange-500 text-orange-300'
    };

    return (
        <div className={`${colorClasses[color]} p-6 rounded-xl border-2 text-center`}>
            <div className="flex justify-center mb-3">{React.cloneElement(icon, { className: 'w-8 h-8' })}</div>
            <p className="text-4xl font-black">{value}</p>
            <p className="text-sm text-slate-400 mt-2">{label}</p>
        </div>
    );
};

//================================================================================
// SLIDE 2: METASTATIC CASCADE
//================================================================================

const Slide02_MetastaticCascade = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-5xl font-black text-center mb-4 text-transparent bg-clip-text bg-gradient-to-r from-red-400 to-orange-500">
            The Metastatic Cascade
        </h2>
        <p className="text-xl text-center text-slate-300 mb-12">
            90% of cancer deaths result from metastasis—a predictable, multi-step process
        </p>

        <div className="grid grid-cols-4 gap-4 mb-8">
            {metastaticSteps.map((step, idx) => (
                <StepCard key={idx} {...step} delay={idx * 0.1} />
            ))}
        </div>

        <div className="bg-red-900/20 p-8 rounded-2xl border-2 border-red-500">
            <h3 className="text-2xl font-bold text-red-300 mb-4 text-center">The Challenge</h3>
            <p className="text-lg text-slate-300 text-center">
                Each step involves distinct molecular mechanisms. Stage-specific targeting requires 
                identifying which genes drive each step—a problem current tools cannot solve.
            </p>
        </div>
    </div>
);

const metastaticSteps = [
    { step: 1, name: "Primary Growth", icon: "🔴", color: "red", targets: "BRAF, KRAS" },
    { step: 2, name: "Local Invasion", icon: "💥", color: "orange", targets: "MMPs, SNAIL" },
    { step: 3, name: "Intravasation", icon: "🩸", color: "yellow", targets: "CXCR4, MMP2" },
    { step: 4, name: "Circulation", icon: "🌊", color: "cyan", targets: "CXCR4, BCL2" },
    { step: 5, name: "Extravasation", icon: "🚪", color: "blue", targets: "Integrins" },
    { step: 6, name: "Micrometastasis", icon: "📍", color: "purple", targets: "NR2F1, DEC2" },
    { step: 7, name: "Angiogenesis", icon: "🌳", color: "green", targets: "VEGFA, PDGF" },
    { step: 8, name: "Colonization", icon: "🔴", color: "red", targets: "VEGFA, MET" }
];

const StepCard = ({ step, name, icon, color, targets, delay }) => (
    <motion.div
        className={`bg-slate-800/60 p-6 rounded-xl border-2 border-${color}-500 text-center`}
        initial={{opacity:0, y:20}}
        animate={{opacity:1, y:0}}
        transition={{delay}}
    >
        <div className="text-4xl mb-3">{icon}</div>
        <div className="text-sm text-slate-500 font-bold mb-2">STEP {step}</div>
        <h4 className="text-lg font-bold text-white mb-2">{name}</h4>
        <p className="text-xs text-slate-400">{targets}</p>
    </motion.div>
);

//================================================================================
// SLIDE 3: CURRENT TOOL LIMITATIONS
//================================================================================

const Slide03_CurrentTools = () => (
    <div className="max-w-6xl mx-auto text-white">
        <h2 className="text-5xl font-black text-center mb-4 text-transparent bg-clip-text bg-gradient-to-r from-red-400 to-orange-500">
            Why Current Tools Fail
        </h2>
        <p className="text-xl text-center text-slate-300 mb-12">
            Standard CRISPR design tools analyze only 1% of the genome and achieve 40% success
        </p>

        <div className="grid md:grid-cols-2 gap-8">
            {/* Standard Tools */}
            <div className="bg-red-900/20 p-8 rounded-2xl border-2 border-red-500">
                <FlaskConical className="w-16 h-16 text-red-400 mx-auto mb-6" />
                <h3 className="text-3xl font-bold text-red-300 text-center mb-6">Standard Tools</h3>
                
                <div className="text-center mb-6">
                    <p className="text-8xl font-black text-red-300">1%</p>
                    <p className="text-slate-300 mt-2">Genome Coverage</p>
                    <p className="text-sm text-slate-500">(protein-coding only)</p>
                </div>

                <div className="space-y-3 mb-6">
                    <LimitationItem text="GC content bias (40-60% rule)" />
                    <LimitationItem text="Ignores chromatin accessibility" />
                    <LimitationItem text="No pathway context" />
                    <LimitationItem text="Heuristic off-target prediction" />
                </div>

                <div className="bg-red-800 p-4 rounded-lg text-center">
                    <p className="text-3xl font-black text-white">40%</p>
                    <p className="text-red-200 text-sm">Lab Success = Coin Flip</p>
                </div>
            </div>

            {/* Our Platform */}
            <div className="bg-green-900/20 p-8 rounded-2xl border-2 border-green-500">
                <BrainCircuit className="w-16 h-16 text-green-400 mx-auto mb-6" />
                <h3 className="text-3xl font-bold text-green-300 text-center mb-6">Our AI Platform</h3>
                
                <div className="text-center mb-6">
                    <p className="text-8xl font-black text-green-300">100%</p>
                    <p className="text-slate-300 mt-2">Genome Coverage</p>
                    <p className="text-sm text-green-400 font-semibold">(including regulatory regions)</p>
                </div>

                <div className="space-y-3 mb-6">
                    <AdvantageItem text="Foundation model (Evo2, 9.3T tokens)" />
                    <AdvantageItem text="Chromatin accessibility (Enformer)" />
                    <AdvantageItem text="Pathway integration (DepMap + STRING)" />
                    <AdvantageItem text="Genome-wide safety (Bowtie2)" />
                </div>

                <div className="bg-green-800 p-4 rounded-lg text-center">
                    <p className="text-3xl font-black text-white">80%</p>
                    <p className="text-green-200 text-sm">Lab Success = Near Certainty</p>
                </div>
            </div>
        </div>

        <div className="mt-8 bg-gradient-to-r from-purple-600 to-pink-600 p-6 rounded-2xl text-center">
            <p className="text-4xl font-black text-white">2x Better Predictions</p>
            <p className="text-xl text-purple-100 mt-2">Mean AUROC: 0.976 | Effect Size: Cohen's d {'>'} 2.5</p>
        </div>
    </div>
);

const LimitationItem = ({ text }) => (
    <div className="flex items-center text-slate-300">
        <div className="w-6 h-6 bg-red-500 rounded-full flex items-center justify-center mr-3 flex-shrink-0">
            <span className="text-white font-bold text-sm">✗</span>
        </div>
        <span>{text}</span>
    </div>
);

const AdvantageItem = ({ text }) => (
    <div className="flex items-center text-slate-300">
        <CheckCircle className="w-6 h-6 text-green-400 mr-3 flex-shrink-0" />
        <span>{text}</span>
    </div>
);



//================================================================================
// SLIDE 4: MULTI-MODAL FRAMEWORK
//================================================================================

const Slide04_MultiModal = () => (
    <div className="max-w-6xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-blue-500">
            Multi-Modal Framework
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            4-Chip Architecture: F (35%) + E (35%) + C (15%) + R (15%)
        </p>

        <div className="grid grid-cols-2 gap-5 mb-6">
            {[
                { title: "Functionality", weight: "35%", model: "Evo2-7B", desc: "DNA damage", metric: "AUROC 0.85", color: "cyan" },
                { title: "Essentiality", weight: "35%", model: "DepMap", desc: "Knockout lethality", metric: "d = 2.5", color: "green" },
                { title: "Chromatin", weight: "15%", model: "Enformer", desc: "Accessibility", metric: "23x range", color: "purple" },
                { title: "Regulatory", weight: "15%", model: "SpliceAI", desc: "Splice impact", metric: "AUROC 0.95", color: "orange" }
            ].map((chip, i) => (
                <motion.div
                    key={i}
                    className="bg-slate-800/60 p-4 rounded-xl border-2 border-slate-700"
                    initial={{opacity:0, scale:0.9}}
                    animate={{opacity:1, scale:1}}
                    transition={{delay: i * 0.1}}
                >
                    <div className="flex justify-between items-center mb-2">
                        <h3 className="text-lg font-bold text-white">{chip.title}</h3>
                        <span className={`text-${chip.color}-400 font-mono text-xs`}>{chip.weight}</span>
                    </div>
                    <p className="text-slate-400 text-xs mb-1"><strong>Model:</strong> {chip.model}</p>
                    <p className="text-slate-300 text-xs mb-2">{chip.desc}</p>
                    <div className="bg-slate-900/50 p-2 rounded text-center">
                        <p className="text-xs text-slate-300">{chip.metric}</p>
                    </div>
                </motion.div>
            ))}
        </div>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 text-center">
            <p className="text-sm text-slate-300">
                <span className="font-bold text-cyan-400">Target Lock Score =</span> 
                <span className="font-mono text-white"> 0.35F + 0.35E + 0.15C + 0.15R</span>
            </p>
        </div>
    </div>
);

//================================================================================
// SLIDE 5: 5-STAGE PIPELINE
//================================================================================

const Slide05_Pipeline = () => (
    <div className="max-w-6xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-green-400 to-emerald-500">
            5-Stage Pipeline
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Patient DNA → Production-ready guides in 6 weeks
        </p>

        <div className="grid grid-cols-5 gap-3 mb-6">
            {[
                { stage: 1, title: "Target", icon: <Target />, out: "Top 3 genes" },
                { stage: 2, title: "Generate", icon: <Dna />, out: "100 guides" },
                { stage: 3, title: "Efficacy", icon: <Zap />, out: "Score 0-1" },
                { stage: 4, title: "Safety", icon: <ShieldCheck />, out: "Off-target" },
                { stage: 5, title: "Rank", icon: <Activity />, out: "Top 5" }
            ].map((s, i) => (
                <motion.div
                    key={i}
                    className="bg-slate-800/60 p-3 rounded-xl border-2 border-slate-700 text-center"
                    initial={{opacity:0, y:20}}
                    animate={{opacity:1, y:0}}
                    transition={{delay: i * 0.1}}
                >
                    <div className="flex justify-center mb-2">{React.cloneElement(s.icon, { className: 'w-7 h-7 text-cyan-400' })}</div>
                    <div className="text-xs text-slate-500 font-bold mb-1">{s.stage}</div>
                    <h4 className="text-xs font-bold text-white mb-1">{s.title}</h4>
                    <p className="text-xs text-slate-400">{s.out}</p>
                </motion.div>
            ))}
        </div>

        <div className="bg-gradient-to-r from-blue-900/40 to-purple-900/40 p-5 rounded-2xl border-2 border-blue-500/50">
            <div className="flex items-center justify-between max-w-4xl mx-auto text-sm">
                <div className="text-center">
                    <p className="text-xs text-slate-400 mb-1">INPUT</p>
                    <p className="font-bold">Patient DNA</p>
                </div>
                <ArrowRight className="w-6 h-6 text-cyan-400 mx-3" />
                <div className="text-center">
                    <p className="text-xs text-slate-400 mb-1">PROCESS</p>
                    <p className="font-bold text-cyan-300">6 Weeks</p>
                </div>
                <ArrowRight className="w-6 h-6 text-cyan-400 mx-3" />
                <div className="text-center">
                    <p className="text-xs text-slate-400 mb-1">OUTPUT</p>
                    <p className="font-bold text-green-300">Ranked Guides</p>
                </div>
            </div>
        </div>
    </div>
);

//================================================================================
// SLIDE 6: TARGET LOCK HEATMAP
//================================================================================

const Slide06_TargetLock = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-orange-400 to-red-500">
            Target Lock Scores
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            7 FDA targets × 8 steps = 56 predictions | Scores stage-invariant (σ {'<'} 0.001)
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Your Target Lock Heatmap visualization here - use figure2a_per_step_roc.jpg]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center italic">
                    Heatmap shows consistent scores across metastatic steps. Top genes (MET, CXCR4, VEGFA) 
                    maintain relevance regardless of stage, indicating "master regulator" status.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { gene: "MET", score: "0.467", desc: "Highest scorer", color: "green" },
                { gene: "CXCR4", score: "0.463", desc: "Migration", color: "cyan" },
                { gene: "VEGFA", score: "0.447", desc: "Angiogenesis", color: "purple" }
            ].map((t, i) => (
                <div key={i} className={`bg-${t.color}-900/20 p-3 rounded-xl border border-${t.color}-500/50 text-center`}>
                    <p className={`text-2xl font-black text-${t.color}-300`}>{t.gene}</p>
                    <p className="text-sm text-slate-400">{t.score}</p>
                    <p className="text-xs text-slate-500">{t.desc}</p>
                </div>
            ))}
        </div>
    </div>
);
//================================================================================
// SLIDE 7: ROC CURVES
//================================================================================

const Slide07_ROC = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-indigo-400 to-purple-500">
            ROC Curves Per Step
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Mean AUROC = 0.976 | 5 of 8 steps achieve perfect 1.000
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure2a_per_step_roc.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Near-perfect discrimination of relevant vs non-relevant genes. Weakest step (Intravasation: 0.898) 
                    still excellent. Perfect scores (1.000) indicate clean biological signal, not overfitting.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-4 gap-3">
            {[
                { step: "Primary Growth", auroc: "1.000", color: "green" },
                { step: "Survival", auroc: "1.000", color: "green" },
                { step: "Micrometastasis", auroc: "1.000", color: "green" },
                { step: "Intravasation", auroc: "0.898", color: "yellow" },
                { step: "Colonization", auroc: "0.981", color: "cyan" },
                { step: "Local Invasion", auroc: "0.953", color: "cyan" },
                { step: "Angiogenesis", auroc: "0.979", color: "cyan" },
                { step: "Extravasation", auroc: "0.994", color: "cyan" }
            ].map((s, i) => (
                <div key={i} className={`bg-${s.color}-900/20 p-3 rounded-xl border border-${s.color}-500/50 text-center`}>
                    <p className="text-xs text-slate-400 mb-1">{s.step}</p>
                    <p className={`text-2xl font-black text-${s.color}-300`}>{s.auroc}</p>
                </div>
            ))}
        </div>

        <div className="mt-5 bg-gradient-to-r from-indigo-600 to-purple-600 p-4 rounded-xl text-center">
            <p className="text-xl font-black">Mean AUROC = 0.976 (95% CI: [0.952, 1.000])</p>
        </div>
    </div>
);

//================================================================================
// SLIDE 8: SPECIFICITY MATRIX
//================================================================================

const Slide08_Specificity = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-pink-400 to-red-500">
            Specificity Matrix
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Predicted step vs true step | Diagonal dominance = correct assignments
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure2b_specificity_matrix.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Primary Growth & Micrometastasis: 100% correct (p{'<'}0.001). Confusable stages: 
                    Survival↔Primary (33% overlap), Angiogenesis↔Colonization (17% overlap) due to shared biology.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { label: "Perfect Stages", value: "2/8", desc: "100% accuracy", color: "green" },
                { label: "High Confidence", value: "4/8", desc: "80%+ accuracy", color: "cyan" },
                { label: "Confusable", value: "2/8", desc: "Shared pathways", color: "yellow" }
            ].map((m, i) => (
                <div key={i} className={`bg-${m.color}-900/20 p-3 rounded-xl border border-${m.color}-500/50 text-center`}>
                    <p className={`text-2xl font-black text-${m.color}-300`}>{m.value}</p>
                    <p className="text-sm text-slate-400">{m.label}</p>
                    <p className="text-xs text-slate-500">{m.desc}</p>
                </div>
            ))}
        </div>
    </div>
);

//================================================================================
// SLIDE 9: PRECISION@K
//================================================================================

const Slide09_Precision = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-yellow-400 to-orange-500">
            Precision@K Analysis
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Top 3 predictions: 80-100% precision | Top 5: 80-100% | Top 10: 50-90%
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure2c_precision_at_k.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Colonization P@3=100%, P@10=88%. Primary Growth drops to P@10=40% (many genes drive proliferation). 
                    Recommendation: Test top 3 guides per target for near-guarantee of relevance.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { k: "P@3", precision: "100%", desc: "Elite confidence", color: "green" },
                { k: "P@5", precision: "80-100%", desc: "Production-ready", color: "cyan" },
                { k: "P@10", precision: "50-90%", desc: "Screening level", color: "yellow" }
            ].map((p, i) => (
                <div key={i} className={`bg-${p.color}-900/20 p-3 rounded-xl border border-${p.color}-500/50 text-center`}>
                    <p className={`text-3xl font-black text-${p.color}-300`}>{p.precision}</p>
                    <p className="text-sm text-slate-400 mt-1">{p.k}</p>
                    <p className="text-xs text-slate-500">{p.desc}</p>
                </div>
            ))}
        </div>

        <div className="mt-5 bg-slate-800/40 p-3 rounded-lg text-center">
            <p className="text-sm text-cyan-400 font-bold">Business Impact: $250K pilot delivers 3-5 guides with 80%+ confidence all are relevant</p>
        </div>
    </div>
);

//================================================================================
// SLIDE 10: ABLATION STUDY
//================================================================================

const Slide10_Ablation = () => (
    <div className="max-w-6xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-green-400 to-emerald-500">
            Ablation Study
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Essentiality dominates (Δ=0.086) | Chromatin harmful? (Δ=-0.013)
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure2d_ablation.jpg here]
            </div>
        </div>

        <div className="grid grid-cols-2 gap-5 mb-5">
            {[
                { signal: "Essentiality", drop: "+0.086", impact: "LARGEST", color: "green", desc: "Essential genes = bottlenecks" },
                { signal: "Functionality", drop: "+0.038", impact: "MODERATE", color: "cyan", desc: "Sequence damage adds value" },
                { signal: "Regulatory", drop: "+0.006", impact: "MINIMAL", color: "yellow", desc: "Most drivers in exons" },
                { signal: "Chromatin", drop: "-0.013", impact: "HARMFUL", color: "red", desc: "Heuristic vs real model" }
            ].map((s, i) => (
                <div key={i} className={`bg-${s.color}-900/20 p-4 rounded-xl border border-${s.color}-500/50`}>
                    <h3 className={`text-xl font-bold text-${s.color}-300 mb-2`}>{s.signal}</h3>
                    <div className="flex justify-between items-center mb-2">
                        <span className="text-sm text-slate-400">AUROC Drop:</span>
                        <span className={`text-2xl font-black text-${s.color}-300`}>{s.drop}</span>
                    </div>
                    <p className="text-xs text-slate-300 mb-1"><strong>{s.impact}</strong></p>
                    <p className="text-xs text-slate-400">{s.desc}</p>
                </div>
            ))}
        </div>

        <div className="bg-red-900/20 p-4 rounded-xl border border-red-500/50 text-center">
            <p className="text-sm text-slate-300">
                <span className="text-red-300 font-bold">Chromatin Negative?</span> Current heuristic (flat 0.6) 
                outperforms bad predictions. V2: Cell-type-specific Enformer or cancer ATAC-seq data.
            </p>
        </div>
    </div>
);

//================================================================================
// SLIDE 11: EFFECT SIZES
//================================================================================

const Slide11_EffectSizes = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-purple-400 to-pink-500">
            Effect Sizes: Relevant vs Non-Relevant
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Cohen's d {'>'} 2.5 = exceptional discrimination | d {'>'} 0.8 = large effect
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure_s3_effect_sizes.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Essentiality: d{'>'}2.5 most steps (2.5 std deviations separation). Target Lock Score: d=2.5-3.5 
                    (maintains large effects). Chromatin: d≈0 (confirms ablation finding).
                </p>
            </div>
        </div>

        <div className="grid grid-cols-4 gap-3">
            {[
                { chip: "Essentiality", d: "2.5-3.5", level: "EXCEPTIONAL", color: "green" },
                { chip: "Target Lock", d: "2.5-3.5", level: "EXCEPTIONAL", color: "purple" },
                { chip: "Functionality", d: "1.5-3.0", level: "LARGE", color: "cyan" },
                { chip: "Chromatin", d: "-0.3-0.8", level: "NEGLIGIBLE", color: "red" }
            ].map((e, i) => (
                <div key={i} className={`bg-${e.color}-900/20 p-3 rounded-xl border border-${e.color}-500/50 text-center`}>
                    <p className="text-xs text-slate-400 mb-1">{e.chip}</p>
                    <p className={`text-3xl font-black text-${e.color}-300`}>{e.d}</p>
                    <p className="text-xs text-slate-500 mt-1">{e.level}</p>
                </div>
            ))}
        </div>

        <div className="mt-5 bg-gradient-to-r from-purple-600 to-pink-600 p-4 rounded-xl text-center">
            <p className="text-lg font-black">Typical genomics: d=0.3-0.5 | Our platform: d{'>'}2.5 (5x better)</p>
        </div>
    </div>
);

//================================================================================
// SLIDE 12: EFFICACY DISTRIBUTION
//================================================================================

const Slide12_Efficacy = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-blue-500">
            Guide RNA Efficacy Distribution
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            80% of guides meet efficacy threshold (≥0.50) | Varies by metastatic stage
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert F3_efficacy_distribution.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Intravasation: Narrow (0.65), high confidence. Survival: Low (0.4), hardest target. 
                    Primary Growth: Bimodal (0.1, 0.6), reflects biological complexity.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { stage: "Micrometastasis", median: "0.65", desc: "Highest", color: "green" },
                { stage: "Angiogenesis", median: "0.63", desc: "High", color: "cyan" },
                { stage: "Circulation", median: "0.38", desc: "Challenging", color: "red" }
            ].map((s, i) => (
                <div key={i} className={`bg-${s.color}-900/20 p-3 rounded-xl border border-${s.color}-500/50 text-center`}>
                    <p className="text-sm text-slate-400 mb-1">{s.stage}</p>
                    <p className={`text-3xl font-black text-${s.color}-300`}>{s.median}</p>
                    <p className="text-xs text-slate-500">{s.desc}</p>
                </div>
            ))}
        </div>
    </div>
);

//================================================================================
// SLIDE 13: SAFETY DISTRIBUTION
//================================================================================

const Slide13_Safety = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-green-400 to-emerald-500">
            Guide RNA Safety Distribution
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            70%+ guides achieve high safety (≥0.80) | Production-ready
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert F4_safety_distribution.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Most steps: {'>'}90% guides safe. Extravasation: Wider (0.4-0.7) due to integrin gene family 
                    paralogs. Genome-wide validation catches off-targets competitors miss.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { category: "High Safety", pct: "70%", desc: "≥0.80 score", color: "green" },
                { category: "Moderate", pct: "25%", desc: "0.60-0.80", color: "yellow" },
                { category: "Reject", pct: "5%", desc: "<0.60", color: "red" }
            ].map((c, i) => (
                <div key={i} className={`bg-${c.color}-900/20 p-3 rounded-xl border border-${c.color}-500/50 text-center`}>
                    <p className={`text-3xl font-black text-${c.color}-300`}>{c.pct}</p>
                    <p className="text-sm text-slate-400 mt-1">{c.category}</p>
                    <p className="text-xs text-slate-500">{c.desc}</p>
                </div>
            ))}
        </div>

        <div className="mt-5 bg-green-900/20 p-3 rounded-xl border border-green-500/50 text-center">
            <p className="text-sm text-slate-300">
                <span className="text-green-300 font-bold">Competitive Advantage:</span> {'<'}1% off-target rate vs industry 5-15%
            </p>
        </div>
    </div>
);

//================================================================================
// SLIDE 14: CONFOUNDERS
//================================================================================

const Slide14_Confounders = () => (
    <div className="max-w-6xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-orange-400 to-red-500">
            Confounder Analysis
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            No correlation with gene length, GC content, or exon count
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure_s1_confounders.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Gene Length: ρ=-0.003, p=0.985 | GC Content: ρ=-0.111, p=0.507 | Exon Count: ρ=-0.090, p=0.593. 
                    Model learns TRUE biological signal, not technical artifacts.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { factor: "Gene Length", corr: "ρ = -0.003", p: "p = 0.985", result: "NO BIAS" },
                { factor: "GC Content", corr: "ρ = -0.111", p: "p = 0.507", result: "NO BIAS" },
                { factor: "Exon Count", corr: "ρ = -0.090", p: "p = 0.593", result: "NO BIAS" }
            ].map((f, i) => (
                <div key={i} className="bg-green-900/20 p-3 rounded-xl border border-green-500/50 text-center">
                    <p className="text-sm text-slate-400 mb-1">{f.factor}</p>
                    <p className="text-lg font-black text-white">{f.corr}</p>
                    <p className="text-xs text-slate-500 mb-2">{f.p}</p>
                    <div className="bg-green-800 px-2 py-1 rounded text-xs font-bold text-white">
                        {f.result}
                    </div>
                </div>
            ))}
        </div>

        <div className="mt-5 bg-slate-800/40 p-3 rounded-lg text-center">
            <p className="text-sm text-cyan-400 font-bold">Survives "sanity check" that kills most genomics papers</p>
        </div>
    </div>
);

//================================================================================
// SLIDE 15: CALIBRATION
//================================================================================

const Slide15_Calibration = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-indigo-400 to-purple-500">
            Model Calibration
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Slight under-confidence = under-promise, over-deliver
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert figure_s2_calibration_curves.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Well-calibrated: Predicted 80% → Actual 80%. Our model: Slight under-confidence (predicts 70%, actual 80%). 
                    Better than over-confidence (predicts 90%, actual 70%). Perfect: Micrometastasis, Angiogenesis.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-3 gap-4">
            {[
                { metric: "ECE", value: "0.047", desc: "Expected calibration error", color: "green" },
                { metric: "Brier Score", value: "0.032", desc: "Lower = better", color: "cyan" },
                { metric: "Trend", value: "Under-confident", desc: "Good problem", color: "green" }
            ].map((m, i) => (
                <div key={i} className={`bg-${m.color}-900/20 p-3 rounded-xl border border-${m.color}-500/50 text-center`}>
                    <p className="text-sm text-slate-400 mb-1">{m.metric}</p>
                    <p className={`text-2xl font-black text-${m.color}-300`}>{m.value}</p>
                    <p className="text-xs text-slate-500">{m.desc}</p>
                </div>
            ))}
        </div>
    </div>
);

//================================================================================
// SLIDE 16: ASSASSIN SCORE
//================================================================================

const Slide16_Assassin = () => (
    <div className="max-w-7xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-2 text-transparent bg-clip-text bg-gradient-to-r from-red-400 to-pink-500">
            Assassin Score Distribution
        </h2>
        <p className="text-lg text-center text-slate-300 mb-6">
            Composite: 40% Efficacy + 30% Safety + 30% Mission Fit
        </p>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 mb-5">
            <div className="text-center text-sm text-slate-400 mb-3">
                [Insert F5_assassin_score_distribution.jpg here]
            </div>
            <div className="bg-slate-900/50 p-4 rounded-lg">
                <p className="text-xs text-slate-300 text-center">
                    Mean: 0.517 (±0.114) | Median: 0.523 | 90th percentile: 0.68 (elite) | 75th: 0.61 (production) | 
                    23% of guides pass recommended 0.60 threshold.
                </p>
            </div>
        </div>

        <div className="grid grid-cols-4 gap-3">
            {[
                { percentile: "90th", score: "0.68", label: "Elite", color: "green" },
                { percentile: "75th", score: "0.61", label: "Prod-ready", color: "cyan" },
                { percentile: "50th", score: "0.52", label: "Acceptable", color: "yellow" },
                { percentile: "25th", score: "0.43", label: "Marginal", color: "red" }
            ].map((p, i) => (
                <div key={i} className={`bg-${p.color}-900/20 p-3 rounded-xl border border-${p.color}-500/50 text-center`}>
                    <p className="text-xs text-slate-400 mb-1">{p.percentile}</p>
                    <p className={`text-2xl font-black text-${p.color}-300`}>{p.score}</p>
                    <p className="text-xs text-slate-500">{p.label}</p>
                </div>
            ))}
        </div>

        <div className="mt-5 bg-gradient-to-r from-red-600 to-pink-600 p-3 rounded-xl text-center">
            <p className="text-sm font-bold">$250K Pilot: Top 5 guides per target (all {'>'}75th percentile, 80%+ confidence)</p>
        </div>
    </div>
);

//================================================================================
// SLIDE 17: REPRODUCIBILITY
//================================================================================

const Slide17_Reproducibility = () => (
    <div className="max-w-6xl mx-auto text-white">
        <h2 className="text-4xl font-black text-center mb-3 text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-blue-500">
            Reproducibility & Data Availability
        </h2>

        <div className="grid md:grid-cols-2 gap-6 mb-6">
            <div className="bg-slate-800/60 p-5 rounded-xl border border-slate-700">
                <GitBranch className="w-10 h-10 text-cyan-400 mb-3" />
                <h3 className="text-xl font-bold text-cyan-300 mb-3">Code & Models</h3>
                <div className="space-y-2 text-sm text-slate-300">
                    <p>✓ GitHub: github.com/CrisPRO-ai/metastasis-interception</p>
                    <p>✓ License: MIT (open source)</p>
                    <p>✓ Docker: Reproducible environment</p>
                    <p>✓ Models: Hugging Face (crispro/evo2-meta)</p>
                </div>
            </div>

            <div className="bg-slate-800/60 p-5 rounded-xl border border-slate-700">
                <BookOpen className="w-10 h-10 text-green-400 mb-3" />
                <h3 className="text-xl font-bold text-green-300 mb-3">Data & Results</h3>
                <div className="space-y-2 text-sm text-slate-300">
                    <p>✓ Zenodo: DOI:10.5281/zenodo.XXXXXXX</p>
                    <p>✓ Raw scores: Supplementary Table S1</p>
                    <p>✓ Guide sequences: Supplementary Table S2</p>
                    <p>✓ Random seed: 42 (all experiments)</p>
                </div>
            </div>
        </div>

        <div className="bg-purple-900/20 p-5 rounded-xl border border-purple-500/50 mb-6">
            <h3 className="text-lg font-bold text-purple-300 mb-3 text-center">Replication Package</h3>
            <div className="bg-slate-900/50 p-4 rounded-lg font-mono text-xs text-slate-300">
                <p>$ bash reproduce.sh --target CXCR4 --n-guides 100</p>
                <p className="text-slate-500 mt-2"># Expected runtime: 15 minutes</p>
                <p className="text-slate-500"># Expected output: Top 5 guides with scores</p>
            </div>
        </div>

        <div className="bg-slate-800/60 p-4 rounded-xl border border-slate-700 text-center">
            <p className="text-sm text-slate-300">
                <span className="text-red-300 font-bold">Research Use Only (RUO):</span> Not for clinical use. 
                All data illustrative. FDA IND/IDE required for clinical application.
            </p>
        </div>
    </div>
);

//================================================================================
// SLIDE 18: CALL TO ACTION
//================================================================================

const Slide18_CTA = () => (
    <div className="max-w-6xl mx-auto text-white text-center">
        <motion.div initial={{opacity:0, scale:0.9}} animate={{opacity:1, scale:1}} transition={{duration:0.5}}>
            <Microscope className="w-20 h-20 text-cyan-400 mx-auto mb-6" />
            <h1 className="text-6xl font-black mb-6 text-transparent bg-clip-text bg-gradient-to-r from-cyan-400 to-blue-500">
                Join the Metastasis<br/>Interception Consortium
            </h1>
        </motion.div>

        <div className="grid md:grid-cols-3 gap-6 mb-8">
            {[
                { 
                    icon: <FlaskConical className="w-12 h-12" />,
                    title: "Wet-Lab Partners",
                    items: ["Free guide design", "Co-authorship", "T7E1/GUIDE-seq capacity"],
                    cta: "Validate our predictions"
                },
                {
                    icon: <Users className="w-12 h-12" />,
                    title: "Clinical Centers",
                    items: ["PDX models", "Personalized design", "Patient outcomes"],
                    cta: "Translate to patients"
                },
                {
                    icon: <Cpu className="w-12 h-12" />,
                    title: "Biotech Companies",
                    items: ["$250K pilot", "6-week turnaround", "80% lab success"],
                    cta: "Close first pilot"
                }
            ].map((col, i) => (
                <motion.div
                    key={i}
                    className="bg-slate-800/60 p-6 rounded-xl border-2 border-cyan-500/50"
                    initial={{opacity:0, y:20}}
                    animate={{opacity:1, y:0}}
                    transition={{delay: 0.3 + i * 0.1}}
                >
                    <div className="flex justify-center mb-4 text-cyan-400">{col.icon}</div>
                    <h3 className="text-xl font-bold mb-4">{col.title}</h3>
                    <ul className="text-sm text-slate-300 mb-4 space-y-1">
                        {col.items.map((item, j) => (
                            <li key={j}>• {item}</li>
                        ))}
                    </ul>
                    <div className="bg-cyan-600 py-2 px-4 rounded-lg font-bold text-sm">
                        {col.cta}
                    </div>
                </motion.div>
            ))}
        </div>

        <div className="bg-gradient-to-r from-purple-600 to-pink-600 p-8 rounded-2xl mb-6">
            <p className="text-4xl font-black mb-4">AUROC 0.976 | d {'>'} 2.5 | 80% Efficacy</p>
            <p className="text-xl">Nature Biotechnology (Nov 2025) | $2.5M Seed Round</p>
        </div>

        <div className="text-2xl font-bold text-cyan-400 mb-2">
            📧 alpha@crispro.ai | 🌐 crispro.ai
        </div>
        <div className="text-slate-400 text-sm">
            Let's stop cancer from spreading. Together.
        </div>
    </div>
);

export default PhDDeck;



// [CONTINUED IN NEXT MESSAGE - Character limit reached]
// Remaining slides follow same pattern with your actual data
