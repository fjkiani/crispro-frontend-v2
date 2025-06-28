/**
 * @file KnowledgeGraphExample.tsx
 * @description This component demonstrates the KnowledgeGraph's capabilities, specifically
 * showcasing how CrisPRO.ai's AI-driven insights and therapeutic context transform
 * standard biological network visualization into an active decision-support tool.
 *
 * CrisPRO.ai Context and Purpose - How CrisPRO Enhances This Visualization:
 * - **AI-Driven Insight Integration:** Illustrates how CrisPRO.ai populates nodes and edges
 *   with AI-generated scores (`aiRelevanceScore`, `predictedImpact`, `confidenceScore`) and
 *   qualitative summaries (`insightSummary`), providing deeper, context-aware understanding
 *   beyond simple relationship mapping.
 * - **Therapeutic Contextualization:** The `therapeuticContext` prop demonstrates CrisPRO.ai's
 *   core ability to filter, prioritize, and annotate graph elements based on specific
 *   user-defined goals (e.g., disease, patient cohort, prophylactic vs. therapeutic intent),
 *   making the visualization directly relevant to the research question.
 * - **Interactive Exploration of AI Outputs:** Allows users to see how CrisPRO.ai's
 *   AI-generated data (e.g., LLM summaries for selected nodes) is presented interactively,
 *   facilitating deeper dives into complex relationships.
 * - **Dynamic Risk & Relevance Assessment:** Shows how AI-weighted edges and nodes can
 *   dynamically shift focus based on CrisPRO.ai's ongoing analysis (e.g., `aiRelevanceScore`
 *   changing based on new literature or simulation results in a real app).
 * - **Foundation for Advanced Features:** This example serves as a testbed for UI patterns
 *   (like the node detail panel) and interaction models (like toggling AI features) that
 *   are central to CrisPRO.ai's goal of making complex CRISPR therapeutic design more
 *   intuitive and data-driven.
 * - **Pathway Illumination Simulation:** The `highlightKeyPathways` toggle simulates how
 *   CrisPRO.ai might proactively guide users by illuminating therapeutically relevant
 *   pathways or chains of evidence within the graph.
 *
 * This example is crucial for understanding how the KnowledgeGraph component will be
 * leveraged within CrisPRO.ai to visualize and interact with the platform's unique
 * AI-generated therapeutic intelligence.
 */
'use client';

import React, { useState, useCallback, useEffect } from 'react';
import KnowledgeGraph, { GraphNode as BaseGraphNode, GraphEdge as BaseGraphEdge } from './KnowledgeGraph'; // Assuming KnowledgeGraph is in the same directory

// Define extended Node and Edge types for CrisPRO.ai specific data
export interface CrisPROGraphNode extends BaseGraphNode {
  baseEvidenceScore?: number; // Renamed from weight
  aiRelevanceScore?: number;
  predictedImpact?: number;
  vusClassification?: 'Pathogenic' | 'Likely Pathogenic' | 'VUS' | 'Likely Benign' | 'Benign';
  crisprTargetabilityScore?: number;
  insightSummary?: string;
  llmGeneratedSummaryTimestamp?: string;
  therapeuticContext?: string; // Context specific to this node
  evidenceStrength?: number; // Kept from original, CrisPRO can also populate this
  prophylacticSuitabilityScore?: number;
  druggabilityIndex?: number;
}

export interface CrisPROGraphEdge extends BaseGraphEdge {
  confidenceScore?: number;
  evidenceSourceType?: 'Literature' | 'ClinicalTrial' | 'CrisPRO_Simulation' | 'CrisPRO_Prediction' | 'User_Annotated';
  mechanismDetails?: string;
  isSynergistic?: boolean;
  isAntagonistic?: boolean;
  inTherapeuticPath?: boolean; // Highlight if edge is part of a key therapeutic pathway
  clinicalRelevance?: number; // Kept from original
  literatureSupport?: number; // Kept from original
}


/**
 * @function KnowledgeGraphExample
 * @description A functional React component that demonstrates the usage of the
 * `KnowledgeGraph` component with sample data and interactive controls,
 * enhanced to showcase CrisPRO.ai's specific data and functionalities.
 */
export default function KnowledgeGraphExample() {
  const initialNodes: CrisPROGraphNode[] = [
    {
      id: '1',
      label: 'Robson et al. 2017',
      type: 'publication',
      description: 'Clinical trial of PARP inhibitors in BRCA-mutated breast cancer',
      baseEvidenceScore: 0.8,
      aiRelevanceScore: 0.92,
      therapeuticContext: 'BRCA-mutated breast cancer',
      evidenceStrength: 0.85,
      llmGeneratedSummaryTimestamp: '2023-10-26T14:00:00Z',
      metadata: { journal: 'NEJM', year: 2017 }
    },
    {
      id: '2',
      label: 'PARP Inhibitors',
      type: 'therapy',
      description: 'Poly (ADP-ribose) polymerase inhibitors used in treatment of BRCA-mutated cancers',
      baseEvidenceScore: 0.9,
      aiRelevanceScore: 0.95,
      insightSummary: 'CrisPRO Agent Insight: Highly relevant for BRCA1/2 mutation carriers. Clinical trials show significant PFS benefits. LLM analysis indicates emerging resistance mechanisms to monitor.',
      therapeuticContext: 'BRCA-mutated breast cancer',
      evidenceStrength: 0.9,
      druggabilityIndex: 0.85,
      llmGeneratedSummaryTimestamp: '2023-10-27T11:00:00Z',
      metadata: { drugClass: 'Enzyme Inhibitor'}
    },
    {
      id: '3',
      label: 'Breast Cancer',
      type: 'outcome', // Consider 'disease' type if more appropriate
      description: 'Malignant breast neoplasm',
      baseEvidenceScore: 1.0,
      aiRelevanceScore: 0.97,
      therapeuticContext: 'BRCA-mutated breast cancer',
      evidenceStrength: 0.98,
      metadata: { prevalence: 'High in BRCA mutation carriers'}
    },
    {
      id: '4',
      label: 'BRCA1 (Variant rs28397696A1)',
      type: 'variant',
      description: 'Pathogenic BRCA1 variant linked to increased breast cancer risk.',
      baseEvidenceScore: 0.85,
      aiRelevanceScore: 0.90,
      predictedImpact: 0.95,
      vusClassification: 'Pathogenic',
      crisprTargetabilityScore: 0.88,
      therapeuticContext: 'BRCA-mutated breast cancer',
      insightSummary: 'CrisPRO Analysis: Pathogenic variant. High predicted impact (0.95). Strong candidate for prophylactic intervention via HDR. Off-target risk assessment by CrisPRO suggests 3 high-priority sites for validation. (LLM-gen v1.3)',
      llmGeneratedSummaryTimestamp: '2023-10-27T10:30:00Z',
      evidenceStrength: 0.92,
      prophylacticSuitabilityScore: 0.9,
      metadata: { gene: 'BRCA1', allele: 'A1'}
    },
    {
      id: '5',
      label: 'Ovarian Cancer',
      type: 'outcome', // Consider 'disease'
      description: 'Malignant ovarian neoplasm, often associated with BRCA mutations.',
      baseEvidenceScore: 0.75,
      aiRelevanceScore: 0.82,
      therapeuticContext: 'BRCA-mutated cancers',
      evidenceStrength: 0.80,
      metadata: { relatedTo: 'BRCA1/2'}
    },
    {
      id: '6',
      label: 'BRCA2 (Variant rs5991231)',
      type: 'variant',
      description: 'Pathogenic BRCA2 variant, common in hereditary breast and ovarian cancer.',
      baseEvidenceScore: 0.80,
      aiRelevanceScore: 0.88,
      predictedImpact: 0.88,
      vusClassification: 'Pathogenic',
      crisprTargetabilityScore: 0.85,
      therapeuticContext: 'BRCA-mutated cancers',
      insightSummary: 'CrisPRO Analysis: Pathogenic. Similar prophylactic potential to BRCA1 variant but may require different gRNA design due to sequence context. CrisPRO Simulation Core predicts 80% correction efficiency with proposed strategy.',
      llmGeneratedSummaryTimestamp: '2023-10-27T12:15:00Z',
      evidenceStrength: 0.85,
      metadata: { gene: 'BRCA2'}
    },
    {
      id: '7',
      label: 'Platinum Chemotherapy',
      type: 'therapy',
      description: 'Platinum-based chemotherapy agents, effective in BRCA-deficient tumors.',
      baseEvidenceScore: 0.7,
      aiRelevanceScore: 0.75,
      therapeuticContext: 'BRCA-mutated cancers',
      insightSummary: 'CrisPRO Agent Insight: Standard of care. LLM analysis of recent literature suggests PARP inhibitor combination significantly improves outcomes.',
      evidenceStrength: 0.78,
      metadata: { mechanism: 'DNA cross-linking' }
    },
    {
      id: '8',
      label: 'Turner et al. 2019',
      type: 'publication',
      description: 'Study on synthetic lethality in BRCA-deficient cells, supporting PARP inhibitor use.',
      baseEvidenceScore: 0.65,
      aiRelevanceScore: 0.70,
      therapeuticContext: 'BRCA-mutated breast cancer',
      evidenceStrength: 0.72,
      metadata: { journal: 'Lancet Oncology', year: 2019 }
    },
    {
      id: '9',
      label: 'DNA Repair Pathway (HR)',
      type: 'pathway',
      description: 'Homologous Recombination pathway, critical for repairing double-strand breaks. Often compromised in BRCA-mutated cells.',
      baseEvidenceScore: 0.9,
      aiRelevanceScore: 0.93,
      therapeuticContext: 'BRCA-mutated cancers',
      insightSummary: 'CrisPRO Agent Insight: Key vulnerability. PARP inhibitors exploit HR deficiency. CrisPRO can model impact of edits on this pathway.',
      metadata: { process: 'Homologous Recombination' }
    },
    {
      id: '10',
      label: 'Cell Cycle Checkpoint (G2/M)',
      type: 'pathway',
      description: 'Regulatory points in the cell cycle ensuring proper cell division. Often dysregulated in cancer.',
      baseEvidenceScore: 0.7,
      aiRelevanceScore: 0.65,
      therapeuticContext: 'General Oncology',
      metadata: { involvedIn: 'Cancer Progression' }
    }
  ];

  const initialEdges: CrisPROGraphEdge[] = [
    { id: 'e1', source: '1', target: '2', type: 'REPORTS_ON_THERAPY', weight: 0.9, confidenceScore: 0.95, label: 'Reports Use Of', clinicalRelevance: 0.9, literatureSupport: 0.88, evidenceSourceType: 'Literature' },
    { id: 'e2', source: '2', target: '3', type: 'TREATS_DISEASE', weight: 0.85, confidenceScore: 0.90, label: 'Treats', clinicalRelevance: 0.92, literatureSupport: 0.9, inTherapeuticPath: true, evidenceSourceType: 'ClinicalTrial' },
    { id: 'e3', source: '4', target: '3', type: 'INCREASES_RISK_FOR_DISEASE', weight: 0.95, confidenceScore: 0.92, label: 'Increases Risk Of', clinicalRelevance: 0.95, literatureSupport: 0.93, inTherapeuticPath: true, evidenceSourceType: 'Literature' },
    { id: 'e4', source: '2', target: '5', type: 'TREATS_DISEASE', weight: 0.80, confidenceScore: 0.85, label: 'Treats', clinicalRelevance: 0.80, literatureSupport: 0.82, evidenceSourceType: 'ClinicalTrial' },
    { id: 'e5', source: '6', target: '5', type: 'INCREASES_RISK_FOR_DISEASE', weight: 0.90, confidenceScore: 0.88, label: 'Increases Risk Of', clinicalRelevance: 0.87, literatureSupport: 0.85, inTherapeuticPath: true, evidenceSourceType: 'Literature'},
    { id: 'e6', source: '7', target: '5', type: 'TREATS_DISEASE', weight: 0.75, confidenceScore: 0.80, label: 'Treats', clinicalRelevance: 0.78, literatureSupport: 0.75, evidenceSourceType: 'ClinicalTrial' },
    { id: 'e7', source: '8', target: '2', type: 'SUPPORTS_THERAPY_USE', weight: 0.70, confidenceScore: 0.75, label: 'Supports Use Of', clinicalRelevance: 0.7, literatureSupport: 0.65, evidenceSourceType: 'Literature' },
    { id: 'e8', source: '7', target: '3', type: 'TREATS_DISEASE', weight: 0.65, confidenceScore: 0.70, label: 'Treats', clinicalRelevance: 0.68, literatureSupport: 0.62, evidenceSourceType: 'ClinicalTrial' },
    { id: 'e9', source: '4', target: '6', type: 'ASSOCIATED_WITH_VARIANT', weight: 0.60, confidenceScore: 0.65, label: 'Co-occurs With', clinicalRelevance: 0.5, literatureSupport: 0.55, evidenceSourceType: 'User_Annotated' }, // Example for co-occurrence
    { id: 'e10', source: '1', target: '8', type: 'RELATED_PUBLICATION', weight: 0.50, confidenceScore: 0.55, label: 'Related Research', clinicalRelevance: 0.4, literatureSupport: 0.45, evidenceSourceType: 'Literature' },
    { id: 'e11', source: '4', target: '9', type: 'IMPACTS_PATHWAY', weight: 0.88, confidenceScore: 0.90, label: 'Impacts HR Pathway', clinicalRelevance: 0.92, literatureSupport: 0.89, inTherapeuticPath: true, mechanismDetails: 'BRCA1 loss impairs Homologous Recombination.', evidenceSourceType: 'CrisPRO_Prediction' },
    { id: 'e12', source: '2', target: '9', type: 'TARGETS_PATHWAY_VULNERABILITY', weight: 0.92, confidenceScore: 0.94, label: 'Targets HR Deficiency', clinicalRelevance: 0.95, literatureSupport: 0.91, inTherapeuticPath: true, mechanismDetails: 'PARP inhibitors exploit synthetic lethality in HR-deficient cells.', evidenceSourceType: 'Literature' },
    { id: 'e13', source: '9', target: '10', type: 'INFLUENCES_PATHWAY', weight: 0.70, confidenceScore: 0.68, label: 'Influences Cell Cycle', clinicalRelevance: 0.60, literatureSupport: 0.65, evidenceSourceType: 'Literature' }
  ];

  const [nodes, setNodes] = useState<CrisPROGraphNode[]>(initialNodes);
  const [edges, setEdges] = useState<CrisPROGraphEdge[]>(initialEdges);
  const [isLoadingLLMInsight, setIsLoadingLLMInsight] = useState(false);
  const [selectedNode, setSelectedNode] = useState<CrisPROGraphNode | null>(null);
  const [selectedEdge, setSelectedEdge] = useState<CrisPROGraphEdge | null>(null);


  const [physicsConfig, setPhysicsConfig] = useState({
    repulsion: 6000,
    springLength: 280,
    damping: 0.85,
    maxVelocity: 35
  });

  const [nodeSpacing, setNodeSpacing] = useState(80);
  const [graphKey, setGraphKey] = useState(Date.now());
  const [simulationActive, setSimulationActive] = useState(true);

  const [enableAIFeatures, setEnableAIFeatures] = useState({
    showCrisPRORelevance: true, // Renamed from aiWeighting
    enableLLMInsights: true,
    highlightTherapeuticPathways: true, // Renamed and default to true
    showPredictedInterventionEffects: false, // New AI Feature
  });

  const [therapeuticContextOptions] = useState([
    'BRCA-mutated Cancers (General)',
    'Prophylactic - BRCA1 High Risk',
    'Therapeutic - Lung Cancer (KRAS G12C)',
    'VUS Investigation - Gene TP53',
    'General Oncology',
  ]);
  const [currentTherapeuticContext, setCurrentTherapeuticContext] = useState(therapeuticContextOptions[0]);

  // Effect to potentially filter/update nodes/edges when therapeutic context changes
  // This is a placeholder for more complex CrisPRO logic
  useEffect(() => {
    console.log("Therapeutic Context Changed:", currentTherapeuticContext);
    // In a real CrisPRO app, this would trigger re-filtering of nodes/edges,
    // re-calculation of aiRelevanceScores, or highlighting based on the context.
    // For this example, we can just re-key the graph to show a visual "refresh"
    // Or, more advanced: dynamically update aiRelevanceScore based on context.
    setNodes(prevNodes => prevNodes.map(n => ({
        ...n,
        // Example: Boost relevance if node's context matches current global context
        aiRelevanceScore: n.therapeuticContext === currentTherapeuticContext || n.therapeuticContext?.includes('General')
                          ? Math.min(1, (n.aiRelevanceScore || 0.5) + 0.1) // Boost
                          : Math.max(0, (n.aiRelevanceScore || 0.5) - 0.1)  // Attenuate
    })));
    setGraphKey(Date.now()); // Force re-render of graph with potentially updated styles
  }, [currentTherapeuticContext]);


  const handleNodeClick = useCallback((node: BaseGraphNode) => { // Use BaseGraphNode for handler type
    const crisproNode = node as CrisPROGraphNode; // Cast to CrisPROGraphNode for use
    setSelectedNode(crisproNode);
    setSelectedEdge(null); // Clear edge selection

    if (enableAIFeatures.enableLLMInsights && crisproNode) {
      setIsLoadingLLMInsight(true);
      setTimeout(() => {
        setSelectedNode(prevNode => {
          if (prevNode && prevNode.id === crisproNode.id) {
            let insight = `(Mock CrisPRO LLM Insight for ${crisproNode.label} | Context: ${currentTherapeuticContext}):
`;
            if (crisproNode.type === 'variant' && crisproNode.predictedImpact && crisproNode.predictedImpact > 0.7) {
              insight += `High predicted pathogenicity (${(crisproNode.predictedImpact * 100).toFixed(0)}%). CrisPRO Therapeutic Strategy Agent suggests exploring HDR for correction if prophylactic context or precise knockout if tumor suppression. Key off-target risks identified in genes X, Y, Z. Consider high-fidelity nucleases.`;
            } else if (crisproNode.type === 'therapy') {
              insight += `This therapy shows high CrisPRO AI Relevance for ${currentTherapeuticContext}. CrisPRO Literature Analysis Module identified 3 new supporting publications this week. CrisPRO Digital Twin simulation predicts 75% efficacy with manageable immunogenicity for this patient profile.`;
            } else if (crisproNode.type === 'pathway') {
                insight += `This pathway is critical in ${currentTherapeuticContext}. CrisPRO Simulation Core can model the effects of interventions targeting this pathway. Current `aiRelevanceScore` is high due to its central role.`;
            }
            else {
              insight += `Standard LLM summary. Further CrisPRO analysis can provide deeper insights into its role within the current therapeutic context.`;
            }
            return { ...prevNode, insightSummary: insight, llmGeneratedSummaryTimestamp: new Date().toISOString() };
          }
          return prevNode;
        });
        setIsLoadingLLMInsight(false);
      }, 1500);
    }
  }, [enableAIFeatures.enableLLMInsights, currentTherapeuticContext]);

  const handleEdgeClick = useCallback((edge: BaseGraphEdge | null) => {
    setSelectedEdge(edge as CrisPROGraphEdge | null);
    setSelectedNode(null); // Clear node selection
     if (enableAIFeatures.enableLLMInsights && edge) {
      setIsLoadingLLMInsight(true);
      const crisproEdge = edge as CrisPROGraphEdge;
      setTimeout(() => {
        setSelectedEdge(prevEdge => {
          if (prevEdge && prevEdge.id === crisproEdge.id) {
            const sourceNode = nodes.find(n => n.id === crisproEdge.source);
            const targetNode = nodes.find(n => n.id === crisproEdge.target);
            let insight = `(Mock CrisPRO LLM Insight for edge: ${sourceNode?.label} -> ${targetNode?.label} | Context: ${currentTherapeuticContext}):
`;
            insight += `This '${crisproEdge.type}' relationship has a CrisPRO confidence score of ${(crisproEdge.confidenceScore!*100).toFixed(0)}%. Evidence type: ${crisproEdge.evidenceSourceType}. `;
            if (crisproEdge.inTherapeuticPath) {
                insight += `Marked as part of a key therapeutic pathway for ${currentTherapeuticContext}. `;
            }
            if (crisproEdge.mechanismDetails) {
                insight += `Mechanism: ${crisproEdge.mechanismDetails}`;
            }
            // You can add more logic here for different edge types or contexts
            return { ...prevEdge, insightSummary: insight } as CrisPROGraphEdge; // Add insightSummary to edge if needed
          }
          return prevEdge;
        });
        setIsLoadingLLMInsight(false);
      }, 1200);
    }
  }, [enableAIFeatures.enableLLMInsights, currentTherapeuticContext, nodes]);


  const handleResetLayout = useCallback(() => {
    setGraphKey(Date.now());
    setSimulationActive(false);
    setSelectedNode(null);
    setSelectedEdge(null);
  }, []);

  const toggleSimulation = useCallback(() => {
    setSimulationActive(prev => !prev);
  }, []);

  const handlePhysicsChange = (param: keyof typeof physicsConfig, value: number) => {
    setPhysicsConfig(prev => ({ ...prev, [param]: value }));
  };

  const handleAIFeatureToggle = (feature: keyof typeof enableAIFeatures) => {
    setEnableAIFeatures(prev => ({ ...prev, [feature]: !prev[feature] }));
     if (feature === 'highlightTherapeuticPathways'){
        setEdges(prevEdges => prevEdges.map(e => ({...e}))); // Trigger re-render for path highlighting
    }
  };

  return (
    <div className="p-4 md:p-6 bg-gray-950 text-white min-h-screen flex flex-col">
      <header className="mb-6">
        <h1 className="text-3xl font-bold text-sky-400">CrisPRO.ai Knowledge Graph Explorer</h1>
        <p className="text-sm text-gray-400">
          Visualizing AI-enhanced biological networks for CRISPR therapeutic design.
        </p>
      </header>

      <div className="mb-6 p-4 bg-slate-800 rounded-lg shadow-lg">
        <h2 className="text-xl font-semibold mb-3 text-sky-300">Graph Configuration & CrisPRO Controls</h2>

        <div className="mb-4">
            <label htmlFor="therapeuticContextInput" className="block text-sm text-gray-300 mb-1">CrisPRO Therapeutic Context:</label>
            <select
              id="therapeuticContextInput"
              value={currentTherapeuticContext}
              onChange={(e) => setCurrentTherapeuticContext(e.target.value)}
              className="w-full p-2 rounded bg-slate-700 border border-slate-600 focus:ring-sky-500 focus:border-sky-500"
            >
              {therapeuticContextOptions.map(option => (
                <option key={option} value={option}>{option}</option>
              ))}
            </select>
            <p className="text-xs text-gray-500 mt-1">
              CrisPRO.ai: Select a context. In the full application, this is set via the 'Therapeutic Context Enabled Mode' and deeply influences all AI analyses and visualizations.
            </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 mb-4">
          <div>
            <h3 className="text-md font-medium mb-2 text-gray-300">Physics Parameters:</h3>
            {Object.entries(physicsConfig).map(([key, value]) => (
              <div key={key} className="mb-2">
                <label className="block text-xs text-gray-400 capitalize">{key.replace(/([A-Z])/g, ' $1')}</label>
                <input
                  type="range"
                  min={key === 'damping' ? 0.5 : (key === 'repulsion' ? 500 : 50)}
                  max={key === 'damping' ? 0.99 : (key === 'repulsion' ? 10000 : (key === 'springLength' ? 500 : 100))}
                  step={key === 'damping' ? 0.01 : (key === 'repulsion' ? 100 : 10)}
                  value={value}
                  onChange={(e) => handlePhysicsChange(key as keyof typeof physicsConfig, parseFloat(e.target.value))}
                  className="w-full h-2 bg-slate-700 rounded-lg appearance-none cursor-pointer accent-sky-500"
                />
                <div className="text-xs text-right text-gray-400">{value.toFixed(key === 'damping' ? 2 : 0)}</div>
              </div>
            ))}
             <div className="mb-2">
                <label className="block text-xs text-gray-400">Node Spacing</label>
                <input
                  type="range" min="20" max="150" step="10" value={nodeSpacing}
                  onChange={(e) => setNodeSpacing(parseInt(e.target.value))}
                  className="w-full h-2 bg-slate-700 rounded-lg appearance-none cursor-pointer accent-sky-500"
                />
                <div className="text-xs text-right text-gray-400">{nodeSpacing}px</div>
              </div>
          </div>

          <div>
            <h3 className="text-md font-medium mb-2 text-gray-300">CrisPRO.ai Features:</h3>
            {Object.entries(enableAIFeatures).map(([key, value]) => (
              <div key={key} className="flex items-center justify-between mb-2">
                <label htmlFor={key} className="text-sm text-gray-400 capitalize">{key.replace(/([A-Z])/g, ' $1')}</label>
                <button
                  id={key}
                  onClick={() => handleAIFeatureToggle(key as keyof typeof enableAIFeatures)}
                  className={`px-3 py-1 text-xs rounded-md transition-colors ${value ? 'bg-green-500 hover:bg-green-600' : 'bg-slate-600 hover:bg-slate-500'}`}
                >
                  {value ? 'Enabled' : 'Disabled'}
                </button>
              </div>
            ))}
            <p className="text-xs text-gray-500 mt-2">
              CrisPRO.ai: Toggle AI-driven enhancements like relevance display, LLM insights, pathway highlighting, and simulated intervention effects.
            </p>
          </div>

          <div>
            <h3 className="text-md font-medium mb-2 text-gray-300">Actions:</h3>
            <div className="flex flex-col space-y-2">
              <button
                onClick={toggleSimulation}
                className={`${simulationActive ? 'bg-orange-500 hover:bg-orange-600' : 'bg-teal-500 hover:bg-teal-600'} text-white px-4 py-2 rounded-md text-sm transition-colors`}
              >
                {simulationActive ? 'Pause Simulation' : 'Start Simulation'}
              </button>
              <button
                onClick={handleResetLayout}
                className="bg-sky-600 hover:bg-sky-700 text-white px-4 py-2 rounded-md text-sm transition-colors"
              >
                Reset Layout
              </button>
            </div>
          </div>
        </div>
      </div>

      <div className="flex flex-col md:flex-row gap-6 flex-grow">
        <div className="flex-grow md:w-2/3 h-[600px] md:h-auto border border-slate-700 rounded-lg shadow-2xl">
          {/* Apply subtle glow if AI features are active */}
          <KnowledgeGraph
            key={graphKey}
            nodes={nodes}
            edges={edges}
            enableDragging={true}
            enableZoom={true}
            showLabels={true}
            showEdgeLabels={true}
            usePhysics={true}
            physicsConfig={physicsConfig}
            nodeSpacing={nodeSpacing}
            onNodeClick={handleNodeClick}
            onEdgeClick={handleEdgeClick} // Added edge click handler
            simulationRunning={simulationActive}
            therapeuticContext={currentTherapeuticContext}
            enableAICustomizations={{ // Pass AI feature toggles
                showCrisPRORelevance: enableAIFeatures.showCrisPRORelevance,
                enableLLMInsights: enableAIFeatures.enableLLMInsights,
                highlightTherapeuticPathways: enableAIFeatures.highlightTherapeuticPathways,
                showPredictedInterventionEffects: enableAIFeatures.showPredictedInterventionEffects,
            }}
            className={`bg-slate-900 rounded-lg w-full h-full transition-all duration-500 
                        ${(enableAIFeatures.showCrisPRORelevance || enableAIFeatures.highlightTherapeuticPathways) ? 'shadow-[0_0_15px_3px_rgba(56,189,248,0.3)]' : ''}`}
          />
        </div>

        <aside className="w-full md:w-1/3 lg:w-1/4 bg-slate-800 p-4 rounded-lg shadow-lg h-full overflow-y-auto">
          <h2 className="text-xl font-semibold mb-3 text-sky-300">
            {selectedNode ? "Node Details" : selectedEdge ? "Edge Details" : "Selection Details"}
          </h2>

          {selectedNode && (
            <div className="space-y-3">
              <h3 className="text-lg font-bold text-sky-400">{selectedNode.label}</h3>
              <p className="text-xs uppercase text-gray-400 tracking-wider">Type: {selectedNode.type}</p>
              {selectedNode.description && (
                <p className="text-sm text-gray-300">{selectedNode.description}</p>
              )}

              <div className="grid grid-cols-2 gap-x-4 mt-2">
                {selectedNode.baseEvidenceScore !== undefined && (
                  <div>
                    <label className="text-xs text-gray-400">Base Evidence Score:</label>
                    <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                      <div
                        className="bg-blue-500 h-2.5 rounded-full transition-all duration-300"
                        style={{ width: `${selectedNode.baseEvidenceScore * 100}%` }}
                        title={`${(selectedNode.baseEvidenceScore * 100).toFixed(0)}%`}
                      />
                    </div>
                  </div>
                )}
                {enableAIFeatures.showCrisPRORelevance && selectedNode.aiRelevanceScore !== undefined && (
                   <div>
                    <label className="text-xs text-gray-400">CrisPRO.ai Relevance:</label>
                    <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                      <div
                        className="bg-green-500 h-2.5 rounded-full transition-all duration-300"
                        style={{ width: `${selectedNode.aiRelevanceScore * 100}%` }}
                        title={`${(selectedNode.aiRelevanceScore * 100).toFixed(0)}%`}
                      />
                    </div>
                  </div>
                )}
              </div>
              
              {enableAIFeatures.showCrisPRORelevance && selectedNode.evidenceStrength !== undefined && (
                 <div className="mt-2">
                  <label className="text-xs text-gray-400">CrisPRO.ai Evidence Strength:</label>
                  <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                    <div
                      className="bg-teal-500 h-2.5 rounded-full transition-all duration-300"
                      style={{ width: `${selectedNode.evidenceStrength * 100}%` }}
                      title={`${(selectedNode.evidenceStrength * 100).toFixed(0)}%`}
                    />
                  </div>
                </div>
              )}
               {selectedNode.vusClassification && (
                <p className="text-sm mt-2">
                    <span className="font-semibold text-gray-300">VUS Classification: </span>
                    <span className="font-bold" style={{color: selectedNode.vusClassification?.includes('Pathogenic') ? '#ef4444' : '#22c55e'}}>{selectedNode.vusClassification}</span>
                </p>
              )}
              {selectedNode.crisprTargetabilityScore !== undefined && (
                 <div className="mt-2">
                  <label className="text-xs text-gray-400">CrisPRO.ai CRISPR Targetability:</label>
                  <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                    <div
                      className="bg-purple-500 h-2.5 rounded-full transition-all duration-300"
                      style={{ width: `${selectedNode.crisprTargetabilityScore * 100}%` }}
                      title={`${(selectedNode.crisprTargetabilityScore * 100).toFixed(0)}%`}
                    />
                  </div>
                </div>
              )}


              {enableAIFeatures.enableLLMInsights && (
                <div className="mt-3 pt-3 border-t border-slate-700">
                  <h4 className="text-sm font-semibold mb-1 text-sky-400">CrisPRO.ai Insight:</h4>
                  {isLoadingLLMInsight ? (
                    <p className="text-sm italic text-sky-400 animate-pulse">CrisPRO processing deeper insight...</p>
                  ) : selectedNode.insightSummary ? (
                    <p className="text-sm italic text-gray-300 whitespace-pre-line">{selectedNode.insightSummary}</p>
                  ) : (
                    <p className="text-sm italic text-gray-500">Click node to attempt LLM insight generation.</p>
                  )}
                  {selectedNode.llmGeneratedSummaryTimestamp && !isLoadingLLMInsight && selectedNode.insightSummary &&(
                    <p className="text-xs text-gray-500 mt-1">Insight as of: {new Date(selectedNode.llmGeneratedSummaryTimestamp).toLocaleString()}</p>
                  )}
                </div>
              )}

              {selectedNode.type === 'variant' && selectedNode.predictedImpact !== undefined && (
                <div className="mt-3 pt-3 border-t border-slate-700">
                  <h4 className="text-sm font-semibold mb-1">CrisPRO.ai Predicted Impact (Variant):</h4>
                  <div className="text-lg font-bold" style={{
                    color: selectedNode.predictedImpact > 0.7 ? '#ef4444'
                         : selectedNode.predictedImpact > 0.4 ? '#f59e0b'
                         : '#22c55e'
                  }}>
                    {selectedNode.predictedImpact > 0.7 ? 'High'
                     : selectedNode.predictedImpact > 0.4 ? 'Medium'
                     : 'Low'}
                    <span className="ml-2 text-sm text-gray-400">
                      (Score: {(selectedNode.predictedImpact * 100).toFixed(0)}%)
                    </span>
                  </div>
                </div>
              )}
               {selectedNode.prophylacticSuitabilityScore !== undefined && (
                 <div className="mt-2 pt-3 border-t border-slate-700">
                  <label className="text-xs text-gray-400">CrisPRO.ai Prophylactic Suitability:</label>
                  <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                    <div
                      className="bg-pink-500 h-2.5 rounded-full transition-all duration-300"
                      style={{ width: `${selectedNode.prophylacticSuitabilityScore * 100}%` }}
                      title={`${(selectedNode.prophylacticSuitabilityScore * 100).toFixed(0)}%`}
                    />
                  </div>
                </div>
              )}


              {selectedNode.metadata && Object.keys(selectedNode.metadata).length > 0 && (
                <div className="mt-3 pt-3 border-t border-slate-700">
                  <h4 className="text-sm font-semibold mb-1 text-gray-300">Additional Details:</h4>
                  {Object.entries(selectedNode.metadata).map(([key, value]) => (
                    <p key={key} className="text-xs text-gray-400">
                      <span className="capitalize font-medium">{key.replace(/([A-Z])/g, ' $1')}:</span> {String(value)}
                    </p>
                  ))}
                </div>
              )}

              {selectedNode.therapeuticContext && (
                 <div className="mt-3 pt-3 border-t border-slate-700">
                  <h4 className="text-sm font-semibold mb-1 text-gray-300">Node Specific Therapeutic Context:</h4>
                  <p className="text-xs text-gray-400 bg-slate-700 px-2 py-1 rounded w-fit">
                    {selectedNode.therapeuticContext}
                  </p>
                </div>
              )}
              <div className="mt-4 pt-3 border-t border-slate-700">
                 <button
                    onClick={() => alert(`Triggering deeper CrisPRO analysis for ${selectedNode.label}... (Conceptual)`)}
                    className="w-full bg-sky-600 hover:bg-sky-700 text-white px-3 py-2 rounded-md text-sm transition-colors"
                  >
                    Run Deeper CrisPRO Analysis
                </button>
              </div>

            </div>
          )}

          {selectedEdge && (
             <div className="space-y-3">
              <h3 className="text-lg font-bold text-sky-400">
                Edge: {nodes.find(n=>n.id === selectedEdge.source)?.label} <span className="text-gray-400">→</span> {nodes.find(n=>n.id === selectedEdge.target)?.label}
              </h3>
              <p className="text-xs uppercase text-gray-400 tracking-wider">Type: {selectedEdge.type}</p>
              {selectedEdge.label && <p className="text-sm text-gray-300">Label: {selectedEdge.label}</p>}

              <div className="grid grid-cols-2 gap-x-4 mt-2">
                {selectedEdge.weight !== undefined && (
                  <div>
                    <label className="text-xs text-gray-400">Base Weight:</label>
                     <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                      <div className="bg-blue-500 h-2.5 rounded-full" style={{ width: `${selectedEdge.weight * 100}%`}} title={`${(selectedEdge.weight * 100).toFixed(0)}%`}/>
                    </div>
                  </div>
                )}
                {enableAIFeatures.showCrisPRORelevance && selectedEdge.confidenceScore !== undefined && (
                  <div>
                    <label className="text-xs text-gray-400">CrisPRO.ai Confidence:</label>
                    <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                      <div className="bg-green-500 h-2.5 rounded-full" style={{ width: `${selectedEdge.confidenceScore * 100}%`}} title={`${(selectedEdge.confidenceScore * 100).toFixed(0)}%`}/>
                    </div>
                  </div>
                )}
              </div>
              {selectedEdge.clinicalRelevance !== undefined && (
                 <div className="mt-2">
                  <label className="text-xs text-gray-400">CrisPRO.ai Clinical Relevance:</label>
                  <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                    <div className="bg-red-500 h-2.5 rounded-full" style={{ width: `${selectedEdge.clinicalRelevance * 100}%`}} title={`${(selectedEdge.clinicalRelevance * 100).toFixed(0)}%`}/>
                  </div>
                </div>
              )}
              {selectedEdge.literatureSupport !== undefined && (
                 <div className="mt-2">
                  <label className="text-xs text-gray-400">CrisPRO.ai Literature Support:</label>
                  <div className="w-full bg-slate-700 h-2.5 rounded-full mt-1">
                    <div className="bg-yellow-500 h-2.5 rounded-full" style={{ width: `${selectedEdge.literatureSupport * 100}%`}} title={`${(selectedEdge.literatureSupport * 100).toFixed(0)}%`}/>
                  </div>
                </div>
              )}

              {selectedEdge.evidenceSourceType && (
                <p className="text-sm mt-2"><span className="font-semibold text-gray-300">Evidence Source:</span> {selectedEdge.evidenceSourceType.replace('_', ' ')}</p>
              )}
              {selectedEdge.mechanismDetails && (
                <p className="text-sm mt-2"><span className="font-semibold text-gray-300">Mechanism Details:</span> {selectedEdge.mechanismDetails}</p>
              )}
               {selectedEdge.inTherapeuticPath && enableAIFeatures.highlightTherapeuticPathways && (
                <p className="text-sm mt-2 px-2 py-1 bg-sky-700 rounded text-sky-200 w-fit">Part of Highlighted Therapeutic Pathway</p>
              )}

              {enableAIFeatures.enableLLMInsights && (
                <div className="mt-3 pt-3 border-t border-slate-700">
                  <h4 className="text-sm font-semibold mb-1 text-sky-400">CrisPRO.ai Insight (Edge):</h4>
                  {isLoadingLLMInsight ? (
                    <p className="text-sm italic text-sky-400 animate-pulse">CrisPRO processing edge insight...</p>
                  ) : (selectedEdge as CrisPROGraphEdge).insightSummary ? ( // Assuming insightSummary can be on edges too
                    <p className="text-sm italic text-gray-300 whitespace-pre-line">{(selectedEdge as CrisPROGraphEdge).insightSummary}</p>
                  ) : (
                    <p className="text-sm italic text-gray-500">Click edge to attempt LLM insight generation.</p>
                  )}
                </div>
              )}


            </div>
          )}

          {!selectedNode && !selectedEdge &&(
            <div className="text-gray-500 text-sm italic">
              Select a node or edge from the graph to view its details and AI-driven insights.
              CrisPRO.ai enhances this data with relevance scores, predictive impacts, and contextual summaries.
            </div>
          )}
        </aside>
      </div>
    </div>
  );
} 