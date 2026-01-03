# Universal Data Pipeline Architecture

## Vision: One Pipeline, Infinite Use Cases

The goal is to transform the myeloma blueprint into a **universal data extraction and analysis pipeline** that can be easily configured for any domain (cancer types, diseases, biological processes) without code changes.

## 1) Current Blueprint Analysis

### What Works Well:
- **Modular Components**: Frontend components are reusable
- **API Orchestration**: Clean proxy pattern to backend services  
- **Provenance Tracking**: Excellent auditability
- **Caching Strategy**: Performance optimization

### What Needs Generalization:
- **Data Sources**: Hardcoded myeloma variants → configurable data extractors
- **Analysis Logic**: Myeloma-specific pathway scoring → domain-agnostic analysis workflows
- **Interpretation Rules**: "Likely Resistant/Sensitive" → configurable decision engines
- **Visualization**: Myeloma-specific charts → generic visualization framework

## 2) Universal Pipeline Architecture

### 2.1 Data Source Layer (Extractors)

#### API Extractor Types:
```python
class DataExtractor:
    def __init__(self, config):
        self.source_type = config['type']  # 'api', 'database', 'file', 'web_scraper'
        self.endpoint = config['endpoint']
        self.auth = config.get('auth', {})
        self.transformers = config.get('transformers', [])
    
    async def extract(self, query_params):
        # Generic extraction logic
        pass

# Specific Implementations
class ClinVarExtractor(DataExtractor):
    """Extract variants from ClinVar API"""
    pass

class EnsemblExtractor(DataExtractor):
    """Extract gene annotations from Ensembl"""
    pass

class PubMedExtractor(DataExtractor):
    """Extract literature data from PubMed"""
    pass

class CosmicExtractor(DataExtractor):
    """Extract cancer variants from COSMIC"""
    pass
```

#### Configuration-Driven Extraction:
```yaml
# config/data_sources.yaml
breast_cancer_panel:
  type: api
  endpoint: https://clinicaltables.nlm.nih.gov/api/ncbi_genes/v3/search
  params:
    q: breast+cancer+gene
    maxList: 100
  transformers:
    - type: filter
      criteria: {evidence_level: high}
    - type: map
      mapping: {gene_symbol: symbol, variant_info: hgvs}

myeloma_variants:
  type: file
  path: data/myeloma_panel.json
  format: json
  schema:
    gene: string
    variant: string
    evidence: string

sarcoma_fusions:
  type: api
  endpoint: https://mitelmandatabase.isb-cgc.org/api/v1/fusions
  params:
    cancer_type: sarcoma
  transformers:
    - type: deduplicate
      key: fusion_genes
```

### 2.2 Analysis Workflow Layer

#### Domain-Agnostic Analysis Engine:
```python
class AnalysisWorkflow:
    def __init__(self, domain_config):
        self.domain = domain_config['domain']
        self.steps = self.build_steps(domain_config['workflow'])
        self.decision_engine = self.build_decision_engine(domain_config['rules'])
    
    def build_steps(self, workflow_config):
        steps = []
        for step_config in workflow_config:
            if step_config['type'] == 'variant_scoring':
                steps.append(VariantScoringStep(step_config))
            elif step_config['type'] == 'pathway_analysis':
                steps.append(PathwayAnalysisStep(step_config))
            elif step_config['type'] == 'literature_mining':
                steps.append(LiteratureMiningStep(step_config))
        return steps
    
    def execute(self, input_data):
        results = {}
        for step in self.steps:
            results[step.name] = step.execute(input_data, results)
        return self.decision_engine.evaluate(results)

# Domain-Specific Workflows
class CancerAnalysisWorkflow(AnalysisWorkflow):
    def __init__(self):
        super().__init__({
            'domain': 'cancer',
            'workflow': [
                {'type': 'variant_scoring', 'endpoint': '/predict_variant_impact'},
                {'type': 'pathway_analysis', 'pathways': ['RAS/MAPK', 'PI3K/AKT', 'TP53']},
                {'type': 'drug_response_prediction', 'method': 'pathway_activity'}
            ],
            'rules': {
                'resistance_threshold': 2.0,
                'pathway_weights': {'RAS/MAPK': 0.7, 'TP53': 0.3}
            }
        })
```

### 2.3 Interpretation & Decision Engine

#### Configurable Decision Rules:
```python
class DecisionEngine:
    def __init__(self, rules_config):
        self.rules = self.parse_rules(rules_config)
        self.templates = rules_config.get('output_templates', {})
    
    def parse_rules(self, config):
        rules = []
        for rule_config in config['rules']:
            if rule_config['type'] == 'threshold':
                rules.append(ThresholdRule(rule_config))
            elif rule_config['type'] == 'pathway_sum':
                rules.append(PathwaySumRule(rule_config))
            elif rule_config['type'] == 'machine_learning':
                rules.append(MLRule(rule_config))
        return rules
    
    def evaluate(self, analysis_results):
        decisions = {}
        for rule in self.rules:
            decisions[rule.name] = rule.evaluate(analysis_results)
        
        # Combine decisions
        final_decision = self.combine_decisions(decisions)
        return {
            'verdict': final_decision,
            'evidence': decisions,
            'confidence': self.calculate_confidence(decisions),
            'rationale': self.generate_rationale(decisions, self.templates)
        }

# Rule Examples
class PathwaySumRule:
    def __init__(self, config):
        self.name = config['name']
        self.pathways = config['pathways']
        self.weights = config['weights']
        self.threshold = config['threshold']
    
    def evaluate(self, results):
        total_score = 0
        for pathway, weight in self.weights.items():
            pathway_score = results.get(pathway, {}).get('activity', 0)
            total_score += pathway_score * weight
        
        return {
            'score': total_score,
            'threshold': self.threshold,
            'meets_criteria': total_score >= self.threshold
        }
```

### 2.4 Visualization Framework

#### Generic Visualization Components:
```python
class VisualizationFactory:
    @staticmethod
    def create_chart(chart_type, data, config):
        if chart_type == 'profile':
            return DeltaProfileChart(data, config)
        elif chart_type == 'heatmap':
            return PathwayHeatmap(data, config)
        elif chart_type == 'network':
            return InteractionNetwork(data, config)
        elif chart_type == 'timeline':
            return EvidenceTimeline(data, config)

class AdaptableDashboard:
    def __init__(self, domain_config):
        self.components = self.build_components(domain_config['dashboard'])
    
    def build_components(self, dashboard_config):
        components = []
        for comp_config in dashboard_config:
            component = {
                'type': comp_config['type'],
                'title': comp_config['title'],
                'data_source': comp_config['data_source'],
                'visualization': VisualizationFactory.create_chart(
                    comp_config['chart_type'],
                    data_source=comp_config['data_source'],
                    config=comp_config.get('config', {})
                )
            }
            components.append(component)
        return components
```

## 3) Implementation Strategy

### 3.1 Configuration-Driven Architecture

#### Domain Configuration Template:
```yaml
# config/domains/breast_cancer.yaml
domain: breast_cancer
name: "Hereditary Breast Cancer Risk Assessment"
description: "AI-powered analysis of hereditary breast cancer variants"

data_sources:
  - name: hereditary_genes
    type: api
    endpoint: https://clinicaltables.nlm.nih.gov/api/ncbi_genes/v3/search
    params:
      q: "hereditary breast cancer gene"
    transformers:
      - type: filter
        criteria: {evidence_level: "high"}
      - type: map
        mapping: {gene_symbol: symbol}

  - name: clinvar_variants
    type: api
    endpoint: https://api.ncbi.nlm.nih.gov/variation/v0/
    params:
      gene: "{gene_symbol}"
    auth:
      type: api_key
      key: NCBI_API_KEY

analysis_workflow:
  - name: variant_impact
    type: evo2_endpoint
    endpoint: /predict_variant_impact
    inputs:
      - source: clinvar_variants
        field: hgvs
    outputs:
      delta_score: result.delta_likelihood_score

  - name: pathway_analysis
    type: custom_analysis
    script: scripts/pathway_analysis.py
    inputs:
      - source: variant_impact
        field: delta_score
    outputs:
      pathway_activity: result.pathway_scores

decision_rules:
  - name: risk_assessment
    type: pathway_sum
    pathways: [BRCA1, BRCA2, TP53, CHEK2]
    weights: {BRCA1: 0.4, BRCA2: 0.4, TP53: 0.1, CHEK2: 0.1}
    threshold: 2.0
    output_template: |
      Risk Level: {risk_category}
      Evidence: {pathway_scores}
      Confidence: {confidence_score}

dashboard:
  - type: summary_card
    title: "Risk Assessment"
    data_source: decision_rules.risk_assessment
    chart_type: gauge

  - type: detail_panel
    title: "Variant Analysis"
    data_source: analysis_workflow.variant_impact
    chart_type: profile

  - type: evidence_panel
    title: "Supporting Evidence"
    data_source: analysis_workflow.variant_impact
    chart_type: timeline
```

### 3.2 Extensible API Extractors

#### Universal API Client:
```python
class UniversalApiClient:
    def __init__(self):
        self.extractors = {}
        self.auth_handlers = {}
        self.transformers = {}
    
    def register_extractor(self, name, extractor_class):
        self.extractors[name] = extractor_class
    
    def register_auth_handler(self, auth_type, handler_class):
        self.auth_handlers[auth_type] = handler_class
    
    def register_transformer(self, transformer_type, transformer_class):
        self.transformers[transformer_type] = transformer_class
    
    async def extract(self, config, query_params=None):
        # Create appropriate extractor
        extractor_class = self.extractors[config['type']]
        extractor = extractor_class(config)
        
        # Apply authentication if needed
        if 'auth' in config:
            auth_handler = self.auth_handlers[config['auth']['type']](config['auth'])
            extractor.set_auth_handler(auth_handler)
        
        # Extract data
        raw_data = await extractor.extract(query_params or config.get('params', {}))
        
        # Apply transformers
        transformed_data = raw_data
        for transformer_config in config.get('transformers', []):
            transformer_class = self.transformers[transformer_config['type']]
            transformer = transformer_class(transformer_config)
            transformed_data = transformer.apply(transformed_data)
        
        return transformed_data

# Usage
client = UniversalApiClient()
await client.extract({
    'type': 'clinvar_api',
    'endpoint': 'https://api.ncbi.nlm.nih.gov/variation/v0/',
    'auth': {'type': 'api_key', 'key': 'NCBI_API_KEY'},
    'params': {'gene': 'BRCA1'},
    'transformers': [
        {'type': 'filter', 'criteria': {'clinical_significance': 'pathogenic'}},
        {'type': 'map', 'mapping': {'variant_name': 'hgvs_p'}}
    ]
})
```

### 3.3 Plugin Architecture for New Domains

#### Domain Plugin Interface:
```python
class DomainPlugin:
    def __init__(self, domain_config):
        self.config = domain_config
        self.name = domain_config['name']
        self.description = domain_config['description']
    
    async def initialize(self):
        """Setup domain-specific resources"""
        pass
    
    def get_data_sources(self):
        """Return domain-specific data source configs"""
        return self.config.get('data_sources', [])
    
    def get_analysis_workflow(self):
        """Return domain-specific analysis steps"""
        return self.config.get('analysis_workflow', [])
    
    def get_decision_rules(self):
        """Return domain-specific decision logic"""
        return self.config.get('decision_rules', [])
    
    def get_dashboard_config(self):
        """Return domain-specific dashboard layout"""
        return self.config.get('dashboard', [])
    
    def get_custom_components(self):
        """Return any domain-specific UI components"""
        return self.config.get('custom_components', [])

# Plugin Registry
class DomainPluginRegistry:
    def __init__(self):
        self.plugins = {}
        self.active_plugins = {}
    
    def register_plugin(self, domain_name, plugin_class, config_path):
        self.plugins[domain_name] = {
            'class': plugin_class,
            'config_path': config_path
        }
    
    async def load_plugin(self, domain_name):
        if domain_name not in self.plugins:
            raise ValueError(f"Plugin {domain_name} not registered")
        
        plugin_info = self.plugins[domain_name]
        config = await self.load_config(plugin_info['config_path'])
        plugin = plugin_info['class'](config)
        await plugin.initialize()
        
        self.active_plugins[domain_name] = plugin
        return plugin
    
    async def load_config(self, config_path):
        # Load YAML/JSON config
        pass
```

## 4) Data Source Capabilities Brainstorm

### 4.1 Biological Databases & APIs

#### Variant Databases:
- **ClinVar**: Pathogenic/likely pathogenic variants
- **COSMIC**: Cancer-specific variants and fusions
- **dbSNP**: Common variants and population frequencies
- **gnomAD**: Population-level variant frequencies
- **DECIPHER**: Rare disease variants

#### Gene & Protein Databases:
- **Ensembl**: Gene annotations, transcripts, regulatory elements
- **UniProt**: Protein sequences and functional annotations
- **HGNC**: Gene nomenclature and aliases
- **Gene Ontology**: Functional annotations
- **STRING**: Protein-protein interactions

#### Clinical & Literature:
- **PubMed**: Literature mining for gene-disease associations
- **ClinicalTrials.gov**: Ongoing clinical trials by gene/target
- **DrugBank**: Drug-target interactions
- **DisGeNET**: Gene-disease associations

### 4.2 Cancer-Specific Resources

#### Cancer Databases:
- **TCGA**: The Cancer Genome Atlas - comprehensive cancer genomics
- **ICGC**: International Cancer Genome Consortium
- **cBioPortal**: Cancer genomics data portal
- **TARGET**: Pediatric cancer genomics
- **PCAWG**: Pan-cancer analysis of whole genomes

#### Cancer Types:
- **Breast Cancer**: BRCA1/2, HER2, hormone receptors
- **Lung Cancer**: EGFR, ALK, KRAS, TP53
- **Colorectal Cancer**: APC, KRAS, TP53, microsatellite instability
- **Melanoma**: BRAF, NRAS, KIT
- **Myeloma**: As currently implemented

### 4.3 External API Integration Patterns

#### Authentication & Rate Limiting:
```python
class ApiRateLimiter:
    def __init__(self, requests_per_minute=60):
        self.requests_per_minute = requests_per_minute
        self.request_times = []
    
    async def wait_if_needed(self):
        now = time.time()
        # Remove requests older than 1 minute
        self.request_times = [t for t in self.request_times if now - t < 60]
        
        if len(self.request_times) >= self.requests_per_minute:
            oldest = min(self.request_times)
            wait_time = 60 - (now - oldest)
            if wait_time > 0:
                await asyncio.sleep(wait_time)
        
        self.request_times.append(now)

class AuthenticatedApiClient:
    def __init__(self, auth_config):
        self.auth_type = auth_config['type']
        if self.auth_type == 'api_key':
            self.api_key = auth_config['key']
        elif self.auth_type == 'oauth':
            self.oauth_config = auth_config
        elif self.auth_type == 'bearer_token':
            self.token = auth_config['token']
    
    def get_headers(self):
        if self.auth_type == 'api_key':
            return {'Authorization': f'Bearer {self.api_key}'}
        elif self.auth_type == 'bearer_token':
            return {'Authorization': f'Bearer {self.token}'}
        return {}
```

#### Error Handling & Retry Logic:
```python
class ResilientApiClient:
    def __init__(self, max_retries=3, backoff_factor=2):
        self.max_retries = max_retries
        self.backoff_factor = backoff_factor
    
    async def make_request(self, url, method='GET', **kwargs):
        for attempt in range(self.max_retries):
            try:
                response = await self._make_request(url, method, **kwargs)
                response.raise_for_status()
                return response
            except Exception as e:
                if attempt == self.max_retries - 1:
                    raise e
                
                wait_time = self.backoff_factor ** attempt
                await asyncio.sleep(wait_time)
    
    async def _make_request(self, url, method, **kwargs):
        # Actual HTTP request implementation
        pass
```

#### Data Transformation Pipeline:
```python
class DataTransformerPipeline:
    def __init__(self, transformers_config):
        self.transformers = []
        for config in transformers_config:
            transformer = self.create_transformer(config)
            self.transformers.append(transformer)
    
    def create_transformer(self, config):
        transformer_type = config['type']
        if transformer_type == 'filter':
            return FilterTransformer(config['criteria'])
        elif transformer_type == 'map':
            return MapTransformer(config['mapping'])
        elif transformer_type == 'deduplicate':
            return DeduplicateTransformer(config['key'])
        elif transformer_type == 'normalize':
            return NormalizeTransformer(config['fields'])
        elif transformer_type == 'join':
            return JoinTransformer(config['join_key'], config['source'])
    
    async def transform(self, data):
        result = data
        for transformer in self.transformers:
            result = await transformer.apply(result)
        return result
```

## 5) Implementation Roadmap

### Phase 1: Core Pipeline Infrastructure (2 weeks)
- [ ] Implement UniversalApiClient with authentication
- [ ] Create DataTransformerPipeline framework
- [ ] Build DomainPlugin interface and registry
- [ ] Set up configuration management system

### Phase 2: Data Source Integration (3 weeks)
- [ ] Implement ClinVar, Ensembl, COSMIC extractors
- [ ] Add authentication handlers (API key, OAuth, Bearer)
- [ ] Create rate limiting and retry logic
- [ ] Build data validation and schema enforcement

### Phase 3: Analysis Workflow Engine (3 weeks)
- [ ] Implement AnalysisWorkflow base class
- [ ] Create domain-agnostic analysis steps
- [ ] Build configurable DecisionEngine
- [ ] Add support for custom analysis scripts

### Phase 4: Visualization Framework (2 weeks)
- [ ] Create generic visualization components
- [ ] Implement AdaptableDashboard system
- [ ] Add chart type registry
- [ ] Support custom visualization plugins

### Phase 5: Domain Templates (2 weeks)
- [ ] Create breast cancer domain template
- [ ] Implement lung cancer configuration
- [ ] Build colorectal cancer setup
- [ ] Add neurological disorder templates

### Phase 6: Testing & Validation (2 weeks)
- [ ] Create comprehensive test suite
- [ ] Validate against known datasets
- [ ] Performance optimization
- [ ] Documentation and examples

## 6) Benefits of This Architecture

### 6.1 For Developers:
- **Rapid Prototyping**: New domains in days, not months
- **Code Reuse**: 90%+ of infrastructure is shared
- **Maintainability**: Centralized configuration management
- **Extensibility**: Plugin architecture for custom logic

### 6.2 For Domain Experts:
- **Self-Service**: Configure new analyses without coding
- **Iterative Development**: Quick feedback loops
- **Consistency**: Standardized workflows across domains
- **Flexibility**: Domain-specific customizations when needed

### 6.3 For Business:
- **Faster Time-to-Market**: New use cases rapidly deployable
- **Scalability**: Handle multiple concurrent domains
- **Cost Efficiency**: Shared infrastructure reduces overhead
- **Innovation**: Easy experimentation with new data sources

## 7) Example: Instantiating a New Domain

```python
# 1. Create domain configuration
config = {
    'domain': 'lung_cancer',
    'name': 'Lung Cancer EGFR Analysis',
    'data_sources': [
        {
            'name': 'egfr_variants',
            'type': 'clinvar_api',
            'endpoint': 'https://api.ncbi.nlm.nih.gov/variation/v0/',
            'params': {'gene': 'EGFR'},
            'transformers': [
                {'type': 'filter', 'criteria': {'clinical_significance': 'pathogenic'}}
            ]
        }
    ],
    'analysis_workflow': [
        {'type': 'variant_scoring', 'endpoint': '/predict_variant_impact'},
        {'type': 'drug_response', 'drugs': ['osimertinib', 'erlotinib']}
    ],
    'decision_rules': [
        {
            'name': 'egfr_sensitivity',
            'type': 'threshold',
            'metric': 'variant_impact',
            'threshold': -0.5,
            'output_template': 'Predicted {drug} sensitivity: {sensitivity_level}'
        }
    ]
}

# 2. Register and load
registry = DomainPluginRegistry()
registry.register_plugin('lung_cancer', CancerAnalysisPlugin, config)
plugin = await registry.load_plugin('lung_cancer')

# 3. Use immediately
results = await plugin.analyze_patient_data(patient_variants)
```

This architecture transforms domain expertise into configuration, enabling rapid deployment of specialized analysis pipelines while maintaining the robustness and features of the original myeloma blueprint.
