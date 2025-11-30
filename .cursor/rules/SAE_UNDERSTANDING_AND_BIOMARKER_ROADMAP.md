# 🧬 SAE UNDERSTANDING & BIOMARKER DISCOVERY ROADMAP

**Date**: January 14, 2025  
**Status**: ✅ **CURRENT STATE ASSESSED** | ✅ **EVO2 SAE NOTEBOOK PATTERN LEARNED** | 🚧 **BIOMARKER DISCOVERY ROADMAP**

---

## 📊 **WHAT I'VE LEARNED ABOUT SAE**

### **1. SAE Theory (From Evo2 Paper)** ✅

**Architecture**:
- **Type**: Batch-TopK Sparse Autoencoder
- **Target Layer**: Layer 26 (Hyena-MR block) - most features of interest
- **Dimensionality**: 
  - Input: `d_model = 4,096` (Evo2 activations)
  - Output: `d_feature = 32,768` (8x overcomplete representation)
  - Batch-TopK: `k = 64` (zeros out all but k largest elements per batch)
- **Training Data**: 1 billion tokens (eukaryotic + prokaryotic genomes)
- **Loss Function**: `L = L_recon + α * L_aux` (reconstruction + auxiliary for dead features)

**Revealed Features** (from paper):
- Exon-intron boundaries
- Transcription factor binding motifs
- Protein structural elements (alpha-helices, beta-sheets)
- Prophage regions & mobile genetic elements
- CRISPR spacer sequences
- Mutation severity signatures (frameshift/stop)

**Key Insight**: SAEs decompose Evo2's black box into 32,768 interpretable biological features that Evo2 learned autonomously.

---

### **2. Current Platform Implementation** ✅

**What We Have NOW** (`sae_feature_service.py`):

**CRITICAL GAP**: We're NOT using actual Evo2 SAE activations. We're computing **proxy SAE features** from:

1. **DNA Repair Capacity** (formula):
   ```python
   dna_repair_capacity = (
       0.6 * pathway_burden_ddr +           # From pathway aggregation
       0.2 * essentiality_hrr_genes +      # From insights bundle
       0.2 * exon_disruption_score          # From insights bundle
   )
   ```

2. **Mechanism Vector** (7D):
   ```python
   mechanism_vector = [
       pathway_burden_ddr,      # From pathway aggregation
       pathway_burden_mapk,      # From pathway aggregation
       pathway_burden_pi3k,      # From pathway aggregation
       pathway_burden_vegf,      # From pathway aggregation
       pathway_burden_her2,      # From pathway aggregation
       1.0 if io_eligible else 0.0,  # From tumor context (TMB/MSI)
       cross_resistance_risk     # From treatment history
   ]
   ```

3. **Pathway Burdens**: Computed from S/P/E pathway aggregation, NOT from SAE features

4. **Resistance Detection**: Uses DNA repair capacity trends, NOT SAE feature activations

**Data Sources** (current):
- ✅ Insights Bundle (functionality, chromatin, essentiality, regulatory)
- ✅ Pathway Scores (from S/P/E framework)
- ✅ Tumor Context (HRD, TMB, MSI, somatic mutations)
- ✅ Treatment History
- ✅ CA-125 Intelligence

**What We DON'T Have**:
- ❌ Access to Evo2 layer 26 activations
- ❌ Trained SAE model (32K features)
- ❌ SAE feature extraction from actual Evo2 activations
- ❌ Correlation analysis infrastructure for biomarker discovery

---

### **3. Evo2 Model Capabilities** ✅

**Good News**: The Evo2 model class (`src/services/genesis_engine/evo2/evo2/models.py`) has:

```python
def forward(
    self,
    input_ids: torch.Tensor,
    return_embeddings: bool = False,
    layer_names=None,
) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
    """Forward pass with optional embedding extraction."""
```

**This means**: We CAN extract layer 26 activations if we:
1. Call `forward()` with `return_embeddings=True`
2. Specify `layer_names=["layer_26"]` (or correct layer name)
3. Get the activations tensor

**BUT**: We still need the trained SAE model to decode these activations into 32K features.

---

## 🎯 **HOW TO ACHIEVE BIOMARKER DISCOVERY**

### **Phase 1: Access Evo2 SAE Features** (Foundation)

#### **Step 1.1: Extract Layer 26 Activations** (2-3 days)

**What**: Modify Evo2 service to return layer 26 activations

**How**:
1. Update `src/services/evo_service/main.py`:
   ```python
   @app.function(...)
   def score_variant_with_activations(
       chrom: str,
       pos: int,
       ref: str,
       alt: str,
       return_activations: bool = False
   ):
       # Get sequence context (8kb window)
       sequence = fetch_sequence(chrom, pos, window=8192)
       
       # Tokenize
       input_ids = evo2_model.tokenizer.encode(sequence)
       
       # Forward pass with activations
       logits, embeddings = evo2_model.forward(
           input_ids,
           return_embeddings=True,
           layer_names=["layer_26"]  # Or correct layer name
       )
       
       # Return activations if requested
       if return_activations:
           return {
               "delta_score": compute_delta(logits, ref, alt),
               "layer_26_activations": embeddings["layer_26"].cpu().numpy()
           }
   ```

2. **New Endpoint**: `POST /api/evo/score_variant_with_activations`

**Acceptance**:
- ✅ Returns layer 26 activations (shape: [batch, seq_len, 4096])
- ✅ Works for all Evo2 models (1B, 7B, 40B)
- ✅ Cached for performance

---

#### **Step 1.2: Load Trained SAE Model** (1-2 weeks)

**What**: Obtain and load a **compatible Evo2 SAE model** using the official Evo2 sparse autoencoder notebook pattern.

**What I Learned (Jan 14, 2025)**:
- The Evo2 team published a mechanistic interpretability notebook at  
  `scripts/evo2/evo2/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb` (now pulled locally).
- The notebook defines three key classes:
  - `ModelScope`: utility for adding/removing PyTorch hooks cleanly.
  - `ObservableEvo2`: wraps an `Evo2` model and exposes `forward`/`generate` with activation caching for specified layers (e.g. `"blocks-26"`).
  - `BatchTopKTiedSAE`: the sparse autoencoder (tied decoder, batch top‑k, 32,768 features) used in the Evo2 paper.
- The notebook downloads SAE weights from Hugging Face (e.g. `Goodfire/Evo-2-Layer-26-Mixed`) and then:
  1. Loads Evo2 (`ObservableEvo2("evo2_7b_262k")`),
  2. Loads `BatchTopKTiedSAE` with the HF checkpoint,
  3. Extracts layer‑26 activations via `forward(..., cache_activations_at=["blocks-26"])`,
  4. Passes those activations through the SAE to get sparse feature activations.

**How (Applied to Our Stack)**:
1. **Download SAE Model Weights (HF)**:
   - Use `hf_hub_download` inside our Modal image, following the notebook’s pattern.
   - Pin the exact repo + filename (e.g. `Goodfire/Evo-2-Layer-26-Mixed`, `sae_weights.pt`) in config for reproducibility.

2. **Create SAE Service** (`api/services/sae_model_service.py`):
   ```python
   class SAEModelService:
       def __init__(self):
           self.sae_model = load_sae_checkpoint("evo2_sae_layer26.pt")
           self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
           self.sae_model.to(self.device)
       
       def extract_features(
           self,
           activations: torch.Tensor  # Shape: [batch, seq_len, 4096]
       ) -> Dict[str, Any]:
           """
           Extract 32K SAE features from layer 26 activations.
           
           Returns:
               {
                   "features": [f0, f1, ..., f32767],  # Activation levels
                   "top_features": [(idx, activation), ...],  # Top k=64
                   "feature_interpretations": {...}  # Known biological meanings
               }
           """
           # Forward through SAE encoder
           with torch.no_grad():
               features = self.sae_model.encoder(activations)
               # Apply Batch-TopK (k=64)
               top_features = self._batch_topk(features, k=64)
           
           return {
               "features": features.cpu().numpy(),
               "top_features": top_features,
               "feature_interpretations": self._interpret_features(top_features)
           }
   ```

3. **Deploy on Modal** (GPU required):
   - Create `src/services/sae_service/main.py`
   - Mount SAE model weights
   - Expose endpoint: `POST /api/sae/extract_features`

**Acceptance**:
- ✅ Loads SAE model successfully
- ✅ Extracts 32K features from layer 26 activations
- ✅ Returns top k=64 features per sequence
- ✅ Performance: <100ms per variant (GPU)

---

#### **Step 1.3: Integrate SAE Extraction into Pipeline** (3-5 days)

**What**: Wire SAE extraction into existing efficacy/resistance workflows

**How**:
1. **Update `sae_feature_service.py`**:
   ```python
   def compute_sae_features(
       self,
       insights_bundle: Dict[str, Any],
       pathway_scores: Dict[str, float],
       tumor_context: Dict[str, Any],
       # NEW: Add actual SAE features
       evo2_activations: Optional[torch.Tensor] = None,
       sae_features: Optional[Dict[str, Any]] = None
   ) -> SAEFeatures:
       # Existing proxy features (keep for backward compatibility)
       dna_repair_capacity = self._compute_dna_repair_capacity(...)
       
       # NEW: Actual SAE features
       if sae_features:
           # Extract known biological features from SAE
           ddr_sae_features = self._extract_ddr_features(sae_features)
           mapk_sae_features = self._extract_mapk_features(sae_features)
           # ... etc
           
           # Enhance proxy features with SAE signals
           dna_repair_capacity = self._enhance_with_sae(
               dna_repair_capacity,
               ddr_sae_features
           )
   ```

2. **Update Evo2 Scorer**:
   - Add `return_activations=True` flag
   - Call SAE service if flag enabled
   - Pass SAE features to `sae_feature_service`

**Acceptance**:
- ✅ SAE features enhance existing proxy features
- ✅ Backward compatible (works without SAE)
- ✅ Provenance tracks SAE vs proxy features

---

### **Phase 2: Biomarker Discovery Infrastructure** (3-4 weeks)

#### **Step 2.1: Large-Scale Patient Dataset** (1-2 weeks)

**What**: Collect patient datasets with outcomes

**Data Requirements**:
- **Genomic**: WGS/WES/targeted panels (somatic + germline)
- **Clinical**: Cancer type, stage, treatment history, response, PFS, OS
- **Size**: 1,000+ patients per cancer type (minimum for statistical power)

**Sources**:
- **TCGA**: Already have ovarian cancer (469 patients with platinum response)
- **cBioPortal**: Extract additional cohorts
- **GDC**: Download public datasets
- **Partnerships**: Clinical trial data (Yale, etc.)

**Infrastructure**:
```python
# New service: api/services/biomarker_dataset_service.py
class BiomarkerDatasetService:
    def load_patient_cohort(
        self,
        study_id: str,
        cancer_type: str
    ) -> List[PatientRecord]:
        """
        Load patient cohort with:
        - Genomic data (variants, CNAs, etc.)
        - Clinical outcomes (response, PFS, OS)
        - Treatment history
        """
        # Extract from cBioPortal/GDC
        # Normalize to common schema
        # Return PatientRecord objects
    ```
```

---

#### **Step 2.2: SAE Feature Extraction Pipeline** (1 week)

**What**: Extract SAE features for all patients in cohort

**How**:
```python
# New script: scripts/extract_sae_features_cohort.py
async def extract_cohort_sae_features(
    cohort: List[PatientRecord],
    evo2_model: str = "evo2_7b"
):
    """
    Extract SAE features for entire cohort.
    
    Process:
    1. For each patient:
       a. Get all somatic mutations
       b. For each mutation:
          - Get 8kb sequence context
          - Call Evo2 with return_activations=True
          - Extract layer 26 activations
          - Call SAE service to get 32K features
          - Store top k=64 features
       c. Aggregate features across all mutations
    2. Store results in database (Supabase/PostgreSQL)
    """
    results = []
    for patient in cohort:
        patient_sae_features = []
        for mutation in patient.somatic_mutations:
            # Call Evo2 + SAE
            sae_result = await extract_sae_features(mutation)
            patient_sae_features.append(sae_result)
        
        # Aggregate (sum, max, mean across mutations)
        aggregated = aggregate_features(patient_sae_features)
        results.append({
            "patient_id": patient.id,
            "sae_features": aggregated,
            "outcomes": patient.outcomes
        })
    
    # Store in database
    store_cohort_sae_features(results)
```

**Performance**:
- **Per-patient**: ~10-50 mutations × 100ms = 1-5 seconds
- **1,000 patients**: ~17-83 minutes (parallelize to 10 workers = 2-8 minutes)

---

#### **Step 2.3: Correlation Analysis Infrastructure** (1 week)

**What**: Statistical modeling to find SAE features correlated with outcomes

**How**:
```python
# New service: api/services/biomarker_correlation_service.py
class BiomarkerCorrelationService:
    def discover_biomarkers(
        self,
        cohort_sae_features: List[Dict],
        outcome_type: str,  # "response", "resistance", "pfs", "os"
        min_correlation: float = 0.3,
        p_value_threshold: float = 0.05
    ) -> List[Biomarker]:
        """
        Discover SAE features correlated with clinical outcomes.
        
        Process:
        1. For each of 32K SAE features:
           a. Extract feature activation levels across cohort
           b. Correlate with outcome (response/resistance/PFS/OS)
           c. Compute correlation coefficient, p-value
           d. Rank by significance
        2. Return top N biomarkers
        """
        biomarkers = []
        for feature_idx in range(32768):
            # Extract feature activations
            activations = [
                patient["sae_features"][feature_idx]
                for patient in cohort_sae_features
            ]
            
            # Extract outcomes
            outcomes = [
                patient["outcomes"][outcome_type]
                for patient in cohort_sae_features
            ]
            
            # Statistical test
            correlation, p_value = pearsonr(activations, outcomes)
            
            if abs(correlation) >= min_correlation and p_value < p_value_threshold:
                biomarkers.append(Biomarker(
                    feature_idx=feature_idx,
                    correlation=correlation,
                    p_value=p_value,
                    interpretation=self._interpret_feature(feature_idx)
                ))
        
        # Rank by absolute correlation
        return sorted(biomarkers, key=lambda x: abs(x.correlation), reverse=True)
```

**Statistical Methods**:
- **Pearson Correlation**: For continuous outcomes (PFS, OS)
- **Chi-square Test**: For binary outcomes (response vs resistance)
- **Cox Proportional Hazards**: For survival analysis (PFS, OS)
- **Logistic Regression**: For multi-class outcomes

**Output**: Top 100-500 SAE features correlated with outcomes (ranked by significance)

---

### **Phase 3: Validation & Application** (4-6 weeks)

#### **Step 3.1: In Silico Validation** (1-2 weeks)

**What**: Cross-reference discovered biomarkers with known biology

**How**:
```python
def validate_biomarkers_in_silico(
    biomarkers: List[Biomarker]
) -> List[ValidatedBiomarker]:
    """
    Validate discovered biomarkers against:
    - TCGA (known cancer drivers)
    - COSMIC (cancer mutations)
    - ClinVar (pathogenic variants)
    - GWAS (genome-wide associations)
    - Literature (PubMed search)
    """
    validated = []
    for biomarker in biomarkers:
        # Check known associations
        known_associations = check_databases(biomarker.feature_idx)
        
        # Literature search
        literature = search_pubmed(
            query=f"SAE feature {biomarker.feature_idx} cancer"
        )
        
        validated.append(ValidatedBiomarker(
            biomarker=biomarker,
            known_associations=known_associations,
            literature_support=literature,
            validation_status="confirmed" if known_associations else "novel"
        ))
    
    return validated
```

**Acceptance**:
- ✅ 50-70% of discovered biomarkers match known biology
- ✅ 20-30% are novel (potential discoveries)
- ✅ Literature support for top 10 biomarkers

---

#### **Step 3.2: Feature Interpretation** (1 week)

**What**: Understand what biological concepts discovered SAE features represent

**How**:
```python
def interpret_sae_features(
    top_biomarkers: List[Biomarker]
) -> Dict[int, FeatureInterpretation]:
    """
    Interpret SAE features using:
    - Contrastive feature search (from Evo2 paper)
    - Sequence enrichment analysis
    - Known feature database (from Evo2 paper)
    """
    interpretations = {}
    for biomarker in top_biomarkers:
        # Contrastive search: Find sequences that activate this feature
        activating_sequences = find_activating_sequences(
            feature_idx=biomarker.feature_idx,
            threshold=0.5
        )
        
        # Analyze sequences for common patterns
        patterns = analyze_sequence_patterns(activating_sequences)
        
        # Match to known biological concepts
        biological_concept = match_to_known_concepts(patterns)
        
        interpretations[biomarker.feature_idx] = FeatureInterpretation(
            feature_idx=biomarker.feature_idx,
            biological_concept=biological_concept,  # e.g., "DNA repair pathway", "Immune evasion"
            sequence_patterns=patterns,
            confidence=compute_confidence(patterns, biological_concept)
        )
    
    return interpretations
```

**Output**: Human-readable interpretations like:
- "SAE Feature 12,345: DNA repair pathway activation (HRD signature)"
- "SAE Feature 23,456: Immune evasion mechanism (PD-L1 upregulation)"
- "SAE Feature 8,901: Drug metabolism site (CYP2D6 variant)"

---

#### **Step 3.3: Clinical Integration** (2-3 weeks)

**What**: Integrate discovered biomarkers into clinical workflows

**How**:
1. **Update SAE Feature Service**:
   ```python
   def compute_sae_features(
       self,
       ...
       discovered_biomarkers: Optional[List[Biomarker]] = None
   ) -> SAEFeatures:
       # Existing features
       dna_repair_capacity = ...
       
       # NEW: Apply discovered biomarkers
       if discovered_biomarkers:
           biomarker_scores = self._apply_biomarkers(
               sae_features,
               discovered_biomarkers
           )
           
           # Enhance predictions
           resistance_risk = self._enhance_with_biomarkers(
               resistance_risk,
               biomarker_scores
           )
   ```

2. **New Endpoint**: `POST /api/biomarkers/predict`
   - Input: Patient genomic data
   - Output: Biomarker-based risk scores (response, resistance, prognosis)

3. **Frontend Integration**:
   - Display biomarker-based predictions
   - Show feature interpretations
   - Provide confidence scores

**Acceptance**:
- ✅ Biomarker predictions improve accuracy by 5-10% over baseline
- ✅ Interpretations are clinically meaningful
- ✅ Confidence scores are calibrated

---

## 📊 **CURRENT STATE vs. BIOMARKER DISCOVERY**

| Component | Current State | Biomarker Discovery Needs |
|-----------|--------------|--------------------------|
| **SAE Features** | Proxy (computed from insights/pathways) | Actual (from Evo2 layer 26 activations) |
| **SAE Model** | ❌ Not loaded | ✅ Trained model (32K features) |
| **Layer 26 Access** | ✅ Possible (Evo2 model supports it) | ✅ Need to implement |
| **Patient Cohorts** | ✅ TCGA ovarian (469 patients) | ✅ Need 1,000+ per cancer type |
| **Correlation Analysis** | ❌ Not implemented | ✅ Statistical modeling infrastructure |
| **Validation** | ❌ Not implemented | ✅ In silico + in vitro pipelines |
| **Clinical Integration** | ✅ SAE features in workflows | ✅ Biomarker-based predictions |

---

## 🎯 **ROADMAP SUMMARY**

### **Phase 1: Foundation** (3-4 weeks)
- ✅ Extract layer 26 activations from Evo2
- ✅ Load trained SAE model
- ✅ Integrate SAE extraction into pipeline

### **Phase 2: Discovery** (3-4 weeks)
- ✅ Large-scale patient datasets
- ✅ SAE feature extraction pipeline
- ✅ Correlation analysis infrastructure

### **Phase 3: Validation** (4-6 weeks)
- ✅ In silico validation
- ✅ Feature interpretation
- ✅ Clinical integration

**Total Timeline**: 10-14 weeks (2.5-3.5 months)

**Key Dependencies**:
1. **SAE Model Availability**: Need access to official Evo2 SAE model weights
2. **Patient Data**: Need partnerships or public datasets (TCGA, cBioPortal)
3. **Compute Resources**: GPU cluster for large-scale SAE extraction (Modal can handle)

---

## 💡 **KEY INSIGHTS**

1. **Current SAE Implementation is a Proxy**: We're computing features from downstream signals (insights, pathways), NOT from actual Evo2 SAE activations. This is still valuable but not true SAE-based biomarker discovery.

2. **Evo2 Model Supports Activation Extraction**: The `forward()` method can return layer embeddings, so we CAN access layer 26 activations. We just need to implement it.

3. **SAE Model is the Missing Piece**: We need the trained SAE model (32K features) to decode activations into interpretable features. This is likely available from the Evo2 repository.

4. **Biomarker Discovery is a Multi-Stage Process**: Not just extracting features, but correlating with outcomes, validating, and interpreting. This requires significant infrastructure.

5. **Current Proxy Features Are Still Valuable**: DNA repair capacity, mechanism vectors, etc. are clinically useful even without true SAE features. Biomarker discovery would ENHANCE these, not replace them.

---

**Confidence Level**: 85% - I understand the theory and current implementation, but biomarker discovery requires infrastructure we don't have yet. The roadmap is realistic if we have access to the SAE model and patient datasets.

