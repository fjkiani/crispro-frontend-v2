# Existing Classifier Training Pipeline Analysis

**Date**: January 21, 2025  
**Purpose**: Document existing classifier training infrastructure and adapt for BRCA classifier

---

## ✅ Found: Complete Adjudicator Training Pipeline

### Pipeline Overview

**Location**: `src/tools/adjudicator_trainer/`

**Pipeline Steps**:
1. **Phase 1.1**: Parse ClinVar data (`01_parse_variant_summary.py`)
2. **Phase 1.3**: Generate embeddings (`02_generate_embeddings.py`)
3. **Phase 2.1**: Train classifier (`03_train_adjudicator.py`)
4. **Phase 2.2**: Visualize performance (`04_visualize_performance.py`)

**Deployment**: `src/services/adjudicator/main.py` (Modal service)

---

## Phase 1.1: ClinVar Data Extraction

**File**: `src/tools/adjudicator_trainer/01_parse_variant_summary.py`

### What It Does

1. **Downloads ClinVar Data**:
   - Source: `variant_summary.txt.gz` from NCBI ClinVar
   - URL: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`

2. **Target Genes** (Includes BRCA1/BRCA2!):
   ```python
   TARGET_GENES = [
       "TP53", "KRAS", "BRAF", "EGFR", "PIK3CA", 
       "PTEN", "BRCA1", "BRCA2", "ATM"  # ✅ BRCA1/BRCA2 included!
   ]
   ```

3. **Filtering Criteria**:
   - Assembly: GRCh38
   - Type: Single nucleotide variant (SNV)
   - Protein change: Must have `(p.` notation (missense variants)
   - Review status: Must contain "criteria provided"
   - Clinical significance: Pathogenic or Benign (mapped to 1 or 0)

4. **Output**:
   - File: `data/adjudicator_training/clinvar_missense_labels.csv`
   - Columns: `gene`, `hgvsp`, `hgvsc`, `clinvar_id`, `significance_raw`, `significance_mapped`, `review_status`

### Key Insights

✅ **BRCA1/BRCA2 Already Included**: The script already extracts BRCA1/BRCA2 variants!  
✅ **Pathogenicity Labels**: Variants are labeled as 1 (Pathogenic) or 0 (Benign)  
✅ **High-Quality Filtering**: Only variants with "criteria provided" review status  
✅ **Reusable**: Can adapt this script for BRCA-only extraction

### Adaptation for BRCA Classifier

**Option 1: Filter Existing Output** (Fastest):
```python
# After running 01_parse_variant_summary.py, filter for BRCA1/BRCA2
df = pd.read_csv("data/adjudicator_training/clinvar_missense_labels.csv")
brca_df = df[df['gene'].isin(['BRCA1', 'BRCA2'])]
brca_df.to_csv("data/brca_training/clinvar_brca_labels.csv", index=False)
```

**Option 2: Modify Script** (More Control):
- Change `TARGET_GENES = ["BRCA1", "BRCA2"]`
- Run script to extract BRCA-only variants
- Output: `data/brca_training/clinvar_brca_labels.csv`

---

## Phase 1.3: Embedding Generation

**File**: `src/tools/adjudicator_trainer/02_generate_embeddings.py`

### What It Does

1. **Uses Zeta Oracle** (Evo2 40B Model):
   - URL: `https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run/invoke`
   - Action: `"embed"`
   - Input: `{ "gene_symbol": gene, "protein_change": hgvsp }`
   - Output: Embedding vector (8192 dimensions)

2. **Embedding Source**:
   - **Model**: Evo2 40B (via Zeta Oracle)
   - **Layer**: `blocks.40.mlp.l3` (layer 40, not layer 26!)
   - **Method**: Mean pooling across sequence length
   - **Dimensions**: 8192 (hidden size of 40B model)

3. **Process**:
   - Async batch processing (100 concurrent requests)
   - Incremental saving (every 500 variants)
   - Resume capability (skips already processed variants)

4. **Output**:
   - File: `data/adjudicator_training/labeled_embeddings.parquet`
   - Columns: `gene`, `hgvsp`, `hgvsc`, `clinvar_id`, `significance_mapped`, `embedding` (JSON string)

### Key Insights

⚠️ **Different from Plan**: 
- **Plan**: Use Evo2 layer 26 activations (as per paper)
- **Adjudicator**: Used Zeta Oracle layer 40 embeddings (8192 dim)

⚠️ **Zeta Oracle Status**: 
- URL points to `zeta-oracle-v2` (may not be available)
- Need to verify if service is still operational

✅ **Alternative Available**: 
- Evo2 Modal service has `/score_variant_with_activations` endpoint
- Can extract layer 26 activations directly (as per plan)

### Adaptation for BRCA Classifier

**Option 1: Use Evo2 Layer 26 Activations** (Matches Paper):
```python
# Use /api/evo/score_variant_with_activations endpoint
# Extract layer 26 activations (32K dimensions for 1B/7B, or 8K for 40B)
# This matches the paper's approach (layer 26, not layer 40)
```

**Option 2: Use Zeta Oracle** (If Available):
```python
# Use existing Zeta Oracle service
# Layer 40 embeddings (8192 dimensions)
# Faster if service is operational
```

**Recommendation**: Use Evo2 layer 26 activations (matches paper, more control)

---

## Phase 2.1: Classifier Training

**File**: `src/tools/adjudicator_trainer/03_train_adjudicator.py`

### What It Does

1. **Model**: GradientBoostingClassifier (scikit-learn)
   ```python
   clf = GradientBoostingClassifier(
       n_estimators=200,
       learning_rate=0.1,
       max_depth=5,
       random_state=42,
       validation_fraction=0.1,  # Early stopping
       n_iter_no_change=10
   )
   ```

2. **Data Split**: 80/20 train/test (stratified)

3. **Input Format**:
   - X: Embedding vectors (numpy array)
   - y: Pathogenicity labels (0 = Benign, 1 = Pathogenic)

4. **Output**:
   - Model: `models/adjudicator/adjudicator_v1.pkl` (joblib format)
   - Evaluation: Accuracy, classification report, confusion matrix

### Key Insights

✅ **Model Type**: GradientBoostingClassifier (not XGBoost, but similar)  
✅ **Training Process**: Standard scikit-learn workflow  
✅ **Early Stopping**: Built-in to prevent overfitting  
✅ **Reusable**: Can adapt for BRCA classifier

### Adaptation for BRCA Classifier

**Option 1: Use XGBoost** (As Per Plan):
```python
from xgboost import XGBClassifier

clf = XGBClassifier(
    n_estimators=200,
    learning_rate=0.1,
    max_depth=5,
    random_state=42,
    eval_metric='logloss'
)
```

**Option 2: Use GradientBoostingClassifier** (As Adjudicator):
```python
from sklearn.ensemble import GradientBoostingClassifier

clf = GradientBoostingClassifier(
    n_estimators=200,
    learning_rate=0.1,
    max_depth=5,
    random_state=42,
    validation_fraction=0.1,
    n_iter_no_change=10
)
```

**Recommendation**: Use XGBoost (as per plan, matches paper approach)

---

## Deployment: Adjudicator Service

**File**: `src/services/adjudicator/main.py`

### What It Does

1. **Modal Service**: Deploys trained classifier as API
2. **Model Loading**: Loads `adjudicator_v1.pkl` on container startup
3. **Endpoint**: `POST /classify`
   - Input: `embedding: list[float]`
   - Output: `{ "prediction": "Pathogenic"|"Benign", "is_pathogenic": bool, "confidence_score": float }`

### Key Insights

✅ **Deployment Pattern**: Modal service with joblib model loading  
✅ **API Format**: Simple embedding → classification  
✅ **Reusable**: Can adapt for BRCA classifier deployment

### Adaptation for BRCA Classifier

**Option 1: New Modal Service** (As Per Plan):
```python
# Create src/services/evo_service/brca_classifier.py
# Add /brca_classifier endpoint to evo_service/main.py
# Load trained BRCA classifier model
```

**Option 2: Extend Adjudicator** (Reuse Existing):
```python
# Add BRCA-specific endpoint to adjudicator service
# Or create separate BRCA adjudicator service
```

**Recommendation**: Add to Evo2 Modal service (as per plan, keeps Evo2 capabilities together)

---

## Re-Training History

### What We Know

1. **Adjudicator Was Trained**:
   - Used ClinVar variants (including BRCA1/BRCA2)
   - Used Zeta Oracle embeddings (Evo2 40B, layer 40)
   - Trained GradientBoostingClassifier
   - Deployed as Modal service

2. **Training Data**:
   - Source: ClinVar `variant_summary.txt.gz`
   - Genes: TP53, KRAS, BRAF, EGFR, PIK3CA, PTEN, **BRCA1, BRCA2**, ATM
   - Variants: Missense SNVs with protein changes
   - Labels: Pathogenic (1) vs. Benign (0)

3. **No Evidence of Re-Training**:
   - Only one model version found: `adjudicator_v1.pkl`
   - No versioning or re-training scripts found
   - May have been re-trained but not documented

### What We Need for BRCA Classifier

1. **BRCA-Only Data**:
   - Filter existing ClinVar extraction for BRCA1/BRCA2 only
   - Or re-run extraction with `TARGET_GENES = ["BRCA1", "BRCA2"]`
   - Target: 1000+ variants (500 pathogenic, 500 benign)

2. **Evo2 Layer 26 Embeddings**:
   - Use `/api/evo/score_variant_with_activations` endpoint
   - Extract layer 26 activations (not layer 40)
   - Dimensions: 32K for 1B/7B, or 8K for 40B

3. **XGBoost Training**:
   - Train on Evo2 layer 26 embeddings
   - Target: AUROC >0.94 (match paper)
   - Save model: `models/brca_classifier/brca_classifier_v1.pkl`

---

## Updated Plan: Adapting Existing Pipeline

### Phase 1: Data Preparation (Week 1, Day 1-2)

**Step 1: Extract BRCA1/BRCA2 Variants**
```bash
# Option 1: Filter existing output
python -c "
import pandas as pd
df = pd.read_csv('data/adjudicator_training/clinvar_missense_labels.csv')
brca_df = df[df['gene'].isin(['BRCA1', 'BRCA2'])]
print(f'Found {len(brca_df)} BRCA variants')
print(brca_df['significance_mapped'].value_counts())
brca_df.to_csv('data/brca_training/clinvar_brca_labels.csv', index=False)
"

# Option 2: Re-run extraction (if need more variants)
# Modify 01_parse_variant_summary.py: TARGET_GENES = ["BRCA1", "BRCA2"]
# python src/tools/adjudicator_trainer/01_parse_variant_summary.py
```

**Step 2: Generate Evo2 Layer 26 Embeddings**
```python
# Create: src/tools/brca_trainer/02_generate_evo2_embeddings.py
# Use: /api/evo/score_variant_with_activations endpoint
# Extract: layer 26 activations (blocks.26)
# Output: data/brca_training/labeled_embeddings.parquet
```

**Key Differences from Adjudicator**:
- ✅ Use Evo2 layer 26 (not Zeta Oracle layer 40)
- ✅ Use `/api/evo/score_variant_with_activations` (not Zeta Oracle)
- ✅ Focus on BRCA1/BRCA2 only (not all genes)

### Phase 2: Model Training (Week 1, Day 3-4)

**Step 3: Train XGBoost Classifier**
```python
# Create: src/tools/brca_trainer/03_train_brca_classifier.py
# Model: XGBoost (not GradientBoostingClassifier)
# Target: AUROC >0.94 on test set
# Output: models/brca_classifier/brca_classifier_v1.pkl
```

**Key Differences from Adjudicator**:
- ✅ Use XGBoost (not GradientBoostingClassifier)
- ✅ Target AUROC >0.94 (match paper)
- ✅ BRCA-specific model (not general adjudicator)

### Phase 3: Deployment (Week 1, Day 5)

**Step 4: Deploy to Evo2 Modal Service**
```python
# Add to: src/services/evo_service/main.py
# Endpoint: /brca_classifier
# Load: models/brca_classifier/brca_classifier_v1.pkl
# Input: Variant (gene, position, ref, alt)
# Output: Pathogenicity probability (0-1), confidence
```

**Key Differences from Adjudicator**:
- ✅ Add to Evo2 service (not separate service)
- ✅ Input: Variant coordinates (not embedding)
- ✅ Extract embedding on-the-fly (not pre-computed)

---

## Critical Updates to Plan

### Gap 1: ClinVar Data Extraction ✅ **RESOLVED**

**Original Issue**: Need to extract BRCA1/BRCA2 variants  
**Solution**: ✅ **EXISTING SCRIPT ALREADY EXTRACTS BRCA1/BRCA2!**

**Action**:
1. Check if `data/adjudicator_training/clinvar_missense_labels.csv` exists
2. Filter for BRCA1/BRCA2 variants
3. If insufficient, re-run extraction with BRCA-only focus

**Time Saved**: 1-2 hours (script already exists and works)

### Gap 2: Embedding Extraction ✅ **PARTIALLY RESOLVED**

**Original Issue**: Need to extract Evo2 embeddings  
**Solution**: ✅ **ENDPOINT EXISTS** (`/score_variant_with_activations`)

**Action**:
1. Verify endpoint works for 1B model (test with sample variant)
2. Create embedding extraction script (adapt from `02_generate_embeddings.py`)
3. Use layer 26 activations (not layer 40 like adjudicator)

**Time Saved**: 30 minutes (pattern exists, just need to adapt)

### Gap 3: Training Infrastructure ✅ **RESOLVED**

**Original Issue**: Where to train classifier  
**Solution**: ✅ **EXISTING TRAINING SCRIPT EXISTS**

**Action**:
1. Adapt `03_train_adjudicator.py` for BRCA classifier
2. Change model to XGBoost (instead of GradientBoostingClassifier)
3. Use same data split and evaluation approach

**Time Saved**: 1 hour (training script exists and works)

---

## Updated Pre-Development Action Items

### P0: Critical (Before Development)

1. ✅ **Verify ClinVar Data** (15 minutes)
   - Check if `data/adjudicator_training/clinvar_missense_labels.csv` exists
   - Filter for BRCA1/BRCA2 variants
   - Count: Need 1000+ variants (500 pathogenic, 500 benign)
   - If insufficient, re-run extraction

2. ✅ **Verify Evo2 Embeddings Extraction** (30 minutes)
   - Test `/api/evo/score_variant_with_activations` with 1B model
   - Verify layer 26 activations are returned
   - Check activation dimensions (should be 32K for 1B/7B)
   - Create extraction script (adapt from `02_generate_embeddings.py`)

3. ⚠️ **Create Base Scorer Interface** (1 hour)
   - Still needed (not found in existing code)
   - Create `api/services/sequence_scorers/base_scorer.py`

### P1: Important (During Week 1)

4. ✅ **Create BRCA Training Scripts** (2-3 hours)
   - Adapt `01_parse_variant_summary.py` for BRCA-only (or filter existing)
   - Adapt `02_generate_embeddings.py` for Evo2 layer 26 (not Zeta Oracle)
   - Adapt `03_train_adjudicator.py` for XGBoost (not GradientBoosting)

5. ✅ **Model Storage** (15 minutes)
   - Use same pattern: `models/brca_classifier/brca_classifier_v1.pkl`
   - Joblib format (same as adjudicator)

---

## Key Learnings

### What Worked Well (Adjudicator)

1. ✅ **ClinVar Extraction**: Script works well, includes BRCA1/BRCA2
2. ✅ **Training Pipeline**: Standard scikit-learn workflow is solid
3. ✅ **Deployment**: Modal service pattern works well
4. ✅ **Incremental Saving**: Good for large datasets

### What to Change (BRCA Classifier)

1. ⚠️ **Embedding Source**: Use Evo2 layer 26 (not Zeta Oracle layer 40)
2. ⚠️ **Model Type**: Use XGBoost (not GradientBoostingClassifier)
3. ⚠️ **Gene Focus**: BRCA1/BRCA2 only (not all genes)
4. ⚠️ **Deployment**: Add to Evo2 service (not separate service)

### What to Keep (BRCA Classifier)

1. ✅ **Data Extraction**: Same ClinVar extraction approach
2. ✅ **Training Workflow**: Same 80/20 split, evaluation approach
3. ✅ **Model Format**: Same joblib format
4. ✅ **Incremental Processing**: Same async batch processing

---

## Updated Timeline

### Week 1: Build (Updated with Existing Infrastructure)

**Day 1 Morning** (2 hours):
- ✅ Verify ClinVar data (15 min) - **FASTER** (script exists)
- ✅ Verify Evo2 embeddings (30 min) - **SAME**
- ✅ Create base scorer interface (1 hour) - **SAME**

**Day 1 Afternoon** (2-3 hours):
- ✅ Filter/extract BRCA variants (30 min) - **FASTER** (can filter existing)
- ✅ Create embedding extraction script (1-2 hours) - **FASTER** (adapt existing)
- ✅ Test embedding extraction (30 min) - **SAME**

**Day 2-3** (4-6 hours):
- ✅ Extract embeddings for all BRCA variants (2-4 hours) - **SAME**
- ✅ Train XGBoost classifier (1-2 hours) - **FASTER** (adapt existing script)
- ✅ Evaluate and validate (1 hour) - **SAME**

**Day 4-5** (3-4 hours):
- ✅ Create Modal endpoint (2 hours) - **SAME**
- ✅ Create backend proxy (1 hour) - **SAME**
- ✅ Test end-to-end (1 hour) - **SAME**

**Total Week 1**: 11-15 hours (was 15-20 hours) - **FASTER** due to existing infrastructure

---

## Summary

### ✅ Good News

1. **ClinVar Extraction**: ✅ Script exists and includes BRCA1/BRCA2
2. **Training Pipeline**: ✅ Complete pipeline exists and works
3. **Deployment Pattern**: ✅ Modal service pattern proven
4. **Time Savings**: ✅ Can save 4-5 hours by adapting existing code

### ⚠️ Adaptations Needed

1. **Embedding Source**: Change from Zeta Oracle (layer 40) to Evo2 (layer 26)
2. **Model Type**: Change from GradientBoosting to XGBoost
3. **Gene Focus**: Filter for BRCA1/BRCA2 only
4. **Deployment**: Add to Evo2 service (not separate)

### 🎯 Updated Recommendation

**Before Starting Development**:
1. ✅ Verify ClinVar data (15 min) - **FASTER**
2. ✅ Verify Evo2 embeddings (30 min) - **SAME**
3. ⚠️ Create base scorer interface (1 hour) - **SAME**

**Total P0 Time**: 1.75 hours (was 2.5-4.5 hours) - **FASTER**

**Then Proceed**: Adapt existing training pipeline for BRCA classifier

---

**ANALYSIS COMPLETE** ✅  
**EXISTING INFRASTRUCTURE FOUND** ✅  
**PLAN UPDATED** ✅  
**READY TO PROCEED** (after P0 gaps, now faster!)

