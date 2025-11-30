# 📖 EXTRACTION PIECE 4.1: Biomarker Correlation Service Implementation
**Source**: Lines 6400-6740 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section documents the complete implementation of the `BiomarkerCorrelationService`, including all statistical methods, correlation computations, feature ranking, and analysis pipeline.

---

## 🔍 KEY FINDINGS

### **Service Architecture**

**Location**: `oncology-coPilot/oncology-backend-minimal/api/services/biomarker_correlation_service.py`

**Main Class**: `BiomarkerCorrelationService`

**Key Components:**
1. Data loading and preparation
2. Statistical correlation computation
3. Feature ranking and filtering
4. Bootstrap confidence intervals
5. Summary generation

---

### **Statistical Methods Implemented**

#### **1. Pearson Correlation**
- **Purpose**: Linear relationship between SAE features and outcomes
- **Implementation**: `scipy.stats.pearsonr` for each feature
- **Output**: Correlation coefficient (r) and p-value
- **Use Case**: Primary ranking metric

#### **2. Spearman Correlation**
- **Purpose**: Non-parametric rank-based correlation
- **Implementation**: `scipy.stats.spearmanr` for each feature
- **Output**: Rank correlation coefficient and p-value
- **Use Case**: Validates Pearson results, handles non-linear relationships

#### **3. Chi-Square Test**
- **Purpose**: Categorical association between features and outcomes
- **Implementation**: Discretizes features into bins, computes chi-square statistic
- **Output**: Chi-square statistic and p-value
- **Use Case**: Detects categorical patterns

#### **4. Cohen's d Effect Size**
- **Purpose**: Standardized mean difference between groups
- **Implementation**: Computes mean difference divided by pooled standard deviation
- **Output**: Effect size (d) for each feature
- **Threshold**: `EFFECT_SIZE_THRESHOLD = 0.3` (small effect minimum)
- **Use Case**: Validates practical significance beyond statistical significance

#### **5. Cross-Validation Stability**
- **Purpose**: Measures feature stability across data splits
- **Implementation**: 5-fold CV, counts how often feature is selected
- **Output**: Stability score (0.0-1.0)
- **Threshold**: `CV_STABILITY_THRESHOLD = 0.6` (selected in 60%+ folds)
- **Use Case**: Ensures features aren't spurious correlations

#### **6. Bootstrap Confidence Intervals**
- **Purpose**: Estimates uncertainty in correlation estimates
- **Implementation**: Resample with replacement, compute correlation distribution
- **Output**: 95% CI (lower, upper) for top features
- **Iterations**: `BOOTSTRAP_ITERATIONS = 1000`
- **Use Case**: Quantifies confidence in top features

#### **7. Multiple Testing Correction**
- **Purpose**: Controls false discovery rate across 32K features
- **Implementation**: Benjamini-Hochberg FDR correction
- **Method**: `scipy.stats.false_discovery_control(method='bh')`
- **Output**: Corrected p-values
- **Use Case**: Prevents false positives from multiple comparisons

---

### **Feature Ranking Criteria**

**Multi-Stage Filtering:**

1. **Statistical Significance**: `pearson_p < P_VALUE_THRESHOLD` (default: 0.05)
2. **Effect Size**: `cohen_d >= EFFECT_SIZE_THRESHOLD` (default: 0.3)
3. **Stability**: `cv_stability >= CV_STABILITY_THRESHOLD` (default: 0.6)
4. **Ranking**: Sort by `|pearson_r|` descending
5. **Top N**: Select top `TOP_N_FEATURES` (default: 20)

**Code:**
```python
significant_features = [
    fc for fc in correlations
    if (fc.pearson_p < P_VALUE_THRESHOLD and
        fc.cohen_d >= EFFECT_SIZE_THRESHOLD and
        fc.cv_stability >= CV_STABILITY_THRESHOLD)
]

significant_features.sort(key=lambda x: abs(x.pearson_r), reverse=True)
top_features = significant_features[:TOP_N_FEATURES]
```

---

### **Data Structures**

#### **FeatureCorrelation**
```python
@dataclass
class FeatureCorrelation:
    feature_index: int
    pearson_r: float
    pearson_p: float  # FDR-corrected
    spearman_r: float
    spearman_p: float
    chi_square: float
    chi_square_p: float
    cohen_d: float
    cv_stability: float
    bootstrap_ci_lower: float
    bootstrap_ci_upper: float
    rank: int = 0
```

#### **BiomarkerSummary**
```python
@dataclass
class BiomarkerSummary:
    top_features: List[Dict]
    total_features_analyzed: int
    significant_features_count: int
    p_value_threshold: float
    effect_size_threshold: float
    cv_stability_threshold: float
    correction_method: str
    cohort_size: int
    outcome_distribution: Dict[str, int]
    provenance: Dict[str, Any]
```

---

### **Analysis Pipeline**

**`run_analysis()` Method:**

1. **Load Data**: `load_data()`
   - Loads SAE cohort JSON
   - Builds feature matrix (patients × 32K features)
   - Extracts outcome vectors

2. **Compute Correlations**: `compute_correlations()`
   - Pearson, Spearman, Chi-square, Cohen's d
   - Cross-validation stability
   - FDR correction

3. **Rank and Filter**: `rank_and_filter_features()`
   - Apply thresholds
   - Sort by correlation strength
   - Select top N

4. **Bootstrap CIs**: `compute_bootstrap_cis_for_top_features()`
   - Compute confidence intervals for top features
   - Update FeatureCorrelation objects

5. **Generate Summary**: `generate_summary()`
   - Create BiomarkerSummary object
   - Include provenance metadata

---

### **Configuration Constants**

```python
# Statistical thresholds
P_VALUE_THRESHOLD = 0.05
EFFECT_SIZE_THRESHOLD = 0.3  # Small effect minimum
CV_STABILITY_THRESHOLD = 0.6  # Selected in 60%+ folds
TOP_N_FEATURES = 20

# Analysis parameters
RANDOM_SEED = 42  # Reproducibility
BOOTSTRAP_ITERATIONS = 1000
CV_FOLDS = 5

# Outcome encoding
OUTCOME_ENCODING = {
    "sensitive": 1.0,
    "resistant": 0.5,
    "refractory": 0.0
}
```

---

### **Data Loading**

**`load_sae_cohort_data()`:**
- Reads JSON from `SAE_COHORT_FILE`
- Handles both old (list) and new (dict with 'patients') formats
- Filters patients with valid variants and SAE features
- Extracts outcome labels (handles both 'platinum_response' and 'outcome' fields)
- Returns: patients list, outcome vector, patient IDs

**`build_feature_matrix()`:**
- Aggregates `top_features` from all variants per patient
- Creates 32K-dimensional vector per patient
- Sums feature activations across variants
- Normalizes by number of variants (mean activation)
- Returns: `[n_patients, 32768]` numpy array

---

### **RUO Endpoint**

**Route**: `POST /api/sae/biomarker_summary`

**Purpose**: Returns biomarker analysis results for internal inspection

**Gating**: `ENABLE_TRUE_SAE` feature flag

**Output**: BiomarkerSummary JSON

**Warning**: RUO/VALIDATION-ONLY, not for clinical use

---

## 📊 KEY INSIGHTS

### **Statistical Rigor**

1. **Multiple Methods**: Uses 5+ statistical tests to validate findings
2. **Effect Size**: Requires practical significance, not just statistical
3. **Stability**: Cross-validation ensures features aren't spurious
4. **Multiple Testing**: FDR correction prevents false discoveries
5. **Uncertainty**: Bootstrap CIs quantify confidence

### **Performance Considerations**

1. **Efficient**: Processes 32K features × 69 patients in ~2 minutes
2. **Memory**: Handles large feature matrices efficiently
3. **Scalable**: Can handle larger cohorts (tested up to 469 patients)
4. **Reproducible**: Fixed random seed ensures consistent results

### **Design Decisions**

1. **Multi-Stage Filtering**: Combines multiple criteria for robustness
2. **Ranking by Correlation**: Prioritizes strongest associations
3. **Bootstrap Only Top Features**: Computationally efficient
4. **Provenance Tracking**: Full metadata for reproducibility

---

## 🔗 CONTEXT & CONNECTIONS

- **Builds on**: Mock data testing (Piece 3.1), Real data extraction (Piece 3.2)
- **Enables**: Biomarker discovery from SAE features
- **Integrates with**: SAE feature extraction pipeline
- **Outputs to**: RUO endpoint for internal inspection

---

## 📝 NOTES

- Service is production-ready and tested
- All statistical methods validated with mock data
- Handles edge cases (missing data, zero variance)
- Full provenance tracking for reproducibility
- Designed for RUO/validation use, not clinical

---

## 🎯 QUESTIONS RESOLVED

- ✅ What statistical methods are used? → Pearson, Spearman, Chi-square, Cohen's d, CV, Bootstrap
- ✅ How are features ranked? → Multi-stage filtering + correlation strength
- ✅ How is multiple testing handled? → FDR Benjamini-Hochberg correction
- ✅ What's the performance? → ~2 minutes for 69 patients × 32K features
- ✅ Is it reproducible? → Yes, fixed random seed and full provenance

