# Validation Improvements Scripts

**Date:** January 21, 2026  
**Purpose:** Address critical issues identified in manuscript audit

---

## 📋 Scripts Overview

### **Task 1.1: Address Overfitting**

1. **`reduce_model_complexity.py`**
   - Reduces 8-pathway model to 3-4 pathways
   - Compares performance (training AUC, CV AUC, EPV ratio)
   - Saves reduced models as JSON

2. **`bootstrap_validation.py`**
   - Performs bootstrap validation (1000 iterations)
   - Computes optimism-corrected AUC
   - Reports overfitting amount

3. **`learning_curve_analysis.py`**
   - Trains models on increasing sample sizes (n=20, 30, 40, 51)
   - Evaluates on held-out test set
   - Generates learning curve plot

### **Task 1.3: Statistical Fixes**

4. **`multiple_testing_correction.py`**
   - Applies Benjamini-Hochberg FDR correction
   - Applies Bonferroni correction
   - Updates Table 1 with corrected p-values

5. **`calibration_analysis.py`**
   - Bins predicted probabilities into 5 quantile groups
   - Computes observed response rate per bin
   - Generates calibration plot
   - Computes Hosmer-Lemeshow test

6. **`improved_cross_validation.py`**
   - Compares 5-fold CV, LOOCV, and bootstrap validation
   - Reports most stable method for n=51
   - Saves comparison results

### **Task 1.4: External Validation**

7. **`external_validation_gse179994.py`**
   - Validates model on independent NSCLC cohort (GSE179994)
   - **Note:** Requires bulk RNA-seq expression data
   - Currently provides template and instructions

### **Task 1.5: Compare Published Signatures**

8. **`compare_published_signatures.py`**
   - Implements simplified versions of TIDE, IMPRES, TMEscore
   - Performs head-to-head comparison on GSE91061
   - Uses DeLong test for statistical comparison
   - **Note:** Full implementations require complete expression matrix

### **Master Script**

9. **`run_all_improvements.py`**
   - Runs all scripts in sequence
   - Provides execution summary
   - Handles errors gracefully

---

## 🚀 Quick Start

### Run All Scripts

```bash
cd publications/06-io-response-prediction
python scripts/validation_improvements/run_all_improvements.py
```

### Run Individual Scripts

```bash
# Task 1.1A: Reduce model complexity
python scripts/validation_improvements/reduce_model_complexity.py

# Task 1.1B: Bootstrap validation
python scripts/validation_improvements/bootstrap_validation.py

# Task 1.1C: Learning curve
python scripts/validation_improvements/learning_curve_analysis.py

# Task 1.3A: Multiple testing correction
python scripts/validation_improvements/multiple_testing_correction.py

# Task 1.3B: Calibration analysis
python scripts/validation_improvements/calibration_analysis.py

# Task 1.3C: Improved CV
python scripts/validation_improvements/improved_cross_validation.py

# Task 1.5: Signature comparison
python scripts/validation_improvements/compare_published_signatures.py
```

---

## 📊 Expected Outputs

### Data Files (in `data/` directory)

- `reduced_model_3pathway.json` - 3-pathway model
- `reduced_model_4pathway.json` - 4-pathway model
- `reduced_model_comparison.csv` - Model comparison table
- `bootstrap_validation_results.json` - Bootstrap results
- `bootstrap_iterations.csv` - Detailed bootstrap iterations
- `learning_curve_results.csv` - Learning curve data
- `pathway_multiple_testing_corrected.csv` - FDR/Bonferroni corrected p-values
- `table1_pathway_performance_corrected.csv` - Publication-ready table
- `calibration_analysis_results.json` - Calibration metrics
- `calibration_bins.csv` - Calibration bin data
- `improved_cv_comparison.json` - CV method comparison
- `improved_cv_comparison.csv` - CV comparison table
- `signature_comparison_gse91061.csv` - Signature comparison table
- `signature_comparison_results.json` - Detailed signature results

### Figures (in `figures/` directory)

- `learning_curve.png` / `.pdf` - Learning curve plot
- `calibration_plot.png` / `.pdf` - Calibration curve

---

## ⚠️ Dependencies

### Required Python Packages

```bash
pip install pandas numpy scipy scikit-learn matplotlib seaborn statsmodels tqdm
```

### Required Data Files

All scripts expect data files in:
```
scripts/data_acquisition/IO/
  - gse91061_analysis_with_composites.csv
  - gse91061_pathway_response_association.csv
  - gse91061_benchmark_comparison.csv
```

### Optional Data Files (for Task 1.4)

- `GSE179994_bulk_expression.tsv` - Bulk RNA-seq expression (or aggregated scRNA-seq)
- `GSE179994_clinical.csv` - Clinical metadata with response labels

---

## 📝 Notes

### Task 1.4 (External Validation)

GSE179994 is **single-cell RNA-seq data**, not bulk RNA-seq. To complete external validation:

1. **Option A:** Aggregate single-cell data to pseudo-bulk expression
   - Sum counts per gene per sample
   - Normalize to TPM
   - Use same gene lists as GSE91061

2. **Option B:** Use alternative bulk RNA-seq dataset
   - GSE168204 (if available)
   - Other melanoma/NSCLC bulk RNA-seq cohorts with IO response data

### Task 1.5 (Published Signatures)

Current implementations are **simplified proxies** using pathway scores. For full comparison:

1. **TIDE:** Download from web server (http://tide.dfci.harvard.edu/) or use R package
2. **IMPRES:** Implement full 15 gene-pair algorithm from paper
3. **TMEscore:** Implement with proper gene lists and weights from paper

Full implementations require:
- Complete expression matrix (not just pathway scores)
- All genes from published gene lists
- Proper weighting schemes and algorithms

---

## ✅ Success Criteria

After running all scripts, you should have:

1. ✅ Reduced model (3-4 pathways) with improved EPV ratio
2. ✅ Bootstrap-corrected AUC (realistic performance estimate)
3. ✅ Learning curve showing sample size adequacy
4. ✅ FDR/Bonferroni corrected p-values for Table 1
5. ✅ Calibration plot and Hosmer-Lemeshow test
6. ✅ Comparison of CV methods (5-fold, LOOCV, bootstrap)
7. ⚠️ External validation (pending data availability)
8. ⚠️ Full signature comparison (pending expression matrix)

---

## 🔧 Troubleshooting

### Issue: "FileNotFoundError: gse91061_analysis_with_composites.csv"

**Solution:** Ensure you're running from the correct directory and data files exist in `scripts/data_acquisition/IO/`

### Issue: "ModuleNotFoundError: statsmodels"

**Solution:** Install missing package: `pip install statsmodels`

### Issue: Bootstrap validation takes too long

**Solution:** Reduce `N_BOOTSTRAP` in `bootstrap_validation.py` (default: 1000)

### Issue: External validation script shows "AWAITING DATA"

**Solution:** This is expected. Provide bulk expression data for GSE179994 or use alternative dataset.

---

## 📚 References

- **TIDE:** Jiang et al. Nature Medicine 2018, PMID: 29910795
- **IMPRES:** Auslander et al. Nature Medicine 2018, PMID: 29910796
- **TMEscore:** Zeng et al. Cancer Immunology Research 2019, PMID: 31088836

---

**Status:** ✅ All scripts created and ready to run  
**Next:** Run `run_all_improvements.py` to execute all validation improvements
