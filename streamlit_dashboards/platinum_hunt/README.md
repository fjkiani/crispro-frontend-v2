# ⚔️ Platinum Response Data Hunt Dashboard

Beautiful, modular Streamlit dashboard showcasing extracted platinum response data.

## 🚀 Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run dashboard
streamlit run app.py
```

## 📁 Structure

```
platinum_hunt/
├── app.py                    # Main entry point
├── config.py                 # Configuration
├── data/                     # Data loading & processing
│   ├── loader.py
│   └── processor.py
├── components/               # Reusable components
│   ├── hero_metrics.py
│   ├── response_charts.py
│   ├── overlap_analysis.py
│   ├── patient_table.py
│   └── validation_status.py
└── pages/                    # Dashboard pages
    ├── overview.py
    ├── patients.py
    ├── overlap.py
    └── validation.py
```

## 📊 Features

- **Overview Dashboard** - Hero metrics, response distribution, validation status
- **Patient Explorer** - Searchable, filterable patient table
- **Overlap Analysis** - Venn diagrams, match rates, sample ID mapping
- **Validation Readiness** - Statistical power analysis, sample size metrics

## 🎯 Data Sources

- `data/validation/tcga_ov_platinum_response_labels.json` - Jr2's 469 patients
- `data/validation/tcga_ov_full_validation_dataset.json` - Zo's 200 patients

## ✅ Status

- ✅ 469 patients extracted
- ✅ 161 patients overlap (34.3% match rate)
- ✅ Validation ready (N=161 exceeds ≥40 threshold by 4x)


