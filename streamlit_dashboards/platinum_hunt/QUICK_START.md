# ⚔️ QUICK START GUIDE

## 🚀 Running the Dashboard

### Option 1: Direct Streamlit Command
```bash
cd streamlit_dashboards/platinum_hunt
streamlit run app.py
```

### Option 2: Using Launch Script
```bash
cd streamlit_dashboards/platinum_hunt
./run_dashboard.sh
```

### Option 3: From Project Root
```bash
streamlit run streamlit_dashboards/platinum_hunt/app.py
```

## 📦 Install Dependencies

```bash
pip install -r streamlit_dashboards/platinum_hunt/requirements.txt
```

## ✅ Verify Data Files Exist

The dashboard expects these files:
- `data/validation/tcga_ov_platinum_response_labels.json` (469 patients)
- `data/validation/tcga_ov_full_validation_dataset.json` (200 patients)

If files are missing, run the extraction first:
```bash
python scripts/platinum_hunt/orchestrator.py
```

## 🎯 Dashboard Pages

1. **Overview** - Hero metrics, response distribution, validation status
2. **Patients** - Searchable patient explorer with filters
3. **Overlap Analysis** - Venn diagrams, match rates, sample ID mapping
4. **Validation Readiness** - Statistical power analysis

## 🐛 Troubleshooting

**Import Errors:**
- Make sure you're running from the correct directory
- Check that all dependencies are installed

**Data Not Loading:**
- Verify data files exist in `data/validation/`
- Check file permissions

**Venn Diagram Not Showing:**
- Install `matplotlib-venn`: `pip install matplotlib-venn`


