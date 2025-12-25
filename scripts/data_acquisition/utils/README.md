# Data Acquisition Utilities

Modular, reusable utilities for TCGA/cBioPortal data acquisition.

## Modules

### `cbioportal_client.py`
Reusable cBioPortal API client with:
- Rate limiting and retry logic
- Caching of mutation profile IDs
- Batch clinical data fetching
- Error handling with exponential backoff

**Usage:**
```python
from cbioportal_client import CBioportalClient

client = CBioportalClient()
patients = client.get_study_patients("ov_tcga_pan_can_atlas_2018")
mutations = client.get_mutations_for_samples(profile_id, sample_ids)
```

### `tmb_calculator.py`
TMB (Tumor Mutational Burden) calculation utilities:
- Extract TMB from clinical data (multiple field name variations)
- Calculate TMB from mutation counts
- Normalize and validate TMB values
- Classify TMB as high/low

**Usage:**
```python
from tmb_calculator import extract_tmb_from_clinical, calculate_tmb_from_mutations

tmb_provided = extract_tmb_from_clinical(clinical_data)
tmb_calculated = calculate_tmb_from_mutations(mutations, exome_size_mb=38.0)
```

### `mutation_processor.py`
Mutation processing utilities:
- Extract and standardize mutation fields
- Find primary tumor samples
- Detect specific mutations (MBD4, TP53, DDR, MMR)

**Usage:**
```python
from mutation_processor import process_mutations, has_mbd4_mutation

mutations = process_mutations(raw_mutations)
has_mbd4 = has_mbd4_mutation(mutations)
```

### `data_quality.py`
Data quality validation and flagging:
- Normalize MSI status values
- Generate data quality flags
- Deduplicate patient records
- Validate patient records

**Usage:**
```python
from data_quality import extract_msi_from_clinical, generate_data_quality_flags

msi_status = extract_msi_from_clinical(clinical_data)
flags = generate_data_quality_flags(has_tmb=True, has_msi=True, ...)
```

## Design Principles

1. **Modularity**: Each module is independent and reusable
2. **No Dependencies**: Utilities don't depend on each other (except imports)
3. **Future-Proof**: Can be reused for other data acquisition tasks
4. **Error Handling**: Graceful degradation when data is missing

## Future Extensions

These utilities can be extended for:
- Other TCGA studies
- Non-TCGA cBioPortal studies
- Custom TMB calculation methods
- Additional mutation flags
- Other biomarker extraction

