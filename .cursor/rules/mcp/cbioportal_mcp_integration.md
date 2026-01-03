# cBioPortal MCP Integration

## Overview

We've installed and tested the cBioPortal MCP server from https://github.com/pickleton89/cbioportal-mcp

## Location
```
tools/cbioportal-mcp/
```

## Capabilities Confirmed

### Working Features
1. **Study Search**: Find cancer studies by keyword (500+ studies)
2. **Molecular Profiles**: Get mutation, CNV, expression profiles
3. **Mutation Data**: Get somatic mutations by gene (BRCA1, BRCA2, TP53, etc.)
4. **Sample Lists**: Get sample cohorts for studies
5. **Clinical Attributes**: List available clinical data fields
6. **Gene Search**: Find genes by symbol/keyword

### Example API Calls
```python
# Get studies
studies = await client.make_api_request('studies', params={'pageSize': 100})

# Get molecular profiles
profiles = await client.make_api_request(f'studies/{study_id}/molecular-profiles')

# Get mutations
muts = await client.make_api_request(
    f'molecular-profiles/{mutation_profile}/mutations',
    params={'entrezGeneId': 672, 'sampleListId': sample_list}  # BRCA1
)

# Get clinical data
clinical = await client.make_api_request(f'studies/{study_id}/clinical-data', params={'pageSize': 10000})
```

## Value for Our Platform

| Capability | cBioPortal MCP Use |
|------------|-------------------|
| **S/P/E Framework** | Somatic mutation data for drug efficacy validation |
| **Trial Matching** | Biomarker eligibility checking |
| **Resistance Prediction** | Mutation profiles |
| **Drug Efficacy** | Gene-drug associations |

## Cursor MCP Configuration

Add to `.cursor/mcp.json`:
```json
{
  "mcpServers": {
    "cbioportal": {
      "command": "uv",
      "args": ["run", "cbioportal-mcp"],
      "cwd": "/Users/fahadkiani/Desktop/development/crispr-assistant-main/tools/cbioportal-mcp"
    }
  }
}
```

## Key Finding: Somatic vs Germline

**cBioPortal = SOMATIC** (tumor mutations)
**PGx Validation needs = GERMLINE** (inherited variants)

For PGx validation (DPYD, TPMT, UGT1A1), we need:
- UK Biobank (germline + medication + outcomes)
- PharmGKB (curated clinical annotations)
- Partner clinical trial data with PGx screening
