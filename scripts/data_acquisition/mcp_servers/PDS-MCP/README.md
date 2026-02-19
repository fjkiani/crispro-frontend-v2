# PDS-MCP Server

MCP server for **Project Data Sphere** — exposes SAS CAS clinical trial data as native AI tools.

## Architecture

```
AI Assistant ──► PDS-MCP Server ──► SAS CAS ──► 102 Caslibs (Trial Data)
     MCP Protocol      │ swat         │ CAS Actions
                       └── .env creds ┘
```

Built on `FastMCP`. Inspired by [SAP Datasphere MCP](https://github.com/MarioDeFelipe/sap-datasphere-mcp) patterns.

## Tools (8)

| Layer | Tool | Description |
|-------|------|-------------|
| Foundation | `test_connection` | Auth to CAS, return status |
| Foundation | `health_check` | Verify env vars + SSL cert |
| Discovery | `list_caslibs` | All 102 caslibs |
| Discovery | `search_caslibs` | Filter by keyword |
| Discovery | `list_tables` | Tables in a caslib |
| Query | `get_table_schema` | Column names, types, row count |
| Query | `preview_table` | First N rows (max 500) |
| Query | `extract_cohort` | Full extraction with filters |

## Quick Start

```bash
cd scripts/data_acquisition/mcp_servers/PDS-MCP
cp .env.example .env
# Fill in SAS_USERNAME, SAS_PASSWORD, SSL_CERT_PATH
pip install -e .
python -m pds_server
```

## Cursor MCP Registration

Add to `.cursor/mcp.json`:

```json
{
  "mcpServers": {
    "pds-mcp": {
      "command": "python",
      "args": ["-m", "pds_server"],
      "cwd": "scripts/data_acquisition/mcp_servers/PDS-MCP"
    }
  }
}
```
