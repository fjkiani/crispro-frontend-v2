"""
PDS-MCP Server — FastMCP tool surface for Project Data Sphere
==============================================================
Modeled after BioMed-MCP. Hardened with SAP Datasphere MCP patterns:
  - Structured response envelopes (status/data/error)
  - Connection caching with auto-reconnect
  - MCP Context error reporting
  - Layered tool taxonomy: Foundation → Discovery → Query
"""

import os
import json
import logging
from typing import Optional
from fastmcp import FastMCP, Context
from dotenv import load_dotenv

load_dotenv()

from .cas_client import CASClient

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pds-mcp-server")

# Initialize FastMCP app
app = FastMCP("pds-mcp-server")

# Lazy singleton (SAP pattern: connection caching)
_client: Optional[CASClient] = None


def get_client() -> CASClient:
    global _client
    if _client is None:
        _client = CASClient()
    return _client


def _envelope(status: str, data=None, error: str = None, tool: str = "") -> str:
    """SAP-inspired structured response envelope."""
    resp = {"tool": tool, "status": status}
    if data is not None:
        resp["data"] = data
    if error:
        resp["error"] = error
    return json.dumps(resp, indent=2, default=str)


# ═══════════════════════════════════════════════════════════════
# FOUNDATION                  (SAP Datasphere: 5 foundation tools)
# ═══════════════════════════════════════════════════════════════

@app.tool(
    annotations={
        "title": "Test PDS Connection",
        "description": "Authenticate to SAS CAS and return server status. Run this first.",
        "readOnlyHint": True,
    }
)
async def test_connection(ctx: Context) -> str:
    """Connect to Project Data Sphere CAS server. Returns connection status."""
    try:
        client = get_client()
        result = client.connect()
        if result["status"] == "connected":
            return _envelope("success", result, tool="test_connection")
        else:
            await ctx.error(f"Connection failed: {result.get('message', 'unknown')}")
            return _envelope("error", error=result.get("message", "unknown"), tool="test_connection")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="test_connection")


@app.tool(
    annotations={
        "title": "PDS Health Check",
        "description": "Verify env vars and SSL cert for PDS access",
        "readOnlyHint": True,
    }
)
async def health_check(ctx: Context) -> str:
    """Check that all required config is present."""
    required = ["SAS_USERNAME", "SAS_PASSWORD"]
    missing = [v for v in required if not os.getenv(v)]
    ssl_path = os.getenv("SSL_CERT_PATH", "")
    ssl_ok = os.path.exists(ssl_path) if ssl_path else False

    status = "healthy" if (not missing and ssl_ok) else "misconfigured"
    data = {
        "missing_env_vars": missing,
        "ssl_cert_path": ssl_path,
        "ssl_cert_exists": ssl_ok,
    }
    if status != "healthy":
        await ctx.error(f"Misconfigured: missing={missing}, ssl_ok={ssl_ok}")
    return _envelope(status, data, tool="health_check")


# ═══════════════════════════════════════════════════════════════
# DISCOVERY               (SAP Datasphere: Space/Catalog discovery)
# ═══════════════════════════════════════════════════════════════

@app.tool(
    annotations={
        "title": "List all PDS Caslibs",
        "description": "Return all caslibs (data libraries) from Project Data Sphere. This is the equivalent of SAP 'list_spaces'.",
        "readOnlyHint": True,
        "openWorldHint": True,
    }
)
async def list_caslibs(ctx: Context) -> str:
    """List all caslibs in Project Data Sphere. Returns name, description, path for each."""
    try:
        client = get_client()
        caslibs = client.list_caslibs()
        return _envelope("success", {"count": len(caslibs), "caslibs": caslibs}, tool="list_caslibs")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="list_caslibs")


@app.tool(
    annotations={
        "title": "Search PDS Caslibs",
        "description": "Filter caslibs by keyword (e.g., 'ovarian', 'PARP', 'platinum'). Equivalent to SAP 'search_catalog'.",
        "readOnlyHint": True,
        "openWorldHint": True,
    }
)
async def search_caslibs(ctx: Context, keyword: str) -> str:
    """
    Search caslibs by keyword (case-insensitive).
    
    Args:
        keyword: Search term to match against caslib names and descriptions
    """
    try:
        client = get_client()
        matches = client.search_caslibs(keyword)
        return _envelope("success", {"keyword": keyword, "matches": len(matches), "caslibs": matches}, tool="search_caslibs")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="search_caslibs")


@app.tool(
    annotations={
        "title": "List tables in a PDS caslib",
        "description": "List all files/tables within a specific caslib. Equivalent to SAP 'list_catalog_assets'.",
        "readOnlyHint": True,
        "openWorldHint": True,
    }
)
async def list_tables(ctx: Context, caslib: str) -> str:
    """
    List files/tables in a caslib.
    
    Args:
        caslib: Name of the caslib to inspect
    """
    try:
        client = get_client()
        tables = client.list_tables(caslib)
        return _envelope("success", {"caslib": caslib, "table_count": len(tables), "tables": tables}, tool="list_tables")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="list_tables")


# ═══════════════════════════════════════════════════════════════
# QUERY / EXTRACTION     (SAP Datasphere: Metadata + ETL tools)
# ═══════════════════════════════════════════════════════════════

@app.tool(
    annotations={
        "title": "Get table schema from PDS",
        "description": "Load a table header and return column names, types, and row count. Equivalent to SAP 'get_table_schema'.",
        "readOnlyHint": True,
        "openWorldHint": True,
    }
)
async def get_table_schema(ctx: Context, caslib: str, table_path: str) -> str:
    """
    Get schema (columns, types, row count) for a table.
    
    Args:
        caslib: Caslib containing the table
        table_path: Path/name of the table file
    """
    try:
        client = get_client()
        schema = client.get_table_schema(caslib, table_path)
        if "error" in schema:
            await ctx.error(schema["error"])
            return _envelope("error", error=schema["error"], tool="get_table_schema")
        return _envelope("success", schema, tool="get_table_schema")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="get_table_schema")


@app.tool(
    annotations={
        "title": "Preview table data from PDS",
        "description": "Load first N rows of a table. Use this for quick inspection before full extraction.",
        "readOnlyHint": True,
        "openWorldHint": True,
    }
)
async def preview_table(ctx: Context, caslib: str, table_path: str, limit: int = 50) -> str:
    """
    Preview first N rows of a table.
    
    Args:
        caslib: Caslib containing the table
        table_path: Path/name of the table file
        limit: Max rows to return (default 50, max 500)
    """
    try:
        limit = min(max(1, limit), 500)  # Guard: cap at 500
        client = get_client()
        data = client.preview_table(caslib, table_path, limit=limit)
        if "error" in data:
            await ctx.error(data["error"])
            return _envelope("error", error=data["error"], tool="preview_table")
        return _envelope("success", data, tool="preview_table")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="preview_table")


@app.tool(
    annotations={
        "title": "Extract cohort from PDS",
        "description": "Full data extraction with optional disease/biomarker filters. Equivalent to SAP 'query_relational_entity'.",
        "readOnlyHint": False,
        "openWorldHint": True,
    }
)
async def extract_cohort(
    ctx: Context,
    caslib: str,
    table_path: str,
    disease_filter: Optional[str] = None,
    biomarker_filter: Optional[str] = None,
    max_rows: int = 5000
) -> str:
    """
    Extract cohort data with optional filters.
    
    Args:
        caslib: Caslib containing the data
        table_path: Path/name of the table
        disease_filter: Optional disease keyword filter (e.g., 'ovarian')
        biomarker_filter: Optional biomarker keyword filter (e.g., 'BRCA')
        max_rows: Maximum rows to extract (default 5000, max 50000)
    """
    try:
        max_rows = min(max(1, max_rows), 50000)  # Guard
        client = get_client()
        result = client.extract_cohort(
            caslib, table_path,
            disease_filter=disease_filter,
            biomarker_filter=biomarker_filter,
            max_rows=max_rows
        )
        if "error" in result:
            await ctx.error(result["error"])
            return _envelope("error", error=result["error"], tool="extract_cohort")
        return _envelope("success", result, tool="extract_cohort")
    except Exception as e:
        await ctx.error(str(e))
        return _envelope("error", error=str(e), tool="extract_cohort")


# ═══════════════════════════════════════════════════════════════
# RESOURCE
# ═══════════════════════════════════════════════════════════════

@app.resource("pds://health")
def pds_health() -> str:
    """Health resource."""
    return "PDS-MCP server is operational"


def main():
    """Run the PDS-MCP server."""
    logger.info("🔱 PDS-MCP Server starting...")
    app.run()


if __name__ == "__main__":
    main()
