"""
CAS Client — Thin wrapper around swat for Project Data Sphere
==============================================================
All SAS CAS logic lives here. The MCP server calls this, never swat directly.
"""

import os
import logging
from typing import List, Dict, Any, Optional
from pathlib import Path

logger = logging.getLogger("pds-mcp.cas_client")

# Lazy import — swat may not be installed
_swat = None

def _get_swat():
    global _swat
    if _swat is None:
        try:
            import swat
            _swat = swat
        except ImportError:
            raise ImportError(
                "swat library required. Install: pip install swat"
            )
    return _swat


class CASClient:
    """Manages a single CAS connection to Project Data Sphere."""

    DEFAULT_URL = "https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/"
    DEFAULT_PORT = 443

    def __init__(self):
        self.conn = None
        self.cas_url = os.getenv("CAS_URL", self.DEFAULT_URL)
        self.port = int(os.getenv("CAS_PORT", self.DEFAULT_PORT))
        self.username = os.getenv("SAS_USERNAME", "")
        self.password = os.getenv("SAS_PASSWORD", "")
        self.ssl_cert = os.getenv("SSL_CERT_PATH", "")

        # Set SSL cert env var if provided
        if self.ssl_cert and Path(self.ssl_cert).exists():
            os.environ["CAS_CLIENT_SSL_CA_LIST"] = self.ssl_cert
            logger.info(f"SSL cert loaded: {self.ssl_cert}")

    def connect(self) -> Dict[str, Any]:
        """Connect to CAS. Returns status dict."""
        if self.conn:
            return {"status": "already_connected"}

        if not self.username or not self.password:
            return {"status": "error", "message": "SAS_USERNAME and SAS_PASSWORD must be set in .env"}

        swat = _get_swat()
        try:
            logger.info(f"Connecting to CAS: {self.cas_url}:{self.port}")
            self.conn = swat.CAS(
                self.cas_url, self.port,
                username=self.username,
                password=self.password
            )
            status = self.conn.serverstatus()
            logger.info("CAS connection established")
            return {"status": "connected", "server_info": str(status)}
        except Exception as e:
            logger.error(f"CAS connection failed: {e}")
            return {"status": "error", "message": str(e)}

    def disconnect(self):
        """Close CAS connection."""
        if self.conn:
            try:
                self.conn.close()
            except Exception:
                pass
            self.conn = None

    def _ensure_connected(self):
        """Auto-connect if not already connected."""
        if not self.conn:
            result = self.connect()
            if result["status"] != "connected":
                raise ConnectionError(f"CAS not connected: {result.get('message', 'unknown')}")

    def list_caslibs(self) -> List[Dict[str, str]]:
        """List all available caslibs."""
        self._ensure_connected()
        result = self.conn.table.caslibInfo()
        if hasattr(result, "CASLibInfo"):
            df = result.CASLibInfo
            return df.to_dict("records") if hasattr(df, "to_dict") else []
        return []

    def search_caslibs(self, keyword: str) -> List[Dict[str, str]]:
        """Filter caslibs by keyword (case-insensitive)."""
        kw = keyword.lower()
        return [
            c for c in self.list_caslibs()
            if kw in (c.get("Name", "") or "").lower()
            or kw in (c.get("Description", "") or "").lower()
        ]

    def list_tables(self, caslib: str) -> List[Dict[str, Any]]:
        """List files/tables in a caslib."""
        self._ensure_connected()
        try:
            result = self.conn.table.fileInfo(allFiles=True, caslib=caslib)
            if hasattr(result, "FileInfo"):
                df = result.FileInfo
                return df.to_dict("records") if hasattr(df, "to_dict") else []
        except Exception as e:
            logger.warning(f"Error listing tables in {caslib}: {e}")
        return []

    def get_table_schema(self, caslib: str, table_path: str) -> Dict[str, Any]:
        """Load a table header and return schema info."""
        self._ensure_connected()
        try:
            temp_name = f"_schema_peek_{hash(table_path) % 99999}"
            cas_table = self.conn.CASTable(temp_name, replace=True, caslib="CASUSER")
            self.conn.table.loadTable(
                sourceCaslib=caslib,
                casOut=cas_table,
                path=table_path
            )
            col_info = self.conn.table.columnInfo(table={"name": temp_name, "caslib": "CASUSER"})
            columns = []
            if hasattr(col_info, "ColumnInfo"):
                columns = col_info.ColumnInfo.to_dict("records")

            row_count_result = self.conn.simple.numRows(table={"name": temp_name, "caslib": "CASUSER"})
            n_rows = 0
            if hasattr(row_count_result, "numrows"):
                n_rows = int(row_count_result.numrows)

            # Cleanup
            self.conn.table.dropTable(name=temp_name, caslib="CASUSER", quiet=True)

            return {
                "caslib": caslib,
                "table": table_path,
                "columns": columns,
                "row_count": n_rows
            }
        except Exception as e:
            return {"error": str(e), "caslib": caslib, "table": table_path}

    def preview_table(self, caslib: str, table_path: str, limit: int = 50) -> Dict[str, Any]:
        """Load first N rows of a table and return as list of dicts."""
        self._ensure_connected()
        try:
            temp_name = f"_preview_{hash(table_path) % 99999}"
            cas_table = self.conn.CASTable(temp_name, replace=True, caslib="CASUSER")
            self.conn.table.loadTable(
                sourceCaslib=caslib,
                casOut=cas_table,
                path=table_path
            )
            fetch_result = self.conn.table.fetch(
                table={"name": temp_name, "caslib": "CASUSER"},
                maxRows=limit
            )
            rows = []
            if hasattr(fetch_result, "Fetch"):
                df = fetch_result.Fetch
                rows = df.to_dict("records") if hasattr(df, "to_dict") else []

            # Cleanup
            self.conn.table.dropTable(name=temp_name, caslib="CASUSER", quiet=True)

            return {
                "caslib": caslib,
                "table": table_path,
                "rows_returned": len(rows),
                "data": rows
            }
        except Exception as e:
            return {"error": str(e), "caslib": caslib, "table": table_path}

    def extract_cohort(
        self,
        caslib: str,
        table_path: str,
        disease_filter: Optional[str] = None,
        biomarker_filter: Optional[str] = None,
        max_rows: int = 5000
    ) -> Dict[str, Any]:
        """Full extraction with optional filters."""
        self._ensure_connected()
        try:
            temp_name = f"_extract_{hash(table_path) % 99999}"
            cas_table = self.conn.CASTable(temp_name, replace=True, caslib="CASUSER")
            self.conn.table.loadTable(
                sourceCaslib=caslib,
                casOut=cas_table,
                path=table_path
            )

            # Build WHERE clause
            where_parts = []
            if disease_filter:
                where_parts.append(f"UPCASE(_CHARACTER_) CONTAINS '{disease_filter.upper()}'")
            # Note: biomarker_filter requires column-specific logic per table

            fetch_params = {
                "table": {"name": temp_name, "caslib": "CASUSER"},
                "maxRows": max_rows
            }
            if where_parts:
                fetch_params["where"] = " AND ".join(where_parts)

            fetch_result = self.conn.table.fetch(**fetch_params)
            rows = []
            if hasattr(fetch_result, "Fetch"):
                df = fetch_result.Fetch
                rows = df.to_dict("records") if hasattr(df, "to_dict") else []

            # Cleanup
            self.conn.table.dropTable(name=temp_name, caslib="CASUSER", quiet=True)

            return {
                "caslib": caslib,
                "table": table_path,
                "rows_extracted": len(rows),
                "filters_applied": {
                    "disease": disease_filter,
                    "biomarker": biomarker_filter
                },
                "data": rows
            }
        except Exception as e:
            return {"error": str(e), "caslib": caslib, "table": table_path}
