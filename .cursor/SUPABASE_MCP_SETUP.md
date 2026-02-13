# Supabase MCP Server Setup

Cursor is configured to use **Supabase’s hosted MCP server** so you can query your Supabase project from the AI assistant.

## Config in use

- **File:** `.cursor/mcp.json`
- **Server:** `https://mcp.supabase.com/mcp`
- **Project:** Scoped to `xfhiwodulrbbtfcqneqt` (from `oncology-coPilot/oncology-backend-minimal/.env`).
- **Mode:** `read_only=true` (no writes; recommended for safety).

## What you need to do

1. **Restart Cursor** (or reload the window) so it picks up `.cursor/mcp.json`.
2. **Sign in when prompted**  
   Cursor will prompt you to log in to Supabase. Use the account that owns the project `xfhiwodulrbbtfcqneqt`.
3. **Confirm the MCP server**  
   In Cursor: **Settings → MCP** and ensure the “supabase” server is enabled.

## Optional: change project or mode

Edit `.cursor/mcp.json` and adjust the `url`:

- **Different project:** replace `xfhiwodulrbbtfcqneqt` with your Project ID (Dashboard → Project Settings → General → Reference ID).
- **Allow writes:** remove `&read_only=true` (not recommended for production DBs).
- **Limit tools:** add `&features=database,docs` (or other [feature groups](https://github.com/supabase-community/supabase-mcp#feature-groups)).

Example (read-only, database + docs only):

```json
"url": "https://mcp.supabase.com/mcp?project_ref=xfhiwodulrbbtfcqneqt&read_only=true&features=database,docs"
```

## If Cursor doesn’t see the config

Cursor can prefer the **global** MCP config. If the project-level one is ignored:

1. Open **Cursor Settings → MCP**.
2. Add the same server manually, or copy the contents of `.cursor/mcp.json` into your global `~/.cursor/mcp.json`.

## Security

- **Read-only** is on by default to avoid accidental writes.
- Use a **dev/staging** project when possible; avoid pointing MCP at production with sensitive data.
- Review any suggested SQL or tool calls before approving them in Cursor.

## Links

- [Supabase MCP setup (official)](https://supabase.com/docs/guides/getting-started/mcp)
- [supabase-mcp repo](https://github.com/supabase-community/supabase-mcp)
- [MCP connection in Dashboard](https://supabase.com/dashboard/project/_?showConnect=true&connectTab=mcp)
