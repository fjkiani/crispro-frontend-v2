# 🎯 MODULE 6: MAIN CLI

## **Purpose**
Orchestrate all modules, handle CLI arguments

## **File Location**
`oncology-coPilot/oncology-backend/scripts/agent_1_seeding/main.py`

## **Key Functions**

### **`parse_args() -> argparse.Namespace`**
CLI argument parser:
- `--limit` (default: 1000)
- `--test-mode` (sets limit=100, skip_embeddings=True)
- `--skip-embeddings`
- `--skip-migration`

### **`main()`**
Main orchestration function:
1. Parse CLI arguments
2. Run schema migration (unless skipped)
3. Fetch trials from API
4. Parse trials
5. Insert into SQLite
6. Embed into ChromaDB (unless skipped)
7. Generate summary report

### **`generate_summary_report(args: argparse.Namespace) -> None`**
Summary report generator:
- Total ovarian trials
- With biomarker requirements (%)
- With locations data (%)
- With NY locations (%)
- ChromaDB embeddings count (if not skipped)

## **Features**
- ✅ CLI arguments: `--limit`, `--test-mode`, `--skip-embeddings`, `--skip-migration`
- ✅ Summary report with stats
- ✅ Test mode support (100 trials, no embeddings)

## **Dependencies**
- All other modules (Config, API, Parsers, Database, Utils)

## **Acceptance Criteria**
- [ ] CLI arguments work correctly
- [ ] Script runs end-to-end without errors
- [ ] Completes in <15 minutes (1000 trials)
- [ ] Summary report generated correctly

## **Time Estimate:** 30 minutes


