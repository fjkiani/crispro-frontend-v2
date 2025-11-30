#!/bin/bash
# ⚔️ ZO AUTONOMOUS MONITORING SCRIPT ⚔️
# Monitors seeding progress and updates JR2 sync files

SYNC_FILE=".cursor/concept/jr2_dossier_mission/00_ZO_JR2_SYNC.json"
ITERATION_LOG=".cursor/concept/jr2_dossier_mission/01_ZO_ITERATION_LOG.md"
SEEDING_LOG="/tmp/seeding_log.txt"

echo "🔥 ZO AUTONOMOUS MONITOR STARTED"
echo "Monitoring seeding progress every 60 seconds..."

while true; do
    # Check if seeding is still running
    if pgrep -f "agent_1_seeding.main" > /dev/null; then
        # Count trials in SQLite
        TRIAL_COUNT=$(sqlite3 /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend/data/clinical_trials.db "SELECT COUNT(*) FROM clinical_trials;" 2>/dev/null || echo "0")
        
        echo "[$(date '+%H:%M:%S')] 🔄 Seeding in progress... Trials: $TRIAL_COUNT"
        
        # Update sync file (simple append for now)
        echo "[$(date '+%H:%M:%S')] Trials seeded: $TRIAL_COUNT" >> .cursor/concept/jr2_dossier_mission/ZO_PROGRESS.log
        
        sleep 60
    else
        echo "[$(date '+%H:%M:%S')] ✅ SEEDING COMPLETE!"
        
        # Get final count
        TRIAL_COUNT=$(sqlite3 /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend/data/clinical_trials.db "SELECT COUNT(*) FROM clinical_trials;" 2>/dev/null || echo "0")
        
        echo "Final trial count: $TRIAL_COUNT"
        echo "Check $SEEDING_LOG for full details"
        
        # Trigger next iteration (AstraDB seeding)
        echo ""
        echo "🎯 READY FOR ITERATION 2: AstraDB Seeding"
        echo "Run: cd oncology-coPilot/oncology-backend-minimal && python3 scripts/seed_astradb_from_sqlite.py --limit $TRIAL_COUNT"
        
        break
    fi
done

echo "✅ MONITOR COMPLETE"

