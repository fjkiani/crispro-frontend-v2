#!/usr/bin/env python3
"""
Remove restored conversations from composerData to see if that makes them visible.
Cursor might only show conversations that are NOT in composerData.
"""

import json
import sqlite3
from pathlib import Path

def remove_restored_from_composerdata():
    """
    Remove restored conversations from composerData.
    Cursor might load conversations directly from bubbles, not from composerData.
    """
    current_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print("🔧 Removing restored conversations from composerData...\n")
    
    if not Path(current_db).exists():
        print(f"❌ Current DB not found!")
        return False
    
    conn = sqlite3.connect(current_db)
    cursor = conn.cursor()
    
    # Get composerData
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    composer_row = cursor.fetchone()
    
    if not composer_row:
        print("❌ No composerData found!")
        conn.close()
        return False
    
    composer_data = json.loads(composer_row[0])
    all_composers = composer_data.get('allComposers', [])
    
    print(f"📊 Total conversations in composerData: {len(all_composers)}")
    
    # Keep only non-restored conversations
    visible_composers = [c for c in all_composers if 'Restored' not in c.get('title', '')]
    removed_count = len(all_composers) - len(visible_composers)
    
    print(f"📊 Keeping {len(visible_composers)} visible conversations")
    print(f"🗑️  Removing {removed_count} restored conversations from composerData\n")
    
    # Update composerData
    composer_data['allComposers'] = visible_composers
    cursor.execute(
        "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
        ('composer.composerData', json.dumps(composer_data))
    )
    conn.commit()
    conn.close()
    
    print("✅ Removed restored conversations from composerData!")
    print("\n🔄 Restart Cursor - conversations should now appear since they only have bubbles (like visible ones)\n")
    return True

if __name__ == "__main__":
    remove_restored_from_composerdata()



