#!/usr/bin/env python3
"""
Quick utility to inspect the format of existing chat bubbles in Cursor's database.
This helps us understand the exact JSON structure Cursor expects.

Usage:
    python inspect_bubble_format.py [--db-path PATH]
"""

import sqlite3
import json
import argparse
from pathlib import Path

def inspect_database(db_path):
    """Inspect existing bubble entries to understand the format."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    print("🔍 Inspecting Cursor database...")
    print(f"Database: {db_path}\n")
    
    # Get a sample bubble entry
    cursor.execute("""
        SELECT key, value FROM cursorDiskKV 
        WHERE key LIKE 'bubbleId:%'
        LIMIT 5
    """)
    
    bubbles = cursor.fetchall()
    
    if not bubbles:
        print("❌ No bubble entries found in database")
        return
    
    print(f"📊 Found {len(bubbles)} sample bubbles\n")
    
    for i, (key, value) in enumerate(bubbles, 1):
        print(f"{'='*80}")
        print(f"BUBBLE #{i}")
        print(f"{'='*80}")
        print(f"Key: {key}")
        print()
        
        # Decode value
        if isinstance(value, bytes):
            value_str = value.decode('utf-8')
        else:
            value_str = str(value)
        
        # Try to parse as JSON
        try:
            bubble_json = json.loads(value_str)
            print("✅ Valid JSON structure:")
            print(json.dumps(bubble_json, indent=2, ensure_ascii=False))
            print()
            print("📋 Field summary:")
            for field, val in bubble_json.items():
                field_type = type(val).__name__
                if isinstance(val, str) and len(val) > 50:
                    val_preview = val[:50] + "..."
                else:
                    val_preview = val
                print(f"  - {field}: {field_type} = {val_preview}")
        except json.JSONDecodeError as e:
            print(f"❌ Invalid JSON: {e}")
            print(f"Raw value (first 500 chars): {value_str[:500]}")
        
        print()
    
    # Check composerData format
    print(f"{'='*80}")
    print("COMPOSER DATA ENTRIES")
    print(f"{'='*80}")
    cursor.execute("""
        SELECT key, value FROM cursorDiskKV 
        WHERE key LIKE 'composerData:%'
        LIMIT 3
    """)
    
    composers = cursor.fetchall()
    
    if composers:
        print(f"📊 Found {len(composers)} composerData entries\n")
        for key, value in composers:
            print(f"Key: {key}")
            if isinstance(value, bytes):
                value_str = value.decode('utf-8')
            else:
                value_str = str(value)
            try:
                composer_json = json.loads(value_str)
                print(json.dumps(composer_json, indent=2))
            except:
                print(f"Raw: {value_str[:200]}")
            print()
    else:
        print("⚠️  No composerData entries found\n")
    
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Inspect Cursor bubble format')
    parser.add_argument(
        '--db-path',
        default=Path.home() / 'Library/Application Support/Cursor/User/globalStorage/state.vscdb',
        help='Path to Cursor database file'
    )
    
    args = parser.parse_args()
    
    db_path = Path(args.db_path)
    if not db_path.exists():
        print(f"❌ Database not found: {db_path}")
        exit(1)
    
    inspect_database(db_path)






