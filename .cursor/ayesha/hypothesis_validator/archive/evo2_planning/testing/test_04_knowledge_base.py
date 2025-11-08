#!/usr/bin/env python3
"""
Test 4: Knowledge Base Extraction

Goal: Verify we can extract targets from food_targets.json
"""

import json
import os
from pathlib import Path

def test_knowledge_base():
    """Test knowledge base target extraction."""
    
    print("🧪 TEST 4: Knowledge Base Extraction\n")
    
    # Find food_targets.json
    base_path = Path(__file__).parent.parent.parent.parent
    food_targets_path = base_path / ".cursor/ayesha/hypothesis_validator/data/food_targets.json"
    
    print(f"📁 Looking for: {food_targets_path}")
    
    if not food_targets_path.exists():
        print(f"   ❌ File not found!")
        print(f"   📝 Need to create food_targets.json first")
        return False
    
    print(f"   ✅ File exists\n")
    
    try:
        with open(food_targets_path) as f:
            food_targets = json.load(f)
        
        print(f"📋 Knowledge Base Contents:")
        print(f"   Compounds: {len(food_targets.get('compounds', []))}")
        print()
        
        # Test extraction for known compounds
        test_compounds = ["Vitamin D", "Curcumin", "NAC"]
        
        for compound in test_compounds:
            print(f"🔍 Looking for: {compound}")
            
            found = False
            for item in food_targets.get('compounds', []):
                if compound.lower() in item.get('compound', '').lower():
                    targets = item.get('targets', [])
                    pathways = item.get('pathways', [])
                    
                    print(f"   ✅ Found!")
                    print(f"   Targets: {targets}")
                    print(f"   Pathways: {pathways}")
                    found = True
                    break
            
            if not found:
                print(f"   ⚠️  Not found in knowledge base")
            
            print()
        
        print("✅ TEST 4 COMPLETE: Knowledge base extraction works")
        return True
        
    except json.JSONDecodeError as e:
        print(f"   ❌ Invalid JSON: {e}")
        return False
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return False

if __name__ == "__main__":
    success = test_knowledge_base()
    if not success:
        print("\n⚠️  TEST 4 FAILED: Need to fix knowledge base before proceeding")

