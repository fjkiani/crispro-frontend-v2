#!/usr/bin/env python3
"""
Real API Test: ChEMBL Target Extraction
Tests the DynamicFoodExtractor service with real ChEMBL API calls.
"""

import sys
from pathlib import Path
import json
import time
from datetime import datetime
import asyncio

# Add backend to path
project_root = Path(__file__).parent.parent.parent.parent.parent
backend_path = project_root / "oncology-coPilot" / "oncology-backend-minimal"
sys.path.insert(0, str(backend_path))
sys.path.insert(0, str(project_root))  # Add project root for module resolution

from api.services.dynamic_food_extraction import DynamicFoodExtractor

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

# Test cases
TEST_CASES = [
    {
        "name": "vitamin_d",
        "compound": "Vitamin D",
        "expected": {
            "min_targets": 1,
            "should_find": ["VDR"],  # Vitamin D Receptor
            "timeout_seconds": 10
        }
    },
    {
        "name": "curcumin",
        "compound": "Curcumin",
        "expected": {
            "min_targets": 2,
            "should_find": ["PTGS2"],  # COX-2
            "timeout_seconds": 10
        }
    },
    {
        "name": "green_tea",
        "compound": "Green Tea Extract",
        "expected": {
            "min_targets": 1,
            "should_find": [],  # May vary
            "timeout_seconds": 10
        }
    },
    {
        "name": "omega3",
        "compound": "Omega-3 Fatty Acids",
        "expected": {
            "min_targets": 1,
            "should_find": [],
            "timeout_seconds": 10
        }
    },
    {
        "name": "unknown_compound",
        "compound": "FakeCompoundXYZ123",
        "expected": {
            "min_targets": 0,
            "should_find": [],
            "timeout_seconds": 10,
            "should_fail": True
        }
    }
]


async def run_test_case(test_case):
    """
    Execute ChEMBL target extraction for a compound.
    """
    print(f"\n{'='*80}")
    print(f"TEST: {test_case['name']}")
    print(f"Compound: {test_case['compound']}")
    print(f"{'='*80}")
    
    result = {
        "test_name": test_case['name'],
        "compound": test_case['compound'],
        "timestamp": datetime.now().isoformat(),
        "targets_found": [],
        "pathways_found": [],
        "mechanisms_found": [],
        "response_time_ms": None,
        "success": False,
        "errors": []
    }
    
    try:
        # Create fresh extractor instance (no cache between tests)
        extractor = DynamicFoodExtractor()
        
        # Extract targets (await async call)
        start_time = time.time()
        extraction_result = await extractor.extract_all(test_case['compound'])
        response_time = (time.time() - start_time) * 1000
        
        result['response_time_ms'] = round(response_time, 2)
        result['targets_found'] = extraction_result.get('targets', [])
        result['pathways_found'] = extraction_result.get('pathways', [])
        result['mechanisms_found'] = extraction_result.get('mechanisms', [])
        
        # Verification
        num_targets = len(result['targets_found'])
        expected = test_case['expected']
        
        print(f"⏱️  Response Time: {response_time:.2f}ms")
        print(f"🎯 Targets Found: {num_targets}")
        
        if num_targets > 0:
            print(f"   Targets: {', '.join(result['targets_found'][:5])}")
        
        # Check expectations
        if expected.get('should_fail'):
            result['success'] = num_targets == 0
            print(f"   {'✅' if result['success'] else '❌'} Expected to fail (no targets)")
        else:
            meets_minimum = num_targets >= expected['min_targets']
            print(f"   {'✅' if meets_minimum else '❌'} Min targets: {expected['min_targets']} (found {num_targets})")
            
            # Check specific targets
            if expected['should_find']:
                found_expected = [
                    target for target in expected['should_find']
                    if target in result['targets_found']
                ]
                print(f"   Expected targets found: {found_expected}")
            
            result['success'] = meets_minimum
        
        # Log full extraction result
        result['full_extraction'] = extraction_result
    
    except Exception as e:
        result['errors'].append(f"Extraction failed: {str(e)}")
        print(f"❌ Error: {e}")
    
    # Save result
    output_file = OUTPUT_DIR / f"chembl_{test_case['name']}.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    
    print(f"💾 Saved to: {output_file}")
    
    return result


async def main():
    """
    Run all ChEMBL extraction tests.
    """
    print("⚔️ REAL API TEST: ChEMBL TARGET EXTRACTION")
    print("=" * 80)
    print(f"Test Cases: {len(TEST_CASES)}")
    print("=" * 80)
    
    results = []
    for test_case in TEST_CASES:
        result = await run_test_case(test_case)
        results.append(result)
        await asyncio.sleep(1)
    
    # Summary
    print(f"\n{'='*80}")
    print("📊 TEST SUMMARY")
    print(f"{'='*80}")
    
    total = len(results)
    successful = sum(1 for r in results if r['success'])
    total_targets = sum(len(r['targets_found']) for r in results)
    
    # Calculate avg response time safely
    times = [r['response_time_ms'] for r in results if r['response_time_ms']]
    avg_response_time = sum(times) / len(times) if times else 0
    
    print(f"Total Tests: {total}")
    print(f"✅ Successful: {successful}")
    print(f"🎯 Total Targets Extracted: {total_targets}")
    print(f"⏱️  Avg Response Time: {avg_response_time:.2f}ms")
    
    summary_file = OUTPUT_DIR / "chembl_summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            "test_suite": "ChEMBL Target Extraction",
            "timestamp": datetime.now().isoformat(),
            "total_tests": total,
            "successful": successful,
            "total_targets_extracted": total_targets,
            "avg_response_time_ms": round(avg_response_time, 2),
            "results": results
        }, f, indent=2, default=str)
    
    print(f"💾 Summary: {summary_file}")
    
    return 0 if successful == total else 1


if __name__ == "__main__":
    exit(asyncio.run(main()))

