#!/usr/bin/env python3
"""
Real API Test: Food Validator
Tests the /api/hypothesis/validate_food_ab_enhanced endpoint with real HTTP calls.
NO MOCKS - captures actual input/output.
"""

import requests
import json
import time
from datetime import datetime
from pathlib import Path

# Backend URL
BASE_URL = "http://127.0.0.1:8000"
ENDPOINT = f"{BASE_URL}/api/hypothesis/validate_food_ab_enhanced"

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

# Test cases
TEST_CASES = [
    {
        "name": "vitamin_d_line3",
        "description": "Vitamin D for Line 3 ovarian cancer with platinum resistance",
        "input": {
            "compound": "Vitamin D",
            "disease": "ovarian_cancer_hgs",
            "treatment_line": 3,
            "prior_therapies": ["carboplatin", "paclitaxel"],
            "use_llm": True
        },
        "expected_structure": {
            "status": str,
            "sae_features": dict,
            "provenance": dict,
            "targets": list,
            "dosage": dict
        }
    },
    {
        "name": "curcumin_line2",
        "description": "Curcumin for Line 2 without LLM",
        "input": {
            "compound": "Curcumin",
            "disease": "ovarian_cancer_hgs",
            "treatment_line": 2,
            "prior_therapies": ["carboplatin"],
            "use_llm": False
        },
        "expected_structure": {
            "status": str,
            "confidence": (float, int),
            "sae_features": dict
        }
    },
    {
        "name": "omega3_line1",
        "description": "Omega-3 for treatment-naive Line 1",
        "input": {
            "compound": "Omega-3 Fatty Acids",
            "disease": "ovarian_cancer_hgs",
            "treatment_line": 1,
            "prior_therapies": [],
            "use_llm": False
        },
        "expected_structure": {
            "status": str,
            "sae_features": dict
        }
    }
]


def verify_structure(response_data, expected_structure):
    """
    Verify the response has expected keys and types.
    Returns dict of verification results.
    """
    results = {}
    
    for key, expected_type in expected_structure.items():
        if key not in response_data:
            results[key] = {
                "present": False,
                "type_correct": False,
                "error": f"Missing key: {key}"
            }
        else:
            value = response_data[key]
            if isinstance(expected_type, tuple):
                type_correct = isinstance(value, expected_type)
            else:
                type_correct = isinstance(value, expected_type)
            
            results[key] = {
                "present": True,
                "type_correct": type_correct,
                "actual_type": type(value).__name__,
                "expected_type": str(expected_type)
            }
    
    return results


def run_test_case(test_case):
    """
    Execute a single test case and capture results.
    """
    print(f"\n{'='*80}")
    print(f"TEST: {test_case['name']}")
    print(f"DESC: {test_case['description']}")
    print(f"{'='*80}")
    
    # Prepare result structure
    result = {
        "test_name": test_case['name'],
        "description": test_case['description'],
        "timestamp": datetime.now().isoformat(),
        "input": test_case['input'],
        "output": None,
        "response_time_ms": None,
        "status_code": None,
        "structure_verification": None,
        "success": False,
        "errors": []
    }
    
    try:
        # Make real API call (using query params, not JSON body)
        start_time = time.time()
        response = requests.post(
            ENDPOINT,
            params=test_case['input'],  # Use params instead of json
            timeout=30
        )
        response_time = (time.time() - start_time) * 1000  # Convert to ms
        
        result['response_time_ms'] = round(response_time, 2)
        result['status_code'] = response.status_code
        
        # Parse response
        if response.status_code == 200:
            response_data = response.json()
            result['output'] = response_data
            
            # Verify structure
            if 'expected_structure' in test_case:
                result['structure_verification'] = verify_structure(
                    response_data,
                    test_case['expected_structure']
                )
            
            # Check if all required fields present
            all_present = all(
                v['present'] for v in result['structure_verification'].values()
            )
            all_types_correct = all(
                v['type_correct'] for v in result['structure_verification'].values()
            )
            
            result['success'] = all_present and all_types_correct
            
            # Log details
            print(f"✅ Status Code: {response.status_code}")
            print(f"⏱️  Response Time: {response_time:.2f}ms")
            print(f"📊 Structure Verification:")
            for key, verification in result['structure_verification'].items():
                status = "✅" if verification['present'] and verification['type_correct'] else "❌"
                print(f"   {status} {key}: {verification}")
            
        else:
            result['errors'].append(f"HTTP {response.status_code}: {response.text}")
            print(f"❌ HTTP Error: {response.status_code}")
            print(f"Response: {response.text}")
    
    except requests.exceptions.Timeout:
        result['errors'].append("Request timeout (>30s)")
        print(f"❌ Timeout")
    
    except requests.exceptions.ConnectionError:
        result['errors'].append("Connection refused - is backend running?")
        print(f"❌ Connection Error - Backend not running?")
    
    except Exception as e:
        result['errors'].append(f"Unexpected error: {str(e)}")
        print(f"❌ Error: {e}")
    
    # Save result to file
    output_file = OUTPUT_DIR / f"food_validator_{test_case['name']}.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    
    print(f"💾 Saved to: {output_file}")
    
    return result


def main():
    """
    Run all test cases and generate summary.
    """
    print("⚔️ REAL API TEST: FOOD VALIDATOR")
    print("=" * 80)
    print(f"Backend: {BASE_URL}")
    print(f"Endpoint: {ENDPOINT}")
    print(f"Test Cases: {len(TEST_CASES)}")
    print("=" * 80)
    
    # Check backend health
    try:
        health_response = requests.get(f"{BASE_URL}/health", timeout=5)
        if health_response.status_code == 200:
            print("✅ Backend is UP")
        else:
            print(f"⚠️  Backend health check returned {health_response.status_code}")
    except:
        print("❌ Backend is DOWN - start it first!")
        print("   cd oncology-coPilot/oncology-backend-minimal")
        print("   source venv/bin/activate")
        print("   uvicorn api.main:app --reload")
        return
    
    # Run all tests
    results = []
    for test_case in TEST_CASES:
        result = run_test_case(test_case)
        results.append(result)
        time.sleep(1)  # Brief pause between tests
    
    # Summary
    print(f"\n{'='*80}")
    print("📊 TEST SUMMARY")
    print(f"{'='*80}")
    
    total = len(results)
    successful = sum(1 for r in results if r['success'])
    failed = total - successful
    avg_response_time = sum(
        r['response_time_ms'] for r in results if r['response_time_ms']
    ) / len([r for r in results if r['response_time_ms']])
    
    print(f"Total Tests: {total}")
    print(f"✅ Successful: {successful}")
    print(f"❌ Failed: {failed}")
    print(f"⏱️  Avg Response Time: {avg_response_time:.2f}ms")
    print(f"\nOutput files in: {OUTPUT_DIR}")
    
    # Save summary
    summary = {
        "test_suite": "Food Validator API",
        "timestamp": datetime.now().isoformat(),
        "total_tests": total,
        "successful": successful,
        "failed": failed,
        "avg_response_time_ms": round(avg_response_time, 2),
        "results": results
    }
    
    summary_file = OUTPUT_DIR / "food_validator_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print(f"💾 Summary saved to: {summary_file}")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    exit(main())

