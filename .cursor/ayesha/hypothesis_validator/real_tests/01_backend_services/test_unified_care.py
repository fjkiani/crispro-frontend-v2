#!/usr/bin/env python3
"""
Real API Test: Unified Care Plan
Tests the /api/ayesha/complete_care_plan endpoint with real HTTP calls.
Verifies parallel execution and integrated confidence.
"""

import requests
import json
import time
from datetime import datetime
from pathlib import Path

# Backend URL
BASE_URL = "http://127.0.0.1:8000"
ENDPOINT = f"{BASE_URL}/api/ayesha/complete_care_plan"

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

# Test cases
TEST_CASES = [
    {
        "name": "ayesha_hrd_positive_line2",
        "description": "HRD+ patient, Line 2, after platinum failure",
        "input": {
            "patient_context": {
                "disease": "ovarian_cancer_hgs",
                "treatment_history": [
                    {
                        "line": 1,
                        "drugs": ["Carboplatin", "Paclitaxel"],
                        "outcome": "partial_response"
                    },
                    {
                        "line": 2,
                        "drugs": ["Olaparib"],
                        "outcome": "progression"
                    }
                ],
                "biomarkers": {
                    "brca1_mutant": False,
                    "brca2_mutant": False,
                    "hrd_positive": True,
                    "tp53_mutant": True
                }
            }
        },
        "verify": [
            "drug_recommendations array present",
            "food_recommendations array present",
            "integrated_confidence is number",
            "provenance tracks both drug and food",
            "timestamp difference < 1 second (parallel execution)"
        ]
    },
    {
        "name": "brca1_germline_treatment_naive",
        "description": "BRCA1 germline positive, treatment-naive",
        "input": {
            "patient_context": {
                "disease": "ovarian_cancer_hgs",
                "treatment_history": [],
                "biomarkers": {
                    "brca1_mutant": True,
                    "brca2_mutant": False,
                    "hrd_positive": False,
                    "tp53_mutant": False
                }
            }
        },
        "verify": [
            "drug recommendations include PARP inhibitors",
            "cross_resistance is LOW (no priors)",
            "food recommendations appropriate for Line 1"
        ]
    }
]


def verify_parallel_execution(provenance):
    """
    Check if drug and food APIs were called in parallel.
    Returns True if timestamp difference < 1 second.
    """
    try:
        if 'drug_analysis' in provenance and 'food_analysis' in provenance:
            drug_ts = datetime.fromisoformat(provenance['drug_analysis'].get('timestamp', ''))
            food_ts = datetime.fromisoformat(provenance['food_analysis'].get('timestamp', ''))
            
            time_diff = abs((drug_ts - food_ts).total_seconds())
            
            return {
                "drug_timestamp": provenance['drug_analysis']['timestamp'],
                "food_timestamp": provenance['food_analysis']['timestamp'],
                "time_difference_seconds": round(time_diff, 3),
                "is_parallel": time_diff < 1.0,
                "execution_mode": "parallel" if time_diff < 1.0 else "sequential"
            }
    except Exception as e:
        return {
            "error": f"Could not verify parallel execution: {str(e)}",
            "is_parallel": None
        }
    
    return {"is_parallel": None, "error": "Missing provenance data"}


def run_test_case(test_case):
    """
    Execute a single test case and capture results.
    """
    print(f"\n{'='*80}")
    print(f"TEST: {test_case['name']}")
    print(f"DESC: {test_case['description']}")
    print(f"{'='*80}")
    
    result = {
        "test_name": test_case['name'],
        "description": test_case['description'],
        "timestamp": datetime.now().isoformat(),
        "input": test_case['input'],
        "output": None,
        "response_time_ms": None,
        "status_code": None,
        "parallel_execution_check": None,
        "success": False,
        "errors": []
    }
    
    try:
        # Make real API call
        start_time = time.time()
        response = requests.post(
            ENDPOINT,
            json=test_case['input'],
            headers={'Content-Type': 'application/json'},
            timeout=60  # Longer timeout for unified call
        )
        response_time = (time.time() - start_time) * 1000
        
        result['response_time_ms'] = round(response_time, 2)
        result['status_code'] = response.status_code
        
        if response.status_code == 200:
            response_data = response.json()
            result['output'] = response_data
            
            # Verify structure
            has_drugs = 'drug_recommendations' in response_data
            has_foods = 'food_recommendations' in response_data
            has_confidence = 'integrated_confidence' in response_data
            has_provenance = 'provenance' in response_data
            
            print(f"✅ Status Code: {response.status_code}")
            print(f"⏱️  Response Time: {response_time:.2f}ms")
            print(f"📊 Structure Check:")
            print(f"   {'✅' if has_drugs else '❌'} drug_recommendations: {len(response_data.get('drug_recommendations', []))} drugs")
            print(f"   {'✅' if has_foods else '❌'} food_recommendations: {len(response_data.get('food_recommendations', []))} foods")
            print(f"   {'✅' if has_confidence else '❌'} integrated_confidence: {response_data.get('integrated_confidence', 'N/A')}")
            print(f"   {'✅' if has_provenance else '❌'} provenance: present")
            
            # Check parallel execution
            if has_provenance:
                parallel_check = verify_parallel_execution(response_data['provenance'])
                result['parallel_execution_check'] = parallel_check
                
                if parallel_check.get('is_parallel'):
                    print(f"⚡ Parallel Execution: YES ({parallel_check['time_difference_seconds']}s)")
                elif parallel_check.get('is_parallel') is False:
                    print(f"⚠️  Sequential Execution: {parallel_check['time_difference_seconds']}s")
                else:
                    print(f"❓ Could not verify parallel execution")
            
            # Check SAE features in drug recommendations
            if has_drugs:
                drugs_with_sae = sum(
                    1 for drug in response_data['drug_recommendations']
                    if 'sae_features' in drug
                )
                print(f"   SAE features in drugs: {drugs_with_sae}/{len(response_data['drug_recommendations'])}")
            
            result['success'] = has_drugs and has_foods and has_confidence
        
        else:
            result['errors'].append(f"HTTP {response.status_code}: {response.text}")
            print(f"❌ HTTP Error: {response.status_code}")
    
    except requests.exceptions.Timeout:
        result['errors'].append("Request timeout (>60s)")
        print(f"❌ Timeout")
    
    except Exception as e:
        result['errors'].append(f"Error: {str(e)}")
        print(f"❌ Error: {e}")
    
    # Save result
    output_file = OUTPUT_DIR / f"unified_care_{test_case['name']}.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    
    print(f"💾 Saved to: {output_file}")
    
    return result


def main():
    """
    Run all test cases.
    """
    print("⚔️ REAL API TEST: UNIFIED CARE PLAN")
    print("=" * 80)
    print(f"Backend: {BASE_URL}")
    print(f"Endpoint: {ENDPOINT}")
    print("=" * 80)
    
    # Check backend health
    try:
        health_response = requests.get(f"{BASE_URL}/health", timeout=5)
        if health_response.status_code == 200:
            print("✅ Backend is UP")
        else:
            print(f"⚠️  Backend health check returned {health_response.status_code}")
    except:
        print("❌ Backend is DOWN")
        return 1
    
    # Run tests
    results = []
    for test_case in TEST_CASES:
        result = run_test_case(test_case)
        results.append(result)
        time.sleep(1)
    
    # Summary
    print(f"\n{'='*80}")
    print("📊 TEST SUMMARY")
    print(f"{'='*80}")
    
    total = len(results)
    successful = sum(1 for r in results if r['success'])
    parallel_confirmed = sum(
        1 for r in results
        if r.get('parallel_execution_check', {}).get('is_parallel')
    )
    
    print(f"Total Tests: {total}")
    print(f"✅ Successful: {successful}")
    print(f"⚡ Parallel Execution Confirmed: {parallel_confirmed}/{total}")
    
    summary_file = OUTPUT_DIR / "unified_care_summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            "test_suite": "Unified Care Plan API",
            "timestamp": datetime.now().isoformat(),
            "total_tests": total,
            "successful": successful,
            "parallel_execution_confirmed": parallel_confirmed,
            "results": results
        }, f, indent=2, default=str)
    
    print(f"💾 Summary: {summary_file}")
    
    return 0 if successful == total else 1


if __name__ == "__main__":
    exit(main())






