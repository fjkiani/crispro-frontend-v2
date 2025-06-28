#!/usr/bin/env python3
"""
Test script for the Zeta Oracle Service deployment.
"""

import os
import requests
import json
import time
import sys

def test_local_deployment():
    """Test the local deployment using modal run"""
    print("🧪 Testing Zeta Oracle Local Deployment...")
    print("Run: modal run deploy_zeta_oracle.py::test_local")
    print("This will test the model locally before deployment.")

def test_api_endpoint(endpoint_url):
    """Test the deployed API endpoint"""
    print(f"🌐 Testing Zeta Oracle API Endpoint: {endpoint_url}")
    
    # Test payload
    test_payload = {
        "reference_sequence": "ATCGATCGATCGATCG",
        "alternate_sequence": "ATCGATCGATCGTTCG"  # Single nucleotide change T->T
    }
    
    headers = {"Content-Type": "application/json"}
    
    try:
        print(f"📤 Sending request to: {endpoint_url}")
        print(f"📦 Payload: {json.dumps(test_payload, indent=2)}")
        
        response = requests.post(
            endpoint_url,
            json=test_payload,
            headers=headers,
            timeout=30
        )
        
        print(f"📥 Response Status: {response.status_code}")
        print(f"📄 Response Headers: {dict(response.headers)}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"✅ Success! Response: {json.dumps(result, indent=2)}")
            
            # Interpret the result
            if "delta_likelihood_score" in result:
                score = result["delta_likelihood_score"]
                if score < -0.1:
                    interpretation = "Likely pathogenic"
                elif score > 0.1:
                    interpretation = "Likely benign"
                else:
                    interpretation = "Uncertain significance"
                    
                print(f"🔬 Interpretation: {interpretation}")
                print(f"📊 Delta Score: {score}")
            
            return True
        else:
            print(f"❌ Error: {response.status_code}")
            print(f"📝 Response: {response.text}")
            return False
            
    except requests.exceptions.Timeout:
        print("⏰ Request timed out - the model might be cold starting")
        return False
    except requests.exceptions.ConnectionError:
        print("🔌 Connection error - check if the endpoint is deployed")
        return False
    except Exception as e:
        print(f"💥 Unexpected error: {e}")
        return False

def get_endpoint_url():
    """Get the endpoint URL from environment or user input"""
    endpoint_url = os.environ.get("ZETA_ORACLE_ENDPOINT")
    
    if not endpoint_url:
        print("🔍 ZETA_ORACLE_ENDPOINT environment variable not found.")
        print("Please provide the endpoint URL:")
        endpoint_url = input("Endpoint URL: ").strip()
    
    if not endpoint_url:
        print("❌ No endpoint URL provided")
        return None
        
    return endpoint_url

def main():
    """Main test function"""
    print("🚀 Zeta Oracle Service Test Suite")
    print("=" * 50)
    
    # Test local deployment first
    print("\n1. Local Deployment Test")
    print("-" * 25)
    test_local_deployment()
    
    # Test API endpoint
    print("\n2. API Endpoint Test")
    print("-" * 20)
    
    endpoint_url = get_endpoint_url()
    if endpoint_url:
        success = test_api_endpoint(endpoint_url)
        if success:
            print("\n🎉 All tests passed!")
        else:
            print("\n❌ API endpoint test failed")
            sys.exit(1)
    else:
        print("⚠️  Skipping API endpoint test - no URL provided")
    
    print("\n" + "=" * 50)
    print("🏁 Test suite completed")

if __name__ == "__main__":
    main() 