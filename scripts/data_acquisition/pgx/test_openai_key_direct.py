#!/usr/bin/env python3
"""
Direct test of OpenAI API key (without BioMed-MCP)

This script tests if the OpenAI API key is valid by making a simple API call.
"""
import os
import sys
from pathlib import Path
from dotenv import load_dotenv
from openai import OpenAI

# Load environment variables
biomed_mcp_path = Path(__file__).parent.parent / "mcp_servers" / "BioMed-MCP"
env_file = biomed_mcp_path / ".env"
if env_file.exists():
    load_dotenv(env_file)
    print(f"✅ Loaded .env from: {env_file}")
else:
    print(f"❌ .env file not found at: {env_file}")
    sys.exit(1)

# Get API key
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    print("❌ OPENAI_API_KEY not found in environment")
    sys.exit(1)

print(f"✅ API key found: {api_key[:20]}...{api_key[-10:]}")

# Test the API key directly
print("\n--- Testing OpenAI API key directly ---")
try:
    client = OpenAI(api_key=api_key)
    
    # Make a simple test call
    response = client.chat.completions.create(
        model="gpt-4o-mini",  # Use cheaper model for testing
        messages=[
            {"role": "user", "content": "Say 'API key is valid' if you can read this."}
        ],
        max_tokens=10
    )
    
    print("✅ API key is VALID!")
    print(f"   Response: {response.choices[0].message.content}")
    print("\n✅ You can now use BioMed-MCP with this key!")
    
except Exception as e:
    print(f"❌ API key test FAILED: {e}")
    print("\n⚠️  Please check:")
    print("   1. API key is correct (no extra spaces/newlines)")
    print("   2. API key hasn't expired")
    print("   3. You have credits/quota on your OpenAI account")
    print("   4. Get a new key from: https://platform.openai.com/account/api-keys")
    sys.exit(1)



