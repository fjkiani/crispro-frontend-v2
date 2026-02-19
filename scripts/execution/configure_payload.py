#!/usr/bin/env python3
"""
PAYLOAD CONFIGURATOR
====================
Extracts credentials from legacy scripts and injects them into the active environment.
Ensures the Ayesha Clinical Intelligence API has access to the Vault.
"""

import os
from pathlib import Path

# Target Credentials (Extracted from explore_project_data_sphere.py)
PDS_USER = "mpm0fxk2"
PDS_PASS = "Thisisfjk12345!"

ENV_PATH = Path(".env")

def inject_credentials():
    print(f"[*] Targeting Environment File: {ENV_PATH.absolute()}")
    
    current_content = ""
    if ENV_PATH.exists():
        current_content = ENV_PATH.read_text()
    
    # Check if already present
    if "PDS_API_KEY" in current_content or "SAS_USERNAME" in current_content:
        print("[!] Credentials already exist in .env. Skipping injection.")
        return

    # Injection Payload
    payload = f"\n# --- MARS DEPLOYMENT CREDENTIALS ---\nSAS_USERNAME={PDS_USER}\nSAS_PASSWORD={PDS_PASS}\nPDS_API_KEY={PDS_PASS} # Mapping pass to key for backward compat\n"
    
    try:
        with open(ENV_PATH, "a") as f:
            f.write(payload)
        print("[+] Credentials successfully injected into .env")
        print(f"    User: {PDS_USER}")
        print("    Pass: ************")
    except Exception as e:
        print(f"[-] Injection Failed: {e}")

if __name__ == "__main__":
    inject_credentials()
