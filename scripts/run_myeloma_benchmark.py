#!/usr/bin/env python3
import json
import os
import sys
import urllib.request

BACKEND = os.environ.get('BACKEND_URL', 'https://crispro-oncology-backend-minimal.vercel.app')
PANEL = os.environ.get('PANEL_PATH', 'docs/benchmarks/myeloma_panel.json')
MODEL = os.environ.get('MODEL_ID', 'evo2_7b')

def post_json(url: str, payload: dict, timeout: int = 600) -> dict:
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode('utf-8')
        return json.loads(body)

def main():
    with open(PANEL, 'r') as f:
        mutations = json.load(f)
    url = f"{BACKEND.rstrip('/')}/api/predict/myeloma_drug_response"
    payload = { 'model_id': MODEL, 'mutations': mutations }
    result = post_json(url, payload)
    print(json.dumps({
        'prediction': result.get('prediction'),
        'pathway_scores': result.get('pathway_scores'),
        'mode': result.get('mode'),
        'upstream_service': result.get('upstream_service')
    }, indent=2))

if __name__ == '__main__':
    main() 