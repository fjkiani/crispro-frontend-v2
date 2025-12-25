#!/usr/bin/env python3
"""
Reusable cBioPortal API Client
===============================
Modular utility for interacting with cBioPortal API.
Can be reused for any future data acquisition tasks.

Author: Agent
Date: January 28, 2025
"""

import httpx
import time
import json
from typing import List, Dict, Optional, Any
from collections import defaultdict

CBIO_BASE = "https://www.cbioportal.org/api"
DEFAULT_TIMEOUT = 120.0
RATE_LIMIT_DELAY = 0.5  # seconds between requests


class CBioportalClient:
    """Reusable cBioPortal API client with rate limiting and error handling."""
    
    def __init__(self, base_url: str = CBIO_BASE, timeout: float = DEFAULT_TIMEOUT, rate_limit_delay: float = RATE_LIMIT_DELAY):
        self.base_url = base_url
        self.timeout = timeout
        self.rate_limit_delay = rate_limit_delay
        self._mutation_profile_cache: Dict[str, str] = {}  # study_id -> profile_id
        self._client = None  # Persistent httpx client
    
    def _get_client(self) -> httpx.Client:
        """Get or create persistent httpx client."""
        if self._client is None:
            self._client = httpx.Client(timeout=self.timeout, follow_redirects=True)
        return self._client
    
    def _make_request(self, method: str, endpoint: str, **kwargs) -> Any:
        """Make HTTP request with retry logic and rate limiting."""
        url = f"{self.base_url}{endpoint}"
        
        # Rate limiting
        time.sleep(self.rate_limit_delay)
        
        max_retries = 3
        for attempt in range(max_retries):
            try:
                client = self._get_client()
                
                if method.upper() == "GET":
                    response = client.get(url, **kwargs)
                elif method.upper() == "POST":
                    response = client.post(url, **kwargs)
                else:
                    raise ValueError(f"Unsupported method: {method}")
                
                response.raise_for_status()
                
                # Check if response is empty
                if not response.text or not response.text.strip():
                    if attempt < max_retries - 1:
                        wait_time = (2 ** attempt) * self.rate_limit_delay
                        print(f"  ⚠️  Empty response. Retrying in {wait_time:.1f}s...")
                        time.sleep(wait_time)
                        continue
                    raise Exception(f"Empty response from {endpoint}")
                
                # Try to parse JSON
                try:
                    return response.json()
                except json.JSONDecodeError as e:
                    # Log first 200 chars of response for debugging
                    preview = response.text[:200] if response.text else "(empty)"
                    if attempt < max_retries - 1:
                        wait_time = (2 ** attempt) * self.rate_limit_delay
                        print(f"  ⚠️  Invalid JSON response (status {response.status_code}). Retrying in {wait_time:.1f}s...")
                        print(f"      Response preview: {preview}")
                        time.sleep(wait_time)
                        continue
                    raise Exception(f"Invalid JSON response from {endpoint}: {preview}")
            
            except httpx.HTTPStatusError as e:
                if e.response.status_code == 429:  # Rate limited
                    wait_time = (2 ** attempt) * self.rate_limit_delay  # Exponential backoff
                    print(f"  ⚠️  Rate limited. Waiting {wait_time:.1f}s before retry {attempt + 1}/{max_retries}...")
                    time.sleep(wait_time)
                    continue
                # Log error response
                error_preview = e.response.text[:200] if e.response.text else "(empty)"
                print(f"  ⚠️  HTTP {e.response.status_code}: {error_preview}")
                raise
            except Exception as e:
                if attempt < max_retries - 1:
                    wait_time = (2 ** attempt) * self.rate_limit_delay
                    print(f"  ⚠️  Request failed: {e}. Retrying in {wait_time:.1f}s...")
                    time.sleep(wait_time)
                    continue
                raise
        
        raise Exception(f"Failed after {max_retries} attempts")
    
    def close(self):
        """Close the persistent client."""
        if self._client:
            self._client.close()
            self._client = None
    
    def get_study_patients(self, study_id: str) -> List[Dict]:
        """Fetch all patients for a study."""
        endpoint = f"/studies/{study_id}/patients"
        params = {"pageSize": 10000}
        return self._make_request("GET", endpoint, params=params)
    
    def get_study_samples(self, study_id: str) -> List[Dict]:
        """Fetch all samples for a study."""
        endpoint = f"/studies/{study_id}/samples"
        params = {"pageSize": 10000}
        return self._make_request("GET", endpoint, params=params)
    
    def get_molecular_profiles(self, study_id: str) -> List[Dict]:
        """Fetch all molecular profiles for a study."""
        endpoint = f"/studies/{study_id}/molecular-profiles"
        return self._make_request("GET", endpoint)
    
    def get_mutation_profile_id(self, study_id: str) -> Optional[str]:
        """Get mutation profile ID for a study (cached)."""
        if study_id in self._mutation_profile_cache:
            return self._mutation_profile_cache[study_id]
        
        profiles = self.get_molecular_profiles(study_id)
        for profile in profiles:
            profile_id = profile.get("molecularProfileId", "")
            if profile_id.endswith("_mutations"):
                self._mutation_profile_cache[study_id] = profile_id
                return profile_id
        
        # Fallback: try pattern matching
        fallback_id = f"{study_id}_mutations"
        print(f"  ⚠️  Mutation profile not found, using fallback: {fallback_id}")
        self._mutation_profile_cache[study_id] = fallback_id
        return fallback_id
    
    def get_mutations_for_samples(self, profile_id: str, sample_ids: List[str]) -> List[Dict]:
        """Fetch mutations for given sample IDs."""
        endpoint = f"/molecular-profiles/{profile_id}/mutations"
        params = {
            "sampleIds": sample_ids,
            "projection": "DETAILED",
            "pageSize": 10000
        }
        result = self._make_request("GET", endpoint, params=params)
        return result if isinstance(result, list) else []
    
    def get_clinical_data(self, study_id: str, entity_type: str = "PATIENT") -> Dict[str, Dict]:
        """
        Fetch clinical data for a study using GET endpoint (more reliable).
        
        Args:
            study_id: cBioPortal study ID
            entity_type: "PATIENT" or "SAMPLE"
        
        Returns:
            Dict mapping entity_id -> {attribute_id: value, ...}
        """
        endpoint = f"/studies/{study_id}/clinical-data"
        params = {
            "clinicalDataType": entity_type,
            "projection": "DETAILED"
        }
        
        # Use longer timeout for large clinical data responses
        original_timeout = self.timeout
        self.timeout = 300.0  # 5 minutes for large responses
        
        try:
            # Rate limiting before large request
            time.sleep(self.rate_limit_delay)
            
            url = f"{self.base_url}{endpoint}"
            client = self._get_client()
            
            print(f"   Making request to {endpoint}...")
            response = client.get(url, params=params)
            response.raise_for_status()
            
            print(f"   Response status: {response.status_code}, length: {len(response.text)} bytes")
            
            # Parse JSON
            if not response.text or not response.text.strip():
                raise Exception("Empty response from clinical data endpoint")
            
            result = response.json()
            print(f"   Parsed {len(result)} clinical data rows")
            
        finally:
            self.timeout = original_timeout
        
        # Convert list of clinical data rows to dict
        all_data = {}
        if isinstance(result, list):
            print(f"   Processing {len(result)} clinical data rows...")
            for row in result:
                entity_id = row.get("entityId")
                if entity_id:
                    if entity_id not in all_data:
                        all_data[entity_id] = {}
                    attr_id = row.get("clinicalAttributeId")
                    value = row.get("value")
                    if attr_id and value:
                        all_data[entity_id][attr_id] = value
        
        return all_data
    
    def get_clinical_attributes(self, study_id: str, clinical_data_type: str = "PATIENT") -> List[Dict]:
        """Get list of available clinical attributes for a study."""
        endpoint = f"/studies/{study_id}/clinical-data"
        params = {
            "clinicalDataType": clinical_data_type,
            "projection": "DETAILED"
        }
        return self._make_request("GET", endpoint, params=params)
