#!/usr/bin/env python3
"""
MOAT Orchestrator Validation Script

Tests the orchestrator endpoints with REAL patient data to validate:
1. API endpoint functionality
2. Data transformation requirements
3. Response format validation
4. Pipeline phase progression
5. Agent output completeness

Uses REAL TCGA data and test patient profiles.
"""

import asyncio
import httpx
import json
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime

API_BASE = "http://localhost:8000"

# Real test patient profiles
TEST_PATIENTS = {
    "ayesha_brca1": {
        "disease": "ovarian",
        "mutations": [
            {
                "gene": "BRCA1",
                "hgvs_p": "p.Arg1835Ter",
                "hgvs_c": "c.5503C>T",
                "chrom": "17",
                "pos": 43044295,
                "ref": "C",
                "alt": "T",
                "consequence": "stop_gained"
            },
            {
                "gene": "TP53",
                "hgvs_p": "p.Arg175His",
                "hgvs_c": "c.524G>A",
                "chrom": "17",
                "pos": 7673802,
                "ref": "G",
                "alt": "A",
                "consequence": "missense_variant"
            }
        ],
        "treatment_line": 1,
        "prior_therapies": []
    },
    "kras_tp53": {
        "disease": "ovarian",
        "mutations": [
            {
                "gene": "KRAS",
                "hgvs_p": "p.G12D",
                "hgvs_c": "c.35G>A",
                "chrom": "12",
                "pos": 25398284,
                "ref": "G",
                "alt": "A",
                "consequence": "missense_variant",
                "vaf": 0.45
            },
            {
                "gene": "TP53",
                "hgvs_p": "p.R175H",
                "hgvs_c": "c.524G>A",
                "chrom": "17",
                "pos": 7673802,
                "ref": "G",
                "alt": "A",
                "consequence": "missense_variant",
                "vaf": 0.32
            }
        ],
        "treatment_line": 1,
        "prior_therapies": []
    },
    "mbd4_tp53_nf1": {
        "disease": "ovarian",
        "mutations": [
            {
                "gene": "MBD4",
                "hgvs_p": "c.1239delA",
                "consequence": "frameshift",
                "zygosity": "homozygous"
            },
            {
                "gene": "TP53",
                "hgvs_p": "p.R273H",
                "consequence": "missense"
            },
            {
                "gene": "NF1",
                "hgvs_p": "p.R1276*",
                "consequence": "stop_gained"
            }
        ],
        "treatment_line": 1,
        "prior_therapies": []
    }
}


class ValidationResult:
    """Stores validation test results."""
    def __init__(self, test_name: str):
        self.test_name = test_name
        self.passed = False
        self.errors = []
        self.warnings = []
        self.data = {}
        self.duration_ms = 0
        
    def add_error(self, error: str):
        self.errors.append(error)
        self.passed = False
        
    def add_warning(self, warning: str):
        self.warnings.append(warning)
        
    def set_data(self, key: str, value: Any):
        self.data[key] = value
        
    def __str__(self):
        status = "✅ PASS" if self.passed and not self.errors else "❌ FAIL"
        result = f"{status} - {self.test_name}\n"
        if self.errors:
            result += f"  Errors: {len(self.errors)}\n"
            for err in self.errors[:3]:  # Show first 3
                result += f"    - {err}\n"
        if self.warnings:
            result += f"  Warnings: {len(self.warnings)}\n"
        if self.duration_ms:
            result += f"  Duration: {self.duration_ms:.0f}ms\n"
        return result


async def test_health_check() -> ValidationResult:
    """Test orchestrator health endpoint."""
    result = ValidationResult("Health Check")
    start = datetime.now()
    
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(f"{API_BASE}/api/orchestrate/health")
            
            if response.status_code == 200:
                data = response.json()
                result.set_data("status", data.get("status"))
                result.set_data("service", data.get("service"))
                result.passed = True
            else:
                result.add_error(f"Health check returned {response.status_code}")
                
    except Exception as e:
        result.add_error(f"Health check failed: {str(e)}")
        
    result.duration_ms = (datetime.now() - start).total_seconds() * 1000
    return result


async def test_orchestrator_full_pipeline(patient_key: str, patient_data: Dict) -> ValidationResult:
    """Test /api/orchestrate/full endpoint with real patient data."""
    result = ValidationResult(f"Full Pipeline - {patient_key}")
    start = datetime.now()
    
    try:
        # Build request according to OrchestratePipelineRequest schema
        request = {
            "disease": patient_data["disease"],
            "mutations": patient_data["mutations"],
            "treatment_line": patient_data.get("treatment_line", 1),
            "prior_therapies": patient_data.get("prior_therapies", []),
            "skip_agents": []
        }
        
        async with httpx.AsyncClient(timeout=180.0) as client:
            print(f"  📤 Sending request for {patient_key}...")
            response = await client.post(
                f"{API_BASE}/api/orchestrate/full",
                json=request
            )
            
            if response.status_code != 200:
                error_text = await response.aread()
                result.add_error(f"API returned {response.status_code}: {error_text.decode()[:200]}")
                return result
                
            data = response.json()
            result.set_data("response", data)
            
            # Validate response structure
            required_fields = [
                "patient_id", "disease", "phase", "progress_percent",
                "completed_agents", "created_at", "updated_at"
            ]
            
            for field in required_fields:
                if field not in data:
                    result.add_error(f"Missing required field: {field}")
                else:
                    result.set_data(field, data[field])
            
            # Validate phase
            valid_phases = ["initialized", "extracting", "analyzing", "ranking", 
                          "matching", "planning", "monitoring", "complete", "error"]
            phase = data.get("phase", "").lower()
            if phase not in valid_phases:
                result.add_warning(f"Unexpected phase: {phase}")
            else:
                result.set_data("phase", phase)
            
            # Validate agent outputs
            agent_outputs = {
                "biomarker_profile": data.get("biomarker_profile"),
                "resistance_prediction": data.get("resistance_prediction"),
                "drug_ranking": data.get("drug_ranking"),
                "trial_matches": data.get("trial_matches"),
                "nutrition_plan": data.get("nutrition_plan"),
                "synthetic_lethality_result": data.get("synthetic_lethality_result"),
                "care_plan": data.get("care_plan"),
                "monitoring_config": data.get("monitoring_config")
            }
            
            result.set_data("agent_outputs", agent_outputs)
            
            # Check which agents completed
            completed = data.get("completed_agents", [])
            result.set_data("completed_agents", completed)
            
            # Validate mechanism vector
            mechanism_vector = data.get("mechanism_vector")
            if mechanism_vector:
                if len(mechanism_vector) != 7:
                    result.add_error(f"Mechanism vector should be 7D, got {len(mechanism_vector)}")
                else:
                    result.set_data("mechanism_vector", mechanism_vector)
            else:
                result.add_warning("Mechanism vector not present")
            
            # Validate mutations
            mutation_count = data.get("mutation_count", 0)
            result.set_data("mutation_count", mutation_count)
            if mutation_count == 0:
                result.add_warning("No mutations extracted")
            
            # Check if pipeline completed
            if phase == "complete":
                result.passed = True
                # Validate all expected agents ran
                expected_agents = ["biomarker", "resistance", "drug_efficacy", 
                                 "trial_matching", "nutrition", "care_plan", "monitoring"]
                missing_agents = [a for a in expected_agents if a not in completed]
                if missing_agents:
                    result.add_warning(f"Expected agents not completed: {missing_agents}")
            elif phase == "error":
                result.add_error(f"Pipeline ended in error phase")
            else:
                result.add_warning(f"Pipeline not complete, phase: {phase}")
                
    except httpx.TimeoutException:
        result.add_error("Request timed out (>180s)")
    except Exception as e:
        result.add_error(f"Request failed: {str(e)}")
        import traceback
        result.set_data("traceback", traceback.format_exc())
        
    result.duration_ms = (datetime.now() - start).total_seconds() * 1000
    return result


async def test_status_endpoint(patient_id: str) -> ValidationResult:
    """Test /api/orchestrate/status/{patient_id} endpoint."""
    result = ValidationResult(f"Status Endpoint - {patient_id}")
    start = datetime.now()
    
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(f"{API_BASE}/api/orchestrate/status/{patient_id}")
            
            if response.status_code == 404:
                result.add_warning(f"Patient {patient_id} not found (may need to run pipeline first)")
                return result
                
            if response.status_code != 200:
                result.add_error(f"Status endpoint returned {response.status_code}")
                return result
                
            data = response.json()
            result.set_data("response", data)
            
            # Validate structure
            required_fields = ["patient_id", "phase", "progress_percent"]
            for field in required_fields:
                if field not in data:
                    result.add_error(f"Missing required field: {field}")
                else:
                    result.set_data(field, data[field])
            
            result.passed = True
            
    except Exception as e:
        result.add_error(f"Status check failed: {str(e)}")
        
    result.duration_ms = (datetime.now() - start).total_seconds() * 1000
    return result


async def test_state_endpoint(patient_id: str) -> ValidationResult:
    """Test /api/patients/{patient_id} endpoint (full state)."""
    result = ValidationResult(f"State Endpoint - {patient_id}")
    start = datetime.now()
    
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(f"{API_BASE}/api/patients/{patient_id}")
            
            if response.status_code == 404:
                result.add_warning(f"Patient {patient_id} not found")
                return result
                
            if response.status_code != 200:
                result.add_error(f"State endpoint returned {response.status_code}")
                return result
                
            data = response.json()
            result.set_data("response", data)
            
            # Validate all agent outputs are present
            agent_outputs = {
                "biomarker_profile": data.get("biomarker_profile"),
                "resistance_prediction": data.get("resistance_prediction"),
                "drug_ranking": data.get("drug_ranking"),
                "trial_matches": data.get("trial_matches"),
                "nutrition_plan": data.get("nutrition_plan"),
                "care_plan": data.get("care_plan"),
                "monitoring_config": data.get("monitoring_config")
            }
            
            result.set_data("agent_outputs", agent_outputs)
            
            # Check completeness
            missing_outputs = [k for k, v in agent_outputs.items() if v is None]
            if missing_outputs:
                result.add_warning(f"Missing agent outputs: {missing_outputs}")
            else:
                result.passed = True
                
    except Exception as e:
        result.add_error(f"State check failed: {str(e)}")
        
    result.duration_ms = (datetime.now() - start).total_seconds() * 1000
    return result


async def validate_data_transformation(orchestrator_response: Dict) -> ValidationResult:
    """Validate that orchestrator response can be transformed to legacy format."""
    result = ValidationResult("Data Transformation Validation")
    
    try:
        # Check if we can map orchestrator response to legacy format
        legacy_mapping = {
            "biomarker_intelligence": orchestrator_response.get("biomarker_profile"),
            "resistance_prediction": orchestrator_response.get("resistance_prediction"),
            "wiwfm": {
                "drugs": orchestrator_response.get("drug_ranking", []),
                "evidence_tier": orchestrator_response.get("drug_ranking", [{}])[0].get("tier") if orchestrator_response.get("drug_ranking") else None,
                "run_signature": orchestrator_response.get("patient_id")
            },
            "trials": {
                "trials": orchestrator_response.get("trial_matches", [])
            },
            "mechanism_map": {
                "vector": orchestrator_response.get("mechanism_vector", [0,0,0,0,0,0,0])
            }
        }
        
        result.set_data("legacy_mapping", legacy_mapping)
        
        # Validate mapping completeness
        if not legacy_mapping["biomarker_intelligence"]:
            result.add_warning("biomarker_intelligence mapping missing")
        if not legacy_mapping["wiwfm"]["drugs"]:
            result.add_warning("wiwfm.drugs mapping missing")
        if not legacy_mapping["trials"]["trials"]:
            result.add_warning("trials.trials mapping missing")
        if not legacy_mapping["mechanism_map"]["vector"]:
            result.add_warning("mechanism_map.vector mapping missing")
            
        # Check mechanism vector format
        mechanism_vector = legacy_mapping["mechanism_map"]["vector"]
        if mechanism_vector and len(mechanism_vector) != 7:
            result.add_error(f"Mechanism vector should be 7D, got {len(mechanism_vector)}")
        else:
            result.passed = True
            
    except Exception as e:
        result.add_error(f"Transformation validation failed: {str(e)}")
        
    return result


async def main():
    """Run all validation tests."""
    print("="*80)
    print("MOAT ORCHESTRATOR VALIDATION - REAL DATA TESTING")
    print("="*80)
    print(f"API Base: {API_BASE}")
    print(f"Timestamp: {datetime.now().isoformat()}\n")
    
    results = []
    
    # Test 1: Health Check
    print("🔍 Test 1: Health Check")
    health_result = await test_health_check()
    print(health_result)
    results.append(health_result)
    
    if not health_result.passed:
        print("\n❌ Health check failed - server may not be running")
        print("   Start server: cd oncology-coPilot/oncology-backend-minimal && uvicorn main:app --reload")
        return
    
    # Test 2-4: Full Pipeline with Real Patients
    patient_ids = []
    for patient_key, patient_data in TEST_PATIENTS.items():
        print(f"\n🔍 Test: Full Pipeline - {patient_key}")
        pipeline_result = await test_orchestrator_full_pipeline(patient_key, patient_data)
        print(pipeline_result)
        results.append(pipeline_result)
        
        if pipeline_result.passed and "patient_id" in pipeline_result.data:
            patient_id = pipeline_result.data["patient_id"]
            patient_ids.append(patient_id)
            
            # Test data transformation
            if "response" in pipeline_result.data:
                print(f"\n  🔍 Validating data transformation for {patient_id}...")
                transform_result = await validate_data_transformation(pipeline_result.data["response"])
                print(f"  {transform_result}")
                results.append(transform_result)
    
    # Test 5-6: Status and State Endpoints
    for patient_id in patient_ids[:2]:  # Test first 2
        print(f"\n🔍 Test: Status Endpoint - {patient_id}")
        status_result = await test_status_endpoint(patient_id)
        print(status_result)
        results.append(status_result)
        
        print(f"\n🔍 Test: State Endpoint - {patient_id}")
        state_result = await test_state_endpoint(patient_id)
        print(state_result)
        results.append(state_result)
    
    # Summary
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    passed = sum(1 for r in results if r.passed and not r.errors)
    failed = sum(1 for r in results if r.errors)
    warnings = sum(1 for r in results if r.warnings)
    
    print(f"Total Tests: {len(results)}")
    print(f"✅ Passed: {passed}")
    print(f"❌ Failed: {failed}")
    print(f"⚠️  Warnings: {warnings}")
    
    if failed > 0:
        print("\n❌ FAILED TESTS:")
        for r in results:
            if r.errors:
                print(f"  - {r.test_name}: {len(r.errors)} errors")
                for err in r.errors[:2]:
                    print(f"    • {err}")
    
    # Save detailed results
    output_file = Path(".cursor/MOAT/validation_results.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    results_dict = {
        "timestamp": datetime.now().isoformat(),
        "api_base": API_BASE,
        "summary": {
            "total": len(results),
            "passed": passed,
            "failed": failed,
            "warnings": warnings
        },
        "results": [
            {
                "test_name": r.test_name,
                "passed": r.passed,
                "errors": r.errors,
                "warnings": r.warnings,
                "duration_ms": r.duration_ms,
                "data_keys": list(r.data.keys())
            }
            for r in results
        ]
    }
    
    with open(output_file, "w") as f:
        json.dump(results_dict, f, indent=2)
    
    print(f"\n📄 Detailed results saved to: {output_file}")
    
    if failed == 0:
        print("\n🎉 All validation tests passed!")
    else:
        print(f"\n⚠️  {failed} test(s) failed. Review results above.")


if __name__ == "__main__":
    asyncio.run(main())

