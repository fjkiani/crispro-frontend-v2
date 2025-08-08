# RUNX1 Integration Test Report
Generated: 2025-07-11T11:53:42.608554

## Summary
- Total Tests: 8
- Passed: 4
- Failed: 4
- Success Rate: 50.0%

## Detailed Results

### RUNX1 Data Loader - ❌ FAIL
- Duration: 0.00s
- Timestamp: 2025-07-11T11:47:44.344667

### RUNX1 Progression Modeler - ❌ FAIL
- Duration: 0.00s
- Timestamp: 2025-07-11T11:47:44.344710

### RUNX1 Integration Analysis - ✅ PASS
- Duration: 0.05s
- Timestamp: 2025-07-11T11:47:44.344755
- Result: {'germline_analysis_complete': True, 'somatic_analysis_complete': True, 'progression_model_complete': True, 'risk_stratification_complete': True, 'clinical_recommendations_complete': True}

### RUNX1 Intervention Design - ✅ PASS
- Duration: 185.01s
- Timestamp: 2025-07-11T11:47:44.391398
- Result: {'guides_found': 5, 'success_probability': 0.1, 'timeline_available': True, 'resources_calculated': True, 'intervention_complete': True}

### RUNX1 Genomic Browser - ❌ FAIL
- Duration: 0.00s
- Timestamp: 2025-07-11T11:50:49.400582

### RUNX1 Demo Scenarios - ❌ FAIL
- Duration: 0.00s
- Timestamp: 2025-07-11T11:50:49.400612

### Clinical Report Generation - ✅ PASS
- Duration: 0.01s
- Timestamp: 2025-07-11T11:50:49.400637
- Result: {'patient_summary_complete': True, 'genetic_profile_complete': True, 'risk_assessment_complete': True, 'clinical_recommendations_complete': True, 'report_functional': True}

### End-to-End Workflow - ✅ PASS
- Duration: 173.20s
- Timestamp: 2025-07-11T11:50:49.408883
- Result: {'analysis_complete': True, 'intervention_designed': True, 'genomic_data_created': True, 'clinical_report_generated': True, 'end_to_end_successful': True, 'workflow_time': 358.263784}
