#!/usr/bin/env python3
"""
RUNX1 Integration Testing Suite
Comprehensive testing for RUNX1-FPD platform integration
"""

import sys
import os
import logging
import traceback
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# Add tools directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# Import RUNX1 components
try:
    from tools.runx1_integration_plan import (
        integrate_runx1_analysis,
        design_runx1_intervention,
        create_runx1_genomic_browser,
        generate_runx1_clinical_report,
        RUNX1DigitalTwinIntegrator
    )
    from tools.runx1_genomic_browser import RUNX1GenomicBrowser
    from tools.runx1_demo_scenarios import RUNX1DemoScenarios
    from tools.runx1_data_loader import RUNX1DataLoader
    from tools.runx1_progression_modeler import RUNX1ProgressionModeler
    from tools.chopchop_integration import ChopChopIntegration
except ImportError as e:
    print(f"Import error: {e}")
    print("Some components may not be available for testing")

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RUNX1IntegrationTester:
    """
    Comprehensive testing suite for RUNX1-FPD platform integration.
    """
    
    def __init__(self):
        """Initialize the testing suite."""
        self.test_results = []
        self.start_time = datetime.now()
        
        # Test data
        self.test_germline_variant = {
            "gene": "RUNX1",
            "protein_change": "p.Arg135fs",
            "variant_type": "germline",
            "id": "RUNX1_p.Arg135fs",
            "position": 36207748,
            "consequence": "frameshift_variant",
            "pathogenicity": "pathogenic",
            "functional_impact": "high"
        }
        
        self.test_somatic_variant = {
            "gene": "ASXL1",
            "protein_change": "p.Gly646fs",
            "variant_type": "somatic",
            "id": "ASXL1_p.Gly646fs",
            "position": 31022441,
            "consequence": "frameshift_variant",
            "pathogenicity": "pathogenic",
            "functional_impact": "high"
        }
    
    def run_test(self, test_name: str, test_function, *args, **kwargs) -> bool:
        """
        Run a single test and record results.
        
        Args:
            test_name: Name of the test
            test_function: Function to test
            *args: Arguments for test function
            **kwargs: Keyword arguments for test function
            
        Returns:
            bool: True if test passed, False otherwise
        """
        logger.info(f"Running test: {test_name}")
        
        try:
            start_time = datetime.now()
            result = test_function(*args, **kwargs)
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            # Determine if test passed
            passed = result is not None and result != {}
            
            test_result = {
                "name": test_name,
                "passed": passed,
                "duration": duration,
                "result": result,
                "error": None,
                "timestamp": start_time.isoformat()
            }
            
            self.test_results.append(test_result)
            
            if passed:
                logger.info(f"✅ {test_name} PASSED ({duration:.2f}s)")
            else:
                logger.warning(f"⚠️ {test_name} FAILED - No result returned ({duration:.2f}s)")
            
            return passed
            
        except Exception as e:
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            test_result = {
                "name": test_name,
                "passed": False,
                "duration": duration,
                "result": None,
                "error": str(e),
                "traceback": traceback.format_exc(),
                "timestamp": start_time.isoformat()
            }
            
            self.test_results.append(test_result)
            logger.error(f"❌ {test_name} FAILED ({duration:.2f}s): {e}")
            
            return False
    
    def test_data_loader(self) -> bool:
        """Test RUNX1 data loader functionality."""
        def test_func():
            try:
                loader = RUNX1DataLoader()
                
                # Test sequence loading
                sequence = loader.load_runx1_sequence()
                if not sequence:
                    return None
                
                # Test variant loading
                variants = loader.load_runx1_variants()
                if not variants:
                    return None
                
                # Test genomic context
                context = loader.get_genomic_context(36207748, 1000)
                if not context:
                    return None
                
                # Test functional domain mapping
                domains = loader.map_functional_domains(36207748)
                if not domains:
                    return None
                
                return {
                    "sequence_length": len(sequence) if sequence else 0,
                    "variants_count": len(variants) if variants else 0,
                    "context_available": bool(context),
                    "domains_mapped": bool(domains)
                }
            except Exception as e:
                logger.error(f"Data loader test failed: {e}")
                return None
        
        return self.run_test("RUNX1 Data Loader", test_func)
    
    def test_progression_modeler(self) -> bool:
        """Test RUNX1 progression modeler functionality."""
        def test_func():
            try:
                modeler = RUNX1ProgressionModeler()
                
                # Test progression modeling
                result = modeler.model_progression(
                    self.test_germline_variant,
                    self.test_somatic_variant,
                    patient_age=35.0
                )
                
                if not result:
                    return None
                
                # Validate result structure
                required_keys = ["transformation_probability", "risk_factors", "intervention_opportunities"]
                for key in required_keys:
                    if key not in result:
                        logger.warning(f"Missing key in progression model: {key}")
                        return None
                
                return {
                    "transformation_probability": result.get("transformation_probability", {}),
                    "risk_factors_count": len(result.get("risk_factors", [])),
                    "intervention_opportunities_count": len(result.get("intervention_opportunities", [])),
                    "model_complete": True
                }
            except Exception as e:
                logger.error(f"Progression modeler test failed: {e}")
                return None
        
        return self.run_test("RUNX1 Progression Modeler", test_func)
    
    def test_integration_analysis(self) -> bool:
        """Test integrated RUNX1 analysis functionality."""
        def test_func():
            try:
                # Test integrated analysis
                result = integrate_runx1_analysis(
                    self.test_germline_variant,
                    self.test_somatic_variant,
                    patient_age=35.0
                )
                
                if not result:
                    return None
                
                # Validate result structure
                required_sections = ["germline_analysis", "somatic_analysis", "progression_model"]
                for section in required_sections:
                    if section not in result:
                        logger.warning(f"Missing section in analysis: {section}")
                        return None
                
                return {
                    "germline_analysis_complete": bool(result.get("germline_analysis")),
                    "somatic_analysis_complete": bool(result.get("somatic_analysis")),
                    "progression_model_complete": bool(result.get("progression_model")),
                    "risk_stratification_complete": bool(result.get("risk_stratification")),
                    "clinical_recommendations_complete": bool(result.get("clinical_recommendations"))
                }
            except Exception as e:
                logger.error(f"Integration analysis test failed: {e}")
                return None
        
        return self.run_test("RUNX1 Integration Analysis", test_func)
    
    def test_intervention_design(self) -> bool:
        """Test RUNX1 intervention design functionality."""
        def test_func():
            try:
                # Test intervention design
                result = design_runx1_intervention(
                    self.test_somatic_variant,
                    intervention_type="crispr_correction"
                )
                
                if not result:
                    return None
                
                # Validate result structure
                required_keys = ["guide_design", "safety_assessment", "experimental_protocol"]
                for key in required_keys:
                    if key not in result:
                        logger.warning(f"Missing key in intervention design: {key}")
                        return None
                
                guide_design = result.get("guide_design", {})
                guides_count = len(guide_design.get("all_guides", []))
                
                return {
                    "guides_found": guides_count,
                    "success_probability": result.get("success_probability", 0.0),
                    "timeline_available": bool(result.get("timeline_estimate")),
                    "resources_calculated": bool(result.get("resource_requirements")),
                    "intervention_complete": True
                }
            except Exception as e:
                logger.error(f"Intervention design test failed: {e}")
                return None
        
        return self.run_test("RUNX1 Intervention Design", test_func)
    
    def test_genomic_browser(self) -> bool:
        """Test RUNX1 genomic browser functionality."""
        def test_func():
            try:
                browser = RUNX1GenomicBrowser()
                
                # Test genomic overview creation
                fig = browser.create_genomic_overview(36207748, 50000)
                if not fig:
                    return None
                
                # Test variant impact visualization
                fig_impact = browser.create_variant_impact_visualization(self.test_germline_variant)
                if not fig_impact:
                    return None
                
                # Test comparison view
                fig_comparison = browser.create_comparison_view([self.test_germline_variant])
                if not fig_comparison:
                    return None
                
                return {
                    "genomic_overview_created": bool(fig),
                    "variant_impact_created": bool(fig_impact),
                    "comparison_view_created": bool(fig_comparison),
                    "browser_functional": True
                }
            except Exception as e:
                logger.error(f"Genomic browser test failed: {e}")
                return None
        
        return self.run_test("RUNX1 Genomic Browser", test_func)
    
    def test_demo_scenarios(self) -> bool:
        """Test RUNX1 demo scenarios functionality."""
        def test_func():
            try:
                scenarios = RUNX1DemoScenarios()
                
                # Test scenario retrieval
                all_scenarios = scenarios.get_all_scenarios()
                if not all_scenarios:
                    return None
                
                # Test specific scenario
                scenario = scenarios.get_scenario("scenario_1_high_risk")
                if not scenario:
                    return None
                
                # Test demo flow creation
                demo_flow = scenarios.create_demo_flow("scenario_1_high_risk")
                if not demo_flow:
                    return None
                
                return {
                    "scenarios_count": len(all_scenarios),
                    "scenario_structure_valid": bool(scenario.get("patient_info")),
                    "demo_flow_created": bool(demo_flow.get("demo_script")),
                    "scenarios_functional": True
                }
            except Exception as e:
                logger.error(f"Demo scenarios test failed: {e}")
                return None
        
        return self.run_test("RUNX1 Demo Scenarios", test_func)
    
    def test_clinical_report_generation(self) -> bool:
        """Test clinical report generation functionality."""
        def test_func():
            try:
                # First run analysis to get results
                analysis_result = integrate_runx1_analysis(
                    self.test_germline_variant,
                    self.test_somatic_variant,
                    patient_age=35.0
                )
                
                if not analysis_result:
                    return None
                
                # Test clinical report generation
                clinical_report = generate_runx1_clinical_report(
                    analysis_result,
                    {"age": 35, "id": "TEST_001"}
                )
                
                if not clinical_report:
                    return None
                
                # Validate report structure
                required_sections = ["patient_summary", "genetic_profile", "risk_assessment"]
                for section in required_sections:
                    if section not in clinical_report:
                        logger.warning(f"Missing section in clinical report: {section}")
                        return None
                
                return {
                    "patient_summary_complete": bool(clinical_report.get("patient_summary")),
                    "genetic_profile_complete": bool(clinical_report.get("genetic_profile")),
                    "risk_assessment_complete": bool(clinical_report.get("risk_assessment")),
                    "clinical_recommendations_complete": bool(clinical_report.get("clinical_recommendations")),
                    "report_functional": True
                }
            except Exception as e:
                logger.error(f"Clinical report generation test failed: {e}")
                return None
        
        return self.run_test("Clinical Report Generation", test_func)
    
    def test_end_to_end_workflow(self) -> bool:
        """Test complete end-to-end workflow."""
        def test_func():
            try:
                # Initialize integrator
                integrator = RUNX1DigitalTwinIntegrator()
                
                # Step 1: Run enhanced two-hit analysis
                analysis_result = integrator.enhanced_two_hit_analysis(
                    self.test_germline_variant,
                    self.test_somatic_variant,
                    patient_age=35.0
                )
                
                if not analysis_result:
                    return None
                
                # Step 2: Design intervention
                intervention_result = integrator.design_precision_intervention(
                    self.test_somatic_variant,
                    intervention_type="crispr_correction"
                )
                
                if not intervention_result:
                    return None
                
                # Step 3: Create genomic browser data
                browser_data = integrator.create_genomic_browser_data(
                    self.test_germline_variant.get("position", 36207748),
                    window=5000
                )
                
                if not browser_data:
                    return None
                
                # Step 4: Generate clinical report
                clinical_report = integrator.generate_clinical_report(
                    analysis_result,
                    {"age": 35, "id": "E2E_TEST"}
                )
                
                if not clinical_report:
                    return None
                
                return {
                    "analysis_complete": bool(analysis_result),
                    "intervention_designed": bool(intervention_result),
                    "genomic_data_created": bool(browser_data),
                    "clinical_report_generated": bool(clinical_report),
                    "end_to_end_successful": True,
                    "workflow_time": (datetime.now() - self.start_time).total_seconds()
                }
            except Exception as e:
                logger.error(f"End-to-end workflow test failed: {e}")
                return None
        
        return self.run_test("End-to-End Workflow", test_func)
    
    def run_all_tests(self) -> Dict:
        """Run all tests and return summary."""
        logger.info("🚀 Starting RUNX1 Integration Test Suite")
        logger.info("=" * 60)
        
        # Run individual component tests
        test_methods = [
            self.test_data_loader,
            self.test_progression_modeler,
            self.test_integration_analysis,
            self.test_intervention_design,
            self.test_genomic_browser,
            self.test_demo_scenarios,
            self.test_clinical_report_generation,
            self.test_end_to_end_workflow
        ]
        
        passed_tests = 0
        total_tests = len(test_methods)
        
        for test_method in test_methods:
            if test_method():
                passed_tests += 1
        
        # Calculate summary
        end_time = datetime.now()
        total_duration = (end_time - self.start_time).total_seconds()
        
        summary = {
            "total_tests": total_tests,
            "passed_tests": passed_tests,
            "failed_tests": total_tests - passed_tests,
            "success_rate": passed_tests / total_tests if total_tests > 0 else 0,
            "total_duration": total_duration,
            "start_time": self.start_time.isoformat(),
            "end_time": end_time.isoformat(),
            "test_results": self.test_results
        }
        
        # Log summary
        logger.info("=" * 60)
        logger.info("🏁 RUNX1 Integration Test Suite Complete")
        logger.info(f"✅ Passed: {passed_tests}/{total_tests} ({summary['success_rate']:.1%})")
        logger.info(f"⏱️ Duration: {total_duration:.2f} seconds")
        
        if passed_tests == total_tests:
            logger.info("🎉 ALL TESTS PASSED! RUNX1 integration is ready for demo!")
        else:
            logger.warning(f"⚠️ {total_tests - passed_tests} tests failed. Review results before demo.")
        
        return summary
    
    def generate_test_report(self) -> str:
        """Generate detailed test report."""
        if not self.test_results:
            return "No test results available."
        
        report = []
        report.append("# RUNX1 Integration Test Report")
        report.append(f"Generated: {datetime.now().isoformat()}")
        report.append("")
        
        # Summary
        total_tests = len(self.test_results)
        passed_tests = sum(1 for result in self.test_results if result["passed"])
        failed_tests = total_tests - passed_tests
        
        report.append("## Summary")
        report.append(f"- Total Tests: {total_tests}")
        report.append(f"- Passed: {passed_tests}")
        report.append(f"- Failed: {failed_tests}")
        report.append(f"- Success Rate: {passed_tests/total_tests:.1%}")
        report.append("")
        
        # Detailed results
        report.append("## Detailed Results")
        report.append("")
        
        for result in self.test_results:
            status = "✅ PASS" if result["passed"] else "❌ FAIL"
            report.append(f"### {result['name']} - {status}")
            report.append(f"- Duration: {result['duration']:.2f}s")
            report.append(f"- Timestamp: {result['timestamp']}")
            
            if result["error"]:
                report.append(f"- Error: {result['error']}")
            
            if result["result"]:
                report.append(f"- Result: {result['result']}")
            
            report.append("")
        
        return "\n".join(report)

def main():
    """Main testing function."""
    print("🚀 RUNX1-FPD Integration Testing Suite")
    print("=" * 50)
    
    # Initialize tester
    tester = RUNX1IntegrationTester()
    
    # Run all tests
    summary = tester.run_all_tests()
    
    # Generate report
    report = tester.generate_test_report()
    
    # Save report
    report_file = f"runx1_integration_test_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
    with open(report_file, 'w') as f:
        f.write(report)
    
    print(f"\n📄 Test report saved to: {report_file}")
    
    # Return summary for programmatic use
    return summary

if __name__ == "__main__":
    main() 