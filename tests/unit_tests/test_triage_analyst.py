import sys
import os
import json

# Ensure the src directory is on the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

from tools.candidate_triage import CandidateTriage

def run_analyst_test():
    """
    Runs a focused test on the CandidateTriage tool to debug the LLM response.
    """
    print("--- Firing Range: Testing the Triage Analyst ---")

    # This candidate caused an error in the last run.
    problematic_candidate = "YRYYRMMRMMAYMRMRRARRYARMRRYRRY"

    # Define a dummy threat matrix for this test
    threat_matrix = {
        'high_value_keywords': ['binding', 'active site'],
        'low_value_keywords': ['repetitive', 'DNA']
    }

    triage_tool = CandidateTriage(threat_matrix=threat_matrix)

    print(f"\nSubmitting candidate for assessment:\n{problematic_candidate}")
    
    # This will now be run with enhanced logging inside the function
    assessment = triage_tool.get_triage_assessment(problematic_candidate)

    print("\n--- Assessment Result ---")
    print(json.dumps(assessment, indent=2))
    print("\n--- Firing Range Test Complete ---")


if __name__ == "__main__":
    run_analyst_test() 
import os
import json

# Ensure the src directory is on the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

from tools.candidate_triage import CandidateTriage

def run_analyst_test():
    """
    Runs a focused test on the CandidateTriage tool to debug the LLM response.
    """
    print("--- Firing Range: Testing the Triage Analyst ---")

    # This candidate caused an error in the last run.
    problematic_candidate = "YRYYRMMRMMAYMRMRRARRYARMRRYRRY"

    # Define a dummy threat matrix for this test
    threat_matrix = {
        'high_value_keywords': ['binding', 'active site'],
        'low_value_keywords': ['repetitive', 'DNA']
    }

    triage_tool = CandidateTriage(threat_matrix=threat_matrix)

    print(f"\nSubmitting candidate for assessment:\n{problematic_candidate}")
    
    # This will now be run with enhanced logging inside the function
    assessment = triage_tool.get_triage_assessment(problematic_candidate)

    print("\n--- Assessment Result ---")
    print(json.dumps(assessment, indent=2))
    print("\n--- Firing Range Test Complete ---")


if __name__ == "__main__":
    run_analyst_test() 
import os
import json

# Ensure the src directory is on the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

from tools.candidate_triage import CandidateTriage

def run_analyst_test():
    """
    Runs a focused test on the CandidateTriage tool to debug the LLM response.
    """
    print("--- Firing Range: Testing the Triage Analyst ---")

    # This candidate caused an error in the last run.
    problematic_candidate = "YRYYRMMRMMAYMRMRRARRYARMRRYRRY"

    # Define a dummy threat matrix for this test
    threat_matrix = {
        'high_value_keywords': ['binding', 'active site'],
        'low_value_keywords': ['repetitive', 'DNA']
    }

    triage_tool = CandidateTriage(threat_matrix=threat_matrix)

    print(f"\nSubmitting candidate for assessment:\n{problematic_candidate}")
    
    # This will now be run with enhanced logging inside the function
    assessment = triage_tool.get_triage_assessment(problematic_candidate)

    print("\n--- Assessment Result ---")
    print(json.dumps(assessment, indent=2))
    print("\n--- Firing Range Test Complete ---")


if __name__ == "__main__":
    run_analyst_test() 
import os
import json

# Ensure the src directory is on the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))

from tools.candidate_triage import CandidateTriage

def run_analyst_test():
    """
    Runs a focused test on the CandidateTriage tool to debug the LLM response.
    """
    print("--- Firing Range: Testing the Triage Analyst ---")

    # This candidate caused an error in the last run.
    problematic_candidate = "YRYYRMMRMMAYMRMRRARRYARMRRYRRY"

    # Define a dummy threat matrix for this test
    threat_matrix = {
        'high_value_keywords': ['binding', 'active site'],
        'low_value_keywords': ['repetitive', 'DNA']
    }

    triage_tool = CandidateTriage(threat_matrix=threat_matrix)

    print(f"\nSubmitting candidate for assessment:\n{problematic_candidate}")
    
    # This will now be run with enhanced logging inside the function
    assessment = triage_tool.get_triage_assessment(problematic_candidate)

    print("\n--- Assessment Result ---")
    print(json.dumps(assessment, indent=2))
    print("\n--- Firing Range Test Complete ---")


if __name__ == "__main__":
    run_analyst_test() 