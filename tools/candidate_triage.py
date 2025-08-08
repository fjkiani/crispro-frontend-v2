
import re
from typing import List, Dict, Any
import json

from tools.llm_api import query_llm

class CandidateTriage:
    """
    A tool to triage and filter candidate protein sequences.
    """

    def __init__(self, threat_matrix: Dict[str, List[str]]):
        """
        Initializes the triage tool with a threat matrix.

        Args:
            threat_matrix: A dictionary containing lists of keywords for triage.
                           Expected keys: 'high_value_keywords', 'low_value_keywords'.
        """
        self.threat_matrix = threat_matrix
        self.llm_provider = "anthropic" 

    def filter_low_complexity(self, candidates: List[str], max_repeat_pct: float = 0.8) -> List[str]:
        """
        Filters out low-complexity sequences based on character repetition.

        Args:
            candidates: A list of candidate sequences.
            max_repeat_pct: The maximum allowed percentage of the most frequent character.

        Returns:
            A list of candidates that pass the complexity filter.
        """
        high_quality_candidates = []
        for cand in candidates:
            if not cand:
                continue
            
            most_common_char_count = cand.count(max(set(cand), key=cand.count))
            if (most_common_char_count / len(cand)) < max_repeat_pct:
                high_quality_candidates.append(cand)
        return high_quality_candidates

    def get_triage_assessment(self, candidate: str) -> Dict[str, Any]:
        """
        Uses an LLM to assess a single candidate sequence based on the threat matrix.

        Args:
            candidate: The sequence to assess.

        Returns:
            A dictionary containing the assessment details (e.g., 'rating', 'reasoning').
        """
        prompt = self._build_triage_prompt(candidate)
        
        try:
            response_text = query_llm(
                prompt=prompt,
                provider=self.llm_provider,
                # In a real scenario, you might add more parameters like temperature
            )
            # The LLM is instructed to return a JSON object, so we parse it.
            # A more robust implementation would have retries and validation.
            assessment_data = json.loads(response_text)
            assessment_data['candidate'] = candidate # Ensure candidate is in the final dict
            return assessment_data

        except Exception as e:
            print(f"ERROR: LLM assessment failed for candidate: {candidate[:30]}... | Error: {e}")
            return {
                "candidate": candidate,
                "rating": "ERROR",
                "reasoning": f"Failed to get a valid assessment from the LLM. Error: {e}",
            }

    def _build_triage_prompt(self, candidate: str) -> str:
        """Builds the LLM prompt for triaging a candidate."""
        
        prompt = f"""
You are a master protein engineer and nanobody design specialist. Your mission is to triage a candidate amino acid sequence and determine if it is a viable candidate for a therapeutic inhibitor.

**Candidate Sequence:**
`{candidate}`

**Triage Keywords (High-Value):**
{self.threat_matrix.get('high_value_keywords', [])}

**Triage Keywords (Low-Value/Red Flags):**
{self.threat_matrix.get('low_value_keywords', [])}

**Analysis Instructions:**

1.  **Analyze the sequence composition.** Does it appear to be a plausible protein sequence? Does it contain unusual patterns or non-amino-acid characters?
2.  **Assess for red flags.** Does the sequence exhibit characteristics of a junk sequence (e.g., extreme repetitiveness, looks like DNA)?
3.  **Assess for high-value potential.** Does the sequence composition suggest it might form a structure capable of binding (e.g., variable regions, some structure)?
4.  **Provide a final rating.** Based on your analysis, rate the candidate as **"VIABLE"** or **"REJECTED"**.

**Output Format (JSON):**
Provide your response as a single JSON object with the following keys:
- `rating`: (string) "VIABLE" or "REJECTED".
- `reasoning`: (string) A brief, 1-2 sentence explanation for your rating.

**Example Reasoning:**
- "VIABLE: The sequence has a mix of hydrophobic and hydrophilic residues, suggesting it could fold into a stable structure. No obvious red flags."
- "REJECTED: The sequence is highly repetitive and contains non-amino acid characters (G, C, T), suggesting it is a DNA artifact, not a protein."

Begin your assessment now.
"""
        return prompt

    def run_triage(self, candidates: List[str]) -> List[Dict[str, Any]]:
        """
        Runs the full triage protocol on a list of candidates.

        Args:
            candidates: A list of candidate sequences.

        Returns:
            A list of triage assessment dictionaries for viable candidates.
        """
        print("--- Starting Candidate Triage Protocol ---")
        
        # Stage 1: Low-Complexity Purge
        print(f"INFO: Received {len(candidates)} candidates for triage.")
        complex_candidates = self.filter_low_complexity(candidates)
        print(f"INFO: {len(complex_candidates)} candidates passed low-complexity filter.")

        # Stage 2: AI Analyst Interrogation
        assessments = []
        for cand in complex_candidates:
            assessment = self.get_triage_assessment(cand)
            assessments.append(assessment)

        print("--- Triage Protocol Complete ---")
        return assessments

if __name__ == '__main__':
    # Example Usage
    dummy_candidates = [
        "AMRAYRAMYARRMAMMMYMMRRYMYRMRYMRRRYRAMAARYAYRMMRYMMARMARMARMYAAYRYYMAAAAGCGTGTTCGGTACAAAAACACGCTTTTTCCCTATTACCATCTTCACTTT",
        "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR",
        "AYYYAARARRRYYYRRRYRRYMYAYRYYRYRAYYRYAGYYRARRAYYMYAARAYRRAAYYRAAAARRAARYRWRYYAAYARYRRRWRYAYARARMRAWARRAARARYARYAARYYAAAAA",
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "YRMYYMRYAMARAARAARRRAAGGAAGCACCCCCCCCACCTGTTGCCTTCAGTATAAACTAAGCTCCAATACCTGACAATGTTTTAGAGAGTTGGGTCTTTTTTTTGCCTACAGAAAGCT",
        "ATGCATGCATGCATGCATGCATGC" # DNA-like
    ]

    # These would normally come from the database/NCBIClient
    dummy_threat_matrix = {
        'high_value_keywords': ['binding', 'active site', 'interface', 'loop', 'variable region'],
        'low_value_keywords': ['repetitive', 'non-standard characters', 'DNA', 'artifact']
    }

    triage_tool = CandidateTriage(threat_matrix=dummy_threat_matrix)
    
    final_report = triage_tool.run_triage(dummy_candidates)

    print("\n--- FINAL TRIAGE REPORT ---")
    import json
    print(json.dumps(final_report, indent=2)) 