#!/usr/bin/env python3
"""
Answer 8 Clinical Questions from MBD4+TP53 Analysis Results

Extracts structured answers to all 8 clinical questions from the complete
analysis results JSON file.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def answer_question_1_variant_impact(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 1: Variant Impact Prediction
    Which mutations are probable drivers?
    """
    mutations = results.get("mutations", [])
    efficacy = results.get("efficacy_prediction", {})
    drugs = efficacy.get("drugs", [])
    
    # Extract variant annotations from efficacy response
    variant_impacts = []
    for mut in mutations:
        gene = mut.get("gene")
        hgvs_p = mut.get("hgvs_p")
        
        # Find drug rankings that mention this gene
        # Handle rationale as list or string
    gene_drugs = []
    for d in drugs:
        rationale = d.get("rationale", "")
        if isinstance(rationale, list):
            rationale_str = " ".join([str(r) for r in rationale])
        else:
            rationale_str = str(rationale)
        if gene.lower() in rationale_str.lower():
            gene_drugs.append(d)
        
        # Get insights for this variant
        insights = results.get("insights_bundle", {})
        functionality = insights.get("functionality", {}).get("functionality_score", 0.0)
        essentiality = insights.get("essentiality", {}).get("essentiality_score", 0.0)
        
        variant_impacts.append({
            "gene": gene,
            "hgvs_p": hgvs_p,
            "consequence": mut.get("consequence"),
            "driver_probability": "HIGH" if functionality > 0.6 or essentiality > 0.7 else "MODERATE",
            "functionality_score": functionality,
            "essentiality_score": essentiality,
            "affects_drugs": len(gene_drugs),
            "proxy_sae_source": "Pathway scores indicate driver pathways (DDR for MBD4, TP53 pathway for TP53)"
        })
    
    return {
        "question": "Variant Impact Prediction",
        "answer": variant_impacts,
        "summary": f"{len([v for v in variant_impacts if v['driver_probability'] == 'HIGH'])} high-probability drivers identified"
    }


def answer_question_2_functional_annotation(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 2: Functional Annotation
    Protein-level effects?
    """
    insights = results.get("insights_bundle", {})
    mutations = results.get("mutations", [])
    
    functional_annotations = []
    for mut in mutations:
        gene = mut.get("gene")
        
        # Extract all 4 insight chips
        functionality = insights.get("functionality", {}).get("functionality_score", 0.0)
        chromatin = insights.get("chromatin", {}).get("chromatin_score", 0.0)
        essentiality = insights.get("essentiality", {}).get("essentiality_score", 0.0)
        regulatory_key = f"regulatory_{gene}"
        regulatory = insights.get(regulatory_key, {}).get("regulatory_impact_score", 0.0)
        
        functional_annotations.append({
            "gene": gene,
            "hgvs_p": mut.get("hgvs_p"),
            "protein_effect": {
                "functionality_change": functionality,
                "chromatin_accessibility": chromatin,
                "gene_essentiality": essentiality,
                "regulatory_impact": regulatory
            },
            "interpretation": {
                "functionality": "HIGH" if functionality > 0.6 else "MODERATE" if functionality > 0.4 else "LOW",
                "essentiality": "HIGH" if essentiality > 0.7 else "MODERATE" if essentiality > 0.5 else "LOW",
                "regulatory": "HIGH" if regulatory > 0.6 else "MODERATE" if regulatory > 0.4 else "LOW"
            },
            "proxy_sae_source": "Insights bundle feeds into SAE features"
        })
    
    return {
        "question": "Functional Annotation",
        "answer": functional_annotations,
        "summary": "Protein-level effects quantified via 4 insight chips"
    }


def answer_question_3_pathway_analysis(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 3: Pathway Analysis
    Dominant pathways and vulnerabilities?
    """
    pathway_scores = results.get("pathway_scores", {})
    sae_features = results.get("sae_features", {})
    
    # Sort pathways by score (descending)
    sorted_pathways = sorted(
        pathway_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )
    
    # Get DNA repair capacity from SAE features
    dna_repair = sae_features.get("dna_repair_capacity", 0.0) if isinstance(sae_features, dict) else 0.0
    
    dominant_pathways = []
    for pathway, score in sorted_pathways[:5]:  # Top 5
        dominant_pathways.append({
            "pathway": pathway.upper(),
            "disruption_score": score,
            "interpretation": "HIGH" if score > 0.7 else "MODERATE" if score > 0.4 else "LOW",
            "vulnerability": "Targetable" if score > 0.6 else "Moderate" if score > 0.4 else "Low"
        })
    
    return {
        "question": "Pathway Analysis",
        "answer": {
            "dominant_pathways": dominant_pathways,
            "dna_repair_capacity": dna_repair,
            "total_pathways": len(pathway_scores),
            "proxy_sae_source": "Pathway scores are the proxy SAE mechanism vector source"
        },
        "summary": f"Top pathway: {dominant_pathways[0]['pathway']} (score: {dominant_pathways[0]['disruption_score']:.2f})"
    }


def answer_question_4_drug_prediction(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 4: Drug and Therapy Prediction
    Most effective drugs?
    """
    efficacy = results.get("efficacy_prediction", {})
    drugs = efficacy.get("drugs", [])
    
    # Sort by efficacy score (descending)
    sorted_drugs = sorted(drugs, key=lambda x: x.get("efficacy_score", 0.0), reverse=True)
    
    top_drugs = []
    for drug in sorted_drugs[:10]:  # Top 10
        top_drugs.append({
            "name": drug.get("name", "Unknown"),
            "efficacy_score": drug.get("efficacy_score", 0.0),
            "confidence": drug.get("confidence", 0.0),
            "evidence_tier": drug.get("evidence_tier", "insufficient"),
            "badges": drug.get("badges", []),
            "rationale": drug.get("rationale", [])[:3],  # First 3 rationale points
            "proxy_sae_source": "Mechanism vector used for drug-pathway alignment"
        })
    
    return {
        "question": "Drug and Therapy Prediction",
        "answer": top_drugs,
        "summary": f"Top drug: {top_drugs[0]['name']} (efficacy: {top_drugs[0]['efficacy_score']:.2f}, confidence: {top_drugs[0]['confidence']:.2f})"
    }


def answer_question_5_trial_matching(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 5: Trial and Biomarker Matching
    Molecular fit trials?
    """
    trials_response = results.get("trial_matching", {})
    trials = trials_response.get("trials", []) if isinstance(trials_response, dict) else []
    
    # Sort by match score if available
    if trials and isinstance(trials[0], dict) and "match_score" in trials[0]:
        sorted_trials = sorted(trials, key=lambda x: x.get("match_score", 0.0), reverse=True)
    else:
        sorted_trials = trials[:10]  # Just take first 10
    
    matched_trials = []
    for trial in sorted_trials[:10]:
        matched_trials.append({
            "nct_id": trial.get("nct_id", "Unknown"),
            "title": trial.get("title", "Unknown"),
            "match_score": trial.get("match_score", 0.0),
            "mechanism_fit": trial.get("mechanism_fit", 0.0),
            "eligibility": trial.get("eligibility_score", 0.0),
            "proxy_sae_source": "Mechanism vector used for trial MoA matching"
        })
    
    return {
        "question": "Trial and Biomarker Matching",
        "answer": matched_trials,
        "summary": f"{len(matched_trials)} trials matched with mechanism fit"
    }


def answer_question_6_metastasis_prediction(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 6: Metastasis Prediction/Surveillance
    Risk profile?
    """
    resistance = results.get("resistance_detection", {})
    sae_features = results.get("sae_features", {})
    
    # Extract resistance signals
    resistance_signals = []
    if isinstance(resistance, dict):
        signals = resistance.get("signals_detected", [])
        for signal in signals:
            if isinstance(signal, dict):
                resistance_signals.append({
                    "signal_type": signal.get("signal_type", "Unknown"),
                    "detected": signal.get("detected", False),
                    "probability": signal.get("probability", 0.0),
                    "confidence": signal.get("confidence", 0.0)
                })
    
    # Get DNA repair capacity trend (if available)
    dna_repair = sae_features.get("dna_repair_capacity", 0.0) if isinstance(sae_features, dict) else 0.0
    
    return {
        "question": "Metastasis Prediction/Surveillance",
        "answer": {
            "resistance_signals": resistance_signals,
            "dna_repair_capacity": dna_repair,
            "risk_level": "HIGH" if len([s for s in resistance_signals if s.get("detected")]) >= 2 else "MODERATE",
            "proxy_sae_source": "DNA repair capacity trends"
        },
        "summary": f"{len([s for s in resistance_signals if s.get('detected')])} resistance signals detected"
    }


def answer_question_7_immunogenicity(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 7: Immunogenicity & Vaccine Targets
    Neoantigens?
    """
    tumor_context = results.get("tumor_context", {})
    sae_features = results.get("sae_features", {})
    
    tmb = tumor_context.get("tmb_score", 0.0)
    msi_status = tumor_context.get("msi_status", "Unknown")
    io_eligible = sae_features.get("io_eligible", False) if isinstance(sae_features, dict) else False
    
    return {
        "question": "Immunogenicity & Vaccine Targets",
        "answer": {
            "tmb_score": tmb,
            "tmb_status": "HIGH" if tmb >= 20 else "MODERATE" if tmb >= 10 else "LOW",
            "msi_status": msi_status,
            "io_eligible": io_eligible,
            "neoantigen_potential": "HIGH" if tmb >= 20 or msi_status == "MSI-H" else "MODERATE",
            "proxy_sae_source": "IO eligibility from tumor context"
        },
        "summary": f"TMB: {tmb} ({'HIGH' if tmb >= 20 else 'MODERATE' if tmb >= 10 else 'LOW'}), IO Eligible: {io_eligible}"
    }


def answer_question_8_nutritional_therapies(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Question 8: Personalized Nutritional/Adjunctive Therapies
    Diet interventions?
    """
    food_results = results.get("nutritional_therapies", {})
    
    validated_compounds = []
    for compound, food_data in food_results.items():
        if isinstance(food_data, dict):
            spe_score = food_data.get("spe_score", 0.0)
            verdict = food_data.get("verdict", "NOT_SUPPORTED")
            
            validated_compounds.append({
                "compound": compound,
                "spe_score": spe_score,
                "verdict": verdict,
                "recommendation": "SUPPORTED" if verdict == "SUPPORTED" else "CONSIDER" if verdict == "WEAK_SUPPORT" else "NOT_SUPPORTED",
                "proxy_sae_source": "Pathway alignment for compound-disease matching"
            })
    
    return {
        "question": "Personalized Nutritional/Adjunctive Therapies",
        "answer": validated_compounds,
        "summary": f"{len([c for c in validated_compounds if c['recommendation'] == 'SUPPORTED'])} compounds supported"
    }


def answer_all_questions(results_file: Path) -> Dict[str, Any]:
    """
    Load analysis results and answer all 8 clinical questions.
    """
    logger.info(f"📖 Loading analysis results from: {results_file}")
    
    with open(results_file, 'r') as f:
        results = json.load(f)
    
    logger.info("🔬 Answering 8 clinical questions...")
    
    answers = {
        "timestamp": results.get("timestamp"),
        "mutations": results.get("mutations"),
        "questions": [
            answer_question_1_variant_impact(results),
            answer_question_2_functional_annotation(results),
            answer_question_3_pathway_analysis(results),
            answer_question_4_drug_prediction(results),
            answer_question_5_trial_matching(results),
            answer_question_6_metastasis_prediction(results),
            answer_question_7_immunogenicity(results),
            answer_question_8_nutritional_therapies(results)
        ]
    }
    
    return answers


def save_answers(answers: Dict[str, Any], output_path: Path):
    """Save question answers to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(answers, f, indent=2, default=str)
    
    logger.info(f"💾 Answers saved to: {output_path}")


def print_summary(answers: Dict[str, Any]):
    """Print summary of all 8 question answers."""
    print("\n" + "=" * 80)
    print("MBD4+TP53 CLINICAL QUESTIONS - ANSWERS SUMMARY")
    print("=" * 80)
    
    for i, q in enumerate(answers.get("questions", []), 1):
        print(f"\n{i}. {q.get('question')}")
        print(f"   {q.get('summary')}")
    
    print("\n" + "=" * 80)


def main():
    """Main entry point."""
    import sys
    
    if len(sys.argv) < 2:
        # Find most recent analysis file
        results_dir = Path("data/validation/mbd4_tp53_analysis")
        if not results_dir.exists():
            logger.error(f"❌ Results directory not found: {results_dir}")
            logger.info("💡 Run scripts/sae/run_mbd4_tp53_analysis.py first")
            return
        
        # Find most recent JSON file
        json_files = list(results_dir.glob("mbd4_tp53_analysis_*.json"))
        if not json_files:
            logger.error(f"❌ No analysis results found in {results_dir}")
            return
        
        results_file = max(json_files, key=lambda p: p.stat().st_mtime)
        logger.info(f"📂 Using most recent results file: {results_file}")
    else:
        results_file = Path(sys.argv[1])
    
    if not results_file.exists():
        logger.error(f"❌ Results file not found: {results_file}")
        return
    
    # Answer all questions
    answers = answer_all_questions(results_file)
    
    # Save answers
    output_dir = Path("data/validation/mbd4_tp53_analysis")
    output_file = output_dir / f"mbd4_tp53_questions_answered_{answers.get('timestamp', 'unknown').replace(':', '-')}.json"
    save_answers(answers, output_file)
    
    # Print summary
    print_summary(answers)
    
    logger.info(f"\n📄 Full answers saved to: {output_file}")


if __name__ == "__main__":
    main()

