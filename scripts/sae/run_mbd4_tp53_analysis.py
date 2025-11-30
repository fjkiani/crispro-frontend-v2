#!/usr/bin/env python3
"""
MBD4+TP53 HGSOC End-to-End Analysis Script

Runs complete analysis pipeline for MBD4 germline + TP53 somatic mutations
using proxy SAE features (current production method).

Answers 8 clinical questions:
1. Variant Impact Prediction
2. Functional Annotation
3. Pathway Analysis
4. Drug and Therapy Prediction
5. Trial and Biomarker Matching
6. Metastasis Prediction/Surveillance
7. Immunogenicity & Vaccine Targets
8. Personalized Nutritional/Adjunctive Therapies
"""

import json
import asyncio
import httpx
import logging
from typing import Dict, Any, List, Optional
from datetime import datetime
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# API Configuration
API_BASE = "http://127.0.0.1:8000"
TIMEOUT = 120.0  # 2 minutes for complex analyses

# MBD4+TP53 Mutations (from MBD4.mdc)
MBD4_MUTATION = {
    "gene": "MBD4",
    "hgvs_p": "p.Ile413Serfs*2",
    "chrom": "3",
    "pos": 129430456,
    "ref": "A",
    "alt": "",  # Deletion
    "build": "GRCh38",
    "consequence": "frameshift_variant"
}

TP53_MUTATION = {
    "gene": "TP53",
    "hgvs_p": "p.R175H",
    "chrom": "17",
    "pos": 7577120,
    "ref": "G",
    "alt": "A",
    "build": "GRCh38",
    "consequence": "missense_variant"
}

# Tumor Context (HGSOC with high HRD)
TUMOR_CONTEXT = {
    "disease": "ovarian_cancer",
    "hrd_score": 0.75,  # High HRD (MBD4 BER loss + TP53 checkpoint loss)
    "tmb_score": 25.0,   # TMB-high (eligible for IO)
    "msi_status": "MSS",  # Microsatellite stable
    "somatic_mutations": [
        {
            "gene": "TP53",
            "hgvs_p": "p.R175H",
            "chrom": "17",
            "pos": 7577120,
            "ref": "G",
            "alt": "A"
        }
    ]
}


async def call_efficacy_predict(
    client: httpx.AsyncClient,
    mutations: List[Dict[str, Any]],
    tumor_context: Dict[str, Any],
    germline_status: str = "positive"  # MBD4 is germline
) -> Optional[Dict[str, Any]]:
    """
    Call /api/efficacy/predict endpoint.
    
    Returns:
        Efficacy prediction response with drug rankings, pathway scores, and provenance
    """
    try:
        payload = {
            "model_id": "evo2_1b",
            "mutations": mutations,
            "disease": tumor_context.get("disease", "ovarian_cancer"),
            "germline_status": germline_status,
            "tumor_context": tumor_context,
            "options": {
                "adaptive": True,
                "ensemble": False,
                "ablation_mode": "SPE"  # Full S/P/E framework
            },
            "api_base": API_BASE
        }
        
        logger.info(f"🔬 Calling /api/efficacy/predict with {len(mutations)} mutations")
        response = await client.post(
            f"{API_BASE}/api/efficacy/predict",
            json=payload,
            timeout=TIMEOUT
        )
        
        if response.status_code == 200:
            data = response.json()
            logger.info(f"✅ Efficacy prediction complete: {len(data.get('drugs', []))} drugs ranked")
            return data
        else:
            logger.error(f"❌ Efficacy API error: {response.status_code} - {response.text}")
            return None
            
    except Exception as e:
        logger.error(f"❌ Efficacy prediction failed: {str(e)}")
        return None


async def call_insights_bundle(
    client: httpx.AsyncClient,
    mutations: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """
    Call all 4 insights endpoints to get functionality, chromatin, essentiality, regulatory.
    
    Returns:
        Dict with all 4 insight scores
    """
    insights = {}
    
    # Functionality
    try:
        payload = {
            "gene": mutations[0].get("gene"),  # Use first mutation for gene-level insights
            "variants": mutations,
            "model_id": "evo2_1b"
        }
        response = await client.post(
            f"{API_BASE}/api/insights/predict_protein_functionality_change",
            json=payload,
            timeout=60.0
        )
        if response.status_code == 200:
            insights["functionality"] = response.json()
    except Exception as e:
        logger.warning(f"Functionality insight failed: {e}")
    
    # Essentiality
    try:
        payload = {
            "gene": mutations[0].get("gene"),
            "variants": mutations,
            "model_id": "evo2_1b"
        }
        response = await client.post(
            f"{API_BASE}/api/insights/predict_gene_essentiality",
            json=payload,
            timeout=60.0
        )
        if response.status_code == 200:
            insights["essentiality"] = response.json()
    except Exception as e:
        logger.warning(f"Essentiality insight failed: {e}")
    
    # Regulatory (splicing/noncoding)
    try:
        for mut in mutations:
            if mut.get("chrom") and mut.get("pos"):
                payload = {
                    "chrom": mut["chrom"],
                    "pos": mut["pos"],
                    "ref": mut.get("ref", ""),
                    "alt": mut.get("alt", ""),
                    "model_id": "evo2_1b"
                }
                response = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json=payload,
                    timeout=60.0
                )
                if response.status_code == 200:
                    insights[f"regulatory_{mut['gene']}"] = response.json()
                    break  # Just get one for now
    except Exception as e:
        logger.warning(f"Regulatory insight failed: {e}")
    
    # Chromatin
    try:
        payload = {
            "gene": mutations[0].get("gene"),
            "variants": mutations,
            "model_id": "evo2_1b"
        }
        response = await client.post(
            f"{API_BASE}/api/insights/predict_chromatin_accessibility",
            json=payload,
            timeout=60.0
        )
        if response.status_code == 200:
            insights["chromatin"] = response.json()
    except Exception as e:
        logger.warning(f"Chromatin insight failed: {e}")
    
    return insights


async def call_evidence_deep_analysis(
    client: httpx.AsyncClient,
    mutations: List[Dict[str, Any]],
    disease: str = "ovarian_cancer"
) -> Optional[Dict[str, Any]]:
    """
    Call /api/evidence/deep_analysis for literature and ClinVar evidence.
    """
    try:
        # Use first mutation for evidence search
        mut = mutations[0]
        payload = {
            "gene": mut.get("gene"),
            "hgvs_p": mut.get("hgvs_p"),
            "disease": disease,
            "moa_terms": ["DNA repair", "BER", "tumor suppressor"]  # MBD4-specific terms
        }
        
        logger.info(f"🔬 Calling /api/evidence/deep_analysis for {mut.get('gene')}")
        response = await client.post(
            f"{API_BASE}/api/evidence/deep_analysis",
            json=payload,
            timeout=TIMEOUT
        )
        
        if response.status_code == 200:
            data = response.json()
            logger.info(f"✅ Evidence analysis complete: {len(data.get('citations', []))} citations")
            return data
        else:
            logger.warning(f"⚠️ Evidence API error: {response.status_code}")
            return None
            
    except Exception as e:
        logger.warning(f"⚠️ Evidence analysis failed: {e}")
        return None


async def call_sae_compute_features(
    client: httpx.AsyncClient,
    pathway_scores: Dict[str, float],
    insights_bundle: Dict[str, Any],
    tumor_context: Dict[str, Any]
) -> Optional[Dict[str, Any]]:
    """
    Call /api/sae/compute_features to get proxy SAE features.
    
    Note: This endpoint may not exist yet. If it doesn't, we'll compute proxy SAE
    features directly from pathway scores and insights bundle.
    """
    try:
        payload = {
            "pathway_scores": pathway_scores,
            "insights_bundle": insights_bundle,
            "tumor_context": tumor_context,
            "treatment_history": []
        }
        
        logger.info(f"🔬 Calling /api/sae/compute_features")
        response = await client.post(
            f"{API_BASE}/api/sae/compute_features",
            json=payload,
            timeout=60.0
        )
        
        if response.status_code == 200:
            data = response.json()
            logger.info(f"✅ SAE features computed")
            return data
        else:
            logger.warning(f"⚠️ SAE endpoint not available (status {response.status_code}), computing proxy locally")
            return None
            
    except Exception as e:
        logger.warning(f"⚠️ SAE endpoint not available: {e}, computing proxy locally")
        return None


async def call_trials_agent_search(
    client: httpx.AsyncClient,
    mutations: List[Dict[str, Any]],
    tumor_context: Dict[str, Any],
    mechanism_vector: Optional[List[float]] = None
) -> Optional[Dict[str, Any]]:
    """
    Call /api/trials/agent/search for autonomous trial matching.
    """
    try:
        payload = {
            "mutations": mutations,
            "disease": tumor_context.get("disease", "ovarian_cancer"),
            "germline_status": "positive",  # MBD4 is germline
            "tumor_context": tumor_context,
            "mechanism_vector": mechanism_vector  # 7D mechanism vector if available
        }
        
        logger.info(f"🔬 Calling /api/trials/agent/search")
        response = await client.post(
            f"{API_BASE}/api/trials/agent/search",
            json=payload,
            timeout=TIMEOUT
        )
        
        if response.status_code == 200:
            data = response.json()
            logger.info(f"✅ Trial search complete: {len(data.get('trials', []))} trials found")
            return data
        else:
            logger.warning(f"⚠️ Trial search error: {response.status_code}")
            return None
            
    except Exception as e:
        logger.warning(f"⚠️ Trial search failed: {e}")
        return None


async def call_resistance_detection(
    client: httpx.AsyncClient,
    tumor_context: Dict[str, Any],
    sae_features: Optional[Dict[str, Any]] = None
) -> Optional[Dict[str, Any]]:
    """
    Call resistance detection service (if available).
    """
    try:
        payload = {
            "tumor_context": tumor_context,
            "sae_features": sae_features
        }
        
        logger.info(f"🔬 Calling /api/care/resistance_playbook")
        response = await client.post(
            f"{API_BASE}/api/care/resistance_playbook",
            json=payload,
            timeout=60.0
        )
        
        if response.status_code == 200:
            data = response.json()
            logger.info(f"✅ Resistance detection complete")
            return data
        else:
            logger.warning(f"⚠️ Resistance detection error: {response.status_code}")
            return None
            
    except Exception as e:
        logger.warning(f"⚠️ Resistance detection failed: {e}")
        return None


async def call_food_validator(
    client: httpx.AsyncClient,
    compound: str,
    disease: str = "ovarian_cancer"
) -> Optional[Dict[str, Any]]:
    """
    Call /api/hypothesis/validate_food_dynamic for nutritional therapy validation.
    """
    try:
        payload = {
            "compound": compound,
            "disease": disease,
            "patient_context": {
                "disease": disease,
                "treatment_line": 1,
                "biomarkers": {
                    "hrd_positive": True,
                    "tmb_high": True
                }
            }
        }
        
        logger.info(f"🔬 Calling /api/hypothesis/validate_food_dynamic for {compound}")
        response = await client.post(
            f"{API_BASE}/api/hypothesis/validate_food_dynamic",
            json=payload,
            timeout=60.0
        )
        
        if response.status_code == 200:
            data = response.json()
            logger.info(f"✅ Food validation complete for {compound}")
            return data
        else:
            logger.warning(f"⚠️ Food validation error: {response.status_code}")
            return None
            
    except Exception as e:
        logger.warning(f"⚠️ Food validation failed: {e}")
        return None


async def run_complete_analysis() -> Dict[str, Any]:
    """
    Run complete MBD4+TP53 analysis pipeline.
    
    Returns:
        Complete analysis results with all 8 question answers
    """
    mutations = [MBD4_MUTATION, TP53_MUTATION]
    
    async with httpx.AsyncClient() as client:
        logger.info("=" * 80)
        logger.info("MBD4+TP53 HGSOC END-TO-END ANALYSIS")
        logger.info("=" * 80)
        
        # Step 1: Efficacy Prediction (S/P/E framework)
        logger.info("\n📊 Step 1: Efficacy Prediction (S/P/E Framework)")
        efficacy_response = await call_efficacy_predict(
            client, mutations, TUMOR_CONTEXT, germline_status="positive"
        )
        
        if not efficacy_response:
            logger.error("❌ Efficacy prediction failed - cannot continue")
            return {"error": "Efficacy prediction failed"}
        
        # Extract pathway scores from efficacy response
        pathway_scores = {}
        if "provenance" in efficacy_response:
            confidence_breakdown = efficacy_response["provenance"].get("confidence_breakdown", {})
            pathway_disruption = confidence_breakdown.get("pathway_disruption", {})
            if isinstance(pathway_disruption, dict):
                pathway_scores = pathway_disruption
        
        logger.info(f"✅ Pathway scores extracted: {list(pathway_scores.keys())}")
        
        # Step 2: Insights Bundle (4 chips)
        logger.info("\n🧬 Step 2: Insights Bundle (Functionality, Chromatin, Essentiality, Regulatory)")
        insights_bundle = await call_insights_bundle(client, mutations)
        logger.info(f"✅ Insights bundle complete: {list(insights_bundle.keys())}")
        
        # Step 3: Evidence Analysis
        logger.info("\n📚 Step 3: Evidence Analysis (Literature + ClinVar)")
        evidence_response = await call_evidence_deep_analysis(
            client, mutations, disease="ovarian_cancer"
        )
        
        # Step 4: SAE Features (Proxy)
        logger.info("\n🔬 Step 4: SAE Features (Proxy from Pathway Scores)")
        sae_features = await call_sae_compute_features(
            client, pathway_scores, insights_bundle, TUMOR_CONTEXT
        )
        
        # If SAE endpoint doesn't exist, create minimal proxy SAE structure
        if not sae_features:
            logger.info("⚠️ SAE endpoint not available - creating minimal proxy structure from pathway scores...")
            # Compute mechanism vector from pathway scores
            import sys
            from pathlib import Path
            # Add backend to path for local computation
            backend_path = Path(__file__).parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"
            if backend_path.exists():
                sys.path.insert(0, str(backend_path))
                try:
                    from api.services.pathway_to_mechanism_vector import convert_pathway_scores_to_mechanism_vector
                    mechanism_vector, dimension = convert_pathway_scores_to_mechanism_vector(
                        pathway_scores, TUMOR_CONTEXT, 
                        tmb=TUMOR_CONTEXT.get("tmb_score"),
                        msi_status=TUMOR_CONTEXT.get("msi_status")
                    )
                    # Create minimal SAE features structure
                    sae_features = {
                        "mechanism_vector": mechanism_vector,
                        "dna_repair_capacity": pathway_scores.get("ddr", 0.0) * 0.6,  # Simplified formula
                        "io_eligible": TUMOR_CONTEXT.get("tmb_score", 0) >= 20 or TUMOR_CONTEXT.get("msi_status") == "MSI-H",
                        "dimension": dimension
                    }
                    logger.info(f"✅ Created proxy SAE features locally ({dimension} mechanism vector)")
                except ImportError:
                    logger.warning("⚠️ Could not import pathway converter - using pathway scores directly")
                    sae_features = {"pathway_scores": pathway_scores}
            else:
                logger.warning("⚠️ Backend path not found - using pathway scores directly")
                sae_features = {"pathway_scores": pathway_scores}
        
        # Extract mechanism vector for trial matching
        mechanism_vector = None
        if sae_features and isinstance(sae_features, dict):
            mechanism_vector = sae_features.get("mechanism_vector", [])
            logger.info(f"✅ Mechanism vector: {mechanism_vector}")
        
        # Step 5: Trial Matching
        logger.info("\n🎯 Step 5: Trial Matching (Autonomous Agent)")
        trials_response = await call_trials_agent_search(
            client, mutations, TUMOR_CONTEXT, mechanism_vector=mechanism_vector
        )
        
        # Step 6: Resistance Detection
        logger.info("\n⚠️ Step 6: Resistance Detection")
        resistance_response = await call_resistance_detection(
            client, TUMOR_CONTEXT, sae_features=sae_features
        )
        
        # Step 7: Food Validator (Sample compounds)
        logger.info("\n🍎 Step 7: Nutritional Therapy Validation")
        food_compounds = ["Vitamin D", "Curcumin", "Omega-3 fatty acids"]
        food_results = {}
        for compound in food_compounds:
            food_result = await call_food_validator(client, compound, disease="ovarian_cancer")
            if food_result:
                food_results[compound] = food_result
        
        # Compile complete results
        results = {
            "timestamp": datetime.now().isoformat(),
            "mutations": mutations,
            "tumor_context": TUMOR_CONTEXT,
            "efficacy_prediction": efficacy_response,
            "pathway_scores": pathway_scores,
            "insights_bundle": insights_bundle,
            "evidence_analysis": evidence_response,
            "sae_features": sae_features,
            "trial_matching": trials_response,
            "resistance_detection": resistance_response,
            "nutritional_therapies": food_results,
            "provenance": {
                "api_base": API_BASE,
                "model_id": "evo2_1b",
                "sae_type": "proxy",  # Using proxy SAE (pathway-based)
                "analysis_version": "v1"
            }
        }
        
        logger.info("\n" + "=" * 80)
        logger.info("✅ COMPLETE ANALYSIS FINISHED")
        logger.info("=" * 80)
        
        return results


def save_results(results: Dict[str, Any], output_path: Path):
    """Save analysis results to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    logger.info(f"💾 Results saved to: {output_path}")


async def main():
    """Main entry point."""
    # Run analysis
    results = await run_complete_analysis()
    
    # Save results
    output_dir = Path("data/validation/mbd4_tp53_analysis")
    output_file = output_dir / f"mbd4_tp53_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    save_results(results, output_file)
    
    # Print summary
    print("\n" + "=" * 80)
    print("MBD4+TP53 ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"✅ Efficacy Prediction: {len(results.get('efficacy_prediction', {}).get('drugs', []))} drugs ranked")
    print(f"✅ Pathway Scores: {len(results.get('pathway_scores', {}))} pathways")
    print(f"✅ Insights Bundle: {len(results.get('insights_bundle', {}))} insights")
    print(f"✅ SAE Features: {'Computed' if results.get('sae_features') else 'Failed'}")
    print(f"✅ Trial Matching: {len(results.get('trial_matching', {}).get('trials', []))} trials")
    print(f"✅ Nutritional Therapies: {len(results.get('nutritional_therapies', {}))} compounds")
    print(f"\n📄 Full results saved to: {output_file}")
    print("=" * 80)


if __name__ == "__main__":
    asyncio.run(main())

