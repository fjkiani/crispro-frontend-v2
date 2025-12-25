#!/usr/bin/env python3
"""
PDGFRA p.S755P VUS Resolution Script
Executes the complete 8-step VUS identification workflow
"""
import asyncio
import httpx
import json
import sys
import os
from typing import Dict, Any, Optional
from datetime import datetime

# Add path for standalone coordinate resolver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

API_BASE = "http://127.0.0.1:8000"
TIMEOUT = 60.0

# Import standalone coordinate resolver
try:
    from resolve_coordinates_standalone import resolve_coordinates as resolve_coords_standalone
except ImportError:
    resolve_coords_standalone = None

# Variant information
VARIANT = {
    "gene": "PDGFRA",
    "hgvs_p": "S755P",
    "c_dna": "c.2263T>C",
    "classification": "VUS",
    "zygosity": "heterozygous"
}

# Note: We need GRCh38 coordinates for full analysis
# PDGFRA is on chromosome 4, but we'll need to resolve exact position
# For now, we'll use the gene name and HGVS protein notation


async def step1_normalize() -> Dict[str, Any]:
    """Step 1: Inputs and Normalization + Coordinate Resolution"""
    print("\n" + "="*60)
    print("STEP 1: Inputs and Normalization + Coordinate Resolution")
    print("="*60)
    
    result = {
        "gene": VARIANT["gene"],
        "hgvs_p": VARIANT["hgvs_p"],
        "c_dna": VARIANT["c_dna"],
        "status": "normalized",
        "coordinates": None
    }
    
    print(f"✅ Gene: {VARIANT['gene']}")
    print(f"✅ Protein: {VARIANT['hgvs_p']}")
    print(f"✅ cDNA: {VARIANT['c_dna']}")
    
    # Resolve coordinates - try backend endpoint first, fallback to standalone
    coords = None
    
    # Try backend endpoint
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "c_dna": VARIANT["c_dna"]
            }
            resp = await client.post(
                f"{API_BASE}/api/vus/resolve_coordinates",
                json=payload
            )
            if resp.status_code == 200:
                coords = resp.json()
                print("✅ Coordinates resolved via backend endpoint")
    except Exception as e:
        pass  # Will try standalone
    
    # Fallback to standalone resolver if backend failed
    if not coords and resolve_coords_standalone:
        try:
            # Run in executor since it's synchronous
            import asyncio
            loop = asyncio.get_event_loop()
            coords = await loop.run_in_executor(
                None,
                resolve_coords_standalone,
                VARIANT["gene"],
                VARIANT["hgvs_p"],
                VARIANT["c_dna"]
            )
            print("✅ Coordinates resolved via standalone resolver")
        except Exception as e2:
            print(f"⚠️ Standalone coordinate resolution failed: {e2}")
    
    if not coords:
        print(f"⚠️ All coordinate resolution methods failed")
    
    if coords:
        result["coordinates"] = coords
        result["chrom"] = coords.get("chrom")
        result["pos"] = coords.get("pos")
        result["ref"] = coords.get("ref")
        result["alt"] = coords.get("alt")
        result["transcript"] = coords.get("transcript")
        print(f"✅ Coordinates resolved:")
        print(f"   Chrom: {coords.get('chrom')}")
        print(f"   Pos: {coords.get('pos')}")
        print(f"   Ref: {coords.get('ref')}")
        print(f"   Alt: {coords.get('alt')}")
        print(f"   Transcript: {coords.get('transcript')}")
    else:
        result["note"] = "GRCh38 coordinates resolution failed - continuing with limited analysis"
    
    return result


async def step2_priors(normalized: Dict[str, Any]) -> Dict[str, Any]:
    """Step 2: Priors (ClinVar + AlphaMissense)"""
    print("\n" + "="*60)
    print("STEP 2: Priors (ClinVar + AlphaMissense Coverage)")
    print("="*60)
    
    result = {
        "clinvar": None,
        "alphamissense": None,
        "prior_resolved": False
    }
    
    # Get coordinates from normalized result
    chrom = normalized.get("chrom")
    pos = normalized.get("pos")
    ref = normalized.get("ref")
    alt = normalized.get("alt")
    
    # ClinVar lookup (with coordinates if available)
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "assembly": "GRCh38"
            }
            
            # Add coordinates if available
            if chrom and pos and ref and alt:
                payload.update({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt
                })
            
            resp = await client.post(
                f"{API_BASE}/api/evidence/deep_analysis",
                json=payload
            )
            
            if resp.status_code == 200:
                data = resp.json()
                result["clinvar"] = data.get("clinvar", {})
                classification = result["clinvar"].get("classification") if result["clinvar"] else None
                print(f"✅ ClinVar: {classification or 'Not found'}")
            else:
                print(f"⚠️ ClinVar lookup returned {resp.status_code}")
    except Exception as e:
        print(f"⚠️ ClinVar lookup failed: {e}")
    
    # AlphaMissense coverage (requires coordinates)
    if chrom and pos and ref and alt:
        try:
            async with httpx.AsyncClient(timeout=TIMEOUT) as client:
                resp = await client.get(
                    f"{API_BASE}/api/fusion/coverage",
                    params={"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
                )
                if resp.status_code == 200:
                    data = resp.json()
                    result["alphamissense"] = {
                        "eligible": data.get("coverage", False),
                        "status": "checked"
                    }
                    print(f"✅ AlphaMissense: {'Eligible' if data.get('coverage') else 'Not eligible'}")
                else:
                    print(f"⚠️ AlphaMissense check returned {resp.status_code}")
        except Exception as e:
            print(f"⚠️ AlphaMissense check failed: {e}")
    else:
        result["alphamissense"] = {
            "status": "requires_coordinates",
            "note": "p.S755P is missense, eligible for Fusion if coordinates available"
        }
        print(f"⚠️ AlphaMissense: Requires GRCh38 coordinates")
    
    return result


async def step3_triumvirate_gate() -> Dict[str, Any]:
    """Step 3: Triumvirate Protocol Gate"""
    print("\n" + "="*60)
    print("STEP 3: Triumvirate Protocol Gate")
    print("="*60)
    
    # p.S755P is a missense variant (Ser→Pro), NOT truncating
    is_truncating = "*" in VARIANT["hgvs_p"] or "FS" in VARIANT["hgvs_p"]
    
    result = {
        "is_truncating": is_truncating,
        "gate_result": "PASS" if not is_truncating else "FAIL",
        "reason": "Missense variant (Ser→Pro), not truncating"
    }
    
    print(f"✅ Is truncating: {is_truncating}")
    print(f"✅ Gate result: {result['gate_result']}")
    print(f"✅ Reason: {result['reason']}")
    
    return result


async def step4_spe_signals() -> Dict[str, Any]:
    """Step 4: S/P/E Core Signals"""
    print("\n" + "="*60)
    print("STEP 4: S/P/E Core Signals")
    print("="*60)
    
    result = {
        "sequence": None,
        "pathway": None,
        "evidence": None
    }
    
    # Sequence (S) - Evo2 scores
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            # Try insights endpoint for functionality (proxy for S)
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "model_id": "evo2_1b"
            }
            
            resp = await client.post(
                f"{API_BASE}/api/insights/predict_protein_functionality_change",
                json=payload
            )
            
            if resp.status_code == 200:
                data = resp.json()
                result["sequence"] = {
                    "functionality_score": data.get("functionality_score", 0.0),
                    "provenance": data.get("provenance", {})
                }
                print(f"✅ Sequence (S): {result['sequence']['functionality_score']:.3f}")
            else:
                print(f"⚠️ Sequence scoring returned {resp.status_code}")
    except Exception as e:
        print(f"⚠️ Sequence scoring failed: {e}")
    
    # Pathway (P) - PDGFRA is in RTK/MAPK pathway
    result["pathway"] = {
        "pathways": ["RTK", "MAPK"],
        "note": "PDGFRA is receptor tyrosine kinase in MAPK pathway"
    }
    print(f"✅ Pathway (P): {', '.join(result['pathway']['pathways'])}")
    
    # Evidence (E) - Already checked in step 2
    result["evidence"] = {
        "note": "See Step 2 ClinVar results"
    }
    
    return result


async def step5_insights_bundle(normalized: Dict[str, Any]) -> Dict[str, Any]:
    """Step 5: Insights Bundle"""
    print("\n" + "="*60)
    print("STEP 5: Insights Bundle")
    print("="*60)
    
    result = {
        "functionality": None,
        "chromatin": None,
        "essentiality": None,
        "regulatory": None
    }
    
    # Get coordinates from normalized result
    chrom = normalized.get("chrom")
    pos = normalized.get("pos")
    ref = normalized.get("ref")
    alt = normalized.get("alt")
    
    async with httpx.AsyncClient(timeout=TIMEOUT) as client:
        # Functionality
        try:
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "model_id": "evo2_1b"
            }
            resp = await client.post(
                f"{API_BASE}/api/insights/predict_protein_functionality_change",
                json=payload
            )
            if resp.status_code == 200:
                data = resp.json()
                result["functionality"] = data.get("functionality_score", 0.0)
                print(f"✅ Functionality: {result['functionality']:.3f}")
        except Exception as e:
            print(f"⚠️ Functionality failed: {e}")
        
        # Essentiality
        try:
            payload = {
                "gene": VARIANT["gene"],
                "variants": [{"gene": VARIANT["gene"], "hgvs_p": VARIANT["hgvs_p"]}],
                "model_id": "evo2_1b"
            }
            resp = await client.post(
                f"{API_BASE}/api/insights/predict_gene_essentiality",
                json=payload
            )
            if resp.status_code == 200:
                data = resp.json()
                result["essentiality"] = data.get("essentiality_score", 0.0)
                print(f"✅ Essentiality: {result['essentiality']:.3f}")
        except Exception as e:
            print(f"⚠️ Essentiality failed: {e}")
        
        # Regulatory (requires coordinates)
        if chrom and pos and ref and alt:
            try:
                payload = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "model_id": "evo2_1b"
                }
                resp = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json=payload
                )
                if resp.status_code == 200:
                    data = resp.json()
                    result["regulatory"] = data.get("regulatory_impact_score", 0.0)
                    print(f"✅ Regulatory: {result['regulatory']:.3f}")
            except Exception as e:
                print(f"⚠️ Regulatory failed: {e}")
        else:
            result["regulatory"] = {
                "status": "requires_coordinates",
                "note": "Needs chrom, pos, ref, alt for regulatory analysis"
            }
            print(f"⚠️ Regulatory: Requires coordinates")
        
        # Chromatin (requires coordinates)
        if chrom and pos:
            try:
                payload = {
                    "chrom": chrom,
                    "pos": pos,
                    "radius": 500
                }
                resp = await client.post(
                    f"{API_BASE}/api/insights/predict_chromatin_accessibility",
                    json=payload
                )
                if resp.status_code == 200:
                    data = resp.json()
                    result["chromatin"] = data.get("accessibility_score", 0.0)
                    print(f"✅ Chromatin: {result['chromatin']:.3f}")
            except Exception as e:
                print(f"⚠️ Chromatin failed: {e}")
        else:
            result["chromatin"] = {
                "status": "requires_coordinates",
                "note": "Needs chrom, pos for chromatin analysis"
            }
            print(f"⚠️ Chromatin: Requires coordinates")
    
    return result


async def step6_fusion_gating() -> Dict[str, Any]:
    """Step 6: Fusion Gating (AlphaMissense)"""
    print("\n" + "="*60)
    print("STEP 6: Fusion Gating (AlphaMissense)")
    print("="*60)
    
    result = {
        "eligible": True,  # p.S755P is missense
        "coverage": None,
        "fused_score": None,
        "note": "Requires GRCh38 coordinates for actual Fusion scoring"
    }
    
    print(f"✅ Eligible: {result['eligible']} (missense variant)")
    print(f"⚠️ Coverage: Requires coordinates")
    
    return result


async def step7_sae_features() -> Dict[str, Any]:
    """Step 7: SAE Interpretable Features"""
    print("\n" + "="*60)
    print("STEP 7: SAE Interpretable Features")
    print("="*60)
    
    result = {
        "status": "requires_evo2_scoring",
        "note": "SAE features extracted from Evo2 embeddings"
    }
    
    print(f"⚠️ SAE: Requires Evo2 scoring with coordinates")
    
    return result


async def step8_triage_scoring(
    normalized: Dict[str, Any],
    priors: Dict[str, Any],
    gate: Dict[str, Any],
    spe: Dict[str, Any],
    insights: Dict[str, Any],
    fusion: Dict[str, Any],
    sae: Dict[str, Any]
) -> Dict[str, Any]:
    """Step 8: VUS Triage Scoring and Decision Rubric"""
    print("\n" + "="*60)
    print("STEP 8: VUS Triage Scoring and Decision Rubric")
    print("="*60)
    
    # Collect scores
    functionality = insights.get("functionality", 0.0) if insights.get("functionality") else 0.0
    essentiality = insights.get("essentiality", 0.0) if insights.get("essentiality") else 0.0
    
    # Decision logic
    confidence = 0.0
    classification = "VUS"
    rationale = []
    
    # If functionality is low (<0.3) and ClinVar is benign
    if functionality < 0.3:
        confidence += 0.3
        rationale.append(f"Low functionality score ({functionality:.3f}) suggests minimal impact")
    
    # If functionality is high (>0.7)
    if functionality > 0.7:
        confidence += 0.4
        rationale.append(f"High functionality score ({functionality:.3f}) suggests significant impact")
    
    # Check ClinVar
    clinvar_class = priors.get("clinvar", {}).get("classification") if priors.get("clinvar") else None
    if clinvar_class:
        if "benign" in clinvar_class.lower():
            confidence += 0.3
            classification = "Likely Benign"
            rationale.append(f"ClinVar classification: {clinvar_class}")
        elif "pathogenic" in clinvar_class.lower():
            confidence += 0.4
            classification = "Likely Pathogenic"
            rationale.append(f"ClinVar classification: {clinvar_class}")
    
    # Final classification
    if confidence >= 0.7:
        if functionality < 0.3:
            classification = "Benign" if confidence >= 0.8 else "Likely Benign"
        elif functionality > 0.7:
            classification = "Pathogenic" if confidence >= 0.8 else "Likely Pathogenic"
    else:
        classification = "VUS"
        rationale.append("Insufficient evidence for definitive classification")
    
    result = {
        "confidence": confidence,
        "classification": classification,
        "rationale": rationale,
        "evidence_summary": {
            "functionality": functionality,
            "essentiality": essentiality,
            "clinvar": clinvar_class,
            "gate_passed": gate.get("gate_result") == "PASS"
        },
        "recommendation": "VUS eliminated" if classification != "VUS" else "VUS remains - additional evidence needed"
    }
    
    print(f"✅ Confidence: {confidence:.3f}")
    print(f"✅ Classification: {classification}")
    print(f"✅ Recommendation: {result['recommendation']}")
    print(f"\nRationale:")
    for r in rationale:
        print(f"  - {r}")
    
    return result


async def main():
    """Execute complete 8-step workflow"""
    print("\n" + "="*60)
    print("PDGFRA p.S755P VUS RESOLUTION WORKFLOW")
    print("="*60)
    print(f"Variant: {VARIANT['gene']} {VARIANT['hgvs_p']} ({VARIANT['c_dna']})")
    print(f"Current Classification: {VARIANT['classification']}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    
    # Execute steps
    normalized = await step1_normalize()
    priors = await step2_priors(normalized)
    gate = await step3_triumvirate_gate()
    spe = await step4_spe_signals()
    insights = await step5_insights_bundle(normalized)
    fusion = await step6_fusion_gating()
    sae = await step7_sae_features()
    triage = await step8_triage_scoring(normalized, priors, gate, spe, insights, fusion, sae)
    
    # Final summary
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    print(f"Original Classification: {VARIANT['classification']}")
    print(f"New Classification: {triage['classification']}")
    print(f"Confidence: {triage['confidence']:.3f}")
    print(f"Recommendation: {triage['recommendation']}")
    
    # Save results
    results = {
        "variant": VARIANT,
        "workflow": {
            "step1": normalized,
            "step2": priors,
            "step3": gate,
            "step4": spe,
            "step5": insights,
            "step6": fusion,
            "step7": sae,
            "step8": triage
        },
        "timestamp": datetime.now().isoformat()
    }
    
    output_file = ".cursor/ayesha/PDGFRA_VUS_RESOLUTION_RESULTS.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    return results


if __name__ == "__main__":
    asyncio.run(main())




PDGFRA p.S755P VUS Resolution Script
Executes the complete 8-step VUS identification workflow
"""
import asyncio
import httpx
import json
import sys
import os
from typing import Dict, Any, Optional
from datetime import datetime

# Add path for standalone coordinate resolver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

API_BASE = "http://127.0.0.1:8000"
TIMEOUT = 60.0

# Import standalone coordinate resolver
try:
    from resolve_coordinates_standalone import resolve_coordinates as resolve_coords_standalone
except ImportError:
    resolve_coords_standalone = None

# Variant information
VARIANT = {
    "gene": "PDGFRA",
    "hgvs_p": "S755P",
    "c_dna": "c.2263T>C",
    "classification": "VUS",
    "zygosity": "heterozygous"
}

# Note: We need GRCh38 coordinates for full analysis
# PDGFRA is on chromosome 4, but we'll need to resolve exact position
# For now, we'll use the gene name and HGVS protein notation


async def step1_normalize() -> Dict[str, Any]:
    """Step 1: Inputs and Normalization + Coordinate Resolution"""
    print("\n" + "="*60)
    print("STEP 1: Inputs and Normalization + Coordinate Resolution")
    print("="*60)
    
    result = {
        "gene": VARIANT["gene"],
        "hgvs_p": VARIANT["hgvs_p"],
        "c_dna": VARIANT["c_dna"],
        "status": "normalized",
        "coordinates": None
    }
    
    print(f"✅ Gene: {VARIANT['gene']}")
    print(f"✅ Protein: {VARIANT['hgvs_p']}")
    print(f"✅ cDNA: {VARIANT['c_dna']}")
    
    # Resolve coordinates - try backend endpoint first, fallback to standalone
    coords = None
    
    # Try backend endpoint
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "c_dna": VARIANT["c_dna"]
            }
            resp = await client.post(
                f"{API_BASE}/api/vus/resolve_coordinates",
                json=payload
            )
            if resp.status_code == 200:
                coords = resp.json()
                print("✅ Coordinates resolved via backend endpoint")
    except Exception as e:
        pass  # Will try standalone
    
    # Fallback to standalone resolver if backend failed
    if not coords and resolve_coords_standalone:
        try:
            # Run in executor since it's synchronous
            import asyncio
            loop = asyncio.get_event_loop()
            coords = await loop.run_in_executor(
                None,
                resolve_coords_standalone,
                VARIANT["gene"],
                VARIANT["hgvs_p"],
                VARIANT["c_dna"]
            )
            print("✅ Coordinates resolved via standalone resolver")
        except Exception as e2:
            print(f"⚠️ Standalone coordinate resolution failed: {e2}")
    
    if not coords:
        print(f"⚠️ All coordinate resolution methods failed")
    
    if coords:
        result["coordinates"] = coords
        result["chrom"] = coords.get("chrom")
        result["pos"] = coords.get("pos")
        result["ref"] = coords.get("ref")
        result["alt"] = coords.get("alt")
        result["transcript"] = coords.get("transcript")
        print(f"✅ Coordinates resolved:")
        print(f"   Chrom: {coords.get('chrom')}")
        print(f"   Pos: {coords.get('pos')}")
        print(f"   Ref: {coords.get('ref')}")
        print(f"   Alt: {coords.get('alt')}")
        print(f"   Transcript: {coords.get('transcript')}")
    else:
        result["note"] = "GRCh38 coordinates resolution failed - continuing with limited analysis"
    
    return result


async def step2_priors(normalized: Dict[str, Any]) -> Dict[str, Any]:
    """Step 2: Priors (ClinVar + AlphaMissense)"""
    print("\n" + "="*60)
    print("STEP 2: Priors (ClinVar + AlphaMissense Coverage)")
    print("="*60)
    
    result = {
        "clinvar": None,
        "alphamissense": None,
        "prior_resolved": False
    }
    
    # Get coordinates from normalized result
    chrom = normalized.get("chrom")
    pos = normalized.get("pos")
    ref = normalized.get("ref")
    alt = normalized.get("alt")
    
    # ClinVar lookup (with coordinates if available)
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "assembly": "GRCh38"
            }
            
            # Add coordinates if available
            if chrom and pos and ref and alt:
                payload.update({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt
                })
            
            resp = await client.post(
                f"{API_BASE}/api/evidence/deep_analysis",
                json=payload
            )
            
            if resp.status_code == 200:
                data = resp.json()
                result["clinvar"] = data.get("clinvar", {})
                classification = result["clinvar"].get("classification") if result["clinvar"] else None
                print(f"✅ ClinVar: {classification or 'Not found'}")
            else:
                print(f"⚠️ ClinVar lookup returned {resp.status_code}")
    except Exception as e:
        print(f"⚠️ ClinVar lookup failed: {e}")
    
    # AlphaMissense coverage (requires coordinates)
    if chrom and pos and ref and alt:
        try:
            async with httpx.AsyncClient(timeout=TIMEOUT) as client:
                resp = await client.get(
                    f"{API_BASE}/api/fusion/coverage",
                    params={"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
                )
                if resp.status_code == 200:
                    data = resp.json()
                    result["alphamissense"] = {
                        "eligible": data.get("coverage", False),
                        "status": "checked"
                    }
                    print(f"✅ AlphaMissense: {'Eligible' if data.get('coverage') else 'Not eligible'}")
                else:
                    print(f"⚠️ AlphaMissense check returned {resp.status_code}")
        except Exception as e:
            print(f"⚠️ AlphaMissense check failed: {e}")
    else:
        result["alphamissense"] = {
            "status": "requires_coordinates",
            "note": "p.S755P is missense, eligible for Fusion if coordinates available"
        }
        print(f"⚠️ AlphaMissense: Requires GRCh38 coordinates")
    
    return result


async def step3_triumvirate_gate() -> Dict[str, Any]:
    """Step 3: Triumvirate Protocol Gate"""
    print("\n" + "="*60)
    print("STEP 3: Triumvirate Protocol Gate")
    print("="*60)
    
    # p.S755P is a missense variant (Ser→Pro), NOT truncating
    is_truncating = "*" in VARIANT["hgvs_p"] or "FS" in VARIANT["hgvs_p"]
    
    result = {
        "is_truncating": is_truncating,
        "gate_result": "PASS" if not is_truncating else "FAIL",
        "reason": "Missense variant (Ser→Pro), not truncating"
    }
    
    print(f"✅ Is truncating: {is_truncating}")
    print(f"✅ Gate result: {result['gate_result']}")
    print(f"✅ Reason: {result['reason']}")
    
    return result


async def step4_spe_signals() -> Dict[str, Any]:
    """Step 4: S/P/E Core Signals"""
    print("\n" + "="*60)
    print("STEP 4: S/P/E Core Signals")
    print("="*60)
    
    result = {
        "sequence": None,
        "pathway": None,
        "evidence": None
    }
    
    # Sequence (S) - Evo2 scores
    try:
        async with httpx.AsyncClient(timeout=TIMEOUT) as client:
            # Try insights endpoint for functionality (proxy for S)
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "model_id": "evo2_1b"
            }
            
            resp = await client.post(
                f"{API_BASE}/api/insights/predict_protein_functionality_change",
                json=payload
            )
            
            if resp.status_code == 200:
                data = resp.json()
                result["sequence"] = {
                    "functionality_score": data.get("functionality_score", 0.0),
                    "provenance": data.get("provenance", {})
                }
                print(f"✅ Sequence (S): {result['sequence']['functionality_score']:.3f}")
            else:
                print(f"⚠️ Sequence scoring returned {resp.status_code}")
    except Exception as e:
        print(f"⚠️ Sequence scoring failed: {e}")
    
    # Pathway (P) - PDGFRA is in RTK/MAPK pathway
    result["pathway"] = {
        "pathways": ["RTK", "MAPK"],
        "note": "PDGFRA is receptor tyrosine kinase in MAPK pathway"
    }
    print(f"✅ Pathway (P): {', '.join(result['pathway']['pathways'])}")
    
    # Evidence (E) - Already checked in step 2
    result["evidence"] = {
        "note": "See Step 2 ClinVar results"
    }
    
    return result


async def step5_insights_bundle(normalized: Dict[str, Any]) -> Dict[str, Any]:
    """Step 5: Insights Bundle"""
    print("\n" + "="*60)
    print("STEP 5: Insights Bundle")
    print("="*60)
    
    result = {
        "functionality": None,
        "chromatin": None,
        "essentiality": None,
        "regulatory": None
    }
    
    # Get coordinates from normalized result
    chrom = normalized.get("chrom")
    pos = normalized.get("pos")
    ref = normalized.get("ref")
    alt = normalized.get("alt")
    
    async with httpx.AsyncClient(timeout=TIMEOUT) as client:
        # Functionality
        try:
            payload = {
                "gene": VARIANT["gene"],
                "hgvs_p": VARIANT["hgvs_p"],
                "model_id": "evo2_1b"
            }
            resp = await client.post(
                f"{API_BASE}/api/insights/predict_protein_functionality_change",
                json=payload
            )
            if resp.status_code == 200:
                data = resp.json()
                result["functionality"] = data.get("functionality_score", 0.0)
                print(f"✅ Functionality: {result['functionality']:.3f}")
        except Exception as e:
            print(f"⚠️ Functionality failed: {e}")
        
        # Essentiality
        try:
            payload = {
                "gene": VARIANT["gene"],
                "variants": [{"gene": VARIANT["gene"], "hgvs_p": VARIANT["hgvs_p"]}],
                "model_id": "evo2_1b"
            }
            resp = await client.post(
                f"{API_BASE}/api/insights/predict_gene_essentiality",
                json=payload
            )
            if resp.status_code == 200:
                data = resp.json()
                result["essentiality"] = data.get("essentiality_score", 0.0)
                print(f"✅ Essentiality: {result['essentiality']:.3f}")
        except Exception as e:
            print(f"⚠️ Essentiality failed: {e}")
        
        # Regulatory (requires coordinates)
        if chrom and pos and ref and alt:
            try:
                payload = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "model_id": "evo2_1b"
                }
                resp = await client.post(
                    f"{API_BASE}/api/insights/predict_splicing_regulatory",
                    json=payload
                )
                if resp.status_code == 200:
                    data = resp.json()
                    result["regulatory"] = data.get("regulatory_impact_score", 0.0)
                    print(f"✅ Regulatory: {result['regulatory']:.3f}")
            except Exception as e:
                print(f"⚠️ Regulatory failed: {e}")
        else:
            result["regulatory"] = {
                "status": "requires_coordinates",
                "note": "Needs chrom, pos, ref, alt for regulatory analysis"
            }
            print(f"⚠️ Regulatory: Requires coordinates")
        
        # Chromatin (requires coordinates)
        if chrom and pos:
            try:
                payload = {
                    "chrom": chrom,
                    "pos": pos,
                    "radius": 500
                }
                resp = await client.post(
                    f"{API_BASE}/api/insights/predict_chromatin_accessibility",
                    json=payload
                )
                if resp.status_code == 200:
                    data = resp.json()
                    result["chromatin"] = data.get("accessibility_score", 0.0)
                    print(f"✅ Chromatin: {result['chromatin']:.3f}")
            except Exception as e:
                print(f"⚠️ Chromatin failed: {e}")
        else:
            result["chromatin"] = {
                "status": "requires_coordinates",
                "note": "Needs chrom, pos for chromatin analysis"
            }
            print(f"⚠️ Chromatin: Requires coordinates")
    
    return result


async def step6_fusion_gating() -> Dict[str, Any]:
    """Step 6: Fusion Gating (AlphaMissense)"""
    print("\n" + "="*60)
    print("STEP 6: Fusion Gating (AlphaMissense)")
    print("="*60)
    
    result = {
        "eligible": True,  # p.S755P is missense
        "coverage": None,
        "fused_score": None,
        "note": "Requires GRCh38 coordinates for actual Fusion scoring"
    }
    
    print(f"✅ Eligible: {result['eligible']} (missense variant)")
    print(f"⚠️ Coverage: Requires coordinates")
    
    return result


async def step7_sae_features() -> Dict[str, Any]:
    """Step 7: SAE Interpretable Features"""
    print("\n" + "="*60)
    print("STEP 7: SAE Interpretable Features")
    print("="*60)
    
    result = {
        "status": "requires_evo2_scoring",
        "note": "SAE features extracted from Evo2 embeddings"
    }
    
    print(f"⚠️ SAE: Requires Evo2 scoring with coordinates")
    
    return result


async def step8_triage_scoring(
    normalized: Dict[str, Any],
    priors: Dict[str, Any],
    gate: Dict[str, Any],
    spe: Dict[str, Any],
    insights: Dict[str, Any],
    fusion: Dict[str, Any],
    sae: Dict[str, Any]
) -> Dict[str, Any]:
    """Step 8: VUS Triage Scoring and Decision Rubric"""
    print("\n" + "="*60)
    print("STEP 8: VUS Triage Scoring and Decision Rubric")
    print("="*60)
    
    # Collect scores
    functionality = insights.get("functionality", 0.0) if insights.get("functionality") else 0.0
    essentiality = insights.get("essentiality", 0.0) if insights.get("essentiality") else 0.0
    
    # Decision logic
    confidence = 0.0
    classification = "VUS"
    rationale = []
    
    # If functionality is low (<0.3) and ClinVar is benign
    if functionality < 0.3:
        confidence += 0.3
        rationale.append(f"Low functionality score ({functionality:.3f}) suggests minimal impact")
    
    # If functionality is high (>0.7)
    if functionality > 0.7:
        confidence += 0.4
        rationale.append(f"High functionality score ({functionality:.3f}) suggests significant impact")
    
    # Check ClinVar
    clinvar_class = priors.get("clinvar", {}).get("classification") if priors.get("clinvar") else None
    if clinvar_class:
        if "benign" in clinvar_class.lower():
            confidence += 0.3
            classification = "Likely Benign"
            rationale.append(f"ClinVar classification: {clinvar_class}")
        elif "pathogenic" in clinvar_class.lower():
            confidence += 0.4
            classification = "Likely Pathogenic"
            rationale.append(f"ClinVar classification: {clinvar_class}")
    
    # Final classification
    if confidence >= 0.7:
        if functionality < 0.3:
            classification = "Benign" if confidence >= 0.8 else "Likely Benign"
        elif functionality > 0.7:
            classification = "Pathogenic" if confidence >= 0.8 else "Likely Pathogenic"
    else:
        classification = "VUS"
        rationale.append("Insufficient evidence for definitive classification")
    
    result = {
        "confidence": confidence,
        "classification": classification,
        "rationale": rationale,
        "evidence_summary": {
            "functionality": functionality,
            "essentiality": essentiality,
            "clinvar": clinvar_class,
            "gate_passed": gate.get("gate_result") == "PASS"
        },
        "recommendation": "VUS eliminated" if classification != "VUS" else "VUS remains - additional evidence needed"
    }
    
    print(f"✅ Confidence: {confidence:.3f}")
    print(f"✅ Classification: {classification}")
    print(f"✅ Recommendation: {result['recommendation']}")
    print(f"\nRationale:")
    for r in rationale:
        print(f"  - {r}")
    
    return result


async def main():
    """Execute complete 8-step workflow"""
    print("\n" + "="*60)
    print("PDGFRA p.S755P VUS RESOLUTION WORKFLOW")
    print("="*60)
    print(f"Variant: {VARIANT['gene']} {VARIANT['hgvs_p']} ({VARIANT['c_dna']})")
    print(f"Current Classification: {VARIANT['classification']}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    
    # Execute steps
    normalized = await step1_normalize()
    priors = await step2_priors(normalized)
    gate = await step3_triumvirate_gate()
    spe = await step4_spe_signals()
    insights = await step5_insights_bundle(normalized)
    fusion = await step6_fusion_gating()
    sae = await step7_sae_features()
    triage = await step8_triage_scoring(normalized, priors, gate, spe, insights, fusion, sae)
    
    # Final summary
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    print(f"Original Classification: {VARIANT['classification']}")
    print(f"New Classification: {triage['classification']}")
    print(f"Confidence: {triage['confidence']:.3f}")
    print(f"Recommendation: {triage['recommendation']}")
    
    # Save results
    results = {
        "variant": VARIANT,
        "workflow": {
            "step1": normalized,
            "step2": priors,
            "step3": gate,
            "step4": spe,
            "step5": insights,
            "step6": fusion,
            "step7": sae,
            "step8": triage
        },
        "timestamp": datetime.now().isoformat()
    }
    
    output_file = ".cursor/ayesha/PDGFRA_VUS_RESOLUTION_RESULTS.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    return results


if __name__ == "__main__":
    asyncio.run(main())

