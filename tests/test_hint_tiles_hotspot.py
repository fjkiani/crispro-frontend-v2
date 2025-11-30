"""
Test Hint Tiles Hotspot Detection Integration (P1.1)
====================================================
Tests that hotspot mutations trigger MEK/RAF hint tiles.

Owner: Zo
Date: January 13, 2025
Manager: SR
"""

import sys
sys.path.insert(0, 'oncology-coPilot/oncology-backend-minimal')

from api.services.hint_tiles_service import get_hint_tiles


def test_hotspot_kras_g12d_hint_tile():
    """Test that KRAS G12D hotspot triggers MEK/RAF hint tile."""
    
    # SAE features with KRAS G12D hotspot
    sae_features = {
        "hotspot_mutation": True,
        "hotspot_details": {
            "gene": "KRAS",
            "mutation": "G12D",
            "pathway": "MAPK",
            "cosmic_id": "COSV55455408"
        },
        "dna_repair_capacity": 0.50,
        "pathway_burden_mapk": 0.72
    }
    
    # Tumor context (post-NGS)
    tumor_context = {
        "hrd_score": 35,
        "somatic_mutations": [
            {"gene": "KRAS", "hgvs_p": "p.G12D"}
        ]
    }
    
    # Generate hint tiles
    result = get_hint_tiles(
        germline_status="negative",
        tumor_context=tumor_context,
        ca125_intelligence={"current_value": 2842, "burden_classification": "EXTENSIVE"},
        next_test_recommendations=[],  # No next tests needed (post-NGS)
        treatment_history=[],
        trials_matched=5,
        sae_features=sae_features
    )
    
    # Verify hotspot hint tile is present
    hint_tiles = result["hint_tiles"]
    
    # Find hotspot tile
    hotspot_tile = None
    for tile in hint_tiles:
        if "Hotspot" in tile.get("title", ""):
            hotspot_tile = tile
            break
    
    assert hotspot_tile is not None, "❌ Hotspot hint tile not generated"
    assert "KRAS" in hotspot_tile["message"], "❌ KRAS not mentioned in message"
    assert "G12D" in hotspot_tile["message"], "❌ G12D not mentioned in message"
    assert "MEK/RAF" in hotspot_tile["message"], "❌ MEK/RAF not mentioned in message"
    assert tile["category"] == "trials_lever", "❌ Wrong category (should be trials_lever)"
    assert any("COSMIC" in reason for reason in hotspot_tile["reasons"]), "❌ COSMIC not mentioned in reasons"
    assert any("RUO" in reason for reason in hotspot_tile["reasons"]), "❌ RUO label missing"
    
    print("✅ KRAS G12D hotspot hint tile: PASSED")
    print(f"   Title: {hotspot_tile['title']}")
    print(f"   Message: {hotspot_tile['message']}")
    print(f"   Reasons: {len(hotspot_tile['reasons'])} reasons")


def test_no_hotspot_no_hint_tile():
    """Test that no hotspot mutation = no hotspot hint tile."""
    
    # SAE features WITHOUT hotspot
    sae_features = {
        "hotspot_mutation": False,
        "hotspot_details": None,
        "dna_repair_capacity": 0.82
    }
    
    # Tumor context (post-NGS, no hotspot)
    tumor_context = {
        "hrd_score": 58,
        "somatic_mutations": [
            {"gene": "BRCA1", "hgvs_p": "p.C61G"}  # Not a hotspot
        ]
    }
    
    # Generate hint tiles
    result = get_hint_tiles(
        germline_status="negative",
        tumor_context=tumor_context,
        ca125_intelligence={"current_value": 2842, "burden_classification": "EXTENSIVE"},
        next_test_recommendations=[],
        treatment_history=[],
        trials_matched=5,
        sae_features=sae_features
    )
    
    # Verify NO hotspot hint tile
    hint_tiles = result["hint_tiles"]
    hotspot_tiles = [t for t in hint_tiles if "Hotspot" in t.get("title", "")]
    
    assert len(hotspot_tiles) == 0, f"❌ Hotspot tile generated when no hotspot present: {hotspot_tiles}"
    
    print("✅ No hotspot = no hint tile: PASSED")


def test_braf_v600e_hotspot_hint_tile():
    """Test that BRAF V600E hotspot triggers MEK/RAF hint tile."""
    
    sae_features = {
        "hotspot_mutation": True,
        "hotspot_details": {
            "gene": "BRAF",
            "mutation": "V600E",
            "pathway": "MAPK",
            "cosmic_id": "COSV55455406"
        }
    }
    
    result = get_hint_tiles(
        germline_status="negative",
        tumor_context={"hrd_score": 40},
        ca125_intelligence={"current_value": 1500, "burden_classification": "SIGNIFICANT"},
        next_test_recommendations=[],
        treatment_history=[],
        trials_matched=3,
        sae_features=sae_features
    )
    
    hotspot_tile = next((t for t in result["hint_tiles"] if "Hotspot" in t.get("title", "")), None)
    
    assert hotspot_tile is not None, "❌ BRAF hotspot hint tile not generated"
    assert "BRAF" in hotspot_tile["message"], "❌ BRAF not mentioned"
    assert "V600E" in hotspot_tile["message"], "❌ V600E not mentioned"
    
    print("✅ BRAF V600E hotspot hint tile: PASSED")


if __name__ == "__main__":
    print("\n⚔️ P1.1: HOTSPOT HINT TILES TEST SUITE")
    print("=" * 80)
    
    test_hotspot_kras_g12d_hint_tile()
    test_no_hotspot_no_hint_tile()
    test_braf_v600e_hotspot_hint_tile()
    
    print("\n✅ ALL HOTSPOT HINT TILE TESTS PASSED!")
    print("=" * 80)







