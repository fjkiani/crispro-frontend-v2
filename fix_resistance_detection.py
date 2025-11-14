#!/usr/bin/env python3
"""
Fix: Add resistance detection logic to guidance system
Based on user requirements for MM drug response prediction
"""

# Resistance markers that need to be implemented:

RESISTANCE_MARKERS = {
    # Proteasome inhibitor resistance
    "proteasome_inhibitor": {
        "PSMB5": {
            # Known bortezomib resistance mutations
            "resistance_variants": ["A49T", "T21A", "S27P", "G223A"],
            "confidence_penalty": -0.3,
            "efficacy_penalty": -0.2,
            "flag": "PSMB5_resistance"
        }
    },
    
    # IMiD resistance  
    "imid": {
        "CRBN": {
            "resistance_variants": ["loss_of_function", "truncating", "splice"],
            "confidence_penalty": -0.4,
            "efficacy_penalty": -0.3,
            "flag": "CRBN_resistance"
        }
    },
    
    # General chemo tolerance
    "chemotherapy": {
        "TP53": {
            # High-risk TP53 mutations confer chemo tolerance
            "resistance_context": ["R248W", "R273H", "R175H", "G245S", "R282W"],
            "confidence_penalty": -0.1,  # Moderate penalty, not absolute
            "efficacy_penalty": -0.05,
            "flag": "TP53_high_risk"
        }
    },
    
    # Platinum resistance
    "platinum": {
        "BRCA1": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"],
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        },
        "BRCA2": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"], 
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        }
    }
}

SENSITIVITY_MARKERS = {
    # HRD positive for platinum sensitivity
    "platinum": {
        "hrd_markers": ["LOH_high", "TAI_high", "LST_high"],
        "confidence_boost": 0.2,
        "efficacy_boost": 0.1,
        "required_genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK2"]
    },
    
    # MAPK dependence for BRAF/MEK
    "mapk_inhibitor": {
        "sensitivity_genes": ["KRAS", "NRAS", "BRAF"],
        "hotspot_variants": ["G12D", "G12V", "G13D", "Q61R", "Q61K", "V600E"],
        "confidence_boost": 0.15,
        "efficacy_boost": 0.1,
        "pathway_requirement": "ras_mapk"
    }
}

def detect_resistance(gene: str, hgvs_p: str, drug_class: str) -> dict:
    """
    Detect resistance markers for a given mutation and drug class
    Returns: {
        "resistance_detected": bool,
        "resistance_type": str,
        "confidence_penalty": float,
        "efficacy_penalty": float,
        "rationale": str
    }
    """
    
    # Extract variant from hgvs_p (e.g., "p.A49T" -> "A49T")
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    
    drug_key = drug_class.lower().replace(" ", "_")
    
    # Check for direct resistance markers
    if drug_key in RESISTANCE_MARKERS:
        drug_resistances = RESISTANCE_MARKERS[drug_key]
        
        if gene in drug_resistances:
            gene_resistance = drug_resistances[gene]
            
            if variant in gene_resistance.get("resistance_variants", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} confers {drug_class} resistance"
                }
    
    # Check for general chemo tolerance (TP53)
    if "chemotherapy" in RESISTANCE_MARKERS:
        chemo_resistances = RESISTANCE_MARKERS["chemotherapy"]
        
        if gene in chemo_resistances:
            gene_resistance = chemo_resistances[gene]
            
            if variant in gene_resistance.get("resistance_context", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} high-risk context for chemotherapy"
                }
    
    return {
        "resistance_detected": False,
        "resistance_type": None,
        "confidence_penalty": 0.0,
        "efficacy_penalty": 0.0,
        "rationale": ""
    }

def detect_sensitivity(gene: str, hgvs_p: str, drug_class: str, fused_s: float = 0.0) -> dict:
    """
    Detect sensitivity markers for a given mutation and drug class
    """
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    drug_key = drug_class.lower().replace(" ", "_")
    
    # MAPK pathway sensitivity
    if "mapk" in drug_key or "braf" in drug_key or "mek" in drug_key:
        if gene in ["KRAS", "NRAS", "BRAF"]:
            mapk_markers = SENSITIVITY_MARKERS["mapk_inhibitor"]
            if variant in mapk_markers["hotspot_variants"] and fused_s >= 0.9:
                return {
                    "sensitivity_detected": True,
                    "sensitivity_type": "MAPK_hotspot",
                    "confidence_boost": mapk_markers["confidence_boost"],
                    "efficacy_boost": mapk_markers["efficacy_boost"],
                    "rationale": f"{gene} {variant} hotspot with high fused S (≥0.9)"
                }
    
    return {
        "sensitivity_detected": False,
        "sensitivity_type": None,
        "confidence_boost": 0.0,
        "efficacy_boost": 0.0,
        "rationale": ""
    }

# Test cases to validate
test_cases = [
    {"gene": "PSMB5", "hgvs_p": "p.A49T", "drug": "proteasome inhibitor", "expected": "resistance"},
    {"gene": "CRBN", "hgvs_p": "p.Q99*", "drug": "IMiD", "expected": "resistance"},
    {"gene": "TP53", "hgvs_p": "p.R248W", "drug": "chemotherapy", "expected": "tolerance"},
    {"gene": "KRAS", "hgvs_p": "p.G12D", "drug": "MEK inhibitor", "expected": "sensitivity"},
    {"gene": "BRAF", "hgvs_p": "p.V600E", "drug": "BRAF inhibitor", "expected": "sensitivity"},
]

if __name__ == "__main__":
    print("🧬 Testing Resistance/Sensitivity Detection")
    print("=" * 50)
    
    for case in test_cases:
        print(f"\n🔬 Testing: {case['gene']} {case['hgvs_p']} → {case['drug']}")
        
        resistance = detect_resistance(case['gene'], case['hgvs_p'], case['drug'])
        sensitivity = detect_sensitivity(case['gene'], case['hgvs_p'], case['drug'], fused_s=0.95)
        
        if resistance['resistance_detected']:
            print(f"  ❌ RESISTANCE: {resistance['resistance_type']}")
            print(f"     Confidence penalty: {resistance['confidence_penalty']}")
            print(f"     Rationale: {resistance['rationale']}")
        elif sensitivity['sensitivity_detected']:
            print(f"  ✅ SENSITIVITY: {sensitivity['sensitivity_type']}")
            print(f"     Confidence boost: {sensitivity['confidence_boost']}")
            print(f"     Rationale: {sensitivity['rationale']}")
        else:
            print(f"  ⚪ No specific resistance/sensitivity markers detected")

Fix: Add resistance detection logic to guidance system
Based on user requirements for MM drug response prediction
"""

# Resistance markers that need to be implemented:

RESISTANCE_MARKERS = {
    # Proteasome inhibitor resistance
    "proteasome_inhibitor": {
        "PSMB5": {
            # Known bortezomib resistance mutations
            "resistance_variants": ["A49T", "T21A", "S27P", "G223A"],
            "confidence_penalty": -0.3,
            "efficacy_penalty": -0.2,
            "flag": "PSMB5_resistance"
        }
    },
    
    # IMiD resistance  
    "imid": {
        "CRBN": {
            "resistance_variants": ["loss_of_function", "truncating", "splice"],
            "confidence_penalty": -0.4,
            "efficacy_penalty": -0.3,
            "flag": "CRBN_resistance"
        }
    },
    
    # General chemo tolerance
    "chemotherapy": {
        "TP53": {
            # High-risk TP53 mutations confer chemo tolerance
            "resistance_context": ["R248W", "R273H", "R175H", "G245S", "R282W"],
            "confidence_penalty": -0.1,  # Moderate penalty, not absolute
            "efficacy_penalty": -0.05,
            "flag": "TP53_high_risk"
        }
    },
    
    # Platinum resistance
    "platinum": {
        "BRCA1": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"],
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        },
        "BRCA2": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"], 
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        }
    }
}

SENSITIVITY_MARKERS = {
    # HRD positive for platinum sensitivity
    "platinum": {
        "hrd_markers": ["LOH_high", "TAI_high", "LST_high"],
        "confidence_boost": 0.2,
        "efficacy_boost": 0.1,
        "required_genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK2"]
    },
    
    # MAPK dependence for BRAF/MEK
    "mapk_inhibitor": {
        "sensitivity_genes": ["KRAS", "NRAS", "BRAF"],
        "hotspot_variants": ["G12D", "G12V", "G13D", "Q61R", "Q61K", "V600E"],
        "confidence_boost": 0.15,
        "efficacy_boost": 0.1,
        "pathway_requirement": "ras_mapk"
    }
}

def detect_resistance(gene: str, hgvs_p: str, drug_class: str) -> dict:
    """
    Detect resistance markers for a given mutation and drug class
    Returns: {
        "resistance_detected": bool,
        "resistance_type": str,
        "confidence_penalty": float,
        "efficacy_penalty": float,
        "rationale": str
    }
    """
    
    # Extract variant from hgvs_p (e.g., "p.A49T" -> "A49T")
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    
    drug_key = drug_class.lower().replace(" ", "_")
    
    # Check for direct resistance markers
    if drug_key in RESISTANCE_MARKERS:
        drug_resistances = RESISTANCE_MARKERS[drug_key]
        
        if gene in drug_resistances:
            gene_resistance = drug_resistances[gene]
            
            if variant in gene_resistance.get("resistance_variants", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} confers {drug_class} resistance"
                }
    
    # Check for general chemo tolerance (TP53)
    if "chemotherapy" in RESISTANCE_MARKERS:
        chemo_resistances = RESISTANCE_MARKERS["chemotherapy"]
        
        if gene in chemo_resistances:
            gene_resistance = chemo_resistances[gene]
            
            if variant in gene_resistance.get("resistance_context", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} high-risk context for chemotherapy"
                }
    
    return {
        "resistance_detected": False,
        "resistance_type": None,
        "confidence_penalty": 0.0,
        "efficacy_penalty": 0.0,
        "rationale": ""
    }

def detect_sensitivity(gene: str, hgvs_p: str, drug_class: str, fused_s: float = 0.0) -> dict:
    """
    Detect sensitivity markers for a given mutation and drug class
    """
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    drug_key = drug_class.lower().replace(" ", "_")
    
    # MAPK pathway sensitivity
    if "mapk" in drug_key or "braf" in drug_key or "mek" in drug_key:
        if gene in ["KRAS", "NRAS", "BRAF"]:
            mapk_markers = SENSITIVITY_MARKERS["mapk_inhibitor"]
            if variant in mapk_markers["hotspot_variants"] and fused_s >= 0.9:
                return {
                    "sensitivity_detected": True,
                    "sensitivity_type": "MAPK_hotspot",
                    "confidence_boost": mapk_markers["confidence_boost"],
                    "efficacy_boost": mapk_markers["efficacy_boost"],
                    "rationale": f"{gene} {variant} hotspot with high fused S (≥0.9)"
                }
    
    return {
        "sensitivity_detected": False,
        "sensitivity_type": None,
        "confidence_boost": 0.0,
        "efficacy_boost": 0.0,
        "rationale": ""
    }

# Test cases to validate
test_cases = [
    {"gene": "PSMB5", "hgvs_p": "p.A49T", "drug": "proteasome inhibitor", "expected": "resistance"},
    {"gene": "CRBN", "hgvs_p": "p.Q99*", "drug": "IMiD", "expected": "resistance"},
    {"gene": "TP53", "hgvs_p": "p.R248W", "drug": "chemotherapy", "expected": "tolerance"},
    {"gene": "KRAS", "hgvs_p": "p.G12D", "drug": "MEK inhibitor", "expected": "sensitivity"},
    {"gene": "BRAF", "hgvs_p": "p.V600E", "drug": "BRAF inhibitor", "expected": "sensitivity"},
]

if __name__ == "__main__":
    print("🧬 Testing Resistance/Sensitivity Detection")
    print("=" * 50)
    
    for case in test_cases:
        print(f"\n🔬 Testing: {case['gene']} {case['hgvs_p']} → {case['drug']}")
        
        resistance = detect_resistance(case['gene'], case['hgvs_p'], case['drug'])
        sensitivity = detect_sensitivity(case['gene'], case['hgvs_p'], case['drug'], fused_s=0.95)
        
        if resistance['resistance_detected']:
            print(f"  ❌ RESISTANCE: {resistance['resistance_type']}")
            print(f"     Confidence penalty: {resistance['confidence_penalty']}")
            print(f"     Rationale: {resistance['rationale']}")
        elif sensitivity['sensitivity_detected']:
            print(f"  ✅ SENSITIVITY: {sensitivity['sensitivity_type']}")
            print(f"     Confidence boost: {sensitivity['confidence_boost']}")
            print(f"     Rationale: {sensitivity['rationale']}")
        else:
            print(f"  ⚪ No specific resistance/sensitivity markers detected")







Fix: Add resistance detection logic to guidance system
Based on user requirements for MM drug response prediction
"""

# Resistance markers that need to be implemented:

RESISTANCE_MARKERS = {
    # Proteasome inhibitor resistance
    "proteasome_inhibitor": {
        "PSMB5": {
            # Known bortezomib resistance mutations
            "resistance_variants": ["A49T", "T21A", "S27P", "G223A"],
            "confidence_penalty": -0.3,
            "efficacy_penalty": -0.2,
            "flag": "PSMB5_resistance"
        }
    },
    
    # IMiD resistance  
    "imid": {
        "CRBN": {
            "resistance_variants": ["loss_of_function", "truncating", "splice"],
            "confidence_penalty": -0.4,
            "efficacy_penalty": -0.3,
            "flag": "CRBN_resistance"
        }
    },
    
    # General chemo tolerance
    "chemotherapy": {
        "TP53": {
            # High-risk TP53 mutations confer chemo tolerance
            "resistance_context": ["R248W", "R273H", "R175H", "G245S", "R282W"],
            "confidence_penalty": -0.1,  # Moderate penalty, not absolute
            "efficacy_penalty": -0.05,
            "flag": "TP53_high_risk"
        }
    },
    
    # Platinum resistance
    "platinum": {
        "BRCA1": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"],
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        },
        "BRCA2": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"], 
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        }
    }
}

SENSITIVITY_MARKERS = {
    # HRD positive for platinum sensitivity
    "platinum": {
        "hrd_markers": ["LOH_high", "TAI_high", "LST_high"],
        "confidence_boost": 0.2,
        "efficacy_boost": 0.1,
        "required_genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK2"]
    },
    
    # MAPK dependence for BRAF/MEK
    "mapk_inhibitor": {
        "sensitivity_genes": ["KRAS", "NRAS", "BRAF"],
        "hotspot_variants": ["G12D", "G12V", "G13D", "Q61R", "Q61K", "V600E"],
        "confidence_boost": 0.15,
        "efficacy_boost": 0.1,
        "pathway_requirement": "ras_mapk"
    }
}

def detect_resistance(gene: str, hgvs_p: str, drug_class: str) -> dict:
    """
    Detect resistance markers for a given mutation and drug class
    Returns: {
        "resistance_detected": bool,
        "resistance_type": str,
        "confidence_penalty": float,
        "efficacy_penalty": float,
        "rationale": str
    }
    """
    
    # Extract variant from hgvs_p (e.g., "p.A49T" -> "A49T")
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    
    drug_key = drug_class.lower().replace(" ", "_")
    
    # Check for direct resistance markers
    if drug_key in RESISTANCE_MARKERS:
        drug_resistances = RESISTANCE_MARKERS[drug_key]
        
        if gene in drug_resistances:
            gene_resistance = drug_resistances[gene]
            
            if variant in gene_resistance.get("resistance_variants", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} confers {drug_class} resistance"
                }
    
    # Check for general chemo tolerance (TP53)
    if "chemotherapy" in RESISTANCE_MARKERS:
        chemo_resistances = RESISTANCE_MARKERS["chemotherapy"]
        
        if gene in chemo_resistances:
            gene_resistance = chemo_resistances[gene]
            
            if variant in gene_resistance.get("resistance_context", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} high-risk context for chemotherapy"
                }
    
    return {
        "resistance_detected": False,
        "resistance_type": None,
        "confidence_penalty": 0.0,
        "efficacy_penalty": 0.0,
        "rationale": ""
    }

def detect_sensitivity(gene: str, hgvs_p: str, drug_class: str, fused_s: float = 0.0) -> dict:
    """
    Detect sensitivity markers for a given mutation and drug class
    """
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    drug_key = drug_class.lower().replace(" ", "_")
    
    # MAPK pathway sensitivity
    if "mapk" in drug_key or "braf" in drug_key or "mek" in drug_key:
        if gene in ["KRAS", "NRAS", "BRAF"]:
            mapk_markers = SENSITIVITY_MARKERS["mapk_inhibitor"]
            if variant in mapk_markers["hotspot_variants"] and fused_s >= 0.9:
                return {
                    "sensitivity_detected": True,
                    "sensitivity_type": "MAPK_hotspot",
                    "confidence_boost": mapk_markers["confidence_boost"],
                    "efficacy_boost": mapk_markers["efficacy_boost"],
                    "rationale": f"{gene} {variant} hotspot with high fused S (≥0.9)"
                }
    
    return {
        "sensitivity_detected": False,
        "sensitivity_type": None,
        "confidence_boost": 0.0,
        "efficacy_boost": 0.0,
        "rationale": ""
    }

# Test cases to validate
test_cases = [
    {"gene": "PSMB5", "hgvs_p": "p.A49T", "drug": "proteasome inhibitor", "expected": "resistance"},
    {"gene": "CRBN", "hgvs_p": "p.Q99*", "drug": "IMiD", "expected": "resistance"},
    {"gene": "TP53", "hgvs_p": "p.R248W", "drug": "chemotherapy", "expected": "tolerance"},
    {"gene": "KRAS", "hgvs_p": "p.G12D", "drug": "MEK inhibitor", "expected": "sensitivity"},
    {"gene": "BRAF", "hgvs_p": "p.V600E", "drug": "BRAF inhibitor", "expected": "sensitivity"},
]

if __name__ == "__main__":
    print("🧬 Testing Resistance/Sensitivity Detection")
    print("=" * 50)
    
    for case in test_cases:
        print(f"\n🔬 Testing: {case['gene']} {case['hgvs_p']} → {case['drug']}")
        
        resistance = detect_resistance(case['gene'], case['hgvs_p'], case['drug'])
        sensitivity = detect_sensitivity(case['gene'], case['hgvs_p'], case['drug'], fused_s=0.95)
        
        if resistance['resistance_detected']:
            print(f"  ❌ RESISTANCE: {resistance['resistance_type']}")
            print(f"     Confidence penalty: {resistance['confidence_penalty']}")
            print(f"     Rationale: {resistance['rationale']}")
        elif sensitivity['sensitivity_detected']:
            print(f"  ✅ SENSITIVITY: {sensitivity['sensitivity_type']}")
            print(f"     Confidence boost: {sensitivity['confidence_boost']}")
            print(f"     Rationale: {sensitivity['rationale']}")
        else:
            print(f"  ⚪ No specific resistance/sensitivity markers detected")

Fix: Add resistance detection logic to guidance system
Based on user requirements for MM drug response prediction
"""

# Resistance markers that need to be implemented:

RESISTANCE_MARKERS = {
    # Proteasome inhibitor resistance
    "proteasome_inhibitor": {
        "PSMB5": {
            # Known bortezomib resistance mutations
            "resistance_variants": ["A49T", "T21A", "S27P", "G223A"],
            "confidence_penalty": -0.3,
            "efficacy_penalty": -0.2,
            "flag": "PSMB5_resistance"
        }
    },
    
    # IMiD resistance  
    "imid": {
        "CRBN": {
            "resistance_variants": ["loss_of_function", "truncating", "splice"],
            "confidence_penalty": -0.4,
            "efficacy_penalty": -0.3,
            "flag": "CRBN_resistance"
        }
    },
    
    # General chemo tolerance
    "chemotherapy": {
        "TP53": {
            # High-risk TP53 mutations confer chemo tolerance
            "resistance_context": ["R248W", "R273H", "R175H", "G245S", "R282W"],
            "confidence_penalty": -0.1,  # Moderate penalty, not absolute
            "efficacy_penalty": -0.05,
            "flag": "TP53_high_risk"
        }
    },
    
    # Platinum resistance
    "platinum": {
        "BRCA1": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"],
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        },
        "BRCA2": {
            "sensitivity_variants": ["truncating", "frameshift", "nonsense"], 
            "resistance_variants": ["reversion", "secondary"],
            "hrd_context_required": True
        }
    }
}

SENSITIVITY_MARKERS = {
    # HRD positive for platinum sensitivity
    "platinum": {
        "hrd_markers": ["LOH_high", "TAI_high", "LST_high"],
        "confidence_boost": 0.2,
        "efficacy_boost": 0.1,
        "required_genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK2"]
    },
    
    # MAPK dependence for BRAF/MEK
    "mapk_inhibitor": {
        "sensitivity_genes": ["KRAS", "NRAS", "BRAF"],
        "hotspot_variants": ["G12D", "G12V", "G13D", "Q61R", "Q61K", "V600E"],
        "confidence_boost": 0.15,
        "efficacy_boost": 0.1,
        "pathway_requirement": "ras_mapk"
    }
}

def detect_resistance(gene: str, hgvs_p: str, drug_class: str) -> dict:
    """
    Detect resistance markers for a given mutation and drug class
    Returns: {
        "resistance_detected": bool,
        "resistance_type": str,
        "confidence_penalty": float,
        "efficacy_penalty": float,
        "rationale": str
    }
    """
    
    # Extract variant from hgvs_p (e.g., "p.A49T" -> "A49T")
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    
    drug_key = drug_class.lower().replace(" ", "_")
    
    # Check for direct resistance markers
    if drug_key in RESISTANCE_MARKERS:
        drug_resistances = RESISTANCE_MARKERS[drug_key]
        
        if gene in drug_resistances:
            gene_resistance = drug_resistances[gene]
            
            if variant in gene_resistance.get("resistance_variants", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} confers {drug_class} resistance"
                }
    
    # Check for general chemo tolerance (TP53)
    if "chemotherapy" in RESISTANCE_MARKERS:
        chemo_resistances = RESISTANCE_MARKERS["chemotherapy"]
        
        if gene in chemo_resistances:
            gene_resistance = chemo_resistances[gene]
            
            if variant in gene_resistance.get("resistance_context", []):
                return {
                    "resistance_detected": True,
                    "resistance_type": gene_resistance["flag"],
                    "confidence_penalty": gene_resistance["confidence_penalty"],
                    "efficacy_penalty": gene_resistance["efficacy_penalty"],
                    "rationale": f"{gene} {variant} high-risk context for chemotherapy"
                }
    
    return {
        "resistance_detected": False,
        "resistance_type": None,
        "confidence_penalty": 0.0,
        "efficacy_penalty": 0.0,
        "rationale": ""
    }

def detect_sensitivity(gene: str, hgvs_p: str, drug_class: str, fused_s: float = 0.0) -> dict:
    """
    Detect sensitivity markers for a given mutation and drug class
    """
    variant = hgvs_p.split(".")[-1] if "." in hgvs_p else hgvs_p
    drug_key = drug_class.lower().replace(" ", "_")
    
    # MAPK pathway sensitivity
    if "mapk" in drug_key or "braf" in drug_key or "mek" in drug_key:
        if gene in ["KRAS", "NRAS", "BRAF"]:
            mapk_markers = SENSITIVITY_MARKERS["mapk_inhibitor"]
            if variant in mapk_markers["hotspot_variants"] and fused_s >= 0.9:
                return {
                    "sensitivity_detected": True,
                    "sensitivity_type": "MAPK_hotspot",
                    "confidence_boost": mapk_markers["confidence_boost"],
                    "efficacy_boost": mapk_markers["efficacy_boost"],
                    "rationale": f"{gene} {variant} hotspot with high fused S (≥0.9)"
                }
    
    return {
        "sensitivity_detected": False,
        "sensitivity_type": None,
        "confidence_boost": 0.0,
        "efficacy_boost": 0.0,
        "rationale": ""
    }

# Test cases to validate
test_cases = [
    {"gene": "PSMB5", "hgvs_p": "p.A49T", "drug": "proteasome inhibitor", "expected": "resistance"},
    {"gene": "CRBN", "hgvs_p": "p.Q99*", "drug": "IMiD", "expected": "resistance"},
    {"gene": "TP53", "hgvs_p": "p.R248W", "drug": "chemotherapy", "expected": "tolerance"},
    {"gene": "KRAS", "hgvs_p": "p.G12D", "drug": "MEK inhibitor", "expected": "sensitivity"},
    {"gene": "BRAF", "hgvs_p": "p.V600E", "drug": "BRAF inhibitor", "expected": "sensitivity"},
]

if __name__ == "__main__":
    print("🧬 Testing Resistance/Sensitivity Detection")
    print("=" * 50)
    
    for case in test_cases:
        print(f"\n🔬 Testing: {case['gene']} {case['hgvs_p']} → {case['drug']}")
        
        resistance = detect_resistance(case['gene'], case['hgvs_p'], case['drug'])
        sensitivity = detect_sensitivity(case['gene'], case['hgvs_p'], case['drug'], fused_s=0.95)
        
        if resistance['resistance_detected']:
            print(f"  ❌ RESISTANCE: {resistance['resistance_type']}")
            print(f"     Confidence penalty: {resistance['confidence_penalty']}")
            print(f"     Rationale: {resistance['rationale']}")
        elif sensitivity['sensitivity_detected']:
            print(f"  ✅ SENSITIVITY: {sensitivity['sensitivity_type']}")
            print(f"     Confidence boost: {sensitivity['confidence_boost']}")
            print(f"     Rationale: {sensitivity['rationale']}")
        else:
            print(f"  ⚪ No specific resistance/sensitivity markers detected")






