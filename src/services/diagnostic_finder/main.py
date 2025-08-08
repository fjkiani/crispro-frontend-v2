from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import logging
import random

# --- Constants based on Jiao et al., 2024, Nature Communications ---
# DOI: 10.1038/s41467-024-50243-x

# Scaffold for the BhCas12b reprogrammed tracrRNA (Rptr).
# The anti-repeat sequence will be inserted into this scaffold.
# This structure is derived from the paper's description and supplementary figures.
# The 'X's represent the long anti-repeat, and 'Y's represent the short anti-repeat.
BHCAS12B_TRACR_SCAFFOLD = {
    "long_anti_repeat": "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",  # 31 nt
    "scaffold_part_1": "GUGAUAGUGU",
    "short_anti_repeat": "YYYYY", # 5 nt
    "scaffold_part_2": "UAAUUCUCUACUUUUGU",
}

# The PAM sequence for the dsDNA activator for BhCas12b
ACTIVATOR_PAM = "ATTG" # Using 'G' for 'N' for simplicity

# The fluorescent reporter
SSDNA_REPORTER = "5'-FAM-TTATT-3'-IABkFQ"

# Optimal Truncation for activator based on the paper's findings (Fig. 4)
NTS_TRUNCATION = 6
TS_TRUNCATION = 2

# --- Utility Functions ---

def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, base) for base in reversed(seq))

# --- FastAPI App ---

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Diagnostic Finder Service",
    description="Designs CRISPR-Cas12-based diagnostic assays using the PUMA system.",
    version="1.1.0",
)

class AssayRequest(BaseModel):
    target_rna_sequence: str
    target_gene_name: str

class AssayResponse(BaseModel):
    reprogrammed_tracrRNA: str
    dsDNA_activator_NTS: str
    dsDNA_activator_TS: str
    ssDNA_reporter: str
    enzyme: str
    protocol: str
    notes: str

@app.post("/design_diagnostic_assay", response_model=AssayResponse)
async def design_diagnostic_assay(request: AssayRequest):
    """
    Designs a complete diagnostic assay package based on the PUMA system
    (Programmable tracrRNAs Unlock Motif-independent detection of ribonucleic Acids).
    
    This endpoint takes a target RNA sequence and returns all the necessary
    components to run a fluorescent-based detection assay using BhCas12b.
    """
    logger.info(f"Received assay design request for gene: {request.target_gene_name}")
    
    # --- Input Validation ---
    # The Rptr needs to bind a 31nt long anti-repeat and have a guide sequence.
    # We'll use a 20nt guide for this implementation. Total required length: 51nt.
    required_len = len(BHCAS12B_TRACR_SCAFFOLD["long_anti_repeat"]) + 20
    if len(request.target_rna_sequence) < required_len:
        raise HTTPException(
            status_code=400,
            detail=f"Target RNA sequence must be at least {required_len} nucleotides long."
        )

    # --- Design Logic based on Jiao et al., 2024 ---

    # 1. Split the target RNA sequence into the guide and anti-repeat regions.
    # The paper is complex; for a robust first version, we'll make some simplifying assumptions.
    # We will use a 20-nt guide sequence. The preceding 36 nt of the target RNA
    # will form the long (31nt) and short (5nt) anti-repeats.
    
    guide_len = 20
    long_ar_len = len(BHCAS12B_TRACR_SCAFFOLD["long_anti_repeat"])
    short_ar_len = len(BHCAS12B_TRACR_SCAFFOLD["short_anti_repeat"])
    
    # The guide is at the 3' end of the target region
    guide_sequence_rna = request.target_rna_sequence[-guide_len:]
    guide_sequence_dna = guide_sequence_rna.replace('U', 'T') # Convert to DNA for the activator

    # The anti-repeats are upstream of the guide
    anti_repeat_rna_region = request.target_rna_sequence[-(guide_len + long_ar_len + short_ar_len):-guide_len]
    
    long_anti_repeat_rna = anti_repeat_rna_region[:long_ar_len]
    short_anti_repeat_rna = anti_repeat_rna_region[long_ar_len:]

    # 2. Construct the reprogrammed tracrRNA (Rptr)
    # We need the complement of the target RNA's anti-repeat regions.
    # Note: These are RNA-RNA pairings.
    rptr_long_ar = reverse_complement(long_anti_repeat_rna).replace('T', 'U')
    rptr_short_ar = reverse_complement(short_anti_repeat_rna).replace('T', 'U')

    reprogrammed_tracrRNA = (
        rptr_long_ar +
        BHCAS12B_TRACR_SCAFFOLD["scaffold_part_1"] +
        rptr_short_ar +
        BHCAS12B_TRACR_SCAFFOLD["scaffold_part_2"]
    )
    
    # 3. Construct the optimized dsDNA activator
    # The target strand (TS) contains the guide sequence.
    # The non-target strand (NTS) is its complement, with the PAM at the 5' end.
    
    activator_target_strand_full = guide_sequence_dna
    activator_non_target_strand_full = ACTIVATOR_PAM + reverse_complement(guide_sequence_dna)
    
    # Apply optimal truncations for enhanced activity
    dsDNA_activator_TS = activator_target_strand_full[:-TS_TRUNCATION]
    dsDNA_activator_NTS = activator_non_target_strand_full[:-NTS_TRUNCATION]

    # 4. Assemble the response
    notes = (
        "Assay components designed based on the PUMA system (Jiao et al., Nat Commun 15, 5909 (2024)) "
        "for direct, PAM-independent RNA detection. Activator is truncated for optimal performance."
    )
    
    return AssayResponse(
        reprogrammed_tracrRNA=reprogrammed_tracrRNA,
        dsDNA_activator_NTS=dsDNA_activator_NTS,
        dsDNA_activator_TS=dsDNA_activator_TS,
        ssDNA_reporter=SSDNA_REPORTER,
        enzyme="BhCas12b",
        protocol="Incubate components at 42°C. Monitor fluorescence at Excitation: 485nm, Emission: 528nm.",
        notes=notes,
    )

if __name__ == "__main__":
    import uvicorn
    # Add a simple health check endpoint
    @app.get("/health")
    def health_check():
        return {"status": "ok"}
    uvicorn.run(app, host="0.0.0.0", port=8006) 
 
from pydantic import BaseModel
import logging
import random

# --- Constants based on Jiao et al., 2024, Nature Communications ---
# DOI: 10.1038/s41467-024-50243-x

# Scaffold for the BhCas12b reprogrammed tracrRNA (Rptr).
# The anti-repeat sequence will be inserted into this scaffold.
# This structure is derived from the paper's description and supplementary figures.
# The 'X's represent the long anti-repeat, and 'Y's represent the short anti-repeat.
BHCAS12B_TRACR_SCAFFOLD = {
    "long_anti_repeat": "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",  # 31 nt
    "scaffold_part_1": "GUGAUAGUGU",
    "short_anti_repeat": "YYYYY", # 5 nt
    "scaffold_part_2": "UAAUUCUCUACUUUUGU",
}

# The PAM sequence for the dsDNA activator for BhCas12b
ACTIVATOR_PAM = "ATTG" # Using 'G' for 'N' for simplicity

# The fluorescent reporter
SSDNA_REPORTER = "5'-FAM-TTATT-3'-IABkFQ"

# Optimal Truncation for activator based on the paper's findings (Fig. 4)
NTS_TRUNCATION = 6
TS_TRUNCATION = 2

# --- Utility Functions ---

def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, base) for base in reversed(seq))

# --- FastAPI App ---

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Diagnostic Finder Service",
    description="Designs CRISPR-Cas12-based diagnostic assays using the PUMA system.",
    version="1.1.0",
)

class AssayRequest(BaseModel):
    target_rna_sequence: str
    target_gene_name: str

class AssayResponse(BaseModel):
    reprogrammed_tracrRNA: str
    dsDNA_activator_NTS: str
    dsDNA_activator_TS: str
    ssDNA_reporter: str
    enzyme: str
    protocol: str
    notes: str

@app.post("/design_diagnostic_assay", response_model=AssayResponse)
async def design_diagnostic_assay(request: AssayRequest):
    """
    Designs a complete diagnostic assay package based on the PUMA system
    (Programmable tracrRNAs Unlock Motif-independent detection of ribonucleic Acids).
    
    This endpoint takes a target RNA sequence and returns all the necessary
    components to run a fluorescent-based detection assay using BhCas12b.
    """
    logger.info(f"Received assay design request for gene: {request.target_gene_name}")
    
    # --- Input Validation ---
    # The Rptr needs to bind a 31nt long anti-repeat and have a guide sequence.
    # We'll use a 20nt guide for this implementation. Total required length: 51nt.
    required_len = len(BHCAS12B_TRACR_SCAFFOLD["long_anti_repeat"]) + 20
    if len(request.target_rna_sequence) < required_len:
        raise HTTPException(
            status_code=400,
            detail=f"Target RNA sequence must be at least {required_len} nucleotides long."
        )

    # --- Design Logic based on Jiao et al., 2024 ---

    # 1. Split the target RNA sequence into the guide and anti-repeat regions.
    # The paper is complex; for a robust first version, we'll make some simplifying assumptions.
    # We will use a 20-nt guide sequence. The preceding 36 nt of the target RNA
    # will form the long (31nt) and short (5nt) anti-repeats.
    
    guide_len = 20
    long_ar_len = len(BHCAS12B_TRACR_SCAFFOLD["long_anti_repeat"])
    short_ar_len = len(BHCAS12B_TRACR_SCAFFOLD["short_anti_repeat"])
    
    # The guide is at the 3' end of the target region
    guide_sequence_rna = request.target_rna_sequence[-guide_len:]
    guide_sequence_dna = guide_sequence_rna.replace('U', 'T') # Convert to DNA for the activator

    # The anti-repeats are upstream of the guide
    anti_repeat_rna_region = request.target_rna_sequence[-(guide_len + long_ar_len + short_ar_len):-guide_len]
    
    long_anti_repeat_rna = anti_repeat_rna_region[:long_ar_len]
    short_anti_repeat_rna = anti_repeat_rna_region[long_ar_len:]

    # 2. Construct the reprogrammed tracrRNA (Rptr)
    # We need the complement of the target RNA's anti-repeat regions.
    # Note: These are RNA-RNA pairings.
    rptr_long_ar = reverse_complement(long_anti_repeat_rna).replace('T', 'U')
    rptr_short_ar = reverse_complement(short_anti_repeat_rna).replace('T', 'U')

    reprogrammed_tracrRNA = (
        rptr_long_ar +
        BHCAS12B_TRACR_SCAFFOLD["scaffold_part_1"] +
        rptr_short_ar +
        BHCAS12B_TRACR_SCAFFOLD["scaffold_part_2"]
    )
    
    # 3. Construct the optimized dsDNA activator
    # The target strand (TS) contains the guide sequence.
    # The non-target strand (NTS) is its complement, with the PAM at the 5' end.
    
    activator_target_strand_full = guide_sequence_dna
    activator_non_target_strand_full = ACTIVATOR_PAM + reverse_complement(guide_sequence_dna)
    
    # Apply optimal truncations for enhanced activity
    dsDNA_activator_TS = activator_target_strand_full[:-TS_TRUNCATION]
    dsDNA_activator_NTS = activator_non_target_strand_full[:-NTS_TRUNCATION]

    # 4. Assemble the response
    notes = (
        "Assay components designed based on the PUMA system (Jiao et al., Nat Commun 15, 5909 (2024)) "
        "for direct, PAM-independent RNA detection. Activator is truncated for optimal performance."
    )
    
    return AssayResponse(
        reprogrammed_tracrRNA=reprogrammed_tracrRNA,
        dsDNA_activator_NTS=dsDNA_activator_NTS,
        dsDNA_activator_TS=dsDNA_activator_TS,
        ssDNA_reporter=SSDNA_REPORTER,
        enzyme="BhCas12b",
        protocol="Incubate components at 42°C. Monitor fluorescence at Excitation: 485nm, Emission: 528nm.",
        notes=notes,
    )

if __name__ == "__main__":
    import uvicorn
    # Add a simple health check endpoint
    @app.get("/health")
    def health_check():
        return {"status": "ok"}
    uvicorn.run(app, host="0.0.0.0", port=8006) 