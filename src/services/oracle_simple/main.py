import modal
import os
import re
import torch
from transformers import T5EncoderModel, T5Tokenizer
from fastapi import FastAPI
from modal import Image, Secret, asgi_app

# --- Image Definition ---
image = (
    Image.from_registry("nvidia/cuda:12.1.1-devel-ubuntu22.04", add_python="3.11")
    .run_commands(
        "apt-get update && apt-get install -y git",
        "python -m pip install modal",
        "pip install pandas==2.2.2",
        "pip install transformers==4.41.2",
        "pip install torch==2.3.0",
        "pip install biopython==1.83",
        "pip install fastapi uvicorn",
    )
    .add_local_dir("tools", remote_path="/root/tools")
    .add_local_dir("data", remote_path="/root/data")
)

app = modal.App(
    "zeta-oracle-simple", 
    image=image,
)

# --- Global Variables for Model and Tokenizer ---
model = None
tokenizer = None
AA_CODES = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Thr': 'T'
}

def _load_model():
    """Function to load the model and tokenizer into memory."""
    global model, tokenizer
    if model is None or tokenizer is None:
        print("Loading model and tokenizer...")
        model_name = "Rostlab/prot_t5_xl_uniref50"
        tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
        model = T5EncoderModel.from_pretrained(model_name).to("cuda")
        print("Model and tokenizer loaded.")

def _apply_hgvs_mutation(sequence: str, hgvsp: str) -> str:
    hgvsp_norm = hgvsp if hgvsp.startswith('p.') else f"p.{hgvsp}"
    match = re.match(r"p\.\(?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|[A-Z])\)?", hgvsp_norm)
    if not match: raise ValueError(f"Invalid HGVS: {hgvsp_norm}")
    
    groups = match.groupdict()
    ref = groups['ref']
    pos = int(groups['pos']) - 1
    alt = groups['alt']
    
    ref_aa_one = AA_CODES.get(ref, ref)
    alt_aa_one = AA_CODES.get(alt, alt)

    if pos >= len(sequence) or sequence[pos] != ref_aa_one:
        raise ValueError("Mismatch or out of bounds")

    mutated_sequence = list(sequence)
    mutated_sequence[pos] = alt_aa_one
    return "".join(mutated_sequence)

def _get_canonical_sequence(gene_symbol: str) -> str:
    from Bio import SeqIO
    fasta_path = f"/root/data/reference/{gene_symbol}.fasta"
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found for {gene_symbol}")
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        return str(record.seq)

# --- Modal Function ---
@app.function(
    gpu="A10G",
    secrets=[Secret.from_name("huggingface-secret")],
    timeout=600,
)
def embed(gene_symbol: str, variant: str) -> list[float]:
    _load_model()
    
    canonical_seq = _get_canonical_sequence(gene_symbol)
    mutated_seq = _apply_hgvs_mutation(canonical_seq, variant)

    sequence_spaced = " ".join(list(mutated_seq))
    
    inputs = tokenizer(sequence_spaced, return_tensors="pt").to("cuda")
    with torch.no_grad():
        outputs = model(**inputs)

    embedding = outputs.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()
    return embedding.tolist()

# --- FastAPI App ---
web_app = FastAPI()

@web_app.post("/embed")
def web_embed(item: dict):
    gene = item.get("gene_symbol")
    variant = item.get("variant")
    if not gene or not variant:
        return {"error": "Missing gene_symbol or variant"}, 400
    
    try:
        embedding_vector = embed.remote(gene, variant)
        return {"embedding": embedding_vector, "status": "success"}
    except Exception as e:
        return {"error": str(e)}, 500

@app.function()
@asgi_app()
def fastapi_app():
    return web_app 