import sys

import modal

from pydantic import BaseModel
from .image import evo2_image

class VariantAnalysisRequest(BaseModel):
    gene_symbol: str
    protein_change: str
    variant_position: int
    alternative: str
    genome: str
    chromosome: str

app = modal.App("genesis-engine", image=evo2_image)

volume = modal.Volume.from_name("hf_cache", create_if_missing=True)
mount_path = "/root/.cache/huggingface"


@app.function(gpu="H100", volumes={mount_path: volume}, timeout=1000)
def run_brca1_analysis():
    import base64
    from io import BytesIO
    from Bio import SeqIO
    import gzip
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import os
    import seaborn as sns
    from sklearn.metrics import roc_auc_score, roc_curve

    from evo2 import Evo2

    WINDOW_SIZE = 8192

    print("Loading evo2 model...")
    model = Evo2('evo2_7b')
    print("Evo2 model loaded")

    brca1_df = pd.read_excel(
        '/evo2/notebooks/brca1/41586_2018_461_MOESM3_ESM.xlsx',
        header=2,
    )
    brca1_df = brca1_df[[
        'chromosome', 'position (hg19)', 'reference', 'alt', 'function.score.mean', 'func.class',
    ]]

    brca1_df.rename(columns={
        'chromosome': 'chrom',
        'position (hg19)': 'pos',
        'reference': 'ref',
        'alt': 'alt',
        'function.score.mean': 'score',
        'func.class': 'class',
    }, inplace=True)

    # Convert to two-class system
    brca1_df['class'] = brca1_df['class'].replace(['FUNC', 'INT'], 'FUNC/INT')

    with gzip.open('/evo2/notebooks/brca1/GRCh37.p13_chr17.fna.gz', "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_chr17 = str(record.seq)
            break

    # Build mappings of unique reference sequences
    ref_seqs = []
    ref_seq_to_index = {}

    # Parse sequences and store indexes
    ref_seq_indexes = []
    var_seqs = []

    brca1_subset = brca1_df.iloc[:500].copy()

    for _, row in brca1_subset.iterrows():
        p = row["pos"] - 1  # Convert to 0-indexed position
        full_seq = seq_chr17

        ref_seq_start = max(0, p - WINDOW_SIZE//2)
        ref_seq_end = min(len(full_seq), p + WINDOW_SIZE//2)
        ref_seq = seq_chr17[ref_seq_start:ref_seq_end]
        snv_pos_in_ref = min(WINDOW_SIZE//2, p)
        var_seq = ref_seq[:snv_pos_in_ref] + \
            row["alt"] + ref_seq[snv_pos_in_ref+1:]

        # Get or create index for reference sequence
        if ref_seq not in ref_seq_to_index:
            ref_seq_to_index[ref_seq] = len(ref_seqs)
            ref_seqs.append(ref_seq)

        ref_seq_indexes.append(ref_seq_to_index[ref_seq])
        var_seqs.append(var_seq)

    ref_seq_indexes = np.array(ref_seq_indexes)

    print(
        f'Scoring likelihoods of {len(ref_seqs)} reference sequences with Evo 2...')
    ref_scores = model.score_sequences(ref_seqs)

    print(
        f'Scoring likelihoods of {len(var_seqs)} variant sequences with Evo 2...')
    var_scores = model.score_sequences(var_seqs)

    # Subtract score of corresponding reference sequences from scores of variant sequences
    delta_scores = np.array(var_scores) - np.array(ref_scores)[ref_seq_indexes]

    # Add delta scores to dataframe
    brca1_subset[f'evo2_delta_score'] = delta_scores

    y_true = (brca1_subset['class'] == 'LOF')
    auroc = roc_auc_score(y_true, -brca1_subset['evo2_delta_score'])

    # --- Calculate threshold START
    y_true = (brca1_subset["class"] == "LOF")

    fpr, tpr, thresholds = roc_curve(y_true, -brca1_subset["evo2_delta_score"])

    optimal_idx = (tpr - fpr).argmax()

    optimal_threshold = -thresholds[optimal_idx]

    lof_scores = brca1_subset.loc[brca1_subset["class"]
                                  == "LOF", "evo2_delta_score"]
    func_scores = brca1_subset.loc[brca1_subset["class"]
                                   == "FUNC/INT", "evo2_delta_score"]

    lof_std = lof_scores.std()
    func_std = func_scores.std()

    confidence_params = {
        "threshold": optimal_threshold,
        "lof_std": lof_std,
        "func_std": func_std
    }

    print("Confidence params:", confidence_params)

    # --- Calculate threshold END

    plt.figure(figsize=(4, 2))

    # Plot stripplot of distributions
    p = sns.stripplot(
        data=brca1_subset,
        x='evo2_delta_score',
        y='class',
        hue='class',
        order=['FUNC/INT', 'LOF'],
        palette=['#777777', 'C3'],
        size=2,
        jitter=0.3,
    )

    # Mark medians from each distribution
    sns.boxplot(showmeans=True,
                meanline=True,
                meanprops={'visible': False},
                medianprops={'color': 'k', 'ls': '-', 'lw': 2},
                whiskerprops={'visible': False},
                zorder=10,
                x="evo2_delta_score",
                y="class",
                data=brca1_subset,
                showfliers=False,
                showbox=False,
                showcaps=False,
                ax=p)
    plt.xlabel('Delta likelihood score, Evo 2')
    plt.ylabel('BRCA1 SNV class')
    plt.tight_layout()

    buffer = BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    plot_data = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return {'variants': brca1_subset.to_dict(orient="records"), "plot": plot_data, "auroc": auroc}


@app.function()
def brca1_example():
    import base64
    from io import BytesIO
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    print("Running BRCA1 variant analysis with Evo2...")

    # Run inference
    result = run_brca1_analysis.remote()

    if "plot" in result:
        plot_data = base64.b64decode(result["plot"])
        with open("brca1_analysis_plot.png", "wb") as f:
            f.write(plot_data)

        img = mpimg.imread(BytesIO(plot_data))
        plt.figure(figsize=(10, 5))
        plt.imshow(img)
        plt.axis("off")
        plt.show()


def get_genome_sequence(position, genome: str, chromosome: str, window_size=8192):
    import requests

    half_window = window_size // 2
    start = max(0, position - 1 - half_window)
    end = position - 1 + half_window + 1

    print(
        f"Fetching {window_size}bp window around position {position} from UCSC API..")
    print(f"Coordinates: {chromosome}:{start}-{end} ({genome})")

    api_url = f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chromosome};start={start};end={end}"
    response = requests.get(api_url)

    if response.status_code != 200:
        raise Exception(
            f"Failed to fetch genome sequence from UCSC API: {response.status_code}")

    genome_data = response.json()

    if "dna" not in genome_data:
        error = genome_data.get("error", "Unknown error")
        raise Exception(f"UCSC API errpr: {error}")

    sequence = genome_data.get("dna", "").upper()
    expected_length = end - start
    if len(sequence) != expected_length:
        print(
            f"Warning: received sequence length ({len(sequence)}) differs from expected ({expected_length})")

    print(
        f"Loaded reference genome sequence window (length: {len(sequence)} bases)")

    return sequence, start


def analyze_variant(relative_pos_in_window, reference, alternative, window_seq, model):
    var_seq = window_seq[:relative_pos_in_window] + \
        alternative + window_seq[relative_pos_in_window+1:]

    ref_score = model.score_sequences([window_seq])[0]
    var_score = model.score_sequences([var_seq])[0]

    delta_score = var_score - ref_score

    threshold = -0.0009178519
    lof_std = 0.0015140239
    func_std = 0.0009016589

    if delta_score < threshold:
        prediction = "Likely pathogenic"
        confidence = min(1.0, abs(delta_score - threshold) / lof_std)
    else:
        prediction = "Likely benign"
        confidence = min(1.0, abs(delta_score - threshold) / func_std)

    return {
        "reference": reference,
        "alternative": alternative,
        "delta_score": float(delta_score),
        "prediction": prediction,
        "classification_confidence": float(confidence)
    }


@app.cls(gpu="H100", volumes={mount_path: volume}, max_containers=3, retries=2, scaledown_window=120)
class Evo2Model:
    @modal.enter()
    def load_evo2_model(self):
        from evo2 import Evo2
        print("Loading evo2 model...")
        self.model = Evo2('evo2_40b')
        print("Evo2 model loaded.")

    # @modal.method()
    @modal.fastapi_endpoint(method="POST")
    def analyze_single_variant(self, request: VariantAnalysisRequest):
        """
        Analyzes a single genomic variant using the Evo2 model and enriches it with external data.
        """
        import requests
        from bs4 import BeautifulSoup
        import json

        try:
            # Step 1: Perform Evo2 analysis (existing logic)
            variant_position = request.variant_position
            alternative = request.alternative
            genome = request.genome
            chromosome = request.chromosome

            print("Genome:", genome)
            print("Chromosome:", chromosome)
            print("Variant position:", variant_position)
            print("Variant alternative:", alternative)

            WINDOW_SIZE = 8192

            window_seq, seq_start = get_genome_sequence(
                position=variant_position,
                genome=genome,
                chromosome=chromosome,
                window_size=WINDOW_SIZE
            )

            print(f"Fetched genome seauence window, first 100: {window_seq[:100]}")

            relative_pos = variant_position - 1 - seq_start
            print(f"Relative position within window: {relative_pos}")

            if relative_pos < 0 or relative_pos >= len(window_seq):
                raise ValueError(
                    f"Variant position {variant_position} is outside the fetched window (start={seq_start+1}, end={seq_start+len(window_seq)})")

            reference = window_seq[relative_pos]
            print("Reference is: " + reference)

            # Analyze the variant
            evo2_result = analyze_variant(
                relative_pos_in_window=relative_pos,
                reference=reference,
                alternative=alternative,
                window_seq=window_seq,
                model=self.model
            )

            # Step 2: Fetch data from ClinVar
            clinvar_data = self._fetch_clinvar_data(request.gene_symbol, request.protein_change)

            # Step 3: (Placeholder) Fetch data from OncoKB
            oncokb_data = {"level": "TBD", "summary": "OncoKB integration pending."}

            # Step 4: Generate AI Summary
            ai_summary = self._generate_ai_summary(
                request.gene_symbol,
                request.protein_change,
                evo2_result,
                clinvar_data
            )

            # Step 5: Combine all results
            final_result = {
                "evo2_analysis": evo2_result,
                "clinvar_data": clinvar_data,
                "oncokb_data": oncokb_data,
                "ai_summary": ai_summary,
                "input_request": request.dict()
            }

            return final_result
        except Exception as e:
            print(f"Error analyzing variant: {e}")
            return {"error": str(e)}

    def _fetch_clinvar_data(self, gene_symbol: str, protein_change: str) -> dict:
        """
        Fetches variant data from the NCBI ClinVar API using E-utils.
        """
        import requests
        import xml.etree.ElementTree as ET

        print(f"Fetching ClinVar data for {gene_symbol} {protein_change}...")
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        # Construct the search term
        # The format "[gene]" AND "[protein change]" is effective for specific variants
        search_term = f"{gene_symbol}[gene] AND {protein_change}[protein change]"
        
        # 1. ESearch: Find the ClinVar record's UID
        esearch_url = f"{base_url}esearch.fcgi?db=clinvar&term={search_term}&retmode=json"
        
        try:
            response = requests.get(esearch_url)
            response.raise_for_status()
            esearch_data = response.json()
            
            ids = esearch_data.get("esearchresult", {}).get("idlist")
            if not ids:
                return {"error": "Variant not found in ClinVar."}
            
            clinvar_id = ids[0] # Use the first ID found

            # 2. ESummary: Get the summary of the record
            esummary_url = f"{base_url}esummary.fcgi?db=clinvar&id={clinvar_id}&retmode=json"
            response = requests.get(esummary_url)
            response.raise_for_status()
            esummary_data = response.json()
            
            # Extract relevant information
            result = esummary_data.get("result", {}).get(clinvar_id, {})
            title = result.get("title", "N/A")
            clinical_significance = result.get("clinical_significance", {}).get("description", "N/A")
            
            print(f"ClinVar data found: Significance - {clinical_significance}")
            return {
                "clinvar_id": clinvar_id,
                "title": title,
                "clinical_significance": clinical_significance
            }

        except requests.exceptions.RequestException as e:
            print(f"Error fetching data from ClinVar API: {e}")
            return {"error": f"ClinVar API request failed: {e}"}
        except Exception as e:
            print(f"An unexpected error occurred during ClinVar fetch: {e}")
            return {"error": f"An unexpected error occurred: {e}"}

    def _generate_ai_summary(self, gene: str, variant: str, evo2_data: dict, clinvar_data: dict) -> str:
        """
        Generates a concise, human-readable summary of the variant analysis.
        """
        import google.generativeai as genai
        import os

        try:
            # Ensure the API key is configured
            # In a real Modal app, this would be handled by modal.Secret
            api_key = os.getenv("GEMINI_API_KEY")
            if not api_key:
                return "Error: GEMINI_API_KEY not configured."
            genai.configure(api_key=api_key)

            model = genai.GenerativeModel('gemini-1.5-flash')

            # Construct a detailed prompt
            prompt = f"""
            You are an expert clinical genomicist. Synthesize the following data for the variant {gene} {variant} into a single, concise clinical significance summary.

            **Data Point 1: Evo2 Deep Learning Model**
            - Predicted Classification: {evo2_data.get('classification', 'N/A')}
            - Delta Likelihood Score: {evo2_data.get('delta_score', 'N/A')}
            - Confidence: {evo2_data.get('confidence', 'N/A')}

            **Data Point 2: ClinVar Public Database**
            - Clinical Significance: {clinvar_data.get('clinical_significance', 'N/A')}
            - ClinVar Title: {clinvar_data.get('title', 'N/A')}

            **Task:**
            Provide a single paragraph summary (2-3 sentences) integrating these findings. Start with the most definitive information (usually ClinVar) and then use the Evo2 data as supporting or clarifying evidence. State the final conclusion clearly.
            Example: "The BRAF V600E variant is classified as Pathogenic by ClinVar and is a well-established oncogenic driver. This is strongly supported by the Evo2 deep learning model, which predicts a pathogenic classification with high confidence. Therefore, this variant is considered clinically actionable."
            """

            response = model.generate_content(prompt)
            summary = response.text.strip()
            print(f"Generated AI summary for {gene} {variant}: {summary}")
            return summary

        except Exception as e:
            print(f"Error generating AI summary: {e}")
            return f"Error during AI summary generation: {e}"


@app.local_entrypoint()
def main():
    # Example of how you'd call the deployed Modal Function from your client
    import requests
    import json    # brca1_example.remote()

    evo2Model = Evo2Model()

    url = evo2Model.analyze_single_variant.web_url

    payload = {
        "variant_position": 43119628,
        "alternative": "G",
        "genome": "hg38",
        "chromosome": "chr17"
    }

    headers = {
        "Content-Type": "application/json"
    }

    response = requests.post(url, json=payload, headers=headers)
    response.raise_for_status()
    result = response.json()
    print(result)
