import torch
import esm
import asyncio
from typing import List, Dict, Any, Tuple
import logging
import re

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# It's good practice to define the model name as a constant
ESM_MODEL_NAME = "esm2_t33_650M_UR50D"

# Amino Acid Code Conversion Dictionary
AA_CODES = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
}

class ESMClient:
    """
    A client to score protein variants using a local ESM-2 model.
    """
    def __init__(self, model_name: str = ESM_MODEL_NAME):
        """
        Initializes the ESMClient.
        """
        self.model_name = model_name
        self.model = None
        self.alphabet = None
        logger.info(f"Initializing ESMClient for model {self.model_name}...")
        
        try:
            # --- DOCTRINE: DELAYED IMPORT & SECURITY OVERRIDE ---
            # We import transformer_engine *inside* this method. This bypasses the
            # local pre-flight check by Modal, as the import only occurs when the
            # container is actually running in the cloud.
            import transformer_engine
            
            # As per the Commander's directive, we are whitelisting a trusted global
            # from the transformer_engine library to resolve a torch.load security error.
            if hasattr(transformer_engine, 'common') and hasattr(transformer_engine.common, 'recipe') and hasattr(transformer_engine.common.recipe, '_OverrideLinearPrecision'):
                torch.serialization.add_safe_globals([transformer_engine.common.recipe._OverrideLinearPrecision])

            self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(self.model_name)
            self.model.eval()
            if torch.cuda.is_available():
                self.model = self.model.cuda()
                logger.info("ESM-2 model loaded to GPU.")
            else:
                logger.info("ESM-2 model loaded to CPU.")
            self.batch_converter = self.alphabet.get_batch_converter()
        except Exception as e:
            logger.error(f"CRITICAL: Failed to load ESM model {self.model_name}: {e}", exc_info=True)
            logger.error("The ESMClient will be non-operational.")
            self.model = None
            self.alphabet = None
            self.batch_converter = None

    def _get_batch_log_likelihood(self, sequences: List[str]) -> List[float]:
        """
        Calculates the log-likelihood of a batch of sequences.
        """
        if not self.model or not self.batch_converter:
            raise RuntimeError("ESMClient is not operational. The ESM model failed to load during initialization.")

        data = [(f"seq_{i}", seq) for i, seq in enumerate(sequences)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        
        if torch.cuda.is_available():
            batch_tokens = batch_tokens.cuda()

        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[self.model.num_layers], return_contacts=False)
        
        log_likelihoods = []
        logits = results["logits"]
        for i, seq in enumerate(sequences):
            seq_tokens = batch_tokens[i, 1 : len(seq) + 1]
            seq_logits = logits[i, 1 : len(seq) + 1]
            log_probs = torch.log_softmax(seq_logits, dim=-1)
            ll = torch.gather(log_probs, 1, seq_tokens.unsqueeze(1)).squeeze(1)
            log_likelihoods.append(ll.sum().item())
            
        return log_likelihoods

    def get_scores(self, protein_sequence: str, variants: List[str]) -> Dict[str, Dict[str, Any]]:
        """
        Retrieves ESM-2 scores for a list of variants on a given protein sequence.
        """
        results = {}
        mutated_sequences = []
        valid_variants = []

        logger.info(f"Scoring {len(variants)} variants for a protein of length {len(protein_sequence)} with ESM-2 (LIVE).")

        for variant in variants:
            try:
                if variant.startswith("p."):
                    variant = variant[2:]
                
                match = re.match(r"([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])", variant)
                if not match:
                    results[variant] = {"esm_score": -999.0, "error": f"Invalid HGVS format: {variant}"}
                    continue

                ref_aa_code, pos_str, alt_aa_code = match.groups()
                pos = int(pos_str) - 1
                ref_aa = AA_CODES.get(ref_aa_code, ref_aa_code)
                alt_aa = AA_CODES.get(alt_aa_code, alt_aa_code)

                if not (0 <= pos < len(protein_sequence)):
                    error_msg = f"Invalid position: {pos+1} is out of bounds for sequence of length {len(protein_sequence)}."
                    logger.warning(f"ESMClient Validation FAIL for '{variant}': {error_msg}")
                    results[variant] = {"esm_score": -999.0, "error": error_msg}
                    continue
                
                sequence_ref_aa = protein_sequence[pos]
                if sequence_ref_aa != ref_aa:
                    error_msg = f"Reference AA mismatch for '{variant}': HGVS says '{ref_aa}', but sequence has '{sequence_ref_aa}' at position {pos+1}."
                    logger.warning(f"ESMClient Validation FAIL: {error_msg}")
                    results[variant] = {"esm_score": -999.0, "error": error_msg}
                    continue

                mutated_seq_list = list(protein_sequence)
                mutated_seq_list[pos] = alt_aa
                mutated_sequences.append("".join(mutated_seq_list))
                valid_variants.append(variant)

            except (ValueError, IndexError) as e:
                results[variant] = {"esm_score": -999.0, "error": f"Invalid variant format or position: {variant}. Details: {e}"}
            except Exception as e:
                logger.error(f"An unexpected error occurred processing variant {variant}: {e}")
                results[variant] = {"esm_score": -999.0, "error": "An internal error occurred."}

        if not valid_variants:
            return results

        all_sequences = [protein_sequence] + mutated_sequences
        
        try:
            log_likelihoods = self._get_batch_log_likelihood(all_sequences)
            wt_ll = log_likelihoods[0]
            mutant_lls = log_likelihoods[1:]

            for i, variant in enumerate(valid_variants):
                delta_ll = mutant_lls[i] - wt_ll
                results[variant] = {"esm_score": delta_ll}

        except Exception as e:
            logger.error(f"Failed during ESM model inference: {e}", exc_info=True)
            for variant in valid_variants:
                results[variant] = {"esm_score": -997.0, "error": "ESM model inference failed."}

        return results

def main():
    # Example protein sequence (portion of human BRCA1)
    sequence = "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCVGENIOLDTVYTACQLSQTFEKREVSSAQVIRFRKKENKNEIHQGSVLTVNVGNQLLSAQRRNQGKEGKIHHQFSWRGENVDIELQTCPNPTKILKGN"
    
    client = ESMClient()
    test_variants = ["M1A", "D2N", "L3P", "S4W", "A5A", "L3X", "INVALID", "M12345A"]
    scores = client.get_scores(sequence, test_variants)
    
    import json
    print(json.dumps(scores, indent=2))

if __name__ == "__main__":
    main() 
    sequence = "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCVGENIOLDTVYTACQLSQTFEKREVSSAQVIRFRKKENKNEIHQGSVLTVNVGNQLLSAQRRNQGKEGKIHHQFSWRGENVDIELQTCPNPTKILKGN"
    
    client = ESMClient()
    test_variants = ["M1A", "D2N", "L3P", "S4W", "A5A", "L3X", "INVALID", "M12345A"]
    scores = client.get_scores(sequence, test_variants)
    
    import json
    print(json.dumps(scores, indent=2))

if __name__ == "__main__":
    main() 
    sequence = "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCVGENIOLDTVYTACQLSQTFEKREVSSAQVIRFRKKENKNEIHQGSVLTVNVGNQLLSAQRRNQGKEGKIHHQFSWRGENVDIELQTCPNPTKILKGN"
    
    client = ESMClient()
    test_variants = ["M1A", "D2N", "L3P", "S4W", "A5A", "L3X", "INVALID", "M12345A"]
    scores = client.get_scores(sequence, test_variants)
    
    import json
    print(json.dumps(scores, indent=2))

if __name__ == "__main__":
    main() 
    sequence = "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCVGENIOLDTVYTACQLSQTFEKREVSSAQVIRFRKKENKNEIHQGSVLTVNVGNQLLSAQRRNQGKEGKIHHQFSWRGENVDIELQTCPNPTKILKGN"
    
    client = ESMClient()
    test_variants = ["M1A", "D2N", "L3P", "S4W", "A5A", "L3X", "INVALID", "M12345A"]
    scores = client.get_scores(sequence, test_variants)
    
    import json
    print(json.dumps(scores, indent=2))

if __name__ == "__main__":
    main() 