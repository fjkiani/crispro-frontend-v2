from Bio.Seq import Seq
from Bio.Data import CodonTable
from loguru import logger

def translate_dna_to_protein(dna_sequence: str) -> str | None:
    """
    Translates a DNA sequence into a protein sequence.
    This function iterates through common genetic code tables to find the one
    that produces the longest translatable sequence with the fewest stop codons,
    making it robust against different biological contexts.
    """
    if not dna_sequence or not isinstance(dna_sequence, str):
        logger.warning("Invalid DNA sequence provided for translation.")
        return None

    best_protein = ""
    min_stop_codons = float('inf')

    # Iterate through a curated list of common codon tables
    # (1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc.)
    for table_id in [1, 2, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16]:
        try:
            # Create a sequence object and translate it
            # We set to_stop=False to translate through stop codons (*)
            # This allows us to count them and assess the quality of the translation.
            protein = str(Seq(dna_sequence).translate(table=table_id, to_stop=False))
            
            # Count the number of stop codons
            stop_codon_count = protein.count('*')

            # We prefer the translation with the fewest stop codons.
            # If counts are equal, we prefer the longer sequence (less likely to be a fragment).
            if stop_codon_count < min_stop_codons:
                min_stop_codons = stop_codon_count
                best_protein = protein
            elif stop_codon_count == min_stop_codons:
                if len(protein) > len(best_protein):
                    best_protein = protein

        except CodonTable.TranslationError as e:
            # This can happen if the sequence length is not a multiple of 3
            logger.trace(f"Translation failed with table {table_id}: {e}")
            continue

    if not best_protein:
        logger.error(f"Failed to translate DNA sequence of length {len(dna_sequence)} using any codon table.")
        return None

    # Before returning, remove the stop codons for a clean protein sequence
    final_protein = best_protein.replace("*", "")
    
    if min_stop_codons > 0:
        logger.warning(f"Translated sequence contained {min_stop_codons} stop codon(s), which were removed.")

    return final_protein 
from Bio.Data import CodonTable
from loguru import logger

def translate_dna_to_protein(dna_sequence: str) -> str | None:
    """
    Translates a DNA sequence into a protein sequence.
    This function iterates through common genetic code tables to find the one
    that produces the longest translatable sequence with the fewest stop codons,
    making it robust against different biological contexts.
    """
    if not dna_sequence or not isinstance(dna_sequence, str):
        logger.warning("Invalid DNA sequence provided for translation.")
        return None

    best_protein = ""
    min_stop_codons = float('inf')

    # Iterate through a curated list of common codon tables
    # (1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc.)
    for table_id in [1, 2, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16]:
        try:
            # Create a sequence object and translate it
            # We set to_stop=False to translate through stop codons (*)
            # This allows us to count them and assess the quality of the translation.
            protein = str(Seq(dna_sequence).translate(table=table_id, to_stop=False))
            
            # Count the number of stop codons
            stop_codon_count = protein.count('*')

            # We prefer the translation with the fewest stop codons.
            # If counts are equal, we prefer the longer sequence (less likely to be a fragment).
            if stop_codon_count < min_stop_codons:
                min_stop_codons = stop_codon_count
                best_protein = protein
            elif stop_codon_count == min_stop_codons:
                if len(protein) > len(best_protein):
                    best_protein = protein

        except CodonTable.TranslationError as e:
            # This can happen if the sequence length is not a multiple of 3
            logger.trace(f"Translation failed with table {table_id}: {e}")
            continue

    if not best_protein:
        logger.error(f"Failed to translate DNA sequence of length {len(dna_sequence)} using any codon table.")
        return None

    # Before returning, remove the stop codons for a clean protein sequence
    final_protein = best_protein.replace("*", "")
    
    if min_stop_codons > 0:
        logger.warning(f"Translated sequence contained {min_stop_codons} stop codon(s), which were removed.")

    return final_protein 
from Bio.Data import CodonTable
from loguru import logger

def translate_dna_to_protein(dna_sequence: str) -> str | None:
    """
    Translates a DNA sequence into a protein sequence.
    This function iterates through common genetic code tables to find the one
    that produces the longest translatable sequence with the fewest stop codons,
    making it robust against different biological contexts.
    """
    if not dna_sequence or not isinstance(dna_sequence, str):
        logger.warning("Invalid DNA sequence provided for translation.")
        return None

    best_protein = ""
    min_stop_codons = float('inf')

    # Iterate through a curated list of common codon tables
    # (1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc.)
    for table_id in [1, 2, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16]:
        try:
            # Create a sequence object and translate it
            # We set to_stop=False to translate through stop codons (*)
            # This allows us to count them and assess the quality of the translation.
            protein = str(Seq(dna_sequence).translate(table=table_id, to_stop=False))
            
            # Count the number of stop codons
            stop_codon_count = protein.count('*')

            # We prefer the translation with the fewest stop codons.
            # If counts are equal, we prefer the longer sequence (less likely to be a fragment).
            if stop_codon_count < min_stop_codons:
                min_stop_codons = stop_codon_count
                best_protein = protein
            elif stop_codon_count == min_stop_codons:
                if len(protein) > len(best_protein):
                    best_protein = protein

        except CodonTable.TranslationError as e:
            # This can happen if the sequence length is not a multiple of 3
            logger.trace(f"Translation failed with table {table_id}: {e}")
            continue

    if not best_protein:
        logger.error(f"Failed to translate DNA sequence of length {len(dna_sequence)} using any codon table.")
        return None

    # Before returning, remove the stop codons for a clean protein sequence
    final_protein = best_protein.replace("*", "")
    
    if min_stop_codons > 0:
        logger.warning(f"Translated sequence contained {min_stop_codons} stop codon(s), which were removed.")

    return final_protein 
from Bio.Data import CodonTable
from loguru import logger

def translate_dna_to_protein(dna_sequence: str) -> str | None:
    """
    Translates a DNA sequence into a protein sequence.
    This function iterates through common genetic code tables to find the one
    that produces the longest translatable sequence with the fewest stop codons,
    making it robust against different biological contexts.
    """
    if not dna_sequence or not isinstance(dna_sequence, str):
        logger.warning("Invalid DNA sequence provided for translation.")
        return None

    best_protein = ""
    min_stop_codons = float('inf')

    # Iterate through a curated list of common codon tables
    # (1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc.)
    for table_id in [1, 2, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16]:
        try:
            # Create a sequence object and translate it
            # We set to_stop=False to translate through stop codons (*)
            # This allows us to count them and assess the quality of the translation.
            protein = str(Seq(dna_sequence).translate(table=table_id, to_stop=False))
            
            # Count the number of stop codons
            stop_codon_count = protein.count('*')

            # We prefer the translation with the fewest stop codons.
            # If counts are equal, we prefer the longer sequence (less likely to be a fragment).
            if stop_codon_count < min_stop_codons:
                min_stop_codons = stop_codon_count
                best_protein = protein
            elif stop_codon_count == min_stop_codons:
                if len(protein) > len(best_protein):
                    best_protein = protein

        except CodonTable.TranslationError as e:
            # This can happen if the sequence length is not a multiple of 3
            logger.trace(f"Translation failed with table {table_id}: {e}")
            continue

    if not best_protein:
        logger.error(f"Failed to translate DNA sequence of length {len(dna_sequence)} using any codon table.")
        return None

    # Before returning, remove the stop codons for a clean protein sequence
    final_protein = best_protein.replace("*", "")
    
    if min_stop_codons > 0:
        logger.warning(f"Translated sequence contained {min_stop_codons} stop codon(s), which were removed.")

    return final_protein 