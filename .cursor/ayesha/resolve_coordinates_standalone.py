#!/usr/bin/env python3
"""
Standalone coordinate resolution for PDGFRA p.S755P
Uses Ensembl VEP API directly (bypasses backend endpoint)
"""
import requests
import json

def resolve_coordinates(gene: str, hgvs_p: str, c_dna: str = None):
    """
    Resolve GRCh38 coordinates from HGVS protein notation using Ensembl VEP.
    
    Returns:
        dict with chrom, pos, ref, alt, transcript
    """
    # Normalize HGVS (remove p. prefix if present)
    if hgvs_p.startswith("p."):
        hgvs_p = hgvs_p[2:]
    
    # Step 1: Get canonical transcript
    lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
    headers = {"Content-Type": "application/json"}
    
    try:
        r = requests.get(lookup_url, headers=headers, timeout=30)
        r.raise_for_status()
        lookup_data = r.json()
        
        transcript_id = lookup_data.get('canonical_transcript')
        if not transcript_id:
            raise ValueError(f"Could not find canonical transcript for {gene}")
        
        print(f"✅ Found canonical transcript: {transcript_id}")
        
    except Exception as e:
        raise ValueError(f"Ensembl lookup failed: {e}")
    
    # Step 2: Query VEP
    vep_url = "https://rest.ensembl.org/vep/human/hgvs"
    vep_headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    vep_hgvs = f"{transcript_id}:p.{hgvs_p}"
    vep_data = {"hgvs_notations": [vep_hgvs]}
    
    try:
        r = requests.post(vep_url, headers=vep_headers, json=vep_data, timeout=30)
        r.raise_for_status()
        vep_result = r.json()
        
        if not vep_result or len(vep_result) == 0:
            raise ValueError(f"VEP returned no data for {gene} {hgvs_p}")
        
        vep_item = vep_result[0]
        
        if 'error' in vep_item:
            raise ValueError(f"VEP error: {vep_item['error']}")
        
        # Extract coordinates
        chrom = vep_item.get('seq_region_name')
        pos = vep_item.get('start')
        end = vep_item.get('end')
        
        # Extract alleles
        allele_string = vep_item.get('allele_string', '')
        alleles = allele_string.split('/') if '/' in allele_string else []
        ref = alleles[0] if len(alleles) > 0 else None
        alt = alleles[1] if len(alleles) > 1 else None
        
        # Fallback: parse from cDNA if provided
        if not ref or not alt:
            if c_dna and '>' in c_dna:
                parts = c_dna.split('>')
                if len(parts) == 2:
                    ref_part = parts[0][-1]
                    alt_part = parts[1]
                    if not ref:
                        ref = ref_part
                    if not alt:
                        alt = alt_part
        
        if not chrom or not pos:
            raise ValueError("Could not extract coordinates from VEP response")
        
        if not ref or not alt:
            raise ValueError("Could not extract ref/alt alleles")
        
        return {
            "chrom": str(chrom).replace('chr', ''),
            "pos": int(pos),
            "ref": ref.upper(),
            "alt": alt.upper(),
            "transcript": transcript_id,
            "assembly": "GRCh38",
            "gene": gene,
            "hgvs_p": f"p.{hgvs_p}",
            "c_dna": c_dna
        }
        
    except Exception as e:
        raise ValueError(f"VEP query failed: {e}")


if __name__ == "__main__":
    # Test with PDGFRA p.S755P
    result = resolve_coordinates("PDGFRA", "S755P", "c.2263T>C")
    print("\n" + "="*60)
    print("COORDINATE RESOLUTION RESULT")
    print("="*60)
    print(json.dumps(result, indent=2))



"""
Standalone coordinate resolution for PDGFRA p.S755P
Uses Ensembl VEP API directly (bypasses backend endpoint)
"""
import requests
import json

def resolve_coordinates(gene: str, hgvs_p: str, c_dna: str = None):
    """
    Resolve GRCh38 coordinates from HGVS protein notation using Ensembl VEP.
    
    Returns:
        dict with chrom, pos, ref, alt, transcript
    """
    # Normalize HGVS (remove p. prefix if present)
    if hgvs_p.startswith("p."):
        hgvs_p = hgvs_p[2:]
    
    # Step 1: Get canonical transcript
    lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
    headers = {"Content-Type": "application/json"}
    
    try:
        r = requests.get(lookup_url, headers=headers, timeout=30)
        r.raise_for_status()
        lookup_data = r.json()
        
        transcript_id = lookup_data.get('canonical_transcript')
        if not transcript_id:
            raise ValueError(f"Could not find canonical transcript for {gene}")
        
        print(f"✅ Found canonical transcript: {transcript_id}")
        
    except Exception as e:
        raise ValueError(f"Ensembl lookup failed: {e}")
    
    # Step 2: Query VEP
    vep_url = "https://rest.ensembl.org/vep/human/hgvs"
    vep_headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    vep_hgvs = f"{transcript_id}:p.{hgvs_p}"
    vep_data = {"hgvs_notations": [vep_hgvs]}
    
    try:
        r = requests.post(vep_url, headers=vep_headers, json=vep_data, timeout=30)
        r.raise_for_status()
        vep_result = r.json()
        
        if not vep_result or len(vep_result) == 0:
            raise ValueError(f"VEP returned no data for {gene} {hgvs_p}")
        
        vep_item = vep_result[0]
        
        if 'error' in vep_item:
            raise ValueError(f"VEP error: {vep_item['error']}")
        
        # Extract coordinates
        chrom = vep_item.get('seq_region_name')
        pos = vep_item.get('start')
        end = vep_item.get('end')
        
        # Extract alleles
        allele_string = vep_item.get('allele_string', '')
        alleles = allele_string.split('/') if '/' in allele_string else []
        ref = alleles[0] if len(alleles) > 0 else None
        alt = alleles[1] if len(alleles) > 1 else None
        
        # Fallback: parse from cDNA if provided
        if not ref or not alt:
            if c_dna and '>' in c_dna:
                parts = c_dna.split('>')
                if len(parts) == 2:
                    ref_part = parts[0][-1]
                    alt_part = parts[1]
                    if not ref:
                        ref = ref_part
                    if not alt:
                        alt = alt_part
        
        if not chrom or not pos:
            raise ValueError("Could not extract coordinates from VEP response")
        
        if not ref or not alt:
            raise ValueError("Could not extract ref/alt alleles")
        
        return {
            "chrom": str(chrom).replace('chr', ''),
            "pos": int(pos),
            "ref": ref.upper(),
            "alt": alt.upper(),
            "transcript": transcript_id,
            "assembly": "GRCh38",
            "gene": gene,
            "hgvs_p": f"p.{hgvs_p}",
            "c_dna": c_dna
        }
        
    except Exception as e:
        raise ValueError(f"VEP query failed: {e}")


if __name__ == "__main__":
    # Test with PDGFRA p.S755P
    result = resolve_coordinates("PDGFRA", "S755P", "c.2263T>C")
    print("\n" + "="*60)
    print("COORDINATE RESOLUTION RESULT")
    print("="*60)
    print(json.dumps(result, indent=2))









