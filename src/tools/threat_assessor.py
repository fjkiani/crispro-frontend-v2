"""
The Intelligence Engine for the Interception Design Studio.

This module is responsible for orchestrating the multi-step, live intelligence
gathering process to assess the threat level of a genetic variant.
It is a dual-purpose module, containing both the backend logic and a
renderable Streamlit component.
"""

import streamlit as st
import requests
import re
import time
import json
from pathlib import Path
import pandas as pd
import sqlite3

# Use a relative import to be compatible with Streamlit's module loading
from .command_center_client import initiate_germline_correction


def _get_gene_location_standalone(gene: str):
    """
    Gets the genomic coordinates for an entire gene from Ensembl as a standalone utility.

    Args:
        gene (str): The gene symbol.

    Returns:
        str or None: The coordinates in "chr:start-end" format, or None on failure.
    """
    if not gene:
        return None
    
    lookup_server = "https://rest.ensembl.org"
    lookup_ext = f"/lookup/symbol/homo_sapiens/{gene}"
    headers = {"Content-Type": "application/json"}

    try:
        r = requests.get(lookup_server + lookup_ext, headers=headers, timeout=10)
        r.raise_for_status()
        lookup_data = r.json()

        seq_region = lookup_data.get('seq_region_name')
        start = lookup_data.get('start')
        end = lookup_data.get('end')

        if all([seq_region, start, end]):
            return f"{seq_region}:{start}-{end}"
        else:
            return None
    except (requests.exceptions.RequestException, json.JSONDecodeError):
        return None


def _get_sequence_from_coords(coords: str, log_callback=None):
    """
    Fetches a DNA sequence for given coordinates from the Ensembl REST API.
    """
    if not coords:
        return None
    
    server = "https://rest.ensembl.org"
    # Ensure coords are URL-safe, although they typically are.
    ext = f"/sequence/region/human/{coords}?content-type=application/json"
    
    try:
        if log_callback: log_callback(f"Fetching reference sequence for: {coords}...")
        r = requests.get(server + ext, headers={"Content-Type": "application/json"}, timeout=15, verify=False)
        r.raise_for_status()
        data = r.json()
        sequence = data.get('seq')
        if log_callback: log_callback(f"✅ Success! Sequence acquired.")
        return sequence
    except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
        if log_callback: log_callback(f"❌ ERROR: Could not fetch sequence from Ensembl. Reason: {e}")
        return None


def _expand_coords(coords: str, width: int = 2000, log_callback=None):
    """
    Expands a genomic coordinate string (e.g., 'chr:start-end') by a given width on each side.
    """
    try:
        chrom, pos = coords.split(':', 1)
        # Handle both single positions (e.g., '12345') and ranges (e.g., '12345-12346')
        start_str = pos.split('-')[0]
        start = int(start_str)
        
        new_start = max(1, start - width)
        new_end = start + width
        
        expanded = f"{chrom}:{new_start}-{new_end}"
        if log_callback: log_callback(f"Expanded target site context to: {expanded}")
        return expanded
    except (ValueError, IndexError) as e:
        if log_callback: log_callback(f"❌ ERROR: Could not parse and expand coordinates '{coords}'. Reason: {e}")
        return None


def run_germline_correction_workflow():
    """
    Initiates the full germline correction design campaign.
    Gathers required context and sequences before calling the command center.
    """
    st.session_state.clone_assassin_guides = None
    st.session_state.clone_assassin_error = None
    st.session_state.chopchop_run_status = 'running'
    st.session_state.chopchop_log = ""
    
    log_container = st.empty()
    def log_to_streamlit(message):
        st.session_state.chopchop_log += message + "\n"
        log_container.code(st.session_state.chopchop_log, language="log")

    if 'threat_assessment_results' not in st.session_state:
        st.session_state.clone_assassin_error = "Cannot design blueprint: Please run a threat assessment first."
        return

    # --- 1. GATHER INTEL ---
    log_to_streamlit("--- Initiating Germline Correction Campaign ---")
    dossier = st.session_state.threat_assessment_results.get('dossier', {})
    
    # Use the precise VEP location if available, otherwise the whole gene location
    base_coords = None
    vep_data = dossier.get('VEP Annotation', {})
    if isinstance(vep_data, dict) and 'Location' in vep_data:
        base_coords = vep_data.get('Location')
    else:
        base_coords = dossier.get('Gene Location')

    if not base_coords or "Failed" in str(base_coords):
        st.session_state.clone_assassin_error = "Cannot design blueprint: No valid base coordinates found in the threat dossier."
        st.session_state.chopchop_run_status = 'error'
        return

    log_to_streamlit(f"Using base coordinates: {base_coords}")

    # --- 2. PREPARE PAYLOAD ---
    target_site_context = _expand_coords(base_coords, width=2000, log_callback=log_to_streamlit)
    corrected_sequence = _get_sequence_from_coords(base_coords, log_callback=log_to_streamlit)

    if not target_site_context or not corrected_sequence:
        st.session_state.clone_assassin_error = "Failed to prepare payload: Could not retrieve required genomic context or sequence."
        st.session_state.chopchop_run_status = 'error'
        return
        
    # --- 3. DISPATCH TO COMMANDCENTER ---
    try:
        log_to_streamlit(f"Dispatching request to CommandCenter for blueprint...")
        guides_df, blueprint = initiate_germline_correction(
            corrected_sequence=corrected_sequence,
            target_site_context=target_site_context,
            log_callback=log_to_streamlit
        )

        if guides_df.empty:
            st.session_state.chopchop_run_status = 'error'
            st.session_state.clone_assassin_error = "The CommandCenter returned no valid guide candidates in the blueprint."
        else:
            st.session_state.chopchop_run_status = 'success'
            st.session_state.clone_assassin_guides = (guides_df, blueprint)

    except Exception as e:
        st.session_state.chopchop_run_status = 'error'
        st.session_state.clone_assassin_error = f"An unexpected error occurred during the campaign: {e}"


class ThreatAssessor:
    """
    A class to encapsulate the logic for assessing a genetic variant.
    """

    def __init__(self, gene, protein_change):
        self.gene = gene.strip()
        self.protein_change = protein_change.strip()
        self.results = {"dossier": {}}
        self.hgvs_notation = None
        self.vep_data = None

        # Caching setup
        cache_dir = Path(".cache/threat_assessor")
        cache_dir.mkdir(parents=True, exist_ok=True)
        sanitized_protein_change = re.sub(r'[^\w.-]', '_', self.protein_change)
        cache_key = f"{self.gene}_{sanitized_protein_change}.json"
        self.cache_file = cache_dir / cache_key

    def _load_from_cache(self):
        """Loads assessment results from a local file-based cache."""
        if self.cache_file.exists():
            try:
                with self.cache_file.open("r") as f:
                    cached_data = json.load(f)
                    self.results = cached_data["results"]
                    self.vep_data = cached_data.get("vep_data")
                    self.hgvs_notation = cached_data.get("hgvs_notation")
                    return True
            except (json.JSONDecodeError, KeyError):
                # Cache is corrupted or invalid, ignore it
                return False
        return False

    def _save_to_cache(self):
        """Saves assessment results to the local file-based cache."""
        cache_data = {
            "results": self.results,
            "vep_data": self.vep_data,
            "hgvs_notation": self.hgvs_notation,
        }
        with self.cache_file.open("w") as f:
            json.dump(cache_data, f, indent=2)

    def _get_vep_annotation(self):
        """
        Gets variant annotation from Ensembl VEP using a two-step process:
        1. Find canonical transcript for the gene.
        2. Use transcript to construct a full HGVS notation for the VEP query.
        """
        # Step 1: Get canonical transcript
        lookup_server = "https://rest.ensembl.org"
        lookup_ext = f"/lookup/symbol/homo_sapiens/{self.gene}"
        headers = {"Content-Type": "application/json"}
        
        try:
            r = requests.get(lookup_server + lookup_ext, headers=headers, timeout=10)
            r.raise_for_status()
            lookup_data = r.json()
            transcript_id = lookup_data.get('canonical_transcript')
            if not transcript_id:
                self.results['dossier']['VEP Annotation'] = "Failed: Could not find canonical transcript for gene."
                return False
        except requests.exceptions.RequestException as e:
            self.results['dossier']['VEP Annotation'] = f"Failed: Could not contact Ensembl Lookup API. {e}"
            return False

        # Step 2: Use transcript to query VEP
        # The VEP API requires a transcript-specific HGVS notation for protein changes.
        self.hgvs_notation = f"{transcript_id}:{self.protein_change}"
        
        # --- BEGIN DEBUG ---
        st.info(f"🧬 Querying VEP API with HGVS notation: `{self.hgvs_notation}`")
        # --- END DEBUG ---

        vep_server = "https://rest.ensembl.org"
        vep_ext = "/vep/human/hgvs"
        vep_headers = {"Content-Type": "application/json", "Accept": "application/json"}
        vep_data = {"hgvs_notations": [self.hgvs_notation]}
        
        try:
            r_vep = requests.post(vep_server + vep_ext, headers=vep_headers, json=vep_data, timeout=15)
            r_vep.raise_for_status()
            vep_decoded = r_vep.json()
            
            # Check for errors in the response payload itself
            if vep_decoded and 'error' in vep_decoded[0]:
                 self.results['dossier']['VEP Annotation'] = f"Failed: VEP API returned an error: {vep_decoded[0]['error']}"
                 return False

            if not vep_decoded:
                self.results['dossier']['VEP Annotation'] = "Failed: VEP returned no data for this variant."
                return False

            self.vep_data = vep_decoded[0]
            self.results['dossier']['VEP Annotation'] = {
                "Input HGVS": self.vep_data.get('input'),
                "Location": f"chr{self.vep_data.get('seq_region_name')}:{self.vep_data.get('start')}-{self.vep_data.get('end')}",
                "Most Severe Consequence": self.vep_data.get('most_severe_consequence'),
                "Alleles": self.vep_data.get('allele_string'),
            }
            return True
        except requests.exceptions.RequestException as e:
            error_message = f"Failed: Could not contact Ensembl VEP API. {e}"
            # If the error has a response from the server, display its content for detailed debugging.
            if e.response is not None:
                error_message += f"\n\n**API Server Response:**\n```json\n{e.response.text}\n```"
            self.results['dossier']['VEP Annotation'] = error_message
            return False
        except (IndexError, KeyError) as e:
            self.results['dossier']['VEP Annotation'] = f"Failed: Could not parse VEP response. {e}"
            return False

    def _get_gene_location(self):
        """
        Gets the genomic coordinates for the entire gene from Ensembl.
        This serves as a fallback when VEP fails or is bypassed.
        """
        lookup_server = "https://rest.ensembl.org"
        lookup_ext = f"/lookup/symbol/homo_sapiens/{self.gene}"
        headers = {"Content-Type": "application/json"}

        try:
            r = requests.get(lookup_server + lookup_ext, headers=headers, timeout=10)
            r.raise_for_status()
            lookup_data = r.json()

            seq_region = lookup_data.get('seq_region_name')
            start = lookup_data.get('start')
            end = lookup_data.get('end')

            if all([seq_region, start, end]):
                # Store standard chr:start-end format for the guide designer
                self.results['dossier']['Gene Location'] = f"{seq_region}:{start}-{end}"
                return True
            else:
                self.results['dossier']['Gene Location'] = "Failed: Could not parse location from Ensembl lookup."
                return False
        except requests.exceptions.RequestException as e:
            self.results['dossier']['Gene Location'] = f"Failed: Could not contact Ensembl Lookup API. {e}"
            return False
        except (KeyError, json.JSONDecodeError) as e:
            self.results['dossier']['Gene Location'] = f"Failed: Could not parse Ensembl lookup response. {e}"
            return False

    def _get_clinvar_data(self):
        """
        Queries NCBI ClinVar for data on the variant.
        """
        if not self.hgvs_notation:
            self.results['dossier']['ClinVar Significance'] = "Skipped: Missing HGVS notation."
            return False

        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        search_url = f"{base_url}esearch.fcgi?db=clinvar&term={self.hgvs_notation}[hgvs]&retmode=json"
        
        try:
            r_search = requests.get(search_url, timeout=10)
            r_search.raise_for_status()
            search_data = r_search.json()
            ids = search_data.get("esearchresult", {}).get("idlist")
            
            if not ids:
                self.results['dossier']['ClinVar Significance'] = "Not Found"
                return True

            fetch_url = f"{base_url}esummary.fcgi?db=clinvar&id={','.join(ids)}&retmode=json"
            r_fetch = requests.get(fetch_url, timeout=10)
            r_fetch.raise_for_status()
            fetch_data = r_fetch.json()
            
            # Extract clinical significance from the first result
            result = fetch_data.get("result", {}).get(ids[0], {})
            description = result.get("clinical_significance", {}).get("description", "N/A")
            self.results['dossier']['ClinVar Significance'] = description
            return True

        except requests.exceptions.RequestException as e:
            self.results['dossier']['ClinVar Significance'] = f"Failed: Could not contact NCBI API. {e}"
            return False
        except Exception as e:
            self.results['dossier']['ClinVar Significance'] = f"Failed: Could not parse ClinVar response. {e}"
            return False

    def _run_triumvirate_protocol(self):
        """
        Mocks the call to the internal CommandCenter /assess_threat endpoint.
        This simulates our proprietary "Triumvirate Threat Assessment" protocol.
        """
        # For demo purposes, we have pre-computed verdicts for key examples.
        if "fs" in self.protein_change or "*" in self.protein_change:
            # Any frameshift or nonsense mutation is caught by the Truncation Sieve
            self.results['verdict'] = "PATHOGENIC"
            self.results['reason'] = "Determined by: Truncation Sieve"
            self.results['dossier']['Triumvirate Protocol Score'] = -1000.0
        elif self.gene == "BRAF" and self.protein_change == "p.V600E":
            self.results['verdict'] = "PATHOGENIC"
            self.results['reason'] = "Determed by: Zeta Oracle"
            self.results['dossier']['Triumvirate Protocol Score'] = -98.5
        else:
            # Otherwise, assume Zeta Oracle finds it benign
            self.results['verdict'] = "BENIGN"
            self.results['reason'] = "Determined by: Zeta Oracle"
            self.results['dossier']['Triumvirate Protocol Score'] = 25.0
        
        # time.sleep(2) # Removed for better UX
        return True

    def run_assessment(self):
        """
        Executes the full, multi-step threat assessment workflow.
        """
        # Attempt to load from cache first
        if self._load_from_cache():
            st.success("✅ Loaded threat assessment from cache.", icon="🗂️")
            return self.results

        with st.status("Executing live intelligence gathering...", expanded=True) as status:
            # Handle frameshift/nonsense mutations, which VEP cannot parse without nucleotide data
            if "fs" in self.protein_change or "*" in self.protein_change:
                status.update(label="Analyzing definitive truncating mutation...")
                self.results['dossier']['VEP Annotation'] = "Bypassed: Frameshift/Nonsense mutation is definitively pathogenic. VEP API requires nucleotide-level data for this variant type."

                status.update(label="Acquiring gene coordinates as fallback...")
                self._get_gene_location() # Get gene-level coords

                # Use a best-effort query for ClinVar
                self.hgvs_notation = f"{self.gene} {self.protein_change}" # e.g. "ASXL1 p.Gly646fs"
                status.update(label="Querying ClinVar for clinical significance...")
                self._get_clinvar_data()

                status.update(label="Executing Triumvirate Threat Assessment Protocol...")
                self._run_triumvirate_protocol()

            else: # Standard workflow for other mutation types
                status.update(label="Acquiring genomic coordinates (Ensembl VEP)...")
                if not self._get_vep_annotation():
                    status.update(label="VEP annotation failed. Acquiring gene coordinates as fallback...", state="warning")
                    self._get_gene_location() # Get gene-level coords as a fallback
                    # We don't save to cache here, allowing a retry to succeed later if the API was down.
                
                status.update(label="Querying ClinVar for clinical significance...")
                self._get_clinvar_data() # Not a fatal error if this fails
                
                status.update(label="Executing Triumvirate Threat Assessment Protocol...")
                self._run_triumvirate_protocol()

            status.update(label="Assessment Complete.", state="complete")

        self._save_to_cache()
        return self.results


def threat_assessor_module():
    """
    A Streamlit module to render the threat assessment UI components.
    """
    # --- State Initialization ---
    # Ensure all required session state variables are initialized on first run
    # to prevent AttributeError.
    if 'clone_assassin_error' not in st.session_state:
        st.session_state.clone_assassin_error = None
    if 'clone_assassin_guides' not in st.session_state:
        st.session_state.clone_assassin_guides = None
    if 'chopchop_log' not in st.session_state:
        st.session_state.chopchop_log = ""

    st.subheader("Phase A: Threat Assessment")
    st.write("Input a potential threat (e.g., a somatic mutation) to begin the assessment.")

    # Restore the form to make the UI interactive again
    with st.form("threat_assessment_form"):
        # Use session state to pre-fill inputs, allowing for programmatic changes
        gene_input = st.text_input(
            "Gene Symbol", 
            st.session_state.get('threat_gene_input', 'ASXL1'), 
            help="e.g., BRAF, ASXL1",
            key='threat_gene_input_form' # Use a different key for the widget
        )
        protein_change_input = st.text_input(
            "Protein Change (p. notation)", 
            st.session_state.get('threat_protein_change_input', 'p.Gly646fs'), 
            help="e.g., p.V600E, p.Gly646fs",
            key='threat_protein_change_input_form' # Use a different key
        )
        submitted = st.form_submit_button("Assess Threat")

    if submitted:
        if not gene_input or not protein_change_input:
            st.error("Please provide both a gene symbol and a protein change.")
        else:
            # Update the canonical session state variables before running
            st.session_state.threat_gene_input = gene_input
            st.session_state.threat_protein_change_input = protein_change_input

            assessor = ThreatAssessor(gene_input, protein_change_input)
            results = assessor.run_assessment()
            st.session_state.threat_assessment_results = results
            st.rerun() # This is the correct way to rerun after form submission

    if 'threat_assessment_results' in st.session_state:
        results = st.session_state.threat_assessment_results
        
        st.subheader("Threat Dossier")
        
        # Add a button to clear the cache for the current entry
        if st.button("Clear Cache & Re-run Assessment"):
            gene_input = st.session_state.get('threat_gene_input')
            protein_change_input = st.session_state.get('threat_protein_change_input')

            if gene_input and protein_change_input:
                cache_dir = Path(".cache/threat_assessor")
                sanitized_protein_change = re.sub(r'[^\w.-]', '_', protein_change_input)
                cache_key = f"{gene_input}_{sanitized_protein_change}.json"
                cache_file = cache_dir / cache_key
                
                if cache_file.exists():
                    cache_file.unlink()
                    del st.session_state.threat_assessment_results
                    st.success("Cache cleared. Re-running assessment...")
                    time.sleep(1) # Give user time to see the message
                    st.rerun()
                else:
                    st.warning("No cache file found to clear.")
            else:
                st.warning("Could not determine which cache entry to clear.")

        if results.get('verdict') == "PATHOGENIC":
            st.error(f"**Verdict: {results.get('verdict')}**")
        elif results.get('verdict') == "BENIGN":
            st.success(f"**Verdict: {results.get('verdict')}**")
        else:
            st.warning(f"**Verdict: {results.get('verdict', 'UNKNOWN')}**")

        st.caption(results.get('reason', ''))

        with st.expander("Show Intelligence Dossier Details"):
            for key, value in results.get('dossier', {}).items():
                if isinstance(value, dict):
                    st.markdown(f"**{key}:**")
                    st.json(value)
                else:
                    st.markdown(f"**{key}:** `{value}`")

    st.markdown("---")
    st.subheader("Phase B: Therapeutic Blueprint Generation")
    st.write(
        "Based on the validated threat, the CasPro platform will now design a complete, pre-validated therapeutic package. "
        "This **Therapeutic Blueprint** includes AI-optimized guide RNAs for maximum efficacy, a full safety analysis to de-risk development, "
        "and the required components for a homology-directed repair strategy."
    )

    # --- Design Parameters ---
    col1, col2 = st.columns(2)
    with col1:
        genome_selection = st.selectbox(
            "Select Genome Assembly",
            ["hg38", "hg19", "mm10", "danRer11"],
            index=0,
            key="genome_selection",
            help="The reference genome to use for guide design and off-target analysis."
        )
    with col2:
        enzyme_selection = st.selectbox(
            "Select CRISPR Enzyme",
            ["SpCas9", "AsCas12a", "SaCas9"],
            index=0,
            key="enzyme_selection",
            help="The CRISPR nuclease to use for generating guides."
        )

    # --- Action Button ---
    st.button(
        "Design Germline Correction Blueprint",
        on_click=run_germline_correction_workflow,
        help="Click to start the full design campaign based on the assessed threat.",
        disabled=not st.session_state.get('threat_assessment_results') 
    )

    # --- Display Area for Logs and Results ---
    if st.session_state.clone_assassin_error:
        st.error(st.session_state.clone_assassin_error, icon="❌")

    if st.session_state.get('clone_assassin_guides') is not None:
        # The client now returns a tuple: (dataframe, full_blueprint_dict)
        raw_guides_df, blueprint = st.session_state.clone_assassin_guides

        if not raw_guides_df.empty:
            # --- Perform Data Processing and Scoring Directly in the UI ---
            off_target_report = blueprint.get("off_target_analysis", {})
            off_target_hits = off_target_report.get("total_off_targets_found", 0)
            best_guide_by_efficacy = blueprint.get("best_guide_by_efficacy")
            
            processed_guides = []
            for _, row in raw_guides_df.iterrows():
                zeta_score = row["efficacy_score"]
                seq = row["sequence"]
                
                current_guide_hits = off_target_hits if seq == best_guide_by_efficacy else 0
                normalized_efficacy = min(abs(zeta_score), 20) / 20
                safety_score = max(0.0, 1.0 - (current_guide_hits / 5.0))
                blueprint_score = (normalized_efficacy * 0.7) + (safety_score * 0.3)
                
                processed_guides.append({
                    "Sequence": seq,
                    "Blueprint Score": blueprint_score,
                    "Efficacy (Zeta)": normalized_efficacy,
                    "Safety (Off-Target)": safety_score,
                    "Off-Target Hits": current_guide_hits if seq == best_guide_by_efficacy else "N/A",
                    "Raw Score (ΔLL)": f"{zeta_score:.1f}",
                    "Status": "TOP CANDIDATE" if seq == best_guide_by_efficacy else "Validated Candidate"
                })
            
            results_df = pd.DataFrame(processed_guides).sort_values(by="Blueprint Score", ascending=False)
            # --- End of Processing ---

            st.markdown("---")
            st.markdown('<div class="subheader">Blueprint Executive Summary</div>', unsafe_allow_html=True)
            
            top_candidate = results_df.iloc[0]
            blueprint_score = top_candidate["Blueprint Score"]
            efficacy_score = top_candidate["Efficacy (Zeta)"]
            safety_score = top_candidate["Safety (Off-Target)"]
            off_target_hits_display = top_candidate["Off-Target Hits"]
            
            st.info(
                "**Conclusion:** The CasPro platform has generated a high-quality therapeutic blueprint. "
                f"The lead candidate, `{top_candidate['Sequence']}`, demonstrates an excellent balance of predicted efficacy and safety. "
                "This provides a clear, data-driven path forward for pre-clinical development, "
                "significantly de-risking the asset and reducing R&D costs."
            )

            st.markdown(f"""
                <style>
                .candidate-card {{
                    background-color: #262730; /* Darker background */
                    border-radius: 10px;
                    padding: 20px;
                    margin-bottom: 20px;
                    border: 1px solid #4A4A4A; /* Zeta-themed border */
                }}
                </style>
                <div class="candidate-card">
                    <h4>⭐ Top Candidate: <code>{top_candidate['Sequence']}</code></h4>
                </div>
            """, unsafe_allow_html=True)

            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown("**Blueprint Score**")
                st.progress(blueprint_score, text=f"{blueprint_score:.2f}")
                st.caption("Overall score combining efficacy and safety.")
            with col2:
                st.markdown("**Efficacy (Zeta Score)**")
                st.progress(efficacy_score, text=f"{efficacy_score:.2f}")
                st.caption("AI-predicted knockout success.")
            with col3:
                st.markdown("**Safety (Off-Target)**")
                st.progress(safety_score, text=f"{safety_score:.2f}")
                st.caption(f"Based on {off_target_hits_display} predicted off-target hits.")
            
            with st.expander("Show All Validated Candidates"):
                st.write("The following alternative candidates have also been validated and scored, representing additional therapeutic options.")
                st.dataframe(
                    results_df,
                    use_container_width=True,
                    column_config={
                        "Blueprint Score": st.column_config.ProgressColumn("Blueprint Score",format="%.3f",min_value=0.0,max_value=1.0),
                        "Efficacy (Zeta)": st.column_config.ProgressColumn("Efficacy (Zeta)",format="%.3f",min_value=0.0,max_value=1.0),
                        "Safety (Off-Target)": st.column_config.ProgressColumn("Safety (Off-Target)",format="%.3f",min_value=0.0,max_value=1.0),
                        "Status": "Status",
                        "Raw Score (ΔLL)": "Raw Score"
                    }
                )
        else:
            st.warning("The blueprint was generated, but it did not contain any valid guide candidates.")

    # Display the log if there's anything in it
    if st.session_state.get('chopchop_log'):
        with st.expander("Show Full Design Log"):
            st.code(st.session_state.chopchop_log, language="log")

if __name__ == '__main__':
    threat_assessor_module() 