# Old code for CRISPResso2 Analysis

# import streamlit as st
# from pathlib import Path
# import os
# import json
# import subprocess
# import shlex
# import tempfile
# import threading
# import queue
# import time
# import pandas as pd
# import re
# from tools.st_utils import (
#     init_session_state,
#     create_educational_sidebar,
#     apply_custom_css,
#     stream_subprocess_output,
#     ask_llm
# )

# # --- Page Config ---
# st.set_page_config(
#     page_title="CRISPResso2 Analysis",
#     page_icon="📊",
#     layout="wide",
# )

# apply_custom_css()
# init_session_state()
# create_educational_sidebar()


# # --- CRISPResso2 Specific Functions ---

# def verify_crispresso_installation(env_mode, conda_env=None):
#     """Verify that CRISPResso2 is installed and available"""
#     cmd = ["conda", "run", "-n", conda_env] if env_mode == "conda" and conda_env else []
#     cmd.extend(["CRISPResso", "--version"])
#     try:
#         result = subprocess.run(cmd, capture_output=True, text=True, check=True)
#         return True, result.stdout.strip()
#     except (subprocess.CalledProcessError, FileNotFoundError) as e:
#         error_msg = e.stderr.strip() if hasattr(e, 'stderr') else "CRISPResso2 not found. Check installation or Conda environment."
#         return False, f"Error: {error_msg}"

# def build_crispresso_command(params):
#     """Build the CRISPResso2 command based on parameters"""
#     cmd = []
#     if params['env_mode'] == 'conda' and params.get('conda_env'):
#         cmd.extend(["conda", "run", "-n", params['conda_env']])
#     cmd.append("CRISPResso")
#     cmd.extend(["-r1", params['fastq_r1']])
#     if params.get('fastq_r2'): cmd.extend(["-r2", params['fastq_r2']])
#     cmd.extend(["-a", params['amplicon_seq'], "-g", params['guide_seq']])
#     if params.get('experiment_name'): cmd.extend(["-n", params['experiment_name']])
#     if params.get('quant_window_size'): cmd.extend(["--quantification_window_size", str(params['quant_window_size'])])
#     if params.get('min_read_quality', 0) > 0: cmd.extend(["--min_average_read_quality", str(params['min_read_quality'])])
#     if params.get('output_dir'): cmd.extend(["--output_folder", params['output_dir']])
#     if params.get('expected_hdr_seq'): cmd.extend(["--expected_hdr_amplicon_seq", params['expected_hdr_seq']])
#     if params.get('base_editor_output'):
#         cmd.append("--base_editor_output")
#         if params.get('conversion_nuc_from') and params.get('conversion_nuc_to'):
#             cmd.extend(["--conversion_nuc_from", params['conversion_nuc_from'], "--conversion_nuc_to", params['conversion_nuc_to']])
#     if params.get('ignore_substitutions'): cmd.append("--ignore_substitutions_at_ends")
#     if params.get('additional_args'): cmd.extend(shlex.split(params['additional_args']))
#     return cmd

# def extract_output_dir(stdout_text):
#     """Extract the CRISPResso2 output directory from stdout"""
#     match = re.search(r"Results will be available in folder:?\s+(\S+)", stdout_text)
#     return match.group(1) if match else None

# def explain_results_with_llm(result_type, data):
#     """Generate an LLM explanation of results data"""
#     prompt = f"Concisely summarize what this {result_type} data from a CRISPResso2 analysis shows about the editing outcomes for a scientist.\n\nDATA:\n{data}"
#     return ask_llm(prompt)

# # --- UI for CRISPResso Page ---

# def crispresso_page():
#     st.markdown('<div class="main-header">CRISPResso2 Analysis</div>', unsafe_allow_html=True)
#     st.markdown("Analyze FASTQ files to quantify editing outcomes from a CRISPR experiment.")
#     tab1, tab2, tab3 = st.tabs(["Input & Parameters", "Run Analysis", "Results"])
#     with tab1: render_crispresso_input_tab()
#     with tab2: render_crispresso_run_tab()
#     with tab3: render_crispresso_results_tab()

# def render_crispresso_input_tab():
#     st.markdown('<div class="subheader">CRISPResso2 Input Data</div>', unsafe_allow_html=True)

#     # FASTQ Inputs
#     st.markdown("**FASTQ Files:**")
#     fastq_input_method = st.radio("Select FASTQ input method:", ["Upload Files", "Specify File Paths"])
#     if fastq_input_method == "Upload Files":
#         col1, col2 = st.columns(2)
#         with col1: fastq_r1_file = st.file_uploader("FASTQ R1 File", type=["fastq", "fq", "gz"])
#         with col2: fastq_r2_file = st.file_uploader("FASTQ R2 File (Paired-end)", type=["fastq", "fq", "gz"])
#         if fastq_r1_file:
#             upload_dir = Path(tempfile.mkdtemp(prefix="crispresso_uploads_"))
#             r1_path = upload_dir / fastq_r1_file.name
#             with open(r1_path, "wb") as f: f.write(fastq_r1_file.getbuffer())
#             st.session_state.crispresso["fastq_r1"] = str(r1_path)
#             if fastq_r2_file:
#                 r2_path = upload_dir / fastq_r2_file.name
#                 with open(r2_path, "wb") as f: f.write(fastq_r2_file.getbuffer())
#                 st.session_state.crispresso["fastq_r2"] = str(r2_path)
#     else:
#         st.session_state.crispresso["fastq_r1"] = st.text_input("FASTQ R1 File Path", st.session_state.crispresso.get("fastq_r1", ""))
#         st.session_state.crispresso["fastq_r2"] = st.text_input("FASTQ R2 File Path", st.session_state.crispresso.get("fastq_r2", ""))

#     # Sequence Inputs
#     st.markdown("**Sequence Information:**")
#     guide_source = st.radio("Guide RNA Source:", ["Direct Input", "From CHOPCHOP Results"]) if st.session_state.chopchop_results.get("guides") else "Direct Input"
#     if guide_source == "From CHOPCHOP Results":
#         guides = st.session_state.chopchop_results["guides"]
#         guide_options = [f"Guide {i+1}: {g.get('seq', '')}" for i, g in enumerate(guides)]
#         selected_guide_text = st.selectbox("Select Guide RNA from CHOPCHOP:", guide_options)
#         selected_guide_index = guide_options.index(selected_guide_text)
#         st.session_state.crispresso["guide_seq"] = guides[selected_guide_index]['seq']
#     else:
#         st.session_state.crispresso["guide_seq"] = st.text_area("Guide RNA Sequence", st.session_state.crispresso.get("guide_seq", ""), help="Enter guide RNA sequence without PAM.")
#     st.session_state.crispresso["amplicon_seq"] = st.text_area("Amplicon Sequence", st.session_state.crispresso.get("amplicon_seq", ""), help="DNA sequence surrounding your target site.")

#     # Parameters & Environment
#     st.markdown("**CRISPResso2 Parameters & Environment:**")
#     col1, col2 = st.columns(2)
#     with col1:
#         st.session_state.crispresso["env_mode"] = "conda" if st.radio("CRISPResso2 Installation:", ["System PATH", "Conda Environment"], index=1) == "Conda Environment" else "path"
#         if st.session_state.crispresso["env_mode"] == "conda":
#             st.session_state.crispresso["conda_env"] = st.text_input("Conda Environment Name", st.session_state.crispresso.get("conda_env", "crispresso2_env"))
#     with col2:
#         st.session_state.crispresso["experiment_name"] = st.text_input("Experiment Name", st.session_state.crispresso.get("experiment_name", "MyAnalysis"))
#         st.session_state.crispresso["output_dir"] = st.text_input("Output Directory", st.session_state.crispresso.get("output_dir", ""), help="Leave empty for automatic.")

#     with st.expander("Advanced Analysis Options"):
#         render_advanced_crispresso_options()

#     if st.button("Validate Inputs & Proceed"):
#         if not all([st.session_state.crispresso.get(k) for k in ["fastq_r1", "amplicon_seq", "guide_seq"]]):
#             st.error("Please provide FASTQ R1, Amplicon Sequence, and Guide RNA sequence.")
#         else:
#             st.success("Inputs valid. Proceed to the 'Run Analysis' tab.")

# def render_advanced_crispresso_options():
#     # To be filled in with more options
#     st.session_state.crispresso["quant_window_size"] = st.number_input("Quantification Window Size", min_value=1, value=st.session_state.crispresso.get("quant_window_size", 10))
#     st.session_state.crispresso["min_read_quality"] = st.number_input("Minimum Read Quality", min_value=0, max_value=42, value=st.session_state.crispresso.get("min_read_quality", 0))

# def render_crispresso_run_tab():
#     st.subheader("Run CRISPResso2 Analysis")
#     if not all([st.session_state.crispresso.get(k) for k in ["fastq_r1", "amplicon_seq", "guide_seq"]]):
#         st.warning("Please complete the input data in the first tab.")
#         return

#     st.markdown("**Summary of Inputs:**")
#     col1, col2 = st.columns(2)
#     with col1:
#         st.write(f"**R1 FASTQ:** `{Path(st.session_state.crispresso['fastq_r1']).name}`")
#         if st.session_state.crispresso.get('fastq_r2'): st.write(f"**R2 FASTQ:** `{Path(st.session_state.crispresso['fastq_r2']).name}`")
#         env = f"Conda: {st.session_state.crispresso['conda_env']}" if st.session_state.crispresso['env_mode'] == 'conda' else 'System PATH'
#         st.write(f"**Environment:** {env}")
#     with col2:
#         st.write(f"**Amplicon Length:** {len(st.session_state.crispresso['amplicon_seq'])} bp")
#         st.write(f"**Guide Length:** {len(st.session_state.crispresso['guide_seq'])} bp")
#         st.write(f"**Experiment Name:** {st.session_state.crispresso.get('experiment_name')}")

#     cmd = build_crispresso_command(st.session_state.crispresso)
#     st.markdown("**Command Preview:**")
#     st.code(" ".join(cmd), language="bash")

#     if st.button("Run CRISPResso2", type="primary"):
#         progress_placeholder = st.empty()
#         progress_placeholder.info("Starting CRISPResso2...")
#         stdout_container = st.empty()
#         stderr_container = st.empty()

#         process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#         stdout_queue = queue.Queue()
#         stderr_queue = queue.Queue()
#         threading.Thread(target=stream_subprocess_output, args=(process, stdout_queue, stderr_queue), daemon=True).start()

#         stdout_text, stderr_text = "", ""
#         progress_bar = st.progress(0)
#         while process.poll() is None or not stdout_queue.empty() or not stderr_queue.empty():
#             try:
#                 while line := stdout_queue.get_nowait():
#                     stdout_text += line
#                     stdout_container.text_area("Standard Output", stdout_text, height=300)
#                     if "Processing" in line:
#                         progress_bar.progress(min(0.9, progress_bar.value + 0.05))
#             except queue.Empty: pass
#             try:
#                 while line := stderr_queue.get_nowait():
#                     stderr_text += line
#                     stderr_container.text_area("Error Output", stderr_text, height=100)
#             except queue.Empty: pass
#             time.sleep(0.1)

#         if process.returncode == 0:
#             progress_placeholder.success("CRISPResso2 completed successfully!")
#             progress_bar.progress(1.0)
#             output_dir = extract_output_dir(stdout_text) or st.session_state.crispresso.get('output_dir')
#             if output_dir:
#                 st.session_state.crispresso_results["path"] = output_dir
#                 st.info("Switch to the 'Results' tab to view the report.")
#         else:
#             progress_placeholder.error(f"CRISPResso2 failed with return code {process.returncode}")

# def render_crispresso_results_tab():
#     st.subheader("CRISPResso2 Results")
#     if not st.session_state.crispresso_results.get("path"):
#         st.info("No CRISPResso2 results available. Run an analysis first.")
#         return

#     results_path = st.session_state.crispresso_results["path"]
#     if not os.path.isdir(results_path):
#         st.error(f"Results directory not found: {results_path}")
#         return

#     st.success(f"Displaying results from: `{results_path}`")
#     report_files = [f for f in os.listdir(results_path) if f.endswith('_report.html')]
#     if report_files:
#         report_path = os.path.join(results_path, report_files[0])
#         # This link won't work in a sandboxed environment but is correct for local execution.
#         st.markdown(f'<a href="file://{os.path.abspath(report_path)}" target="_blank">View Full CRISPResso2 HTML Report</a>', unsafe_allow_html=True)

#     result_tabs = st.tabs(["Summary", "Alleles", "Quantification", "Plots"])
#     with result_tabs[0]: render_summary_results(results_path)
#     with result_tabs[1]: render_dataframe_results(results_path, "Alleles_frequency_table.txt", "Alleles")
#     with result_tabs[2]: render_dataframe_results(results_path, "Quantification_of_editing_frequency.txt", "Quantification")
#     with result_tabs[3]: render_plot_results(results_path)

# def render_summary_results(results_path):
#     st.subheader("Analysis Summary")
#     info_file = os.path.join(results_path, "CRISPResso_info.json")
#     if not os.path.isfile(info_file):
#         st.warning("CRISPResso_info.json not found.")
#         return

#     with open(info_file, 'r') as f:
#         info = json.load(f)

#     col1, col2 = st.columns(2)
#     with col1:
#         st.metric("Total Reads", info.get("total_reads", "N/A"))
#         st.metric("Reads Aligned", info.get("n_reads_aligned", "N/A"))
#     with col2:
#         st.metric("Reads with Modifications", info.get("n_reads_modified", "N/A"))
#         if info.get("total_reads", 0) > 0:
#             mod_rate = (info.get("n_reads_modified", 0) / info["total_reads"]) * 100
#             st.metric("Modification Rate", f"{mod_rate:.2f}%")

#     if "crispresso_summary" not in st.session_state or st.button("Regenerate AI Summary"):
#         with st.spinner("Generating analysis summary with AI..."):
#             summary_data = {k: info.get(k, "N/A") for k in ["total_reads", "n_reads_aligned", "n_reads_modified"]}
#             st.session_state.crispresso_summary = explain_results_with_llm("summary metrics", summary_data)
#     st.markdown("### AI Analysis")
#     st.markdown(st.session_state.crispresso_summary)

# def render_dataframe_results(results_path, filename, title):
#     st.subheader(f"{title} Frequency")
#     file_path = os.path.join(results_path, filename)
#     if not os.path.isfile(file_path):
#         st.warning(f"{filename} not found.")
#         return
#     try:
#         df = pd.read_csv(file_path, sep='\t')
#         st.dataframe(df)
#         st.download_button(f"Download {title} CSV", df.to_csv(index=False).encode('utf-8'), f"{title.lower()}.csv", "text/csv")
#         if st.button(f"Explain {title} Data"):
#             with st.spinner(f"Analyzing {title} data..."):
#                 st.info(explain_results_with_llm(f"{title.lower()} frequency", df.head(5).to_string()))
#     except Exception as e:
#         st.error(f"Error reading {filename}: {e}")

# def render_plot_results(results_path):
#     st.subheader("Result Plots")
#     plot_files = [f for f in os.listdir(results_path) if f.endswith('.png')]
#     if not plot_files:
#         st.warning("No plot images (.png) found in the results directory.")
#         return
#     for plot_file in plot_files:
#         st.image(os.path.join(results_path, plot_file), caption=plot_file)


# if __name__ == "__main__":
#     crispresso_page() 