import streamlit as st
import requests
import time
import json

COMMAND_CENTER_URL = "https://crispro--command-center-v6-orchestrator-web-app.modal.run"

st.set_page_config(
    page_title="Predator Protocol",
    page_icon="🐆",
    layout="wide"
)

st.title("🐆 Predator Protocol: Autonomous Therapeutic Design")
st.warning("This is a weapon of extreme power. Use with discipline and precision.", icon="⚠️")

target_gene = st.text_input(
    "**Enter Target Gene Symbol:**",
    "VEGFA",
    help="Enter the official gene symbol (e.g., 'MMP9', 'VEGFA', 'KRAS')."
)

if st.button("**UNLEASH PREDATOR**", use_container_width=True):
    if not target_gene:
        st.error("A target gene symbol is required to unleash the Predator.")
    else:
        with st.status(f"Unleashing Predator on **{target_gene}**...", expanded=True) as status:
            job_id = None
            try:
                # 1. Initiate the workflow
                status.write("🔥 Dispatching command to CommandCenter...")
                initiate_response = requests.post(
                    f"{COMMAND_CENTER_URL}/workflow/execute_predator_protocol",
                    json={"target_gene_symbol": target_gene},
                    timeout=30
                )
                initiate_response.raise_for_status()
                job_id = initiate_response.json().get("job_id")
                st.info(f"CommandCenter has acknowledged the order. **Job ID: `{job_id}`**")

                # 2. Poll for results
                while True:
                    time.sleep(5)
                    status_response = requests.get(f"{COMMAND_CENTER_URL}/status/{job_id}", timeout=30)
                    status_response.raise_for_status()
                    job_data = status_response.json()
                    
                    current_status = job_data.get("status")
                    message = job_data.get("message")
                    
                    status.update(label=f"Predator campaign against **{target_gene}**: {current_status.upper()}")
                    
                    if current_status == "running":
                        st.write(f"**Log:** {message}")
                    elif current_status == "complete":
                        st.success(f"**MISSION COMPLETE:** {message}")
                        st.subheader("🏆 Final Weapon Blueprints 🏆")
                        st.json(job_data.get("result"))
                        break
                    elif current_status == "failed":
                        st.error(f"**MISSION FAILED:** {message}")
                        break

            except requests.exceptions.RequestException as e:
                st.error(f"A communication error occurred with the CommandCenter: {e}")
            except Exception as e:
                st.error(f"An unexpected error occurred: {e}") 