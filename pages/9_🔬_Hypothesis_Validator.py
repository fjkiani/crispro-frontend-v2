import streamlit as st
import requests

def show_hypothesis_validator():
    """
    Renders the Hypothesis Validator page.
    This page will serve as the command center for OPERATION: DEEP DIVE.
    """
    st.set_page_config(
        page_title="Hypothesis Validator",
        page_icon="🔬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    st.title("🔬 Hypothesis Validator")
    st.markdown("### OPERATION: DEEP DIVE")
    st.info(
        """
        **Objective:** To bypass decades of failed research by using `in silico`
        methods to identify the active anti-angiogenic compound in shark cartilage,
        validate its mechanism against genomic databases, and design a modern
        CRISPR-based experiment to test the hypothesis.
        """
    )

    # --- Phase I: Target Identification ---
    with st.expander("Phase I: Target Identification & Reconnaissance", expanded=True):
        st.subheader("✅ Task 1.1: Identify Putative Active Compound")
        st.success("Target Acquired: **Neovastat (AE-941)**")


        st.subheader("✅ Task 1.2: Uncover Molecular Target")
        st.success("Molecular Target Identified: **VEGF/VEGFR-2 Signaling Axis**")


        st.subheader("Task 1.4: Manual Reconnaissance (Zeta-Stream Ingress Protocol)")
        st.write("Perform targeted web searches using the Tavily AI-powered search API.")
        query = st.text_input("Enter your research query:", "What is the mechanism of action for Neovastat (AE-941)?")

        if st.button("Execute Intelligence Run"):
            with st.spinner("Executing advanced reconnaissance via Tavily API..."):
                try:
                    response = requests.post(
                        "http://localhost:8000/api/intelligence/search",
                        json={"query": query}
                    )
                    response.raise_for_status()
                    data = response.json()
                    results = data.get("results", [])
                    answer = data.get("answer")

                    if answer:
                        st.subheader("Synthesized Answer")
                        st.markdown(answer)
                        st.divider()

                    if results:
                        st.subheader(f"Found {len(results)} Primary Sources:")
                        for result in results:
                            st.markdown(f"#### [{result['title']}]({result['url']})")
                            st.caption(f"Source: {result['url']} | Score: {result['score']:.2f}")
                            st.write(result['content'])
                            st.divider()
                    elif not answer:
                        st.warning("No results found. The target may be obscure or the search terms too narrow.")

                except requests.exceptions.RequestException as e:
                    st.error(f"Failed to connect to intelligence backend: {e}")


    # --- Phase II: In Silico Target Validation (Placeholder) ---
    with st.expander("Phase II: In Silico Target Validation", expanded=False):
        st.write("Awaiting implementation...")

    # --- Phase III: Design CRISPR Kill Vehicle (Placeholder) ---
    with st.expander("Phase III: Design CRISPR-Based Kill Vehicle", expanded=False):
        st.write("Awaiting implementation...")

show_hypothesis_validator() 