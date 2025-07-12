import streamlit as st
import os
import json
import subprocess
import sys
import tempfile
from pathlib import Path
import time
import threading
import queue
import re
import shlex
import shutil
import pandas as pd
import importlib.util
from dotenv import load_dotenv

# --- Dynamic Imports ---

def import_llm_api():
    """Import the LLM API module from tools/llm_api.py"""
    try:
        spec = importlib.util.spec_from_file_location("llm_api", "tools/llm_api.py")
        llm_api = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(llm_api)
        return llm_api
    except Exception as e:
        st.error(f"Error importing LLM API: {str(e)}")
        return None

def import_educational_context():
    """Import the Educational Context module from tools/educational_context.py"""
    try:
        spec = importlib.util.spec_from_file_location("educational_context", "tools/educational_context.py")
        edu_context = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(edu_context)
        return edu_context
    except Exception as e:
        st.warning(f"Educational content not available: {str(e)}")
        return None

def import_guide_finder():
    """Import the intelligent guide finder module"""
    try:
        spec = importlib.util.spec_from_file_location("intelligent_guide_finder", "tools/intelligent_guide_finder.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception as e:
        st.error(f"Error importing Intelligent Guide Finder: {str(e)}")
        return None

def import_next_steps():
    """Import the next steps recommendation module"""
    try:
        spec = importlib.util.spec_from_file_location("next_steps", "tools/next_steps.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception as e:
        st.error(f"Error importing Next Steps Advisor: {str(e)}")
        return None


# --- LLM & API Calls ---

llm_api = import_llm_api()
edu_context_module = import_educational_context()
guide_finder_module = import_guide_finder()
next_steps_module = import_next_steps()


def ask_llm(prompt, provider="gemini"):
    """Call the LLM API with a prompt and return the response"""
    if llm_api is None:
        return "LLM API not available."
    try:
        response = llm_api.query_llm(prompt, provider=provider)
        return response
    except Exception as e:
        st.error(f"Error calling LLM: {str(e)}")
        return f"Error: {str(e)}"

# --- UI & Styling ---

def apply_custom_css():
    """Applies custom CSS styles to the Streamlit app."""
    custom_css = """
    <style>
        .main-header {
            font-size: 28px;
            font-weight: bold;
            color: #2c3e50;
            padding-bottom: 10px;
            border-bottom: 2px solid #3498db;
        }
        .subheader {
            font-size: 22px;
            font-weight: bold;
            color: #34495e;
            margin-top: 15px;
            margin-bottom: 10px;
        }
    </style>
    """
    st.markdown(custom_css, unsafe_allow_html=True)

def page_header(title, subtitle):
    """Creates a standardized page header."""
    st.markdown(f'<div class="main-header">{title}</div>', unsafe_allow_html=True)
    st.markdown(f"<p>{subtitle}</p>")


# --- Session State ---

def init_session_state():
    """Initialize Streamlit session state variables."""
    if 'page' not in st.session_state:
        st.session_state.page = 'Home'
    if 'chopchop_running' not in st.session_state:
        st.session_state.chopchop_running = False
    if 'chopchop_results' not in st.session_state:
        st.session_state.chopchop_results = None
    if 'crispresso_running' not in st.session_state:
        st.session_state.crispresso_running = False
    if 'crispresso_results' not in st.session_state:
        st.session_state.crispresso_results = None
    if 'intelligent_guides_running' not in st.session_state:
        st.session_state.intelligent_guides_running = False
    if 'intelligent_guides_results' not in st.session_state:
        st.session_state.intelligent_guides_results = None
    # This is the single source of truth for the agent's chat history
    if "messages" not in st.session_state:
        st.session_state.messages = [{"role": "assistant", "content": "How can I assist with your campaign?"}]
    # DEPRECATED: Remove old chat history if it exists to prevent conflicts
    if "chat_history" in st.session_state:
        del st.session_state["chat_history"]


# --- Educational Sidebar ---

def create_educational_sidebar():
    """
    Creates the educational sidebar with CRISPR terminology, concepts,
    and the Campaign Advisor.
    """
    if not edu_context_module:
        st.sidebar.warning("Educational content not available.")
        # Still render the agent even if edu content fails
        agent_chat_box()
        return

    with st.sidebar:
        if 'active_mutation' in st.session_state and st.session_state.active_mutation.get("hugo_gene_symbol"):
            hugo = st.session_state.active_mutation["hugo_gene_symbol"]
            prot_change = st.session_state.active_mutation["protein_change"]
            st.info(f"🎯 Actively Designing for: **{hugo} {prot_change}**")
        else:
            st.info("🎯 No active mutation target set.")

        st.session_state.show_therapeutic_context = st.checkbox(
            "Enable Therapeutic Development Context",
            value=st.session_state.get('show_therapeutic_context', False),
            help="Show LLM-generated context on therapeutic challenges.",
            key="therapeutic_context_toggle_sidebar"
        )
        st.markdown("---")

        st.header("CRISPR Education Center")
        ed_tabs = st.tabs(["Terminology", "Concepts"])

        with ed_tabs[0]: # Terminology
            _terminology_tab()
        with ed_tabs[1]: # Concepts
            _concepts_tab()

        st.markdown("---")
        # Integrate the primary agent here
        agent_chat_box()


def _terminology_tab():
    st.subheader("CRISPR Terminology")
    term_query = st.text_input("Search for a term:", key="term_search")
    if term_query:
        term_info = edu_context_module.get_term_definition(term_query, detailed=True)
        if "definition" in term_info:
            st.markdown(f"### {term_info.get('term', term_query)}")
            if "full_name" in term_info: st.markdown(f"*{term_info['full_name']}*")
            st.markdown(f"**Definition:** {term_info['definition']}")
            if "context" in term_info: st.markdown(f"**Context:** {term_info['context']}")
            if "related_terms" in term_info:
                st.markdown("**Related Terms:** " + ", ".join(term_info["related_terms"]))
    else:
        for term in ["CRISPR", "Cas9", "gRNA", "PAM", "Indel"]:
            with st.expander(term):
                st.write(edu_context_module.get_term_definition(term)["definition"])

def _concepts_tab():
    st.subheader("Advanced Concepts")
    concept_options = list(edu_context_module.ADVANCED_CONCEPTS.keys())
    selected_concept = st.selectbox("Select a concept:", concept_options)
    if selected_concept:
        concept_info = edu_context_module.get_advanced_concept(selected_concept)
        st.markdown(f"### {selected_concept}")
        st.markdown(f"**{concept_info['definition']}**")
        for key, value in concept_info.items():
            if key not in ['definition', 'name']:
                st.markdown(f"**{key.replace('_', ' ').capitalize()}:**")
                for item in value:
                    st.markdown(f"- {item if isinstance(item, str) else item['name']}")

# --- Agent Chat Box ---
def agent_chat_box():
    """Renders the context-aware agent chat box in the sidebar."""
    st.subheader("Campaign Advisor")

    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    if prompt := st.chat_input("Ask your advisor...", key="campaign_advisor_input"):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

        with st.chat_message("assistant"):
            with st.spinner("Thinking..."):
                context = "You are 'Zo', a helpful AI assistant and campaign planner for a CRISPR therapeutic design platform. Be concise and helpful. Use the following real-time data to inform your response:\n\n"
                
                if 'threat_assessment_results' in st.session_state and st.session_state.threat_assessment_results:
                    context += f"**Threat Assessment Dossier:** {json.dumps(st.session_state.threat_assessment_results)}\n"
                
                if 'guide_results' in st.session_state and st.session_state.guide_results:
                    context += f"**CRISPR Guide Design Results:** {json.dumps(st.session_state.guide_results)}\n"

                if 'nanobody_result' in st.session_state and st.session_state.nanobody_result:
                    context += f"**Nanobody Design Results:** {json.dumps(st.session_state.nanobody_result)}\n"

                full_prompt = context + "\n\n**User Question:** " + prompt
                
                response = ask_llm(full_prompt, "anthropic")
                st.markdown(response)
        
        st.session_state.messages.append({"role": "assistant", "content": response})


# --- Subprocess & Tooling ---
def stream_subprocess_output(process, output_queue):
    """Stream subprocess output to a queue in a non-blocking way."""
    def reader(pipe, queue):
        try:
            with pipe:
                for line in iter(pipe.readline, ''):
                    queue.put(line)
        finally:
            queue.put(None)

    threading.Thread(target=reader, args=[process.stdout, output_queue], daemon=True).start()


def guide_finder_module():
    """Placeholder for a complex module with its own UI and logic."""
    st.info("Guide Finder Module (Placeholder)")

# --- Educational Context ---
class EducationalContext:
    def __init__(self, file_path="assets/educational_content.json"):
        with open(file_path, "r") as f:
            self.data = json.load(f)

    def get_llm_tooltip(self, content_key):
        return self.data.get(content_key, {}).get("tooltip", "Tooltip not found.")