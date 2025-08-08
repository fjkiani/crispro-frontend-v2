import streamlit as st
import os
import importlib.util
import sys
import re

def main():
    """
    A dynamic Streamlit dashboard that discovers and loads all pages from the 'src/app' directory.
    """
    st.set_page_config(
        page_title="Crispr Assistant - Main Dashboard",
        page_icon="🚀",
        layout="wide"
    )

    st.sidebar.title("Master Control")

    # --- DYNAMIC PAGE DISCOVERY ---
    # The path to the directory containing the page files. Since this script
    # is now inside 'src/app', we can reference it directly.
    pages_dir = os.path.dirname(__file__)
    
    # Add the project root to the python path. The project root is two levels
    # up from the current file's location (src/app -> src -> root).
    project_root = os.path.abspath(os.path.join(pages_dir, '..', '..'))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    
    # Add the 'src' directory to the path as well.
    src_path = os.path.abspath(os.path.join(pages_dir, '..'))
    if src_path not in sys.path:
        sys.path.insert(0, src_path)

    try:
        # Get all .py files, ignoring hidden files and this dashboard itself.
        page_files = [f for f in os.listdir(pages_dir) if f.endswith(".py") and not f.startswith('.') and f != 'main_dashboard.py']
        page_files.sort()

        # Create a mapping from a "pretty" name to the file path
        page_mapping = {}
        for page_file in page_files:
            name_without_ext = os.path.splitext(page_file)[0]
            # Use regex to remove leading numbers and underscores for a cleaner name
            pretty_name = re.sub(r"^\d+[_ ]*", "", name_without_ext).replace('_', ' ')
            if not pretty_name: 
                pretty_name = name_without_ext
            page_mapping[pretty_name] = os.path.join(pages_dir, page_file)

        # Create the selection box in the sidebar
        selected_page_name = st.sidebar.selectbox(
            "Select Page",
            list(page_mapping.keys())
        )

        # --- PAGE RENDERING ---
        if selected_page_name:
            page_path = page_mapping[selected_page_name]
            
            st.title(f"{selected_page_name}")
            st.markdown("---")

            try:
                # Dynamically import and run the selected page's code
                module_name = f"pages.{selected_page_name.replace(' ', '_').replace('🧬', '').replace('⚔️', '')}"
                
                spec = importlib.util.spec_from_file_location(module_name, page_path)
                page_module = importlib.util.module_from_spec(spec)
                
                sys.modules[module_name] = page_module
                spec.loader.exec_module(page_module)

            except Exception as e:
                st.error(f"Error loading page '{selected_page_name}':")
                st.exception(e)

    except FileNotFoundError:
        st.error(f"The directory '{pages_dir}' was not found. Please ensure it exists.")
    except Exception as e:
        st.error("An unexpected error occurred:")
        st.exception(e)


if __name__ == "__main__":
    main() 