"""
Vercel-compatible entry point for the Command Center Backend
This file adapts the main Command Center FastAPI app for Vercel's serverless environment
"""

import sys
import os
from pathlib import Path

# Add the src directory to Python path for imports
current_dir = Path(__file__).parent
src_dir = current_dir.parent.parent.parent  # Go up to src/
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

# Import the FastAPI app from the command center
from services.command_center.main import fastapi_app

# Vercel expects the app to be named 'app'
app = fastapi_app

# Export the app for Vercel
# Vercel will look for a variable named 'app' in this file 