"""EvoService class for Evo2 model serving on Modal."""
import os
import sys
from fastapi import FastAPI
from loguru import logger

# Handle imports for both local and Modal deployment
try:
    from .patches import apply_torch_patches
    from .patches import apply_evo2_runtime_patches
except ImportError:
    # Modal deployment - direct imports
    from patches import apply_torch_patches
    from patches import apply_evo2_runtime_patches


class EvoService:
    """Evo2 service class for sequence scoring and generation."""
    
    def __init__(self, model_id: str = None):
        """Initialize service (model loading happens in load_model_and_api)."""
        self.model = None
        self.fastapi_app = None
        self._model_id = model_id
    
    @property
    def model_id(self) -> str:
        """Get the model ID to load, with fallback to environment variable or default."""
        if self._model_id:
            return self._model_id
        model_id = os.getenv('EVO_MODEL_ID', 'evo2_1b_base')
        # Set default if not present
        if not os.getenv('EVO_MODEL_ID'):
            os.environ['EVO_MODEL_ID'] = model_id
        return model_id
    
    def load_model_and_api(self, job_status_dict):
        """
        Load Evo2 model and initialize FastAPI app.
        
        This method should be called from @modal.enter() in the Modal app definition.
        
        Args:
            job_status_dict: Modal Dict for job status tracking (needed for endpoints)
        """
        # CRITICAL: Apply PyTorch patches BEFORE importing evo2
        # This ensures torch.load() works correctly with HuggingFace checkpoints
        apply_torch_patches()
        
        # Now safe to import evo2
        from evo2 import Evo2
        # Patch Evo2 scoring runtime (handles upstream signature mismatches)
        apply_evo2_runtime_patches()
        
        model_id = self.model_id
        logger.info(f"🚀 Loading Evo2 model: {model_id} ...")
        self.model = Evo2(model_id)
        logger.info("🎉 Evo2 model loaded successfully!")
        
        # Initialize FastAPI app
        self.fastapi_app = FastAPI(title="Evo2 General-Purpose Generation Service")
        
        # Register all endpoints (import here to avoid circular dependency)
        try:
            from .endpoints import register_endpoints
        except ImportError:
            # Modal deployment - direct import
            from endpoints import register_endpoints
        register_endpoints(self.fastapi_app, self, job_status_dict)
        
        return self.fastapi_app
