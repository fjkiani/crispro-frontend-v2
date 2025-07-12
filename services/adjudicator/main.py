import modal
import joblib
import numpy as np
import os

# --- Image Definition ---
# A lightweight image with scikit-learn and other dependencies.
adjudicator_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("scikit-learn", "pandas", "joblib", "numpy", "fastapi")
    .copy_local_file("models/adjudicator/adjudicator_v1.pkl", "/root/adjudicator_v1.pkl")
)

# --- App Definition ---
app = modal.App("adjudicator-service", image=adjudicator_image)

# --- Model Loading & Service Class ---
@app.cls(cpu=1)
class Adjudicator:
    @modal.enter()
    def load_model(self):
        """Load the trained model into memory when the container starts."""
        self.model_path = "/root/adjudicator_v1.pkl"
        if not os.path.exists(self.model_path):
            raise FileNotFoundError(f"Model file not found at {self.model_path}")
        
        self.model = joblib.load(self.model_path)
        print("✅ Adjudicator model loaded successfully.")

    @modal.asgi_app()
    def api(self):
        """A simple FastAPI app to serve the classifier."""
        from fastapi import FastAPI, Body
        from fastapi.responses import JSONResponse

        fastapi_app = FastAPI()

        @fastapi_app.post("/classify")
        def classify_embedding(embedding: list[float] = Body(...)):
            """
            Accepts an embedding vector and returns a pathogenicity classification.
            """
            try:
                if not isinstance(embedding, list) or not all(isinstance(x, (int, float)) for x in embedding):
                    return JSONResponse(content={"error": "Invalid embedding format. Must be a list of numbers."}, status_code=400)

                # The model expects a 2D array, so we reshape the input vector
                embedding_array = np.array(embedding).reshape(1, -1)

                # Get the prediction (0 or 1)
                prediction = self.model.predict(embedding_array)[0]
                
                # Get the probability scores
                probabilities = self.model.predict_proba(embedding_array)[0]
                
                pathogenic_prob = float(probabilities[1])
                interpretation = "Pathogenic" if prediction == 1 else "Benign"

                response = {
                    "prediction": interpretation,
                    "is_pathogenic": bool(prediction == 1),
                    "confidence_score": pathogenic_prob,
                    "status": "success"
                }
                return JSONResponse(content=response, status_code=200)

            except Exception as e:
                return JSONResponse(
                    content={"status": "error", "message": f"Internal server error: {str(e)}"},
                    status_code=500
                )
        
        return fastapi_app

@app.local_entrypoint()
def local_test():
    """A local test function for the Adjudicator."""
    # This is a placeholder for a real embedding.
    # In a real test, you'd get this from the Zeta Oracle.
    mock_embedding = [0.0] * 8192
    
    adjudicator = Adjudicator()
    adjudicator.load_model() # Manually call for local testing

    # Simulate the FastAPI call logic
    embedding_array = np.array(mock_embedding).reshape(1, -1)
    prediction = adjudicator.model.predict(embedding_array)[0]
    probabilities = adjudicator.model.predict_proba(embedding_array)[0]
    
    print("--- 🔬 LOCAL ADJUDICATOR TEST 🔬 ---")
    print(f"  -> Prediction: {'Pathogenic' if prediction == 1 else 'Benign'}")
    print(f"  -> Confidence (Pathogenic): {probabilities[1]:.4f}")
    print("---------------------------------") 