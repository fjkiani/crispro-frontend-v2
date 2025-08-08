import pandas as pd
import numpy as np
import json
import joblib
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import os

# --- Configuration ---
INPUT_FILE = "data/adjudicator_training/labeled_embeddings.parquet"
MODEL_FILE = "models/adjudicator/adjudicator_v1.pkl"
OUTPUT_DIR = "results/adjudicator"
CHART_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "adjudicator_v1_confusion_matrix.png")

# --- Main Logic ---

def load_and_prepare_data(file_path: str):
    """
    Loads the Parquet file and prepares the data for evaluation.
    This is identical to the preparation step in the training script.
    """
    df = pd.read_parquet(file_path)
    df = df.dropna(subset=['embedding'])
    
    # The 'embedding' column is a JSON string, convert it back to a list of floats.
    embeddings = df['embedding'].apply(json.loads).tolist()
    
    X = np.array(embeddings)
    y = df['significance_mapped'].values
    
    return X, y

def generate_confusion_matrix_chart(y_true, y_pred, output_path: str):
    """
    Generates and saves a confusion matrix heatmap.
    """
    cm = confusion_matrix(y_true, y_pred)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Benign', 'Pathogenic'], 
                yticklabels=['Benign', 'Pathogenic'])
    
    plt.title('Adjudicator Model v1 - Confusion Matrix', fontsize=16)
    plt.ylabel('Actual Label', fontsize=12)
    plt.xlabel('Predicted Label', fontsize=12)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    plt.savefig(output_path)
    print(f"✅ Confusion matrix chart saved to {output_path}")

def main():
    """
    Main function to load the model, evaluate it on the test set,
    and generate a confusion matrix visualization.
    """
    print("--- 📊 Generating Adjudicator Performance Visualization 📊 ---")

    # Load data
    if not os.path.exists(INPUT_FILE):
        print(f"FATAL: Input data file not found at {INPUT_FILE}")
        return
    X, y = load_and_prepare_data(INPUT_FILE)
    
    # Split data exactly as in training to get the same test set
    _, X_test, _, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    
    # Load model
    if not os.path.exists(MODEL_FILE):
        print(f"FATAL: Model file not found at {MODEL_FILE}")
        return
    model = joblib.load(MODEL_FILE)
    print(f"Loaded model from {MODEL_FILE}")

    # Make predictions
    y_pred = model.predict(X_test)
    
    # Generate and save chart
    generate_confusion_matrix_chart(y_test, y_pred, CHART_OUTPUT_FILE)
    
    print("--- ✅ Visualization Complete ---")

if __name__ == "__main__":
    main() 