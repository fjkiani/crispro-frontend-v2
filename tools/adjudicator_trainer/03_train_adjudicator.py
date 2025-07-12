import pandas as pd
import numpy as np
import json
import joblib
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import argparse
import os

# --- Configuration ---
INPUT_FILE = "data/adjudicator_training/labeled_embeddings.parquet"
MODEL_OUTPUT_DIR = "models/adjudicator"
MODEL_OUTPUT_FILE = os.path.join(MODEL_OUTPUT_DIR, "adjudicator_v1.pkl")

# --- Main Logic ---

def load_and_prepare_data(file_path: str):
    """
    Loads the Parquet file and prepares the data for training.
    - Extracts embeddings from JSON strings.
    - Separates features (X) and labels (y).
    """
    print(f"🔄 Loading and preparing data from {file_path}...")
    df = pd.read_parquet(file_path)

    # Filter out rows where the embedding is missing
    df = df.dropna(subset=['embedding'])
    print(f"  -> Found {len(df)} variants with embeddings.")

    # The 'embedding' column is a JSON string, convert it back to a list of floats
    # This is a bit slow, but necessary.
    X = np.array(df['embedding'].apply(json.loads).tolist())
    
    # The 'significance_mapped' column contains our 0s (Benign) and 1s (Pathogenic)
    y = df['significance_mapped'].values

    print(f"  -> Data prepared. Feature matrix shape: {X.shape}, Label vector shape: {y.shape}")
    
    # Report class balance
    class_counts = df['significance_mapped'].value_counts()
    print(f"  -> Class balance: {class_counts[1]} Pathogenic (1) vs. {class_counts[0]} Benign (0)")
    
    return X, y

def train_classifier(X_train, y_train):
    """
    Trains a Gradient Boosting Classifier. This model is a good balance of performance and speed.
    """
    print("🤖 Training the Adjudicator model (GradientBoostingClassifier)...")
    
    # Initialize the classifier with some robust starting parameters
    # n_estimators=100 is often a good starting point.
    # We use a validation set to enable early stopping, which prevents overfitting.
    clf = GradientBoostingClassifier(
        n_estimators=200, 
        learning_rate=0.1,
        max_depth=5,
        random_state=42,
        verbose=1,
        validation_fraction=0.1, # Use 10% of training data for early stopping
        n_iter_no_change=10 # Stop if validation score doesn't improve for 10 rounds
    )
    
    clf.fit(X_train, y_train)
    print("  -> Model training complete.")
    return clf

def evaluate_model(clf, X_test, y_test):
    """
    Evaluates the trained model on the test set and prints a report.
    """
    print("\n📊 Evaluating model performance...")
    y_pred = clf.predict(X_test)
    
    accuracy = accuracy_score(y_test, y_pred)
    print(f"  -> Accuracy: {accuracy:.4f}")
    
    print("  -> Classification Report:")
    # Use target names for a more readable report
    report = classification_report(y_test, y_pred, target_names=["Benign (0)", "Pathogenic (1)"])
    print(report)
    
    print("  -> Confusion Matrix:")
    print(confusion_matrix(y_test, y_pred))

def save_model(clf, output_path: str):
    """
    Saves the trained model to a file using joblib for efficient serialization.
    """
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📂 Created model output directory: {output_dir}")

    print(f"\n💾 Saving trained model to {output_path}...")
    joblib.dump(clf, output_path)
    print(f"  -> Model saved successfully.")

def main():
    parser = argparse.ArgumentParser(description="Train the Adjudicator model on variant embeddings.")
    args = parser.parse_args()

    print("\n--- 💥 Operation Adjudicator: Phase 2.1 - Train Classifier 💥 ---")
    
    if not os.path.exists(INPUT_FILE):
        print(f"FATAL: Input data file not found at {INPUT_FILE}")
        return

    # 1. Load and prepare data
    X, y = load_and_prepare_data(INPUT_FILE)

    # 2. Split data into training and testing sets (80/20 split)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    print(f"\nSplitting data into {len(X_train)} training and {len(X_test)} testing samples.")

    # 3. Train the model
    classifier = train_classifier(X_train, y_train)

    # 4. Evaluate the model
    evaluate_model(classifier, X_test, y_test)
    
    # 5. Save the final model
    save_model(classifier, MODEL_OUTPUT_FILE)

    print("\n--- 🔥 Adjudicator Training Complete. The model is ready for deployment. 🔥 ---")

if __name__ == "__main__":
    main() 