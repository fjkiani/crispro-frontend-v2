# 📖 EXTRACTION PIECE 1.3: SAE vs Evo2 Clarification
**Source**: Lines 1077-1125 of `2025-11-19_09-12Z-designing-zero-tolerance-content-filtering-layer.md`  
**Date Extracted**: 2025-01-20  
**Status**: ✅ Complete

---

## 📋 SUMMARY

This section clarifies the critical distinction: SAE is NOT built into Evo2. It's a separate post-hoc interpretability model trained on Evo2 activations. This clarification came from direct analysis of the Evo2 paper.

---

## 🔍 KEY FINDINGS

### **What the Paper Says**

**Section 2.4 (lines 432-434):**
> "To probe what Evo 2 is capturing, we trained SAEs on its representations, or neuron firing patterns, to decompose the model into sparse, high-dimensional representations..."

**Section 4.4.1 (lines 1698-1713):**
> "We trained a BatchTopK sparse autoencoder... on activations from the Evo 2 residual stream following layer 26... Sparse autoencoders take as input model activations 𝑥 and autoencode them into a sparse feature vector 𝑓..."

---

### **What This Means**

1. **SAE is a separate model trained after Evo2**
   - SAE is NOT part of Evo2's architecture
   - SAE is a post-hoc interpretability tool
   - SAE was trained specifically to understand what Evo2 learned

2. **SAE takes Evo2 activations as input**
   - Input: Evo2 layer 26 activations (4096 dimensions)
   - Output: 32,768 interpretable features (8x overcomplete)
   - Process: Autoencodes activations into sparse feature vector

3. **SAE reveals interpretable features**
   - Exons, TF motifs, protein structure, prophage regions
   - Trained on 1B tokens
   - BatchTopK variant (k=64)

---

### **What We Currently Have**

From `SAE_UNDERSTANDING_AND_BIOMARKER_ROADMAP.md`:

> "CRITICAL GAP: We're NOT using actual Evo2 SAE activations. We're computing proxy SAE features from pathway aggregation and insights bundle."

**Current "SAE features":**
- DNA repair capacity: formula from pathway + insights (0.6/0.2/0.2)
- Mechanism vector: 7D from pathway aggregation
- Pathway burdens: from S/P/E, not SAE

**What's missing:**
- Actual Evo2 layer 26 activations
- Trained SAE model weights
- Real SAE feature extraction

---

### **What We Would Need for Real SAE Features**

1. **Extract layer 26 activations from Evo2**
   - The model supports `return_embeddings=True`
   - Need to specify `layer_names=["layer_26"]`
   - Output: 4096-dimensional activation vectors

2. **Load the trained SAE model weights**
   - 32K features, ~500MB–1GB file size
   - Hugging Face model: `Goodfire/Evo-2-Layer-26-Mixed`
   - BatchTopK SAE architecture

3. **Decode activations → SAE features**
   - Forward pass through SAE model
   - Input: [batch, seq_len, 4096] activations
   - Output: [batch, seq_len, 32768] SAE features
   - Extract top-k features per variant

---

### **Paper Resources**

The paper mentions:
- **Evo Mech Interp Visualizer**: https://arcinstitute.org/tools/evo/evo-mech-interp
- **Visualization tool**: For exploring SAE features alongside genomic data
- **Note**: Paper does NOT explicitly state that SAE model weights are released

**What we found:**
- SAE weights available on Hugging Face: `Goodfire/Evo-2-Layer-26-Mixed`
- Used in Phase 1 implementation (Modal SAE service)

---

## 📊 KEY INSIGHTS

### **Critical Distinction**

**SAE ≠ Built into Evo2**
- Evo2 is the foundational model (trained on genomic sequences)
- SAE is the interpretability layer (trained on Evo2's activations)
- SAE helps us understand what Evo2 learned, but it's separate

**Why This Matters:**
- We can't just "turn on" SAE in Evo2
- We need separate infrastructure to:
  1. Extract Evo2 activations
  2. Load SAE model weights
  3. Run SAE forward pass
  4. Extract top features

### **Current Implementation Gap**

**What we have:**
- Proxy SAE features (formulas from pathway/insights)
- Working system using these proxies
- Manager's policy based on proxy features

**What we're building:**
- Real SAE feature extraction (Phase 1 complete)
- Evo2 activations endpoint (Modal service)
- SAE service (Modal service with HF weights)
- Cohort extraction pipeline (in progress)

---

## 🔗 CONTEXT & CONNECTIONS

- **Builds on**: Piece 1.1 (Initial SAE Assessment), Piece 1.2 (Manager Policy)
- **Clarifies**: The fundamental relationship between SAE and Evo2
- **Leads to**: Understanding why Phase 1 implementation was needed
- **Key Insight**: SAE is a separate interpretability model, not part of Evo2

---

## 📝 NOTES

- This clarification was critical for understanding the architecture
- The paper explicitly states "we trained SAEs" - past tense, separate training
- SAE is a post-hoc tool, not built into Evo2's forward pass
- We need both Evo2 (for activations) and SAE (for feature extraction)
- The gap between proxy and real SAE features is now clear

---

## 🎯 QUESTIONS RESOLVED

- ✅ Is SAE built into Evo2? → **NO**, it's a separate model
- ✅ How does SAE relate to Evo2? → SAE is trained on Evo2 activations
- ✅ What do we need for real SAE? → Evo2 activations + SAE model weights
- ✅ Why are we using proxies? → Real SAE infrastructure wasn't built yet

