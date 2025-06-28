# CrisPRO Intelligence Platform

The CrisPRO Intelligence Platform is a comprehensive, AI-powered application for advanced CRISPR research. It moves beyond simple guide design to de-risk and accelerate therapeutic R&D by using sophisticated AI models to intelligently select optimal guide RNA candidates.

The entire platform is orchestrated through a user-friendly Streamlit web application.

## The Challenge: Choosing the Right Guide

In CRISPR-based therapeutic development, choosing the single best guide RNA from hundreds of possibilities is a critical, high-stakes decision. A poor choice leads to failed experiments, wasted resources, and months of lost time. Traditional methods that rely on simple bioinformatics rules often fail to predict *in vivo* success.

## Our Solution: The Intelligent Guide Designer

The core of the platform is the **Intelligent Guide Designer**, which uses a novel AI methodology called **"Generate & Compare"** to predict guide RNA efficacy with unprecedented biological context.

### The "Generate & Compare" Method

Instead of relying on a static list of rules, we use a massive AI model to "think" like an expert biologist.

1.  **The "Expert" AI:** We use a powerful, 1-billion parameter generative DNA model that has learned the "grammar" of functional DNA from hundreds of millions of sequences. It understands the deep, contextual patterns of the genome.
2.  **GENERATE:** We provide a candidate guide sequence and ask the AI to "complete the sentence" by generating what it predicts is the most biologically plausible continuation.
3.  **COMPARE & SCORE:** We then perform a sequence alignment between the AI's ideal continuation and our original guide. A high alignment score indicates the AI found our guide to be "natural" and well-placed in a functional region, making it a high-quality candidate. A low score flags the guide as likely ineffective.

This AI-driven efficacy score is combined with a traditional off-target safety score (from BLAST) to provide a holistic, actionable ranking of all potential guides.

> ### **The Breakthrough: From Counting to Understanding**
>
> Our "Generate & Compare" method is a breakthrough because it's the first approach that approximates a true biological understanding of a guide's genomic context.
>
> 1.  **It Solves for Causation, Not Just Correlation:** The AI model hasn't been programmed with a list of "good guide rules." It has built an internal, implicit model of what *causes* a sequence to be functional. Its score is therefore a proxy for the actual biological plausibility of the guide's target site.
>
> 2.  **It Sees the Invisible Genomic Landscape:** A simple algorithm only sees the 20 letters of the guide. Our AI evaluates those 20 letters against its knowledge of millions of other sequences, implicitly informed by the grammar of promoters, enhancers, and other functional elements. It judges a guide based on how well it "fits in" with the functional landscape of the entire known biological world.
>
> 3.  **It Creates a Virtuous AI Flywheel:** A rule-based system is static. Our AI system learns. As we gather real-world experimental data on guide success, we can use it to fine-tune the model. The more data it gets, the smarter it becomes, and the more valuable the platform becomes. This creates a powerful competitive moat that is nearly impossible for competitors using static rules to overcome.
>
> In summary, we have moved from counting letters to understanding language. We've built an engine that doesn't just analyze a guide sequence; it *reasons* about its potential, providing a glimpse into its future success or failure before a single dollar is spent in the lab.

## How to Use the Platform

The platform is designed to be run as a web application from your local machine.

### 1. Prerequisites
- Python 3.9+
- An environment manager like `venv` or `conda` is recommended.
- Install dependencies:
  ```bash
  pip install -r requirements.txt
  ```

### 2. Configuration
The platform requires API endpoints for its AI models. These are configured in a `.env` file in the project root.

1.  Create a file named `.env`.
2.  Add the following lines, replacing the example URLs with your actual deployed Modal endpoints:
    ```env
    # Endpoint for the AI model that generates DNA sequences for scoring
    EVO2_GENERATIVE_ENDPOINT="https://your-modal-generative-endpoint.modal.run"

    # Endpoint for the AI model that predicts variant effects and drug response
    EVO2_DISCRIMINATIVE_ENDPOINT="https://your-modal-discriminative-endpoint.modal.run"
    ```

### 3. Running the Application
Launch the main Streamlit application with the following command:
```bash
streamlit run streamlit_app.py
```
This will open the CrisPRO platform in your web browser. Navigate to the **"Intelligent Guide Designer"** from the sidebar to use the core feature.

## Technology Stack
*   **Application Framework:** Streamlit
*   **Core Logic:** Python
*   **AI Model Serving:** Modal
*   **Key Libraries:** Biopython, Pandas, Scikit-learn, Requests

## Legacy Tools

The original command-line CHOPCHOP integration scripts still exist within the `tools/chopchop/` directory. However, their functionality has been superseded by the more powerful and user-friendly workflows available in the main Streamlit application.

## License

- This Project: See LICENSE file.
- CHOPCHOP: Apache License 2.0 (see `tools/chopchop/LICENSE`) 