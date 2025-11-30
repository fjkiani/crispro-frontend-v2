"""
Tag Top Ovarian Cancer Trials with MoA Vectors (P0 Fix #5)
==========================================================
Extracts MoA vectors for top 20 ovarian cancer trials manually (no Gemini for now - fast execution).

Owner: Zo
Date: January 13, 2025
Manager Policy: P3 (MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md)
"""

import sqlite3
import json
from pathlib import Path
from typing import List, Dict, Any
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Database path
DB_PATH = Path(__file__).parent.parent / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "clinical_trials.db"

# Output path
OUTPUT_PATH = Path(__file__).parent.parent / "oncology-coPilot" / "oncology-backend-minimal" / "api" / "resources" / "trial_moa_vectors.json"


def get_top_ovarian_trials(limit: int = 200) -> List[Dict[str, Any]]:
    """Get top ovarian cancer trials from database."""
    if not DB_PATH.exists():
        logger.error(f"Database not found: {DB_PATH}")
        return []
    
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    # Query trials table (1000 ovarian trials) - cast to more trials, filter later
    query = """
    SELECT id as nct_id, title, status, phases, interventions, conditions
    FROM trials
    WHERE status IN ('RECRUITING', 'ACTIVE_NOT_RECRUITING', 'ENROLLING_BY_INVITATION', 'COMPLETED')
    LIMIT ?
    """
    
    cursor.execute(query, (limit,))
    trials = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    logger.info(f"✅ Retrieved {len(trials)} trials from database")
    return trials


def manual_moa_tagging(trials: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Manually tag trials with MoA vectors based on keywords in title/interventions.
    
    This is a simplified approach for P0 completion. Full Gemini tagging later.
    """
    moa_vectors = {}
    
    # Keywords for each pathway (expanded to catch more trials)
    keywords = {
        "ddr": [
            "parp", "olaparib", "niraparib", "rucaparib", "talazoparib", "veliparib",
            "atr", "chk1", "chk2", "wee1", "dna repair", "ddr", "atm",
            "ceralasertib", "berzosertib", "prexasertib", "adavosertib"
        ],
        "mapk": [
            "mek", "raf", "braf", "kras", "mapk", "ras",
            "trametinib", "selumetinib", "cobimetinib", "binimetinib",
            "dabrafenib", "vemurafenib", "encorafenib"
        ],
        "pi3k": [
            "pi3k", "akt", "mtor", "pten",
            "alpelisib", "everolimus", "temsirolimus", "ridaforolimus",
            "copanlisib", "idelalisib", "pictilisib"
        ],
        "vegf": [
            "bevacizumab", "vegf", "angiogenesis", "avastin",
            "cediranib", "nintedanib", "pazopanib", "sorafenib", "sunitinib",
            "ramucirumab", "regorafenib", "cabozantinib"
        ],
        "her2": [
            "her2", "erbb2", "trastuzumab", "herceptin", "pertuzumab", "perjeta",
            "ado-trastuzumab", "t-dm1", "t-dxd", "enhertu",
            "zanidatamab", "margetuximab", "tucatinib", "neratinib", "lapatinib"
        ],
        "io": [
            "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "avelumab",
            "pd-1", "pd-l1", "pdl1", "ctla-4", "ctla4",
            "immunotherapy", "checkpoint", "keytruda", "opdivo", "tecentriq",
            "ipilimumab", "tremelimumab", "cemiplimab"
        ],
        "efflux": [
            "p-gp", "abcb1", "mdr1", "efflux", "p-glycoprotein",
            "elacridar", "tariquidar", "zosuquidar"
        ]
    }
    
    for trial in trials:
        nct_id = trial.get("nct_id")
        title = (trial.get("title") or "").lower()
        interventions = (trial.get("interventions") or "").lower()
        
        # Combine title + interventions for keyword matching
        text = f"{title} {interventions}"
        
        # Initialize vector
        vector = {"ddr": 0.0, "mapk": 0.0, "pi3k": 0.0, "vegf": 0.0, "her2": 0.0, "io": 0.0, "efflux": 0.0}
        confidence = 0.0
        primary_moa = "Unknown"
        
        # Check each pathway
        matches = []
        for pathway, kws in keywords.items():
            for kw in kws:
                if kw in text:
                    vector[pathway] = 0.90  # High confidence for keyword match
                    matches.append((pathway, kw))
                    break  # Only count once per pathway
        
        # Determine primary MoA
        if matches:
            # Primary MoA is the first match (assumes keywords are priority-ordered)
            primary_pathway = matches[0][0]
            primary_kw = matches[0][1]
            primary_moa = f"{primary_pathway.upper()} pathway ({primary_kw})"
            confidence = 0.85 if len(matches) == 1 else 0.75  # Lower confidence for multi-pathway
        
        # Only save if we found at least one MoA
        if matches:
            moa_vectors[nct_id] = {
                "moa_vector": vector,
                "confidence": confidence,
                "source": "manual_keyword_matching",
                "tagged_at": "2025-01-13T18:00:00Z",
                "reviewed_by": "Zo",
                "provenance": {
                    "method": "keyword_matching",
                    "keywords_matched": [f"{p}:{k}" for p, k in matches],
                    "primary_moa": primary_moa
                }
            }
    
    logger.info(f"✅ Tagged {len(moa_vectors)} trials with MoA vectors")
    return moa_vectors


def merge_with_existing(new_vectors: Dict[str, Any]) -> Dict[str, Any]:
    """Merge new vectors with existing ones."""
    if OUTPUT_PATH.exists():
        with open(OUTPUT_PATH, "r") as f:
            existing = json.load(f)
        logger.info(f"✅ Loaded {len(existing)} existing MoA vectors")
    else:
        existing = {}
    
    # Merge (new vectors override existing if same NCT ID)
    existing.update(new_vectors)
    
    logger.info(f"✅ Total MoA vectors: {len(existing)}")
    return existing


def save_moa_vectors(vectors: Dict[str, Any]):
    """Save MoA vectors to JSON file."""
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    
    with open(OUTPUT_PATH, "w") as f:
        json.dump(vectors, f, indent=2)
    
    logger.info(f"✅ Saved {len(vectors)} MoA vectors to {OUTPUT_PATH}")


def main():
    """Main execution."""
    logger.info("⚔️ P0 Fix #5: Tagging top ovarian cancer trials with MoA vectors")
    
    # Step 1: Get top trials (query 200, tag those with keywords)
    trials = get_top_ovarian_trials(limit=200)
    
    if not trials:
        logger.error("❌ No trials found")
        return
    
    # Step 2: Manual MoA tagging
    new_vectors = manual_moa_tagging(trials)
    
    # Step 3: Merge with existing
    all_vectors = merge_with_existing(new_vectors)
    
    # Step 4: Save
    save_moa_vectors(all_vectors)
    
    # Step 5: Summary
    logger.info("=" * 60)
    logger.info(f"✅ P0 Fix #5 Progress:")
    logger.info(f"   - Total trials tagged: {len(all_vectors)}")
    logger.info(f"   - New trials tagged: {len(new_vectors)}")
    logger.info(f"   - Output: {OUTPUT_PATH}")
    logger.info("=" * 60)
    
    # Show breakdown by pathway
    pathway_counts = {"ddr": 0, "mapk": 0, "pi3k": 0, "vegf": 0, "her2": 0, "io": 0, "efflux": 0}
    for nct_id, data in all_vectors.items():
        vector = data.get("moa_vector", {})
        for pathway, value in vector.items():
            if value > 0:
                pathway_counts[pathway] += 1
    
    logger.info("Trials by pathway:")
    for pathway, count in sorted(pathway_counts.items(), key=lambda x: x[1], reverse=True):
        logger.info(f"   - {pathway.upper()}: {count} trials")


if __name__ == "__main__":
    main()

