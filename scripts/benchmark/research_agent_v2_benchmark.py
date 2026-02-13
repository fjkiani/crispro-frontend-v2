import sys
import os
import asyncio
import logging
import time
import json

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../oncology-coPilot/oncology-backend')))

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger("BenchmarkV2")

async def run_benchmark():
    logger.info("==================================================")
    logger.info("   RESEARCH PORTAL AGENT BENCHMARK - ITERATION 1")
    logger.info("==================================================\n")

    try:
        # Import the NEW ResearchAgent
        from backend.agents.research.research_agent import ResearchAgent
        agent = ResearchAgent()
        logger.info("✅ ResearchAgent Instantiated.\n")
    except Exception as e:
        logger.error(f"❌ Failed to instantiate ResearchAgent: {e}")
        return

    # Using the same queries as baseline, but expecting improvement on Q1
    queries = [
        {
            "id": "Q1_BROAD_MECHANISM",
            "prompt": "What are the common resistance mechanisms to Cisplatin in Ovarian Cancer?",
            "intent": "research_general", # Updated intent for general research
            "patient_id": "BENCHMARK_PATIENT_001",
            "patient_data": {"patientId": "BENCHMARK_PATIENT_001", "mutations": []}
        },
        {
            "id": "Q2_SPECIFIC_VARIANT",
            "prompt": "Analyze impact of localized KRAS G12C mutation on treatment options.",
            "intent": "analyze_genomic_criterion",
            "patient_id": "BENCHMARK_PATIENT_002",
            "patient_data": {
                "patientId": "BENCHMARK_PATIENT_002", 
                "mutations": [{"hugo_gene_symbol": "KRAS", "protein_change": "G12C", "variant_type": "Missense_Mutation"}]
            }
        },
        # Q3 omitted for brevity as it passed baseline
    ]

    results_log = []

    for q in queries:
        logger.info(f"--- Running Query: {q['id']} ---")
        logger.info(f"Prompt: {q['prompt']}")
        
        start_time = time.time()
        try:
            prompt_details = {
                "prompt": q["prompt"],
                "intent": q["intent"],
                "criterion_id": q["id"]
            }
            
            # Using the simplified run() signature of our ResearchAgent wrapper
            # It expects (patient_data, prompt_details)
            response = await agent.run(q["patient_data"], prompt_details)
            end_time = time.time()
            duration = end_time - start_time
            
            status = response.get("status", "UNKNOWN")
            evidence_len = len(str(response.get("evidence", "")))
            source = response.get("source", "Unknown")
            
            logger.info(f"✅ Completed in {duration:.2f}s")
            logger.info(f"Status: {status}")
            logger.info(f"Source: {source}")
            
            success = status != "ERROR" and status != "UNCLEAR"
            
            results_log.append({
                "query_id": q["id"],
                "duration_seconds": duration,
                "status": status,
                "success": success,
                "source": source
            })
            
            snippet = str(response.get("evidence", ""))[:200].replace("\n", " ")
            logger.info(f"Snippet: {snippet}...\n")

        except Exception as e:
            end_time = time.time()
            logger.error(f"❌ Failed: {e}\n")
            results_log.append({
                "query_id": q["id"],
                "duration_seconds": end_time - start_time,
                "status": "EXCEPTION",
                "success": False
            })

    # Summary
    logger.info("==================================================")
    logger.info("             BENCHMARK V2 SUMMARY")
    logger.info("==================================================")
    for res in results_log:
        icon = "✅" if res["success"] else "❌"
        src = f"[{res['source']}]" if 'source' in res else ""
        logger.info(f"{icon} {res['query_id']}: {res['status']} {src} ({res['duration_seconds']:.2f}s)")
    logger.info("==================================================")

if __name__ == "__main__":
    asyncio.run(run_benchmark())
