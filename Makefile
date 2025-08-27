BACKEND=oncology-coPilot/oncology-backend-minimal

.PHONY: backend chromatin enformer borzoi all

backend:
	cd $(BACKEND) && source venv/bin/activate && PYTHONPATH=. uvicorn api.main:app --host 127.0.0.1 --port 8000 --log-level info

enformer:
	uvicorn tools.chromatin.enformer_server:app --host 127.0.0.1 --port 7001 --log-level warning

borzoi:
	uvicorn tools.chromatin.borzoi_server:app --host 127.0.0.1 --port 7002 --log-level warning

chromatin:
	@echo "Starting Enformer and Borzoi proxies..."
	@$(MAKE) -j2 enformer borzoi

all:
	@echo "Start backend and chromatin proxies in separate terminals:"
	@echo "  make backend"
	@echo "  make chromatin"

.PHONY: bench_variant
bench_variant:
	python3 tools/benchmarks/variant_auroc.py --download --n_pos 200 --n_neg 200 --allow_grch37 --api_base http://127.0.0.1:8000 --model_id evo2_7b --out tools/benchmarks/variant_auroc_results.json


