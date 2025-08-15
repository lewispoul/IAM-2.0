.PHONY: env run test lint

# Detect platform python
PYTHON ?= python3

# Create or update local venv if not using conda
env:
	@if [ -f chem-env.yaml ]; then \
		echo "Conda environment detected. Use: conda env create -f chem-env.yaml && conda activate chem-env"; \
	else \
		$(PYTHON) -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt; \
	fi
	@echo "Run: make verify"

verify:
	@$(PYTHON) scripts/verify_env.py || true

run:
	uvicorn iam.backend.app:app --reload --port 8010

test:
	@echo "Run unit tests here (e.g., pytest -q)."

lint:
	@echo "Run linters here (e.g., ruff . || flake8 .)."
