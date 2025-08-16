.PHONY: help env env-update test test-cov lint format type-check security clean docker-build docker-dev docker-clean install-dev run dev docs

# Configuration
PYTHON := python3
CONDA_ENV := iam2
DOCKER_TAG := iam2:latest
DEV_COMPOSE := docker-compose.dev.yml

help: ## Show this help message
	@echo "IAM2.0 Migration - Available Commands:"
	@echo "======================================"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

# Environment Management
env: ## Create conda environment from chem-env.yaml
	@if [ -f chem-env.yaml ]; then \
		echo "Creating conda environment: $(CONDA_ENV)"; \
		conda env create -f chem-env.yaml; \
		echo "Environment created! Activate with: conda activate $(CONDA_ENV)"; \
	else \
		echo "Creating Python virtual environment..."; \
		$(PYTHON) -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt; \
		echo "Virtual environment created! Activate with: source .venv/bin/activate"; \
	fi

env-update: ## Update conda environment
	@echo "Updating conda environment: $(CONDA_ENV)"
	conda env update -f chem-env.yaml --prune
	@echo "Environment updated!"

verify: ## Verify environment setup
	@$(PYTHON) scripts/verify_env.py || true

install-dev: ## Install development dependencies
	@echo "Installing development dependencies..."
	pip install -r requirements.txt
	pip install -e . || true
	pre-commit install || echo "pre-commit not available"
	@echo "Development dependencies installed!"

# Testing
test: ## Run test suite
	@echo "Running test suite..."
	pytest tests/ -v

test-cov: ## Run tests with coverage report
	@echo "Running tests with coverage..."
	pytest tests/ -v --cov=backend --cov-report=html --cov-report=term-missing

test-integration: ## Run integration tests only
	@echo "Running integration tests..."
	pytest tests/ -v -m integration

test-smoke: ## Run smoke tests only
	@echo "Running smoke tests..."
	pytest tests/test_smoke.py -v

# Code Quality
lint: ## Run linting with ruff
	@echo "Running ruff linter..."
	ruff check backend/ tests/ --show-source --statistics

lint-fix: ## Run linting with auto-fix
	@echo "Running ruff with auto-fix..."
	ruff check backend/ tests/ --fix

format: ## Format code with black
	@echo "Formatting code with black..."
	black backend/ tests/

format-check: ## Check code formatting
	@echo "Checking code formatting..."
	black --check backend/ tests/

type-check: ## Run type checking with mypy
	@echo "Running type checker..."
	mypy backend/ --ignore-missing-imports

security: ## Run security checks
	@echo "Running security checks..."
	bandit -r backend/ -f json -o bandit-report.json || true
	safety check --json --output safety-report.json || true
	@echo "Security reports generated: bandit-report.json, safety-report.json"

pre-commit: ## Run all pre-commit hooks
	@echo "Running pre-commit hooks..."
	pre-commit run --all-files

# Application Management
run: ## Run the backend server
	@echo "Starting IAM2.0 backend server..."
	uvicorn iam.backend.app:app --host 0.0.0.0 --port 8000 --reload

dev: ## Run development server with auto-reload
	@echo "Starting development server..."
	./scripts/run_backend.sh || uvicorn iam.backend.app:app --host 0.0.0.0 --port 8000 --reload

worker: ## Run background worker
	@echo "Starting background worker..."
	cd backend && python -m backend.worker

# Docker Management
docker-build: ## Build Docker image
	@echo "Building Docker image: $(DOCKER_TAG)"
	docker build -t $(DOCKER_TAG) .

docker-dev: ## Start development environment with Docker Compose
	@echo "Starting development environment..."
	docker-compose -f $(DEV_COMPOSE) up --build

docker-dev-bg: ## Start development environment in background
	@echo "Starting development environment in background..."
	docker-compose -f $(DEV_COMPOSE) up --build -d

docker-stop: ## Stop Docker Compose services
	@echo "Stopping development environment..."
	docker-compose -f $(DEV_COMPOSE) down

docker-clean: ## Clean Docker resources
	@echo "Cleaning Docker resources..."
	docker-compose -f $(DEV_COMPOSE) down -v --remove-orphans
	docker system prune -f

# Utility Commands
clean: ## Clean temporary files and caches
	@echo "Cleaning temporary files..."
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + || true
	find . -type f -name ".coverage" -delete || true
	find . -type d -name "htmlcov" -exec rm -rf {} + || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + || true

status: ## Show project status
	@echo "IAM2.0 Migration Status:"
	@echo "========================"
	@echo "Environment: $(CONDA_ENV)"
	@echo "Python: $(shell which python || echo 'Not found')"
	@echo "Git branch: $(shell git branch --show-current 2>/dev/null || echo 'No git')"
	@echo ""

health: ## Check system health
	@echo "System Health Check:"
	@echo "==================="
	@echo "Python: $(shell python --version 2>/dev/null || echo 'Not found')"
	@echo "Conda: $(shell conda --version 2>/dev/null || echo 'Not found')"
	@echo "Docker: $(shell docker --version 2>/dev/null || echo 'Not found')"

# Quick commands for common workflows
setup: env install-dev ## Complete setup: create environment and install dependencies
	@echo "Setup complete! Next steps:"
	@echo "1. conda activate $(CONDA_ENV) (if using conda)"
	@echo "2. make test"
	@echo "3. make dev"

ci: lint format-check type-check test ## Run CI checks locally
	@echo "All CI checks passed!"
