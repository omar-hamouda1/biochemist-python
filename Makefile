# ─── Biochemist Python — Makefile ─────────────────────────────────────
# Shortcuts for common development tasks
#
# Usage:
#   make setup      — Create conda environment
#   make test       — Run all tests
#   make lint       — Check code style
#   make clean      — Remove generated files
#   make merge      — Merge docking + ADMET results
# ──────────────────────────────────────────────────────────────────────

.PHONY: setup test lint clean merge help

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

setup: ## Create conda environment from environment.yml
	conda env create -f environment.yml
	@echo "✅ Environment created. Activate with: conda activate biochem"

test: ## Run all unit tests with pytest
	pytest tests/ -v --tb=short

lint: ## Check code style with flake8
	flake8 src/ scripts/ --max-line-length=100 --ignore=E501,W503

clean: ## Remove generated/temporary files
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} + 2>/dev/null || true
	rm -f temp_*.pdb 2>/dev/null || true
	@echo "✅ Cleaned up generated files"

merge: ## Merge all pipeline results into final summary
	python scripts/merge_final_results.py

lab: ## Launch JupyterLab
	jupyter lab --notebook-dir=notebooks
