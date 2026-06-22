# Makefile -- reproducibility entry points for BRCA-PathwayML
#
# Editor requirement E3 (BMC Artificial Intelligence revision): one-command
# reproduction of the headline results from committed weights/intermediates.
#
# Use the pinned environment described in requirements.txt. The `python` below
# should resolve to a Python 3.12 interpreter with `pip install -r
# requirements.txt` already applied (e.g. an activated venv). The Dockerfile in
# this repo builds that environment (including R + genefu) if you prefer.
#
# Two reproducibility tiers (see data/processed/REGENERATION.md):
#   Tier 1 (default): reproduce from COMMITTED intermediates -- no raw download,
#                     no R required for the modeling steps.
#   Tier 2 (all-from-raw): rebuild everything from raw expression. Needs the
#                     multi-GB raw inputs and R/genefu; partly notebook-based.

PYTHON ?= python3

.DEFAULT_GOAL := help

.PHONY: help reproduce-headline smoke test lint consistency phase3 figures \
        all all-from-raw migrate-pickles fetch

help: ## Show this help
	@echo "BRCA-PathwayML -- make targets:"
	@echo ""
	@echo "  reproduce-headline  Reproduce headline tables + TNBC endpoint, then verify"
	@echo "  smoke               Fast recompute smoke + primary-endpoint tests"
	@echo "  test                Full pytest suite (PAM50-ROR test needs R/genefu)"
	@echo "  lint                ruff check + mypy (mypy non-blocking)"
	@echo "  consistency         Recompute-vs-committed submission consistency gate"
	@echo "  phase3              Run the four Phase-3 exploratory analyses"
	@echo "  figures             Build manuscript figures/tables"
	@echo "  all                 From committed intermediates: train -> external -> phase3 -> figures -> consistency"
	@echo "  all-from-raw        Tier-2 full rebuild from raw data (needs raw inputs + R; see REGENERATION.md)"
	@echo "  migrate-pickles     Make committed pickles load cold (no __main__ shim)"
	@echo ""
	@echo "Set PYTHON=... to point at your pinned 3.12 interpreter (default: python3)."

reproduce-headline: ## Reproduce headline external-validation tables + endpoint, then verify
	$(PYTHON) scripts/external_validation.py
	$(PYTHON) scripts/check_submission_consistency.py

smoke: ## Fast recompute smoke test + primary-endpoint test
	$(PYTHON) -m pytest tests/test_smoke.py tests/test_primary_endpoint.py -q

test: ## Full test suite (test_pam50_official.py auto-skips without R/genefu)
	$(PYTHON) -m pytest tests/ -q

lint: ## Lint with ruff; type-check with mypy (mypy is non-blocking)
	ruff check .
	mypy src scripts || true

consistency: ## Recompute-vs-committed submission consistency gate (exit 0 = all pass)
	$(PYTHON) scripts/check_submission_consistency.py

phase3: ## Phase-3 exploratory analyses (read external-validation results)
	$(PYTHON) scripts/within_subtype_analysis.py
	$(PYTHON) scripts/calibration_analysis.py
	$(PYTHON) scripts/dca_analysis.py
	$(PYTHON) scripts/stability_analysis.py

figures: ## Build manuscript figures and tables
	$(PYTHON) scripts/build_table1.py
	$(PYTHON) scripts/build_cohort_flow.py
	$(PYTHON) scripts/build_forest_plot.py
	$(PYTHON) scripts/compute_model_attributions.py

all: ## Full Tier-1 chain from committed intermediates
	$(PYTHON) scripts/train_ml_zoo.py
	$(PYTHON) scripts/external_validation.py
	$(MAKE) phase3
	$(MAKE) figures
	$(PYTHON) scripts/check_submission_consistency.py

# Tier 2: rebuild from raw expression. REQUIRES the multi-GB raw inputs that are
# NOT committed (see data/processed/REGENERATION.md and data/raw/README.md) and
# an R install with genefu/GSVA/mice (see R/install.R).
#
# NOTE: the raw-expression -> expression-parquet parse step currently lives in
# notebooks/01_data_preprocessing.ipynb (run top-to-bottom) and is not yet a
# headless script, so this target is NOT fully turnkey -- run that notebook
# between `fetch` and the pipeline steps. Tier 1 (`make all`) is the supported
# review path.
fetch: ## Download raw expression + gene sets (Tier 2 only; see data/raw/README.md)
	cd data/raw && bash download_data.sh
	cd data/raw && bash download_gene_sets.sh

all-from-raw: fetch ## Tier-2 full rebuild from raw (see NOTE above re: notebook parse step)
	$(PYTHON) scripts/harmonize_cohorts.py
	$(PYTHON) scripts/score_pathways.py
	$(PYTHON) scripts/compute_clinical_baselines.py
	$(MAKE) all

migrate-pickles: ## Rewrite committed pickles so they load without a __main__ shim
	$(PYTHON) scripts/migrate_pickles.py
