# BRCA-PathwayML

[![CI](https://github.com/suhaanthayyil/BioInformed-BRCA-Risk/actions/workflows/ci.yml/badge.svg)](https://github.com/suhaanthayyil/BioInformed-BRCA-Risk/actions/workflows/ci.yml)

Interpretable pathway-based machine learning for breast cancer prognosis: a four-cohort benchmark against PAM50-ROR.

**Authors:** Suhaan Thayyil, Eshaan Nidee

## Overview

This repository contains the code, analysis pipelines, and results for a four-cohort study comparing a pathway-based survival machine learning model (Gradient Boosted Survival Analysis) against official genefu PAM50-ROR for breast cancer prognosis.

The pre-registered primary endpoint (TNBC meta-analyzed delta C-index >= +0.03, p < 0.05) was **not met**. The observed result was delta = +0.0144, 95% CI [-0.0472, +0.0760], p = 0.6466. All secondary analyses are exploratory.

## Key Results

- **Harmonized database:** 4,532 samples across four cohorts (TCGA-BRCA, GSE96058/SCAN-B, METABRIC, GSE20685).
- **Analysis sets:** 4,003 patients in patient-characteristics summaries; 4,002 metric-evaluable patients in the official PAM50-ROR discrimination analysis.
- **Headline model:** Gradient Boosted Survival (internal TCGA CV C-index = 0.642, delta vs Cox = +0.042).
- **TNBC primary endpoint:** NOT MET (meta delta = +0.0144, p = 0.6466).
- **Overall discrimination:** comparable to PAM50-ROR across external cohorts.

## Reproducibility traceability

Helpers for verifying the central reported values:

- `scripts/check_submission_consistency.py` verifies the central submitted values against repository CSV/JSON outputs.
- `results/Table_S2_internal_cv_manuscript_crosswalk.csv` maps internal cross-validation values between the manuscript summary and exact reproducible output.
- `results/Table_1_denominator_crosswalk.csv` documents patient-characteristics versus metric-evaluable denominators.

## Repository Structure

```
src/              survival.py, meta.py, pathways.py (authoritative pathway defs),
                  baselines.py; src/ml/ (deepsurv.py, wrappers.py = pickleable
                  model wrappers). NOTE: features.py, models.py, data_loader.py,
                  and preprocessing.py are DEPRECATED v1 classification modules
                  (notebook-only); see their headers.
src/ml/           DeepSurv + the survival-model wrappers used by the saved weights
R/                Official genefu PAM50-ROR scoring (pam50.R); install.R for R deps
scripts/          Analysis pipeline + reviewer helpers (see inventory below)
notebooks/        Jupyter notebooks (01-14: original analyses)
results/          All results tables (Table_1..6, S1..S17) and QC reports
figures/          Generated figures
data/             Raw, processed, and clinical data (committed intermediates +
                  download/regeneration scripts; large raw expression not redistributed)
models/           Trained model artifacts (committed; see models/README.md + SHA256)
tests/            Pytest suite (incl. recompute smoke + primary-endpoint tests)
docs/             Pre-registration endpoint, cohort notes, data dictionary/use,
                  harmonization, hyperparameters, model SHA256 manifest
Makefile, Dockerfile, .github/workflows/ci.yml, pyproject.toml   reproducibility infra
```

### Scripts inventory

Pipeline (in order): `harmonize_cohorts.py` -> `score_pathways.py` ->
`compute_clinical_baselines.py` (R/genefu) -> `train_ml_zoo.py` ->
`external_validation.py`. Phase-3 exploratory: `within_subtype_analysis.py`,
`calibration_analysis.py`, `dca_analysis.py`, `stability_analysis.py`
(`phase3_common.py` is a shared helper). Figures/tables: `build_table1.py`,
`build_cohort_flow.py`, `build_forest_plot.py`, `compute_model_attributions.py`
(feature attributions backing R1.6/R3.4). Sensitivity: `feature_ablation.py`,
`pathway_scoring_sensitivity.py`, `external_only_tnbc_sensitivity.py`,
`recurrence_sensitivity.py` (R2.1/R3.1, Table S16), `learning_curve.py`
(R1.2/R2.2/R3.2, Table S17). Reviewer helpers: `predict.py` (load a committed
weight and score), `migrate_pickles.py` (make pickles load cold),
`check_submission_consistency.py` (recompute-vs-committed gate),
`build_manuscript_docx.py`. `run_rescue_analysis.py` / `screen_rescue_top_candidates.py`
are exploratory rescue analyses (over-fit; not a result claim).

## Manuscript

The manuscript and supplementary materials are maintained separately from this
code repository and are not included here. This repository is the code, data,
and reproducibility artifacts that back the submitted results.

## Data Sources

- **TCGA-BRCA:** GDC Data Portal (https://portal.gdc.cancer.gov)
- **GSE96058/SCAN-B:** NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058)
- **METABRIC:** cBioPortal (https://www.cbioportal.org/study/summary?id=brca_metabric)
- **GSE20685:** NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685)

Large raw expression matrices are not redistributed (size + source terms; see
`docs/DATA_USE.md`). The committed clinical CSVs and derived pathway/baseline
intermediates are de-identified open-access data. Source URLs and SHA256s for
every raw input are in `data/processed/*.metadata.json`; download with
`data/raw/download_data.sh` and `data/raw/download_gene_sets.sh`.

## Pre-registration

The primary endpoint was locked before any external validation was run. See `docs/PRIMARY_ENDPOINT.md` (committed at `7783f4f0bfaaa6bdc611c78d33ccda621c6b243d`, before any results commits).

## Reproduction

The trained model weights and the pipeline intermediates they consume are
committed, so the headline results reproduce **without** downloading the
multi-GB raw expression data (two tiers documented in
`data/processed/REGENERATION.md`).

```bash
pip install -r requirements.txt          # pinned; or use the Dockerfile (adds R/genefu)

make smoke                # recompute the headline C-index + TNBC endpoint from weights+data
make reproduce-headline   # regenerate the external-validation tables + endpoint, then verify
make test                 # full pytest suite (PAM50-ROR test needs R/genefu; auto-skips otherwise)
make consistency          # recompute-vs-committed gate
make all                  # Tier-1: train -> external_validation -> phase3 -> figures -> consistency
```

Load a committed weight directly:

```bash
python3 scripts/predict.py --from-features data/processed/04_pathway_features.parquet --cohort TCGA-BRCA
```

Full pipeline order (Tier-2, needs raw data + R/genefu): `harmonize_cohorts.py`
-> `score_pathways.py` -> `compute_clinical_baselines.py` -> `train_ml_zoo.py`
-> `external_validation.py` -> phase-3/figures (`make all-from-raw`). Install R
deps with `Rscript R/install.R`.

## License

MIT License (code). See `LICENSE`. Third-party data remain under their source
terms — see `docs/DATA_USE.md`.
