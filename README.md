# BRCA-PathwayML

[![CI](https://github.com/suhaanthayyil/BioInformed-BRCA-Risk/actions/workflows/ci.yml/badge.svg)](https://github.com/suhaanthayyil/BioInformed-BRCA-Risk/actions/workflows/ci.yml)

Pathway-based survival machine learning for breast cancer prognosis: a four-cohort benchmark against genefu PAM50 ROR-S.

**Authors:** Suhaan Thayyil, Eshaan Nidee

## Overview

Survival models trained on seven MSigDB pathway scores in TCGA-BRCA and applied without refitting to GSE96058/SCAN-B, METABRIC and GSE20685, compared with genefu::rorS, a research implementation of PAM50 ROR-S (not the clinical Prosigna assay).

The primary endpoint was registered on 2026-05-16 (`docs/PRIMARY_ENDPOINT.md`). TCGA-BRCA and GSE96058 had already been analysed in the April 2026 preprint; METABRIC and GSE20685 were new at registration. Every deviation from the registered plan is listed in `docs/REGISTRATION_DEVIATIONS.md` (Table S20).

## Revision (October 2026)

The analysis was re-run after the journal's pre-publication checks:

- **Cohorts rebuilt from source.** The submitted analysis reused v1 cohorts built for a binary 5-year outcome: patients censored before 5 years were excluded, TCGA (n = 213) contained 37 normal-tissue and 3 metastatic samples, and GSE96058 (n = 1,483) contained 60 technical replicates. All eligible patients are now included with no follow-up restriction (TCGA-BRCA 1,064; GSE96058 3,273; METABRIC 1,980; GSE20685 327). See `docs/COHORT_NOTES.md` and Figure 1.
- **Stage coding fixed.** The old parser only read Roman numerals, so stage was missing for every external patient. `src/clinical.py` now reads TCGA, METABRIC and GSE20685 (TNM) stage; GSE96058 releases none.
- **Model inputs as registered.** Seven pathway scores only. Age and stage are sensitivity analyses (Table 6, Table S5).
- **Headline rule as registered.** Best of RSF, GBSA and DeepSurv by internal CV; Cox and elastic-net Cox are comparators. DeepSurv early stopping no longer sees the held-out fold.
- **TCGA in pooled analyses** uses out-of-fold predictions, not in-sample ones.
- **All registered endpoints reported** (Table S18), and the 2026-05-16 locked weights are evaluated on the full cohorts (Table S19, `models/locked_v1/`).

## Key Results

- **Headline model:** DeepSurv, chosen from RSF, GBSA and DeepSurv on TCGA only (expanded grid, top 5 per model re-scored with 5x5 repeated CV; C-index 0.593; Cox PH 0.614).
- **TNBC primary endpoint:** NOT MET. Delta C vs PAM50 ROR-S = +0.020 (95% CI -0.040 to +0.081), p = 0.51; external cohorts only +0.013 (p = 0.68). Comparable to PAM50 ROR-S in TNBC, not superior.
- **Overall discrimination:** lower than PAM50 ROR-S in GSE96058 (-0.077) and METABRIC (-0.043); registered non-inferiority (margin -0.02) not met.
- **Age confounding:** the endpoint is all-cause overall survival. Age alone exceeds PAM50 ROR-S in GSE96058, and adding age to the pathways reproduces the apparent advantage of the submitted model (Table S5). With age given to both sides, PAM50 ROR-S still adds more than the pathway model (Table S21). The post-hoc rescue search selected a pathways + age + stage Cox model (Tables S4-S7); its advantage comes from the clinical covariates.

## Reproducibility traceability

- `scripts/check_submission_consistency.py` checks the central values against the CSV/JSON outputs.
- `results/Table_1_denominator_crosswalk.csv` documents patient-characteristics versus metric-evaluable denominators.
- `data/processed/ingest_flow.json` records every inclusion step and the SHA256 of every raw input.

## Repository Structure

```
src/              survival.py, meta.py, pathways.py (authoritative pathway defs),
                  clinical.py (feature lists, stage parsing), baselines.py;
                  src/ml/ (deepsurv.py, wrappers.py, registry.py = model loading
                  and retraining). NOTE: features.py, models.py, data_loader.py,
                  and preprocessing.py are DEPRECATED v1 classification modules
                  (notebook-only); see their headers.
src/ml/           DeepSurv + the survival-model wrappers used by the saved weights
R/                genefu PAM50 ROR-S (research implementation, not the Prosigna assay) scoring (pam50.R); install.R for R deps
scripts/          Analysis pipeline + reviewer helpers (see inventory below)
notebooks/        Jupyter notebooks (01-14: original analyses)
results/          Results tables (Table_1..7, S1..S21), QC reports; legacy_v1_preprint/ and
                  archive_2026-05_submission/ hold superseded outputs
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

Pipeline (in order): `ingest_raw.py` -> `harmonize_cohorts.py` -> `score_pathways.py` ->
`compute_clinical_baselines.py` (R/genefu) -> `train_ml_zoo.py` ->
`external_validation.py`. Registered endpoints: `registered_endpoints.py`
(Table S18); locked-model sensitivity: `locked_model_sensitivity.py` (Table S19).
Phase-3 exploratory: `within_subtype_analysis.py`,
`calibration_analysis.py`, `dca_analysis.py`, `stability_analysis.py`
(`phase3_common.py` is a shared helper). Figures/tables: `build_table1.py`,
`build_figure1_consort.py`, `build_forest_plot.py`, `compute_model_attributions.py`
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

## Registration

The primary endpoint was committed at `7783f4f` on 2026-05-16, before METABRIC and GSE20685 were obtained but after the TCGA-BRCA to GSE96058 analysis of the April 2026 preprint (first commit, 2026-04-01). The corrected timeline is in `docs/PRIMARY_ENDPOINT.md`; deviations are in `docs/REGISTRATION_DEVIATIONS.md`.

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

Full pipeline order (Tier-2, needs raw data + R/genefu): `ingest_raw.py` -> `harmonize_cohorts.py`
-> `score_pathways.py` -> `compute_clinical_baselines.py` -> `train_ml_zoo.py`
-> `external_validation.py` -> phase-3/figures (`make all-from-raw`). Install R
deps with `Rscript R/install.R`.

## License

MIT License (code). See `LICENSE`. Third-party data remain under their source
terms — see `docs/DATA_USE.md`.
