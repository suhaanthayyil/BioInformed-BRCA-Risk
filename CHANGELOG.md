# Changelog

All notable changes to this repository. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [1.1.0] — Revision 1 (BMC Artificial Intelligence major revision)

Repository changes made for the major-revision resubmission. The scientific
result is unchanged: the pre-registered TNBC superiority endpoint remains **not
met** (meta ΔC-index +0.0144, 95% CI [-0.0472, +0.0760], p = 0.6466).

### Added (reproducibility — editor E3)
- Trained model weights committed under `models/` (6 artifacts) with
  `models/README.md` and `docs/MODEL_SHA256.txt`.
- `src/ml/wrappers.py` (pickleable model wrappers) and `scripts/migrate_pickles.py`
  so the weights load cold (no `__main__` shim); `scripts/predict.py` reviewer
  loader/scorer.
- Pipeline intermediates committed under `data/processed/`
  (`04_pathway_features.parquet`, `pathway_scores_all.parquet`,
  `baselines_pam50.parquet`, `unified_cohorts.duckdb`).
- `tests/test_smoke.py` (recompute headline C-index from weights+data) and
  `tests/test_primary_endpoint.py` (recompute TNBC endpoint).
- `Makefile`, `Dockerfile`, `.github/workflows/ci.yml`, `pyproject.toml`,
  `conftest.py`, `.python-version`, `R/install.R`.
- `data/raw/download_gene_sets.sh`; `download_data.sh` extended to METABRIC +
  GSE20685 with SHA256 verification; `data/processed/REGENERATION.md`.

### Added (reviewer-driven analyses)
- `scripts/recurrence_sensitivity.py` → `results/Table_S16_metabric_recurrence_sensitivity.csv`
  (R2.1/R3.1: METABRIC DFS/recurrence endpoint).
- `scripts/learning_curve.py` → `results/Table_S17_learning_curve.csv`,
  `results/Table_S17_train_cohort_crosscheck.csv`, `figures/fig_learning_curve.*`
  (R1.2/R2.2/R3.2: training-size sensitivity + train-on-each-cohort cross-check).

### Added (documentation)
- `docs/RESPONSE_TO_REVIEWERS.md`, `docs/DATA_DICTIONARY.md`, `docs/DATA_USE.md`,
  `docs/HYPERPARAMETERS.md`, `docs/HARMONIZATION.md` (+
  `results/reports/cohort_gene_overlap.csv`).

### Changed
- Regenerated `data/processed/02_tcga_feature_matrix.csv` with the seven
  authoritative pathways (was a stale gene-level schema); provenance clarified in
  `docs/COHORT_NOTES.md`.
- Secondary-analysis tables labelled `analysis_type=exploratory_secondary`
  (Table_S8/S9/S10/S11); `Table_S13` labelled
  `interpretation=descriptive_signal_not_mechanistic` (R1.1/R3.3/R1.6/R3.4).
- Forest plots annotated with per-cohort event counts; pooled row relabelled
  (R1.3/R1-minor).
- Deprecation headers on the orphaned v1 `src/` modules
  (`features/models/data_loader/preprocessing.py`); `src/pathways.py` noted as
  the authoritative pathway definition.
- `CITATION.cff` completed (version, date-released, repository-code).
- Fixed stale `render_docx.py` references in the audit docs.

### Notes
- Manuscript prose, the two R2.4 bibliography DOIs, and the R2.5 typo are
  authored separately (see `docs/RESPONSE_TO_REVIEWERS.md`).
- No archival DOI (Zenodo) minted for this revision.
