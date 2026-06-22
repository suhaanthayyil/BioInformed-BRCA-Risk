# Submitted Manuscript to Repository Crosswalk

Date: 2026-05-24

Submitted manuscript reviewed: `bmcmanu.docx` (submitted manuscript file)

Submitted manuscript SHA256: `b01ef9151a5780b118edcda0841d848432ae2aec0892472261a9d4ec4cd461ac`

Repository baseline commit before this reconciliation: `91878fdb72ec6d400232357dcacfbf5a9da44dfe`

## Purpose

This document maps the submitted BMC Artificial Intelligence manuscript to the exact repository files that support each central result. It exists so reviewers can trace the submitted claims without ambiguity.

The analytic endpoint in the original pre-registration file is intentionally left unchanged. The venue line reflects the current submitted journal.

## Authorship and Contributions

The repository-facing authors are:

- Suhaan Thayyil
- Eshaan Nidee

No software tool, automation, or third-party system is listed as an author, contributor, or acknowledgement.

## Manuscript Framing

The manuscript title is:

> Pathway-based machine learning is competitive with but does not exceed PAM50-ROR for breast cancer prognosis: a pre-registered multi-cohort benchmark in 4,532 patients

The repository supports this framing. The pre-registered TNBC endpoint was not met, and the model should not be described as superior to PAM50-ROR.

## Primary Endpoint Trace

Primary endpoint source files:

- `docs/PRIMARY_ENDPOINT.md`
- `data/processed/phase6_primary_endpoint.json`
- `results/Table_3_head_to_head.csv`

Locked result:

- Model: `Gradient_Boosted_Survival`
- Comparator: official `PAM50_ROR_official`
- Cohorts in TNBC meta-analysis: TCGA-BRCA, GSE96058, METABRIC
- TNBC delta C-index: `+0.0144`
- 95% CI: `[-0.0472, +0.0760]`
- p-value: `0.6466`
- I2: `10.5%`
- Endpoint status: not met

## Pre-registration and Model-Selection Wording

The locked pre-registration file specifies:

- TNBC random-effects delta C-index versus PAM50-ROR as the primary endpoint.
- Success threshold of delta C-index `>= +0.03` and `p < 0.05`.
- Selection of the best non-linear machine-learning model by internal cross-validation.

The submitted manuscript's wording about a locked model-selection rule should be read as the implementation layer of the internal cross-validation selection process. The original pre-registration document remains the controlling historical record.

The 50/50 overall/TNBC composite score appears in the repository as part of the transparent post-primary rescue/sensitivity analysis, not as a change to the locked primary endpoint:

- `results/reports/FINAL_SUMMARY.md`
- `data/processed/rescue_analysis.metadata.json`
- `results/Table_S4_rescue_internal_cv.csv`

## Internal Cross-Validation Trace

Primary source file:

- `results/Table_S2_ml_internal_cv.csv`

Crosswalk file:

- `results/Table_S2_internal_cv_manuscript_crosswalk.csv`

The headline model value in the manuscript, GBSA C-index `0.642`, matches the repository exact value `0.641826` after rounding. The repository exact values should be treated as authoritative for reproduction.

## External Discrimination Trace

Primary source files:

- `results/Table_2_external_validation.csv`
- `results/Table_3_head_to_head.csv`

Submitted manuscript external C-index values trace as follows:

- GSE20685: GBSA `0.565`, PAM50-ROR `0.635`, delta `-0.070`
- GSE96058: GBSA `0.663`, PAM50-ROR `0.633`, delta `+0.030`
- METABRIC: GBSA `0.555`, PAM50-ROR `0.595`, delta `-0.040`
- TCGA-BRCA training reference: GBSA `0.673`, PAM50-ROR `0.504`, delta `+0.169`

These values match the repository after rounding.

## Denominator Conventions

Primary source files:

- `results/Table_1_patient_characteristics.csv`
- `results/Table_1_denominator_crosswalk.csv`
- `results/Table_2_external_validation.csv`

The manuscript uses two related denominators:

- Harmonized/evaluable cohort counts for patient characteristics.
- Metric-evaluable counts for discrimination models, which exclude non-positive survival time and model-specific missingness.

METABRIC has one patient with non-positive survival time in the discrimination dataset. This is why patient-characteristics tables can show `1,980` evaluable samples while model C-index tables can show `1,979`.

Receptor-status percentages are calculated among samples with known receptor status.

## Calibration, DCA, and Reclassification Trace

Calibration:

- Source file: `results/Table_S10_calibration.csv`
- Submitted 5-year mean ICI: GBSA `0.125`, PAM50-ROR `0.155`

Decision curve analysis:

- Source file: `results/Table_S11_dca.csv`
- Submitted average net-benefit delta: `+0.016`

NRI/IDI:

- Source file: `results/Table_S12_nri_idi.csv`
- IDI is positive in every cohort.
- The thresholds are 5-year mortality-risk thresholds, not survival-probability thresholds.

## Stability Trace

Primary source file:

- `results/Table_S13_stability.csv`

Implementation source file:

- `scripts/stability_analysis.py`

The analysis used 100 stratified 80% TCGA subsamples. The submitted manuscript's stability result matches the repository:

- Mean pairwise Spearman rho: `0.395`
- 95% empirical interval: `[-0.356, +0.917]`
- Top-3 frequency: Pathway_Hormone `92%`, age_at_dx `75%`, Pathway_DNA_Repair `62%`

## Sensitivity Analyses Trace

External-only TNBC sensitivity:

- `results/Table_S15_external_only_tnbc.csv`
- Meta delta `-0.0052`
- p-value `0.8557`

Feature ablation:

- `results/Table_6_feature_ablation_summary.csv`
- `results/Table_S5_feature_ablation.csv`

Pathway-scoring sensitivity:

- `results/Table_S14_pathway_scoring_sensitivity.csv`

These analyses support the manuscript's conclusion that no tested sensitivity analysis establishes external superiority over PAM50-ROR.

## Supplementary Table Index

Within-subtype analyses (cited in the main text):

- `results/Table_S8_within_subtype_external.csv` -- per-cohort within-subtype external discrimination (overall vs PAM50-ROR), now labelled `analysis_type=exploratory_secondary`.
- `results/Table_S9_within_subtype_meta.csv` -- random-effects meta of the within-subtype deltas.

Numbering note: there is no `Table_S3`. Supplementary tables S4-S7 are the
post-hoc *rescue* analyses (disavowed as over-fit in `docs/STORY.md`); the gap
at S3 is intentional and is not a missing file.

Revision-1 additions (BMC Artificial Intelligence major revision):

- `results/Table_S16_metabric_recurrence_sensitivity.csv` -- locked GBSA vs official PAM50-ROR on the METABRIC DFS/recurrence endpoint (reviewers R2.1/R3.1).
- `results/Table_S17_learning_curve.csv` and `results/Table_S17_train_cohort_crosscheck.csv` -- training-size sensitivity and train-on-each-cohort cross-check (reviewers R1.2/R2.2/R3.2).

## Clinical Covariate Availability

The model feature schema contains age and ordinal stage. Age is available for nearly all included samples. Stage is stored as a harmonized field where source cohorts provide it and is handled by the preprocessing imputer where missing.

This is the intended reading of manuscript language describing shared clinical covariates.

## Repository QC

Current checks (all green at the revision-1 tip):

- `python3 -m pytest tests/` -- includes the recompute smoke test
  (`tests/test_smoke.py`, reproduces the headline per-cohort C-index from the
  committed weights + data) and the TNBC primary-endpoint recompute
  (`tests/test_primary_endpoint.py`).
- `python3 -m ruff check .` (config in `pyproject.toml`).
- `python3 scripts/check_submission_consistency.py` -- recompute-vs-committed
  gate; exits non-zero on any mismatch.
- `make smoke` / `make reproduce-headline` reproduce the headline results from a
  fresh clone without raw data; CI (`.github/workflows/ci.yml`) runs all of the
  above on every push.

## Reviewer-Facing Summary

The repository is aligned with the submitted manuscript's central claim:

> Pathway-based survival machine learning was competitive with but did not significantly exceed official PAM50-ROR, including in the pre-registered TNBC primary endpoint.

The locked primary endpoint is negative and is preserved transparently.
