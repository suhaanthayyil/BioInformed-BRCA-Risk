# Results Log

## Submitted Manuscript Crosswalk

The submitted BMC Artificial Intelligence manuscript is mapped to repository source files in `docs/SUBMITTED_MANUSCRIPT_REPOSITORY_CROSSWALK.md`.

Additional traceability files:

- `scripts/check_submission_consistency.py`: executable consistency check for central manuscript values.
- `results/Table_S2_internal_cv_manuscript_crosswalk.csv`: internal CV manuscript summary versus exact repository output.
- `results/Table_1_denominator_crosswalk.csv`: cohort denominator conventions.

Important denominator note: METABRIC has 1,980 patients in the patient-characteristics set but 1,979 in C-index calculations because one sample has non-positive survival time.

## Sensitivity Analyses

### Table S15: External-only TNBC Sensitivity

Restricts the TNBC delta C-index analysis to external cohorts only (GSE96058, METABRIC), excluding TCGA-BRCA (training set).

| Cohort | n | Events | GBSA C | PAM50 C | Delta | 95% CI | p |
|--------|---|--------|--------|---------|-------|--------|---|
| GSE96058 | 55 | 27 | 0.6055 | 0.5759 | +0.0295 | [-0.1313, 0.1947] | 0.7291 |
| METABRIC | 319 | 168 | 0.5150 | 0.5250 | -0.0100 | [-0.0696, 0.0508] | 0.7435 |
| **Meta RE** | 374 | 195 | -- | -- | **-0.0052** | [-0.0617, 0.0512] | 0.8557 |
| **Meta FE** | 374 | 195 | -- | -- | **-0.0052** | [-0.0617, 0.0512] | 0.8557 |

tau2 = 0.0000, I2 = 0.0%.

### Table 6 / Table S5: Feature Ablation

Eight feature-set combinations trained with locked GBSA architecture on TCGA, evaluated on three external cohorts.

| Feature Set | n_features | Mean Harrell C | Mean Delta vs PAM50 |
|-------------|-----------|----------------|---------------------|
| Pathways+Clinical | 9 | 0.5942 | -0.0267 |
| Clinical_only | 2 | 0.5740 | -0.0470 |
| Leave_one_out_Hormone | 8 | 0.5481 | -0.0728 |
| Hormone_only | 1 | 0.5438 | -0.0771 |
| Pathways_only | 7 | 0.5156 | -0.1054 |
| Top3_pathways | 3 | 0.5080 | -0.1130 |
| Immune_only | 1 | 0.4935 | -0.1275 |
| Proliferation_only | 1 | 0.4698 | -0.1511 |

### Table S14: Pathway Scoring Sensitivity

Three aggregation methods with locked GBSA architecture:

| Method | Mean Harrell C | Mean Delta vs PAM50 |
|--------|----------------|---------------------|
| rank_percentile (production) | 0.5942 | -0.0267 |
| ssgsea (gseapy) | 0.5931 | -0.0278 |
| mean_z | 0.5673 | -0.0537 |

Rank-percentile and ssGSEA perform comparably; mean-z shows modest degradation.
