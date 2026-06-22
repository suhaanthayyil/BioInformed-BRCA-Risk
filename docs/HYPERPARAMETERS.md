# Hyperparameter Search Spaces and Model Selection

Addresses reviewer comment **R2.3**.

This document surfaces the hyperparameter search spaces, selection rule, and
selected configurations that are already coded and recorded in the repository.
It introduces no new analysis — every search space below is transcribed from
`scripts/train_ml_zoo.py` (`model_spaces()`), and every selected configuration
is taken verbatim from `results/Table_S2_ml_internal_cv.csv` (`config` column)
and `data/processed/ml_model_zoo.metadata.json`.

## Selection rule (all models)

- **Training cohort:** TCGA-BRCA only (`n_training_samples = 213`,
  `n_events = 132`), after excluding rows with non-positive or missing
  `os_days`/`os_event`.
- **Cross-validation:** 5-fold **stratified** CV (`StratifiedKFold`,
  `n_splits = 5`, `shuffle = True`), stratified on the event indicator.
- **Selection metric:** **Harrell's concordance index (C-index)**,
  `concordance_index_censored` from scikit-survival, computed out-of-fold; the
  configuration with the highest mean 5-fold C-index is selected per model
  (`best = max(candidates, key=mean_cindex)`).
- **Seed:** `SEED = 42` (NumPy, PyTorch, and all estimator `random_state`).
- **Preprocessing inside each fold:** median imputation then standard scaling
  (`SimpleImputer(strategy="median")` -> `StandardScaler`), fit on the training
  fold only to avoid leakage.
- **Feature set (fixed, 9 columns):** `Pathway_Immune`, `Pathway_Proliferation`,
  `Pathway_DNA_Repair`, `Pathway_Metabolism`, `Pathway_Stromal_EMT`,
  `Pathway_Apoptosis_Stress`, `Pathway_Hormone`, `age_at_dx`, `stage_ordinal`.

### Grid-mode caveat (recorded in metadata)

`data/processed/ml_model_zoo.metadata.json` records `"grid_mode":
"compact_local"` with the note:

> "Compact grid used for local tractability; XGBoost Cox disabled by default
> after repeated native process termination. Set ENABLE_XGBOOST_SURVIVAL=1 to
> retry that optional cross-check."

The QC report (`results/reports/ml_model_zoo_qc.md`, written by the same script)
states: *"Grid mode: compact local grid. This is explicitly recorded and should
not be described as the full protocol grid."* These grids are the
compact local grids actually executed; they should be reported as such.

## Per-model search spaces and selected configurations

C-index values below are the exact repository values from
`results/Table_S2_ml_internal_cv.csv` (mean of 5 out-of-fold folds).

### Cox_PH (penalized Cox, lifelines `CoxPHFitter`)
- **Search space (3 configs):** `penalizer in {0.01, 0.1, 1.0}`.
- **Selected:** `{"penalizer": 0.01}`.
- **Selected mean CV C-index:** 0.599862 (reference baseline; `delta_vs_cox = 0`).

### Elastic_Net_Cox (scikit-survival `CoxnetSurvivalAnalysis`)
- **Search space (12 configs):** `alpha in {0.001, 0.01, 0.1, 1.0}` x
  `l1_ratio in {0.1, 0.5, 0.9}` (one alpha per fit, `fit_baseline_model=True`,
  `max_iter=100000`).
- **Selected:** `{"alpha": 0.001, "l1_ratio": 0.1}`.
- **Selected mean CV C-index:** 0.601334 (delta vs Cox +0.001472).

### Random_Survival_Forest (scikit-survival `RandomSurvivalForest`)
- **Search space (18 configs):** `n_estimators in {500}` x
  `max_depth in {None, 4, 8}` x `min_samples_leaf in {1, 5, 10}` x
  `max_features in {"sqrt", None}` (`random_state=42`, `n_jobs=2`).
- **Selected:** `{"n_estimators": 500, "max_depth": 4, "min_samples_leaf": 1, "max_features": "sqrt"}`.
- **Selected mean CV C-index:** 0.600705 (delta vs Cox +0.000843).

### Gradient_Boosted_Survival (scikit-survival `GradientBoostingSurvivalAnalysis`) — HEADLINE
- **Search space (16 configs):** `n_estimators in {100, 300}` x
  `learning_rate in {0.03, 0.1}` x `max_depth in {1, 3}` x
  `subsample in {0.8, 1.0}` (`random_state=42`).
- **Selected:** `{"n_estimators": 100, "learning_rate": 0.03, "max_depth": 1, "subsample": 1.0}`.
- **Selected mean CV C-index:** 0.641826 (delta vs Cox +0.041964). This is the
  headline model (`headline_model = "Gradient_Boosted_Survival"` in the
  metadata).

### DeepSurv (PyTorch MLP, `src/ml/deepsurv.py`)
- **Search space (4 configs):**
  - `{hidden_dims:(32,),    dropout:0.1, lr:1e-3, weight_decay:1e-4, epochs:120, patience:15}`
  - `{hidden_dims:(64,32),  dropout:0.1, lr:1e-3, weight_decay:1e-4, epochs:150, patience:20}`
  - `{hidden_dims:(64,32),  dropout:0.2, lr:5e-4, weight_decay:1e-4, epochs:150, patience:20}`
  - `{hidden_dims:(128,64), dropout:0.2, lr:5e-4, weight_decay:1e-4, epochs:150, patience:20}`
- **Selected:** `{"hidden_dims": [128, 64], "dropout": 0.2, "lr": 0.0005, "weight_decay": 0.0001, "epochs": 150, "patience": 20}`.
- **Selected mean CV C-index:** 0.609042 (delta vs Cox +0.009181).

### Stacked_Ensemble (Cox meta-learner over top-3 base models)
- **Construction:** not a grid search. The top-3 models by mean CV C-index are
  selected, their out-of-fold risk predictions are stacked, and a penalized Cox
  meta-model (`CoxPHFitter(penalizer=0.1)`) is fit on those OOF predictions.
- **Selected (top-3 inputs):** `Gradient_Boosted_Survival`, `DeepSurv`,
  `Elastic_Net_Cox`.
- **Mean CV C-index (on OOF stack):** 0.626274 (delta vs Cox +0.026412).

### XGBoost_Cox (optional, disabled by default)
- **Status:** Not run in the default pipeline. Gated behind environment variable
  `ENABLE_XGBOOST_SURVIVAL=1`, because the local native XGBoost library
  repeatedly terminated the process without a Python traceback (recorded in the
  metadata note and QC report).
- **Search space when enabled (4 configs):** `n_estimators in {100}` x
  `learning_rate in {0.03, 0.1}` x `max_depth in {1, 3}` x `subsample in {0.8}`,
  with `colsample_bytree = 0.9` (objective `survival:cox`).

## Manuscript rounding crosswalk

`results/Table_S2_internal_cv_manuscript_crosswalk.csv` maps the rounded summary
values reported in the submitted manuscript to the exact reproducible repository
values. The headline GBSA value (manuscript 0.642 vs repository 0.641826) and
the Cox value (0.600 vs 0.599862) match after rounding; the repository exact
values are authoritative for reproduction.
