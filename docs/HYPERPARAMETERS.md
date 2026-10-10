# Hyperparameter Search Spaces and Model Selection

Addresses reviewer comment **R2.3**.

This document surfaces the hyperparameter search spaces, selection rule, and
selected configurations that are already coded and recorded in the repository.
It introduces no new analysis — every search space below is transcribed from
`scripts/train_ml_zoo.py` (`model_spaces()`), and every selected configuration
is taken verbatim from `results/Table_S2_ml_internal_cv.csv` (`config` column)
and `data/processed/ml_model_zoo.metadata.json`.

## Selection rule (all models)

- **Training cohort:** TCGA-BRCA only (`n_training_samples = 1,064`,
  `n_events = 150`): all eligible female primary-tumour patients with TCGA-CDR
  overall survival, no follow-up restriction.
- **Cross-validation:** 5-fold **stratified** CV (`StratifiedKFold`,
  `n_splits = 5`, `shuffle = True`), stratified on the event indicator.
- **Selection metric:** **Harrell's concordance index (C-index)**,
  `concordance_index_censored` from scikit-survival, computed out-of-fold;
  every configuration is scored with 5-fold CV, the top 5 per model are
  re-scored with 5x5 repeated stratified CV (25 folds), and the highest
  25-fold mean is selected (`refine()` in `scripts/train_ml_zoo.py`). The same
  rule applies to every model, including the Cox comparators.
- **Seed:** `SEED = 42` (NumPy, PyTorch, and all estimator `random_state`).
- **Preprocessing inside each fold:** median imputation then standard scaling
  (`SimpleImputer(strategy="median")` -> `StandardScaler`), fit on the training
  fold only to avoid leakage.
- **Feature set (fixed, 7 columns, as registered):** `Pathway_Immune`,
  `Pathway_Proliferation`, `Pathway_DNA_Repair`, `Pathway_Metabolism`,
  `Pathway_Stromal_EMT`, `Pathway_Apoptosis_Stress`, `Pathway_Hormone`.
- **Headline model:** the registered candidates (RSF, GBSA, DeepSurv) only; the
  one with the highest mean internal CV C-index is the headline. Cox PH,
  elastic-net Cox and the stacked ensemble are comparators.
- **DeepSurv early stopping:** a stratified 20% split of the training fold
  (never the held-out fold); CPU, single thread, seed 42.

### Tuning history (all rounds disclosed)

| Round | When | Grid (RSF / GBSA / DeepSurv) | Selection | Headline | Internal CV C |
|---|---|---|---|---|---|
| Submitted | 2026-05-16 | 18 / 16 / 4, 9 features (pathways + age + broken stage), TCGA n = 213 | best 5-fold mean over all six models | GBSA | 0.642 |
| Revision 0 (not used) | 2026-10-01 | 18 / 16 / 4, seven pathways + age, TCGA n = 1,064 | best registered non-linear model, 5-fold | DeepSurv | 0.656 |
| Revision 1 | 2026-10-01 | 18 / 16 / 4, seven pathways, TCGA n = 1,064 | best registered non-linear model, 5-fold | DeepSurv | 0.604 |
| Revision 2 | 2026-10-02 | 75 / 81 / 40 | same rule, 5-fold | DeepSurv | 0.610 |
| Revision 3 (final) | 2026-10-03 | 75 / 81 / 40 | top 5 per model re-scored with 5x5 repeated CV | DeepSurv | 0.593 |

Revision 0 added age to the pathways. Its TNBC primary delta was +0.074
(p = 0.019), but age confounds a comparison with PAM50 ROR-S on all-cause
survival, so the analysis returned to the registered seven-pathway inputs
(revision 1 onward); the pathways + age model is reported as a sensitivity
analysis (Table 6, Table S5). Every round selected models with TCGA-BRCA
internal cross-validation only. The
external cohorts were not blind during the revision: their results under the
submitted model, and under each earlier round, had been seen before the next
round was run. Revision 3 is final. External results for rounds 1-3 were
similar (TNBC primary delta +0.008, +0.020, +0.020; overall delta in GSE96058
-0.083, -0.077, -0.077).

## Per-model search spaces and selected configurations

C-index values below are the exact repository values from
`results/Table_S2_ml_internal_cv.csv` (mean over the 25 repeated-CV folds).

### Cox_PH (penalized Cox, lifelines `CoxPHFitter`)
- **Search space (3 configs):** `penalizer in {0.01, 0.1, 1.0}`.
- **Selected:** `{"penalizer": 0.1}`.
- **Selected mean CV C-index:** 0.614014 (delta vs Cox +0.000000).

### Elastic_Net_Cox (scikit-survival `CoxnetSurvivalAnalysis`)
- **Search space (12 configs):** `alpha in {0.001, 0.01, 0.1, 1.0}` x
  `l1_ratio in {0.1, 0.5, 0.9}` (one alpha per fit, `fit_baseline_model=True`,
  `max_iter=100000`).
- **Selected:** `{"alpha": 0.01, "l1_ratio": 0.1}`.
- **Selected mean CV C-index:** 0.613778 (delta vs Cox -0.000236).

### Random_Survival_Forest (scikit-survival `RandomSurvivalForest`)
- **Search space (75 configs):** `n_estimators = 500` x
  `max_depth in {None, 3, 5, 8, 12}` x `min_samples_leaf in {3, 5, 10, 20, 40}` x
  `max_features in {"sqrt", 0.5, None}` (`random_state=42`).
- **Selected:** `{"n_estimators": 500, "max_depth": 5, "min_samples_leaf": 40, "max_features": 0.5}`.
- **Selected mean CV C-index:** 0.577450 (delta vs Cox -0.036564).

### Gradient_Boosted_Survival (scikit-survival `GradientBoostingSurvivalAnalysis`)
- **Search space (81 configs):** `n_estimators in {100, 300, 600}` x
  `learning_rate in {0.01, 0.03, 0.1}` x `max_depth in {1, 2, 3}` x
  `subsample in {0.5, 0.8, 1.0}` (`random_state=42`).
- **Selected:** `{"n_estimators": 100, "learning_rate": 0.01, "max_depth": 1, "subsample": 0.8}`.
- **Selected mean CV C-index:** 0.573880 (delta vs Cox -0.040135).

### DeepSurv (PyTorch MLP, `src/ml/deepsurv.py`) — HEADLINE
- **Search space (40 configs):** `hidden_dims in {(16,), (32,), (64,), (64, 32), (128, 64)}` x
  `dropout in {0.1, 0.3}` x `lr in {1e-3, 3e-4}` x `weight_decay in {1e-4, 1e-2}`;
  `epochs = 300`, `patience = 30`, early stopping on a 20% split of the training fold.
- **Selected:** `{"hidden_dims": [128, 64], "dropout": 0.3, "lr": 0.0003, "weight_decay": 0.0001, "epochs": 300, "patience": 30}`.
- **Selected mean CV C-index:** 0.592974 (delta vs Cox -0.021041). This is the headline model.

### Stacked_Ensemble (Cox meta-learner over top-3 base models)
- **Construction:** not a grid search. The top-3 models by mean CV C-index are
  selected, their out-of-fold risk predictions are stacked, and a penalized Cox
  meta-model (`CoxPHFitter(penalizer=0.1)`) is fit on those OOF predictions.
- **Selected (top-3 inputs):** `Cox_PH`, `Elastic_Net_Cox`, `DeepSurv`.
- **Mean CV C-index (on OOF stack):** 0.612761 (delta vs Cox -0.001253).

### XGBoost_Cox (optional, disabled by default)
- **Status:** Not run in the default pipeline. Gated behind environment variable
  `ENABLE_XGBOOST_SURVIVAL=1`, because the local native XGBoost library
  repeatedly terminated the process without a Python traceback (recorded in the
  metadata note and QC report).
- **Search space when enabled (4 configs):** `n_estimators in {100}` x
  `learning_rate in {0.03, 0.1}` x `max_depth in {1, 3}` x `subsample in {0.8}`,
  with `colsample_bytree = 0.9` (objective `survival:cox`).
