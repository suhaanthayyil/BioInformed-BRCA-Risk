# Trained model artifacts

Trained survival-model weights for the BRCA-PathwayML benchmark. All models were
trained on **TCGA-BRCA only** (n = 1,064 patients,
150 deaths, seed = 42, 5-fold stratified internal CV) using the
seven registered pathway features:

```
Pathway_Immune, Pathway_Proliferation, Pathway_DNA_Repair, Pathway_Metabolism, Pathway_Stromal_EMT, Pathway_Apoptosis_Stress, Pathway_Hormone
```

The **headline model** is the best registered non-linear candidate (RSF, GBSA or
DeepSurv) by internal CV: `deepsurv.pt` (DeepSurv).

`locked_v1/` holds, unchanged, the weights locked on 2026-05-16 (TCGA n = 213,
nine features including a stage covariate that was missing for every external
patient), the pathway features and genefu scores they were evaluated on, and
the harmonized database of that version. The weights are used only by
`scripts/locked_model_sensitivity.py`; the data copies are kept for reference.

## Loading

The artifacts load cold (no import shim) — the wrapper classes live in
`src/ml/wrappers.py`:

```bash
python scripts/predict.py --from-features data/processed/04_pathway_features.parquet --cohort TCGA-BRCA
python scripts/predict.py --model random_survival_forest --input my_features.csv -o risk.csv
```

Each `.pkl` is a compressed joblib file (`joblib.load`, or `src/ml/registry.load_artifact`) holding a dict: `model_name`, `config`, `preprocessor`
(sklearn `SimpleImputer(median)` → `StandardScaler`, fit on TCGA-BRCA only —
**no leakage**: external cohorts are only `.transform`-ed), `model`, `features`,
`training_cohort`, `cv_mean_cindex`. `deepsurv.pt` is a `torch` checkpoint
(`config`, `preprocessor`, `state_dict`, `features`, `input_dim`); rebuild with
`CoxMLP` from `src/ml/deepsurv.py` (see `scripts/predict.py`).

If a fresh clone ever fails to load a pickle, run `python scripts/migrate_pickles.py` once.

## Inventory (SHA256 over the committed artifacts)

| File | Model | Internal CV C-index (5x5 repeated; ensemble 5-fold) | Selected config | SHA256 |
|------|-------|--------------------:|-----------------|--------|
| `cox_ph.pkl` | Cox PH | 0.6140 | `{"penalizer": 0.1}` | `b00c90d8…b65fb4` |
| `elastic_net_cox.pkl` | Elastic Net Cox | 0.6138 | `{"alpha": 0.01, "l1_ratio": 0.1}` | `ab9748b1…b479fa` |
| `stacked_ensemble.pkl` | Stacked Ensemble | 0.6128 | `{"top_models": ["Cox_PH", "Elastic_Net_Cox", "DeepSurv"]}` | `ea647938…27b847` |
| `deepsurv.pt` | DeepSurv **(headline)** | 0.5930 | `{"hidden_dims": [128, 64], "dropout": 0.3, "lr": 0.0003, "weight_decay": 0.0001, "epochs": 300, "patience": 30}` | `d0a87cd7…c4c7ad` |
| `random_survival_forest.pkl` | Random Survival Forest | 0.5775 | `{"n_estimators": 500, "max_depth": 5, "min_samples_leaf": 40, "max_features": 0.5}` | `c0b31011…0df8f5` |
| `gradient_boosted_survival.pkl` | Gradient Boosted Survival | 0.5739 | `{"n_estimators": 100, "learning_rate": 0.01, "max_depth": 1, "subsample": 0.8}` | `0af5f432…527524` |

Full 64-char digests: `shasum -a 256 models/*.pkl models/*.pt` (cross-check against
`docs/MODEL_SHA256.txt`). Per-model fold C-indexes, AUC, and Brier scores are in
`results/Table_S2_ml_internal_cv.csv`.

## Reference environment

Trained/serialized under the pinned stack in `requirements.txt` — notably
`scikit-survival==0.23.0`, `scikit-learn==1.5.2`, `torch==2.4.1`,
`lifelines==0.30.0`, `xgboost==2.1.4`. The pickles also load under newer minor
versions; use the pinned versions (or the `Dockerfile`) for exact numerical
reproduction.

## Not shipped

`rescue_headline.pkl` is **intentionally excluded** (gitignored). It is a
post-hoc rescue exploration that is over-fit; it is
not a result of the manuscript and must not be used as a prognostic model.
