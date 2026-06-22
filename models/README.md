# Trained model artifacts

Trained survival-model weights for the BRCA-PathwayML benchmark, provided for
reproducibility (BMC Artificial Intelligence revision, editor requirement E3).
All models were trained on **TCGA-BRCA only** (n = 213 samples, 132 events,
seed = 42, 5-fold internal stratified CV) using the nine features below, then
locked before external validation.

```
Features (ordered): Pathway_Immune, Pathway_Proliferation, Pathway_DNA_Repair,
Pathway_Metabolism, Pathway_Stromal_EMT, Pathway_Apoptosis_Stress,
Pathway_Hormone, age_at_dx, stage_ordinal
```

The **headline model** is `gradient_boosted_survival.pkl` (Gradient Boosted
Survival Analysis), internal CV Harrell C = 0.642, Δ vs Cox = +0.042.

## Loading

The artifacts load cold (no import shim) — the wrapper classes live in
`src/ml/wrappers.py`:

```bash
python scripts/predict.py --from-features data/processed/04_pathway_features.parquet --cohort TCGA-BRCA
python scripts/predict.py --model random_survival_forest --input my_features.csv -o risk.csv
```

Each `.pkl` unpickles to a dict: `model_name`, `config`, `preprocessor`
(sklearn `SimpleImputer(median)` → `StandardScaler`, fit on TCGA-BRCA only —
**no leakage**: external cohorts are only `.transform`-ed), `model`, `features`,
`training_cohort`, `cv_mean_cindex`. `deepsurv.pt` is a `torch` checkpoint
(`config`, `preprocessor`, `state_dict`, `features`, `input_dim`); rebuild with
`CoxMLP` from `src/ml/deepsurv.py` (see `scripts/predict.py`).

If a fresh clone ever fails to load a pickle, run `python scripts/migrate_pickles.py` once.

## Inventory (SHA256 over the committed artifacts)

| File | Manuscript model | Internal CV C-index | Selected config | SHA256 |
|------|------------------|--------------------:|-----------------|--------|
| `gradient_boosted_survival.pkl` | Gradient Boosted Survival **(headline)** | 0.6418 | n_estimators=100, lr=0.03, max_depth=1, subsample=1.0 | `ef89e35d…931a7d` |
| `stacked_ensemble.pkl` | Stacked Ensemble (GBSA+DeepSurv+ElasticNet) | 0.6263 | meta-Cox over top-3 | `46987bbc…d39411` |
| `deepsurv.pt` | DeepSurv | 0.6090 | hidden=[128,64], dropout=0.2, lr=5e-4, wd=1e-4, epochs=150 | `f8d06703…b29124` |
| `elastic_net_cox.pkl` | Elastic-Net Cox | 0.6013 | alpha=0.001, l1_ratio=0.1 | `80d21ed1…7af061` |
| `random_survival_forest.pkl` | Random Survival Forest | 0.6007 | n_estimators=500, max_depth=4, max_features=sqrt | `325767e2…0dcedf` |
| `cox_ph.pkl` | Cox PH baseline | 0.5999 | penalizer=0.01 | `b2ee6d96…9728d8` |

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
post-hoc rescue exploration that `docs/STORY.md` documents as over-fit; it is
not a result of the manuscript and must not be used as a prognostic model.
