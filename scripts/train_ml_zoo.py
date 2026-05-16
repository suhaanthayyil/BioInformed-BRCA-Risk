#!/usr/bin/env python3
"""Phase 5: train the survival ML model zoo on TCGA only."""

from __future__ import annotations

import json
import os
import pickle
import sys
from copy import deepcopy
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import duckdb
import joblib
import numpy as np
import pandas as pd
import torch
import xgboost as xgb
from lifelines import CoxPHFitter
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import GradientBoostingSurvivalAnalysis, RandomSurvivalForest
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.ml.deepsurv import DeepSurvConfig, DeepSurvEstimator  # noqa: E402
from src.survival import brier_score_at, time_dependent_auc  # noqa: E402


DATA = REPO_ROOT / "data"
PROCESSED = DATA / "processed"
RESULTS = REPO_ROOT / "results"
REPORTS = RESULTS / "reports"
MODELS = REPO_ROOT / "models"
LOGS = REPO_ROOT / "logs"
SEED = 42
N_SPLITS = 5

PATHWAY_FEATURES = [
    "Pathway_Immune",
    "Pathway_Proliferation",
    "Pathway_DNA_Repair",
    "Pathway_Metabolism",
    "Pathway_Stromal_EMT",
    "Pathway_Apoptosis_Stress",
    "Pathway_Hormone",
]


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "agent_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def stage_to_ordinal(value) -> float:
    text = "" if pd.isna(value) else str(value).upper()
    if "IV" in text:
        return 4.0
    if "III" in text:
        return 3.0
    if "II" in text:
        return 2.0
    if "I" in text:
        return 1.0
    return np.nan


def load_modeling_frame() -> pd.DataFrame:
    features = pd.read_parquet(PROCESSED / "04_pathway_features.parquet")
    con = duckdb.connect(str(PROCESSED / "unified_cohorts.duckdb"), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
    finally:
        con.close()
    frame = (
        features.merge(samples, on=["sample_id", "cohort"], how="inner")
        .merge(survival, on="sample_id", how="inner")
    )
    frame["stage_ordinal"] = frame["stage_tnm"].map(stage_to_ordinal)
    frame["age_at_dx"] = pd.to_numeric(frame["age_at_dx"], errors="coerce")
    return frame


def feature_columns(frame: pd.DataFrame) -> list[str]:
    clinical = ["age_at_dx", "stage_ordinal"]
    return [c for c in [*PATHWAY_FEATURES, *clinical] if c in frame.columns]


def make_preprocessor() -> Pipeline:
    return Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())])


def cindex(time, event, risk) -> float:
    if len(time) < 5 or np.sum(event) < 2:
        return np.nan
    return float(concordance_index_censored(event.astype(bool), time.astype(float), risk.astype(float))[0])


class LifelinesCoxWrapper:
    def __init__(self, penalizer: float):
        self.penalizer = penalizer
        self.model = CoxPHFitter(penalizer=penalizer)

    def fit(self, x, time, event):
        df = pd.DataFrame(x, columns=[f"x{i}" for i in range(x.shape[1])])
        df["T"] = time
        df["E"] = event.astype(int)
        self.model.fit(df, duration_col="T", event_col="E", show_progress=False)
        return self

    def predict(self, x):
        df = pd.DataFrame(x, columns=[f"x{i}" for i in range(x.shape[1])])
        return self.model.predict_partial_hazard(df).to_numpy(float)


class XGBCoxWrapper:
    def __init__(self, **params):
        self.params = params
        self.model = xgb.XGBRegressor(
            objective="survival:cox",
            eval_metric="cox-nloglik",
            random_state=SEED,
            tree_method="hist",
            n_jobs=1,
            verbosity=0,
            **params,
        )

    def fit(self, x, time, event):
        signed_time = np.where(event.astype(bool), time, -time)
        self.model.fit(x, signed_time, verbose=False)
        return self

    def predict(self, x):
        return self.model.predict(x)


class SkSurvAdapter:
    """Pickleable adapter giving sksurv estimators a common fit signature."""

    def __init__(self, estimator):
        self.inner = estimator

    def fit(self, x, time, event):
        self.inner.fit(x, Surv.from_arrays(event.astype(bool), time.astype(float)))
        return self

    def predict(self, x):
        return np.asarray(self.inner.predict(x), dtype=float)


def evaluate_config(model_name, make_model, config, x_df, time, event, folds):
    oof = np.full(len(x_df), np.nan, dtype=float)
    fold_scores = []
    for fold, (train_idx, test_idx) in enumerate(folds, start=1):
        pre = make_preprocessor()
        x_train = pre.fit_transform(x_df.iloc[train_idx])
        x_test = pre.transform(x_df.iloc[test_idx])
        t_train, e_train = time[train_idx], event[train_idx]
        t_test, e_test = time[test_idx], event[test_idx]
        try:
            model = make_model(config)
            if model_name == "DeepSurv":
                model.fit(x_train, t_train, e_train, x_test, t_test, e_test)
            else:
                model.fit(x_train, t_train, e_train)
            risk = np.asarray(model.predict(x_test), dtype=float)
            oof[test_idx] = risk
            fold_scores.append(cindex(t_test, e_test, risk))
        except Exception as exc:
            log(f"Phase 5 {model_name} config failed on fold {fold}: {config} :: {exc}")
            return None
    return {
        "model": model_name,
        "config": config,
        "mean_cindex": float(np.nanmean(fold_scores)),
        "sd_cindex": float(np.nanstd(fold_scores, ddof=1)),
        "fold_cindexes": [float(x) for x in fold_scores],
        "oof_risk": oof,
    }


def fit_final_model(model_name, make_model, config, x_df, time, event):
    pre = make_preprocessor()
    x = pre.fit_transform(x_df)
    model = make_model(config)
    model.fit(x, time, event)
    return {"model_name": model_name, "config": config, "preprocessor": pre, "model": model}


def save_model_artifact(model_name: str, artifact: dict) -> str:
    MODELS.mkdir(exist_ok=True)
    safe = model_name.lower().replace(" ", "_").replace("-", "_")
    if model_name == "DeepSurv":
        path = MODELS / "deepsurv.pt"
        torch.save(
            {
                "config": artifact["config"],
                "preprocessor": artifact["preprocessor"],
                "state_dict": artifact["model"].state_dict(),
                "features": artifact["features"],
                "input_dim": len(artifact["features"]),
            },
            path,
        )
        return str(path.relative_to(REPO_ROOT))
    path = MODELS / f"{safe}.pkl"
    with path.open("wb") as fh:
        pickle.dump(artifact, fh)
    return str(path.relative_to(REPO_ROOT))


def model_spaces():
    spaces = {
        "Cox_PH": (
            lambda cfg: LifelinesCoxWrapper(**cfg),
            [{"penalizer": p} for p in [0.01, 0.1, 1.0]],
        ),
        "Elastic_Net_Cox": (
            lambda cfg: CoxnetSurvivalAnalysis(
                l1_ratio=cfg["l1_ratio"],
                alphas=[cfg["alpha"]],
                fit_baseline_model=True,
                max_iter=100000,
            ),
            [{"alpha": a, "l1_ratio": l} for a in [0.001, 0.01, 0.1, 1.0] for l in [0.1, 0.5, 0.9]],
        ),
        "Random_Survival_Forest": (
            lambda cfg: RandomSurvivalForest(
                n_estimators=cfg["n_estimators"],
                max_depth=cfg["max_depth"],
                min_samples_leaf=cfg["min_samples_leaf"],
                max_features=cfg["max_features"],
                random_state=SEED,
                n_jobs=2,
            ),
            [
                {"n_estimators": n, "max_depth": d, "min_samples_leaf": leaf, "max_features": mf}
                for n in [500]
                for d in [None, 4, 8]
                for leaf in [1, 5, 10]
                for mf in ["sqrt", None]
            ],
        ),
        "Gradient_Boosted_Survival": (
            lambda cfg: GradientBoostingSurvivalAnalysis(random_state=SEED, **cfg),
            [
                {"n_estimators": n, "learning_rate": lr, "max_depth": d, "subsample": ss}
                for n in [100, 300]
                for lr in [0.03, 0.1]
                for d in [1, 3]
                for ss in [0.8, 1.0]
            ],
        ),
        "DeepSurv": (
            lambda cfg: DeepSurvEstimator(DeepSurvConfig(**cfg)),
            [
                {"hidden_dims": (32,), "dropout": 0.1, "lr": 1e-3, "weight_decay": 1e-4, "epochs": 120, "patience": 15},
                {"hidden_dims": (64, 32), "dropout": 0.1, "lr": 1e-3, "weight_decay": 1e-4, "epochs": 150, "patience": 20},
                {"hidden_dims": (64, 32), "dropout": 0.2, "lr": 5e-4, "weight_decay": 1e-4, "epochs": 150, "patience": 20},
                {"hidden_dims": (128, 64), "dropout": 0.2, "lr": 5e-4, "weight_decay": 1e-4, "epochs": 150, "patience": 20},
            ],
        ),
    }
    if os.environ.get("ENABLE_XGBOOST_SURVIVAL") == "1":
        spaces["XGBoost_Cox"] = (
            lambda cfg: XGBCoxWrapper(**cfg),
            [
                {
                    "n_estimators": n,
                    "learning_rate": lr,
                    "max_depth": d,
                    "subsample": ss,
                    "colsample_bytree": 0.9,
                }
                for n in [100]
                for lr in [0.03, 0.1]
                for d in [1, 3]
                for ss in [0.8]
            ],
        )
    return spaces


def train_stacked_ensemble(results, artifacts, x_df, time, event):
    ranked = sorted(results, key=lambda r: r["mean_cindex"], reverse=True)
    top = ranked[:3]
    oof = pd.DataFrame({item["model"]: item["oof_risk"] for item in top})
    df = oof.copy()
    df["T"] = time
    df["E"] = event.astype(int)
    meta = CoxPHFitter(penalizer=0.1)
    meta.fit(df, duration_col="T", event_col="E")
    meta_risk = meta.predict_partial_hazard(oof).to_numpy(float)
    artifact = {
        "model_name": "Stacked_Ensemble",
        "top_models": [item["model"] for item in top],
        "meta_model": meta,
        "features": list(x_df.columns),
    }
    path = MODELS / "stacked_ensemble.pkl"
    with path.open("wb") as fh:
        pickle.dump(artifact, fh)
    return {
        "model": "Stacked_Ensemble",
        "config": {"top_models": [item["model"] for item in top]},
        "mean_cindex": cindex(time, event, meta_risk),
        "sd_cindex": np.nan,
        "fold_cindexes": [],
        "oof_risk": meta_risk,
        "artifact_path": str(path.relative_to(REPO_ROOT)),
    }


def main() -> None:
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    log("Phase 5 ML model zoo training started.")

    frame = load_modeling_frame()
    train = frame[frame["cohort"].eq("TCGA-BRCA")].copy()
    train = train[train["os_days"].notna() & train["os_event"].notna() & (train["os_days"] > 0)]
    features = feature_columns(train)
    x_df = train[features].copy()
    time = pd.to_numeric(train["os_days"], errors="coerce").to_numpy(float)
    event = pd.to_numeric(train["os_event"], errors="coerce").fillna(0).to_numpy(int).astype(bool)
    y_surv = Surv.from_arrays(event, time)
    folds = list(StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED).split(x_df, event.astype(int)))

    results = []
    artifacts = {}
    for model_name, (make_model, configs) in model_spaces().items():
        log(f"Phase 5 training {model_name} with {len(configs)} compact-grid configs.")
        candidates = []
        for config in configs:
            # scikit-survival models consume structured y, handled by adapter below.
            if model_name in {"Elastic_Net_Cox", "Random_Survival_Forest", "Gradient_Boosted_Survival"}:
                def maker(cfg, base=make_model):
                    return SkSurvAdapter(base(cfg))

                candidate = evaluate_config(model_name, maker, config, x_df, time, event, folds)
            else:
                candidate = evaluate_config(model_name, make_model, config, x_df, time, event, folds)
            if candidate is not None and np.isfinite(candidate["mean_cindex"]):
                candidates.append(candidate)
        if not candidates:
            log(f"Phase 5 no valid candidates for {model_name}.")
            continue
        best = max(candidates, key=lambda item: item["mean_cindex"])
        results.append(best)

        if model_name in {"Elastic_Net_Cox", "Random_Survival_Forest", "Gradient_Boosted_Survival"}:
            def final_maker(cfg, base=make_model):
                return SkSurvAdapter(base(cfg))

            artifact = fit_final_model(model_name, final_maker, best["config"], x_df, time, event)
        else:
            artifact = fit_final_model(model_name, make_model, best["config"], x_df, time, event)
        artifact["features"] = features
        artifact["training_cohort"] = "TCGA-BRCA"
        artifact["cv_mean_cindex"] = best["mean_cindex"]
        artifact_path = save_model_artifact(model_name, artifact)
        best["artifact_path"] = artifact_path
        artifacts[model_name] = artifact

    if len(results) >= 3:
        stack = train_stacked_ensemble(results, artifacts, x_df, time, event)
        results.append(stack)

    rows = []
    for item in results:
        risk = item["oof_risk"]
        rows.append(
            {
                "model": item["model"],
                "mean_cv_cindex": item["mean_cindex"],
                "sd_cv_cindex": item["sd_cindex"],
                "fold_cindexes": json.dumps(item["fold_cindexes"]),
                "config": json.dumps(item["config"], default=list),
                "artifact_path": item.get("artifact_path", ""),
                **time_dependent_auc(time, event, risk, years=(5, 10)),
                **brier_score_at(time, event, risk, years=(5,)),
            }
        )

    table = pd.DataFrame(rows).sort_values("mean_cv_cindex", ascending=False)
    cox = table.loc[table["model"].eq("Cox_PH"), "mean_cv_cindex"]
    cox_value = float(cox.iloc[0]) if len(cox) else np.nan
    table["delta_vs_cox"] = table["mean_cv_cindex"] - cox_value
    headline = table.iloc[0].to_dict()

    RESULTS.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    table.to_csv(RESULTS / "Table_S2_ml_internal_cv.csv", index=False)
    pd.DataFrame(
        {
            "sample_id": train["sample_id"].to_numpy(),
            "time": time,
            "event": event.astype(int),
            **{item["model"]: item["oof_risk"] for item in results},
        }
    ).to_csv(RESULTS / "05_ml_oof_predictions.csv", index=False)

    metadata = {
        "generated_at": timestamp(),
        "seed": SEED,
        "n_splits": N_SPLITS,
        "training_cohort": "TCGA-BRCA",
        "n_training_samples": int(len(train)),
        "n_events": int(event.sum()),
        "features": features,
        "grid_mode": "compact_local",
        "note": (
            "Compact grid used for local tractability; XGBoost Cox disabled by default "
            "after repeated native process termination. Set ENABLE_XGBOOST_SURVIVAL=1 "
            "to retry that optional cross-check."
        ),
        "headline_model": headline["model"],
        "headline_cv_cindex": headline["mean_cv_cindex"],
        "headline_delta_vs_cox": headline["delta_vs_cox"],
    }
    (PROCESSED / "ml_model_zoo.metadata.json").write_text(json.dumps(metadata, indent=2))
    report = [
        "# Phase 5 ML Model Zoo QC",
        "",
        f"Generated: {timestamp()}",
        "",
        f"Training cohort: TCGA-BRCA, n={len(train)}, events={int(event.sum())}",
        f"Features: {', '.join(features)}",
        "",
        "## Internal 5-fold CV",
        "",
        table.to_markdown(index=False, floatfmt=".3f"),
        "",
        f"Headline model by TCGA CV: {headline['model']} (C-index {headline['mean_cv_cindex']:.3f}, delta vs Cox {headline['delta_vs_cox']:.3f}).",
        "",
        "Grid mode: compact local grid. This is explicitly recorded and should not be described as the full protocol grid.",
        "XGBoost Cox is disabled by default because the local native library repeatedly terminated the process without a Python traceback.",
        "",
    ]
    (REPORTS / "ml_model_zoo_qc.md").write_text("\n".join(report))
    log(f"Phase 5 ML model zoo completed. Headline={headline['model']}.")


if __name__ == "__main__":
    main()
