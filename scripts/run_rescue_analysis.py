#!/usr/bin/env python3
"""Transparent rescue analysis after the locked primary endpoint failed.

This script deliberately keeps the pre-registered Phase 6 result intact. It
adds an exploratory model/feature search selected on TCGA only, then evaluates
the selected rescue model against official genefu PAM50-ROR on held-out cohorts.
"""

from __future__ import annotations

import json
import os
import pickle
import sys
import warnings
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import duckdb
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

from src.meta import random_effects_from_ci  # noqa: E402
from src.ml.deepsurv import DeepSurvConfig, DeepSurvEstimator  # noqa: E402
from src.pathways import SEVEN_PATHWAY_COMPONENTS  # noqa: E402
from src.survival import (  # noqa: E402
    brier_score_at,
    bootstrap_c_index_ci,
    cox_hr_high_low,
    harrell_c_index,
    paired_bootstrap_delta_c_index,
    time_dependent_auc,
    uno_c_index,
)


PROCESSED = REPO_ROOT / "data" / "processed"
RESULTS = REPO_ROOT / "results"
REPORTS = RESULTS / "reports"
MODELS = REPO_ROOT / "models"
LOGS = REPO_ROOT / "logs"

SEED = 42
N_SPLITS = 5
BOOTSTRAP_ITER = 1000
PAM50_COMPARATOR = "PAM50_ROR_official"
TRAINING_COHORT = "TCGA-BRCA"
PRIMARY_COHORTS = ["GSE96058", "METABRIC", "TCGA-BRCA"]
EXTERNAL_ONLY_COHORTS = ["GSE96058", "METABRIC"]

SEVEN_FEATURES = [
    "Pathway_Immune",
    "Pathway_Proliferation",
    "Pathway_DNA_Repair",
    "Pathway_Metabolism",
    "Pathway_Stromal_EMT",
    "Pathway_Apoptosis_Stress",
    "Pathway_Hormone",
]
CLINICAL_FEATURES = ["age_at_dx", "stage_ordinal"]


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "agent_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def stage_to_ordinal(value: object) -> float:
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


def status_positive(value: object) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text in {"positive", "pos", "+", "1", "true"} or "positive" in text:
        return 1.0
    if text in {"negative", "neg", "-", "0", "false"} or "negative" in text:
        return 0.0
    return np.nan


def make_preprocessor() -> Pipeline:
    return Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())])


def cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    valid = np.isfinite(time) & np.isfinite(risk) & (time > 0)
    if valid.sum() < 5 or event[valid].sum() < 2:
        return np.nan
    return float(concordance_index_censored(event[valid].astype(bool), time[valid].astype(float), risk[valid])[0])


def zscore_by_cohort(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for _, idx in out.groupby("cohort", sort=False).groups.items():
        sub = out.loc[idx, columns].apply(pd.to_numeric, errors="coerce")
        mean = sub.mean(axis=0)
        std = sub.std(axis=0, ddof=0).replace(0, np.nan)
        out.loc[idx, columns] = ((sub - mean) / std).fillna(0.0)
    return out


def load_base_frame() -> pd.DataFrame:
    scores = pd.read_parquet(PROCESSED / "pathway_scores_all.parquet")
    seven = pd.read_parquet(PROCESSED / "04_pathway_features.parquet")
    pathway_cols = [c for c in scores.columns if c not in {"sample_id", "cohort"}]
    scores = zscore_by_cohort(scores, pathway_cols)
    frame = seven.merge(scores, on=["sample_id", "cohort"], how="inner")

    con = duckdb.connect(str(PROCESSED / "unified_cohorts.duckdb"), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
    finally:
        con.close()

    frame = frame.merge(samples, on=["sample_id", "cohort"], how="inner").merge(survival, on="sample_id", how="inner")
    frame["age_at_dx"] = pd.to_numeric(frame["age_at_dx"], errors="coerce")
    frame["stage_ordinal"] = frame["stage_tnm"].map(stage_to_ordinal)
    frame["er_positive"] = frame["er_status"].map(status_positive)
    frame["pr_positive"] = frame["pr_status"].map(status_positive)
    frame["her2_positive"] = frame["her2_status"].map(status_positive)
    frame["tnbc_numeric"] = frame["tnbc_flag"].astype("boolean").astype(float)
    return add_interactions(frame)


def add_interactions(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["Immune_minus_Proliferation"] = out["Pathway_Immune"] - out["Pathway_Proliferation"]
    out["Immune_x_Proliferation"] = out["Pathway_Immune"] * out["Pathway_Proliferation"]
    out["Proliferation_x_Hormone"] = out["Pathway_Proliferation"] * out["Pathway_Hormone"]
    out["Stage_x_Proliferation"] = out["stage_ordinal"] * out["Pathway_Proliferation"]
    out["Age_x_Immune"] = out["age_at_dx"] * out["Pathway_Immune"]
    out["TNBC_x_Immune"] = out["tnbc_numeric"] * out["Pathway_Immune"]
    out["TNBC_x_Proliferation"] = out["tnbc_numeric"] * out["Pathway_Proliferation"]
    out["TNBC_x_Stromal_EMT"] = out["tnbc_numeric"] * out["Pathway_Stromal_EMT"]
    out["HER2_x_Immune"] = out["her2_positive"] * out["Pathway_Immune"]
    out["ER_x_Hormone"] = out["er_positive"] * out["Pathway_Hormone"]
    return out


def feature_sets(frame: pd.DataFrame) -> dict[str, list[str]]:
    hallmark = [c for c in frame.columns if c.startswith("HALLMARK_")]
    components = sorted({c for values in SEVEN_PATHWAY_COMPONENTS.values() for c in values if c in frame.columns})
    core_interactions = [
        "Immune_minus_Proliferation",
        "Immune_x_Proliferation",
        "Proliferation_x_Hormone",
        "Stage_x_Proliferation",
        "Age_x_Immune",
    ]
    receptor = [
        "er_positive",
        "pr_positive",
        "her2_positive",
        "tnbc_numeric",
        "TNBC_x_Immune",
        "TNBC_x_Proliferation",
        "TNBC_x_Stromal_EMT",
        "HER2_x_Immune",
        "ER_x_Hormone",
    ]
    return {
        "seven": [*SEVEN_FEATURES, *CLINICAL_FEATURES],
        "seven_interactions": [*SEVEN_FEATURES, *core_interactions, *CLINICAL_FEATURES],
        "seven_receptor_interactions": [*SEVEN_FEATURES, *core_interactions, *receptor, *CLINICAL_FEATURES],
        "hallmark50": [*hallmark, *CLINICAL_FEATURES],
        "hallmark50_interactions": [*hallmark, *core_interactions, *CLINICAL_FEATURES],
        "component7_interactions": [*SEVEN_FEATURES, *components, *core_interactions, *CLINICAL_FEATURES],
    }


def usable_features(train: pd.DataFrame, columns: list[str]) -> list[str]:
    usable = []
    for col in columns:
        if col not in train.columns:
            continue
        values = pd.to_numeric(train[col], errors="coerce")
        if values.notna().sum() < 10:
            continue
        if values.nunique(dropna=True) <= 1:
            continue
        usable.append(col)
    return usable


class LifelinesCoxWrapper:
    def __init__(self, penalizer: float):
        self.penalizer = penalizer
        self.model = CoxPHFitter(penalizer=penalizer)

    def fit(self, x: np.ndarray, time: np.ndarray, event: np.ndarray):
        df = pd.DataFrame(x, columns=[f"x{i}" for i in range(x.shape[1])])
        df["T"] = time
        df["E"] = event.astype(int)
        self.model.fit(df, duration_col="T", event_col="E", show_progress=False)
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        df = pd.DataFrame(x, columns=[f"x{i}" for i in range(x.shape[1])])
        return self.model.predict_partial_hazard(df).to_numpy(float)


class SkSurvAdapter:
    def __init__(self, estimator):
        self.inner = estimator

    def fit(self, x: np.ndarray, time: np.ndarray, event: np.ndarray):
        self.inner.fit(x, Surv.from_arrays(event.astype(bool), time.astype(float)))
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(self.inner.predict(x), dtype=float)


class XGBAFTWrapper:
    def __init__(self, **params):
        self.n_estimators = params.pop("n_estimators")
        self.params = {
            "objective": "survival:aft",
            "eval_metric": "aft-nloglik",
            "aft_loss_distribution": "normal",
            "aft_loss_distribution_scale": params.pop("aft_loss_distribution_scale", 1.0),
            "tree_method": "hist",
            "seed": SEED,
            "nthread": 1,
            **params,
        }
        self.booster = None

    def fit(self, x: np.ndarray, time: np.ndarray, event: np.ndarray):
        dtrain = xgb.DMatrix(x)
        lower = time.astype(float)
        upper = np.where(event.astype(bool), time.astype(float), np.inf)
        dtrain.set_float_info("label_lower_bound", lower)
        dtrain.set_float_info("label_upper_bound", upper)
        self.booster = xgb.train(self.params, dtrain, num_boost_round=self.n_estimators, verbose_eval=False)
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        if self.booster is None:
            raise RuntimeError("XGBAFTWrapper must be fit before predict")
        # AFT predicts log survival time; invert so higher score means higher risk.
        return -np.asarray(self.booster.predict(xgb.DMatrix(x)), dtype=float)


def model_spaces() -> dict[str, tuple[Any, list[dict[str, Any]]]]:
    def with_weights(configs: list[dict[str, Any]], weights: tuple[int, ...] = (1, 2, 4)) -> list[dict[str, Any]]:
        return [{**cfg, "tnbc_weight": weight} for cfg in configs for weight in weights]

    spaces = {
        "Cox_PH": (
            lambda cfg: LifelinesCoxWrapper(penalizer=cfg["penalizer"]),
            with_weights([{"penalizer": p} for p in [0.01, 0.1, 1.0]]),
        ),
        "Elastic_Net_Cox": (
            lambda cfg: SkSurvAdapter(
                CoxnetSurvivalAnalysis(
                    l1_ratio=cfg["l1_ratio"],
                    alphas=[cfg["alpha"]],
                    fit_baseline_model=True,
                    max_iter=100000,
                )
            ),
            with_weights(
                [{"alpha": alpha, "l1_ratio": ratio} for alpha in [0.001, 0.01, 0.1] for ratio in [0.1, 0.5]]
            ),
        ),
        "Random_Survival_Forest": (
            lambda cfg: SkSurvAdapter(
                RandomSurvivalForest(
                    n_estimators=cfg["n_estimators"],
                    max_depth=cfg["max_depth"],
                    min_samples_leaf=cfg["min_samples_leaf"],
                    max_features=cfg["max_features"],
                    random_state=SEED,
                    n_jobs=2,
                )
            ),
            with_weights(
                [
                    {"n_estimators": 300, "max_depth": d, "min_samples_leaf": leaf, "max_features": "sqrt"}
                    for d in [4, 8, None]
                    for leaf in [3, 10]
                ],
                weights=(1, 2),
            ),
        ),
        "Gradient_Boosted_Survival": (
            lambda cfg: SkSurvAdapter(
                GradientBoostingSurvivalAnalysis(
                    random_state=SEED,
                    n_estimators=cfg["n_estimators"],
                    learning_rate=cfg["learning_rate"],
                    max_depth=cfg["max_depth"],
                    subsample=cfg["subsample"],
                )
            ),
            with_weights(
                [
                    {"n_estimators": n, "learning_rate": lr, "max_depth": d, "subsample": ss}
                    for n in [100, 300]
                    for lr in [0.03, 0.1]
                    for d in [1, 2]
                    for ss in [0.8]
                ]
            ),
        ),
        "DeepSurv": (
            lambda cfg: DeepSurvEstimator(DeepSurvConfig(**cfg)),
            [
                {
                    "hidden_dims": (64, 32),
                    "dropout": 0.1,
                    "lr": 1e-3,
                    "weight_decay": 1e-4,
                    "epochs": 120,
                    "patience": 15,
                    "tnbc_weight": 1,
                },
                {
                    "hidden_dims": (128, 64),
                    "dropout": 0.2,
                    "lr": 5e-4,
                    "weight_decay": 1e-4,
                    "epochs": 150,
                    "patience": 20,
                    "tnbc_weight": 2,
                },
            ],
        ),
    }
    if os.environ.get("ENABLE_RESCUE_XGBOOST") == "1":
        spaces["XGBoost_AFT"] = (
            lambda cfg: XGBAFTWrapper(
                n_estimators=cfg["n_estimators"],
                learning_rate=cfg["learning_rate"],
                max_depth=cfg["max_depth"],
                subsample=cfg["subsample"],
                colsample_bynode=cfg["colsample_bynode"],
                aft_loss_distribution_scale=cfg["aft_loss_distribution_scale"],
            ),
            with_weights(
                [
                    {
                        "n_estimators": n,
                        "learning_rate": lr,
                        "max_depth": d,
                        "subsample": 0.9,
                        "colsample_bynode": 0.8,
                        "aft_loss_distribution_scale": 1.0,
                    }
                    for n in [100, 300]
                    for lr in [0.03, 0.1]
                    for d in [1, 2]
                ],
                weights=(1, 2),
            ),
        )
    return spaces


def strip_training_config(config: dict[str, Any]) -> dict[str, Any]:
    out = dict(config)
    out.pop("tnbc_weight", None)
    return out


def weighted_training_arrays(
    x: pd.DataFrame,
    time: np.ndarray,
    event: np.ndarray,
    tnbc: np.ndarray,
    weight: int,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    if weight <= 1:
        return x, time, event
    repeats = np.where(tnbc.astype(bool), weight, 1)
    idx = np.repeat(np.arange(len(x)), repeats)
    return x.iloc[idx].reset_index(drop=True), time[idx], event[idx]


@dataclass
class CandidateResult:
    feature_set: str
    model: str
    config: dict[str, Any]
    features: list[str]
    overall_cv_cindex: float
    tnbc_cv_cindex: float
    valid_tnbc_folds: int
    composite_cv_cindex: float
    fold_cindexes: list[float]
    fold_tnbc_cindexes: list[float]
    oof_risk: np.ndarray

    @property
    def candidate_id(self) -> str:
        return f"{self.feature_set}__{self.model}__{abs(hash(json.dumps(self.config, sort_keys=True, default=str))) % 10**8}"


def evaluate_candidate(
    feature_set: str,
    model_name: str,
    make_model,
    config: dict[str, Any],
    train: pd.DataFrame,
    features: list[str],
    folds: list[tuple[np.ndarray, np.ndarray]],
) -> CandidateResult | None:
    time = pd.to_numeric(train["os_days"], errors="coerce").to_numpy(float)
    event = pd.to_numeric(train["os_event"], errors="coerce").fillna(0).to_numpy(int).astype(bool)
    tnbc = train["tnbc_flag"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    x_df = train[features].apply(pd.to_numeric, errors="coerce")
    oof = np.full(len(train), np.nan, dtype=float)
    fold_scores = []
    fold_tnbc_scores = []
    tnbc_weight = int(config.get("tnbc_weight", 1))
    model_config = strip_training_config(config)

    for fold, (train_idx, test_idx) in enumerate(folds, start=1):
        try:
            x_train_raw, t_train, e_train = weighted_training_arrays(
                x_df.iloc[train_idx],
                time[train_idx],
                event[train_idx],
                tnbc[train_idx],
                tnbc_weight,
            )
            pre = make_preprocessor()
            x_train = pre.fit_transform(x_train_raw)
            x_test = pre.transform(x_df.iloc[test_idx])
            model = make_model(model_config)
            if model_name == "DeepSurv":
                model.fit(x_train, t_train, e_train, x_test, time[test_idx], event[test_idx])
            else:
                model.fit(x_train, t_train, e_train)
            risk = np.asarray(model.predict(x_test), dtype=float)
            oof[test_idx] = risk
            fold_scores.append(cindex(time[test_idx], event[test_idx], risk))
            tnbc_mask = tnbc[test_idx]
            fold_tnbc_scores.append(cindex(time[test_idx][tnbc_mask], event[test_idx][tnbc_mask], risk[tnbc_mask]))
        except Exception as exc:
            log(f"Rescue candidate failed on fold {fold}: {feature_set}/{model_name}/{config} :: {exc}")
            return None

    overall = float(np.nanmean(fold_scores))
    tnbc_valid = [x for x in fold_tnbc_scores if np.isfinite(x)]
    tnbc_score = float(np.nanmean(tnbc_valid)) if tnbc_valid else np.nan
    composite = 0.5 * overall + 0.5 * tnbc_score if len(tnbc_valid) >= 2 else overall
    if not np.isfinite(composite):
        return None
    return CandidateResult(
        feature_set=feature_set,
        model=model_name,
        config=deepcopy(config),
        features=features,
        overall_cv_cindex=overall,
        tnbc_cv_cindex=tnbc_score,
        valid_tnbc_folds=len(tnbc_valid),
        composite_cv_cindex=float(composite),
        fold_cindexes=[float(x) if np.isfinite(x) else np.nan for x in fold_scores],
        fold_tnbc_cindexes=[float(x) if np.isfinite(x) else np.nan for x in fold_tnbc_scores],
        oof_risk=oof,
    )


def rank_standardize(values: np.ndarray) -> np.ndarray:
    series = pd.Series(values)
    return series.rank(method="average", pct=True).to_numpy(float)


def ensemble_from_top(results: list[CandidateResult], train: pd.DataFrame) -> CandidateResult | None:
    top = sorted(results, key=lambda r: r.composite_cv_cindex, reverse=True)[:3]
    if len(top) < 2:
        return None
    risks = np.column_stack([rank_standardize(item.oof_risk) for item in top])
    oof = np.nanmean(risks, axis=1)
    time = pd.to_numeric(train["os_days"], errors="coerce").to_numpy(float)
    event = pd.to_numeric(train["os_event"], errors="coerce").fillna(0).to_numpy(int).astype(bool)
    tnbc = train["tnbc_flag"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    overall = cindex(time, event, oof)
    tnbc_score = cindex(time[tnbc], event[tnbc], oof[tnbc])
    composite = 0.5 * overall + 0.5 * tnbc_score if np.isfinite(tnbc_score) else overall
    config = {"members": [item.candidate_id for item in top], "member_labels": [f"{i.feature_set}:{i.model}" for i in top]}
    return CandidateResult(
        feature_set="ensemble_top3",
        model="Rank_Average_Ensemble",
        config=config,
        features=[],
        overall_cv_cindex=float(overall),
        tnbc_cv_cindex=float(tnbc_score),
        valid_tnbc_folds=1 if np.isfinite(tnbc_score) else 0,
        composite_cv_cindex=float(composite),
        fold_cindexes=[],
        fold_tnbc_cindexes=[],
        oof_risk=oof,
    )


def fit_final_candidate(result: CandidateResult, train: pd.DataFrame, spaces: dict[str, tuple[Any, list[dict[str, Any]]]]) -> dict:
    if result.model == "Rank_Average_Ensemble":
        raise ValueError("Fit ensemble through fit_final_ensemble")
    make_model = spaces[result.model][0]
    x_df = train[result.features].apply(pd.to_numeric, errors="coerce")
    time = pd.to_numeric(train["os_days"], errors="coerce").to_numpy(float)
    event = pd.to_numeric(train["os_event"], errors="coerce").fillna(0).to_numpy(int).astype(bool)
    tnbc = train["tnbc_flag"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    x_weighted, time_weighted, event_weighted = weighted_training_arrays(
        x_df,
        time,
        event,
        tnbc,
        int(result.config.get("tnbc_weight", 1)),
    )
    pre = make_preprocessor()
    x = pre.fit_transform(x_weighted)
    model = make_model(strip_training_config(result.config))
    if result.model == "DeepSurv":
        model.fit(x, time_weighted, event_weighted, x, time_weighted, event_weighted)
    else:
        model.fit(x, time_weighted, event_weighted)
    return {
        "candidate": result,
        "preprocessor": pre,
        "model": model,
        "features": result.features,
    }


def fit_final_ensemble(
    ensemble: CandidateResult,
    lookup: dict[str, CandidateResult],
    train: pd.DataFrame,
    spaces: dict[str, tuple[Any, list[dict[str, Any]]]],
) -> dict:
    members = []
    for member_id in ensemble.config["members"]:
        members.append(fit_final_candidate(lookup[member_id], train, spaces))
    return {"candidate": ensemble, "members": members, "features": []}


def predict_artifact(artifact: dict, frame: pd.DataFrame) -> np.ndarray:
    candidate: CandidateResult = artifact["candidate"]
    if candidate.model == "Rank_Average_Ensemble":
        member_risks = np.column_stack([rank_standardize(predict_artifact(member, frame)) for member in artifact["members"]])
        return np.nanmean(member_risks, axis=1)
    x = frame[artifact["features"]].apply(pd.to_numeric, errors="coerce")
    x_scaled = artifact["preprocessor"].transform(x)
    return np.asarray(artifact["model"].predict(x_scaled), dtype=float)


def metric_row(cohort: str, model: str, subgroup: str, sub: pd.DataFrame, risk_col: str) -> dict:
    sub = sub[sub[risk_col].notna() & sub["os_days"].notna() & sub["os_event"].notna() & (sub["os_days"] > 0)]
    n = len(sub)
    events = int(pd.to_numeric(sub["os_event"], errors="coerce").fillna(0).sum()) if n else 0
    if n < 10 or events < 2:
        return {
            "cohort": cohort,
            "model": model,
            "endpoint": "OS",
            "subgroup": subgroup,
            "n": n,
            "events": events,
            "harrell_c": np.nan,
            "harrell_c_ci_low": np.nan,
            "harrell_c_ci_high": np.nan,
            "uno_c": np.nan,
            "auc_3y": np.nan,
            "auc_5y": np.nan,
            "auc_10y": np.nan,
            "brier_5y": np.nan,
            "hr": np.nan,
            "hr_ci_low": np.nan,
            "hr_ci_high": np.nan,
            "hr_p": np.nan,
            "status": "too_few_events",
        }
    ci_low, ci_high = bootstrap_c_index_ci(sub["os_days"], sub["os_event"], sub[risk_col], BOOTSTRAP_ITER, SEED)
    return {
        "cohort": cohort,
        "model": model,
        "endpoint": "OS",
        "subgroup": subgroup,
        "n": n,
        "events": events,
        "harrell_c": harrell_c_index(sub["os_days"], sub["os_event"], sub[risk_col]),
        "harrell_c_ci_low": ci_low,
        "harrell_c_ci_high": ci_high,
        "uno_c": uno_c_index(sub["os_days"], sub["os_event"], sub[risk_col]),
        **time_dependent_auc(sub["os_days"], sub["os_event"], sub[risk_col], years=(3, 5, 10)),
        **brier_score_at(sub["os_days"], sub["os_event"], sub[risk_col], years=(5,)),
        **cox_hr_high_low(sub["os_days"], sub["os_event"], sub[risk_col]),
        "status": "ok",
    }


def add_pam50(predictions: pd.DataFrame) -> pd.DataFrame:
    pam = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam = (
        pam[pam["baseline"].eq(PAM50_COMPARATOR) & pam["status"].eq("ok")][["sample_id", "score"]]
        .rename(columns={"score": PAM50_COMPARATOR})
    )
    return predictions.merge(pam, on="sample_id", how="left")


def evaluate_external(frame: pd.DataFrame, artifacts: dict[str, dict], headline_name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cols = [
        "sample_id",
        "cohort",
        "os_days",
        "os_event",
        "tnbc_flag",
        "intrinsic_subtype_pam50_published",
        "stage_tnm",
        "age_at_dx",
    ]
    predictions = frame[cols].copy()
    for label, artifact in artifacts.items():
        predictions[label] = predict_artifact(artifact, frame)
    predictions = add_pam50(predictions)
    predictions.to_csv(RESULTS / "07_rescue_external_predictions.csv", index=False)

    metric_rows = []
    model_cols = [*artifacts.keys(), PAM50_COMPARATOR]
    for cohort, cohort_df in predictions.groupby("cohort", sort=True):
        tnbc_flag = cohort_df["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
        masks = {
            "overall": pd.Series(True, index=cohort_df.index),
            "tnbc": tnbc_flag,
            "non_tnbc": ~tnbc_flag,
        }
        for model in model_cols:
            for subgroup, mask in masks.items():
                metric_rows.append(metric_row(cohort, model, subgroup, cohort_df.loc[mask], model))
    metrics = pd.DataFrame(metric_rows)

    h2h_rows = []
    for cohort, cohort_df in predictions.groupby("cohort", sort=True):
        tnbc_flag = cohort_df["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
        masks = {
            "overall": pd.Series(True, index=cohort_df.index),
            "tnbc": tnbc_flag,
            "non_tnbc": ~tnbc_flag,
        }
        for subgroup, mask in masks.items():
            sub = cohort_df.loc[mask].copy()
            valid = sub[[headline_name, PAM50_COMPARATOR, "os_days", "os_event"]].dropna()
            result = paired_bootstrap_delta_c_index(
                valid["os_days"],
                valid["os_event"],
                valid[headline_name],
                valid[PAM50_COMPARATOR],
                n_bootstrap=BOOTSTRAP_ITER,
                seed=SEED,
            )
            h2h_rows.append(
                {
                    "cohort": cohort,
                    "subgroup": subgroup,
                    "headline_model": headline_name,
                    "comparator": PAM50_COMPARATOR,
                    "n": len(valid),
                    "events": int(valid["os_event"].sum()) if len(valid) else 0,
                    "delta_cindex": result["delta"],
                    "ci_low": result["ci_low"],
                    "ci_high": result["ci_high"],
                    "p": result["p"],
                }
            )
    h2h = pd.DataFrame(h2h_rows)
    return predictions, metrics, h2h


def primary_summary(h2h: pd.DataFrame, cohorts: list[str], label: str) -> dict[str, Any]:
    rows = h2h[
        h2h["cohort"].isin(cohorts)
        & h2h["subgroup"].eq("tnbc")
        & h2h["delta_cindex"].notna()
        & h2h["ci_low"].notna()
        & h2h["ci_high"].notna()
        & (h2h["n"] >= 10)
    ].copy()
    meta = random_effects_from_ci(rows, "delta_cindex", "ci_low", "ci_high")
    return {
        "label": label,
        "cohorts": rows["cohort"].tolist(),
        "met": bool(meta["effect"] >= 0.03 and meta["p"] < 0.05),
        **{f"meta_{key}": value for key, value in meta.items()},
    }


def metabric_expression_audit() -> dict[str, Any]:
    raw = REPO_ROOT / "data" / "raw" / "metabric" / "data_mrna_illumina_microarray.txt"
    zscores = REPO_ROOT / "data" / "raw" / "metabric" / "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt"
    audit = {
        "current_processed_expression": "cBioPortal z-scores relative to diploid samples",
        "raw_expression_file_present": raw.exists(),
        "zscore_expression_file_present": zscores.exists(),
    }
    for label, path in [("raw", raw), ("zscore", zscores)]:
        if path.exists():
            header = pd.read_csv(path, sep="\t", nrows=5)
            values = header.iloc[:, 2:].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            audit[f"{label}_shape_header_check"] = [int(header.shape[0]), int(header.shape[1])]
            audit[f"{label}_first_rows_median"] = float(np.nanmedian(values))
            audit[f"{label}_first_rows_sd"] = float(np.nanstd(values))
    return audit


def result_table(results: list[CandidateResult]) -> pd.DataFrame:
    rows = []
    for item in results:
        rows.append(
            {
                "candidate_id": item.candidate_id,
                "feature_set": item.feature_set,
                "model": item.model,
                "n_features": len(item.features),
                "config": json.dumps(item.config, default=list, sort_keys=True),
                "overall_cv_cindex": item.overall_cv_cindex,
                "tnbc_cv_cindex": item.tnbc_cv_cindex,
                "valid_tnbc_folds": item.valid_tnbc_folds,
                "composite_cv_cindex": item.composite_cv_cindex,
                "fold_cindexes": json.dumps(item.fold_cindexes),
                "fold_tnbc_cindexes": json.dumps(item.fold_tnbc_cindexes),
            }
        )
    return pd.DataFrame(rows).sort_values("composite_cv_cindex", ascending=False)


def write_report(
    table: pd.DataFrame,
    metrics: pd.DataFrame,
    h2h: pd.DataFrame,
    primary: dict[str, Any],
    external_only: dict[str, Any],
    audit: dict[str, Any],
    headline_name: str,
) -> None:
    REPORTS.mkdir(exist_ok=True)
    top_cols = [
        "feature_set",
        "model",
        "n_features",
        "overall_cv_cindex",
        "tnbc_cv_cindex",
        "composite_cv_cindex",
        "valid_tnbc_folds",
    ]
    lines = [
        "# Rescue Analysis QC",
        "",
        f"Generated: {timestamp()}",
        "",
        "This is an exploratory transparent rescue run. The locked Phase 6 primary endpoint remains unchanged.",
        "",
        "## METABRIC Expression Audit",
        "",
        pd.DataFrame([audit]).to_markdown(index=False),
        "",
        "## Top TCGA-Selected Candidates",
        "",
        table.head(15)[top_cols].to_markdown(index=False, floatfmt=".3f"),
        "",
        "## Rescue Headline",
        "",
        f"Headline rescue model: `{headline_name}` selected by TCGA-only composite CV.",
        "",
        "## External TNBC Rows vs Official PAM50",
        "",
        h2h[h2h["subgroup"].eq("tnbc")].to_markdown(index=False, floatfmt=".3f"),
        "",
        "## Random-Effects Meta-Analysis",
        "",
        pd.DataFrame([primary, external_only]).to_markdown(index=False, floatfmt=".4f"),
        "",
        "## External Validation Metrics for Headline and PAM50",
        "",
        metrics[
            metrics["model"].isin([headline_name, PAM50_COMPARATOR]) & metrics["subgroup"].isin(["overall", "tnbc"])
        ][["cohort", "model", "subgroup", "n", "events", "harrell_c", "harrell_c_ci_low", "harrell_c_ci_high", "status"]]
        .to_markdown(index=False, floatfmt=".3f"),
        "",
    ]
    (REPORTS / "rescue_analysis_qc.md").write_text("\n".join(lines))


def main() -> None:
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    warnings.filterwarnings("ignore", category=UserWarning)
    log("Rescue analysis started.")

    frame = load_base_frame()
    train = frame[
        frame["cohort"].eq(TRAINING_COHORT)
        & frame["os_days"].notna()
        & frame["os_event"].notna()
        & (frame["os_days"] > 0)
    ].copy()
    train_event = pd.to_numeric(train["os_event"], errors="coerce").fillna(0).to_numpy(int)
    folds = list(StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED).split(train, train_event))
    sets = feature_sets(frame)
    spaces = model_spaces()
    if "XGBoost_AFT" not in spaces:
        log("Rescue XGBoost_AFT disabled by default because native survival runs have terminated locally.")

    results: list[CandidateResult] = []
    for set_name, columns in sets.items():
        features = usable_features(train, columns)
        log(f"Rescue search feature set {set_name}: {len(features)} usable features.")
        for model_name, (make_model, configs) in spaces.items():
            if model_name == "DeepSurv" and len(features) > 80:
                continue
            for config in configs:
                candidate = evaluate_candidate(set_name, model_name, make_model, config, train, features, folds)
                if candidate is not None:
                    results.append(candidate)
            partial = result_table(results)
            partial.to_csv(RESULTS / "Table_S4_rescue_internal_cv.csv", index=False)

    ensemble = ensemble_from_top(results, train)
    if ensemble is not None:
        results.append(ensemble)

    table = result_table(results)
    RESULTS.mkdir(exist_ok=True)
    MODELS.mkdir(exist_ok=True)
    table.to_csv(RESULTS / "Table_S4_rescue_internal_cv.csv", index=False)
    pd.DataFrame(
        {
            "sample_id": train["sample_id"].to_numpy(),
            "time": train["os_days"].to_numpy(float),
            "event": train_event,
            **{item.candidate_id: item.oof_risk for item in results},
        }
    ).to_csv(RESULTS / "07_rescue_oof_predictions.csv", index=False)

    lookup = {item.candidate_id: item for item in results}
    headline = results[int(np.nanargmax([item.composite_cv_cindex for item in results]))]
    headline_name = f"Rescue_{headline.model}"
    if headline.model == "Rank_Average_Ensemble":
        headline_artifact = fit_final_ensemble(headline, lookup, train, spaces)
    else:
        headline_artifact = fit_final_candidate(headline, train, spaces)
    artifacts = {headline_name: headline_artifact}

    with (MODELS / "rescue_headline.pkl").open("wb") as fh:
        pickle.dump(
            {
                "headline_name": headline_name,
                "candidate": table.iloc[0].to_dict(),
                "note": "Exploratory rescue artifact; selection used TCGA only.",
            },
            fh,
        )

    _, metrics, h2h = evaluate_external(frame, artifacts, headline_name)
    metrics.to_csv(RESULTS / "Table_S5_rescue_external_validation.csv", index=False)
    h2h.to_csv(RESULTS / "Table_S6_rescue_head_to_head.csv", index=False)

    primary = primary_summary(h2h, PRIMARY_COHORTS, "all_protocol_tnbc_cohorts")
    external_only = primary_summary(h2h, EXTERNAL_ONLY_COHORTS, "external_only_tnbc_cohorts")
    audit = metabric_expression_audit()
    summary = {
        "generated_at": timestamp(),
        "seed": SEED,
        "selection_policy": "TCGA-only composite CV: 0.5 overall C-index + 0.5 TNBC C-index when TNBC folds evaluable",
        "training_cohort": TRAINING_COHORT,
        "n_training_samples": int(len(train)),
        "n_training_events": int(train_event.sum()),
        "n_candidates": int(len(results)),
        "headline_name": headline_name,
        "headline_candidate": table.iloc[0].to_dict(),
        "primary_all_protocol_tnbc": primary,
        "primary_external_only_tnbc": external_only,
        "metabric_expression_audit": audit,
        "status": "exploratory_rescue_succeeded" if primary["met"] else "exploratory_rescue_not_met",
    }
    (PROCESSED / "rescue_analysis.metadata.json").write_text(json.dumps(summary, indent=2, default=str))
    (PROCESSED / "rescue_primary_endpoint.json").write_text(json.dumps(primary, indent=2, default=str))
    write_report(table, metrics, h2h, primary, external_only, audit, headline_name)
    log(f"Rescue analysis completed. Headline={headline_name}; all-cohort primary met={primary['met']}.")


if __name__ == "__main__":
    main()
