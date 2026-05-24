"""Shared helpers for Phase Three BMC Artificial Intelligence submission preparation."""

from __future__ import annotations

import json
import pickle
import sys
import __main__
from datetime import datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.train_ml_zoo import SkSurvAdapter, stage_to_ordinal  # noqa: E402,F401
from src.meta import random_effects_from_ci  # noqa: E402
from src.survival import (  # noqa: E402
    bootstrap_c_index_ci,
    harrell_c_index,
    paired_bootstrap_delta_c_index,
)

DATA = REPO_ROOT / "data"
PROCESSED = DATA / "processed"
RESULTS = REPO_ROOT / "results"
REPORTS = RESULTS / "reports"
FIGURES = REPO_ROOT / "figures"
MANUSCRIPT = REPO_ROOT / "manuscript"
MODELS = REPO_ROOT / "models"
DOCS = REPO_ROOT / "docs"
LOGS = REPO_ROOT / "logs"

SEED = 42
BOOTSTRAP_ITER = 1000
HEADLINE_MODEL = "Gradient_Boosted_Survival"
PAM50_MODEL = "PAM50_ROR_official"
COX_MODEL = "Cox_PH"
EXTERNAL_COHORTS = ["GSE96058", "METABRIC", "GSE20685"]
ALL_COHORTS = ["TCGA-BRCA", "GSE96058", "METABRIC", "GSE20685"]

# Phase 5 pickles were created by running scripts/train_ml_zoo.py as __main__.
# This alias lets pickle resolve the saved adapter class without touching scores.
setattr(__main__, "SkSurvAdapter", SkSurvAdapter)


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log_phase3(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "phase3_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def ensure_dirs() -> None:
    for path in [
        RESULTS,
        REPORTS,
        FIGURES,
        MANUSCRIPT / "figures",
        MANUSCRIPT / "tables",
        MANUSCRIPT / "additional_files",
        DOCS,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def load_samples_survival() -> pd.DataFrame:
    con = duckdb.connect(str(PROCESSED / "unified_cohorts.duckdb"), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
        qc = con.execute("select * from qc").fetchdf()
    finally:
        con.close()
    out = samples.merge(survival, on="sample_id", how="left").merge(qc, on="sample_id", how="left")
    out["age_at_dx"] = pd.to_numeric(out["age_at_dx"], errors="coerce")
    out["os_days"] = pd.to_numeric(out["os_days"], errors="coerce")
    out["os_event"] = pd.to_numeric(out["os_event"], errors="coerce").fillna(0).astype(int)
    out["stage_ordinal"] = out["stage_tnm"].map(stage_to_ordinal)
    return out


def load_predictions() -> pd.DataFrame:
    pred = pd.read_csv(RESULTS / "06_external_model_predictions.csv")
    pam = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam = pam[pam["baseline"].eq(PAM50_MODEL) & pam["status"].eq("ok")][
        ["sample_id", "pam50_subtype_genefu", "rorS_risk", "n_pam50_genes_input"]
    ]
    pred = pred.merge(pam, on="sample_id", how="left")
    pred["os_days"] = pd.to_numeric(pred["os_days"], errors="coerce")
    pred["os_event"] = pd.to_numeric(pred["os_event"], errors="coerce").fillna(0).astype(int)
    pred["tnbc_flag"] = pred["tnbc_flag"].astype("boolean")
    return pred


def load_modeling_frame() -> pd.DataFrame:
    features = pd.read_parquet(PROCESSED / "04_pathway_features.parquet")
    frame = features.merge(load_samples_survival(), on=["sample_id", "cohort"], how="inner")
    return frame


def load_headline_artifact() -> dict:
    with (MODELS / "gradient_boosted_survival.pkl").open("rb") as fh:
        return pickle.load(fh)


def predict_headline(frame: pd.DataFrame, artifact: dict | None = None) -> np.ndarray:
    if artifact is None:
        artifact = load_headline_artifact()
    x = frame[artifact["features"]]
    x_scaled = artifact["preprocessor"].transform(x)
    return np.asarray(artifact["model"].predict(x_scaled), dtype=float)


def valid_survival_subset(df: pd.DataFrame, risk_cols: list[str]) -> pd.DataFrame:
    cols = ["os_days", "os_event", *risk_cols]
    out = df.dropna(subset=cols).copy()
    out = out[out["os_days"] > 0]
    return out


def metric_cindex(df: pd.DataFrame, risk_col: str) -> tuple[float, float, float]:
    cval = harrell_c_index(df["os_days"], df["os_event"], df[risk_col])
    low, high = bootstrap_c_index_ci(df["os_days"], df["os_event"], df[risk_col], BOOTSTRAP_ITER, SEED)
    return cval, low, high


def paired_delta(df: pd.DataFrame, model_col: str, comparator_col: str) -> dict:
    return paired_bootstrap_delta_c_index(
        df["os_days"],
        df["os_event"],
        df[model_col],
        df[comparator_col],
        n_bootstrap=BOOTSTRAP_ITER,
        seed=SEED,
    )


def meta_from_rows(rows: pd.DataFrame) -> dict:
    if rows.empty:
        return {
            "n_studies": 0,
            "effect": np.nan,
            "se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "z": np.nan,
            "p": np.nan,
            "tau2": np.nan,
            "i2": np.nan,
        }
    return random_effects_from_ci(rows, "delta_vs_pam50", "delta_ci_low", "delta_ci_high")


def binary_horizon_outcome(df: pd.DataFrame, horizon_days: int) -> pd.DataFrame:
    out = df.copy()
    event_before = (out["os_event"].astype(int).eq(1)) & (out["os_days"] <= horizon_days)
    censored_before = (out["os_event"].astype(int).eq(0)) & (out["os_days"] < horizon_days)
    out = out[~censored_before].copy()
    out["horizon_event"] = event_before.loc[out.index].astype(int)
    return out


def ci(values: list[float] | np.ndarray) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return np.nan, np.nan
    return tuple(np.quantile(arr, [0.025, 0.975]).astype(float))


def spearman_pairwise(rank_matrix: np.ndarray) -> np.ndarray:
    values = []
    for i in range(rank_matrix.shape[0]):
        for j in range(i + 1, rank_matrix.shape[0]):
            rho = stats.spearmanr(rank_matrix[i], rank_matrix[j], nan_policy="omit").statistic
            if np.isfinite(rho):
                values.append(float(rho))
    return np.asarray(values, dtype=float)


def save_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, default=str))
