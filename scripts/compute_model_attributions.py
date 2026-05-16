#!/usr/bin/env python3
"""Phase 5 model attribution artifacts."""

from __future__ import annotations

import pickle
import sys
from datetime import datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import shap
import torch
from captum.attr import IntegratedGradients
from sksurv.metrics import concordance_index_censored

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.ml.deepsurv import CoxMLP  # noqa: E402
from scripts.train_ml_zoo import SkSurvAdapter, feature_columns, stage_to_ordinal  # noqa: E402,F401


PROCESSED = REPO_ROOT / "data" / "processed"
RESULTS = REPO_ROOT / "results"
MODELS = REPO_ROOT / "models"
LOGS = REPO_ROOT / "logs"
SEED = 42


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "agent_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def load_training_data():
    features = pd.read_parquet(PROCESSED / "04_pathway_features.parquet")
    con = duckdb.connect(str(PROCESSED / "unified_cohorts.duckdb"), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
    finally:
        con.close()
    frame = features.merge(samples, on=["sample_id", "cohort"], how="inner").merge(survival, on="sample_id", how="inner")
    frame["stage_ordinal"] = frame["stage_tnm"].map(stage_to_ordinal)
    frame["age_at_dx"] = pd.to_numeric(frame["age_at_dx"], errors="coerce")
    train = frame[frame["cohort"].eq("TCGA-BRCA") & frame["os_days"].notna() & frame["os_event"].notna()].copy()
    cols = feature_columns(train)
    return train, cols


def cindex(time, event, risk):
    return float(concordance_index_censored(event.astype(bool), time.astype(float), risk.astype(float))[0])


def load_pickle(path: Path):
    with path.open("rb") as fh:
        return pickle.load(fh)


def rsf_permutation_importance(train: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    artifact = load_pickle(MODELS / "random_survival_forest.pkl")
    x = train[cols]
    x_scaled = artifact["preprocessor"].transform(x)
    model = artifact["model"]
    base = cindex(train["os_days"].to_numpy(float), train["os_event"].to_numpy(bool), model.predict(x_scaled))
    rng = np.random.default_rng(SEED)
    rows = []
    for idx, feature in enumerate(cols):
        drops = []
        for _ in range(50):
            x_perm = x_scaled.copy()
            x_perm[:, idx] = rng.permutation(x_perm[:, idx])
            perm = cindex(train["os_days"].to_numpy(float), train["os_event"].to_numpy(bool), model.predict(x_perm))
            drops.append(base - perm)
        rows.append(
            {
                "feature": feature,
                "base_cindex": base,
                "importance_mean_delta_cindex": float(np.mean(drops)),
                "importance_sd": float(np.std(drops, ddof=1)),
            }
        )
    return pd.DataFrame(rows).sort_values("importance_mean_delta_cindex", ascending=False)


def gbsa_shap(train: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    artifact = load_pickle(MODELS / "gradient_boosted_survival.pkl")
    x = train[cols]
    x_scaled = artifact["preprocessor"].transform(x)
    model = artifact["model"]
    background = shap.sample(pd.DataFrame(x_scaled, columns=cols), min(50, len(x_scaled)), random_state=SEED)
    explainer = shap.Explainer(lambda z: model.predict(np.asarray(z, dtype=float)), background, algorithm="permutation")
    values = explainer(pd.DataFrame(x_scaled, columns=cols), max_evals=2 * len(cols) + 1)
    rows = []
    for sample_id, arr in zip(train["sample_id"].astype(str), values.values):
        for feature, value in zip(cols, arr):
            rows.append({"sample_id": sample_id, "feature": feature, "shap_value": float(value)})
    return pd.DataFrame(rows)


def deepsurv_integrated_gradients(train: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    payload = torch.load(MODELS / "deepsurv.pt", map_location="cpu", weights_only=False)
    config = payload["config"]
    hidden_dims = tuple(config["hidden_dims"])
    model = CoxMLP(payload["input_dim"], hidden_dims, config["dropout"])
    model.load_state_dict(payload["state_dict"])
    model.eval()

    x_scaled = payload["preprocessor"].transform(train[cols]).astype("float32")
    x_t = torch.tensor(x_scaled, dtype=torch.float32)
    baseline = torch.zeros_like(x_t)
    ig = IntegratedGradients(model)
    attrs = ig.attribute(x_t, baselines=baseline, n_steps=32).detach().numpy()
    rows = []
    for sample_id, arr in zip(train["sample_id"].astype(str), attrs):
        for feature, value in zip(cols, arr):
            rows.append({"sample_id": sample_id, "feature": feature, "integrated_gradient": float(value)})
    return pd.DataFrame(rows)


def main() -> None:
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    log("Phase 5 model attribution started.")
    train, cols = load_training_data()
    rsf = rsf_permutation_importance(train, cols)
    rsf.to_csv(RESULTS / "05_rsf_feature_importance.csv", index=False)
    gbsa = gbsa_shap(train, cols)
    gbsa.to_csv(RESULTS / "05_gbsa_shap_values.csv", index=False)
    ig = deepsurv_integrated_gradients(train, cols)
    ig.to_csv(RESULTS / "05_deepsurv_integrated_gradients.csv", index=False)
    log("Phase 5 model attribution completed.")


if __name__ == "__main__":
    main()
