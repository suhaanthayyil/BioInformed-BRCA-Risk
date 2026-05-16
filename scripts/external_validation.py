#!/usr/bin/env python3
"""Phase 6 external validation and pre-registered TNBC primary endpoint."""

from __future__ import annotations

import json
import pickle
import sys
from datetime import datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import torch
from lifelines import CoxPHFitter

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.train_ml_zoo import (  # noqa: E402,F401
    LifelinesCoxWrapper,
    SkSurvAdapter,
    XGBCoxWrapper,
    feature_columns,
    stage_to_ordinal,
)
from src.meta import random_effects_from_ci  # noqa: E402
from src.ml.deepsurv import CoxMLP  # noqa: E402
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
BOOTSTRAP_ITER = 1000
PAM50_COMPARATOR = "PAM50_ROR_official"


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    LOGS.mkdir(exist_ok=True)
    with (LOGS / "agent_log.md").open("a") as fh:
        fh.write(f"\n- {timestamp()} {message}\n")


def load_frame() -> pd.DataFrame:
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
    return frame


def load_pickle(path: Path):
    with path.open("rb") as fh:
        return pickle.load(fh)


def load_model_artifacts() -> dict[str, dict]:
    artifacts = {}
    for name, filename in {
        "Cox_PH": "cox_ph.pkl",
        "Elastic_Net_Cox": "elastic_net_cox.pkl",
        "Random_Survival_Forest": "random_survival_forest.pkl",
        "Gradient_Boosted_Survival": "gradient_boosted_survival.pkl",
    }.items():
        artifacts[name] = load_pickle(MODELS / filename)
    payload = torch.load(MODELS / "deepsurv.pt", map_location="cpu", weights_only=False)
    config = payload["config"]
    model = CoxMLP(payload["input_dim"], tuple(config["hidden_dims"]), config["dropout"])
    model.load_state_dict(payload["state_dict"])
    model.eval()
    artifacts["DeepSurv"] = {
        "model_name": "DeepSurv",
        "config": config,
        "preprocessor": payload["preprocessor"],
        "model": model,
        "features": payload["features"],
    }
    artifacts["Stacked_Ensemble"] = load_pickle(MODELS / "stacked_ensemble.pkl")
    return artifacts


def predict_artifact(name: str, artifact: dict, x: pd.DataFrame, base_predictions: dict[str, np.ndarray] | None = None) -> np.ndarray:
    if name == "Stacked_Ensemble":
        top = artifact["top_models"]
        if base_predictions is None:
            raise ValueError("Stacked ensemble requires base predictions")
        meta_x = pd.DataFrame({model: base_predictions[model] for model in top})
        return artifact["meta_model"].predict_partial_hazard(meta_x).to_numpy(float)

    pre = artifact["preprocessor"]
    x_scaled = pre.transform(x[artifact["features"]])
    if name == "DeepSurv":
        with torch.no_grad():
            risk = artifact["model"](torch.tensor(x_scaled, dtype=torch.float32)).detach().numpy()
        return risk.astype(float)
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


def add_combined_model(predictions: pd.DataFrame, headline: str) -> pd.DataFrame:
    pam_all = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam = (
        pam_all[pam_all["baseline"].eq(PAM50_COMPARATOR) & pam_all["status"].eq("ok")][["sample_id", "score"]]
        .rename(columns={"score": "pam50_score"})
    )
    merged = predictions.merge(pam, on="sample_id", how="left")
    train = merged[merged["cohort"].eq("TCGA-BRCA") & merged[headline].notna() & merged["pam50_score"].notna()].copy()
    if len(train) < 20 or train["os_event"].sum() < 3:
        merged["Combined_ML_PAM50"] = np.nan
        return merged
    df = train[[headline, "pam50_score", "os_days", "os_event"]].rename(columns={"os_days": "T", "os_event": "E"})
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(df, "T", "E")
    valid = merged[headline].notna() & merged["pam50_score"].notna()
    merged.loc[valid, "Combined_ML_PAM50"] = cph.predict_partial_hazard(merged.loc[valid, [headline, "pam50_score"]]).to_numpy(float)
    return merged


def main() -> None:
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    log("Phase 6 external validation started.")
    metadata = json.loads((PROCESSED / "ml_model_zoo.metadata.json").read_text())
    headline = metadata["headline_model"]
    frame = load_frame()
    features = metadata["features"]
    artifacts = load_model_artifacts()

    rows = []
    for cohort, cohort_df in frame.groupby("cohort", sort=True):
        x = cohort_df[features]
        base_predictions = {}
        cohort_pred = cohort_df[
            [
                "sample_id",
                "cohort",
                "os_days",
                "os_event",
                "tnbc_flag",
                "intrinsic_subtype_pam50_published",
                "stage_tnm",
                "age_at_dx",
            ]
        ].copy()
        for name, artifact in artifacts.items():
            if name == "Stacked_Ensemble":
                continue
            risk = predict_artifact(name, artifact, x)
            base_predictions[name] = risk
            cohort_pred[name] = risk
        if "Stacked_Ensemble" in artifacts:
            cohort_pred["Stacked_Ensemble"] = predict_artifact("Stacked_Ensemble", artifacts["Stacked_Ensemble"], x, base_predictions)
        rows.append(cohort_pred)

    predictions = pd.concat(rows, ignore_index=True)
    pam50_all = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam50 = (
        pam50_all[pam50_all["baseline"].eq(PAM50_COMPARATOR) & pam50_all["status"].eq("ok")][["sample_id", "score"]]
        .rename(columns={"score": PAM50_COMPARATOR})
    )
    predictions = predictions.merge(pam50, on="sample_id", how="left")
    predictions = add_combined_model(predictions, headline)
    predictions.to_csv(RESULTS / "06_external_model_predictions.csv", index=False)

    model_cols = [
        "Cox_PH",
        "Elastic_Net_Cox",
        "Random_Survival_Forest",
        "Gradient_Boosted_Survival",
        "DeepSurv",
        "Stacked_Ensemble",
        "Combined_ML_PAM50",
        PAM50_COMPARATOR,
    ]
    metric_rows = []
    for cohort, cohort_df in predictions.groupby("cohort", sort=True):
        tnbc_flag = cohort_df["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
        masks = {
            "overall": pd.Series(True, index=cohort_df.index),
            "tnbc": tnbc_flag,
            "non_tnbc": ~tnbc_flag,
        }
        for model in model_cols:
            if model not in cohort_df:
                continue
            for subgroup, mask in masks.items():
                metric_rows.append(metric_row(cohort, model, subgroup, cohort_df.loc[mask], model))

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(RESULTS / "Table_1_ml_external_validation.csv", index=False)

    h2h_rows = []
    comparators = [PAM50_COMPARATOR, "Cox_PH", "Combined_ML_PAM50"]
    for cohort, cohort_df in predictions.groupby("cohort", sort=True):
        tnbc_flag = cohort_df["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
        for subgroup, mask in {
            "overall": pd.Series(True, index=cohort_df.index),
            "tnbc": tnbc_flag,
            "non_tnbc": ~tnbc_flag,
        }.items():
            sub = cohort_df.loc[mask].copy()
            for comparator in comparators:
                if comparator not in sub or headline not in sub:
                    continue
                valid = sub[[headline, comparator, "os_days", "os_event"]].dropna()
                result = paired_bootstrap_delta_c_index(
                    valid["os_days"],
                    valid["os_event"],
                    valid[headline],
                    valid[comparator],
                    n_bootstrap=BOOTSTRAP_ITER,
                    seed=SEED,
                )
                h2h_rows.append(
                    {
                        "cohort": cohort,
                        "subgroup": subgroup,
                        "headline_model": headline,
                        "comparator": comparator,
                        "n": len(valid),
                        "events": int(valid["os_event"].sum()) if len(valid) else 0,
                        "delta_cindex": result["delta"],
                        "ci_low": result["ci_low"],
                        "ci_high": result["ci_high"],
                        "p": result["p"],
                    }
                )

    h2h = pd.DataFrame(h2h_rows)
    h2h.to_csv(RESULTS / "Table_2_head_to_head.csv", index=False)
    primary_rows = h2h[
        h2h["subgroup"].eq("tnbc")
        & h2h["comparator"].eq(PAM50_COMPARATOR)
        & h2h["delta_cindex"].notna()
        & h2h["ci_low"].notna()
        & h2h["ci_high"].notna()
        & (h2h["n"] >= 10)
    ].copy()
    meta = random_effects_from_ci(primary_rows, "delta_cindex", "ci_low", "ci_high")
    primary_met = bool(meta["effect"] >= 0.03 and meta["p"] < 0.05)
    primary = {
        "generated_at": timestamp(),
        "headline_model": headline,
        "primary_endpoint": f"TNBC random-effects delta C-index vs {PAM50_COMPARATOR}",
        "pam50_comparator": PAM50_COMPARATOR,
        "threshold_delta": 0.03,
        "threshold_p": 0.05,
        "met": primary_met,
        **{f"meta_{k}": v for k, v in meta.items()},
        "cohorts_in_meta": primary_rows["cohort"].tolist(),
    }
    (PROCESSED / "phase6_primary_endpoint.json").write_text(json.dumps(primary, indent=2))

    REPORTS.mkdir(exist_ok=True)
    report_lines = [
        "# Phase 6 External Validation QC",
        "",
        f"Generated: {timestamp()}",
        f"Headline model: {headline}",
        "",
        "## External Validation Metrics",
        "",
        metrics[
            (metrics["subgroup"].eq("overall"))
            & (metrics["model"].isin([headline, "Cox_PH", PAM50_COMPARATOR]))
        ][["cohort", "model", "n", "events", "harrell_c", "harrell_c_ci_low", "harrell_c_ci_high", "status"]].to_markdown(
            index=False, floatfmt=".3f"
        ),
        "",
        "## Primary Endpoint Rows",
        "",
        primary_rows.to_markdown(index=False, floatfmt=".3f"),
        "",
        "## Random-Effects Primary Endpoint",
        "",
        pd.DataFrame([primary]).to_markdown(index=False, floatfmt=".4f"),
        "",
    ]
    (REPORTS / "external_validation_qc.md").write_text("\n".join(report_lines))
    log(f"Phase 6 external validation completed. Primary endpoint met={primary_met}.")


if __name__ == "__main__":
    main()
