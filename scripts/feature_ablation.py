#!/usr/bin/env python3
"""Task 2: Feature-ablation and incremental-value analysis (Tables 6, S5).

Eight feature set combinations trained on TCGA with locked GBSA architecture,
applied to three external cohorts. Reports Harrell C, Uno C, 5-year Brier,
ICI, and paired delta C vs PAM50-ROR-only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import GradientBoostingSurvivalAnalysis
from sksurv.util import Surv

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.survival import (  # noqa: E402
    brier_score_at,
    bootstrap_c_index_ci,
    harrell_c_index,
    paired_bootstrap_delta_c_index,
    uno_c_index,
)

DATA = REPO_ROOT / "data"
PROCESSED = DATA / "processed"
RESULTS = REPO_ROOT / "results"
SEED = 42
BOOTSTRAP_ITER = 1000
N_SPLITS = 5
PAM50 = "PAM50_ROR_official"
EXTERNAL_COHORTS = ["GSE96058", "METABRIC", "GSE20685"]

PATHWAY_FEATURES = [
    "Pathway_Immune",
    "Pathway_Proliferation",
    "Pathway_DNA_Repair",
    "Pathway_Metabolism",
    "Pathway_Stromal_EMT",
    "Pathway_Apoptosis_Stress",
    "Pathway_Hormone",
]


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


# Eight feature-set combinations
FEATURE_SETS = {
    "Pathways_only": PATHWAY_FEATURES,
    "Clinical_only": ["age_at_dx", "stage_ordinal"],
    "Pathways+Clinical": PATHWAY_FEATURES + ["age_at_dx", "stage_ordinal"],
    "Proliferation_only": ["Pathway_Proliferation"],
    "Hormone_only": ["Pathway_Hormone"],
    "Immune_only": ["Pathway_Immune"],
    "Top3_pathways": ["Pathway_Hormone", "Pathway_Proliferation", "Pathway_Immune"],
    "Leave_one_out_Hormone": [p for p in PATHWAY_FEATURES if p != "Pathway_Hormone"]
        + ["age_at_dx", "stage_ordinal"],
}


def load_frame() -> pd.DataFrame:
    features = pd.read_parquet(PROCESSED / "04_pathway_features.parquet")
    con = duckdb.connect(str(PROCESSED / "unified_cohorts.duckdb"), read_only=True)
    try:
        samples = con.execute("select * from samples").fetchdf()
        survival = con.execute("select * from survival").fetchdf()
    finally:
        con.close()
    frame = features.merge(samples, on=["sample_id", "cohort"], how="inner").merge(
        survival, on="sample_id", how="inner"
    )
    frame["stage_ordinal"] = frame["stage_tnm"].map(stage_to_ordinal)
    frame["age_at_dx"] = pd.to_numeric(frame["age_at_dx"], errors="coerce")
    frame["os_days"] = pd.to_numeric(frame["os_days"], errors="coerce")
    frame["os_event"] = pd.to_numeric(frame["os_event"], errors="coerce").fillna(0).astype(int)
    return frame


def make_preprocessor() -> Pipeline:
    return Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())])


def get_locked_gbsa_config() -> dict:
    # Load the best config from the internal CV table
    cv_table = pd.read_csv(RESULTS / "Table_S2_ml_internal_cv.csv")
    gbsa_row = cv_table[cv_table["model"] == "Gradient_Boosted_Survival"].iloc[0]
    return json.loads(gbsa_row["config"])


def train_and_predict(feature_names: list[str], train_df: pd.DataFrame,
                      test_dfs: dict[str, pd.DataFrame], config: dict) -> dict[str, np.ndarray]:
    """Train GBSA on TCGA with given features, predict on external cohorts."""
    available = [f for f in feature_names if f in train_df.columns]
    if not available:
        return {c: np.full(len(df), np.nan) for c, df in test_dfs.items()}

    pre = make_preprocessor()
    x_train = pre.fit_transform(train_df[available])
    t_train = train_df["os_days"].to_numpy(float)
    e_train = train_df["os_event"].to_numpy(bool)

    y_train = Surv.from_arrays(e_train, t_train)
    gbsa = GradientBoostingSurvivalAnalysis(random_state=SEED, **config)
    gbsa.fit(x_train, y_train)

    predictions = {}
    for cohort, cdf in test_dfs.items():
        x_test = pre.transform(cdf[available])
        predictions[cohort] = np.asarray(gbsa.predict(x_test), dtype=float)
    return predictions


def compute_ici(time, event, risk) -> float:
    """Integrated Calibration Index - mean absolute calibration error."""
    time_arr = pd.to_numeric(pd.Series(time), errors="coerce").to_numpy(float)
    event_arr = pd.to_numeric(pd.Series(event), errors="coerce").fillna(0).to_numpy(int)
    risk_arr = pd.to_numeric(pd.Series(risk), errors="coerce").to_numpy(float)
    valid = np.isfinite(time_arr) & (time_arr > 0) & np.isfinite(risk_arr)
    if valid.sum() < 20 or event_arr[valid].sum() < 3:
        return np.nan

    z = (risk_arr[valid] - np.mean(risk_arr[valid])) / (np.std(risk_arr[valid]) or 1.0)
    pred_prob = 1.0 / (1.0 + np.exp(-z))

    # Use deciles of predicted risk
    try:
        deciles = pd.qcut(pred_prob, 10, duplicates="drop")
        grouped = pd.DataFrame({
            "pred": pred_prob,
            "event": event_arr[valid],
            "time": time_arr[valid],
            "group": deciles,
        }).groupby("group", observed=True)

        ici_values = []
        for _, grp in grouped:
            if len(grp) < 5:
                continue
            obs_rate = grp["event"].mean()
            pred_rate = grp["pred"].mean()
            ici_values.append(abs(obs_rate - pred_rate))
        return float(np.mean(ici_values)) if ici_values else np.nan
    except Exception:
        return np.nan


def main() -> None:
    np.random.seed(SEED)
    frame = load_frame()
    config = get_locked_gbsa_config()
    print(f"Locked GBSA config: {config}")

    # Training set
    train = frame[
        frame["cohort"].eq("TCGA-BRCA")
        & frame["os_days"].notna()
        & (frame["os_days"] > 0)
        & frame["os_event"].notna()
    ].copy()

    # External cohorts
    ext_dfs = {}
    for cohort in EXTERNAL_COHORTS:
        cdf = frame[frame["cohort"].eq(cohort)].copy()
        ext_dfs[cohort] = cdf

    # Load PAM50 scores
    pam50_scores = pd.read_csv(RESULTS / "06_external_model_predictions.csv")[
        ["sample_id", PAM50]
    ]

    rows = []
    for set_name, feature_list in FEATURE_SETS.items():
        print(f"  Training feature set: {set_name} ({len(feature_list)} features)")
        predictions = train_and_predict(feature_list, train, ext_dfs, config)

        for cohort in EXTERNAL_COHORTS:
            cdf = ext_dfs[cohort].copy()
            cdf["ablation_risk"] = predictions[cohort]
            cdf = cdf.merge(pam50_scores, on="sample_id", how="left")
            valid = cdf[
                cdf["ablation_risk"].notna()
                & cdf[PAM50].notna()
                & cdf["os_days"].notna()
                & (cdf["os_days"] > 0)
                & cdf["os_event"].notna()
            ].copy()
            n = len(valid)
            events = int(valid["os_event"].sum())

            if n < 10 or events < 2:
                rows.append({
                    "feature_set": set_name,
                    "n_features": len(feature_list),
                    "cohort": cohort,
                    "n": n,
                    "events": events,
                    "harrell_c": np.nan,
                    "harrell_c_ci_low": np.nan,
                    "harrell_c_ci_high": np.nan,
                    "uno_c": np.nan,
                    "brier_5y": np.nan,
                    "ici": np.nan,
                    "delta_vs_pam50": np.nan,
                    "delta_ci_low": np.nan,
                    "delta_ci_high": np.nan,
                    "delta_p": np.nan,
                })
                continue

            hc = harrell_c_index(valid["os_days"], valid["os_event"], valid["ablation_risk"])
            ci_low, ci_high = bootstrap_c_index_ci(
                valid["os_days"], valid["os_event"], valid["ablation_risk"],
                BOOTSTRAP_ITER, SEED,
            )
            uc = uno_c_index(valid["os_days"], valid["os_event"], valid["ablation_risk"])
            brier = brier_score_at(valid["os_days"], valid["os_event"], valid["ablation_risk"], years=(5,))
            ici = compute_ici(valid["os_days"], valid["os_event"], valid["ablation_risk"])

            delta = paired_bootstrap_delta_c_index(
                valid["os_days"], valid["os_event"],
                valid["ablation_risk"], valid[PAM50],
                n_bootstrap=BOOTSTRAP_ITER, seed=SEED,
            )

            rows.append({
                "feature_set": set_name,
                "n_features": len(feature_list),
                "cohort": cohort,
                "n": n,
                "events": events,
                "harrell_c": hc,
                "harrell_c_ci_low": ci_low,
                "harrell_c_ci_high": ci_high,
                "uno_c": uc,
                "brier_5y": brier.get("brier_5y", np.nan),
                "ici": ici,
                "delta_vs_pam50": delta["delta"],
                "delta_ci_low": delta["ci_low"],
                "delta_ci_high": delta["ci_high"],
                "delta_p": delta["p"],
            })

    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "Table_S5_feature_ablation.csv", index=False)

    # Summary table (Table 6): mean across external cohorts per feature set
    summary_rows = []
    for set_name, grp in detail.groupby("feature_set", sort=False):
        summary_rows.append({
            "feature_set": set_name,
            "n_features": int(grp["n_features"].iloc[0]),
            "mean_harrell_c": grp["harrell_c"].mean(),
            "mean_uno_c": grp["uno_c"].mean(),
            "mean_brier_5y": grp["brier_5y"].mean(),
            "mean_ici": grp["ici"].mean(),
            "mean_delta_vs_pam50": grp["delta_vs_pam50"].mean(),
            "n_cohorts": len(grp),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS / "Table_6_feature_ablation_summary.csv", index=False)

    # Print results
    print("\n" + "=" * 70)
    print("Table S5: Feature Ablation Detail (per cohort)")
    print("=" * 70)
    for _, row in detail.iterrows():
        print(
            f"  {row['feature_set']:30s} | {row['cohort']:10s} | "
            f"C={row['harrell_c']:.4f} | delta={row['delta_vs_pam50']:+.4f} "
            f"[{row['delta_ci_low']:.4f}, {row['delta_ci_high']:.4f}] p={row['delta_p']:.4f}"
        )

    print("\n" + "=" * 70)
    print("Table 6: Feature Ablation Summary (mean across external cohorts)")
    print("=" * 70)
    for _, row in summary.iterrows():
        print(
            f"  {row['feature_set']:30s} | nfeat={row['n_features']:2d} | "
            f"mean C={row['mean_harrell_c']:.4f} | "
            f"mean delta vs PAM50={row['mean_delta_vs_pam50']:+.4f}"
        )


if __name__ == "__main__":
    main()
