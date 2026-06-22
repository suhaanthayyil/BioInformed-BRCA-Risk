#!/usr/bin/env python3
"""Stage 3: per-cohort calibration analysis."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    ALL_COHORTS,
    BOOTSTRAP_ITER,
    COX_MODEL,
    FIGURES,
    HEADLINE_MODEL,
    PAM50_MODEL,
    RESULTS,
    SEED,
    binary_horizon_outcome,
    ci,
    ensure_dirs,
    load_predictions,
    log_phase3,
)

MODELS = [HEADLINE_MODEL, COX_MODEL, PAM50_MODEL]
HORIZONS = {5: 5 * 365, 10: 10 * 365}


def fit_probability_mapping(train: pd.DataFrame, risk_col: str) -> LogisticRegression | None:
    df = train.dropna(subset=[risk_col, "horizon_event"])
    if len(df) < 30 or df["horizon_event"].nunique() < 2:
        return None
    model = LogisticRegression(random_state=SEED)
    model.fit(df[[risk_col]], df["horizon_event"])
    return model


def predict_probability(model: LogisticRegression | None, df: pd.DataFrame, risk_col: str) -> np.ndarray:
    if model is None:
        return np.full(len(df), np.nan)
    return model.predict_proba(df[[risk_col]])[:, 1]


def logit(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, 1e-4, 1 - 1e-4)
    return np.log(clipped / (1 - clipped))


def calibration_stats(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    valid = np.isfinite(y) & np.isfinite(p)
    y = y[valid].astype(int)
    p = p[valid].astype(float)
    if len(y) < 30 or len(np.unique(y)) < 2:
        return {"slope": np.nan, "intercept": np.nan, "calibration_in_large": np.nan, "ici": np.nan}
    x = logit(p).reshape(-1, 1)
    cal = LogisticRegression(random_state=SEED)
    cal.fit(x, y)
    slope = float(cal.coef_[0][0])
    intercept = float(cal.intercept_[0])
    calibration_in_large = float(y.mean() - p.mean())
    tmp = pd.DataFrame({"y": y, "p": p})
    tmp["bin"] = pd.qcut(tmp["p"], q=min(10, tmp["p"].nunique()), duplicates="drop")
    grouped = tmp.groupby("bin", observed=False).agg(obs=("y", "mean"), pred=("p", "mean"), n=("y", "size"))
    ici = float(np.average(np.abs(grouped["obs"] - grouped["pred"]), weights=grouped["n"]))
    return {"slope": slope, "intercept": intercept, "calibration_in_large": calibration_in_large, "ici": ici}


def bootstrap_stats(y: np.ndarray, p: np.ndarray) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(SEED)
    values = {"slope": [], "intercept": [], "calibration_in_large": [], "ici": []}
    n = len(y)
    if n < 30:
        return {key: (np.nan, np.nan) for key in values}
    for _ in range(BOOTSTRAP_ITER):
        idx = rng.integers(0, n, n)
        stat = calibration_stats(y[idx], p[idx])
        for key in values:
            values[key].append(stat[key])
    return {key: ci(vals) for key, vals in values.items()}


def calibration_curve_points(y: np.ndarray, p: np.ndarray) -> pd.DataFrame:
    tmp = pd.DataFrame({"y": y, "p": p}).dropna()
    if len(tmp) < 30:
        return pd.DataFrame(columns=["pred", "obs"])
    tmp["bin"] = pd.qcut(tmp["p"], q=min(10, tmp["p"].nunique()), duplicates="drop")
    return tmp.groupby("bin", observed=False).agg(pred=("p", "mean"), obs=("y", "mean")).reset_index(drop=True)


def run_analysis() -> tuple[pd.DataFrame, dict[tuple[str, int, str], pd.DataFrame]]:
    pred = load_predictions()
    rows = []
    curves = {}
    for years, days in HORIZONS.items():
        horizon = binary_horizon_outcome(pred, days)
        train = horizon[horizon["cohort"].eq("TCGA-BRCA")].copy()
        mappings = {model: fit_probability_mapping(train, model) for model in MODELS}
        for cohort in ALL_COHORTS:
            cohort_df = horizon[horizon["cohort"].eq(cohort)].copy()
            for model in MODELS:
                probs = predict_probability(mappings[model], cohort_df, model)
                y = cohort_df["horizon_event"].to_numpy(int)
                stat = calibration_stats(y, probs)
                boot = bootstrap_stats(y, probs)
                curves[(cohort, years, model)] = calibration_curve_points(y, probs)
                rows.append(
                    {
                        "cohort": cohort,
                        "model": model,
                        "time_horizon": f"{years}y",
                        "n": len(cohort_df),
                        "events": int(cohort_df["horizon_event"].sum()) if len(cohort_df) else 0,
                        "slope": stat["slope"],
                        "slope_ci_low": boot["slope"][0],
                        "slope_ci_high": boot["slope"][1],
                        "intercept": stat["intercept"],
                        "intercept_ci_low": boot["intercept"][0],
                        "intercept_ci_high": boot["intercept"][1],
                        "calibration_in_large": stat["calibration_in_large"],
                        "calibration_in_large_ci_low": boot["calibration_in_large"][0],
                        "calibration_in_large_ci_high": boot["calibration_in_large"][1],
                        "ici": stat["ici"],
                        "ici_ci_low": boot["ici"][0],
                        "ici_ci_high": boot["ici"][1],
                    }
                )
    return pd.DataFrame(rows), curves


def plot_calibration(curves: dict[tuple[str, int, str], pd.DataFrame]) -> None:
    colors = {HEADLINE_MODEL: "#31688e", COX_MODEL: "#35b779", PAM50_MODEL: "#440154"}
    fig, axes = plt.subplots(2, 4, figsize=(12, 6), sharex=True, sharey=True)
    for row, years in enumerate(HORIZONS):
        for col, cohort in enumerate(ALL_COHORTS):
            ax = axes[row, col]
            ax.plot([0, 1], [0, 1], color="black", lw=1, ls=":")
            for model in MODELS:
                pts = curves.get((cohort, years, model), pd.DataFrame())
                if not pts.empty:
                    ax.plot(pts["pred"], pts["obs"], marker="o", ms=3, lw=1, label=model, color=colors[model])
            ax.set_title(f"{cohort}, {years}y", fontsize=9)
            ax.grid(alpha=0.2)
            if row == 1:
                ax.set_xlabel("Predicted mortality")
            if col == 0:
                ax.set_ylabel("Observed mortality")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(FIGURES / "fig_calibration_plots.pdf")
    fig.savefig(FIGURES / "fig_calibration_plots.tiff", dpi=300)
    plt.close(fig)


def update_story(table: pd.DataFrame) -> None:
    story = REPO_ROOT / "docs" / "STORY.md"
    if not story.exists():
        return  # STORY.md was removed from the repo; nothing to update
    text = story.read_text()
    mean_ici = table.groupby("model")["ici"].mean(numeric_only=True)
    if HEADLINE_MODEL in mean_ici and PAM50_MODEL in mean_ici:
        if mean_ici[HEADLINE_MODEL] < mean_ici[PAM50_MODEL]:
            verdict = f"headline ML showed lower mean ICI than PAM50-ROR ({mean_ici[HEADLINE_MODEL]:.3f} vs {mean_ici[PAM50_MODEL]:.3f})."
        else:
            verdict = f"headline ML did not improve mean ICI versus PAM50-ROR ({mean_ici[HEADLINE_MODEL]:.3f} vs {mean_ici[PAM50_MODEL]:.3f})."
    else:
        verdict = "not evaluable in all cohorts."
    text = text.replace("- Calibration: pending Stage 3.", f"- Calibration: {verdict}")
    story.write_text(text)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 3 calibration analysis started.")
    table, curves = run_analysis()
    table["analysis_type"] = "exploratory_secondary"  # secondary metric; not the primary endpoint
    table.to_csv(RESULTS / "Table_S10_calibration.csv", index=False)
    plot_calibration(curves)
    update_story(table)
    log_phase3("Stage 3 calibration analysis completed.")


if __name__ == "__main__":
    main()
