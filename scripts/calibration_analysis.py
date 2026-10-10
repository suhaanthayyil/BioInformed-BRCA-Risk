#!/usr/bin/env python3
"""Stage 3: per-cohort calibration analysis."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import patsy
import pandas as pd
from lifelines import CoxPHFitter
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
    ci,
    fit_risk_mapping,
    horizon_reached,
    ipcw_horizon,
    km_risk,
    predicted_risk,
    valid_survival_subset,
    ensure_dirs,
    load_predictions,
    log_phase3,
)

MODELS = [HEADLINE_MODEL, COX_MODEL, PAM50_MODEL]
HORIZONS = {5: 5 * 365, 10: 10 * 365}


def logit(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, 1e-4, 1 - 1e-4)
    return np.log(clipped / (1 - clipped))


def smoothed_ici(df: pd.DataFrame, p: np.ndarray, days: float) -> tuple[float, float, float]:
    x = np.log(-np.log(np.clip(1 - p, 1e-6, 1 - 1e-6)))
    basis = np.asarray(patsy.dmatrix("cr(x, df=3) - 1", {"x": x}))
    data = pd.DataFrame(basis, columns=[f"s{i}" for i in range(basis.shape[1])])
    data["T"] = df["os_days"].to_numpy(float)
    data["E"] = df["os_event"].to_numpy(int)
    try:
        cox = CoxPHFitter(penalizer=1e-4).fit(data, "T", "E")
    except Exception:
        return np.nan, np.nan, np.nan
    observed = 1 - cox.predict_survival_function(data.drop(columns=["T", "E"]), times=[days]).to_numpy().ravel()
    gap = np.abs(observed - p)
    return float(gap.mean()), float(np.median(gap)), float(np.quantile(gap, 0.9))


def calibration_stats(df: pd.DataFrame, p: np.ndarray, days: float) -> dict[str, float]:
    keep = np.isfinite(p)
    df, p = df[keep], p[keep]
    y, w = ipcw_horizon(df, days)
    use = w > 0
    empty = {"slope": np.nan, "intercept": np.nan, "calibration_in_large": np.nan, "ici": np.nan, "e50": np.nan, "e90": np.nan, "decile_error": np.nan}
    if use.sum() < 30 or len(np.unique(y[use])) < 2:
        return empty
    cal = LogisticRegression(random_state=SEED)
    cal.fit(logit(p[use]).reshape(-1, 1), y[use], sample_weight=w[use])
    bins = pd.qcut(p, q=min(10, len(np.unique(p))), labels=False, duplicates="drop")
    gaps, sizes = [], []
    for b in np.unique(bins):
        part = df[bins == b]
        gaps.append(abs(km_risk(part, days) - p[bins == b].mean()))
        sizes.append(len(part))
    ici, e50, e90 = smoothed_ici(df, p, days)
    return {
        "slope": float(cal.coef_[0][0]),
        "intercept": float(cal.intercept_[0]),
        "calibration_in_large": float(km_risk(df, days) - p.mean()),
        "ici": ici,
        "e50": e50,
        "e90": e90,
        "decile_error": float(np.average(gaps, weights=sizes)),
    }


def bootstrap_stats(df: pd.DataFrame, p: np.ndarray, days: float) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(SEED)
    values = {"slope": [], "intercept": [], "calibration_in_large": [], "ici": [], "e50": [], "e90": [], "decile_error": []}
    n = len(df)
    if n < 30:
        return {key: (np.nan, np.nan) for key in values}
    for _ in range(BOOTSTRAP_ITER):
        idx = rng.integers(0, n, n)
        stat = calibration_stats(df.iloc[idx], p[idx], days)
        for key in values:
            values[key].append(stat[key])
    return {key: ci(vals) for key, vals in values.items()}


def calibration_curve_points(df: pd.DataFrame, p: np.ndarray, days: float) -> pd.DataFrame:
    keep = np.isfinite(p)
    df, p = df[keep], p[keep]
    if len(df) < 30:
        return pd.DataFrame(columns=["pred", "obs"])
    bins = pd.qcut(p, q=min(10, len(np.unique(p))), labels=False, duplicates="drop")
    return pd.DataFrame(
        [{"pred": p[bins == b].mean(), "obs": km_risk(df[bins == b], days)} for b in np.unique(bins)]
    )


def run_analysis() -> tuple[pd.DataFrame, dict[tuple[str, int, str], pd.DataFrame]]:
    pred = valid_survival_subset(load_predictions(), [])
    rows = []
    curves = {}
    train = pred[pred["cohort"].eq("TCGA-BRCA")]
    mappings = {model: fit_risk_mapping(train, model) for model in MODELS}
    for years, days in HORIZONS.items():
        for cohort in ALL_COHORTS:
            cohort_df = pred[pred["cohort"].eq(cohort)].copy()
            if not horizon_reached(cohort_df, days):
                log_phase3(f"Calibration skipped {cohort} {years}y: follow-up does not reach the horizon.")
                continue
            for model in MODELS:
                sub = cohort_df[cohort_df[model].notna()]
                probs = predicted_risk(mappings[model], sub, model, days)
                stat = calibration_stats(sub, probs, days)
                boot = bootstrap_stats(sub, probs, days)
                curves[(cohort, years, model)] = calibration_curve_points(sub, probs, days)
                rows.append(
                    {
                        "cohort": cohort,
                        "model": model,
                        "time_horizon": f"{years}y",
                        "n": len(sub),
                        "events_by_horizon": int(((sub["os_event"] == 1) & (sub["os_days"] <= days)).sum()),
                        "observed_km_risk": km_risk(sub, days),
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
                        "e50": stat["e50"],
                        "e90": stat["e90"],
                        "decile_error": stat["decile_error"],
                        "decile_error_ci_low": boot["decile_error"][0],
                        "decile_error_ci_high": boot["decile_error"][1],
                        "ici_method": "Austin 2020: Cox with restricted cubic spline of cloglog(predicted risk)",
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
            if not any((cohort, years, m) in curves for m in MODELS):
                ax.text(0.5, 0.5, "follow-up does not\nreach horizon", ha="center", va="center", fontsize=8)
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
