#!/usr/bin/env python3
"""Stage 4: exploratory decision curve analysis and NRI/IDI."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.calibration_analysis import fit_probability_mapping, predict_probability  # noqa: E402
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
MODEL_LABELS = {
    HEADLINE_MODEL: "Pathway ML",
    COX_MODEL: "Cox baseline",
    PAM50_MODEL: "PAM50-ROR",
}
THRESHOLDS = np.round(np.arange(0.05, 0.5001, 0.025), 3)
HORIZON_DAYS = 5 * 365


def net_benefit(y: np.ndarray, p: np.ndarray, threshold: float) -> float:
    valid = np.isfinite(y) & np.isfinite(p)
    y = y[valid].astype(int)
    p = p[valid].astype(float)
    if len(y) == 0:
        return np.nan
    positive = p >= threshold
    tp = int(((y == 1) & positive).sum())
    fp = int(((y == 0) & positive).sum())
    odds = threshold / (1 - threshold)
    return float(tp / len(y) - fp / len(y) * odds)


def treat_all_net_benefit(y: np.ndarray, threshold: float) -> float:
    valid = np.isfinite(y)
    y = y[valid].astype(int)
    if len(y) == 0:
        return np.nan
    prevalence = y.mean()
    return float(prevalence - (1 - prevalence) * threshold / (1 - threshold))


def nri_idi_stats(y: np.ndarray, p_model: np.ndarray, p_ref: np.ndarray, threshold: float) -> dict[str, float]:
    valid = np.isfinite(y) & np.isfinite(p_model) & np.isfinite(p_ref)
    y = y[valid].astype(int)
    p_model = p_model[valid].astype(float)
    p_ref = p_ref[valid].astype(float)
    if len(y) < 30 or y.sum() == 0 or (y == 0).sum() == 0:
        return {"nri": np.nan, "idi": np.nan}
    cat_model = (p_model >= threshold).astype(int)
    cat_ref = (p_ref >= threshold).astype(int)
    event = y == 1
    nonevent = y == 0
    event_up = ((cat_model > cat_ref) & event).sum() / event.sum()
    event_down = ((cat_model < cat_ref) & event).sum() / event.sum()
    nonevent_down = ((cat_model < cat_ref) & nonevent).sum() / nonevent.sum()
    nonevent_up = ((cat_model > cat_ref) & nonevent).sum() / nonevent.sum()
    nri = float((event_up - event_down) + (nonevent_down - nonevent_up))
    slope_model = p_model[event].mean() - p_model[nonevent].mean()
    slope_ref = p_ref[event].mean() - p_ref[nonevent].mean()
    idi = float(slope_model - slope_ref)
    return {"nri": nri, "idi": idi}


def bootstrap_nri_idi(y: np.ndarray, p_model: np.ndarray, p_ref: np.ndarray, threshold: float) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(SEED)
    values = {"nri": [], "idi": []}
    n = len(y)
    if n < 30:
        return {key: (np.nan, np.nan) for key in values}
    for _ in range(BOOTSTRAP_ITER):
        idx = rng.integers(0, n, n)
        stat = nri_idi_stats(y[idx], p_model[idx], p_ref[idx], threshold)
        for key in values:
            values[key].append(stat[key])
    return {key: ci(vals) for key, vals in values.items()}


def probability_table() -> pd.DataFrame:
    pred = binary_horizon_outcome(load_predictions(), HORIZON_DAYS)
    train = pred[pred["cohort"].eq("TCGA-BRCA")].copy()
    mappings = {model: fit_probability_mapping(train, model) for model in MODELS}
    parts = []
    for cohort in ALL_COHORTS:
        cohort_df = pred[pred["cohort"].eq(cohort)].copy()
        keep = ["sample_id", "cohort", "os_days", "os_event", "horizon_event"]
        out = cohort_df[keep].copy()
        for model in MODELS:
            out[f"{model}_prob_5y"] = predict_probability(mappings[model], cohort_df, model)
        parts.append(out)
    return pd.concat(parts, ignore_index=True)


def run_dca(probs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort in ALL_COHORTS:
        sub = probs[probs["cohort"].eq(cohort)].copy()
        y = sub["horizon_event"].to_numpy(int)
        for threshold in THRESHOLDS:
            rows.append(
                {
                    "cohort": cohort,
                    "threshold": threshold,
                    "strategy": "Treat all",
                    "net_benefit": treat_all_net_benefit(y, threshold),
                    "n": len(sub),
                    "events": int(y.sum()),
                }
            )
            rows.append(
                {
                    "cohort": cohort,
                    "threshold": threshold,
                    "strategy": "Treat none",
                    "net_benefit": 0.0,
                    "n": len(sub),
                    "events": int(y.sum()),
                }
            )
            for model in MODELS:
                rows.append(
                    {
                        "cohort": cohort,
                        "threshold": threshold,
                        "strategy": MODEL_LABELS[model],
                        "net_benefit": net_benefit(y, sub[f"{model}_prob_5y"].to_numpy(float), threshold),
                        "n": len(sub),
                        "events": int(y.sum()),
                    }
                )
    return pd.DataFrame(rows)


def run_nri_idi(probs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort in ALL_COHORTS:
        sub = probs[probs["cohort"].eq(cohort)].copy()
        y = sub["horizon_event"].to_numpy(int)
        p_model = sub[f"{HEADLINE_MODEL}_prob_5y"].to_numpy(float)
        p_ref = sub[f"{PAM50_MODEL}_prob_5y"].to_numpy(float)
        for threshold in (0.10, 0.20):
            stat = nri_idi_stats(y, p_model, p_ref, threshold)
            boot = bootstrap_nri_idi(y, p_model, p_ref, threshold)
            rows.append(
                {
                    "cohort": cohort,
                    "threshold": threshold,
                    "n": len(sub),
                    "events": int(y.sum()),
                    "nri": stat["nri"],
                    "nri_ci_low": boot["nri"][0],
                    "nri_ci_high": boot["nri"][1],
                    "idi": stat["idi"],
                    "idi_ci_low": boot["idi"][0],
                    "idi_ci_high": boot["idi"][1],
                    "comparison": "Gradient_Boosted_Survival minus PAM50_ROR_official",
                    "analysis_type": "exploratory_secondary",
                }
            )
    return pd.DataFrame(rows)


def plot_dca(dca: pd.DataFrame) -> None:
    colors = {
        "Treat all": "#999999",
        "Treat none": "#000000",
        "Pathway ML": "#31688e",
        "Cox baseline": "#35b779",
        "PAM50-ROR": "#440154",
    }
    linestyles = {"Treat all": "--", "Treat none": ":", "Pathway ML": "-", "Cox baseline": "-", "PAM50-ROR": "-"}
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
    for ax, cohort in zip(axes.ravel(), ALL_COHORTS, strict=False):
        sub = dca[dca["cohort"].eq(cohort)]
        for strategy, group in sub.groupby("strategy"):
            ax.plot(
                group["threshold"],
                group["net_benefit"],
                label=strategy,
                color=colors.get(strategy, "#333333"),
                ls=linestyles.get(strategy, "-"),
                lw=1.4,
            )
        ax.set_title(cohort, fontsize=9)
        ax.set_xlabel("5-year mortality threshold")
        ax.set_ylabel("Net benefit")
        ax.grid(alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(FIGURES / "fig_dca_per_cohort.pdf")
    fig.savefig(FIGURES / "fig_dca_per_cohort.tiff", dpi=300)
    plt.close(fig)


def update_story(dca: pd.DataFrame) -> None:
    story = REPO_ROOT / "docs" / "STORY.md"
    text = story.read_text()
    wide = dca.pivot_table(index=["cohort", "threshold"], columns="strategy", values="net_benefit")
    if {"Pathway ML", "PAM50-ROR"}.issubset(wide.columns):
        delta = (wide["Pathway ML"] - wide["PAM50-ROR"]).dropna()
        mean_delta = float(delta.mean()) if len(delta) else np.nan
        verdict = (
            f"mean 5-year net benefit delta versus PAM50-ROR was {mean_delta:+.4f} across cohorts and thresholds."
        )
    else:
        verdict = "not evaluable for all strategies."
    text = text.replace("- Net benefit (DCA): pending Stage 4.", f"- Net benefit (DCA): {verdict}")
    story.write_text(text)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 4 DCA and NRI/IDI analysis started.")
    probs = probability_table()
    dca = run_dca(probs)
    nri_idi = run_nri_idi(probs)
    dca["analysis_type"] = "exploratory_secondary"  # secondary metric; not the primary endpoint
    dca.to_csv(RESULTS / "Table_S11_dca.csv", index=False)
    nri_idi.to_csv(RESULTS / "Table_S12_nri_idi.csv", index=False)
    plot_dca(dca)
    update_story(dca)
    log_phase3("Stage 4 DCA and NRI/IDI analysis completed.")


if __name__ == "__main__":
    main()
