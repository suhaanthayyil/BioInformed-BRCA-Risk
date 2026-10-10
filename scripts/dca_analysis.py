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
    ipcw_horizon,
    km_risk,
    predicted_risk,
    valid_survival_subset,
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


def net_benefit(df: pd.DataFrame, p: np.ndarray, threshold: float, days: float) -> float:
    keep = np.isfinite(p)
    df, p = df[keep], p[keep]
    if len(df) == 0:
        return np.nan
    positive = p >= threshold
    share = positive.mean()
    if share == 0:
        return 0.0
    risk = km_risk(df[positive], days)
    odds = threshold / (1 - threshold)
    return float(share * risk - share * (1 - risk) * odds)


def treat_all_net_benefit(df: pd.DataFrame, threshold: float, days: float) -> float:
    risk = km_risk(df, days)
    return float(risk - (1 - risk) * threshold / (1 - threshold))


def nri_idi_stats(y: np.ndarray, w: np.ndarray, p_model: np.ndarray, p_ref: np.ndarray, threshold: float) -> dict[str, float]:
    valid = np.isfinite(p_model) & np.isfinite(p_ref) & (w > 0)
    y, w, p_model, p_ref = y[valid], w[valid], p_model[valid], p_ref[valid]
    ev, ne = w * (y == 1), w * (y == 0)
    if len(y) < 30 or ev.sum() == 0 or ne.sum() == 0:
        return {"nri": np.nan, "idi": np.nan}
    up = (p_model >= threshold) & (p_ref < threshold)
    down = (p_model < threshold) & (p_ref >= threshold)
    nri = (ev[up].sum() - ev[down].sum()) / ev.sum() + (ne[down].sum() - ne[up].sum()) / ne.sum()
    slope_model = np.average(p_model, weights=ev) - np.average(p_model, weights=ne)
    slope_ref = np.average(p_ref, weights=ev) - np.average(p_ref, weights=ne)
    return {"nri": float(nri), "idi": float(slope_model - slope_ref)}


def bootstrap_nri_idi(y, w, p_model, p_ref, threshold) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(SEED)
    values = {"nri": [], "idi": []}
    n = len(y)
    if n < 30:
        return {key: (np.nan, np.nan) for key in values}
    for _ in range(BOOTSTRAP_ITER):
        idx = rng.integers(0, n, n)
        stat = nri_idi_stats(y[idx], w[idx], p_model[idx], p_ref[idx], threshold)
        for key in values:
            values[key].append(stat[key])
    return {key: ci(vals) for key, vals in values.items()}


def probability_table() -> pd.DataFrame:
    pred = valid_survival_subset(load_predictions(), MODELS)
    train = pred[pred["cohort"].eq("TCGA-BRCA")]
    mappings = {model: fit_risk_mapping(train, model) for model in MODELS}
    for model in MODELS:
        pred[f"{model}_prob_5y"] = predicted_risk(mappings[model], pred, model, HORIZON_DAYS)
    return pred


def run_dca(probs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort in ALL_COHORTS:
        sub = probs[probs["cohort"].eq(cohort)].copy()
        n = len(sub)
        events = int(((sub["os_event"] == 1) & (sub["os_days"] <= HORIZON_DAYS)).sum())
        for threshold in THRESHOLDS:
            base = {"cohort": cohort, "threshold": threshold, "n": n, "events_by_horizon": events}
            rows.append({**base, "strategy": "Treat all", "net_benefit": treat_all_net_benefit(sub, threshold, HORIZON_DAYS)})
            rows.append({**base, "strategy": "Treat none", "net_benefit": 0.0})
            for model in MODELS:
                nb = net_benefit(sub, sub[f"{model}_prob_5y"].to_numpy(float), threshold, HORIZON_DAYS)
                rows.append({**base, "strategy": MODEL_LABELS[model], "net_benefit": nb})
    return pd.DataFrame(rows)


def run_nri_idi(probs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort in ALL_COHORTS:
        sub = probs[probs["cohort"].eq(cohort)].copy()
        y, w = ipcw_horizon(sub, HORIZON_DAYS)
        p_model = sub[f"{HEADLINE_MODEL}_prob_5y"].to_numpy(float)
        p_ref = sub[f"{PAM50_MODEL}_prob_5y"].to_numpy(float)
        for threshold in (0.10, 0.20):
            stat = nri_idi_stats(y, w, p_model, p_ref, threshold)
            boot = bootstrap_nri_idi(y, w, p_model, p_ref, threshold)
            rows.append(
                {
                    "cohort": cohort,
                    "threshold": threshold,
                    "n": len(sub),
                    "events_by_horizon": int(y.sum()),
                    "nri": stat["nri"],
                    "nri_ci_low": boot["nri"][0],
                    "nri_ci_high": boot["nri"][1],
                    "idi": stat["idi"],
                    "idi_ci_low": boot["idi"][0],
                    "idi_ci_high": boot["idi"][1],
                    "comparison": f"{HEADLINE_MODEL} minus {PAM50_MODEL}",
                    "method": "IPCW-weighted",
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
    if not story.exists():
        return  # STORY.md was removed from the repo; nothing to update
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
