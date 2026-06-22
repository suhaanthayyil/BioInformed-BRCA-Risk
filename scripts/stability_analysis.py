#!/usr/bin/env python3
"""Stage 5: exploratory subsampling stability for the locked headline model."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit
from sksurv.ensemble import GradientBoostingSurvivalAnalysis
from sksurv.util import Surv

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    FIGURES,
    HEADLINE_MODEL,
    RESULTS,
    SEED,
    ensure_dirs,
    load_modeling_frame,
    log_phase3,
    spearman_pairwise,
)
from scripts.train_ml_zoo import cindex, make_preprocessor  # noqa: E402

N_SUBSAMPLES = 100


def load_config() -> dict:
    table = pd.read_csv(RESULTS / "Table_S2_ml_internal_cv.csv")
    row = table[table["model"].eq(HEADLINE_MODEL)].iloc[0]
    return json.loads(row["config"])


def permutation_importance(
    model: GradientBoostingSurvivalAnalysis,
    preprocessor,
    x_test: pd.DataFrame,
    time: np.ndarray,
    event: np.ndarray,
    features: list[str],
    rng: np.random.Generator,
) -> dict[str, float]:
    x_scaled = preprocessor.transform(x_test)
    baseline = cindex(time, event, model.predict(x_scaled))
    values = {}
    for feature in features:
        shuffled = x_test.copy()
        shuffled[feature] = rng.permutation(shuffled[feature].to_numpy())
        risk = model.predict(preprocessor.transform(shuffled))
        values[feature] = float(baseline - cindex(time, event, risk))
    return values


def run_analysis() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    frame = load_modeling_frame()
    train = frame[frame["cohort"].eq("TCGA-BRCA")].copy()
    train = train[train["os_days"].notna() & train["os_event"].notna() & (train["os_days"] > 0)]
    metadata = json.loads((REPO_ROOT / "data" / "processed" / "ml_model_zoo.metadata.json").read_text())
    features = metadata["features"]
    x = train[features].copy()
    time = pd.to_numeric(train["os_days"], errors="coerce").to_numpy(float)
    event = pd.to_numeric(train["os_event"], errors="coerce").fillna(0).to_numpy(int).astype(bool)
    config = load_config()
    splitter = StratifiedShuffleSplit(n_splits=N_SUBSAMPLES, train_size=0.8, random_state=SEED)
    rng = np.random.default_rng(SEED)
    importance_rows = []
    rank_rows = []
    for run, (sub_idx, oob_idx) in enumerate(splitter.split(x, event.astype(int)), start=1):
        x_sub = x.iloc[sub_idx]
        x_oob = x.iloc[oob_idx]
        time_sub, event_sub = time[sub_idx], event[sub_idx]
        time_oob, event_oob = time[oob_idx], event[oob_idx]
        if event_oob.sum() < 5:
            log_phase3(f"Stage 5 stability run {run} skipped because out-of-bag events={int(event_oob.sum())}.")
            continue
        preprocessor = make_preprocessor()
        x_sub_scaled = preprocessor.fit_transform(x_sub)
        model = GradientBoostingSurvivalAnalysis(random_state=SEED + run, **config)
        model.fit(x_sub_scaled, Surv.from_arrays(event_sub, time_sub))
        values = permutation_importance(model, preprocessor, x_oob, time_oob, event_oob, features, rng)
        ranked = pd.Series(values).rank(ascending=False, method="average")
        for feature in features:
            importance_rows.append({"run": run, "feature": feature, "importance": values[feature], "rank": ranked[feature]})
        rank_rows.append([float(ranked[feature]) for feature in features])
    importance = pd.DataFrame(importance_rows)
    rank_matrix = np.asarray(rank_rows, dtype=float)
    rho = spearman_pairwise(rank_matrix)
    summary = (
        importance.groupby("feature")
        .agg(
            mean_importance=("importance", "mean"),
            sd_importance=("importance", "std"),
            mean_rank=("rank", "mean"),
            sd_rank=("rank", "std"),
            top3_rate=("rank", lambda s: float((s <= 3).mean())),
        )
        .reset_index()
        .sort_values("mean_rank")
    )
    summary["mean_pairwise_spearman"] = float(np.nanmean(rho)) if len(rho) else np.nan
    summary["spearman_ci_low"] = float(np.nanquantile(rho, 0.025)) if len(rho) else np.nan
    summary["spearman_ci_high"] = float(np.nanquantile(rho, 0.975)) if len(rho) else np.nan
    summary["n_subsamples"] = rank_matrix.shape[0]
    summary["n_pairwise_comparisons"] = len(rho)
    return summary, importance, rank_matrix


def plot_heatmap(rank_matrix: np.ndarray, features: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    im = ax.imshow(rank_matrix.T, aspect="auto", cmap="viridis_r", interpolation="nearest")
    ax.set_yticks(np.arange(len(features)))
    ax.set_yticklabels(features, fontsize=7)
    ax.set_xlabel("Subsample run")
    ax.set_title("Headline model permutation-importance rank stability")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Rank, lower is more important")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_stability_heatmap.pdf")
    fig.savefig(FIGURES / "fig_stability_heatmap.tiff", dpi=300)
    plt.close(fig)


def update_story(summary: pd.DataFrame) -> None:
    story = REPO_ROOT / "docs" / "STORY.md"
    text = story.read_text()
    mean_rho = float(summary["mean_pairwise_spearman"].iloc[0]) if not summary.empty else np.nan
    low = float(summary["spearman_ci_low"].iloc[0]) if not summary.empty else np.nan
    high = float(summary["spearman_ci_high"].iloc[0]) if not summary.empty else np.nan
    top = ", ".join(summary.head(3)["feature"].tolist()) if not summary.empty else "not evaluable"
    verdict = f"mean pairwise Spearman rho={mean_rho:.3f} (95% empirical interval {low:.3f} to {high:.3f}); top-ranked features: {top}."
    text = text.replace("- Stability of pathway attributions: pending Stage 5.", f"- Stability of pathway attributions: {verdict}")
    story.write_text(text)
    if np.isfinite(mean_rho) and mean_rho < 0.6:
        log_phase3(f"Stage 5 stability flagged low mean Spearman rho={mean_rho:.3f}.")


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 5 stability analysis started.")
    summary, importance, rank_matrix = run_analysis()
    # R1.6/R3.4: these rankings are descriptive signals from this modelling setup,
    # not mechanistic biological conclusions. Mark it on the artifact itself.
    summary["interpretation"] = "descriptive_signal_not_mechanistic"
    importance["interpretation"] = "descriptive_signal_not_mechanistic"
    summary.to_csv(RESULTS / "Table_S13_stability.csv", index=False)
    importance.to_csv(RESULTS / "Table_S13_stability_ranks.csv", index=False)
    features = summary["feature"].tolist() if not summary.empty else []
    plot_heatmap(rank_matrix, features)
    update_story(summary)
    log_phase3("Stage 5 stability analysis completed.")


if __name__ == "__main__":
    main()
