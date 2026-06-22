#!/usr/bin/env python3
"""Figure 6: permutation-importance rank stability across 100 stratified subsamples.

Panel A: per-subsample rank of each feature (1 = most important).
Panel B: mean rank +/- SD with the fraction of subsamples in which the feature
ranked in the top 3.

Subsamples are stratified 80% draws of TCGA-BRCA without replacement (not
bootstrap resamples), so the x-axis is labelled accordingly.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES = REPO_ROOT / "figures"
RESULTS = REPO_ROOT / "results"

DISPLAY = {
    "Pathway_Hormone": "Pathway-Hormone",
    "age_at_dx": "age at diagnosis",
    "Pathway_DNA_Repair": "Pathway-DNA-Repair",
    "stage_ordinal": "stage_ordinal",
    "Pathway_Proliferation": "Pathway-Proliferation",
    "Pathway_Apoptosis_Stress": "Pathway-Apoptosis-Stress",
    "Pathway_Metabolism": "Pathway-Metabolism",
    "Pathway_Immune": "Pathway-Immune",
    "Pathway_Stromal_EMT": "Pathway-Stromal-EMT",
}


def main() -> None:
    df = pd.read_csv(RESULTS / "Table_S13_stability_ranks.csv")
    piv = df.pivot(index="run", columns="feature", values="rank").sort_index()
    order = piv.mean().sort_values().index.tolist()  # most important (lowest mean rank) first
    piv = piv[order]
    labels = [DISPLAY.get(f, f) for f in order]
    mat = piv.to_numpy().T  # features x runs
    n_runs = piv.shape[0]
    mean_rank = piv.mean().to_numpy()
    sd_rank = piv.std().to_numpy()
    top3 = (piv <= 3).mean().to_numpy() * 100.0

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.2), gridspec_kw={"width_ratios": [1.55, 1.0]})

    # Panel A: rank heatmap
    im = axA.imshow(mat, aspect="auto", cmap="viridis_r", vmin=1, vmax=mat.shape[0], interpolation="nearest")
    axA.set_yticks(range(len(labels)))
    axA.set_yticklabels(labels, fontsize=9)
    xt = [t - 1 for t in (1, 20, 40, 60, 80, 100) if t <= n_runs]
    axA.set_xticks(xt)
    axA.set_xticklabels([t + 1 for t in xt])
    axA.set_xlabel("Stratified 80% subsample", fontsize=11)
    axA.set_title("A", loc="left", fontsize=14, fontweight="bold")
    cbar = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.02)
    cbar.set_label("Rank (1 = most important)", fontsize=9)

    # Panel B: mean rank +/- SD with top-3 frequency
    y = np.arange(len(labels))
    axB.barh(y, mean_rank, color="#5b8db8", edgecolor="none")
    axB.errorbar(mean_rank, y, xerr=sd_rank, fmt="none", ecolor="#1a1a1a", capsize=3, lw=1.0)
    for yi, m, s, t in zip(y, mean_rank, sd_rank, top3, strict=False):
        axB.text(min(m + s + 0.25, mat.shape[0] + 0.6), yi, f"top-3: {t:.0f}%",
                 va="center", ha="left", fontsize=8.5)
    axB.set_yticks(y)
    axB.set_yticklabels(labels, fontsize=9)
    axB.invert_yaxis()
    axB.set_xlim(0, mat.shape[0] + 1)
    axB.set_xlabel("Mean rank ± SD", fontsize=11)
    axB.set_title("B", loc="left", fontsize=14, fontweight="bold")
    axB.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIGURES / "fig_figure6_stability.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure6_stability.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure6_stability.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure 6 written: runs =", n_runs, "features =", len(labels))


if __name__ == "__main__":
    main()
