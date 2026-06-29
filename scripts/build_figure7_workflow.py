#!/usr/bin/env python3
"""Figure 7 / graphical abstract: study design and data-processing workflow.

A single schematic summarising cohort harmonization, pathway scoring, locked-model
training, external validation, and the pre-registered primary result.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

FIGURES = Path(__file__).resolve().parents[1] / "figures"

BLUE_F, BLUE_E = "#dbe9f6", "#2c7fb8"
GREY_F, GREY_E = "#f2f2f2", "#555555"
GREEN_F, GREEN_E = "#e5f5e0", "#31a354"
GOLD_F, GOLD_E = "#fde9c8", "#e6a23c"
RED = "#c0392b"


def box(ax, x, y, w, h, text, fc, ec, fs=9, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.004,rounding_size=0.015",
                                fc=fc, ec=ec, lw=1.2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal")


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=13, lw=1.1, color="#444444"))


def main() -> None:
    fig, ax = plt.subplots(figsize=(9.5, 8.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.965, "Four public breast cancer cohorts  (harmonized N = 4,532)",
            ha="center", fontsize=10.5, fontweight="bold")
    cohorts = [("TCGA-BRCA\nn = 213\n(training)", 0.16),
               ("GSE96058 / SCAN-B\nn = 1,483", 0.385),
               ("METABRIC\nn = 1,980 eval.", 0.62),
               ("GSE20685\nn = 327", 0.85)]
    for txt, x in cohorts:
        box(ax, x, 0.875, 0.205, 0.078, txt, GREY_F, GREY_E, fs=8)
        arrow(ax, x, 0.836, x, 0.781)

    box(ax, 0.5, 0.742, 0.84, 0.072,
        "Per-cohort gene-level z-scoring (no scaling transfer)  →  7 MSigDB pathway activity scores",
        BLUE_F, BLUE_E, fs=9)
    arrow(ax, 0.5, 0.706, 0.5, 0.671)
    box(ax, 0.5, 0.635, 0.84, 0.062,
        "9 features = 7 pathway scores + 2 clinical covariates (age, ordinal stage)",
        BLUE_F, BLUE_E, fs=9)
    arrow(ax, 0.5, 0.604, 0.30, 0.548)

    box(ax, 0.21, 0.505, 0.36, 0.085,
        "Train 6 survival models on\nTCGA-BRCA (5-fold CV)", GREY_F, GREY_E, fs=8.5)
    arrow(ax, 0.39, 0.505, 0.45, 0.505)
    box(ax, 0.63, 0.505, 0.36, 0.085,
        "Lock best model\n(Gradient-Boosted Survival)\nby internal CV C-index", GOLD_F, GOLD_E, fs=8.3, bold=True)
    arrow(ax, 0.63, 0.4625, 0.5, 0.410)

    box(ax, 0.5, 0.375, 0.84, 0.062,
        "Apply locked model WITHOUT refitting to 3 external cohorts vs official genefu PAM50-ROR",
        GREEN_F, GREEN_E, fs=9)
    arrow(ax, 0.5, 0.344, 0.5, 0.284)

    box(ax, 0.5, 0.247, 0.84, 0.074,
        "Pre-registered primary endpoint (locked 2026-05-16):\nTNBC meta-analyzed ΔC-index ≥ +0.03, p < 0.05",
        GREY_F, GREY_E, fs=9)
    arrow(ax, 0.5, 0.210, 0.5, 0.158)
    box(ax, 0.5, 0.117, 0.84, 0.082,
        "Primary endpoint NOT MET:  ΔC = +0.0144 (95% CI -0.047 to +0.076), p = 0.65",
        "#fdecea", RED, fs=9.6, bold=True)
    ax.text(0.5, 0.028,
            "Take-home: pathway-level transcriptomic ML is competitive with but does not exceed PAM50-ROR.",
            ha="center", fontsize=9.2, fontstyle="italic")

    fig.savefig(FIGURES / "fig_figure7_workflow.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure7_workflow.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure7_workflow.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure 7 written")


if __name__ == "__main__":
    main()
