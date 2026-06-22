#!/usr/bin/env python3
"""Figure 1: four-cohort CONSORT-style inclusion flow.

Two evaluable denominators are reported separately and must not be conflated:
  - Characteristics-evaluable set (patient-characteristics tables): METABRIC 1,980 -> N = 4,003.
  - PAM50-ROR head-to-head C-index set: one METABRIC sample with a non-positive
    survival time is not C-index-evaluable, so METABRIC 1,979 -> N = 4,002.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

FIGURES = Path(__file__).resolve().parents[1] / "figures"

COHORTS = ["TCGA-BRCA", "GSE96058", "METABRIC", "GSE20685"]
INITIAL = [213, 1483, 2509, 327]          # harmonized samples; sum = 4,532
EXCLUDED = [0, 0, 529, 0]                  # missing OS / harmonization filters
CHARACTERISTICS = [213, 1483, 1980, 327]  # characteristics-evaluable; sum = 4,003
HEAD_TO_HEAD = [213, 1483, 1979, 327]      # C-index-evaluable; sum = 4,002


def n(v: int) -> str:
    return f"{v:,}"


def box(ax, x, y, text, fc, ec, w=0.205, h=0.135, bold=False):
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        fc=fc, ec=ec, lw=1.1,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=8.5,
            fontweight="bold" if bold else "normal")


def arrow(ax, x, y1, y2):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
                                 mutation_scale=12, lw=1.0, color="#555555"))


def main() -> None:
    fig, ax = plt.subplots(figsize=(13, 7.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    xs = [0.14, 0.38, 0.62, 0.86]
    ys = [0.88, 0.66, 0.44, 0.22]
    gray_fc, gray_ec = "#f2f2f2", "#444444"
    blue_fc, blue_ec = "#dbe9f6", "#2c7fb8"
    green_fc, green_ec = "#e5f5e0", "#31a354"

    for x, cohort, ini, exc, char, h2h in zip(
        xs, COHORTS, INITIAL, EXCLUDED, CHARACTERISTICS, HEAD_TO_HEAD, strict=False
    ):
        box(ax, x, ys[0], f"{cohort}\nInitial harmonized samples\nn = {n(ini)}", gray_fc, gray_ec)
        box(ax, x, ys[1], f"Excluded for missing OS\nor harmonization filters\nn = {n(exc)}", gray_fc, gray_ec)
        box(ax, x, ys[2], f"Final evaluable\n(patient characteristics)\nn = {n(char)}", blue_fc, blue_ec)
        box(ax, x, ys[3], f"Evaluable in PAM50-ROR\nhead-to-head (C-index)\nn = {n(h2h)}", green_fc, green_ec)
        arrow(ax, x, ys[0] - 0.068, ys[1] + 0.068)
        arrow(ax, x, ys[1] - 0.068, ys[2] + 0.068)
        arrow(ax, x, ys[2] - 0.068, ys[3] + 0.068)

    footer = (
        f"Total harmonized N = {n(sum(INITIAL))}     |     "
        f"Characteristics-evaluable N = {n(sum(CHARACTERISTICS))}     |     "
        f"PAM50-ROR head-to-head N = {n(sum(HEAD_TO_HEAD))}"
    )
    ax.text(0.5, 0.045, footer, ha="center", va="center", fontsize=12, fontweight="bold")

    fig.tight_layout()
    fig.savefig(FIGURES / "fig_figure1_consort.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure1_consort.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure1_consort.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure 1 written:", sum(INITIAL), sum(CHARACTERISTICS), sum(HEAD_TO_HEAD))


if __name__ == "__main__":
    main()
