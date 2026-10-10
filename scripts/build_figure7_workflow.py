#!/usr/bin/env python3
"""Figure 7 / graphical abstract: study design and data-processing workflow.

A single schematic summarising cohort harmonization, pathway scoring, model
training, external validation, and the registered primary result.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"

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
    flow = json.loads((ROOT / "results" / "Figure1_cohort_flow_counts.json").read_text())
    primary = json.loads((ROOT / "data" / "processed" / "phase6_primary_endpoint.json").read_text())
    total = sum(v["analysed"] for v in flow.values())
    fig, ax = plt.subplots(figsize=(9.5, 8.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.965, f"Four public breast cancer cohorts  (analysed N = {total:,})",
            ha="center", fontsize=10.5, fontweight="bold")
    cohorts = [(f"TCGA-BRCA\nn = {flow['TCGA-BRCA']['analysed']:,}\n(training)", 0.16),
               (f"GSE96058 / SCAN-B\nn = {flow['GSE96058']['analysed']:,}", 0.385),
               (f"METABRIC\nn = {flow['METABRIC']['analysed']:,}", 0.62),
               (f"GSE20685\nn = {flow['GSE20685']['analysed']:,}", 0.85)]
    for txt, x in cohorts:
        box(ax, x, 0.875, 0.205, 0.078, txt, GREY_F, GREY_E, fs=8)
        arrow(ax, x, 0.836, x, 0.781)

    box(ax, 0.5, 0.742, 0.84, 0.072,
        "Per-cohort gene-level z-scoring (no scaling transfer)  →  7 MSigDB pathway activity scores",
        BLUE_F, BLUE_E, fs=9)
    arrow(ax, 0.5, 0.706, 0.5, 0.671)
    box(ax, 0.5, 0.635, 0.84, 0.062,
        "Model inputs: the 7 pathway scores only (as registered); age and stage in sensitivity analyses",
        BLUE_F, BLUE_E, fs=9)
    arrow(ax, 0.5, 0.604, 0.30, 0.548)

    box(ax, 0.21, 0.505, 0.36, 0.085,
        "Train 6 survival models on\nTCGA-BRCA (5-fold CV)", GREY_F, GREY_E, fs=8.5)
    arrow(ax, 0.39, 0.505, 0.45, 0.505)
    box(ax, 0.63, 0.505, 0.36, 0.085,
        f"Fix best non-linear model\n({primary['headline_model'].replace('_', ' ')})\nby internal CV C-index", GOLD_F, GOLD_E, fs=8.3, bold=True)
    arrow(ax, 0.63, 0.4625, 0.5, 0.410)

    box(ax, 0.5, 0.375, 0.84, 0.062,
        "Apply fixed model WITHOUT refitting to 3 external cohorts\nvs genefu PAM50 ROR-S (research implementation); TCGA scored out-of-fold",
        GREEN_F, GREEN_E, fs=9)
    arrow(ax, 0.5, 0.344, 0.5, 0.284)

    box(ax, 0.5, 0.247, 0.84, 0.074,
        "Primary endpoint (registered 2026-05-16; analysis revised, see Table S20):\nTNBC meta-analyzed ΔC-index ≥ +0.03, p < 0.05",
        GREY_F, GREY_E, fs=9)
    arrow(ax, 0.5, 0.210, 0.5, 0.158)
    status = "MET" if primary["met"] else "NOT MET"
    box(ax, 0.5, 0.117, 0.84, 0.082,
        f"Primary endpoint {status}:  ΔC = {primary['meta_effect']:+.3f} "
        f"(95% CI {primary['meta_ci_low']:+.3f} to {primary['meta_ci_high']:+.3f}), p = {primary['meta_p']:.2f}",
        "#fdecea", RED, fs=9.6, bold=True)

    fig.savefig(FIGURES / "fig_figure7_workflow.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure7_workflow.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure7_workflow.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure 7 written")


if __name__ == "__main__":
    main()
