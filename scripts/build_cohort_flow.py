#!/usr/bin/env python3
"""Stage 7: cohort flow diagram."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    ALL_COHORTS,
    FIGURES,
    PAM50_MODEL,
    PROCESSED,
    ensure_dirs,
    load_samples_survival,
    log_phase3,
)


def flow_counts() -> dict[str, dict[str, int]]:
    frame = load_samples_survival()
    pam = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam_ok = set(pam[pam["baseline"].eq(PAM50_MODEL) & pam["status"].eq("ok")]["sample_id"])
    counts = {}
    for cohort in ALL_COHORTS:
        sub = frame[frame["cohort"].eq(cohort)].copy()
        survival_valid = sub["os_days"].notna() & (sub["os_days"] > 0)
        pam50_valid = sub["sample_id"].isin(pam_ok)
        missing_survival = int((~survival_valid).sum())
        after_survival = sub[survival_valid]
        missing_pam50 = int((~after_survival["sample_id"].isin(pam_ok)).sum())
        final_n = int((survival_valid & pam50_valid).sum())
        counts[cohort] = {
            "initial": len(sub),
            "missing_survival": missing_survival,
            "missing_pam50": missing_pam50,
            "short_followup": 0,
            "final": final_n,
        }
    return counts


def draw_box(ax, x: float, y: float, text: str, width: float = 0.20, height: float = 0.12) -> None:
    rect = Rectangle((x - width / 2, y - height / 2), width, height, fc="#f7f7f7", ec="#333333", lw=0.8)
    ax.add_patch(rect)
    ax.text(x, y, text, ha="center", va="center", fontsize=7, wrap=True)


def draw_arrow(ax, x: float, y1: float, y2: float) -> None:
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=9, lw=0.8, color="#333333"))


def plot_flow(counts: dict[str, dict[str, int]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    xs = [0.125, 0.375, 0.625, 0.875]
    ys = [0.84, 0.62, 0.40, 0.18]
    for x, cohort in zip(xs, ALL_COHORTS, strict=False):
        c = counts[cohort]
        labels = [
            f"{cohort}\nInitial harmonized samples\nn={c['initial']}",
            f"Excluded missing OS\nn={c['missing_survival']}",
            f"Excluded missing official PAM50\nn={c['missing_pam50']}",
            f"Final analysis set\nn={c['final']}",
        ]
        for y, label in zip(ys, labels, strict=False):
            draw_box(ax, x, y, label)
        draw_arrow(ax, x, ys[0] - 0.07, ys[1] + 0.07)
        draw_arrow(ax, x, ys[1] - 0.07, ys[2] + 0.07)
        draw_arrow(ax, x, ys[2] - 0.07, ys[3] + 0.07)
    ax.text(0.5, 0.97, "Cohort inclusion flow for Phase Three BMC Artificial Intelligence analysis", ha="center", va="center", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_cohort_flow.pdf")
    fig.savefig(FIGURES / "fig_cohort_flow.tiff", dpi=300)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 7 cohort flow diagram started.")
    counts = flow_counts()
    plot_flow(counts)
    log_phase3("Stage 7 cohort flow diagram completed.")


if __name__ == "__main__":
    main()
