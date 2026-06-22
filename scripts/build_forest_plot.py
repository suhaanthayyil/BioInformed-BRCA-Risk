#!/usr/bin/env python3
"""Stage 8: central BMC forest plot for overall and TNBC delta C-index."""

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
    FIGURES,
    HEADLINE_MODEL,
    PAM50_MODEL,
    RESULTS,
    ensure_dirs,
    log_phase3,
)
from src.meta import random_effects_from_ci  # noqa: E402

LOCKED_TNBC = {"delta": 0.0144, "ci_low": -0.0472, "ci_high": 0.0760, "p": 0.6466}


def subset_rows(table: pd.DataFrame, subgroup: str) -> tuple[pd.DataFrame, dict]:
    rows = table[
        table["headline_model"].eq(HEADLINE_MODEL)
        & table["comparator"].eq(PAM50_MODEL)
        & table["subgroup"].eq(subgroup)
    ].copy()
    rows = rows.rename(columns={"delta_cindex": "delta_vs_pam50", "ci_low": "delta_ci_low", "ci_high": "delta_ci_high"})
    if subgroup == "tnbc":
        meta = {
            "n_studies": len(rows),
            "effect": LOCKED_TNBC["delta"],
            "ci_low": LOCKED_TNBC["ci_low"],
            "ci_high": LOCKED_TNBC["ci_high"],
            "p": LOCKED_TNBC["p"],
            "i2": np.nan,
            "tau2": np.nan,
        }
    else:
        meta = random_effects_from_ci(rows, "delta_vs_pam50", "delta_ci_low", "delta_ci_high")
    return rows, meta


def draw_panel(ax, rows: pd.DataFrame, meta: dict, title: str, primary_note: str | None = None) -> None:
    rows = rows.sort_values("cohort")
    y_positions = np.arange(len(rows), 0, -1)
    for y, (_, row) in zip(y_positions, rows.iterrows(), strict=False):
        x = row["delta_vs_pam50"]
        low = row["delta_ci_low"]
        high = row["delta_ci_high"]
        ax.errorbar(x, y, xerr=[[x - low], [high - x]], fmt="o", color="#31688e", capsize=3)
        events = int(row["events"]) if "events" in row and np.isfinite(row["events"]) else 0
        ax.text(-0.185, y, f"{row['cohort']} (n={int(row['n'])}, e={events})", ha="left", va="center", fontsize=8)
        ax.text(0.185, y, f"{x:+.3f} [{low:+.3f}, {high:+.3f}]", ha="right", va="center", fontsize=8)
    total_events = int(rows["events"].fillna(0).sum()) if "events" in rows else 0
    meta_y = 0
    effect = meta["effect"]
    low = meta["ci_low"]
    high = meta["ci_high"]
    ax.errorbar(effect, meta_y, xerr=[[effect - low], [high - effect]], fmt="D", color="#440154", capsize=4)
    ax.text(-0.185, meta_y, f"Pooled (random-effects, e={total_events})", ha="left", va="center", fontsize=8, fontweight="bold")
    ax.text(0.185, meta_y, f"{effect:+.3f} [{low:+.3f}, {high:+.3f}]", ha="right", va="center", fontsize=8)
    ax.axvline(0, color="black", lw=1)
    ax.set_xlim(-0.20, 0.20)
    ax.set_ylim(-0.8, len(rows) + 1)
    ax.set_yticks([])
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.set_xlabel("Delta C-index, pathway ML minus PAM50-ROR")
    p_text = f"p={meta['p']:.4f}" if np.isfinite(meta.get("p", np.nan)) else "p=NA"
    i2 = meta.get("i2", np.nan)
    tau2 = meta.get("tau2", np.nan)
    heterogeneity = (
        f"I2={100 * i2:.1f}%, tau2={tau2:.4f}" if np.isfinite(i2) and np.isfinite(tau2) else "I2=NA, tau2=NA"
    )
    note = f"{p_text}; {heterogeneity}"
    if primary_note:
        note = f"{note}\n{primary_note}"
    ax.text(0.02, 0.96, note, transform=ax.transAxes, ha="left", va="top", fontsize=8)
    ax.grid(axis="x", alpha=0.2)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 8 main forest plot started.")
    table = pd.read_csv(RESULTS / "Table_2_head_to_head.csv")
    overall_rows, overall_meta = subset_rows(table, "overall")
    tnbc_rows, tnbc_meta = subset_rows(table, "tnbc")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
    draw_panel(axes[0], overall_rows, overall_meta, "A. Overall cohort")
    draw_panel(
        axes[1],
        tnbc_rows,
        tnbc_meta,
        "B. TNBC subgroup",
        "Pre-registered primary endpoint. Result: NOT MET.",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_main_forest.pdf")
    fig.savefig(FIGURES / "fig_main_forest.tiff", dpi=300)
    plt.close(fig)
    log_phase3("Stage 8 main forest plot completed.")


if __name__ == "__main__":
    main()
