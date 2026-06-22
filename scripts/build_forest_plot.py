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

LOCKED_TNBC = {"delta": 0.0144, "ci_low": -0.0472, "ci_high": 0.0760, "p": 0.6466, "i2": 0.1053, "tau2": 0.00045}


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
            "i2": LOCKED_TNBC["i2"],
            "tau2": LOCKED_TNBC["tau2"],
        }
    else:
        meta = random_effects_from_ci(rows, "delta_vs_pam50", "delta_ci_low", "delta_ci_high")
    return rows, meta


def draw_panel(ax, rows: pd.DataFrame, meta: dict, title: str, primary_note: str | None = None) -> None:
    # Plot only evaluable cohorts (drop rows without a computable delta/CI, e.g.
    # GSE20685 in the TNBC panel, which lacks receptor-status annotation).
    rows = rows.dropna(subset=["delta_vs_pam50", "delta_ci_low", "delta_ci_high"]).sort_values("cohort")
    y_positions = np.arange(len(rows), 0, -1)
    # Adaptive x-range: span the full data (incl. the pooled estimate) so no CI
    # whisker is clipped; reserve left/right margins for the label columns.
    lows = [*rows["delta_ci_low"].tolist(), meta["ci_low"]]
    highs = [*rows["delta_ci_high"].tolist(), meta["ci_high"]]
    dmin, dmax = min(lows), max(highs)
    center, half = (dmin + dmax) / 2.0, max((dmax - dmin) / 2.0, 0.05)
    xlo, xhi = center - 2.0 * half, center + 2.0 * half  # data occupies the middle ~50%
    label_x, value_x = xlo + 0.02 * (xhi - xlo), xhi - 0.02 * (xhi - xlo)
    for y, (_, row) in zip(y_positions, rows.iterrows(), strict=False):
        x, low, high = row["delta_vs_pam50"], row["delta_ci_low"], row["delta_ci_high"]
        ax.errorbar(x, y, xerr=[[x - low], [high - x]], fmt="o", color="#31688e", capsize=3)
        events = int(row["events"]) if "events" in row and np.isfinite(row["events"]) else 0
        ax.text(label_x, y, f"{row['cohort']} (n={int(row['n'])}, e={events})", ha="left", va="center", fontsize=8)
        ax.text(value_x, y, f"{x:+.3f} [{low:+.3f}, {high:+.3f}]", ha="right", va="center", fontsize=8)
    total_events = int(rows["events"].fillna(0).sum()) if "events" in rows else 0
    effect, low, high = meta["effect"], meta["ci_low"], meta["ci_high"]
    ax.errorbar(effect, 0, xerr=[[effect - low], [high - effect]], fmt="D", color="#440154", capsize=4)
    ax.text(label_x, 0, f"Pooled (random-effects, e={total_events})", ha="left", va="center", fontsize=8, fontweight="bold")
    ax.text(value_x, 0, f"{effect:+.3f} [{low:+.3f}, {high:+.3f}]", ha="right", va="center", fontsize=8)
    ax.axvline(0, color="black", lw=1)
    ax.axvline(0.03, color="#d62728", ls="--", lw=0.9)  # pre-registered +0.03 superiority margin
    ax.text(0.03, len(rows) + 0.55, "+0.03 margin", color="#d62728", ha="center", va="bottom", fontsize=7)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-0.8, len(rows) + 1.1)
    ax.set_yticks([])
    ax.set_xticks([t for t in np.round(np.arange(-0.3, 0.41, 0.1), 1) if dmin - 0.02 <= t <= dmax + 0.02])
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.set_xlabel("Delta C-index, pathway ML minus PAM50-ROR")
    p_text = f"p={meta['p']:.4f}" if np.isfinite(meta.get("p", np.nan)) else "p=NA"
    i2, tau2 = meta.get("i2", np.nan), meta.get("tau2", np.nan)
    heterogeneity = (
        f"I2={100 * i2:.1f}%, tau2={tau2:.4f}" if np.isfinite(i2) and np.isfinite(tau2) else "I2=NA, tau2=NA"
    )
    note = f"{p_text}; {heterogeneity}"
    if primary_note:
        note = f"{note}\n{primary_note}"
    ax.text(0.02, 0.97, note, transform=ax.transAxes, ha="left", va="top", fontsize=8)
    ax.grid(axis="x", alpha=0.2)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 8 main forest plot started.")
    table = pd.read_csv(RESULTS / "Table_2_head_to_head.csv")
    overall_rows, overall_meta = subset_rows(table, "overall")
    tnbc_rows, tnbc_meta = subset_rows(table, "tnbc")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    draw_panel(axes[0], overall_rows, overall_meta, "A. Overall cohort")
    draw_panel(
        axes[1],
        tnbc_rows,
        tnbc_meta,
        "B. TNBC subgroup",
        "Pre-registered primary endpoint. Result: NOT MET.\nGSE20685 lacked receptor status and is not evaluable for TNBC.",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_main_forest.pdf")
    fig.savefig(FIGURES / "fig_main_forest.tiff", dpi=300)
    plt.close(fig)
    log_phase3("Stage 8 main forest plot completed.")


if __name__ == "__main__":
    main()
