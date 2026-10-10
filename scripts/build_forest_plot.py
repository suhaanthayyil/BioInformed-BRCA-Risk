#!/usr/bin/env python3
"""Stage 8: central BMC forest plot for overall and TNBC delta C-index."""

from __future__ import annotations

import json
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
    PROCESSED,
    RESULTS,
    ensure_dirs,
    log_phase3,
)
from src.meta import random_effects_from_ci  # noqa: E402

def subset_rows(table: pd.DataFrame, subgroup: str) -> tuple[pd.DataFrame, dict]:
    rows = table[
        table["headline_model"].eq(HEADLINE_MODEL)
        & table["comparator"].eq(PAM50_MODEL)
        & table["subgroup"].eq(subgroup)
    ].copy()
    rows = rows.rename(columns={"delta_cindex": "delta_vs_pam50", "ci_low": "delta_ci_low", "ci_high": "delta_ci_high"})
    if subgroup == "tnbc":
        primary = json.loads((PROCESSED / "phase6_primary_endpoint.json").read_text())
        meta = {k.removeprefix("meta_"): v for k, v in primary.items() if k.startswith("meta_")}
    else:
        meta = random_effects_from_ci(rows, "delta_vs_pam50", "delta_ci_low", "delta_ci_high")
    return rows, meta


def draw_panel(
    ax,
    rows: pd.DataFrame,
    meta: dict,
    title: str,
    primary_note: str | None = None,
) -> None:
    # Plot only evaluable cohorts (drop rows without a computable delta/CI, e.g.
    # GSE20685 in the TNBC panel, which lacks receptor-status annotation).
    rows = rows.dropna(subset=["delta_vs_pam50", "delta_ci_low", "delta_ci_high"]).sort_values("cohort")
    y_positions = np.arange(len(rows), 0, -1)
    lows = [*rows["delta_ci_low"].tolist(), meta["ci_low"]]
    highs = [*rows["delta_ci_high"].tolist(), meta["ci_high"]]
    dmin, dmax = min(lows), max(highs)
    pad = 0.1 * (dmax - dmin)
    labels, values = [], []
    for y, (_, row) in zip(y_positions, rows.iterrows(), strict=False):
        x, low, high = row["delta_vs_pam50"], row["delta_ci_low"], row["delta_ci_high"]
        ax.errorbar(x, y, xerr=[[x - low], [high - x]], fmt="o", color="#31688e", capsize=3)
        events = int(row["events"]) if np.isfinite(row["events"]) else 0
        name = f"{row['cohort']} (out-of-fold)" if row["cohort"] == "TCGA-BRCA" else row["cohort"]
        labels.append(f"{name}\nn={int(row['n']):,}, deaths={events:,}")
        values.append(f"{x:+.3f} [{low:+.3f}, {high:+.3f}]")
    total_events = int(rows["events"].fillna(0).sum())
    effect, low, high = meta["effect"], meta["ci_low"], meta["ci_high"]
    ax.errorbar(effect, 0, xerr=[[effect - low], [high - effect]], fmt="D", color="#440154", capsize=4)
    labels.append(f"Pooled (random effects)\ndeaths={total_events:,}")
    values.append(f"{effect:+.3f} [{low:+.3f}, {high:+.3f}]")
    ax.axvline(0, color="black", lw=1)
    ax.axvline(0.03, color="#d62728", ls="--", lw=0.9, label="+0.03 registered margin")
    ax.set_xlim(dmin - pad, dmax + pad)
    ax.set_ylim(-0.7, len(rows) + 0.7)
    ticks = [*y_positions, 0]
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=8)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    for y, text in zip(ticks, values, strict=True):
        ax.text(1.02, y, text, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
    ax.set_xlabel("Delta C-index, pathway ML minus genefu PAM50 ROR-S")
    p_text = f"p = {meta['p']:.3f}" if np.isfinite(meta.get("p", np.nan)) else "p = NA"
    i2, tau2 = meta.get("i2", np.nan), meta.get("tau2", np.nan)
    het = f"I2 = {100 * i2:.0f}%, tau2 = {tau2:.4f}" if np.isfinite(i2) and np.isfinite(tau2) else "I2 = NA"
    subtitle = f"{p_text}; {het}" + (f"\n{primary_note}" if primary_note else "")
    ax.set_title(f"{title}\n{subtitle}", loc="left", fontsize=9)
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    ax.grid(axis="x", alpha=0.2)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 8 main forest plot started.")
    table = pd.read_csv(RESULTS / "Table_2_head_to_head.csv")
    overall_rows, overall_meta = subset_rows(table, "overall")
    tnbc_rows, tnbc_meta = subset_rows(table, "tnbc")
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4))
    primary = json.loads((PROCESSED / "phase6_primary_endpoint.json").read_text())
    status = "MET" if primary["met"] else "NOT MET"
    draw_panel(axes[0], overall_rows, overall_meta, "A. Overall cohort")
    draw_panel(
        axes[1],
        tnbc_rows,
        tnbc_meta,
        "B. TNBC subgroup",
        f"Primary endpoint: {status}. GSE20685 has no receptor status.",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_main_forest.pdf")
    fig.savefig(FIGURES / "fig_main_forest.png", dpi=200)
    fig.savefig(FIGURES / "fig_main_forest.tiff", dpi=300)
    plt.close(fig)
    log_phase3("Stage 8 main forest plot completed.")


if __name__ == "__main__":
    main()
