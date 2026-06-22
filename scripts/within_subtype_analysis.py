#!/usr/bin/env python3
"""Stage 2: external within-subtype analysis."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    EXTERNAL_COHORTS,
    FIGURES,
    HEADLINE_MODEL,
    PAM50_MODEL,
    RESULTS,
    ensure_dirs,
    load_predictions,
    log_phase3,
    meta_from_rows,
    metric_cindex,
    paired_delta,
    valid_survival_subset,
)

SUBTYPES = ["LumA", "LumB", "Her2", "Basal", "Normal"]


def logrank_median_split(df: pd.DataFrame, risk_col: str) -> float:
    median = df[risk_col].median()
    high = df[risk_col] >= median
    if high.nunique() < 2:
        return np.nan
    result = logrank_test(
        df.loc[high, "os_days"],
        df.loc[~high, "os_days"],
        event_observed_A=df.loc[high, "os_event"],
        event_observed_B=df.loc[~high, "os_event"],
    )
    return float(result.p_value)


def plot_km(df: pd.DataFrame, cohort: str, subtype: str, risk_col: str) -> None:
    outdir = FIGURES / "km_within_subtype"
    outdir.mkdir(parents=True, exist_ok=True)
    median = df[risk_col].median()
    groups = {"High ML risk": df[risk_col] >= median, "Low ML risk": df[risk_col] < median}
    if any(mask.sum() == 0 for mask in groups.values()):
        log_phase3(f"Within-subtype KM skipped {cohort} {subtype}: median split produced an empty group.")
        return
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    kmf = KaplanMeierFitter()
    for label, mask in groups.items():
        kmf.fit(df.loc[mask, "os_days"] / 365.25, df.loc[mask, "os_event"], label=label)
        kmf.plot_survival_function(ax=ax, ci_show=False)
    ax.set_title(f"{cohort} {subtype}")
    ax.set_xlabel("Years")
    ax.set_ylabel("Overall survival")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    stem = outdir / f"km_{cohort}_{subtype}".replace("/", "_")
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".tiff"), dpi=300)
    plt.close(fig)


def run_analysis() -> tuple[pd.DataFrame, pd.DataFrame]:
    pred = load_predictions()
    rows = []
    for cohort in EXTERNAL_COHORTS:
        cohort_df = pred[pred["cohort"].eq(cohort)].copy()
        for subtype in SUBTYPES:
            sub = cohort_df[cohort_df["pam50_subtype_genefu"].eq(subtype)].copy()
            sub = valid_survival_subset(sub, [HEADLINE_MODEL, PAM50_MODEL])
            n = len(sub)
            events = int(sub["os_event"].sum()) if n else 0
            if n < 30 or events < 10:
                rows.append(
                    {
                        "subtype": subtype,
                        "cohort": cohort,
                        "n": n,
                        "events": events,
                        "ml_cindex": np.nan,
                        "ml_ci_low": np.nan,
                        "ml_ci_high": np.nan,
                        "pam50_cindex": np.nan,
                        "delta_vs_pam50": np.nan,
                        "delta_ci_low": np.nan,
                        "delta_ci_high": np.nan,
                        "logrank_p": np.nan,
                        "status": "skipped_low_n_or_events",
                    }
                )
                log_phase3(f"Within-subtype skipped {cohort} {subtype}: n={n}, events={events}.")
                continue
            ml_c, ml_low, ml_high = metric_cindex(sub, HEADLINE_MODEL)
            pam_c, _, _ = metric_cindex(sub, PAM50_MODEL)
            delta = paired_delta(sub, HEADLINE_MODEL, PAM50_MODEL)
            logrank_p = logrank_median_split(sub, HEADLINE_MODEL)
            plot_km(sub, cohort, subtype, HEADLINE_MODEL)
            rows.append(
                {
                    "subtype": subtype,
                    "cohort": cohort,
                    "n": n,
                    "events": events,
                    "ml_cindex": ml_c,
                    "ml_ci_low": ml_low,
                    "ml_ci_high": ml_high,
                    "pam50_cindex": pam_c,
                    "delta_vs_pam50": delta["delta"],
                    "delta_ci_low": delta["ci_low"],
                    "delta_ci_high": delta["ci_high"],
                    "logrank_p": logrank_p,
                    "status": "ok",
                }
            )
    table = pd.DataFrame(rows)
    meta_rows = []
    for subtype, sub in table[table["status"].eq("ok")].groupby("subtype"):
        meta = meta_from_rows(sub)
        directions = np.sign(sub["delta_vs_pam50"].dropna())
        consistent = bool(len(directions) > 0 and (directions.gt(0).all() or directions.lt(0).all()))
        meta_rows.append(
            {
                "subtype": subtype,
                "n_cohorts": meta["n_studies"],
                "meta_delta": meta["effect"],
                "meta_ci_low": meta["ci_low"],
                "meta_ci_high": meta["ci_high"],
                "meta_p": meta["p"],
                "i2": meta["i2"],
                "consistent_direction": consistent,
                "verdict": "positive_signal"
                if meta["effect"] > 0.02 and meta["p"] < 0.10 and consistent
                else "comparable_no_consistent_advantage",
            }
        )
    meta_table = pd.DataFrame(meta_rows)
    return table, meta_table


def plot_forest(table: pd.DataFrame, meta_table: pd.DataFrame) -> None:
    plot_rows = table[table["status"].eq("ok")].copy()
    fig, ax = plt.subplots(figsize=(7.2, max(4, 0.32 * (len(plot_rows) + len(meta_table) + 2))))
    y = 0
    yticks = []
    labels = []
    for subtype in SUBTYPES:
        sub = plot_rows[plot_rows["subtype"].eq(subtype)]
        if sub.empty:
            continue
        ax.text(-0.32, y, subtype, fontsize=9, fontweight="bold", va="center")
        yticks.append(y)
        labels.append(subtype)
        y += 1
        for _, row in sub.iterrows():
            ax.errorbar(
                row["delta_vs_pam50"],
                y,
                xerr=[[row["delta_vs_pam50"] - row["delta_ci_low"]], [row["delta_ci_high"] - row["delta_vs_pam50"]]],
                fmt="o",
                color="#31688e",
                capsize=2,
            )
            events = int(row["events"]) if "events" in row and pd.notna(row["events"]) else 0
            ax.text(-0.30, y, f"{row['cohort']} (n={int(row['n'])}, e={events})", fontsize=8, va="center")
            y += 1
        meta = meta_table[meta_table["subtype"].eq(subtype)]
        if not meta.empty:
            row = meta.iloc[0]
            ax.errorbar(
                row["meta_delta"],
                y,
                xerr=[[row["meta_delta"] - row["meta_ci_low"]], [row["meta_ci_high"] - row["meta_delta"]]],
                fmt="D",
                color="#440154",
                capsize=3,
            )
            ax.text(-0.30, y, f"Meta p={row['meta_p']:.3f}", fontsize=8, va="center")
            y += 1
        y += 0.5
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("Delta C-index: ML minus PAM50-ROR")
    ax.set_yticks([])
    ax.set_title("Exploratory within-subtype external validation")
    ax.set_xlim(-0.35, 0.35)
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "fig_within_subtype_forest.pdf")
    fig.savefig(FIGURES / "fig_within_subtype_forest.tiff", dpi=300)
    plt.close(fig)


def update_story(meta_table: pd.DataFrame) -> None:
    story = REPO_ROOT / "docs" / "STORY.md"
    text = story.read_text()
    if meta_table.empty or not meta_table["verdict"].eq("positive_signal").any():
        verdict = "comparable, with no consistent within-subtype advantage over PAM50-ROR."
    else:
        positives = ", ".join(meta_table.loc[meta_table["verdict"].eq("positive_signal"), "subtype"].tolist())
        verdict = f"exploratory positive signal in {positives}, requiring cautious interpretation."
    text = text.replace("- Within-subtype risk stratification: pending Stage 2.", f"- Within-subtype risk stratification: {verdict}")
    story.write_text(text)


def main() -> None:
    ensure_dirs()
    log_phase3("Stage 2 within-subtype analysis started.")
    table, meta_table = run_analysis()
    table["analysis_type"] = "exploratory_secondary"  # secondary subgroup analysis; not the primary endpoint
    meta_table["analysis_type"] = "exploratory_secondary"
    table.to_csv(RESULTS / "Table_S8_within_subtype_external.csv", index=False)
    meta_table.to_csv(RESULTS / "Table_S9_within_subtype_meta.csv", index=False)
    plot_forest(table, meta_table)
    update_story(meta_table)
    log_phase3("Stage 2 within-subtype analysis completed.")


if __name__ == "__main__":
    main()
