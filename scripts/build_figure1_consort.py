#!/usr/bin/env python3
"""Figure 1: four-cohort CONSORT-style inclusion flow, built from the data."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    FIGURES,
    HEADLINE_MODEL,
    PAM50_MODEL,
    PROCESSED,
    RESULTS,
    load_modeling_frame,
    load_samples_survival,
)

COHORTS = ["TCGA-BRCA", "GSE96058", "METABRIC", "GSE20685"]


def n(v: int) -> str:
    return f"{v:,}"


def flow() -> dict[str, dict]:
    ingest = json.loads((PROCESSED / "ingest_flow.json").read_text())
    harmonized = load_samples_survival()
    features = load_modeling_frame()
    pred = pd.read_csv(RESULTS / "06_external_model_predictions.csv")
    evaluable = pred[pred[HEADLINE_MODEL].notna() & pred[PAM50_MODEL].notna() & pred["os_days"].gt(0)]

    out = {}
    t = ingest["TCGA-BRCA"]
    out["TCGA-BRCA"] = {
        "source": f"Primary-tumour RNA-seq\n(UCSC Xena / TCGA-CDR)\nn = {n(t['unique_patients'])} patients",
        "excluded": [
            ("male or sex missing", t["excluded_male_or_sex_missing"]),
            ("TCGA-CDR redacted", t["excluded_redacted"]),
            ("missing or zero OS time", t["excluded_missing_or_zero_os"]),
        ],
    }
    g = ingest["GSE96058"]
    out["GSE96058"] = {
        "source": f"GEO samples\nn = {n(g['geo_samples'])}",
        "excluded": [
            ("technical replicates", g["excluded_technical_replicates"]),
            ("missing OS", g["excluded_missing_or_zero_os"]),
        ],
    }
    m = harmonized[harmonized["cohort"].eq("METABRIC")]
    m_expr = features[features["cohort"].eq("METABRIC")]
    out["METABRIC"] = {
        "source": f"cBioPortal clinical records\nn = {n(len(m))}",
        "excluded": [("no expression profile", len(m) - len(m_expr))],
    }
    e = harmonized[harmonized["cohort"].eq("GSE20685")]
    out["GSE20685"] = {"source": f"GEO samples\nn = {n(len(e))}", "excluded": [("none", 0)]}

    for cohort in COHORTS:
        analysed = features[features["cohort"].eq(cohort)]
        ev = evaluable[evaluable["cohort"].eq(cohort)]
        out[cohort]["analysed"] = len(analysed)
        out[cohort]["analysed_events"] = int(analysed["os_event"].sum())
        out[cohort]["evaluable"] = len(ev)
        out[cohort]["evaluable_events"] = int(ev["os_event"].sum())
        out[cohort]["evaluable_excluded"] = len(analysed) - len(ev)
    (RESULTS / "Figure1_cohort_flow_counts.json").write_text(json.dumps(out, indent=2))
    return out


def box(ax, x, y, text, fc, ec, w=0.215, h=0.15):
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.006,rounding_size=0.012", fc=fc, ec=ec, lw=1.1
        )
    )
    ax.text(x, y, text, ha="center", va="center", fontsize=8.2)


def arrow(ax, x, y1, y2):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=12, lw=1.0, color="#555555"))


def main() -> None:
    counts = flow()
    fig, ax = plt.subplots(figsize=(13, 8.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    xs = [0.13, 0.38, 0.63, 0.88]
    ys = [0.89, 0.67, 0.45, 0.23]
    gray, blue, green = ("#f2f2f2", "#444444"), ("#dbe9f6", "#2c7fb8"), ("#e5f5e0", "#31a354")

    for x, cohort in zip(xs, COHORTS, strict=True):
        c = counts[cohort]
        excluded = "\n".join(f"{label}: {n(v)}" for label, v in c["excluded"])
        box(ax, x, ys[0], f"{cohort}\n{c['source']}", *gray)
        box(ax, x, ys[1], f"Excluded\n{excluded}", *gray)
        box(ax, x, ys[2], f"Analysed (pathway scores + OS)\nn = {n(c['analysed'])}\ndeaths = {n(c['analysed_events'])}", *blue)
        box(
            ax,
            x,
            ys[3],
            f"Evaluable vs PAM50 ROR-S\n(C-index; OS time > 0)\nn = {n(c['evaluable'])}, deaths = {n(c['evaluable_events'])}",
            *green,
        )
        for a, b in zip(ys[:-1], ys[1:], strict=True):
            arrow(ax, x, a - 0.075, b + 0.075)

    total = sum(counts[c]["analysed"] for c in COHORTS)
    total_eval = sum(counts[c]["evaluable"] for c in COHORTS)
    ax.text(
        0.5,
        0.06,
        f"Analysed N = {n(total)}     |     Evaluable vs PAM50 ROR-S N = {n(total_eval)}\n"
        "No patient was excluded on the basis of follow-up length; censored patients are retained at their last follow-up.",
        ha="center",
        va="center",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_figure1_consort.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure1_consort.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "fig_figure1_consort.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "source"} for k, v in counts.items()}, indent=1))


if __name__ == "__main__":
    main()
