#!/usr/bin/env python3
"""Score the 2026-05-16 locked model on all eligible patients (Table S19)."""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    HEADLINE_MODEL,
    MODELS,
    PAM50_MODEL,
    RESULTS,
    load_modeling_frame,
    load_predictions,
    metric_cindex,
    paired_delta,
    valid_survival_subset,
)
from src.meta import random_effects_from_ci  # noqa: E402
from src.ml.wrappers import SkSurvAdapter  # noqa: E402,F401

EXTERNAL = ["GSE96058", "METABRIC", "GSE20685"]
LOCKED = MODELS / "locked_v1" / "gradient_boosted_survival.pkl"


def locked_risk(frame: pd.DataFrame) -> np.ndarray:
    with LOCKED.open("rb") as fh:
        artifact = pickle.load(fh)
    x = frame.copy()
    x["stage_ordinal"] = np.nan
    return np.asarray(artifact["model"].predict(artifact["preprocessor"].transform(x[artifact["features"]])), dtype=float)


def evaluate(df: pd.DataFrame, risk_col: str, label: str, cohort: str, population: str, subgroup: str) -> dict | None:
    sub = valid_survival_subset(df, [risk_col, PAM50_MODEL])
    if len(sub) < 10 or sub["os_event"].sum() < 2:
        return None
    c, low, high = metric_cindex(sub, risk_col)
    pam, _, _ = metric_cindex(sub, PAM50_MODEL)
    d = paired_delta(sub, risk_col, PAM50_MODEL)
    return {
        "model": label,
        "cohort": cohort,
        "population": population,
        "subgroup": subgroup,
        "n": len(sub),
        "events": int(sub["os_event"].sum()),
        "model_c": c,
        "model_c_ci_low": low,
        "model_c_ci_high": high,
        "pam50_rors_c": pam,
        "delta": d["delta"],
        "ci_low": d["ci_low"],
        "ci_high": d["ci_high"],
        "p": d["p"],
    }


def main() -> None:
    pred = load_predictions()
    frame = load_modeling_frame()
    frame = frame[frame["cohort"].isin(EXTERNAL)].copy()
    frame["locked_v1_gbsa"] = locked_risk(frame)
    pred = pred.merge(frame[["sample_id", "locked_v1_gbsa", "in_legacy_v1_cohort"]], on="sample_id", how="inner")

    rows = []
    for cohort in EXTERNAL:
        cdf = pred[pred["cohort"].eq(cohort)]
        populations = {"all_eligible": cdf}
        if cohort == "GSE96058":
            populations["legacy_v1_subset"] = cdf[cdf["in_legacy_v1_cohort"].astype(bool)]
            populations["excluded_from_legacy"] = cdf[~cdf["in_legacy_v1_cohort"].astype(bool)]
        for population, pdf in populations.items():
            tnbc = pdf["tnbc_flag"].fillna(False).astype(bool)
            for subgroup, sdf in {"overall": pdf, "tnbc": pdf[tnbc]}.items():
                for col, label in [("locked_v1_gbsa", "Locked v1 GBSA (TCGA n=213, 2026-05-16)"), (HEADLINE_MODEL, f"Revised {HEADLINE_MODEL} (all eligible TCGA)")]:
                    row = evaluate(sdf, col, label, cohort, population, subgroup)
                    if row:
                        rows.append(row)
    table = pd.DataFrame(rows)

    meta_rows = []
    for label, g in table[table["subgroup"].eq("tnbc") & table["population"].eq("all_eligible")].groupby("model"):
        m = random_effects_from_ci(g, "delta", "ci_low", "ci_high")
        meta_rows.append(
            {
                "model": label,
                "cohort": "Meta_RE_external",
                "population": "all_eligible",
                "subgroup": "tnbc",
                "n": int(g["n"].sum()),
                "events": int(g["events"].sum()),
                "delta": m["effect"],
                "ci_low": m["ci_low"],
                "ci_high": m["ci_high"],
                "p": m["p"],
            }
        )
    table = pd.concat([table, pd.DataFrame(meta_rows)], ignore_index=True)
    table.to_csv(RESULTS / "Table_S19_locked_model_full_cohorts.csv", index=False)
    print(table.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
