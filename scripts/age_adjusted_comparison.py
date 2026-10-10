#!/usr/bin/env python3
"""Exploratory: pathway ML vs PAM50 ROR-S when both are given age (Table S21)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from lifelines import CoxPHFitter

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.phase3_common import (  # noqa: E402
    HEADLINE_MODEL,
    PAM50_MODEL,
    RESULTS,
    load_predictions,
    meta_from_rows,
    metric_cindex,
    paired_delta,
    valid_survival_subset,
)

EXTERNAL = ["GSE96058", "METABRIC", "GSE20685"]
ML = "Age_plus_ML"
PAM = "Age_plus_PAM50_RORS"
AGE = "Age_only"


def fit(train: pd.DataFrame, cols: list[str]) -> CoxPHFitter:
    return CoxPHFitter(penalizer=0.01).fit(train[[*cols, "os_days", "os_event"]], "os_days", "os_event")


def main() -> None:
    pred = load_predictions()
    pred = pred.dropna(subset=[HEADLINE_MODEL, PAM50_MODEL, "age_at_dx", "os_days"])
    pred = pred[pred["os_days"] > 0].copy()
    train = pred[pred["cohort"].eq("TCGA-BRCA")]
    models = {
        ML: fit(train, ["age_at_dx", HEADLINE_MODEL]),
        PAM: fit(train, ["age_at_dx", PAM50_MODEL]),
        AGE: fit(train, ["age_at_dx"]),
    }
    cols = {ML: ["age_at_dx", HEADLINE_MODEL], PAM: ["age_at_dx", PAM50_MODEL], AGE: ["age_at_dx"]}
    for name, model in models.items():
        pred[name] = model.predict_partial_hazard(pred[cols[name]]).to_numpy(float)

    rows = []
    for cohort in EXTERNAL:
        cdf = pred[pred["cohort"].eq(cohort)]
        tnbc = cdf["tnbc_flag"].fillna(False).astype(bool)
        for subgroup, sub in {"overall": cdf, "tnbc": cdf[tnbc]}.items():
            sub = valid_survival_subset(sub, [ML, PAM, AGE])
            if len(sub) < 30 or sub["os_event"].sum() < 10:
                continue
            ml_c, _, _ = metric_cindex(sub, ML)
            pam_c, _, _ = metric_cindex(sub, PAM)
            age_c, _, _ = metric_cindex(sub, AGE)
            d = paired_delta(sub, ML, PAM)
            rows.append(
                {
                    "cohort": cohort,
                    "subgroup": subgroup,
                    "n": len(sub),
                    "events": int(sub["os_event"].sum()),
                    "age_only_c": age_c,
                    "age_plus_pam50_c": pam_c,
                    "age_plus_ml_c": ml_c,
                    "delta_vs_pam50": d["delta"],
                    "delta_ci_low": d["ci_low"],
                    "delta_ci_high": d["ci_high"],
                    "p": d["p"],
                }
            )
    table = pd.DataFrame(rows)
    pooled = []
    for subgroup, g in table.groupby("subgroup"):
        m = meta_from_rows(g)
        pooled.append(
            {
                "cohort": "Meta_RE_external",
                "subgroup": subgroup,
                "n": int(g["n"].sum()),
                "events": int(g["events"].sum()),
                "delta_vs_pam50": m["effect"],
                "delta_ci_low": m["ci_low"],
                "delta_ci_high": m["ci_high"],
                "p": m["p"],
            }
        )
    table = pd.concat([table, pd.DataFrame(pooled)], ignore_index=True)
    table["analysis_type"] = "exploratory_post_hoc"
    table.to_csv(RESULTS / "Table_S21_age_adjusted_comparison.csv", index=False)
    print(table.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
