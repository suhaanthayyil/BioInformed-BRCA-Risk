#!/usr/bin/env python3
"""Recurrence-endpoint sensitivity analysis (reviewers R2.1 / R3.1).

PAM50-ROR is a risk-of-recurrence tool, while the primary analysis benchmarks
discrimination on overall survival (OS) because OS is the endpoint available
across all four public cohorts. METABRIC is the one cohort that carries a
disease-free-survival (DFS / recurrence) endpoint for the same samples that
have pathway features and an official PAM50-ROR score (n = 1,980).

This script re-runs the locked Gradient-Boosted-Survival model versus official
PAM50-ROR on METABRIC using the DFS endpoint (and OS on the same samples, for a
like-for-like contrast), overall and in the TNBC subgroup. It is an exploratory,
secondary sensitivity analysis — not the pre-registered primary endpoint.

Output: results/Table_S16_metabric_recurrence_sensitivity.csv

GSE96058 and GSE20685 do not provide recurrence fields in the parsed metadata,
so this analysis is METABRIC-only by necessity.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.external_validation import (  # noqa: E402
    PAM50_COMPARATOR,
    load_frame,
    load_pickle,
    predict_artifact,
)
from src.survival import (  # noqa: E402
    bootstrap_c_index_ci,
    harrell_c_index,
    paired_bootstrap_delta_c_index,
)

RESULTS = REPO_ROOT / "results"
MODELS = REPO_ROOT / "models"
PROCESSED = REPO_ROOT / "data" / "processed"
SEED = 42
BOOTSTRAP_ITER = 1000
HEADLINE = "Gradient_Boosted_Survival"
COHORT = "METABRIC"
ENDPOINTS = [("OS", "os_days", "os_event"), ("DFS", "dfs_days", "dfs_event")]


def metric(df: pd.DataFrame, time_col: str, event_col: str, risk_col: str) -> dict:
    sub = df[
        df[risk_col].notna()
        & df[time_col].notna()
        & df[event_col].notna()
        & (df[time_col] > 0)
    ]
    n = len(sub)
    events = int(pd.to_numeric(sub[event_col], errors="coerce").fillna(0).sum()) if n else 0
    if n < 10 or events < 2:
        return {"n": n, "events": events, "c": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    c = harrell_c_index(sub[time_col], sub[event_col], sub[risk_col])
    lo, hi = bootstrap_c_index_ci(sub[time_col], sub[event_col], sub[risk_col], BOOTSTRAP_ITER, SEED)
    return {"n": n, "events": events, "c": float(c), "ci_low": float(lo), "ci_high": float(hi)}


def main() -> None:
    frame = load_frame()
    mb = frame[frame["cohort"].eq(COHORT)].copy()
    if mb.empty:
        raise SystemExit("No METABRIC samples in the harmonized frame.")

    import json

    meta = json.loads((PROCESSED / "ml_model_zoo.metadata.json").read_text())
    features = meta["features"]

    artifact = load_pickle(MODELS / "gradient_boosted_survival.pkl")
    mb[HEADLINE] = predict_artifact(HEADLINE, artifact, mb[features])

    pam = pd.read_parquet(PROCESSED / "baselines_pam50.parquet")
    pam = (
        pam[pam["baseline"].eq(PAM50_COMPARATOR) & pam["status"].eq("ok")][["sample_id", "score"]]
        .rename(columns={"score": PAM50_COMPARATOR})
    )
    mb = mb.merge(pam, on="sample_id", how="left")

    tnbc = mb["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
    subgroups = {"overall": pd.Series(True, index=mb.index), "tnbc": tnbc}

    rows = []
    for endpoint, time_col, event_col in ENDPOINTS:
        for subgroup, mask in subgroups.items():
            sub = mb.loc[mask].copy()
            gb = metric(sub, time_col, event_col, HEADLINE)
            pm = metric(sub, time_col, event_col, PAM50_COMPARATOR)
            valid = sub[[HEADLINE, PAM50_COMPARATOR, time_col, event_col]].dropna()
            valid = valid[valid[time_col] > 0]
            if len(valid) >= 10 and int(valid[event_col].sum()) >= 2:
                delta = paired_bootstrap_delta_c_index(
                    valid[time_col], valid[event_col], valid[HEADLINE], valid[PAM50_COMPARATOR],
                    n_bootstrap=BOOTSTRAP_ITER, seed=SEED,
                )
            else:
                delta = {"delta": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
            rows.append(
                {
                    "cohort": COHORT,
                    "endpoint": endpoint,
                    "subgroup": subgroup,
                    "n": gb["n"],
                    "events": gb["events"],
                    "gbsa_harrell_c": gb["c"],
                    "gbsa_ci_low": gb["ci_low"],
                    "gbsa_ci_high": gb["ci_high"],
                    "pam50_ror_harrell_c": pm["c"],
                    "pam50_ror_ci_low": pm["ci_low"],
                    "pam50_ror_ci_high": pm["ci_high"],
                    "delta_cindex": delta["delta"],
                    "delta_ci_low": delta["ci_low"],
                    "delta_ci_high": delta["ci_high"],
                    "delta_p": delta["p"],
                    "comparison": f"{HEADLINE} minus {PAM50_COMPARATOR}",
                    "analysis_type": "exploratory_secondary",
                }
            )

    out = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    dest = RESULTS / "Table_S16_metabric_recurrence_sensitivity.csv"
    out.to_csv(dest, index=False)
    print(f"Wrote {dest} ({len(out)} rows)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
