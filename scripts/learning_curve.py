#!/usr/bin/env python3
"""Training-size sensitivity / learning curve (reviewers R1.2 / R2.2 / R3.2).

Reviewers asked how the small TCGA-BRCA training cohort (n = 213) affects model
stability and external generalisation, and why TCGA was the sole training set.
This script:

1. Retrains the locked headline configuration (Gradient Boosted Survival,
   n_estimators=100, learning_rate=0.03, max_depth=1, subsample=1.0) on growing
   random fractions of TCGA-BRCA (20-100%), reporting internal 5-fold CV C-index
   and pooled external C-index as a function of training size, with bootstrap
   bands over repeated subsamples.
2. Runs a train-on-METABRIC (the largest cohort) cross-check, validating on the
   remaining cohorts, so the TCGA-as-training choice is justified empirically
   rather than only argued.

Outputs:
  results/Table_S17_learning_curve.csv
  results/Table_S17_train_cohort_crosscheck.csv
  figures/fig_learning_curve.png / .pdf

Exploratory, secondary analysis.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sksurv.ensemble import GradientBoostingSurvivalAnalysis  # noqa: E402
from sksurv.metrics import concordance_index_censored  # noqa: E402
from sksurv.util import Surv  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.train_ml_zoo import load_modeling_frame, make_preprocessor  # noqa: E402

RESULTS = REPO_ROOT / "results"
FIGURES = REPO_ROOT / "figures"
PROCESSED = REPO_ROOT / "data" / "processed"
SEED = 42
N_REPEATS = 25
N_SPLITS = 5
FRACTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
GBSA_CONFIG = {"n_estimators": 100, "learning_rate": 0.03, "max_depth": 1, "subsample": 1.0}
EXTERNAL_COHORTS = ["GSE96058", "METABRIC", "GSE20685"]


def valid_rows(frame: pd.DataFrame) -> pd.DataFrame:
    f = frame[frame["os_days"].notna() & frame["os_event"].notna() & (frame["os_days"] > 0)].copy()
    f["os_event"] = pd.to_numeric(f["os_event"], errors="coerce").fillna(0).astype(int)
    f["os_days"] = pd.to_numeric(f["os_days"], errors="coerce").astype(float)
    return f


def fit_gbsa(x: pd.DataFrame, time: np.ndarray, event: np.ndarray):
    pre = make_preprocessor().fit(x)
    model = GradientBoostingSurvivalAnalysis(random_state=SEED, **GBSA_CONFIG)
    model.fit(pre.transform(x), Surv.from_arrays(event.astype(bool), time.astype(float)))
    return pre, model


def cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    if len(time) < 5 or int(event.sum()) < 2:
        return np.nan
    return float(concordance_index_censored(event.astype(bool), time.astype(float), risk.astype(float))[0])


def internal_cv_c(x: pd.DataFrame, time: np.ndarray, event: np.ndarray, rng: np.random.Generator) -> float:
    """5-fold OOF Harrell C on a (sub)sample."""
    from sklearn.model_selection import StratifiedKFold

    seed = int(rng.integers(0, 2**31 - 1))
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
    oof = np.full(len(x), np.nan)
    for tr, te in skf.split(x, event):
        pre, model = fit_gbsa(x.iloc[tr], time[tr], event[tr])
        oof[te] = model.predict(pre.transform(x.iloc[te]))
    return cindex(time, event, oof)


def external_pooled_c(pre, model, externals: dict[str, pd.DataFrame], features: list[str]) -> float:
    cs, weights = [], []
    for _cohort, df in externals.items():
        risk = model.predict(pre.transform(df[features]))
        c = cindex(df["os_days"].to_numpy(), df["os_event"].to_numpy(), risk)
        if np.isfinite(c):
            cs.append(c)
            weights.append(int(df["os_event"].sum()))
    if not cs:
        return np.nan
    return float(np.average(cs, weights=weights))  # event-weighted pooled external C


def main() -> None:
    meta = json.loads((PROCESSED / "ml_model_zoo.metadata.json").read_text())
    features = meta["features"]
    frame = valid_rows(load_modeling_frame())

    tcga = frame[frame["cohort"].eq("TCGA-BRCA")].reset_index(drop=True)
    externals = {c: frame[frame["cohort"].eq(c)].reset_index(drop=True) for c in EXTERNAL_COHORTS}
    n_full = len(tcga)
    rng = np.random.default_rng(SEED)

    rows = []
    for frac in FRACTIONS:
        n = max(N_SPLITS * 2, int(round(frac * n_full)))
        internal_scores, external_scores = [], []
        for _ in range(N_REPEATS):
            idx = rng.choice(n_full, size=min(n, n_full), replace=False)
            sub = tcga.iloc[idx]
            x = sub[features]
            t = sub["os_days"].to_numpy()
            e = sub["os_event"].to_numpy()
            if int(e.sum()) < N_SPLITS:
                continue
            internal_scores.append(internal_cv_c(x.reset_index(drop=True), t, e, rng))
            pre, model = fit_gbsa(x, t, e)
            external_scores.append(external_pooled_c(pre, model, externals, features))
        internal_scores = [s for s in internal_scores if np.isfinite(s)]
        external_scores = [s for s in external_scores if np.isfinite(s)]
        rows.append(
            {
                "fraction": frac,
                "n_train": min(n, n_full),
                "internal_cv_c_mean": float(np.mean(internal_scores)) if internal_scores else np.nan,
                "internal_cv_c_lo": float(np.quantile(internal_scores, 0.025)) if internal_scores else np.nan,
                "internal_cv_c_hi": float(np.quantile(internal_scores, 0.975)) if internal_scores else np.nan,
                "external_pooled_c_mean": float(np.mean(external_scores)) if external_scores else np.nan,
                "external_pooled_c_lo": float(np.quantile(external_scores, 0.025)) if external_scores else np.nan,
                "external_pooled_c_hi": float(np.quantile(external_scores, 0.975)) if external_scores else np.nan,
                "n_repeats": len(external_scores),
                "analysis_type": "exploratory_secondary",
            }
        )

    table = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    table.to_csv(RESULTS / "Table_S17_learning_curve.csv", index=False)
    print("Learning curve:")
    print(table.to_string(index=False))

    # --- train-on-each-cohort cross-check ---
    cross_rows = []
    all_cohorts = ["TCGA-BRCA", *EXTERNAL_COHORTS]
    for train_cohort in all_cohorts:
        tr = frame[frame["cohort"].eq(train_cohort)]
        if int(tr["os_event"].sum()) < N_SPLITS:
            continue
        pre, model = fit_gbsa(tr[features], tr["os_days"].to_numpy(), tr["os_event"].to_numpy())
        val = {c: frame[frame["cohort"].eq(c)] for c in all_cohorts if c != train_cohort}
        per = {}
        for c, df in val.items():
            risk = model.predict(pre.transform(df[features]))
            per[c] = cindex(df["os_days"].to_numpy(), df["os_event"].to_numpy(), risk)
        pooled = external_pooled_c(pre, model, val, features)
        cross_rows.append(
            {
                "train_cohort": train_cohort,
                "n_train": len(tr),
                "events_train": int(tr["os_event"].sum()),
                "validation_pooled_c": pooled,
                **{f"val_c_{c}": per.get(c, np.nan) for c in all_cohorts},
                "analysis_type": "exploratory_secondary",
            }
        )
    cross = pd.DataFrame(cross_rows)
    cross.to_csv(RESULTS / "Table_S17_train_cohort_crosscheck.csv", index=False)
    print("\nTrain-cohort cross-check:")
    print(cross.to_string(index=False))

    # --- figure ---
    FIGURES.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(table["n_train"], table["internal_cv_c_mean"], "-o", color="#1f77b4", label="Internal 5-fold CV (TCGA)")
    ax.fill_between(table["n_train"], table["internal_cv_c_lo"], table["internal_cv_c_hi"], color="#1f77b4", alpha=0.15)
    ax.plot(table["n_train"], table["external_pooled_c_mean"], "-s", color="#d62728", label="External pooled (event-weighted)")
    ax.fill_between(table["n_train"], table["external_pooled_c_lo"], table["external_pooled_c_hi"], color="#d62728", alpha=0.15)
    ax.axhline(0.5, color="gray", ls=":", lw=1, label="Chance")
    ax.set_xlabel("TCGA-BRCA training samples (n)")
    ax.set_ylabel("Harrell C-index")
    ax.set_title("Learning curve: headline GBSA vs training size")
    ax.set_ylim(0.45, 0.75)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_learning_curve.png", dpi=200)
    fig.savefig(FIGURES / "fig_learning_curve.pdf")
    print(f"\nWrote {FIGURES / 'fig_learning_curve.png'} and .pdf")


if __name__ == "__main__":
    main()
