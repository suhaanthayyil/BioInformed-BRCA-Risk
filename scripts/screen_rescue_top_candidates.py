#!/usr/bin/env python3
"""Post-hoc external screen for near-top rescue candidates.

This script does not choose the rescue headline. It exists only to document
whether any near-top TCGA-selected candidate looks externally promising.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.run_rescue_analysis import (  # noqa: E402
    PAM50_COMPARATOR,
    PRIMARY_COHORTS,
    TRAINING_COHORT,
    CandidateResult,
    add_pam50,
    feature_sets,
    fit_final_candidate,
    load_base_frame,
    model_spaces,
    predict_artifact,
    usable_features,
)
from src.survival import harrell_c_index  # noqa: E402


RESULTS = REPO_ROOT / "results"


def rebuild_candidate(row: pd.Series, train: pd.DataFrame, sets: dict[str, list[str]]) -> CandidateResult:
    config = json.loads(row["config"])
    features = usable_features(train, sets[row["feature_set"]])
    return CandidateResult(
        feature_set=row["feature_set"],
        model=row["model"],
        config=config,
        features=features,
        overall_cv_cindex=float(row["overall_cv_cindex"]),
        tnbc_cv_cindex=float(row["tnbc_cv_cindex"]),
        valid_tnbc_folds=int(row["valid_tnbc_folds"]),
        composite_cv_cindex=float(row["composite_cv_cindex"]),
        fold_cindexes=[],
        fold_tnbc_cindexes=[],
        oof_risk=np.array([]),
    )


def screen_candidate(label: str, candidate: CandidateResult, frame: pd.DataFrame, train: pd.DataFrame, spaces: dict) -> dict:
    artifact = fit_final_candidate(candidate, train, spaces)
    pred = frame[["sample_id", "cohort", "os_days", "os_event", "tnbc_flag"]].copy()
    pred[label] = predict_artifact(artifact, frame)
    pred = add_pam50(pred)

    rows = []
    for cohort, sub in pred[pred["cohort"].isin(PRIMARY_COHORTS)].groupby("cohort"):
        tnbc = sub["tnbc_flag"].astype("boolean").fillna(False).astype(bool)
        scored = sub[tnbc].dropna(subset=[label, PAM50_COMPARATOR, "os_days", "os_event"])
        if len(scored) < 10 or scored["os_event"].sum() < 2:
            continue
        ml_c = harrell_c_index(scored["os_days"], scored["os_event"], scored[label])
        pam50_c = harrell_c_index(scored["os_days"], scored["os_event"], scored[PAM50_COMPARATOR])
        rows.append(
            {
                "cohort": cohort,
                "n": len(scored),
                "events": int(scored["os_event"].sum()),
                "delta_cindex": ml_c - pam50_c,
                "ml_c": ml_c,
                "pam50_c": pam50_c,
            }
        )

    all_delta = float(np.mean([row["delta_cindex"] for row in rows])) if rows else np.nan
    external = [row["delta_cindex"] for row in rows if row["cohort"] != TRAINING_COHORT]
    external_delta = float(np.mean(external)) if external else np.nan
    return {
        "label": label,
        "feature_set": candidate.feature_set,
        "model": candidate.model,
        "overall_cv_cindex": candidate.overall_cv_cindex,
        "tnbc_cv_cindex": candidate.tnbc_cv_cindex,
        "composite_cv_cindex": candidate.composite_cv_cindex,
        "mean_delta_all_protocol_tnbc": all_delta,
        "mean_delta_external_only_tnbc": external_delta,
        "rows": json.dumps(rows),
    }


def main() -> None:
    frame = load_base_frame()
    train = frame[
        frame["cohort"].eq(TRAINING_COHORT)
        & frame["os_days"].notna()
        & frame["os_event"].notna()
        & (frame["os_days"] > 0)
    ].copy()
    sets = feature_sets(frame)
    spaces = model_spaces()
    cv = pd.read_csv(RESULTS / "Table_S4_rescue_internal_cv.csv").head(12)

    rows = []
    for i, (_, row) in enumerate(cv.iterrows()):
        if row["model"] == "Rank_Average_Ensemble":
            continue
        candidate = rebuild_candidate(row, train, sets)
        label = f"Posthoc_{row['model']}_{row['feature_set'][:10]}_{i}"
        rows.append(screen_candidate(label, candidate, frame, train, spaces))

    out = pd.DataFrame(rows).sort_values("mean_delta_external_only_tnbc", ascending=False)
    out.to_csv(RESULTS / "Table_S7_rescue_posthoc_top_candidate_screen.csv", index=False)


if __name__ == "__main__":
    main()
